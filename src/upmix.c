#include "upmix.h"
#include "fx.h"
#include "thread.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

/* Samples live as -1..1 floats everywhere else in this program, but stage
 * boundaries round on a 32 bit integer grid, so the stage-round path has to
 * work in those units. The conversion is exact in both directions. */
static void lift(const float *src, size_t n, unsigned stride, unsigned off,
                 double *dst) {
    for (size_t i = 0; i < n; i++) dst[i] = (double)src[i * stride + off] * PS_FULL32;
}
/* Going back to floats needs care. The working values span the full 32 bit
 * integer range, and a float carries only 24 bits of mantissa, so the bottom
 * byte would be dropped on the way out and could tip the final rounding either
 * way. Putting the value on the 24 bit grid first makes the conversion exact,
 * since a 24 bit number does fit a float exactly. */
static void drop(const double *src, size_t n, float *dst, int q24) {
    if (q24) {
        for (size_t i = 0; i < n; i++)
            dst[i] = (float)(ps_regrid24(src[i]) / PS_FULL32);
    } else {
        for (size_t i = 0; i < n; i++) dst[i] = (float)(src[i] / PS_FULL32);
    }
}

/* One surround side: high band out of the stem, then blended onto the
 * surround. Both sides are identical work on different channels, so they run
 * side by side. */
typedef struct {
    const ps_audio *stem;
    ps_audio       *sur[2];      /* Ls, Rs - mono, overwritten in place */
    double          hz, gain_db, c_sur, c_high;
    int             flags;
    int             q24;
    int             failed;
} fs_job;

static void fs_side(int idx, void *ctx) {
    fs_job *j = (fs_job *)ctx;
    const ps_audio *stem = j->stem;
    ps_audio *sur = j->sur[idx];
    size_t n = sur->frames;
    if (stem->frames < n) n = stem->frames;

    double *hi = (double *)malloc(n * sizeof(double));
    double *lo = (double *)malloc(n * sizeof(double));
    if (!hi || !lo) { free(hi); free(lo); j->failed = 1; return; }

    /* The high band, taken from the matching channel of the original stem.
     * Two identical highpasses in series give a 12 dB per octave slope. */
    lift(stem->data, n, stem->channels, (unsigned)idx, hi);
    ps_biquad hp = ps_biquad_highpass(j->hz, (double)stem->rate);
    ps_biquad_state st = {0, 0, 0, 0};
    ps_biquad_run(&hp, &st, hi, n, j->flags);
    ps_biquad_state st2 = {0, 0, 0, 0};
    ps_biquad_run(&hp, &st2, hi, n, j->flags);
    ps_fx_gain(hi, n, j->gain_db, j->flags);

    /* Run as separate programs this band would sit in a 24 bit file before
     * the blend, so the comparison path loses the same bits at that point. */
    if (j->flags & PS_FX_STAGE_ROUND)
        for (size_t i = 0; i < n; i++) hi[i] = ps_regrid24(hi[i]);

    lift(sur->data, n, sur->channels, 0, lo);
    for (size_t i = 0; i < n; i++) {
        double v = lo[i] * j->c_sur + hi[i] * j->c_high;
        if (j->flags & PS_FX_STAGE_ROUND) v = ps_round32(v);
        lo[i] = v;
    }
    drop(lo, n, sur->data, j->q24);
    sur->frames = n;

    free(hi); free(lo);
}

int ps_freqsplit(const char *stemPath, const char *lsPath, const char *rsPath,
                 double hz, double gain_db, int blend_full,
                 int flags, ps_sample_fmt outfmt, const char **msg) {
    ps_audio stem = {0}, ls = {0}, rs = {0};
    int rc = 1;

    if (ps_wav_read(stemPath, &stem, msg)) goto done;
    if (ps_wav_read(lsPath, &ls, msg))     goto done;
    if (ps_wav_read(rsPath, &rs, msg))     goto done;
    if (stem.channels < 2) { *msg = "the stem must be stereo"; goto done; }
    if (ls.channels != 1 || rs.channels != 1) {
        *msg = "the surround channels must be mono"; goto done;
    }

    fs_job j;
    j.stem = &stem;
    j.sur[0] = &ls; j.sur[1] = &rs;
    j.hz = hz; j.gain_db = gain_db;
    j.c_sur  = blend_full ? 1.0 : 0.5;
    j.c_high = 0.5;
    j.flags = flags;
    j.q24 = (outfmt == PS_FMT_PCM24);
    j.failed = 0;

    ps_parallel_for(2, fs_side, &j, 2);
    if (j.failed) { *msg = "out of memory"; goto done; }

    if (ps_wav_write(lsPath, &ls, outfmt, msg)) goto done;
    if (ps_wav_write(rsPath, &rs, outfmt, msg)) goto done;
    rc = 0;
done:
    ps_audio_free(&stem); ps_audio_free(&ls); ps_audio_free(&rs);
    return rc;
}

int ps_decorrelate(const char *path, double delay_sec,
                   double ap1_hz, double ap2_hz, double ap_width,
                   int flags, ps_sample_fmt outfmt, const char **msg) {
    ps_audio a = {0};
    int rc = 1;
    double *buf = NULL;

    if (ps_wav_read(path, &a, msg)) return 1;
    if (a.channels != 1) { *msg = "decorrelate expects a mono file"; goto done; }

    size_t pad = ps_fx_delay_samples(delay_sec, a.rate);
    size_t n   = a.frames + pad;

    buf = (double *)malloc(n * sizeof(double));
    if (!buf) { *msg = "out of memory"; goto done; }
    for (size_t i = 0; i < pad; i++) buf[i] = 0.0;
    lift(a.data, a.frames, 1, 0, buf + pad);

    /* Both allpasses see the delayed signal, so the padding is in place first.
     * They leave the frequency balance alone and only move phase about, which
     * is what stops the two surrounds folding back together into mono. */
    /* A frequency of zero means "no filter here". Without this guard the
     * coefficients come out of a division by zero and the channel is filled
     * with noise rather than left alone - which is worse than an error,
     * because the result still plays. */
    if (ap1_hz > 0.0) {
        ps_biquad q1 = ps_biquad_allpass(ap1_hz, ap_width, (double)a.rate);
        ps_biquad_state s1 = {0, 0, 0, 0};
        ps_biquad_run(&q1, &s1, buf, n, flags);
    }
    if (ap2_hz > 0.0) {
        ps_biquad q2 = ps_biquad_allpass(ap2_hz, ap_width, (double)a.rate);
        ps_biquad_state s2 = {0, 0, 0, 0};
        ps_biquad_run(&q2, &s2, buf, n, flags);
    }

    float *out = (float *)malloc(n * sizeof(float));
    if (!out) { *msg = "out of memory"; goto done; }
    drop(buf, n, out, outfmt == PS_FMT_PCM24);
    free(a.data);
    a.data = out; a.frames = n;

    if (ps_wav_write(path, &a, outfmt, msg)) goto done;
    rc = 0;
done:
    free(buf);
    ps_audio_free(&a);
    return rc;
}

/* ---------------------------------------------------------------------------
 * Weighted sums of input channels.
 * ------------------------------------------------------------------------- */

int ps_mix_matrix(const char *map, int nch, double *m, int maxOut, const char **msg) {
    if (nch < 1 || nch > PS_MIX_MAX_CH) { *msg = "too many input channels"; return -1; }
    for (int i = 0; i < maxOut * nch; i++) m[i] = 0.0;

    int no = 0;
    const char *p = map;
    while (*p) {
        if (no >= maxOut) { *msg = "too many output channels"; return -1; }
        double *row = m + (size_t)no * nch;
        for (;;) {
            while (*p == ' ') p++;
            char *end;
            long src = strtol(p, &end, 10);
            if (end == p) { *msg = "expected a channel number"; return -1; }
            if (src < 1 || src > nch) { *msg = "the map names a channel that does not exist"; return -1; }
            p = end;
            double g = 1.0;
            while (*p == ' ') p++;
            if (*p == '*') {
                p++;
                g = strtod(p, &end);
                if (end == p) { *msg = "expected a number after '*'"; return -1; }
                p = end;
            }
            /* Adding rather than assigning, so naming a channel twice sums it
             * twice instead of quietly keeping only the last mention. */
            row[src - 1] += g;
            while (*p == ' ') p++;
            if (*p == '+') { p++; continue; }
            break;
        }
        no++;
        if (*p == ';') {
            p++;
            /* A map ending in ';' is almost certainly one that got cut short.
             * Taking it as "that many outputs and no more" would quietly drop a
             * channel, so say so instead. */
            if (!*p) { *msg = "the map ends with a ';' and no channel after it"; return -1; }
            continue;
        }
        if (*p) { *msg = "unexpected character in the map"; return -1; }
    }
    if (no == 0) { *msg = "the map is empty"; return -1; }
    return no;
}

typedef struct {
    const ps_audio *in;
    const int      *file;    /* per global channel: which file it came from */
    const int      *local;   /* per global channel: its index within that file */
    int             nch;
    const double   *m;       /* nOut rows of nch weights */
    int             nOut;
    size_t          frames;
    float          *dst;
    double          gain;
    const ps_filt  *filts;
    int             nFilt;
    unsigned        rate;
    int             flags, q24;
    int             failed;
} mix_job;

static void mix_one(int oc, void *ctx) {
    mix_job *j = (mix_job *)ctx;
    const double *row = j->m + (size_t)oc * j->nch;

    /* Only the channels with a weight are worth touching, so gather them once
     * rather than walking the whole row for every sample. */
    int idx[PS_MIX_MAX_CH], n = 0;
    for (int c = 0; c < j->nch; c++) if (row[c] != 0.0) idx[n++] = c;

    /* Filters have to see the channel in order, so the matrix result is built
     * out in full before they run. Without filters there is nothing to hold
     * and the samples can go straight out. */
    double *tmp = NULL;
    if (j->nFilt > 0) {
        tmp = (double *)malloc(j->frames * sizeof(double));
        if (!tmp) { j->failed = 1; return; }
    }

    for (size_t i = 0; i < j->frames; i++) {
        double acc = 0.0;
        for (int k = 0; k < n; k++) {
            int c = idx[k];
            const ps_audio *a = &j->in[j->file[c]];
            /* A file shorter than the longest counts as silence past its end,
             * so one short input cannot cut the result off. */
            double v = (i < a->frames)
                     ? (double)a->data[i * a->channels + j->local[c]] * PS_FULL32
                     : 0.0;
            acc += v * row[c];
        }
        if (j->flags & PS_FX_STAGE_ROUND) acc = ps_round32(acc);
        if (tmp) tmp[i] = acc;
        else {
            if (j->gain != 1.0) {
                acc *= j->gain;
                if (j->flags & PS_FX_STAGE_ROUND) acc = ps_round32(acc);
            }
            if (j->q24) acc = ps_regrid24(acc);
            j->dst[i * (size_t)j->nOut + oc] = (float)(acc / PS_FULL32);
        }
    }

    if (tmp) {
        for (int f = 0; f < j->nFilt; f++) {
            ps_biquad q;
            switch (j->filts[f].kind) {
                case PS_FILT_LOWPASS:
                    q = ps_biquad_lowpass(j->filts[f].hz, (double)j->rate); break;
                case PS_FILT_HIGHSHELF:
                    q = ps_biquad_highshelf(j->filts[f].db, j->filts[f].hz,
                                            (double)j->rate); break;
                default:
                    q = ps_biquad_highpass(j->filts[f].hz, (double)j->rate); break;
            }
            ps_biquad_state st = {0, 0, 0, 0};
            ps_biquad_run(&q, &st, tmp, j->frames, j->flags);
        }
        for (size_t i = 0; i < j->frames; i++) {
            double v = tmp[i];
            if (j->gain != 1.0) {
                v *= j->gain;
                if (j->flags & PS_FX_STAGE_ROUND) v = ps_round32(v);
            }
            if (j->q24) v = ps_regrid24(v);
            j->dst[i * (size_t)j->nOut + oc] = (float)(v / PS_FULL32);
        }
        free(tmp);
    }
}

int ps_mix(const char **inPaths, int nIn, const char *outPath,
           const char *map, double gain_db,
           const ps_filt *filts, int nFilt, size_t trim,
           int flags, ps_sample_fmt outfmt, const char **msg) {
    ps_audio *in = NULL;
    double *m = NULL;
    float *dst = NULL;
    int file[PS_MIX_MAX_CH], local[PS_MIX_MAX_CH];
    int rc = 1, loaded = 0;

    if (nIn < 1) { *msg = "no input files"; return 1; }
    in = (ps_audio *)calloc((size_t)nIn, sizeof(ps_audio));
    m  = (double *)malloc((size_t)PS_MIX_MAX_OUT * PS_MIX_MAX_CH * sizeof(double));
    if (!in || !m) { *msg = "out of memory"; goto done; }

    size_t frames = 0;
    unsigned rate = 0;
    int nch = 0;
    for (int i = 0; i < nIn; i++) {
        if (ps_wav_read(inPaths[i], &in[i], msg)) goto done;
        loaded = i + 1;
        if (i == 0) rate = in[i].rate;
        else if (in[i].rate != rate) { *msg = "the inputs have different sample rates"; goto done; }
        for (unsigned c = 0; c < in[i].channels; c++) {
            if (nch >= PS_MIX_MAX_CH) { *msg = "too many input channels in total"; goto done; }
            file[nch] = i; local[nch] = (int)c; nch++;
        }
        if (in[i].frames > frames) frames = in[i].frames;
    }

    int nOut = ps_mix_matrix(map, nch, m, PS_MIX_MAX_OUT, msg);
    if (nOut < 0) goto done;

    dst = (float *)malloc(frames * (size_t)nOut * sizeof(float));
    if (!dst) { *msg = "out of memory"; goto done; }

    mix_job j;
    j.in = in; j.file = file; j.local = local; j.nch = nch;
    j.m = m; j.nOut = nOut;
    j.frames = frames; j.dst = dst;
    j.gain = pow(10.0, gain_db / 20.0);
    j.filts = filts; j.nFilt = nFilt;
    j.rate = rate;
    j.flags = flags;
    j.q24 = (outfmt == PS_FMT_PCM24);
    j.failed = 0;

    ps_parallel_for(nOut, mix_one, &j, nOut);
    if (j.failed) { *msg = "out of memory"; goto done; }

    /* Trimming is only ever used to take back off the padding that was added
     * so the overlapping transform had something to chew on at the ends. */
    if (trim > 0 && trim < frames) frames = trim;

    ps_audio out;
    out.data = dst; out.frames = frames;
    out.channels = (unsigned)nOut; out.rate = rate;
    out.src_fmt = in[0].src_fmt;
    if (ps_wav_write(outPath, &out, outfmt, msg)) goto done;
    rc = 0;
done:
    for (int i = 0; i < loaded; i++) ps_audio_free(&in[i]);
    free(in); free(m); free(dst);
    return rc;
}
