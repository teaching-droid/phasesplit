#include "split.h"
#include "fft.h"
#include "dsp.h"
#include "thread.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Scratch space for one worker. Each thread needs its own set, since all of
 * this is written while a frame is being processed. */
typedef struct {
    ps_fft *fft;
    float  *block;                 /* one frame of interleaved input */
    float  *frameL, *frameR;
    float  *Lre, *Lim, *Rre, *Rim;
    float  *Cre, *Cim;
    float  *outC;
} ps_worker;

struct ps_split {
    ps_split_cfg cfg;
    int    hop;
    int    bins;
    int    nworkers;

    float *win;        /* analysis window                                    */
    float *wout;       /* synthesis window with the overlap weight folded in */

    ps_worker *w;
};

void ps_split_default_cfg(ps_split_cfg *cfg) {
    if (!cfg) return;
    cfg->window      = 4096;
    cfg->overlap     = 4;
    cfg->strength    = 1.0f;
    cfg->dual_centre = 0;
    cfg->threads     = 0;          /* 0 = one per hardware thread */
}

static void worker_free(ps_worker *w) {
    if (!w) return;
    ps_fft_destroy(w->fft);
    free(w->block); free(w->frameL); free(w->frameR);
    free(w->Lre); free(w->Lim); free(w->Rre); free(w->Rim);
    free(w->Cre); free(w->Cim); free(w->outC);
    memset(w, 0, sizeof(*w));
}

static int worker_init(ps_worker *w, int N, int bins) {
    memset(w, 0, sizeof(*w));
    w->fft    = ps_fft_create(N);
    w->block  = (float *)malloc((size_t)N * 2 * sizeof(float));
    w->frameL = (float *)malloc((size_t)N * sizeof(float));
    w->frameR = (float *)malloc((size_t)N * sizeof(float));
    w->outC   = (float *)malloc((size_t)N * sizeof(float));
    w->Lre = (float *)malloc((size_t)bins * sizeof(float));
    w->Lim = (float *)malloc((size_t)bins * sizeof(float));
    w->Rre = (float *)malloc((size_t)bins * sizeof(float));
    w->Rim = (float *)malloc((size_t)bins * sizeof(float));
    w->Cre = (float *)malloc((size_t)bins * sizeof(float));
    w->Cim = (float *)malloc((size_t)bins * sizeof(float));
    if (!w->fft || !w->block || !w->frameL || !w->frameR || !w->outC ||
        !w->Lre || !w->Lim || !w->Rre || !w->Rim || !w->Cre || !w->Cim) {
        worker_free(w); return 1;
    }
    return 0;
}

ps_split *ps_split_create(const ps_split_cfg *cfg) {
    if (!cfg) return NULL;
    if (cfg->window < 64 || (cfg->window & (cfg->window - 1)) != 0) return NULL;
    if (cfg->overlap < 2 || (cfg->overlap & (cfg->overlap - 1)) != 0) return NULL;
    if (cfg->overlap > cfg->window) return NULL;

    /* The caller decides the instruction set before getting here, so this only
     * builds the table if nothing has chosen one yet. */
    ps_dsp_init_once();

    ps_split *s = (ps_split *)calloc(1, sizeof(ps_split));
    if (!s) return NULL;
    s->cfg  = *cfg;
    s->hop  = cfg->window / cfg->overlap;
    s->bins = cfg->window / 2 + 1;

    int nt = cfg->threads > 0 ? cfg->threads : ps_cpu_threads();
    if (nt < 1)   nt = 1;
    if (nt > 256) nt = 256;
    s->nworkers = nt;

    const int N = cfg->window;
    s->win  = (float *)malloc((size_t)N * sizeof(float));
    s->wout = (float *)malloc((size_t)N * sizeof(float));
    s->w    = (ps_worker *)calloc((size_t)nt, sizeof(ps_worker));
    if (!s->win || !s->wout || !s->w) { ps_split_destroy(s); return NULL; }

    for (int i = 0; i < nt; i++)
        if (worker_init(&s->w[i], N, s->bins) != 0) { ps_split_destroy(s); return NULL; }

    for (int i = 0; i < N; i++)
        s->win[i] = (float)(0.5 * (1.0 - cos(2.0 * M_PI * (double)i / (double)N)));

    /* The window goes on twice, once in and once out, so what has to come back
     * flat is the sum of w^2 over the overlapping frames. For a Hann window
     * that holds at 75% overlap but not at 50%, so it is worked out per
     * position rather than assumed, then folded into the synthesis window to
     * keep it out of the inner loop. */
    for (int i = 0; i < s->hop; i++) {
        double acc = 0.0;
        for (int k = 0; i + k * s->hop < N; k++) {
            double w = s->win[i + k * s->hop];
            acc += w * w;
        }
        const float inv = (acc > 1e-12) ? (float)(1.0 / acc) : 1.0f;
        for (int k = 0; i + k * s->hop < N; k++)
            s->wout[i + k * s->hop] = s->win[i + k * s->hop] * inv;
    }

    return s;
}

void ps_split_destroy(ps_split *s) {
    if (!s) return;
    if (s->w) { for (int i = 0; i < s->nworkers; i++) worker_free(&s->w[i]); free(s->w); }
    free(s->win); free(s->wout);
    free(s);
}

int ps_split_latency(const ps_split *s) { return s ? s->cfg.window : 0; }
int ps_split_threads(const ps_split *s) { return s ? s->nworkers : 0; }

typedef struct {
    ps_split    *s;
    const float *in;
    size_t       frames;
    float       *accC, *accL, *accR;
    size_t       chunk;        /* output samples per slice */
    int          wantSides;
} job;

/* Each slice owns a stretch of the output and nothing else writes there, so
 * the workers never have to synchronise. A frame that straddles a boundary is
 * computed by both neighbours, each keeping only the part inside its own
 * slice - that repeats a little work but removes locking from the hot path. */
static void do_chunk(int idx, void *vctx) {
    job *j = (job *)vctx;
    ps_split *s = j->s;
    ps_worker *w = &s->w[idx];

    const int N   = s->cfg.window;
    const int hop = s->hop;

    const long long c0 = (long long)idx * (long long)j->chunk;
    long long c1 = c0 + (long long)j->chunk;
    if (c1 > (long long)j->frames) c1 = (long long)j->frames;
    if (c0 >= c1) return;

    /* Frame k covers output [-N + k*hop, k*hop), so these are the frames that
     * reach into this slice at all. */
    const long long klo = c0 / hop;
    const long long khi = (c1 + N) / hop + 1;

    for (long long k = klo; k <= khi; k++) {
        const long long start = -(long long)N + k * (long long)hop;
        if (start >= (long long)j->frames) break;

        /* the part of the frame that is inside the file */
        const long long lo = start < 0 ? -start : 0;
        const long long hi = (start + N > (long long)j->frames)
                           ? ((long long)j->frames - start) : N;
        if (hi <= lo) continue;

        /* ...and the part of that which belongs to this slice */
        long long wlo = lo, whi = hi;
        if (start + wlo < c0) wlo = c0 - start;
        if (start + whi > c1) whi = c1 - start;
        if (whi <= wlo) continue;

        if (lo > 0) memset(w->block, 0, (size_t)lo * 2 * sizeof(float));
        if (hi < N) memset(w->block + hi * 2, 0, (size_t)(N - hi) * 2 * sizeof(float));
        memcpy(w->block + lo * 2, j->in + (size_t)(start + lo) * 2,
               (size_t)(hi - lo) * 2 * sizeof(float));

        ps_window_split(w->block, s->win, w->frameL, w->frameR, N);
        ps_fft_forward(w->fft, w->frameL, w->Lre, w->Lim);
        ps_fft_forward(w->fft, w->frameR, w->Rre, w->Rim);
        ps_estimate_centre(w->Lre, w->Lim, w->Rre, w->Rim, w->Cre, w->Cim,
                           s->bins, s->cfg.strength);

        /* One inverse transform, not three. The sides are the frame minus the
         * centre and the inverse transform is linear, so that subtraction can
         * be done in the time domain against the windowed frame already held. */
        ps_fft_inverse(w->fft, w->Cre, w->Cim, w->outC);

        const int    n   = (int)(whi - wlo);
        const size_t dst = (size_t)(start + wlo);
        if (j->accC) ps_scale_add(j->accC + dst, w->outC + wlo, s->wout + wlo, 1.0f, n);
        if (j->wantSides) {
            ps_scale_add_diff(j->accL + dst, w->frameL + wlo, w->outC + wlo, s->wout + wlo, n);
            ps_scale_add_diff(j->accR + dst, w->frameR + wlo, w->outC + wlo, s->wout + wlo, n);
        }
    }
}

int ps_split_process(ps_split *s, const float *in, size_t frames,
                     float *centre, float *sides) {
    if (!s || !in || frames == 0) return 1;

    const int cch = s->cfg.dual_centre ? 2 : 1;
    float *accC = NULL, *accL = NULL, *accR = NULL;
    int    ownC = 0;

    if (centre) {
        if (cch == 1) { accC = centre; memset(accC, 0, frames * sizeof(float)); }
        else { accC = (float *)calloc(frames, sizeof(float)); ownC = 1; if (!accC) return 1; }
    }
    if (sides) {
        accL = (float *)calloc(frames, sizeof(float));
        accR = (float *)calloc(frames, sizeof(float));
        if (!accL || !accR) { if (ownC) free(accC); free(accL); free(accR); return 1; }
    }

    /* Keep each slice comfortably longer than a window, or the shared boundary
     * work starts to cost more than the extra thread is worth. */
    int nchunks = s->nworkers;
    const size_t minChunk = (size_t)s->cfg.window * 8;
    while (nchunks > 1 && frames / (size_t)nchunks < minChunk) nchunks--;

    job j;
    j.s = s; j.in = in; j.frames = frames;
    j.accC = accC; j.accL = accL; j.accR = accR;
    j.wantSides = sides ? 1 : 0;
    j.chunk = (frames + (size_t)nchunks - 1) / (size_t)nchunks;

    ps_parallel_for(nchunks, do_chunk, &j, nchunks);

    if (sides) {
        for (size_t t = 0; t < frames; t++) {
            sides[t * 2]     = accL[t];
            sides[t * 2 + 1] = accR[t];
        }
        free(accL); free(accR);
    }
    if (centre && cch == 2)
        for (size_t t = 0; t < frames; t++) centre[t * 2] = centre[t * 2 + 1] = accC[t];
    if (ownC) free(accC);
    return 0;
}
