/* Built-in checks. These run from the binary itself (phasesplit --selftest)
 * so the maths can be verified on the machine that will actually run it,
 * without needing a build of the test tooling.
 */
#include "fft.h"
#include "split.h"
#include "cpu.h"
#include "dsp.h"
#include "fx.h"
#include "upmix.h"
#include "version.h"
#include "selftest.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static int g_pass, g_fail;

static void ok(const char *name, int cond, const char *detail) {
    if (cond) { g_pass++; printf("  PASS  %-46s %s\n", name, detail ? detail : ""); }
    else      { g_fail++; printf("  FAIL  %-46s %s\n", name, detail ? detail : ""); }
}

/* Straight O(n^2) DFT in double, used only as the yardstick the fast
 * transform is checked against. Slow on purpose: it is obviously correct. */
static void naive_dft(const float *x, int n, int bin, double *re, double *im) {
    double sr = 0.0, si = 0.0;
    for (int t = 0; t < n; t++) {
        double a = -2.0 * M_PI * (double)bin * (double)t / (double)n;
        sr += (double)x[t] * cos(a);
        si += (double)x[t] * sin(a);
    }
    *re = sr; *im = si;
}

static int test_fft(int n) {
    char label[96];
    ps_fft *p = ps_fft_create(n);
    if (!p) { ok("fft create", 0, "returned NULL"); return 1; }

    const int bins = ps_fft_bins(p);
    float *x   = (float *)malloc((size_t)n * sizeof(float));
    float *re  = (float *)malloc((size_t)bins * sizeof(float));
    float *im  = (float *)malloc((size_t)bins * sizeof(float));
    float *rt  = (float *)malloc((size_t)n * sizeof(float));
    if (!x || !re || !im || !rt) { ok("fft alloc", 0, NULL); return 1; }

    /* 1. an impulse has a flat spectrum of magnitude one */
    memset(x, 0, (size_t)n * sizeof(float));
    x[0] = 1.0f;
    ps_fft_forward(p, x, re, im);
    double worst = 0.0;
    for (int k = 0; k < bins; k++) {
        double m = sqrt((double)re[k] * re[k] + (double)im[k] * im[k]);
        double e = fabs(m - 1.0);
        if (e > worst) worst = e;
    }
    snprintf(label, sizeof(label), "n=%d impulse gives a flat spectrum", n);
    { char d[64]; snprintf(d, sizeof(d), "max error %.2e", worst); ok(label, worst < 1e-5, d); }

    /* 2. a constant signal puts everything in bin zero */
    for (int t = 0; t < n; t++) x[t] = 1.0f;
    ps_fft_forward(p, x, re, im);
    double leak = 0.0;
    for (int k = 1; k < bins; k++) {
        double m = sqrt((double)re[k] * re[k] + (double)im[k] * im[k]);
        if (m > leak) leak = m;
    }
    snprintf(label, sizeof(label), "n=%d DC lands only in bin 0", n);
    { char d[80]; snprintf(d, sizeof(d), "bin0 %.1f, worst other %.2e", (double)re[0], leak);
      ok(label, fabs((double)re[0] - n) < 1e-2 * n && leak < 1e-3, d); }

    /* 3. a sine exactly on a bin centre lands in that bin */
    const int kb = n / 8;
    for (int t = 0; t < n; t++) x[t] = (float)sin(2.0 * M_PI * kb * t / n);
    ps_fft_forward(p, x, re, im);
    double onbin = sqrt((double)re[kb] * re[kb] + (double)im[kb] * im[kb]);
    double offbin = 0.0;
    for (int k = 0; k < bins; k++) {
        if (k == kb) continue;
        double m = sqrt((double)re[k] * re[k] + (double)im[k] * im[k]);
        if (m > offbin) offbin = m;
    }
    snprintf(label, sizeof(label), "n=%d on-bin sine is isolated", n);
    { char d[80]; snprintf(d, sizeof(d), "in bin %.1f, worst leak %.2e", onbin, offbin);
      ok(label, onbin > 0.4 * n && offbin < 1e-2, d); }

    /* 4. against the slow reference on an awkward signal */
    srand(12345);
    for (int t = 0; t < n; t++) x[t] = (float)((rand() / (double)RAND_MAX) * 2.0 - 1.0);
    ps_fft_forward(p, x, re, im);
    double maxrel = 0.0;
    for (int k = 0; k <= 6; k++) {
        int bin = (k * (bins - 1)) / 6;
        double dr, di;
        naive_dft(x, n, bin, &dr, &di);
        double dm = sqrt(dr * dr + di * di);
        double er = fabs(dr - (double)re[bin]);
        double ei = fabs(di - (double)im[bin]);
        double rel = (er + ei) / (dm > 1e-9 ? dm : 1.0);
        if (rel > maxrel) maxrel = rel;
    }
    snprintf(label, sizeof(label), "n=%d matches a direct DFT", n);
    { char d[64]; snprintf(d, sizeof(d), "worst relative %.2e", maxrel); ok(label, maxrel < 1e-4, d); }

    /* 5. forward then inverse returns the original */
    ps_fft_forward(p, x, re, im);
    ps_fft_inverse(p, re, im, rt);
    double rterr = 0.0;
    for (int t = 0; t < n; t++) {
        double e = fabs((double)rt[t] - (double)x[t]);
        if (e > rterr) rterr = e;
    }
    snprintf(label, sizeof(label), "n=%d round trip is lossless", n);
    { char d[64]; snprintf(d, sizeof(d), "max error %.2e", rterr); ok(label, rterr < 1e-4, d); }

    /* 6. Parseval: energy is conserved between the two domains */
    double et = 0.0, ef = 0.0;
    for (int t = 0; t < n; t++) et += (double)x[t] * x[t];
    for (int k = 0; k < bins; k++) {
        double m2 = (double)re[k] * re[k] + (double)im[k] * im[k];
        ef += (k == 0 || k == bins - 1) ? m2 : 2.0 * m2;
    }
    ef /= (double)n;
    double perr = fabs(et - ef) / (et > 1e-9 ? et : 1.0);
    snprintf(label, sizeof(label), "n=%d conserves energy", n);
    { char d[64]; snprintf(d, sizeof(d), "mismatch %.2e", perr); ok(label, perr < 1e-4, d); }

    free(x); free(re); free(im); free(rt);
    ps_fft_destroy(p);
    return 0;
}


/* Build a stereo test signal: a centred component plus two hard panned ones,
 * all on different frequencies so the result can be checked bin by bin. */
static void make_mix(float *st, int frames, int rate,
                     double fC, double fL, double fR,
                     double aC, double aL, double aR) {
    for (int t = 0; t < frames; t++) {
        double c = aC * sin(2.0 * M_PI * fC * t / rate);
        double l = aL * sin(2.0 * M_PI * fL * t / rate + 0.7);
        double r = aR * sin(2.0 * M_PI * fR * t / rate + 1.9);
        st[t * 2]     = (float)(c + l);
        st[t * 2 + 1] = (float)(c + r);
    }
}

static double band_rms(const float *x, int n, int stride, int rate, double f0) {
    /* Energy near one frequency, found by correlating against a complex tone.
     * Enough to say how much of a known component ended up in a channel. */
    double sr = 0.0, si = 0.0;
    for (int t = 0; t < n; t++) {
        double a = 2.0 * M_PI * f0 * t / rate;
        sr += (double)x[t * stride] * cos(a);
        si += (double)x[t * stride] * sin(a);
    }
    return 2.0 * sqrt(sr * sr + si * si) / (double)n;
}

static void test_split_reconstruction(int window, int overlap) {
    char label[96], detail[96];
    const int rate = 44100, frames = 44100;

    ps_split_cfg cfg;
    ps_split_default_cfg(&cfg);
    cfg.window = window; cfg.overlap = overlap;

    ps_split *s = ps_split_create(&cfg);
    if (!s) { ok("split create", 0, "returned NULL"); return; }

    float *in  = (float *)malloc((size_t)frames * 2 * sizeof(float));
    float *cen = (float *)malloc((size_t)frames * sizeof(float));
    float *sid = (float *)malloc((size_t)frames * 2 * sizeof(float));
    if (!in || !cen || !sid) { ok("split alloc", 0, NULL); return; }

    make_mix(in, frames, rate, 700.0, 2200.0, 3300.0, 0.30, 0.25, 0.25);
    ps_split_process(s, in, (size_t)frames, cen, sid);

    /* The whole point of taking the sides as a remainder: putting the parts
     * back together has to give the input again, sample for sample. Any error
     * here means the windowing or the overlap-add is wrong. */
    double worst = 0.0;
    for (int t = window; t < frames - window; t++) {
        double l = (double)cen[t] + (double)sid[t * 2];
        double r = (double)cen[t] + (double)sid[t * 2 + 1];
        double el = fabs(l - (double)in[t * 2]);
        double er = fabs(r - (double)in[t * 2 + 1]);
        if (el > worst) worst = el;
        if (er > worst) worst = er;
    }
    snprintf(label, sizeof(label), "w=%d o=%dx centre+sides rebuilds the input", window, overlap);
    snprintf(detail, sizeof(detail), "max error %.2e", worst);
    ok(label, worst < 1e-4, detail);

    /* Separation: the centred tone should be in the centre output and the
     * panned tones should not. */
    int mid = frames / 4, span = frames / 2;
    double cIn  = band_rms(in  + (size_t)mid * 2, span, 2, rate, 700.0);
    double cOut = band_rms(cen + mid,             span, 1, rate, 700.0);
    double lLeak = band_rms(cen + mid,            span, 1, rate, 2200.0);
    double cKept = 20.0 * log10((cOut > 1e-12 ? cOut : 1e-12) / (cIn > 1e-12 ? cIn : 1e-12));
    double lRej  = 20.0 * log10((lLeak > 1e-12 ? lLeak : 1e-12) / (cIn > 1e-12 ? cIn : 1e-12));

    snprintf(label, sizeof(label), "w=%d o=%dx keeps the centred tone", window, overlap);
    snprintf(detail, sizeof(detail), "%.1f dB of it", cKept);
    ok(label, cKept > -3.0, detail);

    snprintf(label, sizeof(label), "w=%d o=%dx rejects a panned tone", window, overlap);
    snprintf(detail, sizeof(detail), "%.1f dB leak", lRej);
    ok(label, lRej < -20.0, detail);

    free(in); free(cen); free(sid);
    ps_split_destroy(s);
}


/* ---------------------------------------------------------------------------
 * Parameter sweep.
 *
 * A long transform separates steady tones better, because each bin is
 * narrower. The same length smears anything sudden, because one frame covers
 * more time and the overlap-add spreads a click across the whole window - some
 * of it landing before the click happened, which is what pre-echo means.
 * Neither number alone picks a window length, so both get measured here and
 * the defaults are chosen from the table rather than from habit.
 * ------------------------------------------------------------------------- */

static double measure_rejection(int window, int overlap) {
    /* Four seconds: a 16384 window is 0.37 s long, so a one second
     * signal would be measuring its own edges more than its separation. */
    const int rate = 44100, frames = 44100 * 4;
    ps_split_cfg cfg; ps_split_default_cfg(&cfg);
    cfg.window = window; cfg.overlap = overlap;
    ps_split *s = ps_split_create(&cfg);
    if (!s) return 0.0;

    float *in  = (float *)malloc((size_t)frames * 2 * sizeof(float));
    float *cen = (float *)malloc((size_t)frames * sizeof(float));
    if (!in || !cen) { ps_split_destroy(s); free(in); free(cen); return 0.0; }

    make_mix(in, frames, rate, 700.0, 2200.0, 3300.0, 0.30, 0.25, 0.25);
    ps_split_process(s, in, (size_t)frames, cen, NULL);

    int mid = frames / 4, span = frames / 2;
    double ref  = band_rms(in + (size_t)mid * 2, span, 2, rate, 700.0);
    double leak = band_rms(cen + mid,            span, 1, rate, 2200.0);
    (void)window;
    double db = 20.0 * log10((leak > 1e-14 ? leak : 1e-14) / (ref > 1e-14 ? ref : 1e-14));

    free(in); free(cen); ps_split_destroy(s);
    return db;
}

static double measure_smear(int window, int overlap, double *peak_db) {
    /* Pre-onset leakage, which is the honest cost of a long transform.
     *
     * The signal holds one frequency throughout. Before the switch it sits in
     * the left channel only, so none of it belongs in the centre. After the
     * switch the same frequency is centred, so all of it does. One gain is
     * worked out per frame and applied across that whole frame, so the frame
     * straddling the switch applies some of the "keep it" gain to the part
     * that came before - and centre output appears where there should be
     * none. The longer the window, the further back that reaches.
     *
     * Returned as how far back the leak stays audible, in ms, with its level
     * relative to the wanted signal alongside. */
    const int rate = 44100, frames = rate * 3;
    const int sw = frames / 2;
    const double f = 1000.0;

    ps_split_cfg cfg; ps_split_default_cfg(&cfg);
    cfg.window = window; cfg.overlap = overlap;
    ps_split *s = ps_split_create(&cfg);
    if (!s) { if (peak_db) *peak_db = 0.0; return 0.0; }

    float *in  = (float *)malloc((size_t)frames * 2 * sizeof(float));
    float *cen = (float *)malloc((size_t)frames * sizeof(float));
    if (!in || !cen) { ps_split_destroy(s); free(in); free(cen); if (peak_db) *peak_db = 0.0; return 0.0; }

    for (int t = 0; t < frames; t++) {
        double v = 0.4 * sin(2.0 * M_PI * f * t / rate);
        if (t < sw) { in[t * 2] = (float)v; in[t * 2 + 1] = 0.0f; }   /* left only */
        else        { in[t * 2] = (float)v; in[t * 2 + 1] = (float)v; } /* centred  */
    }
    ps_split_process(s, in, (size_t)frames, cen, NULL);

    /* Reference: how loud the centre is once the signal really is centred. */
    double post = 0.0;
    int pn = 0;
    for (int t = sw + rate / 4; t < sw + rate / 2; t++) { post += (double)cen[t] * cen[t]; pn++; }
    post = sqrt(post / (pn ? pn : 1));

    /* Walk backwards from the switch and find how far the leak carries. */
    const double thr = post * 0.01;          /* one hundredth, about -40 dB */
    int reach = 0;
    const int win = rate / 200;              /* 5 ms measuring blocks */
    for (int b = 1; b * win < sw; b++) {
        double e = 0.0;
        for (int t = sw - b * win; t < sw - (b - 1) * win; t++) e += (double)cen[t] * cen[t];
        e = sqrt(e / win);
        if (e >= thr) reach = b; else if (b > 4) break;
    }
    double ms = 1000.0 * (double)(reach * win) / (double)rate;

    /* Loudest leak in the 100 ms before the switch, against the wanted level */
    double peak = 0.0;
    for (int t = sw - rate / 10; t < sw; t++) {
        double a = fabs((double)cen[t]);
        if (a > peak) peak = a;
    }
    if (peak_db) *peak_db = 20.0 * log10((peak > 1e-14 ? peak : 1e-14) / (post > 1e-14 ? post : 1e-14));

    free(in); free(cen); ps_split_destroy(s);
    return ms;
}

static double measure_speed(int window, int overlap) {
    const int rate = 44100, frames = rate * 10;   /* ten seconds of audio */
    ps_split_cfg cfg; ps_split_default_cfg(&cfg);
    cfg.window = window; cfg.overlap = overlap;
    ps_split *s = ps_split_create(&cfg);
    if (!s) return 0.0;

    float *in  = (float *)malloc((size_t)frames * 2 * sizeof(float));
    float *cen = (float *)malloc((size_t)frames * sizeof(float));
    float *sid = (float *)malloc((size_t)frames * 2 * sizeof(float));
    if (!in || !cen || !sid) { ps_split_destroy(s); free(in); free(cen); free(sid); return 0.0; }

    srand(999);
    for (int t = 0; t < frames * 2; t++) in[t] = (float)((rand() / (double)RAND_MAX) - 0.5);

    clock_t t0 = clock();
    ps_split_process(s, in, (size_t)frames, cen, sid);
    clock_t t1 = clock();

    free(in); free(cen); free(sid); ps_split_destroy(s);
    return 1000.0 * (double)(t1 - t0) / (double)CLOCKS_PER_SEC;
}

int ps_sweep(void) {
    static const int windows[]  = { 1024, 2048, 4096, 8192, 16384 };
    static const int overlaps[] = { 2, 4, 8 };

    printf("phasesplit parameter sweep\n\n");
    printf("  rejection  how far a steady panned tone is pushed out of the centre (lower is better)\n");
    printf("  smear      how long a panned click leaks into the centre for (lower is better)\n");
    printf("  leak       loudest part of that leak, against the click itself\n");
    printf("  time       to process 10 s of stereo, centre and sides\n\n");
    printf("  %8s %8s %12s %10s %10s %9s\n",
           "window", "overlap", "rejection", "smear", "leak", "time");
    printf("  %8s %8s %12s %10s %10s %9s\n",
           "------", "-------", "---------", "-----", "----", "----");

    for (size_t wi = 0; wi < sizeof(windows) / sizeof(windows[0]); wi++) {
        for (size_t oi = 0; oi < sizeof(overlaps) / sizeof(overlaps[0]); oi++) {
            int w = windows[wi], o = overlaps[oi];
            double leak = 0.0;
            double rej  = measure_rejection(w, o);
            double sm   = measure_smear(w, o, &leak);
            double ms   = measure_speed(w, o);
            printf("  %8d %7dx %10.1f dB %7.0f ms %7.1f dB %6.0f ms\n",
                   w, o, rej, sm, leak, ms);
        }
    }
    printf("\n  one window covers 23 ms at 1024, 93 ms at 4096, 186 ms at 8192, 372 ms at 16384\n");
    return 0;
}


/* Every instruction set path has to give the same answer. A fast path that
 * quietly differs is worse than no fast path, because the difference only
 * turns up on someone else's machine. Each available path processes the same
 * signal and the results are compared sample for sample against scalar. */
static void test_simd_agreement(void) {
    const int rate = 44100, frames = rate / 2;
    const ps_isa native = ps_cpu_isa();

    float *in   = (float *)malloc((size_t)frames * 2 * sizeof(float));
    float *cenS = (float *)malloc((size_t)frames * sizeof(float));
    float *sidS = (float *)malloc((size_t)frames * 2 * sizeof(float));
    float *cenX = (float *)malloc((size_t)frames * sizeof(float));
    float *sidX = (float *)malloc((size_t)frames * 2 * sizeof(float));
    if (!in || !cenS || !sidS || !cenX || !sidX) { ok("simd alloc", 0, NULL); return; }

    /* Awkward on purpose: a centred tone, two panned tones and some noise, so
     * the bin arithmetic hits its clamps and its near-silent guard as well. */
    srand(4242);
    for (int t = 0; t < frames; t++) {
        double cc = 0.30 * sin(2.0 * M_PI *  700.0 * t / rate);
        double ll = 0.20 * sin(2.0 * M_PI * 2200.0 * t / rate + 0.7);
        double rr = 0.20 * sin(2.0 * M_PI * 3300.0 * t / rate + 1.9);
        double nz = 0.02 * ((rand() / (double)RAND_MAX) - 0.5);
        in[t * 2]     = (float)(cc + ll + nz);
        in[t * 2 + 1] = (float)(cc + rr - nz);
    }

    ps_split_cfg cfg; ps_split_default_cfg(&cfg);

    ps_cpu_force(PS_ISA_SCALAR);
    ps_dsp_init();
    {
        ps_split *s = ps_split_create(&cfg);
        ps_split_process(s, in, (size_t)frames, cenS, sidS);
        ps_split_destroy(s);
    }

    static const ps_isa paths[] = { PS_ISA_SSE2, PS_ISA_AVX2 };
    for (size_t i = 0; i < sizeof(paths) / sizeof(paths[0]); i++) {
        char label[96], detail[96];
        if (paths[i] > native) {
            snprintf(label, sizeof(label), "%s path", ps_isa_name(paths[i]));
            ok(label, 1, "not supported by this CPU, skipped");
            continue;
        }
        ps_cpu_force(paths[i]);
        ps_dsp_init();
        ps_split *s = ps_split_create(&cfg);
        ps_split_process(s, in, (size_t)frames, cenX, sidX);
        ps_split_destroy(s);

        double worst = 0.0;
        for (int t = 0; t < frames; t++) {
            double a = fabs((double)cenX[t] - (double)cenS[t]);
            double b = fabs((double)sidX[t * 2] - (double)sidS[t * 2]);
            double d = fabs((double)sidX[t * 2 + 1] - (double)sidS[t * 2 + 1]);
            if (a > worst) worst = a;
            if (b > worst) worst = b;
            if (d > worst) worst = d;
        }
        snprintf(label, sizeof(label), "%s agrees with scalar", ps_isa_name(paths[i]));
        snprintf(detail, sizeof(detail), "max difference %.2e", worst);
        if (paths[i] == PS_ISA_SSE2) {
            /* This one has to match exactly, not merely closely. SSE2 is the
             * default precisely so the same input gives the same file on any
             * machine, and a tolerance here would let that quietly lapse - as
             * it does if the compiler is allowed to fuse multiply and add in
             * the scalar code, which some do by default. */
            ok(label, worst == 0.0, detail);
        } else {
            /* Not bit identical by requirement: fused multiply-add rounds once
             * where scalar rounds twice. It does have to be inaudibly close. */
            ok(label, worst < 1e-5, detail);
        }
    }

    ps_cpu_force((ps_isa)-1);
    ps_dsp_init();
    free(in); free(cenS); free(sidS); free(cenX); free(sidX);
}

static void test_simd_speed(void) {
    const int rate = 44100, frames = rate * 10;
    const ps_isa native = ps_cpu_isa();

    float *in  = (float *)malloc((size_t)frames * 2 * sizeof(float));
    float *cen = (float *)malloc((size_t)frames * sizeof(float));
    float *sid = (float *)malloc((size_t)frames * 2 * sizeof(float));
    if (!in || !cen || !sid) return;
    srand(7);
    for (int t = 0; t < frames * 2; t++) in[t] = (float)((rand() / (double)RAND_MAX) - 0.5);

    ps_split_cfg cfg; ps_split_default_cfg(&cfg);
    double base = 0.0;

    static const ps_isa paths[] = { PS_ISA_SCALAR, PS_ISA_SSE2, PS_ISA_AVX2 };
    for (size_t i = 0; i < sizeof(paths) / sizeof(paths[0]); i++) {
        if (paths[i] > native) continue;
        ps_cpu_force(paths[i]);
        ps_dsp_init();
        ps_split *s = ps_split_create(&cfg);
        clock_t t0 = clock();
        ps_split_process(s, in, (size_t)frames, cen, sid);
        clock_t t1 = clock();
        ps_split_destroy(s);
        double ms = 1000.0 * (double)(t1 - t0) / (double)CLOCKS_PER_SEC;
        if (i == 0) base = ms;
        printf("        %-10s %6.0f ms for 10 s of audio   %5.1fx real time   %4.2fx scalar\n",
               ps_isa_name(paths[i]), ms, ms > 0 ? 10000.0 / ms : 0.0,
               ms > 0 ? base / ms : 0.0);
    }
    ps_cpu_force((ps_isa)-1);
    ps_dsp_init();
    free(in); free(cen); free(sid);
}

/* Magnitude of a biquad at one frequency, for checking filter shapes. */
static double biquad_mag(const ps_biquad *q, double f, double rate) {
    double w = 2.0 * M_PI * f / rate;
    double cr = q->b0 + q->b1 * cos(w) + q->b2 * cos(2 * w);
    double ci =       -q->b1 * sin(w) - q->b2 * sin(2 * w);
    double dr = 1.0   + q->a1 * cos(w) + q->a2 * cos(2 * w);
    double di =       -q->a1 * sin(w) - q->a2 * sin(2 * w);
    return sqrt((cr * cr + ci * ci) / (dr * dr + di * di));
}

/* The delay lengths that matter, taken from the chain this replaces.
 * 0.015 and 0.005 are the interesting ones: both look like an exact half
 * sample once multiplied in double, but the true products fall on opposite
 * sides of the half way point, so they must round in opposite directions. Get
 * this wrong and a surround channel sits one sample out from the others. */
static void test_delay_rounding(void) {
    static const struct { double sec; size_t want; } cases[] = {
        { 0.015,   661 }, { 0.005,   221 }, { 0.0125,  551 }, { 0.02,    882 },
        { 0.017,   750 }, { 0.0102,  450 }, { 0.009,   397 }, { 0.0255, 1125 },
    };
    int bad = 0;
    size_t first_bad_got = 0, first_bad_want = 0; double first_bad_sec = 0;
    for (size_t i = 0; i < sizeof(cases) / sizeof(cases[0]); i++) {
        size_t got = ps_fx_delay_samples(cases[i].sec, 44100);
        if (got != cases[i].want) {
            if (!bad) { first_bad_sec = cases[i].sec; first_bad_got = got;
                        first_bad_want = cases[i].want; }
            bad++;
        }
    }
    char d[96];
    if (bad) snprintf(d, sizeof(d), "%g s gave %u, wanted %u",
                      first_bad_sec, (unsigned)first_bad_got, (unsigned)first_bad_want);
    else     snprintf(d, sizeof(d), "all 8 lengths agree");
    ok("delay lengths round the way they must", bad == 0, d);
}

static void test_fx_filters(void) {
    const double sr = 44100.0;
    ps_biquad hp = ps_biquad_highpass(7000.0, sr);

    /* A highpass must stop DC completely and pass the top end untouched. */
    double dc = biquad_mag(&hp, 0.0, sr), top = biquad_mag(&hp, 22049.0, sr);
    { char d[64]; snprintf(d, sizeof(d), "DC %.2e, Nyquist %.6f", dc, top);
      ok("highpass blocks DC and passes the top", dc < 1e-12 && fabs(top - 1.0) < 1e-6, d); }

    /* At the corner a Butterworth sits at -3 dB. */
    double corner = 20.0 * log10(biquad_mag(&hp, 7000.0, sr));
    { char d[64]; snprintf(d, sizeof(d), "%.4f dB at the corner", corner);
      ok("highpass is -3 dB at its corner", fabs(corner + 3.0103) < 0.01, d); }

    /* The same shape upside down: passes DC, stops the top, -3 dB at the corner. */
    ps_biquad lp = ps_biquad_lowpass(90.0, sr);
    double ldc = biquad_mag(&lp, 0.0, sr), ltop = biquad_mag(&lp, 22049.0, sr);
    double lcorner = 20.0 * log10(biquad_mag(&lp, 90.0, sr));
    { char d[80]; snprintf(d, sizeof(d), "DC %.6f, Nyquist %.2e, corner %.4f dB",
                           ldc, ltop, lcorner);
      ok("lowpass passes DC and blocks the top",
         fabs(ldc - 1.0) < 1e-9 && ltop < 1e-9 && fabs(lcorner + 3.0103) < 0.01, d); }

    /* The shelf: flat below, the asked-for gain above, and only one pole. The
     * order matters - a second order shelf was tried first and never fitted,
     * so if a2 or b2 ever stops being zero the shape has drifted. */
    ps_biquad sh = ps_biquad_highshelf(-3.0, 7000.0, sr);
    double sdc = 20.0 * log10(biquad_mag(&sh, 1.0, sr));
    double stop = 20.0 * log10(biquad_mag(&sh, 22049.0, sr));
    /* The shelf reaches exactly the gain it was given, so this compares against
     * 3.0 dB. The half power figure of 3.0103 belongs to a filter corner, which
     * is a different thing entirely. */
    { char d[80]; snprintf(d, sizeof(d), "%.4f dB low, %.4f dB high", sdc, stop);
      ok("high shelf is flat below and cut above",
         fabs(sdc) < 0.001 && fabs(stop + 3.0) < 0.001, d); }
    ok("the shelf really is first order", sh.b2 == 0.0 && sh.a2 == 0.0,
       "no second order terms");

    /* Gains of opposite sign have to undo one another, which they only do if
     * the two corners really are mirrored about the nominal frequency. */
    ps_biquad up = ps_biquad_highshelf(6.0, 5000.0, sr);
    ps_biquad dn = ps_biquad_highshelf(-6.0, 5000.0, sr);
    double worstpair = 0.0;
    for (double f = 20.0; f < 22050.0; f *= 1.1) {
        double e = fabs(biquad_mag(&up, f, sr) * biquad_mag(&dn, f, sr) - 1.0);
        if (e > worstpair) worstpair = e;
    }
    { char d[64]; snprintf(d, sizeof(d), "worst %.2e", worstpair);
      ok("boost and cut of the same size cancel", worstpair < 1e-12, d); }

    /* An allpass earns its name by leaving every magnitude alone; only phase
     * moves. If this drifts the surrounds would be recoloured, not widened. */
    ps_biquad ap = ps_biquad_allpass(200.0, 1.0, sr);
    double worst = 0.0;
    for (double f = 10.0; f < 22050.0; f *= 1.05) {
        double e = fabs(biquad_mag(&ap, f, sr) - 1.0);
        if (e > worst) worst = e;
    }
    { char d[64]; snprintf(d, sizeof(d), "worst deviation %.2e", worst);
      ok("allpass leaves every magnitude alone", worst < 1e-9, d); }

    /* Rounding goes away from zero and stops at the 32 bit rails. */
    int r = ps_round32(0.5) == 1.0 && ps_round32(-0.5) == -1.0
         && ps_round32(1.4) == 1.0 && ps_round32(-1.4) == -1.0
         && ps_round32(4e9) == 2147483647.0;
    ok("sample rounding goes away from zero", r, r ? "and clips at the rails" : "");

    /* A value already on the 24 bit grid must survive untouched. */
    int keep = ps_regrid24(256.0) == 256.0 && ps_regrid24(-256.0) == -256.0
            && ps_regrid24(128.0) == 256.0;
    ok("24-bit regrid keeps what already fits", keep, "");
}

/* The mixing map. The important one is that "1+2" means both at full weight.
 * A spec language that silently divides by the number of terms is how a
 * surround channel can sit 6 dB low for years without anyone spotting it, so
 * that behaviour is checked for explicitly and must not appear. */
static void test_mix_map(void) {
    double m[PS_MIX_MAX_OUT * PS_MIX_MAX_CH];
    const char *msg = NULL;
    char d[96];

    int n = ps_mix_matrix("1;2;3", 3, m, PS_MIX_MAX_OUT, &msg);
    int straight = (n == 3) && m[0] == 1.0 && m[1] == 0.0 && m[2] == 0.0
                            && m[3 + 1] == 1.0 && m[6 + 2] == 1.0;
    ok("a plain map passes channels through", straight, "");

    n = ps_mix_matrix("1*0.5+2*0.5", 2, m, PS_MIX_MAX_OUT, &msg);
    snprintf(d, sizeof(d), "got %g and %g", n == 1 ? m[0] : -1.0, n == 1 ? m[1] : -1.0);
    ok("explicit weights are used as written", n == 1 && m[0] == 0.5 && m[1] == 0.5, d);

    n = ps_mix_matrix("1+2", 2, m, PS_MIX_MAX_OUT, &msg);
    snprintf(d, sizeof(d), "got %g and %g, wanted 1 and 1",
             n == 1 ? m[0] : -1.0, n == 1 ? m[1] : -1.0);
    ok("nothing is scaled implicitly", n == 1 && m[0] == 1.0 && m[1] == 1.0, d);

    n = ps_mix_matrix("1+1", 1, m, PS_MIX_MAX_OUT, &msg);
    snprintf(d, sizeof(d), "got %g, wanted 2", n == 1 ? m[0] : -1.0);
    ok("naming a channel twice adds it twice", n == 1 && m[0] == 2.0, d);

    n = ps_mix_matrix("1*-0.5", 1, m, PS_MIX_MAX_OUT, &msg);
    ok("a negative weight is allowed", n == 1 && m[0] == -0.5, "for cancelling");

    struct { const char *map; int nch; } bad[] = {
        { "", 2 }, { "9", 2 }, { "0", 2 }, { "1*", 2 },
        { "1;", 2 }, { "1+", 2 }, { "1x2", 2 }, { "abc", 2 },
    };
    int caught = 0;
    for (size_t i = 0; i < sizeof(bad) / sizeof(bad[0]); i++) {
        msg = NULL;
        if (ps_mix_matrix(bad[i].map, bad[i].nch, m, PS_MIX_MAX_OUT, &msg) < 0 && msg)
            caught++;
    }
    snprintf(d, sizeof(d), "%d of 8 rejected with a reason", caught);
    ok("a broken map is refused, not guessed at", caught == 8, d);
}

/* An unbalanced pair: the same tone in both channels, but one of them 40 dB
 * down. Whatever the two share cannot be louder than the quieter of them, so
 * the centre has to come out at about the level of the quiet channel.
 *
 * The estimate before the bound was added is a geometric mean of the two, which
 * for a 40 dB imbalance lands 20 dB too high - a tenfold overshoot. That cost
 * several dB on dense material, where hard panned parts are everywhere, so this
 * checks the bound is still doing its job. */
static void test_centre_bound(void) {
    const int rate = 44100, frames = 44100;
    const double quiet = 0.01;            /* -40 dB against the loud channel */

    ps_split_cfg cfg;
    ps_split_default_cfg(&cfg);
    ps_split *s = ps_split_create(&cfg);
    if (!s) { ok("split create", 0, "returned NULL"); return; }

    float *in  = (float *)malloc((size_t)frames * 2 * sizeof(float));
    float *cen = (float *)malloc((size_t)frames * sizeof(float));
    float *sid = (float *)malloc((size_t)frames * 2 * sizeof(float));
    if (!in || !cen || !sid) { ok("bound alloc", 0, NULL); return; }

    for (int t = 0; t < frames; t++) {
        double v = 0.5 * sin(2.0 * M_PI * 1000.0 * t / rate);
        in[t * 2]     = (float)v;
        in[t * 2 + 1] = (float)(v * quiet);
    }
    ps_split_process(s, in, (size_t)frames, cen, sid);

    /* Measured over the middle, away from the start-up and run-out ramps. */
    double e = 0.0;
    int lo = frames / 4, hi = frames - frames / 4;
    for (int t = lo; t < hi; t++) e += (double)cen[t] * cen[t];
    double rms = sqrt(e / (hi - lo));
    double want = 0.5 * quiet / sqrt(2.0);         /* the quiet channel's level */
    double over = 20.0 * log10(rms / want);

    char d[96];
    snprintf(d, sizeof(d), "%.1f dB against the quiet channel", over);
    /* Allow a little slack for windowing, but nothing like the 20 dB the
     * unbounded estimate would give. */
    ok("centre stays under the quieter channel", over < 3.0, d);

    free(in); free(cen); free(sid);
    ps_split_destroy(s);
}

int ps_selftest(void) {
    g_pass = g_fail = 0;
    printf("phasesplit %s self test\n\n", PS_VERSION);
    printf("-- transform\n");
    test_fft(64); test_fft(512); test_fft(4096); test_fft(8192);

    printf("\n-- framing and separation\n");
    test_split_reconstruction(1024, 2);
    test_split_reconstruction(1024, 4);
    test_split_reconstruction(2048, 2);
    test_split_reconstruction(4096, 4);
    test_split_reconstruction(4096, 8);
    test_split_reconstruction(8192, 4);
    test_split_reconstruction(16384, 4);
    test_centre_bound();

    printf("\n-- surround stage effects\n");
    test_delay_rounding();
    test_fx_filters();
    test_mix_map();

    printf("\n-- instruction set paths (this CPU: %s)\n", ps_isa_name(ps_cpu_isa()));
    test_simd_agreement();
    printf("\n-- throughput\n");
    test_simd_speed();

    printf("\n  %d passed, %d failed\n", g_pass, g_fail);
    return g_fail == 0 ? 0 : 1;
}
