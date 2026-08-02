#include "fft.h"

#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

struct ps_fft {
    int       nreal;     /* real transform length            */
    int       n;         /* complex transform length, nreal/2 */
    int       log2n;
    unsigned *rev;       /* bit reversal permutation for n    */
    float    *tw_re;     /* butterfly twiddles, n/2 of them   */
    float    *tw_im;
    float    *sp_re;     /* split-step twiddles, n/2          */
    float    *sp_im;
    float    *zr;        /* working complex buffer            */
    float    *zi;
};

static int ilog2(int v) {
    int r = 0;
    while ((1 << r) < v) r++;
    return r;
}

ps_fft *ps_fft_create(int nreal) {
    if (nreal < 8 || (nreal & (nreal - 1)) != 0) return NULL;

    ps_fft *p = (ps_fft *)calloc(1, sizeof(ps_fft));
    if (!p) return NULL;

    p->nreal = nreal;
    p->n     = nreal / 2;
    p->log2n = ilog2(p->n);

    const int n = p->n;
    p->rev   = (unsigned *)malloc((size_t)n * sizeof(unsigned));
    p->tw_re = (float *)malloc((size_t)(n / 2 + 1) * sizeof(float));
    p->tw_im = (float *)malloc((size_t)(n / 2 + 1) * sizeof(float));
    p->sp_re = (float *)malloc((size_t)(n / 2 + 1) * sizeof(float));
    p->sp_im = (float *)malloc((size_t)(n / 2 + 1) * sizeof(float));
    p->zr    = (float *)malloc((size_t)n * sizeof(float));
    p->zi    = (float *)malloc((size_t)n * sizeof(float));
    if (!p->rev || !p->tw_re || !p->tw_im || !p->sp_re || !p->sp_im || !p->zr || !p->zi) {
        ps_fft_destroy(p);
        return NULL;
    }

    /* Bit reversal table for the complex transform. */
    for (int i = 0; i < n; i++) {
        unsigned r = 0;
        for (int b = 0; b < p->log2n; b++)
            if (i & (1 << b)) r |= 1u << (p->log2n - 1 - b);
        p->rev[i] = r;
    }

    /* Twiddles are computed in double and stored as float: the error in the
     * table itself would otherwise dominate everything downstream. */
    for (int k = 0; k <= n / 2; k++) {
        double a = -2.0 * M_PI * (double)k / (double)n;
        p->tw_re[k] = (float)cos(a);
        p->tw_im[k] = (float)sin(a);

        double b = -2.0 * M_PI * (double)k / (double)nreal;
        p->sp_re[k] = (float)cos(b);
        p->sp_im[k] = (float)sin(b);
    }
    return p;
}

void ps_fft_destroy(ps_fft *p) {
    if (!p) return;
    free(p->rev); free(p->tw_re); free(p->tw_im);
    free(p->sp_re); free(p->sp_im); free(p->zr); free(p->zi);
    free(p);
}

int ps_fft_size(const ps_fft *p) { return p ? p->nreal : 0; }
int ps_fft_bins(const ps_fft *p) { return p ? p->n + 1 : 0; }

/* In-place iterative radix-2 complex transform on split arrays.
 * sign is -1 for forward, +1 for inverse (no scaling applied here). */
static void cfft(const ps_fft *p, float *re, float *im, int sign) {
    const int n = p->n;

    for (int i = 0; i < n; i++) {
        unsigned j = p->rev[i];
        if ((unsigned)i < j) {
            float t;
            t = re[i]; re[i] = re[j]; re[j] = t;
            t = im[i]; im[i] = im[j]; im[j] = t;
        }
    }

    for (int len = 2; len <= n; len <<= 1) {
        const int half = len >> 1;
        const int step = n / len;          /* twiddle stride */
        for (int i = 0; i < n; i += len) {
            for (int k = 0; k < half; k++) {
                const int   ti = k * step;
                const float wr = p->tw_re[ti];
                const float wi = (sign < 0) ? p->tw_im[ti] : -p->tw_im[ti];
                const int   a  = i + k;
                const int   b  = a + half;
                const float xr = re[b] * wr - im[b] * wi;
                const float xi = re[b] * wi + im[b] * wr;
                re[b] = re[a] - xr;  im[b] = im[a] - xi;
                re[a] = re[a] + xr;  im[a] = im[a] + xi;
            }
        }
    }
}

void ps_fft_forward(ps_fft *p, const float *in, float *out_re, float *out_im) {
    const int n = p->n;

    /* Pack the real signal as n complex values: z[k] = x[2k] + i x[2k+1] */
    for (int k = 0; k < n; k++) {
        p->zr[k] = in[2 * k];
        p->zi[k] = in[2 * k + 1];
    }

    cfft(p, p->zr, p->zi, -1);

    /* Split step: recover the spectrum of the real signal from Z.
     * Even and odd parts are separated, then recombined with a half-bin
     * twiddle. Bin 0 and the Nyquist bin are both real. */
    const float z0r = p->zr[0], z0i = p->zi[0];
    out_re[0] = z0r + z0i;  out_im[0] = 0.0f;
    out_re[n] = z0r - z0i;  out_im[n] = 0.0f;

    for (int k = 1; k <= n / 2; k++) {
        const int   nk = n - k;
        const float ar = p->zr[k],  ai = p->zi[k];
        const float br = p->zr[nk], bi = -p->zi[nk];   /* conj(Z[n-k]) */

        const float er = 0.5f * (ar + br), ei = 0.5f * (ai + bi);   /* even part  */
        const float or_ = 0.5f * (ai - bi), oi = -0.5f * (ar - br); /* odd part   */

        const float wr = p->sp_re[k], wi = p->sp_im[k];
        const float tr = or_ * wr - oi * wi;
        const float ti = or_ * wi + oi * wr;

        out_re[k] = er + tr;  out_im[k] = ei + ti;
        if (k != nk) {                       /* mirror bin, conjugate symmetric */
            out_re[nk] = er - tr;  out_im[nk] = -(ei - ti);
        }
    }
}

void ps_fft_inverse(ps_fft *p, const float *in_re, const float *in_im, float *out) {
    const int n = p->n;

    /* Undo the split step to get Z back. */
    const float x0 = in_re[0], xn = in_re[n];
    p->zr[0] = 0.5f * (x0 + xn);
    p->zi[0] = 0.5f * (x0 - xn);

    for (int k = 1; k <= n / 2; k++) {
        const int   nk = n - k;
        const float ar = in_re[k],  ai = in_im[k];
        const float br = in_re[nk], bi = -in_im[nk];    /* conj(X[n-k]) */

        const float er = 0.5f * (ar + br), ei = 0.5f * (ai + bi);
        const float dr = 0.5f * (ar - br), di = 0.5f * (ai - bi);

        /* multiply the difference by conj(w) to undo the forward twiddle */
        const float wr = p->sp_re[k], wi = -p->sp_im[k];
        const float fr = dr * wr - di * wi;
        const float fi = dr * wi + di * wr;

        /* z = even + i * odd */
        p->zr[k]  = er - fi;
        p->zi[k]  = ei + fr;
        if (k != nk) {
            p->zr[nk] = er + fi;
            p->zi[nk] = -(ei - fr);
        }
    }

    cfft(p, p->zr, p->zi, +1);

    const float s = 1.0f / (float)n;
    for (int k = 0; k < n; k++) {
        out[2 * k]     = p->zr[k] * s;
        out[2 * k + 1] = p->zi[k] * s;
    }
}
