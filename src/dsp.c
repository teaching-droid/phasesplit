#include "dsp.h"
#include "cpu.h"

#include <math.h>

#if defined(_MSC_VER)
  #include <intrin.h>
  #define PS_TARGET_AVX2
#elif defined(__x86_64__) || defined(__i386__)
  #include <immintrin.h>
  #define PS_TARGET_AVX2 __attribute__((target("avx2,fma")))
#endif

#if defined(__x86_64__) || defined(_M_X64) || defined(__i386__) || defined(_M_IX86)
  #define PS_X86 1
#else
  #define PS_X86 0
#endif

/* ---------------------------------------------------------------- scalar -- */

static void centre_scalar(const float *Lre, const float *Lim,
                          const float *Rre, const float *Rim,
                          float *Cre, float *Cim, int bins, float strength) {
    for (int k = 0; k < bins; k++) {
        const float sr = Lre[k] + Rre[k], si = Lim[k] + Rim[k];
        const float dr = Lre[k] - Rre[k], di = Lim[k] - Rim[k];
        const float ss = sr * sr + si * si;
        const float dd = dr * dr + di * di;
        float diff = ss - dd;
        if (diff < 0.0f) diff = 0.0f;
        float g = 0.0f;
        if (ss > 1e-20f) {
            /* |S|^2 - |D|^2 is 4*Re(L conj(R)), so this estimate is the
             * geometric mean of the two channels. Whatever is common to both
             * cannot be louder than the quieter of them, and the geometric
             * mean is above the smaller value whenever the two are unbalanced.
             * Without this bound every strongly panned bin hands too much to
             * the centre - which is exactly the material dense productions are
             * made of, and it costs several dB on them. */
            const float ll = Lre[k] * Lre[k] + Lim[k] * Lim[k];
            const float rr = Rre[k] * Rre[k] + Rim[k] * Rim[k];
            const float lim2 = (ll < rr) ? ll : rr;
            float num = 0.5f * sqrtf(diff);
            if (num * num > lim2) num = sqrtf(lim2);
            g = num / sqrtf(ss) * strength;
            if (g > 1.0f) g = 1.0f;
        }
        Cre[k] = sr * g;
        Cim[k] = si * g;
    }
}

static void subtract_scalar(float *Lre, float *Lim, float *Rre, float *Rim,
                            const float *Cre, const float *Cim, int bins) {
    for (int k = 0; k < bins; k++) {
        Lre[k] -= Cre[k]; Lim[k] -= Cim[k];
        Rre[k] -= Cre[k]; Rim[k] -= Cim[k];
    }
}

static void scaleadd_scalar(float *dst, const float *src, const float *wa, float wb, int n) {
    for (int i = 0; i < n; i++) dst[i] += src[i] * wa[i] * wb;
}

static void scaleadddiff_scalar(float *dst, const float *a, const float *b,
                                const float *wa, int n) {
    for (int i = 0; i < n; i++) dst[i] += (a[i] - b[i]) * wa[i];
}

static void winsplit_scalar(const float *in, const float *win,
                            float *outL, float *outR, int n) {
    for (int i = 0; i < n; i++) {
        const float w = win[i];
        outL[i] = in[i * 2]     * w;
        outR[i] = in[i * 2 + 1] * w;
    }
}

#if PS_X86
/* ------------------------------------------------------------------ SSE2 -- */

static void centre_sse2(const float *Lre, const float *Lim,
                        const float *Rre, const float *Rim,
                        float *Cre, float *Cim, int bins, float strength) {
    const __m128 half = _mm_set1_ps(0.5f);
    const __m128 one  = _mm_set1_ps(1.0f);
    const __m128 zero = _mm_setzero_ps();
    const __m128 eps  = _mm_set1_ps(1e-20f);
    const __m128 str  = _mm_set1_ps(strength);

    int k = 0;
    for (; k + 4 <= bins; k += 4) {
        __m128 lr = _mm_loadu_ps(Lre + k), li = _mm_loadu_ps(Lim + k);
        __m128 rr = _mm_loadu_ps(Rre + k), ri = _mm_loadu_ps(Rim + k);

        __m128 sr = _mm_add_ps(lr, rr), si = _mm_add_ps(li, ri);
        __m128 dr = _mm_sub_ps(lr, rr), di = _mm_sub_ps(li, ri);

        __m128 ss = _mm_add_ps(_mm_mul_ps(sr, sr), _mm_mul_ps(si, si));
        __m128 dd = _mm_add_ps(_mm_mul_ps(dr, dr), _mm_mul_ps(di, di));

        __m128 diff = _mm_max_ps(_mm_sub_ps(ss, dd), zero);
        __m128 live = _mm_cmpgt_ps(ss, eps);          /* mask: is there a sum at all */
        __m128 safe = _mm_or_ps(_mm_and_ps(live, ss), _mm_andnot_ps(live, one));

        /* Hold the centre down to the quieter channel - see the scalar version
         * for why the estimate runs above it on strongly panned material. */
        __m128 ll   = _mm_add_ps(_mm_mul_ps(lr, lr), _mm_mul_ps(li, li));
        __m128 rr2  = _mm_add_ps(_mm_mul_ps(rr, rr), _mm_mul_ps(ri, ri));
        __m128 lim2 = _mm_min_ps(ll, rr2);
        __m128 num  = _mm_mul_ps(half, _mm_sqrt_ps(diff));
        num = _mm_min_ps(num, _mm_sqrt_ps(lim2));

        __m128 g = _mm_mul_ps(_mm_div_ps(num, _mm_sqrt_ps(safe)), str);
        g = _mm_min_ps(g, one);
        g = _mm_and_ps(g, live);

        _mm_storeu_ps(Cre + k, _mm_mul_ps(sr, g));
        _mm_storeu_ps(Cim + k, _mm_mul_ps(si, g));
    }
    if (k < bins)
        centre_scalar(Lre + k, Lim + k, Rre + k, Rim + k, Cre + k, Cim + k, bins - k, strength);
}

static void subtract_sse2(float *Lre, float *Lim, float *Rre, float *Rim,
                          const float *Cre, const float *Cim, int bins) {
    int k = 0;
    for (; k + 4 <= bins; k += 4) {
        __m128 cr = _mm_loadu_ps(Cre + k), ci = _mm_loadu_ps(Cim + k);
        _mm_storeu_ps(Lre + k, _mm_sub_ps(_mm_loadu_ps(Lre + k), cr));
        _mm_storeu_ps(Lim + k, _mm_sub_ps(_mm_loadu_ps(Lim + k), ci));
        _mm_storeu_ps(Rre + k, _mm_sub_ps(_mm_loadu_ps(Rre + k), cr));
        _mm_storeu_ps(Rim + k, _mm_sub_ps(_mm_loadu_ps(Rim + k), ci));
    }
    if (k < bins) subtract_scalar(Lre + k, Lim + k, Rre + k, Rim + k, Cre + k, Cim + k, bins - k);
}

static void scaleadd_sse2(float *dst, const float *src, const float *wa, float wb, int n) {
    const __m128 b = _mm_set1_ps(wb);
    int i = 0;
    for (; i + 4 <= n; i += 4) {
        __m128 v = _mm_mul_ps(_mm_mul_ps(_mm_loadu_ps(src + i), _mm_loadu_ps(wa + i)), b);
        _mm_storeu_ps(dst + i, _mm_add_ps(_mm_loadu_ps(dst + i), v));
    }
    if (i < n) scaleadd_scalar(dst + i, src + i, wa + i, wb, n - i);
}

static void scaleadddiff_sse2(float *dst, const float *a, const float *b,
                              const float *wa, int n) {
    int i = 0;
    for (; i + 4 <= n; i += 4) {
        __m128 v = _mm_mul_ps(_mm_sub_ps(_mm_loadu_ps(a + i), _mm_loadu_ps(b + i)),
                              _mm_loadu_ps(wa + i));
        _mm_storeu_ps(dst + i, _mm_add_ps(_mm_loadu_ps(dst + i), v));
    }
    if (i < n) scaleadddiff_scalar(dst + i, a + i, b + i, wa + i, n - i);
}

static void winsplit_sse2(const float *in, const float *win,
                          float *outL, float *outR, int n) {
    int i = 0;
    for (; i + 4 <= n; i += 4) {
        __m128 a = _mm_loadu_ps(in + i * 2);       /* L0 R0 L1 R1 */
        __m128 b = _mm_loadu_ps(in + i * 2 + 4);   /* L2 R2 L3 R3 */
        __m128 l = _mm_shuffle_ps(a, b, _MM_SHUFFLE(2, 0, 2, 0));
        __m128 r = _mm_shuffle_ps(a, b, _MM_SHUFFLE(3, 1, 3, 1));
        __m128 w = _mm_loadu_ps(win + i);
        _mm_storeu_ps(outL + i, _mm_mul_ps(l, w));
        _mm_storeu_ps(outR + i, _mm_mul_ps(r, w));
    }
    if (i < n) winsplit_scalar(in + i * 2, win + i, outL + i, outR + i, n - i);
}

/* ------------------------------------------------------------------ AVX2 -- */

PS_TARGET_AVX2
static void centre_avx2(const float *Lre, const float *Lim,
                        const float *Rre, const float *Rim,
                        float *Cre, float *Cim, int bins, float strength) {
    const __m256 half = _mm256_set1_ps(0.5f);
    const __m256 one  = _mm256_set1_ps(1.0f);
    const __m256 zero = _mm256_setzero_ps();
    const __m256 eps  = _mm256_set1_ps(1e-20f);
    const __m256 str  = _mm256_set1_ps(strength);

    int k = 0;
    for (; k + 8 <= bins; k += 8) {
        __m256 lr = _mm256_loadu_ps(Lre + k), li = _mm256_loadu_ps(Lim + k);
        __m256 rr = _mm256_loadu_ps(Rre + k), ri = _mm256_loadu_ps(Rim + k);

        __m256 sr = _mm256_add_ps(lr, rr), si = _mm256_add_ps(li, ri);
        __m256 dr = _mm256_sub_ps(lr, rr), di = _mm256_sub_ps(li, ri);

        __m256 ss = _mm256_fmadd_ps(si, si, _mm256_mul_ps(sr, sr));
        __m256 dd = _mm256_fmadd_ps(di, di, _mm256_mul_ps(dr, dr));

        __m256 diff = _mm256_max_ps(_mm256_sub_ps(ss, dd), zero);
        __m256 live = _mm256_cmp_ps(ss, eps, _CMP_GT_OQ);
        __m256 safe = _mm256_blendv_ps(one, ss, live);

        /* Same bound as the other two paths: the centre cannot be louder than
         * the quieter channel. */
        __m256 ll   = _mm256_fmadd_ps(li, li, _mm256_mul_ps(lr, lr));
        __m256 rr2  = _mm256_fmadd_ps(ri, ri, _mm256_mul_ps(rr, rr));
        __m256 lim2 = _mm256_min_ps(ll, rr2);
        __m256 num  = _mm256_mul_ps(half, _mm256_sqrt_ps(diff));
        num = _mm256_min_ps(num, _mm256_sqrt_ps(lim2));

        __m256 g = _mm256_mul_ps(_mm256_div_ps(num, _mm256_sqrt_ps(safe)), str);
        g = _mm256_min_ps(g, one);
        g = _mm256_and_ps(g, live);

        _mm256_storeu_ps(Cre + k, _mm256_mul_ps(sr, g));
        _mm256_storeu_ps(Cim + k, _mm256_mul_ps(si, g));
    }
    if (k < bins)
        centre_sse2(Lre + k, Lim + k, Rre + k, Rim + k, Cre + k, Cim + k, bins - k, strength);
}

PS_TARGET_AVX2
static void subtract_avx2(float *Lre, float *Lim, float *Rre, float *Rim,
                          const float *Cre, const float *Cim, int bins) {
    int k = 0;
    for (; k + 8 <= bins; k += 8) {
        __m256 cr = _mm256_loadu_ps(Cre + k), ci = _mm256_loadu_ps(Cim + k);
        _mm256_storeu_ps(Lre + k, _mm256_sub_ps(_mm256_loadu_ps(Lre + k), cr));
        _mm256_storeu_ps(Lim + k, _mm256_sub_ps(_mm256_loadu_ps(Lim + k), ci));
        _mm256_storeu_ps(Rre + k, _mm256_sub_ps(_mm256_loadu_ps(Rre + k), cr));
        _mm256_storeu_ps(Rim + k, _mm256_sub_ps(_mm256_loadu_ps(Rim + k), ci));
    }
    if (k < bins) subtract_sse2(Lre + k, Lim + k, Rre + k, Rim + k, Cre + k, Cim + k, bins - k);
}

PS_TARGET_AVX2
static void scaleadd_avx2(float *dst, const float *src, const float *wa, float wb, int n) {
    const __m256 b = _mm256_set1_ps(wb);
    int i = 0;
    for (; i + 8 <= n; i += 8) {
        __m256 v = _mm256_mul_ps(_mm256_mul_ps(_mm256_loadu_ps(src + i), _mm256_loadu_ps(wa + i)), b);
        _mm256_storeu_ps(dst + i, _mm256_add_ps(_mm256_loadu_ps(dst + i), v));
    }
    if (i < n) scaleadd_sse2(dst + i, src + i, wa + i, wb, n - i);
}

PS_TARGET_AVX2
static void scaleadddiff_avx2(float *dst, const float *a, const float *b,
                              const float *wa, int n) {
    int i = 0;
    for (; i + 8 <= n; i += 8) {
        __m256 v = _mm256_mul_ps(_mm256_sub_ps(_mm256_loadu_ps(a + i), _mm256_loadu_ps(b + i)),
                                 _mm256_loadu_ps(wa + i));
        _mm256_storeu_ps(dst + i, _mm256_add_ps(_mm256_loadu_ps(dst + i), v));
    }
    if (i < n) scaleadddiff_sse2(dst + i, a + i, b + i, wa + i, n - i);
}

PS_TARGET_AVX2
static void winsplit_avx2(const float *in, const float *win,
                          float *outL, float *outR, int n) {
    int i = 0;
    for (; i + 8 <= n; i += 8) {
        __m256 a = _mm256_loadu_ps(in + i * 2);
        __m256 b = _mm256_loadu_ps(in + i * 2 + 8);
        /* Pull the even and odd slots apart, then put the halves back in order. */
        __m256 e = _mm256_shuffle_ps(a, b, _MM_SHUFFLE(2, 0, 2, 0));
        __m256 o = _mm256_shuffle_ps(a, b, _MM_SHUFFLE(3, 1, 3, 1));
        __m256 l = _mm256_permutevar8x32_ps(e, _mm256_setr_epi32(0, 1, 4, 5, 2, 3, 6, 7));
        __m256 r = _mm256_permutevar8x32_ps(o, _mm256_setr_epi32(0, 1, 4, 5, 2, 3, 6, 7));
        __m256 w = _mm256_loadu_ps(win + i);
        _mm256_storeu_ps(outL + i, _mm256_mul_ps(l, w));
        _mm256_storeu_ps(outR + i, _mm256_mul_ps(r, w));
    }
    if (i < n) winsplit_sse2(in + i * 2, win + i, outL + i, outR + i, n - i);
}
#endif /* PS_X86 */

/* -------------------------------------------------------------- dispatch -- */

static void (*fn_centre)(const float *, const float *, const float *, const float *,
                         float *, float *, int, float) = centre_scalar;
static void (*fn_subtract)(float *, float *, float *, float *,
                           const float *, const float *, int) = subtract_scalar;
static void (*fn_scaleadd)(float *, const float *, const float *, float, int) = scaleadd_scalar;
static void (*fn_scaleadddiff)(float *, const float *, const float *, const float *, int) = scaleadddiff_scalar;
static void (*fn_winsplit)(const float *, const float *, float *, float *, int) = winsplit_scalar;
static const char *g_path = "scalar";
static int g_ready = 0;

void ps_dsp_init_once(void) {
    if (!g_ready) ps_dsp_init();
}

void ps_dsp_init(void) {
    g_ready = 1;
    const ps_isa isa = ps_cpu_isa();
    fn_centre = centre_scalar; fn_subtract = subtract_scalar;
    fn_scaleadd = scaleadd_scalar; fn_winsplit = winsplit_scalar;
    fn_scaleadddiff = scaleadddiff_scalar;
    g_path = "scalar";
#if PS_X86
    if (isa >= PS_ISA_SSE2) {
        fn_centre = centre_sse2; fn_subtract = subtract_sse2;
        fn_scaleadd = scaleadd_sse2; fn_winsplit = winsplit_sse2;
        fn_scaleadddiff = scaleadddiff_sse2;
        g_path = "SSE2";
    }
    if (isa >= PS_ISA_AVX2) {
        fn_centre = centre_avx2; fn_subtract = subtract_avx2;
        fn_scaleadd = scaleadd_avx2; fn_winsplit = winsplit_avx2;
        fn_scaleadddiff = scaleadddiff_avx2;
        g_path = "AVX2+FMA";
    }
#else
    (void)isa;
#endif
}

const char *ps_dsp_path(void) { return g_path; }

void ps_estimate_centre(const float *Lre, const float *Lim,
                        const float *Rre, const float *Rim,
                        float *Cre, float *Cim, int bins, float strength) {
    fn_centre(Lre, Lim, Rre, Rim, Cre, Cim, bins, strength);
}
void ps_subtract_centre(float *Lre, float *Lim, float *Rre, float *Rim,
                        const float *Cre, const float *Cim, int bins) {
    fn_subtract(Lre, Lim, Rre, Rim, Cre, Cim, bins);
}
void ps_scale_add(float *dst, const float *src, const float *wa, float wb, int n) {
    fn_scaleadd(dst, src, wa, wb, n);
}
void ps_scale_add_diff(float *dst, const float *a, const float *b, const float *wa, int n) {
    fn_scaleadddiff(dst, a, b, wa, n);
}
void ps_window_split(const float *in, const float *win, float *outL, float *outR, int n) {
    fn_winsplit(in, win, outL, outR, n);
}
