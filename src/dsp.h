/* The inner loops, in scalar, SSE2 and AVX2 forms.
 *
 * Which one runs is decided once at start-up. Every version has to produce the
 * same numbers as the scalar one - there is a check for exactly that in the
 * self test, because a fast path that quietly differs is worse than no fast
 * path at all.
 */
#ifndef PS_DSP_H
#define PS_DSP_H

void ps_dsp_init(void);
/* Sets the table up only if it has not been chosen already. */
void ps_dsp_init_once(void);
const char *ps_dsp_path(void);      /* which one ended up selected */

/* Per bin: work out the shared component of the two channels.
 * See split.h for where the arithmetic comes from. */
void ps_estimate_centre(const float *Lre, const float *Lim,
                        const float *Rre, const float *Rim,
                        float *Cre, float *Cim,
                        int bins, float strength);

/* Take the centre back out of both channels, leaving the sides. */
void ps_subtract_centre(float *Lre, float *Lim, float *Rre, float *Rim,
                        const float *Cre, const float *Cim, int bins);

/* dst[i] += src[i] * wa[i] * wb  (the overlap-add step) */
void ps_scale_add(float *dst, const float *src, const float *wa, float wb, int n);

/* dst[i] += (a[i] - b[i]) * wa[i]
 * The sides are the frame minus the centre, and taking the difference here
 * rather than in the frequency domain removes two inverse transforms. */
void ps_scale_add_diff(float *dst, const float *a, const float *b, const float *wa, int n);

/* Deinterleave a stereo block and apply the window in one pass. */
void ps_window_split(const float *interleaved, const float *win,
                     float *outL, float *outR, int n);

#endif
