/* Real-input FFT, written for this project so the whole tool stays free of
 * third party code.
 *
 * An N-point real transform is done as an N/2-point complex transform plus a
 * split step, which is the standard trick and roughly halves the work.
 *
 * Complex data is held as separate real and imaginary arrays rather than
 * interleaved pairs. That costs nothing here and makes the per-bin arithmetic
 * straightforward to vectorise later, since eight real parts and eight
 * imaginary parts load as two clean vectors.
 */
#ifndef PS_FFT_H
#define PS_FFT_H

typedef struct ps_fft ps_fft;

/* nreal must be a power of two and at least 8. */
ps_fft *ps_fft_create(int nreal);
void    ps_fft_destroy(ps_fft *p);

int     ps_fft_size(const ps_fft *p);      /* nreal */
int     ps_fft_bins(const ps_fft *p);      /* nreal/2 + 1 */

/* time -> frequency. in holds nreal samples, out_re/out_im hold nreal/2+1. */
void    ps_fft_forward(ps_fft *p, const float *in, float *out_re, float *out_im);

/* frequency -> time, scaled so that forward followed by inverse is identity. */
void    ps_fft_inverse(ps_fft *p, const float *in_re, const float *in_im, float *out);

#endif
