/* The handful of effects the surround chain actually needs.
 *
 * These exist so the whole chain can run in one pass. Handing a file to an
 * external tool a dozen times costs far more in reading and writing it than
 * the arithmetic ever does, and each hand-off also rounds the samples back to
 * 24 bits, so the rounding error piles up a dozen times over.
 *
 * The filter shapes are the standard published audio EQ cookbook biquads. The
 * exact parameter readings were established by measurement - feeding an
 * impulse through the earlier processing and solving for the coefficients -
 * rather than by copying anyone's source, so nothing here is a derived work.
 * What that measurement established:
 *
 *   highpass f        two pole, Q = 1/sqrt(2)  (Butterworth)
 *   allpass  f w      Q = f/w, so the width is a bandwidth in Hz, not a Q
 *   treble   g f      FIRST order shelf, corners at tan(pi f/fs) times and
 *                     divided by sqrt of the gain - not a two pole shelf
 *   delay    t        round(t * rate) samples of silence in front, in double
 *                     arithmetic - 0.015 s at 44100 gives 661 rather than 662
 *                     because 0.015 has no exact binary form
 *
 * PS_FX_STAGE_ROUND puts every stage's output back on the 32 bit sample
 * grid before the next one sees it. That is what happens unavoidably when
 * each stage is a separate program handing on a file, so switching it on
 * makes the two directly comparable. It is not wanted in itself: without it
 * everything stays in double and is rounded once, at the end.
 */
#ifndef PS_FX_H
#define PS_FX_H

#include <stddef.h>

#define PS_FX_STAGE_ROUND 1

typedef struct { double b0, b1, b2, a1, a2; } ps_biquad;
typedef struct { double x1, x2, y1, y2; } ps_biquad_state;

/* Full scale as a 32 bit integer, which is the grid stage boundaries round
 * on. Samples are held as -1..1 elsewhere, so this is the conversion factor. */
#define PS_FULL32 2147483648.0

ps_biquad ps_biquad_highpass(double f0, double rate);
ps_biquad ps_biquad_lowpass(double f0, double rate);
ps_biquad ps_biquad_allpass(double f0, double width_hz, double rate);
ps_biquad ps_biquad_highshelf(double gain_db, double f0, double rate);

/* Round and clip one sample to the 32 bit integer grid, away from zero, which
 * is what a stage boundary costs when the stages are separate programs. */
double ps_round32(double v);

/* Run a biquad over n samples in place. State carries across calls so a signal
 * can be filtered in blocks. With flags & PS_FX_STAGE_ROUND each output sample is
 * put on the integer grid before it feeds back into the filter. */
void ps_biquad_run(const ps_biquad *c, ps_biquad_state *st,
                   double *buf, size_t n, int flags);

/* Multiply by a gain in dB. */
void ps_fx_gain(double *buf, size_t n, double db, int flags);

/* How many samples of silence "delay t" puts in front. */
size_t ps_fx_delay_samples(double seconds, unsigned rate);

/* A 24 bit file between two stages costs a rounding on the way out and gives
 * nothing back on the way in. Reproducing it is needed to compare like with
 * like; in one pass it simply never happens. */
double ps_regrid24(double v);

#endif
