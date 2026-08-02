/* Spectral centre/side separation.
 *
 * The idea, in one paragraph. In each frequency bin, write the two channels as
 * L = C + A and R = C + B, where C is whatever is common to both and A and B
 * are what belongs to only one of them. Take the sum and difference,
 * S = L + R = 2C + A + B and D = L - R = A - B. If A, B and C are not
 * correlated with each other, their powers add, so |S|^2 = 4|C|^2 + |A|^2 +
 * |B|^2 and |D|^2 = |A|^2 + |B|^2. Subtracting gives the common power
 * directly: |C|^2 = (|S|^2 - |D|^2) / 4. The phase of C follows S, since the
 * uncorrelated parts contribute no consistent phase. The sides are then simply
 * what is left over, L - C and R - C, which means centre plus sides always
 * reconstructs the input exactly.
 */
#ifndef PS_SPLIT_H
#define PS_SPLIT_H

#include <stddef.h>

typedef struct ps_split ps_split;

typedef struct {
    int   window;      /* transform length in samples, power of two */
    int   overlap;     /* number of overlapping frames, e.g. 4 = 75% */
    float strength;    /* 0..2, how much of the estimated centre to pull out */
    int   dual_centre; /* non-zero: write the centre as two channels          */
    int   threads;     /* 0 = one per hardware thread                          */
} ps_split_cfg;

void ps_split_default_cfg(ps_split_cfg *cfg);

ps_split *ps_split_create(const ps_split_cfg *cfg);
void      ps_split_destroy(ps_split *s);

/* Offline processing of a whole stereo buffer.
 * in       interleaved stereo, frames * 2 samples
 * centre   frames * (dual_centre ? 2 : 1), may be NULL if not wanted
 * sides    frames * 2, may be NULL if not wanted
 * Returns 0 on success. */
int ps_split_process(ps_split *s, const float *in, size_t frames,
                     float *centre, float *sides);

/* Latency in samples: the output is delayed by one window. */
int ps_split_latency(const ps_split *s);

/* How many workers were actually created. */
int ps_split_threads(const ps_split *s);

#endif
