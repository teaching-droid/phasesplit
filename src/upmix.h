/* The surround stages, each collapsed into a single pass.
 *
 * These replace groups of separate tool invocations that each read the whole
 * file, did one small thing, and wrote it back out. On a four minute file the
 * reading and writing alone costs around 150 ms per pass while the arithmetic
 * costs almost nothing, so removing the passes is worth far more than making
 * any individual effect faster.
 *
 *   freqsplit   was 4 passes: two to pull the high band out of the stem, two
 *               more to blend it onto the surrounds
 *   decorrelate was 1 pass per channel, each a delay followed by two allpasses
 *
 * With PS_FX_STAGE_ROUND the result matches running the stages separately,
 * bit for bit, down to the rounding each hand-off cost. Without it
 * everything stays in double until the single write at the end.
 */
#ifndef PS_UPMIX_H
#define PS_UPMIX_H

#include "wav.h"

/* Pull everything above hz out of the stem, lift it by gain_db, and blend it
 * onto the two surround channels. blend_full keeps the surround at its own
 * level; otherwise both halves are halved, which is the older behaviour.
 * lsPath and rsPath are read and then overwritten. */
int ps_freqsplit(const char *stemPath, const char *lsPath, const char *rsPath,
                 double hz, double gain_db, int blend_full,
                 int flags, ps_sample_fmt outfmt, const char **msg);

/* Delay a surround channel and run it through two allpasses, in place.
 * The delay lengthens the signal, and the allpasses run over the result. */
int ps_decorrelate(const char *path, double delay_sec,
                   double ap1_hz, double ap2_hz, double ap_width,
                   int flags, ps_sample_fmt outfmt, const char **msg);

/* Build output channels as weighted sums of input channels.
 *
 * Nearly every remaining step of the surround chain is this and nothing else:
 * gather some channels, scale them, add them up, write the result. Done as
 * separate programs that is thirty-odd passes over the audio; here it is one.
 *
 * The inputs are treated as one long list of channels laid end to end, so with
 * a stereo file followed by a mono one, channels 1 and 2 are the stereo pair
 * and channel 3 is the mono. Numbering starts at 1.
 *
 * The map is output channels separated by ';', each a sum of terms separated
 * by '+', each term either "N" or "N*gain":
 *
 *     "1*0.5+3*0.5;3*0.5+2*0.5"    two outputs, each mixing two inputs
 *     "1;2;3;4;5;6"                six outputs, straight through
 *
 * Every weight has to be written down. Leaving them implicit and having the
 * program quietly divide by the number of terms is how a surround channel can
 * end up 6 dB down without anyone noticing.
 *
 * Each output channel can then be filtered, gained and cut to length, in that
 * order. Together that covers nearly every step the chain still does outside
 * the two stages above: the low frequency channel is a matrix down to mono
 * followed by two lowpasses and a gain, and trimming is how the padding added
 * for the overlap gets taken back off.
 *
 * gain_db is applied after the filters; pass 0 for none. trim is a length in
 * samples, or 0 to keep everything. */
typedef enum { PS_FILT_LOWPASS, PS_FILT_HIGHPASS, PS_FILT_HIGHSHELF } ps_filt_kind;
typedef struct { ps_filt_kind kind; double hz; double db; } ps_filt;

#define PS_MIX_MAX_FILT 8

int ps_mix(const char **inPaths, int nIn, const char *outPath,
           const char *map, double gain_db,
           const ps_filt *filts, int nFilt, size_t trim,
           int flags, ps_sample_fmt outfmt, const char **msg);

#define PS_MIX_MAX_OUT 32
#define PS_MIX_MAX_CH  64

/* Turn a map into a plain table of weights, m[out * nch + in]. Returns the
 * number of output channels, or -1 with a reason in msg. Separate from ps_mix
 * so the self test can check the weights without going near a file - in
 * particular that nothing is scaled implicitly. */
int ps_mix_matrix(const char *map, int nch, double *m, int maxOut, const char **msg);

#endif
