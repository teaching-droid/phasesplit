/* Minimal WAV reader/writer for the formats this tool has to deal with:
 * PCM 16/24/32 bit integer and 32/64 bit float, plain RIFF or WAVE_FORMAT_EXTENSIBLE.
 *
 * Samples are always handed to the caller as interleaved 32 bit float in the
 * range -1..1, whatever was on disk, and written back in whatever format the
 * caller asks for.
 */
#ifndef PS_WAV_H
#define PS_WAV_H

#include <stddef.h>
#include <stdint.h>

typedef enum {
    PS_FMT_PCM16 = 0,
    PS_FMT_PCM24,
    PS_FMT_PCM32,
    PS_FMT_F32,
    PS_FMT_F64
} ps_sample_fmt;

typedef struct {
    float        *data;        /* interleaved, frames * channels           */
    size_t        frames;      /* sample frames (per channel)              */
    unsigned      channels;
    unsigned      rate;        /* sample rate in Hz                        */
    ps_sample_fmt src_fmt;     /* what it was on disk, so it can round-trip */
} ps_audio;

/* Returns 0 on success, non-zero on failure; msg gets a short reason. */
int  ps_wav_read(const char *path, ps_audio *out, const char **msg);
int  ps_wav_write(const char *path, const ps_audio *in, ps_sample_fmt fmt, const char **msg);
void ps_audio_free(ps_audio *a);

const char *ps_fmt_name(ps_sample_fmt f);

#endif
