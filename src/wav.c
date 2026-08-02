#include "wav.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ---- little endian readers, so the code does not depend on host layout ---- */

static uint16_t rd16(const unsigned char *p) { return (uint16_t)(p[0] | (p[1] << 8)); }
static uint32_t rd32(const unsigned char *p) {
    return (uint32_t)p[0] | ((uint32_t)p[1] << 8) | ((uint32_t)p[2] << 16) | ((uint32_t)p[3] << 24);
}
static void wr16(unsigned char *p, uint16_t v) { p[0] = (unsigned char)(v & 0xff); p[1] = (unsigned char)(v >> 8); }
static void wr32(unsigned char *p, uint32_t v) {
    p[0] = (unsigned char)(v & 0xff);        p[1] = (unsigned char)((v >> 8) & 0xff);
    p[2] = (unsigned char)((v >> 16) & 0xff); p[3] = (unsigned char)((v >> 24) & 0xff);
}

const char *ps_fmt_name(ps_sample_fmt f) {
    switch (f) {
        case PS_FMT_PCM16: return "16 bit PCM";
        case PS_FMT_PCM24: return "24 bit PCM";
        case PS_FMT_PCM32: return "32 bit PCM";
        case PS_FMT_F32:   return "32 bit float";
        case PS_FMT_F64:   return "64 bit float";
    }
    return "?";
}

void ps_audio_free(ps_audio *a) {
    if (a && a->data) { free(a->data); a->data = NULL; a->frames = 0; }
}

/* 24 bit signed little endian -> int32, sign extended */
static int32_t rd24s(const unsigned char *p) {
    int32_t v = (int32_t)((uint32_t)p[0] | ((uint32_t)p[1] << 8) | ((uint32_t)p[2] << 16));
    if (v & 0x800000) v |= ~0xffffff;   /* sign extend */
    return v;
}

int ps_wav_read(const char *path, ps_audio *out, const char **msg) {
    const char *dummy = NULL;
    if (!msg) msg = &dummy;
    *msg = NULL;
    memset(out, 0, sizeof(*out));

    FILE *f = fopen(path, "rb");
    if (!f) { *msg = "cannot open file"; return 1; }

    unsigned char hdr[12];
    if (fread(hdr, 1, 12, f) != 12 || memcmp(hdr, "RIFF", 4) || memcmp(hdr + 8, "WAVE", 4)) {
        fclose(f); *msg = "not a RIFF/WAVE file"; return 1;
    }

    unsigned      channels = 0, bits = 0;
    unsigned      rate = 0;
    int           isFloat = 0;
    long          dataPos = -1;
    uint32_t      dataLen = 0;
    int           haveFmt = 0;

    /* Walk the chunks. Anything that is not fmt or data is skipped, which is
     * what lets files with LIST/INFO or broadcast extensions through. */
    for (;;) {
        unsigned char ch[8];
        size_t got = fread(ch, 1, 8, f);
        if (got != 8) break;
        uint32_t id  = rd32(ch);
        uint32_t len = rd32(ch + 4);
        long     pos = ftell(f);

        if (!memcmp(ch, "fmt ", 4)) {
            unsigned char fmt[40];
            uint32_t n = len < sizeof(fmt) ? len : (uint32_t)sizeof(fmt);
            if (fread(fmt, 1, n, f) != n) { fclose(f); *msg = "truncated fmt chunk"; return 1; }
            uint16_t tag = rd16(fmt);
            channels = rd16(fmt + 2);
            rate     = rd32(fmt + 4);
            bits     = rd16(fmt + 14);
            if (tag == 0xFFFE && n >= 40) {          /* WAVE_FORMAT_EXTENSIBLE */
                uint16_t sub = rd16(fmt + 24);       /* first 2 bytes of the GUID */
                isFloat = (sub == 3);
            } else {
                isFloat = (tag == 3);
                if (tag != 1 && tag != 3 && tag != 0xFFFE) {
                    fclose(f); *msg = "unsupported WAV encoding (not PCM or float)"; return 1;
                }
            }
            haveFmt = 1;
        } else if (!memcmp(ch, "data", 4)) {
            dataPos = pos;
            dataLen = len;
        }

        (void)id;
        if (fseek(f, pos + (long)len + ((long)len & 1), SEEK_SET) != 0) break;  /* chunks are word aligned */
        if (dataPos >= 0 && haveFmt) break;
    }

    if (!haveFmt || dataPos < 0) { fclose(f); *msg = "missing fmt or data chunk"; return 1; }
    if (channels == 0)           { fclose(f); *msg = "zero channels"; return 1; }

    ps_sample_fmt sf;
    unsigned bytes = bits / 8;
    if      (isFloat && bits == 32) sf = PS_FMT_F32;
    else if (isFloat && bits == 64) sf = PS_FMT_F64;
    else if (!isFloat && bits == 16) sf = PS_FMT_PCM16;
    else if (!isFloat && bits == 24) sf = PS_FMT_PCM24;
    else if (!isFloat && bits == 32) sf = PS_FMT_PCM32;
    else { fclose(f); *msg = "unsupported bit depth"; return 1; }

    size_t frames = dataLen / (bytes * channels);
    if (frames == 0) { fclose(f); *msg = "no sample data"; return 1; }

    float *pcm = (float *)malloc(frames * channels * sizeof(float));
    if (!pcm) { fclose(f); *msg = "out of memory"; return 1; }

    unsigned char *raw = (unsigned char *)malloc((size_t)bytes * channels * 4096);
    if (!raw) { free(pcm); fclose(f); *msg = "out of memory"; return 1; }

    if (fseek(f, dataPos, SEEK_SET) != 0) { free(raw); free(pcm); fclose(f); *msg = "seek failed"; return 1; }

    const float inv16 = 1.0f / 32768.0f;
    const float inv24 = 1.0f / 8388608.0f;
    const float inv32 = 1.0f / 2147483648.0f;

    size_t done = 0;
    while (done < frames) {
        size_t want = frames - done; if (want > 4096) want = 4096;
        size_t rdF  = fread(raw, (size_t)bytes * channels, want, f);
        if (rdF == 0) break;
        size_t n = rdF * channels;
        float *dst = pcm + done * channels;
        switch (sf) {
            case PS_FMT_PCM16: for (size_t i = 0; i < n; i++) dst[i] = (float)(int16_t)rd16(raw + i * 2) * inv16; break;
            case PS_FMT_PCM24: for (size_t i = 0; i < n; i++) dst[i] = (float)rd24s(raw + i * 3) * inv24; break;
            case PS_FMT_PCM32: for (size_t i = 0; i < n; i++) dst[i] = (float)(int32_t)rd32(raw + i * 4) * inv32; break;
            case PS_FMT_F32:   memcpy(dst, raw, n * 4); break;
            case PS_FMT_F64:   for (size_t i = 0; i < n; i++) { double d; memcpy(&d, raw + i * 8, 8); dst[i] = (float)d; } break;
        }
        done += rdF;
    }
    free(raw);
    fclose(f);

    if (done == 0) { free(pcm); *msg = "could not read sample data"; return 1; }

    out->data     = pcm;
    out->frames   = done;
    out->channels = channels;
    out->rate     = rate;
    out->src_fmt  = sf;
    return 0;
}

static int32_t clamp32(double v, double lo, double hi) {
    /* Round to nearest, not toward zero. Truncating here would pull every
     * sample a little way toward silence instead of scattering the error
     * either side of it, which is a bias rather than noise and costs about
     * half a bit of accuracy on every write. */
    v = floor(v + 0.5);
    if (v < lo) v = lo;
    if (v > hi) v = hi;
    return (int32_t)v;
}

int ps_wav_write(const char *path, const ps_audio *in, ps_sample_fmt fmt, const char **msg) {
    const char *dummy = NULL;
    if (!msg) msg = &dummy;
    *msg = NULL;

    unsigned bytes, bits;
    int isFloat = 0;
    switch (fmt) {
        case PS_FMT_PCM16: bytes = 2; bits = 16; break;
        case PS_FMT_PCM24: bytes = 3; bits = 24; break;
        case PS_FMT_PCM32: bytes = 4; bits = 32; break;
        case PS_FMT_F32:   bytes = 4; bits = 32; isFloat = 1; break;
        case PS_FMT_F64:   bytes = 8; bits = 64; isFloat = 1; break;
        default: *msg = "bad output format"; return 1;
    }

    FILE *f = fopen(path, "wb");
    if (!f) { *msg = "cannot create file"; return 1; }

    const unsigned ch    = in->channels;
    const uint64_t nData = (uint64_t)in->frames * ch * bytes;
    /* Anything at or over 4 GB will not fit a RIFF size field. Say so rather
     * than writing a header that lies about the length. */
    if (nData + 44 > 0xFFFFFFFFull) { fclose(f); *msg = "output would exceed the 4 GB WAV limit"; return 1; }

    unsigned char h[44];
    memcpy(h, "RIFF", 4);          wr32(h + 4, (uint32_t)(36 + nData));
    memcpy(h + 8, "WAVEfmt ", 8);  wr32(h + 16, 16);
    wr16(h + 20, (uint16_t)(isFloat ? 3 : 1));
    wr16(h + 22, (uint16_t)ch);
    wr32(h + 24, in->rate);
    wr32(h + 28, in->rate * ch * bytes);          /* byte rate   */
    wr16(h + 32, (uint16_t)(ch * bytes));         /* block align */
    wr16(h + 34, (uint16_t)bits);
    memcpy(h + 36, "data", 4);     wr32(h + 40, (uint32_t)nData);
    if (fwrite(h, 1, 44, f) != 44) { fclose(f); *msg = "write failed"; return 1; }

    unsigned char *raw = (unsigned char *)malloc((size_t)bytes * ch * 4096);
    if (!raw) { fclose(f); *msg = "out of memory"; return 1; }

    size_t done = 0;
    while (done < in->frames) {
        size_t want = in->frames - done; if (want > 4096) want = 4096;
        size_t n = want * ch;
        const float *src = in->data + done * ch;
        for (size_t i = 0; i < n; i++) {
            double v = (double)src[i];
            switch (fmt) {
                case PS_FMT_PCM16: { int32_t s = clamp32(v * 32768.0, -32768.0, 32767.0);         wr16(raw + i * 2, (uint16_t)s); break; }
                case PS_FMT_PCM24: { int32_t s = clamp32(v * 8388608.0, -8388608.0, 8388607.0);
                                     raw[i*3] = (unsigned char)(s & 0xff);
                                     raw[i*3+1] = (unsigned char)((s >> 8) & 0xff);
                                     raw[i*3+2] = (unsigned char)((s >> 16) & 0xff); break; }
                case PS_FMT_PCM32: { int32_t s = clamp32(v * 2147483648.0, -2147483648.0, 2147483647.0); wr32(raw + i * 4, (uint32_t)s); break; }
                case PS_FMT_F32:   { float  s = (float)v;  memcpy(raw + i * 4, &s, 4); break; }
                case PS_FMT_F64:   {                       memcpy(raw + i * 8, &v, 8); break; }
            }
        }
        if (fwrite(raw, (size_t)bytes * ch, want, f) != want) { free(raw); fclose(f); *msg = "write failed"; return 1; }
        done += want;
    }

    free(raw);
    fclose(f);
    return 0;
}
