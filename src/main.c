/* phasesplit - spectral centre/side separation for stereo audio.
 *
 * Splits a stereo file into the part that is common to both channels and the
 * part that is not. The two outputs always add back up to the input.
 */
#include "wav.h"
#include "split.h"
#include "selftest.h"
#include "cpu.h"
#include "dsp.h"
#include "fx.h"
#include "upmix.h"
#include "version.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static void usage(void) {
    printf("phasesplit %s - spectral centre/side separation\n\n", PS_VERSION);
    printf("usage: phasesplit <input.wav> [options]\n\n");
    printf("  -c <file>        write the centre channel\n");
    printf("  -s <file>        write the side channels\n");
    printf("  -d, --dual       write the centre as two channels instead of one\n");
    printf("  -w, --window N   transform length, power of two (default 4096)\n");
    printf("  -o, --overlap N  overlapping frames, 2 4 or 8 (default 4)\n");
    printf("      --strength X level of the extracted centre, 0..2 (default 1).\n");
    printf("                   Only the level: whatever it takes out of the\n");
    printf("                   centre is what the sides keep, so the two still\n");
    printf("                   add back to the input. It does not change how\n");
    printf("                   the centre is picked out.\n");
    printf("  -b, --bits N     output depth: 16, 24, 32 or f32 (default: as input)\n");
    printf("  -q, --quiet      no progress output\n");
    printf("      --threads N  worker threads (default: one per core)\n");
    printf("      --isa NAME   scalar, sse2, avx2 or auto (default sse2)\n");
    printf("      --selftest   run the built-in checks and exit\n");
    printf("      --sweep      measure the window/overlap trade-off and exit\n");
    printf("      --cpuinfo    show what this machine supports and exit\n");
    printf("      --copy <f>   write the input straight back out (file layer check)\n");
    printf("  -v, --version    print the version and exit\n");
    printf("  -h, --help       this help\n\n");
    printf("The centre and the sides always sum back to the input, so nothing\n");
    printf("is invented and nothing is lost.\n\n");
    printf("surround stages, each one pass instead of several:\n\n");
    printf("  --freqsplit <stem.wav> --ls <L.wav> --rs <R.wav> --hz N --gain DB\n");
    printf("                   blend the band above N Hz from the stem onto the\n");
    printf("                   surrounds. --blend full keeps the surround level;\n");
    printf("                   the default halves both halves.\n");
    printf("  --decorrelate <file.wav> --delay SEC --ap1 HZ --ap2 HZ [--apwidth HZ]\n");
    printf("                   delay a surround channel and turn its phase.\n");
    printf("  --mix <out.wav> --in <f> [--in <f> ...] --map \"SPEC\" [--gain DB]\n");
    printf("                   build output channels as weighted sums of input\n");
    printf("                   ones. Inputs are laid end to end and numbered\n");
    printf("                   from 1. Outputs are separated by ';', terms\n");
    printf("                   within one by '+', each term N or N*gain:\n");
    printf("                     --map \"1*0.5+3*0.5;2\"\n");
    printf("                   Every weight is written out in full; nothing is\n");
    printf("                   scaled behind your back. Each output can then be\n");
    printf("                   shaped, in this order:\n");
    printf("                     --lowpass HZ / --highpass HZ  (repeatable)\n");
    printf("                     --highshelf DB HZ   tilt the top of the band\n");
    printf("                     --gain DB\n");
    printf("                     --trim N      keep the first N samples\n");
    printf("  --stage-round    round to the sample grid after every stage, not\n");
    printf("                   just at the end. Only useful for comparing against\n");
    printf("                   running the stages as separate programs. Off by\n");
    printf("                   default: rounding once is more accurate.\n");
}

/* Returns 0 if this run was one of the surround stages, so main can stop. */
static int surround_mode(int argc, char **argv, int *rc) {
    const char *stem = NULL, *ls = NULL, *rs = NULL, *dec = NULL;
    const char *mixOut = NULL, *map = NULL;
    const char *mixIn[64]; int nMixIn = 0;
    ps_filt filts[PS_MIX_MAX_FILT]; int nFilt = 0;
    long trim = 0;
    double hz = 7000.0, gain = 0.0, delay = 0.0;
    double ap1 = 0.0, ap2 = 0.0, apw = 1.0;
    int full = 0, flags = 0, bits = 24, quiet = 0;

    for (int i = 1; i < argc; i++) {
        const char *a = argv[i];
        if      (!strcmp(a, "--mix")         && i + 1 < argc) mixOut = argv[++i];
        else if (!strcmp(a, "--map")         && i + 1 < argc) map    = argv[++i];
        else if (!strcmp(a, "--in")          && i + 1 < argc) {
            if (nMixIn >= 64) { fprintf(stderr, "error: too many --in files\n"); *rc = 1; return 0; }
            mixIn[nMixIn++] = argv[++i];
        }
        else if ((!strcmp(a, "--lowpass") || !strcmp(a, "--highpass")) && i + 1 < argc) {
            if (nFilt >= PS_MIX_MAX_FILT) {
                fprintf(stderr, "error: too many filters\n"); *rc = 1; return 0;
            }
            filts[nFilt].kind = !strcmp(a, "--lowpass") ? PS_FILT_LOWPASS : PS_FILT_HIGHPASS;
            filts[nFilt].hz = atof(argv[++i]);
            filts[nFilt].db = 0.0;
            nFilt++;
        }
        else if (!strcmp(a, "--highshelf") && i + 2 < argc) {
            if (nFilt >= PS_MIX_MAX_FILT) {
                fprintf(stderr, "error: too many filters\n"); *rc = 1; return 0;
            }
            filts[nFilt].kind = PS_FILT_HIGHSHELF;
            filts[nFilt].db = atof(argv[++i]);
            filts[nFilt].hz = atof(argv[++i]);
            nFilt++;
        }
        else if (!strcmp(a, "--trim")        && i + 1 < argc) trim = atol(argv[++i]);
        else if (!strcmp(a, "--freqsplit")   && i + 1 < argc) stem  = argv[++i];
        else if (!strcmp(a, "--decorrelate") && i + 1 < argc) dec   = argv[++i];
        else if (!strcmp(a, "--ls")          && i + 1 < argc) ls    = argv[++i];
        else if (!strcmp(a, "--rs")          && i + 1 < argc) rs    = argv[++i];
        else if (!strcmp(a, "--hz")          && i + 1 < argc) hz    = atof(argv[++i]);
        else if (!strcmp(a, "--gain")        && i + 1 < argc) gain  = atof(argv[++i]);
        else if (!strcmp(a, "--delay")       && i + 1 < argc) delay = atof(argv[++i]);
        else if (!strcmp(a, "--ap1")         && i + 1 < argc) ap1   = atof(argv[++i]);
        else if (!strcmp(a, "--ap2")         && i + 1 < argc) ap2   = atof(argv[++i]);
        else if (!strcmp(a, "--apwidth")     && i + 1 < argc) apw   = atof(argv[++i]);
        else if (!strcmp(a, "--bits")        && i + 1 < argc) bits  = atoi(argv[++i]);
        else if (!strcmp(a, "--blend")       && i + 1 < argc) full  = !strcmp(argv[++i], "full");
        else if (!strcmp(a, "--stage-round"))  flags |= PS_FX_STAGE_ROUND;
        else if (!strcmp(a, "-q") || !strcmp(a, "--quiet")) quiet = 1;
    }
    if (!stem && !dec && !mixOut) return 1;

    ps_sample_fmt fmt = (bits == 16) ? PS_FMT_PCM16
                      : (bits == 32) ? PS_FMT_PCM32 : PS_FMT_PCM24;
    const char *msg = "unknown error";

    if (mixOut) {
        if (!map || nMixIn == 0) {
            fprintf(stderr, "error: --mix needs at least one --in and a --map\n");
            *rc = 1; return 0;
        }
        if (ps_mix(mixIn, nMixIn, mixOut, map, gain, filts, nFilt,
                   trim > 0 ? (size_t)trim : 0, flags, fmt, &msg)) {
            fprintf(stderr, "error: %s\n", msg); *rc = 1; return 0;
        }
        if (!quiet)
            printf("  mix: %d file%s -> %s%s%s\n", nMixIn, nMixIn == 1 ? "" : "s", mixOut,
                   nFilt ? ", filtered" : "",
                   (flags & PS_FX_STAGE_ROUND) ? ", rounding every stage" : "");
    }
    if (stem) {
        if (!ls || !rs) {
            fprintf(stderr, "error: --freqsplit needs both --ls and --rs\n");
            *rc = 1; return 0;
        }
        if (ps_freqsplit(stem, ls, rs, hz, gain, full, flags, fmt, &msg)) {
            fprintf(stderr, "error: %s\n", msg); *rc = 1; return 0;
        }
        if (!quiet)
            printf("  freqsplit: %g Hz, %+g dB, %s blend%s\n",
                   hz, gain, full ? "full" : "classic",
                   (flags & PS_FX_STAGE_ROUND) ? ", rounding every stage" : "");
    }
    if (dec) {
        if (ps_decorrelate(dec, delay, ap1, ap2, apw, flags, fmt, &msg)) {
            fprintf(stderr, "error: %s\n", msg); *rc = 1; return 0;
        }
        if (!quiet)
            printf("  decorrelate: %g s delay, allpass %g and %g Hz%s\n",
                   delay, ap1, ap2,
                   (flags & PS_FX_STAGE_ROUND) ? ", rounding every stage" : "");
    }
    *rc = 0;
    return 0;
}

int main(int argc, char **argv) {
    if (argc < 2) { usage(); return 1; }
    if (!strcmp(argv[1], "-h") || !strcmp(argv[1], "--help")) { usage(); return 0; }
    if (!strcmp(argv[1], "-v") || !strcmp(argv[1], "--version")) {
        printf("phasesplit %s\n", PS_VERSION);
        return 0;
    }
    if (!strcmp(argv[1], "--selftest")) return ps_selftest();
    if (!strcmp(argv[1], "--sweep"))    return ps_sweep();
    if (!strcmp(argv[1], "--cpuinfo"))  { ps_cpu_report(); return 0; }
    { int rc; if (!surround_mode(argc, argv, &rc)) return rc; }

    const char   *inPath  = argv[1];
    const char   *cenPath = NULL, *sidPath = NULL, *copyTo = NULL;
    int           quiet = 0, bits = -1;
    const char   *isaName = NULL;
    ps_split_cfg  cfg;
    ps_split_default_cfg(&cfg);

    for (int i = 2; i < argc; i++) {
        const char *a = argv[i];
        if      (!strcmp(a, "-c")                             && i + 1 < argc) cenPath = argv[++i];
        else if (!strcmp(a, "-s")                             && i + 1 < argc) sidPath = argv[++i];
        else if ((!strcmp(a, "-w") || !strcmp(a, "--window")) && i + 1 < argc) cfg.window = atoi(argv[++i]);
        else if ((!strcmp(a, "-o") || !strcmp(a, "--overlap")) && i + 1 < argc) cfg.overlap = atoi(argv[++i]);
        else if (!strcmp(a, "--strength")                     && i + 1 < argc) cfg.strength = (float)atof(argv[++i]);
        else if ((!strcmp(a, "-b") || !strcmp(a, "--bits"))   && i + 1 < argc) {
            const char *v = argv[++i];
            if (!strcmp(v, "f32") || !strcmp(v, "float")) bits = 0; else bits = atoi(v);
        }
        else if (!strcmp(a, "-d") || !strcmp(a, "--dual"))  cfg.dual_centre = 1;
        else if (!strcmp(a, "-q") || !strcmp(a, "--quiet")) quiet = 1;
        else if (!strcmp(a, "--threads") && i + 1 < argc)    cfg.threads = atoi(argv[++i]);
        else if (!strcmp(a, "--isa") && i + 1 < argc)       isaName = argv[++i];
        else if (!strcmp(a, "--copy") && i + 1 < argc)      copyTo = argv[++i];
        else { fprintf(stderr, "unknown or incomplete option: %s\n", a); return 1; }
    }

    /* SSE2 by default. It matches the scalar path bit for bit, so a file
     * processed here comes out the same as one processed anywhere else. AVX2
     * is a little quicker but its fused multiply-add rounds differently, which
     * can move the last bit of a 24 bit sample - inaudible, but it would mean
     * two machines no longer produce identical files. */
    if (!isaName || !strcmp(isaName, "sse2"))   ps_cpu_force(PS_ISA_SSE2);
    else if (!strcmp(isaName, "scalar"))        ps_cpu_force(PS_ISA_SCALAR);
    else if (!strcmp(isaName, "avx2"))          ps_cpu_force(PS_ISA_AVX2);
    else if (!strcmp(isaName, "auto"))          ps_cpu_force((ps_isa)-1);
    else { fprintf(stderr, "unknown --isa: %s (want scalar, sse2, avx2 or auto)\n", isaName); return 1; }
    ps_dsp_init();

    if (!cenPath && !sidPath && !copyTo) {
        fprintf(stderr, "nothing to do: give -c and/or -s (or --copy)\n");
        return 1;
    }

    ps_audio a;
    const char *why = NULL;
    if (ps_wav_read(inPath, &a, &why) != 0) {
        fprintf(stderr, "error: %s: %s\n", inPath, why ? why : "read failed");
        return 1;
    }

    ps_sample_fmt outFmt = a.src_fmt;
    if      (bits == 16) outFmt = PS_FMT_PCM16;
    else if (bits == 24) outFmt = PS_FMT_PCM24;
    else if (bits == 32) outFmt = PS_FMT_PCM32;
    else if (bits == 0)  outFmt = PS_FMT_F32;

    if (!quiet) {
        printf("  in     : %s\n", inPath);
        printf("           %u ch, %u Hz, %s, %.2f s\n", a.channels, a.rate,
               ps_fmt_name(a.src_fmt), (double)a.frames / (a.rate ? a.rate : 1));
    }

    int rc = 0;

    if (copyTo) {
        if (ps_wav_write(copyTo, &a, a.src_fmt, &why) != 0) {
            fprintf(stderr, "error: %s: %s\n", copyTo, why ? why : "write failed");
            rc = 1;
        }
        ps_audio_free(&a);
        return rc;
    }

    if (a.channels != 2) {
        fprintf(stderr, "error: %s has %u channels, this needs stereo\n", inPath, a.channels);
        ps_audio_free(&a);
        return 1;
    }

    ps_split *sp = ps_split_create(&cfg);
    if (!sp) {
        fprintf(stderr, "error: bad settings (window must be a power of two >= 64, overlap 2/4/8)\n");
        ps_audio_free(&a);
        return 1;
    }

    const int cch = cfg.dual_centre ? 2 : 1;
    float *cen = cenPath ? (float *)malloc(a.frames * (size_t)cch * sizeof(float)) : NULL;
    float *sid = sidPath ? (float *)malloc(a.frames * 2 * sizeof(float)) : NULL;
    if ((cenPath && !cen) || (sidPath && !sid)) {
        fprintf(stderr, "error: out of memory for the output buffers\n");
        free(cen); free(sid); ps_split_destroy(sp); ps_audio_free(&a);
        return 1;
    }

    if (!quiet) printf("  split  : window %d, overlap %dx, strength %.2f, %s\n",
                       cfg.window, cfg.overlap, (double)cfg.strength, ps_dsp_path());
    if (!quiet) printf("  threads: %d\n", ps_split_threads(sp));

    clock_t t0 = clock();
    ps_split_process(sp, a.data, a.frames, cen, sid);
    clock_t t1 = clock();

    if (!quiet) {
        double secs = (double)(t1 - t0) / (double)CLOCKS_PER_SEC;
        double dur  = (double)a.frames / (a.rate ? a.rate : 1);
        printf("  done   : %.2f s (%.0fx real time)\n", secs, secs > 0.0 ? dur / secs : 0.0);
    }

    if (cenPath) {
        ps_audio o = a;
        o.data = cen; o.channels = (unsigned)cch;
        if (ps_wav_write(cenPath, &o, outFmt, &why) != 0) {
            fprintf(stderr, "error: %s: %s\n", cenPath, why ? why : "write failed"); rc = 1;
        } else if (!quiet) printf("  centre : %s\n", cenPath);
    }
    if (sidPath) {
        ps_audio o = a;
        o.data = sid; o.channels = 2;
        if (ps_wav_write(sidPath, &o, outFmt, &why) != 0) {
            fprintf(stderr, "error: %s: %s\n", sidPath, why ? why : "write failed"); rc = 1;
        } else if (!quiet) printf("  sides  : %s\n", sidPath);
    }

    free(cen); free(sid);
    ps_split_destroy(sp);
    ps_audio_free(&a);
    return rc;
}
