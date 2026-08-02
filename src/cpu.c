#include "cpu.h"

#include <stdio.h>
#include <string.h>

#if defined(_MSC_VER)
  #include <intrin.h>
#elif defined(__x86_64__) || defined(__i386__)
  #include <cpuid.h>
  #include <immintrin.h>
#endif

static int    g_probed = 0;
static ps_isa g_native = PS_ISA_SCALAR;
static ps_isa g_forced = (ps_isa)-1;

#if defined(__x86_64__) || defined(_M_X64) || defined(__i386__) || defined(_M_IX86)

static void cpuid_regs(int leaf, int sub, unsigned r[4]) {
#if defined(_MSC_VER)
    int t[4];
    __cpuidex(t, leaf, sub);
    r[0] = (unsigned)t[0]; r[1] = (unsigned)t[1];
    r[2] = (unsigned)t[2]; r[3] = (unsigned)t[3];
#else
    unsigned a, b, c, d;
    __cpuid_count(leaf, sub, a, b, c, d);
    r[0] = a; r[1] = b; r[2] = c; r[3] = d;
#endif
}

static unsigned long long read_xcr0(void) {
#if defined(_MSC_VER)
    return _xgetbv(0);
#else
    unsigned lo, hi;
    __asm__ __volatile__("xgetbv" : "=a"(lo), "=d"(hi) : "c"(0));
    return ((unsigned long long)hi << 32) | lo;
#endif
}

static void probe(void) {
    unsigned r[4];
    cpuid_regs(0, 0, r);
    const unsigned maxLeaf = r[0];

    g_native = PS_ISA_SCALAR;
    if (maxLeaf < 1) return;

    cpuid_regs(1, 0, r);
    const int hasSSE2   = (r[3] & (1u << 26)) != 0;
    const int hasOSXSAVE= (r[2] & (1u << 27)) != 0;
    const int hasAVX    = (r[2] & (1u << 28)) != 0;
    const int hasFMA    = (r[2] & (1u << 12)) != 0;

    if (hasSSE2) g_native = PS_ISA_SSE2;

    /* AVX registers are only usable if the operating system saves them, which
     * it advertises through XCR0. Skipping this check is how programs crash on
     * otherwise capable machines. */
    int osAvx = 0;
    if (hasOSXSAVE && hasAVX) {
        unsigned long long x = read_xcr0();
        osAvx = ((x & 0x6) == 0x6);      /* XMM and YMM state saved */
    }

    if (osAvx && hasFMA && maxLeaf >= 7) {
        cpuid_regs(7, 0, r);
        if (r[1] & (1u << 5)) g_native = PS_ISA_AVX2;
    }
}

/* AVX-512 is reported separately: it is detected and printed, but the DSP does
 * not use it. See the note in ps_cpu_report for why. */
static int have_avx512f(void) {
    unsigned r[4];
    cpuid_regs(0, 0, r);
    if (r[0] < 7) return 0;
    cpuid_regs(1, 0, r);
    if (!(r[2] & (1u << 27))) return 0;              /* OSXSAVE */
    unsigned long long x = read_xcr0();
    if ((x & 0xE6) != 0xE6) return 0;                /* opmask + ZMM state saved */
    cpuid_regs(7, 0, r);
    return (r[1] & (1u << 16)) != 0;                 /* AVX512F */
}

void ps_cpu_report(void) {
    unsigned r[4];
    char brand[49];
    memset(brand, 0, sizeof(brand));
    cpuid_regs(0x80000000u, 0, r);
    if (r[0] >= 0x80000004u) {
        for (int i = 0; i < 3; i++) {
            cpuid_regs(0x80000002u + (unsigned)i, 0, r);
            memcpy(brand + i * 16, r, 16);
        }
    }
    printf("  cpu      : %s\n", brand[0] ? brand : "unknown");
    printf("  detected : %s\n", ps_isa_name(ps_cpu_isa()));
    printf("  avx512f  : %s\n", have_avx512f() ? "yes (present, not used - see below)" : "no");
    printf("\n");
    printf("  The wide paths are picked at run time, so this same binary also\n");
    printf("  starts on a machine that only has SSE2.\n");
}

#else   /* not x86 */

static void probe(void) { g_native = PS_ISA_SCALAR; }

void ps_cpu_report(void) {
    printf("  cpu      : not x86\n");
    printf("  detected : %s\n", ps_isa_name(ps_cpu_isa()));
}

#endif

ps_isa ps_cpu_isa(void) {
    if (!g_probed) { probe(); g_probed = 1; }
    if (g_forced != (ps_isa)-1 && g_forced <= g_native) return g_forced;
    return g_native;
}

void ps_cpu_force(ps_isa i) {
    if (!g_probed) { probe(); g_probed = 1; }
    g_forced = i;
}

const char *ps_isa_name(ps_isa i) {
    switch (i) {
        case PS_ISA_SCALAR: return "scalar";
        case PS_ISA_SSE2:   return "SSE2";
        case PS_ISA_AVX2:   return "AVX2+FMA";
    }
    return "?";
}
