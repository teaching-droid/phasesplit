/* Runtime CPU feature detection.
 *
 * The binary is built without any -arch/-m switch, so it starts on anything
 * with plain SSE2 (which every x86-64 machine has). The wider paths are picked
 * here at start-up instead, which means one executable runs everywhere rather
 * than either crashing on older hardware or leaving performance unused on new.
 */
#ifndef PS_CPU_H
#define PS_CPU_H

typedef enum {
    PS_ISA_SCALAR = 0,
    PS_ISA_SSE2,
    PS_ISA_AVX2      /* AVX2 with FMA */
} ps_isa;

/* Highest usable instruction set, taking into account that the OS must also
 * agree to preserve the wide registers across a context switch. */
ps_isa      ps_cpu_isa(void);
const char *ps_isa_name(ps_isa i);

/* Force a lower path than detected; used by the tests to prove every path
 * agrees. Passing something the CPU cannot do is ignored. */
void        ps_cpu_force(ps_isa i);

/* Human readable feature list, for --cpuinfo. */
void        ps_cpu_report(void);

#endif
