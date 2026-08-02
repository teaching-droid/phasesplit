#include "fx.h"
#include <math.h>

/* The rest of the program is built with fast floating point, which lets the
 * compiler reassociate arithmetic. That is fine for the spectral work, where
 * the last bit does not matter, but not here: a biquad feeds its own output
 * back in, so any reordering compounds sample after sample and the claim that
 * an exact match would quietly stop being true. */
#ifdef _MSC_VER
#pragma float_control(precise, on)
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* Cookbook two pole highpass. Q is fixed at 1/sqrt(2), which is the
 * Butterworth case - flat in the passband with no peak at the corner. */
ps_biquad ps_biquad_highpass(double f0, double rate) {
    const double Q  = 0.70710678118654752440;
    double w  = 2.0 * M_PI * f0 / rate;
    double c  = cos(w), s = sin(w);
    double al = s / (2.0 * Q);
    double a0 = 1.0 + al;
    ps_biquad q;
    q.b0 = ((1.0 + c) / 2.0) / a0;
    q.b1 = (-(1.0 + c))      / a0;
    q.b2 = ((1.0 + c) / 2.0) / a0;
    q.a1 = (-2.0 * c)        / a0;
    q.a2 = (1.0 - al)        / a0;
    return q;
}

/* The same shape the other way up, and the same Q. Used to build the low
 * frequency channel, where two in series give a 24 dB per octave slope. */
ps_biquad ps_biquad_lowpass(double f0, double rate) {
    const double Q  = 0.70710678118654752440;
    double w  = 2.0 * M_PI * f0 / rate;
    double c  = cos(w), s = sin(w);
    double al = s / (2.0 * Q);
    double a0 = 1.0 + al;
    ps_biquad q;
    q.b0 = ((1.0 - c) / 2.0) / a0;
    q.b1 = (1.0 - c)         / a0;
    q.b2 = ((1.0 - c) / 2.0) / a0;
    q.a1 = (-2.0 * c)        / a0;
    q.a2 = (1.0 - al)        / a0;
    return q;
}

/* Cookbook allpass: flat magnitude, phase turned around f0. The width is a
 * bandwidth in Hz, so Q is the centre divided by it. A width of 1 Hz at 200 Hz
 * therefore means Q = 200, which is an extremely narrow phase twist - that is
 * what the earlier processing asked for, and it is deliberate. */
ps_biquad ps_biquad_allpass(double f0, double width_hz, double rate) {
    double w  = 2.0 * M_PI * f0 / rate;
    double c  = cos(w), s = sin(w);
    /* alpha = sin(w) / (2Q) with Q = f0/width, written out in one step. Going
     * via Q rounds twice - once forming f0/width, once dividing by it - and at
     * a Q in the hundreds those two roundings are enough to move the last bit
     * of the output. Same algebra, one rounding. */
    double al = (width_hz > 0.0)
              ? s * width_hz / (2.0 * f0)
              : s / (2.0 * 0.70710678118654752440);
    double a0 = 1.0 + al;
    ps_biquad q;
    q.b0 = (1.0 - al) / a0;
    q.b1 = (-2.0 * c) / a0;
    q.b2 = 1.0;
    q.a1 = (-2.0 * c) / a0;
    q.a2 = (1.0 - al) / a0;
    return q;
}

/* High shelf: unity at DC, the asked-for gain at the top, one corner.
 *
 * This one is FIRST order, which took some finding. Fitting it as a two pole
 * shelf never converged because there is no second order in it at all - the
 * measured response has a pole and a zero sitting on the same spot, cancelling
 * out, with a single pole and zero left doing the work.
 *
 * Written in the bilinear variable k = tan(corner/2), the two corners sit
 * either side of the nominal frequency by the square root of the gain: the
 * pole at K*sqrt(g) and the zero at K/sqrt(g). Their product is K squared, so
 * the shelf is centred on the frequency asked for however much gain is
 * applied. Checked against nine gain and frequency combinations, and the
 * corners land on those values to every digit that could be measured. */
ps_biquad ps_biquad_highshelf(double gain_db, double f0, double rate) {
    double g  = pow(10.0, gain_db / 20.0);
    double K  = tan(M_PI * f0 / rate);
    double rg = sqrt(g);
    double kp = K * rg;         /* pole corner  */
    double kz = K / rg;         /* zero corner  */
    double P  = (1.0 - kp) / (1.0 + kp);
    double Z  = (1.0 - kz) / (1.0 + kz);
    /* Scaled so DC comes through untouched; the top then lands on g by
     * construction rather than by adjustment. */
    double b0 = (1.0 - P) / (1.0 - Z);
    ps_biquad q;
    q.b0 = b0;
    q.b1 = -b0 * Z;
    q.b2 = 0.0;
    q.a1 = -P;
    q.a2 = 0.0;
    return q;
}

double ps_round32(double v) {
    /* Away from zero, then clipped to the signed 32 bit range. */
    double r = (v < 0.0) ? v - 0.5 : v + 0.5;
    if (r >= 2147483647.0) return 2147483647.0;
    if (r <= -2147483648.0) return -2147483648.0;
    return (double)(long long)r;   /* truncation completes the rounding */
}

void ps_biquad_run(const ps_biquad *c, ps_biquad_state *st,
                   double *buf, size_t n, int flags) {
    double b0 = c->b0, b1 = c->b1, b2 = c->b2, a1 = c->a1, a2 = c->a2;
    double x1 = st->x1, x2 = st->x2, y1 = st->y1, y2 = st->y2;
    if (flags & PS_FX_STAGE_ROUND) {
        /* The rounding happens on the way out only. What the filter feeds back
         * into itself stays at full precision, and it has to: these allpasses
         * run at a Q in the hundreds, so a rounded value in the feedback path
         * does not stay a half-bit error, it is amplified every sample until
         * the output is visibly wrong. */
        for (size_t i = 0; i < n; i++) {
            double x = buf[i];
            double y = b0 * x + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2;
            buf[i] = ps_round32(y);
            x2 = x1; x1 = x;
            y2 = y1; y1 = y;
        }
    } else {
        for (size_t i = 0; i < n; i++) {
            double x = buf[i];
            double y = b0 * x + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2;
            buf[i] = y;
            x2 = x1; x1 = x;
            y2 = y1; y1 = y;
        }
    }
    st->x1 = x1; st->x2 = x2; st->y1 = y1; st->y2 = y2;
}

void ps_fx_gain(double *buf, size_t n, double db, int flags) {
    double g = pow(10.0, db / 20.0);
    if (flags & PS_FX_STAGE_ROUND) {
        for (size_t i = 0; i < n; i++) buf[i] = ps_round32(buf[i] * g);
    } else {
        for (size_t i = 0; i < n; i++) buf[i] *= g;
    }
}

double ps_regrid24(double v) {
    /* 32 bit grid down to 24 bit and back, rounding to nearest. */
    double q = floor(v / 256.0 + 0.5);
    if (q >  8388607.0) q =  8388607.0;
    if (q < -8388608.0) q = -8388608.0;
    return q * 256.0;
}

size_t ps_fx_delay_samples(double seconds, unsigned rate) {
    /* Rounding to nearest, but the ties have to be broken on the true product
     * rather than the rounded one, and the two disagree.
     *
     * 0.015 s at 44100 is 661.49999999999997 exactly, which is just under the
     * half way point, so it belongs at 661. Multiplied in double it lands on
     * 661.5 on the nose and would go to 662. Meanwhile 0.005 s is genuinely a
     * shade over 220.5 and does belong at 221. Both look like exact halves
     * after the multiply, so the rounded product alone cannot tell them apart.
     *
     * fma gives the part of the product that the multiply discarded, which is
     * what settles it. Being one sample out here would shift a whole surround
     * channel against the others. */
    double v = seconds * (double)rate;
    if (v <= 0.0) return 0;
    double err  = fma(seconds, (double)rate, -v);   /* exact remainder */
    double base = floor(v);
    double frac = v - base;
    if (frac > 0.5) return (size_t)base + 1;
    if (frac < 0.5) return (size_t)base;
    return (err > 0.0) ? (size_t)base + 1 : (size_t)base;
}
