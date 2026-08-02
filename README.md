# phasesplit

[![build](https://github.com/teaching-droid/phasesplit/actions/workflows/build.yml/badge.svg)](https://github.com/teaching-droid/phasesplit/actions/workflows/build.yml)

Splits a stereo file into the part that is common to both channels and the part
that is not, and does the channel mixing and filtering that usually goes with
it. Written in C with no dependencies beyond the standard library.

Prebuilt binaries for Linux, macOS and Windows are on the
[releases page](https://github.com/teaching-droid/phasesplit/releases).

The two outputs always add back up to the input, sample for sample. The centre
is estimated, and the sides are then taken as the remainder rather than
estimated separately, so nothing is invented and nothing is lost. Content
panned hard to one side yields a centre of exactly zero.

## Building

Windows, with the Visual Studio C++ tools installed:

    build.bat

Anywhere else:

    cmake -B build -DCMAKE_BUILD_TYPE=Release
    cmake --build build
    ctest --test-dir build

There is no `-march` or `/arch` switch in either build. The wide instruction
paths are chosen at run time after asking the processor what it supports, so
one binary runs anywhere from an old SSE2 machine upwards.

## Use

Split a file:

    phasesplit input.wav -c centre.wav -s sides.wav

Useful options:

    -w, --window N    transform length, power of two (default 4096)
    -o, --overlap N   overlapping frames, 2 4 or 8 (default 4)
    -b, --bits N      output depth: 16, 24, 32 or f32
        --isa NAME    scalar, sse2, avx2 or auto (default sse2)
        --threads N   worker threads (default one per core)
        --selftest    run the built-in checks
        --sweep       measure the window and overlap trade-off
        --cpuinfo     show what this machine supports

It also builds output channels as weighted sums of input channels, with
optional filtering, in one pass:

    phasesplit --mix out.wav --in a.wav --in b.wav --map "1*0.5+3*0.5;2"

Inputs are laid end to end and numbered from 1, so with a stereo file followed
by a mono one, channels 1 and 2 are the stereo pair and channel 3 is the mono.
Output channels are separated by `;`, terms within one by `+`, and each term is
either `N` or `N*weight`.

Every weight is written out in full. Specification languages that let weights
be left off and then quietly divide by the number of terms are a good way to
end up with a channel several dB down without noticing.

Each output can then be shaped, in this order:

    --lowpass HZ / --highpass HZ    two pole, repeatable
    --highshelf DB HZ               tilt the top of the band
    --gain DB
    --trim N                        keep the first N samples

## Notes

A few things that were settled by measurement rather than assumption, and are
worth knowing if you change the code.

**The centre cannot be louder than the quieter channel.** The estimate works
out to the square root of the real part of the cross spectrum, which is a
geometric mean of the two channels, and a geometric mean sits above the smaller
of two unequal values. Without a bound, every strongly panned bin hands too
much to the centre. Adding it was worth 1.5 to 2.3 dB on dense material, where
hard panned parts are everywhere, and helped on sparse material too.

**Window length barely matters on real material.** Across a sixteen fold range
the difference is around 0.2 dB. A sweep using steady tones suggests otherwise,
but a steady tone is the easiest possible case for a long transform and nothing
in music resembles one. 4096 with fourfold overlap is the default because it is
as good as anything longer while costing a quarter of the latency.

**SSE2 is the default rather than AVX2.** The SSE2 path is bit for bit
identical to the scalar one. AVX2 is a little quicker but its fused multiply
add rounds once where the others round twice, which can move the last bit of a
24 bit sample. Inaudible, but it would mean two machines no longer produce
identical files. Pass `--isa avx2` to opt in.

**Filters keep full precision in their feedback.** Rounding an output sample
before it is fed back sounds harmless and is not: at the high Q values used for
phase work the error is amplified every sample rather than staying put.

**`--strength` is a level control only.** It scales the extracted centre, and
the sides keep whatever it leaves, so the two still sum to the input. It does
not change how the centre is picked out.

## Self test

    phasesplit --selftest

64 checks covering the transform against a direct reference, energy
conservation, exact reconstruction at every supported window and overlap,
agreement between the scalar and vector paths, filter shapes, and the parsing
of channel maps.

## Licence

MIT. See [LICENSE](LICENSE).
