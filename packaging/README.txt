phasesplit
==========

Splits a stereo file into the part that is common to both channels and the
part that is not, and builds output channels as weighted sums of input ones.

The two outputs always add back up to the input, sample for sample. The
centre is worked out and the sides are then taken as the remainder rather
than worked out separately, so nothing is invented and nothing is lost.
Anything panned hard to one side gives a centre of exactly zero.

Run phasesplit --version to see which build this is.


Installing
----------

There is nothing to install. The program is a single file and needs no
libraries beyond the ones the system already has. Put it wherever suits and
run it from a terminal.


Splitting a file
----------------

    phasesplit input.wav -c centre.wav -s sides.wav

Either output can be left out if only one is wanted. The most useful
settings:

    -w, --window N    transform length, a power of two (default 4096)
    -o, --overlap N   overlapping frames, 2 4 or 8 (default 4)
    -b, --bits N      output depth: 16, 24, 32 or f32
        --strength X  level of the extracted centre, 0 to 2 (default 1)
        --threads N   worker threads (default one per core)

The defaults are chosen to be sensible on real material and there is rarely
a reason to change them.


Combining channels
------------------

Output channels can be built as weighted sums of input channels, filtered
and trimmed, in a single pass:

    phasesplit --mix out.wav --in a.wav --in b.wav --map "1*0.5+3*0.5;2"

Inputs are laid end to end and numbered from 1, so with a stereo file
followed by a mono one, channels 1 and 2 are the stereo pair and channel 3
is the mono. Output channels are separated by ';', the terms within one by
'+', and each term is either N or N*weight.

Every weight is written out in full, deliberately. Leaving them off and
having the program quietly divide by the number of terms is an easy way to
end up with a channel several dB down without noticing.

Each output can then be shaped, in this order:

    --lowpass HZ / --highpass HZ    two pole, may be given more than once
    --highshelf DB HZ               tilt the top of the band
    --gain DB
    --trim N                        keep the first N samples


Everything else
---------------

    phasesplit --help

lists every option, including the surround stages.


Checking the build
------------------

    phasesplit --selftest

runs 64 checks on this machine and prints what it found: the transform
against a direct reference, energy conservation, exact reconstruction at
every supported window and overlap, agreement between the plain and the
vectorised code paths, filter shapes, and channel map parsing.

    phasesplit --cpuinfo

shows which instruction sets this processor offers. The program picks one at
run time, so the same file works on an old machine and a new one.


Licence
-------

MIT. The full text is in the LICENSE file next to this one.

Source, issues and newer versions:
https://github.com/teaching-droid/phasesplit
