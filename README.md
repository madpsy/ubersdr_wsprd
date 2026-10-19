# wsprd and jt9

WSJT-X's WSPR decoder, `wsprd`, and its `jt9` for FST4W, FT8 and FT4 (and
FT2, as WSJT-X improved decodes it), rewritten in C++ to take a fraction of
the CPU.

**A drop-in replacement for WSJT-X's `wsprd`.** It takes exactly the same
arguments, plus one more (`-n`), and gives the same output and files, byte
for byte. Put it where WSJT-X's `wsprd` is and whatever runs `wsprd` runs
this one, with nothing to change.

- **The same decodes, to the bit.** The decoder is WSJT-X's own `wsprd.c`
  (with its OSD, `osdwspr.f90`) at WSJT-X 3.2.0-rc1, commit `567ad29ce`,
  kept line for line where it could be, so the two read side by side. Every
  decode's SNR, DT, frequency, drift and sync are the reference's from the
  same samples, and so is everything it prints and writes: stdout,
  `ALL_WSPR.TXT`, `wspr_spots.txt`, `hashtable.txt` and the `-c` file.
- **A tenth of the CPU on x86-64, a sixth on a Raspberry Pi 5.** The same
  arithmetic, done in a cheaper order and many lanes at once (AVX2 on x86-64,
  NEON on arm64).
- **One option more, `-n`**: no drift search, and each signal's spread.
- Linux (x86-64, arm64) and Windows.
- The original works are the WSJT-X project's, and this is released under
  the same licence, the GNU GPL v3 (see [Credits and licence](#credits-and-licence)).

It decodes what WSJT-X's `wsprd` decodes, WSPR-2 and WSPR-15. FST4W is
`jt9`'s (below).

**And `jt9` for FST4W, FT8, FT4 and FT2.** WSJT-X decodes FST4W (FST4W-120
to -1800), FT8 and FT4 with its `jt9`; this `jt9` does the same from WAV
files, as `jt9 -p 120 --fst4w ... file.wav`, `jt9 -8 ... file.wav` and
`jt9 -5 ... file.wav` do: the same command line, the same output and files
(`fst4_decodes.dat`, `decoded.txt`, `fst4w_calls.txt`), byte for byte, in
about a seventh of the CPU for FST4W, under half for FT8 and about a third
for FT4. FT2, which only WSJT-X improved has, is `jt9 --ft2 ... file.wav`
as improved's `jt9` runs it: improved's FT2 code on the FT4 decoder above,
in a quarter to a sixth of improved's CPU on x86-64 (about a third on a
Pi 5). FST4W, FT8, FT4 and FT2 so far, more
of `jt9`'s modes to follow. Linux (x86-64, arm64) and Windows. For FT8,
FT4 and FT2 it also runs as a back end over shared memory (`jt9 -s KEY`),
as WSJT-X 2.7 and WSJT-X improved run their `jt9`, so a program driving
`jt9` that way (jt9-decode, for one) can use it: see
[jt9 -s](#jt9--s-shared-memory). It does not replace `jt9` inside WSJT-X
itself, which runs it for every mode: see [jt9: FST4W](#jt9-fst4w),
[jt9: FT8](#jt9-ft8), [jt9: FT4](#jt9-ft4) and [jt9: FT2](#jt9-ft2).

Originally designed for UberSDR Skimmer.

## Download

Nothing to build: each release has a ready-made binary. Pick the file for
your operating system and CPU from the
[latest release](https://github.com/madpsy/ubersdr_wsprd/releases/tag/latest):

| File | For |
|---|---|
| `wsprd-linux-x86_64` | Linux on an x86-64 CPU with AVX2 |
| `wsprd-linux-aarch64` | Linux on 64-bit ARM (NEON): a Raspberry Pi 3, 4, 5 or Zero 2 W on a 64-bit OS, for one |
| `wsprd-windows-x86_64.exe` | 64-bit Windows 10 or 11 on an x86-64 CPU with AVX2; nothing else to install |
| `jt9-linux-x86_64` | `jt9` (FST4W, FT8, FT4 and FT2) for Linux on an x86-64 CPU with AVX2 |
| `jt9-linux-aarch64` | `jt9` (FST4W, FT8, FT4 and FT2) for Linux on 64-bit ARM |
| `jt9-windows-x86_64.exe` | `jt9` (FST4W, FT8, FT4 and FT2) for 64-bit Windows 10 or 11 on an x86-64 CPU with AVX2 |

**The CPU: the prebuilt binaries need AVX2 on x86-64 (Linux and Windows)
and NEON on ARM.** The decoder's searches run eight lanes at a time with
these vector instructions:

- **x86-64: AVX2** (with FMA and BMI2, the x86-64-v3 level). Intel: Core
  i3/i5/i7 from the 4th generation (Haswell, 2013) on, and Xeons from E3/E5
  v3; Pentium, Celeron and Atom only from Alder Lake-N (N100, N200, N305,
  2023), not the earlier ones (Gemini Lake, Jasper Lake and the like). AMD:
  every Ryzen, EPYC and Threadripper (Zen, 2017 on), and Excavator APUs
  (2015). On Linux, `grep -o -m1 avx2 /proc/cpuinfo` prints `avx2` if yours
  has it. Without it the program stops at once ("Illegal instruction"):
  build it yourself instead, without AVX2 (see [Building it](#building-it)).
- **ARM: NEON**, which every 64-bit ARM CPU has: on a Raspberry Pi 3, 4, 5,
  400, 500, Compute Module 3, 4 or 5, or Zero 2 W, running a 64-bit OS
  (Raspberry Pi OS 64-bit, `uname -m` says `aarch64`). Not the original Pi,
  Pi 2 (v1.1) or Pi Zero / Zero W, nor a 32-bit OS (`armv7l`): build it
  yourself there (see [Building it](#building-it)), without the vector code.

Built without the vector code, it decodes the same, only slower.

The release `latest` always holds the newest, so its links stay the same,
e.g. `https://github.com/madpsy/ubersdr_wsprd/releases/download/latest/wsprd-linux-x86_64`;
each version is also kept as a release of its own (`vYYYY.MM.DD`), and
`wsprd --version` says which you have.

Each is one file with FFTW built in: nothing else to install (`jt9` on Linux also
uses glibc's vector maths library, `libmvec`, part of the C library itself).

Then rename it `wsprd` (`wsprd.exe` on Windows) and put it in place of
WSJT-X's: whatever ran WSJT-X's `wsprd` runs this one, with nothing else to
change. Keep WSJT-X's as `wsprd.orig` if you like, to go back. For example,
where WSJT-X is installed from a package:

```
curl -LO https://github.com/madpsy/ubersdr_wsprd/releases/download/latest/wsprd-linux-x86_64
chmod +x wsprd-linux-x86_64
sudo mv /usr/bin/wsprd /usr/bin/wsprd.orig
sudo cp wsprd-linux-x86_64 /usr/bin/wsprd
wsprd --version
```

On Windows, WSJT-X's is `wsprd.exe` in its `bin` folder (by default
`C:\WSJT\wsjtx\bin`). Or leave WSJT-X's where it is and point your program
at this one.

An update of WSJT-X puts its own `wsprd` back: copy this one in again after.

## Using it

Exactly as WSJT-X's `wsprd`: the same arguments, and `-n` besides.

```
wsprd [options...] infile        (infile: a .wav or .c2)

  -a <path>  path to writeable data files, default="."
  -B         disable block demodulation - use single-symbol noncoherent demod
  -c         write .c2 file at the end of the first pass
  -C         maximum number of decoder cycles per bit, default 10000
  -d         deeper search. Slower, a few more decodes
  -e x       transceiver dial frequency error in Hz
  -f x       transceiver dial frequency in MHz
  -H         do not use (or update) the hash table
  -J         use the stack decoder instead of Fano decoder
  -m         decode wspr-15 .wav file
  -o n       (0<=n<=5), decoding depth for OSD, default is disabled
  -q         quick mode - doesn't dig deep for weak signals
  -s         single pass mode, no subtraction (same as original wsprd)
  -v         verbose mode (shows dupes)
  -w         wideband mode - decode signals within +/- 150 Hz of center
  -z x       fano metric table bias, default is 0.45
  -n         no drift, does not attempt to follow a drifting frequency
  -0         streaming mode: PCM audio on stdin, decodes as NDJSON on stdout
```

For example, a 20 m period, decoded deep:

```
./wsprd -c -C 500 -o 4 -d -f 14.0956 261011_1200.wav
```

The date and time of each decode come from the file's name (`YYMMDD_HHMM`),
as with WSJT-X's. WSPR-2 and WSPR-15 (`-m`) only, as WSJT-X's `wsprd`:
FST4W's periods go to WSJT-X's `jt9`.

**32-bit float `.wav` files too.** A mono file of 32-bit IEEE float samples
(`WAVE_FORMAT_IEEE_FLOAT`, or `WAVE_FORMAT_EXTENSIBLE` with that subformat,
its chunks read wherever they are, as sox and ka9q-radio's `radiod` write
them) is decoded from its samples as they are, at full scale 1.0, with no
rounding to 16 bits and no need to normalize it first: a 16-bit sample `s`
is `s/32768` of that scale, so a float file of the same samples decodes byte
for byte as the 16-bit file does (`scripts/float-input-check.sh`). The level
does not matter: an SDR's unscaled output, its noise at -135 dBFS, decodes
as it would at any other level. Any other file is read as WSJT-X's `wsprd`
reads it, 16-bit samples after a 44-byte header, whatever its header says.
`-0` (streaming) stays 16-bit.

### `-n`: no drift, and the signal's spread

With `-n` the decoder does not search for a drifting signal (each candidate's
drift is held at 0 in every pass), and each line of `ALL_WSPR.TXT` ends in
one more column, the signal's spread in Hz (`%6.3f`): how wide it arrived,
which is what the path did to it. Without `-n`, nothing differs from WSJT-X's
`wsprd`.

The spread is measured as WSJT-X's `jt9` measures FST4W's (`dopspread` in
`lib/fst4_decode.f90`, its w50): the samples times the decoded signal's
conjugate, transformed, and the width holding the middle half of the power
within 1 Hz, less the noise beside it. Its bins are FST4W-120's (1/240 Hz;
WSPR-15's 1/1800 Hz). A steady signal reads about 0.008 Hz; a Gaussian
Doppler spread of 0.02 to 0.5 Hz reads within its noise of the 1.35 sigma
expected, at SNRs down to -25 dB. It costs about 0.7 ms a decode.

`wsprd --version` (alone) prints the release it came from, where WSJT-X's
would refuse `--version`.

## jt9: FST4W

Run it as WSJT-X's `jt9` is run on WAV files, for FST4W:

```
jt9 -a DIR -p 120 --fst4w -f 1500 -F 100 -d 3 YYMMDD_HHMM.wav [more.wav ...]
```

`-p` the period in seconds (FST4W-120 to -1800: 120, 300, 900, 1800; and
15, 30 or 60), `--fst4w` (or `-W`; `-Y` to print the hash22 values), `-f`
and `-F` the frequency and tolerance (at most 100 Hz for FST4W), `-d` the
depth (1 to 3), `-a` and `-t` the data and temporary directories, `-c` and
`-x` the calls fst4_decode is given (`K1ABC` and `W9XYZ` unless given), `-X`
the noise blanker's flags, `-q` no `<DecodeFinished>`. Several files in one
run share one hash table and calls list, as with WSJT-X's. Everything is
WSJT-X's: the decodes, every field of `fst4_decodes.dat` (with a `decdata`
file in the data directory), `decoded.txt`, `fst4w_calls.txt`, and its
messages, down to its command line's quirks (`jt9 -h` lists the options).

One option more, `--spread fast|exact|off`: the spread (w50, with a
`plotspec` file in the directory it runs in). `fast`, the default, measures
it from the same product decimated by 120 before a 120 times smaller FFT,
within 0.0001 Hz of WSJT-X's (the printed w50 the same on 463 of 465 test
decodes) and several times cheaper; `exact` is WSJT-X's own, to the bit; `off`
measures none.

`-c` and `-x` take standard calls only for FST4W.

`jt9 --fst4w` takes a mono 32-bit float `.wav` as well, as `wsprd` does
(above): its samples are decoded as they are, scaled to the 16-bit full
scale, with nothing rounded or clipped, and a float file of the same samples
as a 16-bit one decodes byte for byte as it does. Two things differ for a
float file because its level is its own: only an all-zero file counts as
silent (a 16-bit file whose first 15 s have an RMS under 0.5 is skipped, but
an unscaled SDR's quiet band sits far below that), and the noise blanker's
histogram puts |x| of 32768 and over in its top bin, which no 16-bit sample
reaches. FT8, FT4 and FT2 take float files too (below).

## jt9: FT8

Run it as WSJT-X's `jt9` is run on WAV files, for FT8:

```
jt9 -a DIR -8 -d 2 YYMMDD_HHMMSS.wav [more.wav ...]
```

`-8` (or `--ft8`) FT8, `-d` the depth (1 to 3, default 1, as WSJT-X's),
`-L` and `-H` the lowest and highest frequency searched (200 and 4000 Hz),
`-f` and `-F` the QSO frequency and its tolerance (1500 and 20 Hz), `-c`
and `-x` my call and his call (`K1ABC` and `W9XYZ` unless given, `b` for
none), `-g` his grid, `-G` my grid, `-Q` the QSO progress (0 to 5, for the
a-priori passes on my call), `-X` the contest mode (NA VHF, EU VHF, Field
Day, RTTY Roundup, WW Digi, Hound), `-p` how many seconds of each file are
read (at most an FT8 period's), `-a` and `-t` the data and temporary
directories, `-q` no `<DecodeFinished>` (`jt9 -h` lists the options). Several files in
one run are a sequence of periods, as WSJT-X's: its a-priori decoding of
calls heard in earlier periods (a7) and its hash tables run across them.

Everything is WSJT-X's, byte for byte: the decodes and their order, SNR,
DT, frequency and the `a1`..., `?` annotations on stdout, `decoded.txt`,
the messages (`EOF on input file` for a short WAV, `STOP 2`...), the exit
status, down to its command line's quirks. The decoder is WSJT-X's
`ft8_decode.f90` and what it calls (`sync8`, `ft8b`, `decode174_91`,
`osd174_91`, `subtractft8`, `ft8_a7`, `ft8_a8`, `packjt77`...), as WSJT-X
3.2.0-rc1, commit `567ad29ce`, builds it; the maths library's functions it
calls (glibc's, and its vector `cosf`, `sinf` and `log10f`) are carried
in `src/ft8_maths.cpp`, so it decodes the same on every platform.

On arm64 FFTW's NEON code rounds its transforms differently from the x86
code WSJT-X's build uses, and a decode near a threshold can go the other
way: on our test set 6 of 179 cases print a line differently there (a
decode more or less, or a DT or SNR one step apart). On x86-64 Linux and
Windows every case is the reference's to the bit.

FT8's Fox mode (`-X` with Fox's bit) is refused.

`jt9 -8`, `jt9 -5` and `jt9 --ft2` (FT4 and FT2, below) take a mono 32-bit float `.wav` as well,
as `jt9 --fst4w` does: its samples are decoded as they are, scaled to the
16-bit full scale, nothing rounded or clipped, and a float file of the same
samples as a 16-bit one decodes byte for byte as it does. One thing differs
for a float file, as for FST4W's: only an all-zero period counts as silent
(WSJT-X skips a period whose RMS is under 0.5 of a 16-bit step, and an
unscaled SDR's quiet band sits far below that).

### --feed: one jt9 for a stream of files

`--feed` (ours, FT8, FT4 and FT2) reads the files to decode from stdin, one a
line, and decodes each as it arrives, with one decoder for all of them, so
the a-priori decoding of calls heard in earlier periods (a7) and the hash
tables run across them as they do in WSJT-X's own live jt9:

```
inotifywait -m -q -e moved_to --format '%w%f' /var/lib/ka9q-radio/ft8 | jt9 -a DIR -8 -d 2 --feed
```

A line is `PATH`, or `PATH HHMMSS` to give the period's UTC; without it the
UTC is jt9's from the name (`..._HHMMSS.wav`), else ka9q-radio's `pcmrecord`'s
(`YYYYMMDDTHHMMSSZ_...`), else 0. FT8 needs it: its a7 passes tell even
periods from odd ones by it. Each file's decodes and `<DecodeFinished>` are
printed as `jt9 -8 FILE` (or `-5`, `--ft2`) prints them, then `<FileDone> PATH N` (N its
decodes), and stdout is flushed, so a caller feeding many files knows where
each ends. A file that cannot be read is said on stderr, gets `N = 0`, and
the feed goes on. It ends when stdin does.

This is what jt9-decode's `-F` does over WSJT-X's shared memory (`jt9 -s`,
below), without the shared memory, and so it takes float files as well,
which that memory's 16-bit `d2` cannot carry.

## jt9: FT4

Run it as WSJT-X's `jt9` is run on WAV files, for FT4:

```
jt9 -a DIR -5 -d 2 YYMMDD_HHMMSS.wav [more.wav ...]
```

`-5` (or `--ft4`) FT4, with the options FT8 has that WSJT-X's FT4 reads:
`-d` the depth (1 to 3), `-L` and `-H` the frequencies searched, `-f` the
QSO frequency, `-c` and `-x` my call and his call (`b` for none), `-Q` the
QSO progress (the a-priori passes on my and his call), `-X` the contest
mode, `-a`, `-t` and `-q` as for FT8; `-F`, `-g` and `-G` are taken but, as
in WSJT-X's, not used. Each file is one 7.5 s period (21 blocks of 3456
samples, whatever `-p` says, `EOF on input file` when a file is shorter);
several files are a sequence of periods, the hash tables kept across them.

Everything is WSJT-X's, byte for byte, as for FT8: stdout and its `+`
lines, `decoded.txt`, the messages, the exit status. The decoder is WSJT-X's
`ft4_decode.f90` and what it calls (`getcandidates4`, `sync4d`,
`get_ft4_bitmetrics`, `decode174_91`, `osd174_91`, `subtractft4`,
`packjt77`...), as WSJT-X 3.2.0-rc1, commit `567ad29ce`, builds it, sharing
FT8's LDPC decoder, `packjt77` and carried maths functions. Every case is
the reference's on x86-64 Linux and Windows; on arm64 FFTW's NEON code
moves the last bits of its sums, but no printed line of the decoder's 82
test cases differed.

## jt9: FT2

Run it as WSJT-X improved's `jt9` is run on WAV files, for FT2:

```
jt9 -a DIR --ft2 -d 2 YYMMDD_HHMMSS.wav [more.wav ...]
```

`--ft2` FT2 (only the long option: improved's `jt9` lists `-2` but does
not take it, and neither does this one), with FT4's options above. Each
file is one 3.75 s period (10 blocks of 3456 samples, 34,560 in all,
whatever `-p` says, `EOF on input file` when a file is shorter); several
files are a sequence of periods, the hash tables kept across them. It
takes float files and `--feed` as FT8 does (above).

WSJT-X itself has no FT2; WSJT-X improved 3.2.0 decodes it with its FT4
decoder on the audio stretched to twice its length, with FT2's own sizes,
sync threshold, frequency range and output line. That FT2 code is
improved's, line for line; the FT4 decoder under it is the one above
(WSJT-X 3.2.0-rc1, commit `567ad29ce`), not improved's, which added to its
FT4 (an extra a-priori pass, more bit metrics, a slightly different front
end). So the reference is a hybrid: improved's FT2 code built on
`567ad29ce`'s FT4, as a `jt9` of its own (`ft2-jt9`). Against it,
everything is byte for byte on x86-64 Linux and Windows (stdout and its
`*` lines, `decoded.txt`, the messages, the exit status), and no printed
line differed on arm64. Against improved 3.2.0's own `jt9 --ft2`, most
output is the same (every recorded period's); where it differs, it is
improved's FT4 or its older command line and WAV reader: a decode more or
less, `-0.0` for `0.0`, 1 Hz or 1 dB apart, runtime errors where this one
checks its options, and odd WAV files it reads that this one refuses as
WSJT-X's does. No FT2 was heard on the air while this was made, so its
tests are generated FT2 periods and recorded 40 m periods without any.

The other modes of WSJT-X's `jt9` (JT65, Q65, MSK144...) and its stdin
stream are not here yet: they are unknown options to this one, so do not put
it in place of WSJT-X's `jt9` inside WSJT-X. Run it where `jt9 --fst4w`,
`jt9 -8`, `jt9 -5` or `jt9 --ft2` is run on files, or where `jt9 -s` is run
for FT8, FT4 or FT2 (below). `jt9_wisdom.dat` is kept as WSJT-X's `jt9`
keeps it (FFTW's plans, read from `-a`'s directory at the start and written
back at a normal end, so a second run starts at once); `timer.out` is not
written.

## jt9 -s: shared memory

```
jt9 -s KEY [-w patience] [-m threads] [-e path] [-a path] [-t path] [-r path]
```

As WSJT-X 2.7 and WSJT-X improved run their `jt9`: the program that starts
it makes Qt's shared memory under `KEY` (Qt 5's System V memory and
semaphore on Linux, named objects on Windows), laid out as their `dec_data`,
puts each period's 12 kHz samples in `d2` and the decode's parameters in
`params`, and sets `ipc[1]` to 1. `jt9` prints each decode's line and
`<DecodeFinished>` on stdout (and writes `decoded.txt` in `-t`'s directory),
then waits for `ipc[2]` to be 1 before the next; `ipc[1]` 999 ends it. FT8
with `ndiskdat` set gets jt9's early passes (at 41 and 47 half-symbols, then
50); without it, one pass at `nzhsym`, as WSJT-X's live decodes. `-w`, `-m`,
`-e` and `-r` are read but change nothing (FFTW's planning is the
decoders' own).

FT8 (`nmode` 8), FT4 (5) and FT2 (52) only. A request for another mode,
or for FT8's multithreaded decoder (`lmultift8`), its Fox or SuperFox mode or
an EME delay, is said once on stderr and answered with `<DecodeFinished>`
and no decodes, so the caller never waits. WSJT-X 3.x's newer
shared-memory protocol is not taken.

Run by jt9-decode, against WSJT-X 3.1.0's own `jt9` run by it on the same
recorded periods (FT8 and FT4 at every depth: 27 periods, 1,135 decodes), every
decode line was the same.

## Its CPU

Against WSJT-X's own `wsprd`, built from the same commit with WSJT-X's
Release flags on each machine. User CPU over six real 2-minute periods (five
on 40 m, WSJT-X's 20 m sample, all in `recordings/`) with `-C 500 -o 4 -d`,
the median of three runs; every decode the same:

| CPU for 6 periods | WSJT-X's wsprd | this wsprd | |
|---|---:|---:|---:|
| x86-64 (Ryzen 7 5700U) | 13.7 s | 1.5 s | 9.1x |
| Raspberry Pi 5 | 19.6 s | 3.4 s | 5.7x |

Period by period, side by side (user CPU, seconds; decodes the number of
`ALL_WSPR.TXT` lines, identical from both):

| Period | Decodes | x86-64: WSJT-X's | this | | Pi 5: WSJT-X's | this | |
|---|---:|---:|---:|---:|---:|---:|---:|
| 40 m, 19:14 | 3 | 2.13 | 0.22 | 9.5x | 3.03 | 0.52 | 5.8x |
| 40 m, 19:16 | 8 | 2.81 | 0.29 | 9.7x | 3.95 | 0.70 | 5.6x |
| 40 m, 19:18 | 2 | 2.00 | 0.23 | 8.6x | 2.93 | 0.50 | 5.8x |
| 40 m, 19:20 | 4 | 2.06 | 0.23 | 9.0x | 2.98 | 0.52 | 5.7x |
| 40 m, 19:22 | 4 | 1.72 | 0.19 | 9.1x | 2.49 | 0.47 | 5.3x |
| 20 m, WSJT-X's sample | 8 | 3.00 | 0.34 | 8.9x | 4.19 | 0.71 | 5.9x |
| **all six** | 29 | **13.7** | **1.5** | **9.1x** | **19.6** | **3.4** | **5.7x** |

Without AVX2 (built with `-DWSPRD_X86_ARCH=`, for an x86-64 CPU that lacks
it) the decodes are the same, and it still takes half WSJT-X's CPU; with AVX2
(the prebuilt binaries) four times less again. On the same x86-64 machine,
the same six periods:

| Period | WSJT-X's | this, without AVX2 | this, with AVX2 | AVX2's gain |
|---|---:|---:|---:|---:|
| 40 m, 19:14 | 2.13 | 0.96 | 0.25 | 3.9x |
| 40 m, 19:16 | 2.81 | 1.51 | 0.31 | 4.9x |
| 40 m, 19:18 | 2.00 | 0.96 | 0.23 | 4.1x |
| 40 m, 19:20 | 2.06 | 1.01 | 0.24 | 4.2x |
| 40 m, 19:22 | 1.72 | 0.93 | 0.22 | 4.3x |
| 20 m, WSJT-X's sample | 3.00 | 1.34 | 0.35 | 3.8x |
| **all six** | **13.7** | **6.7** | **1.6** | **4.2x** |

About 21 MB while a WSPR period decodes (140 MB for WSPR-15), freed after.
To time it yourself, on any of `recordings/`:

```
time ./build/wsprd -C 500 -o 4 -d -f 7.0386 recordings/40m-wspr-20261010-1914/261010_1914.wav
```

### jt9

Against WSJT-X's `jt9` built from the same commit with WSJT-X's Release
flags on each machine (`-p P --fst4w -f 1500 -F 100 -d 3`, with `plotspec`
and `decdata`; user and system CPU, best of three; WSJT-X's run once a
period, as it is run on files): nine periods, five real 40 m ones (busy with
WSPR, which FST4W tries and fails on), a busy generated FST4W-120 and one
each of FST4W-300, -900 and -1800; every decode the same.

| CPU for 9 periods | WSJT-X's jt9 | this jt9 | |
|---|---:|---:|---:|
| x86-64 (Ryzen 7 5700U) | 32.8 s | 4.4 s | 7.5x |
| Raspberry Pi 5 | 56.4 s | 8.5 s | 6.6x |

With `--spread exact` (WSJT-X's own spread, an FFT of twice the period's
length a decode) 10.8 s and 28.5 s. An FST4W-1800 period takes about 0.5 GB
while it decodes, WSJT-X's 1.1 GB.

FT8, against WSJT-X's FT8 decoder from the same commit with WSJT-X's
Release flags on each machine, period by period as WSJT-X decodes them
(three passes a period), CPU of the decoding alone, best of three, every
decode the same: six busy 40 m periods, six quieter 20 m ones, four of a
busy generated band and one weak generated period, at each depth:

| CPU for 17 periods | WSJT-X's | this | |
|---|---:|---:|---:|
| x86-64 (Ryzen 7 5700U), depth 1 | 7.9 s | 3.9 s | 2.1x |
| x86-64, depth 2 | 24.5 s | 9.7 s | 2.5x |
| x86-64, depth 3 | 62.8 s | 21.8 s | 2.9x |
| Raspberry Pi 5, depth 1 | 13.6 s | 8.5 s | 1.6x |
| Raspberry Pi 5, depth 2 | 46.5 s | 21.8 s | 2.1x |
| Raspberry Pi 5, depth 3 | 125.2 s | 50.4 s | 2.5x |

A busy 40 m period at depth 2 takes 0.67 s on the x86-64 machine (WSJT-X's
1.62 s), 1.53 s on the Pi 5 (3.08 s).

FT4 the same way (one decoder call a 7.5 s period): eight recorded 40 m
periods, eight 20 m ones, two of a crowded generated band and two weak or
empty generated ones, at each depth (a different x86-64 machine from the
tables above):

| CPU for 20 periods | WSJT-X's | this | |
|---|---:|---:|---:|
| x86-64 (Ryzen 9 3900X), depth 1 | 0.50 s | 0.14 s | 3.7x |
| x86-64, depth 2 | 2.08 s | 0.72 s | 2.9x |
| x86-64, depth 3 | 3.70 s | 0.86 s | 4.3x |
| Raspberry Pi 5, depth 1 | 0.90 s | 0.34 s | 2.6x |
| Raspberry Pi 5, depth 2 | 3.97 s | 1.90 s | 2.1x |
| Raspberry Pi 5, depth 3 | 7.26 s | 2.22 s | 3.3x |

Like WSJT-X's, `jt9 -5` reads its FFTs' plans back from `jt9_wisdom.dat`
(in `-a`'s directory); its first run there plans them afresh, about half a
second on x86-64, and every run after starts at once.

FT2 the same way (one decoder call a 3.75 s period, as improved's GUI
makes it), against WSJT-X improved 3.2.0's own `jt9` built with its
Release flags on each machine: sixteen recorded 40 m periods (with no FT2
in them), four of a crowded generated band and five weak, faded or empty
generated ones, at each depth (improved's FT4 finds a few decodes more or
less on them; the hybrid reference above, every decode the same, takes
3.6x, 3.4x and 5.6x this one's CPU on x86-64, 2.6x, 2.4x and 4.4x on the
Pi):

| CPU for 25 periods | improved's jt9 | this | |
|---|---:|---:|---:|
| x86-64 (Ryzen 9 3900X), depth 1 | 1.34 s | 0.32 s | 4.2x |
| x86-64, depth 2 | 2.29 s | 0.54 s | 4.2x |
| x86-64, depth 3 | 4.22 s | 0.70 s | 6.0x |
| Raspberry Pi 5, depth 1 | 2.52 s | 0.85 s | 3.0x |
| Raspberry Pi 5, depth 2 | 4.29 s | 1.48 s | 2.9x |
| Raspberry Pi 5, depth 3 | 8.68 s | 1.84 s | 4.7x |

Run on files (`jt9 --ft2`, each program's whole run, no wisdom file for
either), 239 recorded 40 m periods in one run take 3.5 s against
improved's 14.2 s at depth 2 (4.0x) and 4.4 s against 28.6 s at depth 3
(6.5x); a single file, mostly FFTW's planning, 0.51 s against 0.69 s.

## Building it

Most users only need a prebuilt binary from
[the latest release](https://github.com/madpsy/ubersdr_wsprd/releases/tag/latest) (see
[Download](#download)); building is for another platform, or to change or
check the code.

On Debian, Ubuntu or Raspberry Pi OS:

```
sudo apt install build-essential cmake pkg-config libfftw3-dev git
git clone https://github.com/madpsy/ubersdr_wsprd.git
cd ubersdr_wsprd
./build.sh
```

That builds for the machine it runs on, x86-64 or arm64 (a Pi builds its own
in a few minutes): `build/wsprd` and `build/jt9`, and stripped, in
`build/dist`. Or by hand:

```
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

FFTW is linked into the binary (`libfftw3-dev`'s static library); add
`-DWSPRD_STATIC_FFTW=OFF` to link the shared one instead.

On an x86-64 CPU without AVX2, build for the compiler's baseline instead
(the same decodes, slower):

```
cmake -B build -DCMAKE_BUILD_TYPE=Release -DWSPRD_X86_ARCH=
cmake --build build
```

A 32-bit ARM system (`armv7l`) builds the plain code by itself.

**Windows**, cross-compiled on Linux with Docker (mingw-w64, FFTW built in,
so `wsprd.exe` and `jt9.exe` are one file each that need nothing beside them):

```
./build.sh --windows
```

They land in `build/dist/wsprd-windows-x86_64.exe` and
`jt9-windows-x86_64.exe`, and when `wine` is installed they are checked:
`wsprd.exe` decodes as the Linux build does, and `jt9.exe` prints and writes
what WSJT-X's `jt9` does on every case of `tests/data/fst4w`,
`tests/data/ft8`, `tests/data/ft4` and `tests/data/ft2` (the hybrid's, for
FT2) whose inputs are there (`scripts/jt9-wine-check.sh`; `fort.52`, the
spread's debug spectrum, to within its last printed digit: Windows' maths
library and FFTW's build round differently). FT8's, FT4's and FT2's are
byte for byte; FFTW for Windows is built
with Ubuntu's SIMD set (SSE2 and AVX) for that.

Build with the settings in `CMakeLists.txt`: the decoder's files are
compiled with no fused multiply-adds (`-ffp-contract=off`), as WSJT-X's own
build rounds them, on every platform; that is part of what keeps the decodes
the reference's.

## Checking it against WSJT-X's own

```
./build.sh --check
```

checks the build against WSJT-X's own outputs stored in `tests/data`,
nothing of WSJT-X built, in a few minutes: `scripts/wsprd-stored-check.sh`
runs `wsprd` on every case in `tests/data/wsprd/cases` as WSJT-X's was run
(stdout, the files written and the exit statuses byte for byte), and
`scripts/jt9-wine-check.sh` runs `jt9` on every stored FST4W, FT8, FT4 and
FT2 case whose inputs are here (the same script checks `jt9.exe` under
wine). `scripts/float-input-check.sh` writes each 16-bit recording
again as 32-bit float, the same samples, and requires `wsprd` (under four
option sets) and `jt9 --fst4w` to leave the same files from it, byte for byte.
`scripts/ft-float-feed-check.sh` decodes each FT8, FT4 and FT2 recording set as
one run (one decoder for its files, in order) under three option sets, from
the 16-bit files, from their float copies, and through `--feed` from each,
and requires every way to print what jt9 prints from the 16-bit files, and
leave the same `decoded.txt`, byte for byte. On 64-bit ARM a few of jt9's outputs differ in known ways, FFTW's
NEON code rounding differently from its AVX code the stored outputs were
made with (the spread's debug spectrum, `fort.52`, in its last digits; one
weak FT8 decode in three short periods): each is listed with its exact ARM
output in `tests/data/jt9-arm64-accepted.txt`, and the check counts those
as accepted, not as failures; any other difference still fails it. The deep check builds WSJT-X's own programs and runs them beside
these on every option and recording:

```
sudo apt install gfortran
git clone https://github.com/wsjtx/wsjtx ~/repos/wsjtx && git -C ~/repos/wsjtx checkout 567ad29ce
./build.sh --deep-check
```

`scripts/wsprd-oracle.sh --build-only` builds WSJT-X's `wsprd` from its
source with WSJT-X's Release flags (`WSJTX=` another checkout), and
`scripts/wsprd-compare.sh` runs both side by side, each in fresh
directories, comparing everything they print and write byte for byte:

- every case in `tests/data/wsprd/cases`: the WSJT-X sample, the five 40 m
  periods, synthetic periods (SNR sweeps down to where Fano gives up and the
  OSD takes over, crowds past the 50-decode cap, drift, hashed and compound
  calls, WSPR-15), each under many option sets, seeded hash tables, runs in
  sequence in one data directory, read-only files (see its README);
- the command line: no arguments, bad and missing options, clustered flags
  and attached values, operands anywhere, `--`, missing and misnamed files,
  `-a` directories, odd values, `-m`, `-0`;
- every recording under each option alone and the common sets.

802 cases, all identical (only `wspr_timer.out`'s times, which are the
clock's, are compared by their shape).

For `jt9` (Docker needed): `scripts/fst4w-oracle.sh --sim-only` builds
WSJT-X's `jt9` with WSJT-X's own CMake in a container
(`docker/wsjtx-oracle.Dockerfile`) and makes the generated inputs
(`recordings/fst4w-sim`: FST4W of every period and type, hashed calls, weak
and faded signals, with WSPR beside them; and WSJT-X's own FST4W-1800 sample),
and `scripts/jt9-compare.sh` runs both, comparing stdout, stderr, the exit
status and every file written, byte for byte: every case in
`tests/data/fst4w/cases`, the command line (bad, missing and attached
values, long options, `--`, odd and broken WAV files, every option's
values), and every recording. For FT8 the same on the cases of
`tests/data/ft8` and every FT8 recording in `recordings/` (generated
periods: QSOs, crowded and weak bands, every message type, signals late
and at the band's edges; WSJT-X's own FT8 sample and a recorded period;
and the file-mode inputs: a-priori decoding on my and his call at every
QSO progress, contest modes, short and odd WAVs, sequences), at several
depths and options, and FT8's command line (every option's values, odd
and broken WAV files, calls and grids as Fortran reads them, several files
in a run). For FT4 likewise (`tests/data/ft4` and every FT4 recording in
`recordings/`: generated periods, WSJT-X's own FT4 sample and eight
recorded 40 m periods; a-priori decoding at every QSO progress, the
contest modes, short WAVs, sequences, the command line with `-5`). FT2's
reference is not WSJT-X's `jt9` but the hybrid `ft2-jt9` (improved
3.2.0's FT2 code built into the `567ad29ce` tree, [jt9: FT2](#jt9-ft2)),
which `--deep-check` does not build: the 19 `jt9 --ft2` cases of
`tests/data/ft2` (generated FT2 periods and twenty recorded 40 m
periods: sequences, a-priori decoding with calls, the defaults, short,
late and quiet periods, the band's edges) are checked against its outputs
stored there, and FT2's other cases (every generated period at each
depth, odd WAVs and cuts around the 34,560 samples read, the command line
with `--ft2`) run only with that `jt9` given as `scripts/jt9-compare.sh`'s
second program, and are skipped otherwise (230 of them). Over 650
cases run by `--deep-check`, every one identical on x86-64 Linux (25 FT8 cases
that need recorded periods not included here are skipped). Against the
hybrid itself all 468 of FT2's cases, with 239 recorded 40 m periods, were
identical. On
arm64 FFTW's NEON code rounds the spread's debug spectrum (`fort.52`)
differently in its last printed digit, and a few FT8 decodes near a
threshold (see [jt9: FT8](#jt9-ft8)); everything else is the same.

On another platform than x86-64 Linux, FFTW's code paths differ, so a .wav's
first transforms differ in their last bits; no WSPR or FST4W decode differed
in any case.

## What is here

| | |
|---|---|
| `src/wsprd.cpp`, `wsprd_fano.cpp`, `wsprd_unpack.cpp`, `wsprd_osd.cpp` | the decoder: WSJT-X's `wsprd.c`, `fano.c`, `wsprd_utils.c`/`wsprsim_utils.c` and `osdwspr.f90`/`indexx.f90` (each file's head says what changed and why) |
| `src/wsprd_files.cpp` | its command line, files and output: `wsprd.c`'s `main()` |
| `src/wsprd_main.cpp` | the program |
| `src/fft.cpp`, `fft_wisdom.cpp` | a small FFTW wrapper, for the spread |
| `src/fst4w*.cpp` | jt9's FST4W decoder: WSJT-X's `fst4_decode.f90` and what it calls (`get_fst4_bitmetrics`, `decode240_74`, `osd240_74`, `packjt77`'s unpacking, `gen_fst4wave`...) |
| `src/ft8*.cpp`, `packjt77.cpp` | jt9's FT8 decoder: WSJT-X's `ft8_decode.f90` and what it calls (`sync8`, `ft8b`, `decode174_91`, `osd174_91`, `subtractft8`, `ft8_a7`, `ft8_a8`, `packjt77`), and the maths library's functions it uses (`ft8_maths.cpp`) |
| `src/ft4*.cpp` | jt9's FT4 decoder: WSJT-X's `ft4_decode.f90` and what it calls (`getcandidates4`, `sync4d`, `get_ft4_bitmetrics`, `subtractft4`, `gen_ft4wave`), on FT8's LDPC decoder, `packjt77` and maths |
| `src/ft2*.cpp` | WSJT-X improved's FT2: its decoder's FT2 branch (the audio stretched 2x, FT2's sizes and threshold) and `ft2_decoded`'s output, on the FT4 decoder |
| `src/jt9_cli.cpp`, `fst4w_jt9.cpp`, `ft8_jt9.cpp`, `ft4_jt9.cpp`, `ft2_jt9.cpp`, `jt9_main.cpp` | jt9's command line, files and output: `jt9.f90`, `decoder.f90`, `options.f90` |
| `scripts/` | the comparisons with WSJT-X's `wsprd` and `jt9`, and the synthetic periods' generators |
| `tests/data/wsprd`, `tests/data/fst4w`, `tests/data/ft8`, `tests/data/ft4`, `tests/data/ft2` | WSJT-X's `wsprd`'s and `jt9`'s outputs (the hybrid's for FT2), case by case, and wsprd's synthetic inputs |
| `recordings/` | the six real WSPR periods, and FT8's, FT4's and FT2's inputs |

## Credits and licence

The original works are the WSJT-X project's, and the credit for the decoder
is theirs: <https://wsjt.sourceforge.io/wsjtx.html>. In particular
`wsprd.c`, Copyright 2001-2018 Joe Taylor, K1JT, with much of it by Steven
Franke, K9AN, 2014-2018; `fano.c`, Phil Karn, KA9Q, 1994; `osdwspr.f90`,
Steven Franke and Joe Taylor; FST4W's `fst4_decode.f90` and its routines,
Steven Franke, K9AN, and Joe Taylor, K1JT; FT8's `ft8_decode.f90`,
`sync8.f90`, `ft8b.f90`, `decode174_91.f90`, `osd174_91.f90` and
`subtractft8.f90`, and FT4's `ft4_decode.f90` and its routines, Joe
Taylor, K1JT, and Steven Franke, K9AN;
`packjt77.f90`, Joe Taylor, Steven Franke and the WSJT-X developers; FT2's
code, WSJT-X improved 3.2.0, Uwe Risse, DG2YCB, and the WSJT-X developers. The
maths functions in `src/ft8_maths.cpp` are the GNU C Library's (LGPL 2.1
or later).

This code is released under the same licence as WSJT-X, the GNU General
Public License, version 3 (`LICENSE`).
