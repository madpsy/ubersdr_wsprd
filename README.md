# wsprd

WSJT-X's WSPR decoder, `wsprd`, rewritten in C++ to take a fraction of the
CPU.

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

It decodes what WSJT-X's `wsprd` decodes, WSPR-2 and WSPR-15, and nothing
else: **not FST4W**, which is decoded by WSJT-X's `jt9` (`jt9 --fst4w`), a
separate program this does not replace. Keep `jt9` installed for FST4W.

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

On Linux it needs FFTW's single-precision library:
`sudo apt install libfftw3-single3` on Debian, Ubuntu or Raspberry Pi OS.

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

About 21 MB while a WSPR period decodes (140 MB for WSPR-15), freed after.
To time it yourself, on any of `recordings/`:

```
time ./build/wsprd -C 500 -o 4 -d -f 7.0386 recordings/40m-wspr-20261010-1914/261010_1914.wav
```

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
in a minute or two): `build/wsprd`, and stripped, in `build/dist`. Or by hand:

```
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
```

On an x86-64 CPU without AVX2, build for the compiler's baseline instead
(the same decodes, slower):

```
cmake -B build -DCMAKE_BUILD_TYPE=Release -DWSPRD_X86_ARCH=
cmake --build build
```

A 32-bit ARM system (`armv7l`) builds the plain code by itself.

**Windows**, cross-compiled on Linux with Docker (mingw-w64, FFTW built in,
so `wsprd.exe` is one file that needs nothing beside it):

```
./build.sh --windows
```

It lands in `build/dist/wsprd-windows-x86_64.exe`, and when `wine` is
installed it is checked to decode as the Linux build does.

Build with the settings in `CMakeLists.txt`: the decoder's files are
compiled with no fused multiply-adds (`-ffp-contract=off`), as WSJT-X's own
build rounds them, on every platform; that is part of what keeps the decodes
the reference's.

## Checking it against WSJT-X's own

```
sudo apt install gfortran
git clone https://github.com/wsjtx/wsjtx ~/repos/wsjtx && git -C ~/repos/wsjtx checkout 567ad29ce
./build.sh --check
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

On another platform than x86-64 Linux, FFTW's code paths differ, so a .wav's
first transforms differ in their last bits; no decode differed in any case.

## What is here

| | |
|---|---|
| `src/wsprd.cpp`, `wsprd_fano.cpp`, `wsprd_unpack.cpp`, `wsprd_osd.cpp` | the decoder: WSJT-X's `wsprd.c`, `fano.c`, `wsprd_utils.c`/`wsprsim_utils.c` and `osdwspr.f90`/`indexx.f90` (each file's head says what changed and why) |
| `src/wsprd_files.cpp` | its command line, files and output: `wsprd.c`'s `main()` |
| `src/wsprd_main.cpp` | the program |
| `src/fft.cpp`, `fft_wisdom.cpp` | a small FFTW wrapper, for the spread |
| `scripts/` | the comparison with WSJT-X's `wsprd`, and the synthetic periods' generator |
| `tests/data/wsprd` | WSJT-X's `wsprd`'s outputs, case by case, and the synthetic inputs |
| `recordings/` | the six real periods |

## Credits and licence

The original works are the WSJT-X project's, and the credit for the decoder
is theirs: <https://wsjt.sourceforge.io/wsjtx.html>. In particular
`wsprd.c`, Copyright 2001-2018 Joe Taylor, K1JT, with much of it by Steven
Franke, K9AN, 2014-2018; `fano.c`, Phil Karn, KA9Q, 1994; `osdwspr.f90`,
Steven Franke and Joe Taylor.

This code is released under the same licence as WSJT-X, the GNU General
Public License, version 3 (`LICENSE`).
