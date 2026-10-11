# wsprd oracle

What WSJT-X's own `wsprd` prints and writes, case by case, for the in-process
WSPR decoder to match byte for byte. Made by `scripts/wsprd-oracle.sh` from
WSJT-X at commit 567ad29ce (`lib/wsprd`, `lib/indexx.f90`, `lib/lookup3.c`,
`lib/init_random_seed.c`: CMake/Sources.cmake's `wsprd_CSRCS` plus the two
Fortran files), built by hand with gcc/gfortran 13.3 against the system FFTW
3.3.10 (x86_64, AVX2), twice:

- with WSJT-X's Release flags (C: `-Wall -Wextra -fopenmp -pthread -O3
  -DNDEBUG -fdata-sections -ffunction-sections`; Fortran: `-O3 -DNDEBUG
  -funroll-loops -fno-f2c -fno-second-underscore ...`; link `--gc-sections`),
- with `-O2 -ffp-contract=off`.

The two agree on every case, to the byte (stdout, ALL_WSPR.TXT,
wspr_spots.txt, hashtable.txt and the `-c` file); the files here are the
Release build's.

## Inputs

`inputs/` (in git): four synthetic WSPR-2 `.c2` files, 360026 bytes each,
and two `.wav` files cut short from `sim-wav2`.
`recordings/wsprd-sim/` (not in git; the script makes it): those again, the
other synthetic files and WSJT-X's sample `samples/WSPR/150426_0918.wav`.
`SHA256SUMS` in each. Cases on a missing input are to be skipped.

| file | what |
| --- | --- |
| `sim-types-261010_1200.c2` | 12 signals, -15 to -24 dB: type 1, type 2 (`PJ4/K1ABC`, `G4XYZ/P`, `K9AN/12`, `VE3/W1AW`), type 3 (`<PJ4/K1ABC> FK52UD`, a hash nobody sent: `<...>`, a 4-character grid), powers 0 and 60 |
| `sim-sweep-261010_1202.c2` | 12 calls, -22 to -35 dB: Fano's threshold, pass 3's block demodulation, OSD when the calls are known |
| `cut-odd-261010_1218.wav` | `sim-wav2`'s first 10000 samples and the low byte of the next (an odd data length: glibc's fread keeps that byte in wsprd's buffer) |
| `cut-tiny-261010_1218.wav` | `sim-wav2`'s first 5 samples (fewer than 22: the header's shorts stay in wsprd's buffer after them, as samples) |
| `sim-crowd-261010_1210.c2` | 27 signals 2-7 Hz apart, -12 to -28 dB, drift -1/0/1; one call twice 3 Hz apart (the dupe check) |
| `sim-drift-261010_1206.c2` | 10 signals, drift -4 to +4 Hz, start -1 to 3.5 s |
| `sim-sweep2-261010_1204.c2` | the sweep calls, -26.5 to -32 dB in 0.5 dB steps, three start times |
| `sim-wide-261010_1212.c2` | signals at +/-107 to 149 Hz: outside +/-110 only `-w` decodes |
| `sim-edge-261010_1214.c2` | starts at -3, -1.5, 4.5, 6.5 and 8 s (partial symbols), +5 and -3 dB |
| `sim-noise-261010_1216.c2` | noise only |
| `sim-wspr15-261010_1215.c2` | a .c2 whose header says WSPR-15 (its SNR scale and f/8, DT*8 printing) |
| `sim-strong-261010_1158.c2` | 21 calls at -10 dB: its `hashtable.txt` seeds the `-seeded` cases |
| `sim-wav2-261010_1218.wav` | 12 kHz WSPR-2 audio, 6 signals (readwavfile's FFTs) |
| `sim-wav15-261010_1230.wav` | 12 kHz WSPR-15 audio, 15 min, 3 signals (`-m`) |
| `150426_0918.wav` | WSJT-X's sample, dial 14.0956 |
| `recordings/40m-wspr-20261010-1914/261010_1914.wav` .. `_1922.wav` | five consecutive periods off the air on 40 m, dial 7.0386 (cut from an UberSDR IQ recording by `scripts/wspr-cut.py`): `<...>`, `<IU0JJD>`, `SV1/DF3LZ` |

The synthetic files come from `scripts/wsprd-sim.c` (wsprsim with a seeded
generator, several signals, drift, start time, WSPR-15 and .wav output). The
noise is seeded (xorshift64*, not `rand()`), so another machine gives the
same bytes as long as its libm's `cos`, `sin`, `log`, `sqrt` and `pow` agree;
the checksums say if not.

## Cases

`cases/<input>/<options>/`, 222 runs in 200 cases. Options, each on every
input (`-m` added for `sim-wav15`):

| name | wsprd options |
| --- | --- |
| `skim` | `-c -C 500 -o 4 -d` (deep) |
| `skim-nodeep` | `-c -C 500 -o 4` |
| `default` | none |
| `wide` | `-w` |
| `quick` | `-q` |
| `single` | `-s` |
| `osd0` | `-o 0` |
| `nohash` | `-H -o 4` |
| `noblock` | `-B` |
| `stack` | `-J` |
| `bias` | `-z 0.6` |
| `dialerr` | `-e 2.5 -d` (wsprd adds the error to fmin/fmax again on each pass) |
| `verbose` | `-v -d` |

all with `-f <dial>` (ignored for a .c2, which carries its own). Then
`sim-sweep/{skim,osd0,nohash}-seeded`, `sim-sweep2/{skim,stack}-seeded`
(`-J -o 4 -d`) and `sim-types/skim-seeded`, which start with
`sim-strong/skim`'s hashtable (OSD only accepts a call already in the
hashtable), and `sequence/*`: several runs in one data dir, as a station's
data dir carries hashtable.txt, ALL_WSPR.TXT and wspr_wisdom.dat from period
to period. On the 40 m recording: `40m-HHMM/skim-o3` (the deep options with
OSD depth 3, `-c -C 500 -o 3 -d`) for each period,
`sequence/40m-{skim-o3,default,stack-o2}` (all five periods in one data dir,
with those options, none, and `-J -o 2 -d -w`), `sequence/40m-skim-o4` (all
five with the deep options, `-C 500 -o 4 -d`, no `-c`),
`40m-1916/single` (`-s`) and
`40m-1920/nohash-o2` (`-H -o 2 -q`). `cut-odd/skim-o3` and `cut-tiny/osd0`
on the cut `.wav` files. `names/short`: `sim-types` as `x.c2`, `a.c2`,
`1200.c2` and `x1200.c2`, names too short for the date (11 characters before
the suffix) or the time (4), where wsprd reads on into the arguments before
(argv's strings lie end to end). `readonly/{hashtable,all-wspr,spots,timer}`:
that file read-only, so wsprd cannot open it to write and crashes (status
139) where it first uses it, losing what stdio had not yet written;
`readonly/{all-wspr,spots}-none` with no decodes (`cut-tiny`), where the
crash comes at closing the files after the timings, and
`readonly/all-wspr-none-hash` with `hashtable-odd.txt` (an odd line) as the
seed, which wsprd then never rewrites. `sim-sweep/stack-z{20,18}`: the stack
decoder with a bias so large its bucket index wraps, a crash in the first
pass and after it. New
cases keep OSD at depth 3 or less: deeper is too slow to run in bulk.

Each case directory:

- `args.txt`: one line per run, wsprd's argv after the program name, without
  `-a`. Inputs are relative to the repository root.
- `hashtable.in` (seeded cases only): copy it to `hashtable.txt` in the data
  dir before the first run.
- `files.txt` (some cases): files to put in the data dir before the first
  run, one a line: `SRC NAME` copies SRC (relative to the repository root;
  `-`: an empty file) to NAME, `SRC NAME ro` makes it read-only too. Where
  read-only does not stop writing (root), the case is skipped.
- For run N (from 1): `stdout.N.txt`; `wspr_spots.N.txt` (wsprd rewrites it
  each run); `c2.N.sha256`, the sha256 of the `000000_0001.c2` that `-c`
  writes to the current directory (not the data dir) after pass 1.
- `status.txt`: each run's exit status, one per line (0, or 139 for a
  crash).
- After the last run: `ALL_WSPR.TXT` (appended to by each run) and
  `hashtable.txt` (absent with `-H`, which neither reads nor writes it; a
  seeded `-H` case leaves the seed as it was).

Each case runs in a fresh, empty data dir (`-a .`) that is also the current
directory. `wspr_timer.out` (timings) and `wspr_wisdom.dat` (FFTW's) are not
kept. No run writes to stderr.

## FFTW

wsprd plans with `FFTW_ESTIMATE` and first imports `wspr_wisdom.dat` from
its data dir. ESTIMATE uses wisdom when there is some, so wisdom changes the
plan and with it the bits: with MEASURE wisdom for wsprd's own transforms
(`fftwf-wisdom -m cof512 cob46080 rof1474560`) the .wav cases' `-c` files
change (their decodes did not, this time); the .c2 cases do not change. Other
wisdom in the process (in-place forward transforms of 1024-8192 points)
changes nothing. The plans differ by CPU and FFTW build too (AVX
codelets here), so the .wav cases are only exact against an FFTW like this
one. `scripts/wsprd-oracle.sh --wisdom FILE` repeats every case with FILE
as the wisdom and lists what changes. `wsprd-measure-wisdom.txt` is that
MEASURE wisdom: the decoder's tests import it into the process's planner
and checks the decoder's results do not change (the decoder plans with
only its own wisdom, wspr_wisdom.dat's, in the planner).

`osd_vectors.txt`: soft-symbol vectors and the reference Fortran osdwspr's
answers, for the OSD's tests.

## Regenerating

    scripts/wsprd-oracle.sh             # WSJTX=~/repos/wsjtx by default

rebuilds both references, the inputs and every case, and stops if the two
builds disagree.
