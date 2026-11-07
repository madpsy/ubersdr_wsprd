# MSK144 test corpus

The reference outputs for `jt9 -k` (`--msk144`): WSJT-X 3.2.0-rc1, commit
567ad29ce, built in the `wsjtx-oracle` image by `scripts/msk144-oracle.sh`.
Every piece of the port checks itself against these files, or against a small
Fortran oracle built from the same tree.

## Regenerate

    git clone https://github.com/wsjtx/wsjtx ~/repos/wsjtx
    git -C ~/repos/wsjtx checkout 567ad29ce
    scripts/msk144-oracle.sh all            # builds, then makes the corpus
    scripts/msk144-oracle.sh --print-env    # paths of the binaries

`DOCKER=podman` picks the container runner; without it the script uses docker
if present, else podman. The script builds jt9 and `msk144sim` three times:

| build     | what                                                              |
|-----------|-------------------------------------------------------------------|
| `release` | WSJT-X's own Release flags: the reference, `$MSK144_REF`          |
| `o2`      | `-O2 -ffp-contract=off`; must agree with `release` on every run   |
| `patched` | `release` with `bshmsg` and `bswl` read from `MSK144_SHMSG=1` and `MSK144_SWL=1` in the environment (`$MSK144_REF_PATCHED`) |

WSJT-X's `decode_msk144.f90` sets both switches false and offers no way to
turn them on. The sed patch in the script (recorded in `BUILD.txt`) replaces
the two initial values by reads of those variables. With both set, both are
forced true. With one set, the build gives what `--msk-shmsg` or `--msk-swl`
alone must give. With neither, it decodes as `release` (the script checks this
on the two WSJT-X samples).

`scripts/msk144-oracle.sh fortran FILE.f90` compiles a Fortran oracle against
the reference's `libwsjt_fort.a` and prints the binary's path. `exec CMD`
runs a command in the image with the work dir and repository mounted at their
own paths.

## Layout

    tests/data/msk144/cases/CASE/jt9-VARIANT/{args.txt,stdout.txt,status.txt[,stderr.txt],decoded.txt}
    tests/data/msk144/inputs/SHA256SUMS      sha256sum -c, run from the repository root
    recordings/msk144-samples/               WSJT-X's samples
    recordings/msk144-sim/                   the generated grid (+ MANIFEST: case, T/R)
    recordings/msk144-real/                  a real slot
    recordings/msk144-uacnj/                 20 real slots from UACNJ, each with its JSON sidecar

`args.txt` holds `jt9` and its arguments, paths from the repository root.
`decoded.txt` is the empty file jt9 opens: MSK144 writes nothing to it.
`scripts/jt9-compare.sh` names a run `msk144-CASE/VARIANT`.

## Variants

| variant    | arguments (T/R of the case in place of 15)                         | from    |
|------------|--------------------------------------------------------------------|---------|
| `d1`       | `-a . -k -p 15 -f 1500 -F 200 -d 1`                                | release |
| `d3`       | `-a . -k -p 15 -f 1500 -F 200 -d 3`                                | release |
| `defaults` | `-a . -k -p 15` (jt9's own `-f`, `-F`, `-d`)                       | release |
| `ms`       | `-Y --msk144 -p 15 -f 1500 -F 200 -d 3 -a .`: the line meteor-scatter runs (`build_jt9_cmd` in its `core/decoder.py`), with `-a .` for its work dir | release |
| `bare`     | `-a . -k` alone (jt9's 60 s period); the four samples only         | release |
| `shmsg`    | as `d3`, plus `--msk-shmsg`                                        | patched, `MSK144_SHMSG=1` |
| `swl`      | as `d3`, plus `--msk-swl`                                          | patched, `MSK144_SWL=1` |

`-Y` is jt9's FST4W switch; `--msk144` after it wins, as the last mode does.
`shmsg` and `swl` exist for the samples, the real slot, the SNR sweep, the
message types, the pings, the MSK40 cases, the two `msk144sim` cases, silence,
noise and the clipped file. Our jt9 takes `--msk-shmsg` and `--msk-swl` (long
options only, off by default; Task 7): the reference has no such option to run
beside, so `jt9-compare.sh` and `jt9-wine-check.sh` check those two variants
against the stored outputs, and run them against any jt9 whose usage lists
both options.

## Task 6's cases (`scripts/msk144-oracle.sh focus`)

Behaviours that each look wrong and must stay as the reference has them. Made
one reference process per file with `OMP_DYNAMIC` and `OMP_THREAD_LIMIT` unset:
after its first file the reference's `decoder.f90` calls `omp_set_dynamic(.true.)`,
so a second file's answer in the same process depends on the load. A
multi-file case's stored output is the single-file runs end to end, never a
multi-file reference run. `jt9-compare.sh` checks those cases (more than one
`.wav` in `args.txt`, or a `stdin.txt`) against the stored files alone, and
runs `--feed` with `stdin.txt` (paths from the repository's root) as stdin.

| case                  | variants and what they show                                                                 |
|-----------------------|---------------------------------------------------------------------------------------------|
| `repeat-across-files` | `files`: one command line, `feed`: `--feed`. Six files: a decode in the first block twice (its SNR carries any state from the file before), the 181211_120500 sample twice, a file with a hashed call that the earlier message resolves (`<W9XYZ> PJ4/K1ABC RRR`), a file with the hashed message alone (`<...>`), both again. `feed` prints `<FileDone> PATH N` after each file, N its decode lines |
| `lastmode-Y`          | `Ymsk144` (`-Y --msk144`: MSK144), `msk144Y` (`--msk144 -Y`: FST4W)                        |
| `lastmode-8`          | `8msk144` (MSK144), `msk1448` (`--msk144 -8`: FT8)                                         |
| `ntol-wide`           | `-F 300`, `1000`, `1001`, `5000` on a -200 Hz offset (jt9.f90 limits `-F` to 1000); the `f2000` runs of `freq-edge` cover `-F` above 1000 from the other side |
| `freq-edge`           | `-f 0`, `200`, `800`, `1801`, `1802`, `1803`, `2000`, `6000`, `9000` with `-F` 200 to 5000, on an 800 Hz carrier. `-f 1801`, `1802` and `1803` (each with `-F 5000`) sit 1001, 1002 and 1003 Hz from the carrier, where a limit of 1001 would answer differently; `-f 2000` (with `-F 5000` and `-F 1000`) is 1200 Hz away, past what the limit of 1000 reaches |
| `offgrid-sm06`, `offgrid-sm04` | `d1`, `d3`: a 1 s ping at -6 dB and at -4 dB on a 1537 Hz carrier, which lies on none of the averaging searches' frequency grids: `-d 1` (short pings only) misses it and `-d 3` (the averages) finds it, so msk144sync's step of 10/navg Hz shows |
| `swl-state` (Task 7) | `files`, `feed`: `--msk-swl` on three files: three pings of the MSK40 frame `<N0ABC W1AW> RRR`, a standard message from `N0ABC` to `W1AW`, the MSK40 file again. The reference (one process a file, patched build with `MSK144_SWL=1`) prints `<1591> RRR` for both MSK40 files: its table of recent-call hashes is empty at the start of each. A port that kept the table of hashes (but not the calls, which each file starts without) from the middle file into the last would print `< > RRR` there |
| `swl-pair` (Task 7) | `swl`: `--msk-swl` on one file: the standard message at 5 s, which gives the decoder `N0ABC` and `W1AW`, and the MSK40 pings of `<N0ABC W1AW> RRR` around it. The ping before it prints as `<1591> RRR`, the one after it by the pair's calls: the table of recent-call hashes is made when the standard message decodes |

Their inputs are listed in `recordings/msk144-sim/MANIFEST.focus` (made by
`scripts/msk144-sim-grid.sh`, which gives `corpus` no variants on them).

## The final fix wave's cases (`scripts/msk144-oracle.sh extra`)

Made one reference process per file, `OMP_DYNAMIC` and `OMP_THREAD_LIMIT` unset.

| case | variants and what they show |
|------|-----------------------------|
| `odd-trunc`  | `d1`, `d3`, `defaults`, `ms`: a WAV whose header says 15 s but whose file stops at 9 s (cut from `181211_120500`). The reference and the port print `WAV chunk extends beyond end of file` and exit 2. 3.0.2 decoded the 9 s that were there and exited 0 (`AB-3.0.2.md`) |
| `odd-stream` | the same four: a WAV with the streaming writers' placeholder size 0xFFFFFFFF in both headers, the full 15 s of data. The reference and the port print `WAV chunk has an invalid size` and exit 2; 3.0.2 decoded the file and exited 0 |
| `real-uacnj-*` | `ms` only: 20 real 30 s slots from UACNJ's meteor-scatter archive (50.260 MHz, slots of 2026-11-04 to 2026-11-06), run as UACNJ's station runs them, `-Y --msk144 -p 30 -f 1500 -F 200 -d 3 -a .` (its `tr_period_sec` is 30) |

The two odd WAVs are `@I/` inputs, made by `scripts/jt9-compare.sh`'s generator
like the others. The `real-uacnj-*` inputs are in `recordings/msk144-uacnj/`,
12 kHz s16 mono, each beside the JSON sidecar meteor-scatter wrote (its decodes,
the frontend status). Station data of this kind is public.

## The float cases (ours: no reference)

WSJT-X's `jt9` refuses a float WAV for MSK144, so these cases have no
reference output. Their stored files come from our port (scripts/msk144-float-inputs.py
makes the two inputs in `recordings/msk144-float/` from `181211_120500`), and
`jt9-compare.sh` and `jt9-wine-check.sh` check them against those files alone.
Run `scripts/float-input-check.sh` for the other direction: every 16-bit case run
again on its float copy must leave the same files.

| case          | variants | what it shows |
|---------------|----------|---------------|
| `float-nan`   | `ms`, `swl` | a float file with NaN, +inf and -inf written over it: 10 NaN, one +inf and one -inf in the first 0.6 s (the first non-silent block, where a NaN would become the noise level for the whole file and print every SNR as 0), then more away from the decode. The port exits 0 at the normal speed, prints the decode at 8.7 s with its true SNR (8) and `<DecodeFinished>`. A block with a non-finite sample counts as silent. It writes `decoded.txt` (empty), because a NaN in the first 180000 samples makes the rms NaN and `!(rms < 0.5)` holds, as in `decoder.f90`'s semantics |
| `float-level` | `ms` | the same samples times 2^-20 (an rms far below one 16-bit count). It decodes the same line as `sample-181211_120500`. It leaves no `decoded.txt`: that file's rms test (`decoder.f90`) (`decoder.f90:104-107`: `!(rms < 0.5)` over the first 180000 samples) runs on the scaled values, and finds a quiet file |

## Inputs

**WSJT-X's samples.** `181211_120500.wav` and `181211_120800.wav` come from
`samples/MSK144/` in the 567ad29ce tree. `160915_113100.wav` and
`160915_113230.wav` are not in that tree or in 3.0.2's. They came from an
earlier session's scratch directory, where someone had fetched them from
WSJT-X's sample download (https://sourceforge.net/projects/wsjt/files/samples/MSK144/);
the session kept no record of the exact URL. Both are 12 kHz mono 16-bit,
15 s. Neither decodes with the reference (the 2016 recordings predate the
128/90 LDPC format of 2019), so they test the no-decode path on real audio.

**Generated grid.** `scripts/msk144-sim-grid.sh` drives two generators:
WSJT-X's `msk144sim` (two cases, `stock-*`, its usage examples) and
`scripts/msk144-sim.f90` (`msk144-sim`), which takes the choices `msk144sim`
fixes: pings per slot and their times, seeded noise, MSK40 frames. Both run
against WSJT-X's own `genmsk_128_90` (which calls `genmsk40`), `makepings`'
ping shape and `gran`'s noise scaling. File *n* has seed 1000+*n*. Axes, each
swept with the others at their base (`K1ABC W9XYZ EN37`, 1500 Hz, one ping
0.2 s wide at 5 s, T/R 15 s):

- `snr-*`: -8 to +10 dB in 2 dB steps
- `freq-*`: offsets 0, +-100, +-200 Hz, at +4 and -2 dB
- `width-*`: 0.05, 0.2, 1, 5 s and `slot` (the whole slot at peak amplitude), at +4 and 0 dB
- `pings1-*`, `pings3-*`: one and three pings a slot, at 0 and +6 dB
- `tr30-*`: T/R 30 s, two pings, at -2, +4, +10 dB; `tr30-wide-s04` has width 2.5 s
- `msg-*`: CQ, report, RRR, prefix, suffix, free text, at +6 dB
- `msk40-*`: four MSK40 frames (`<K1ABC W9XYZ> +03` and three more), for Task 7's short-message flag
- `odd-*`: `silence` (zeros), `noise` (noise alone), `tr30-run-p15` (30 s file with `-p 15`), `short7` (7 s file, `-p 15`), `clipped` (`181211_120800` at 40x gain, clipped); and `stereo`, `rate11025`, `bits8`, `header`, `tiny`, made by `scripts/jt9-compare.sh`'s own generator (`@I/` in `args.txt`)

**Real.** `msk144_261023_105030.wav`: 30 s of 6 m (50.260 MHz) audio from
AC0G/B4's meteor-scatter recorder, slot start 2026-10-23 10:50:30 UTC. Its
sidecar JSON records one decode, `C67MIP/P BI7CEG/P NN73` at -8 dB, and a
frontend status of "partial" (radiod status packets missing for part of the
slot). Two cases run it: `real-261023-p30`
(T/R 30 s) and `real-261023-p15` (its first 15 s). Copied from
`proofs/agu-2026/spectrum/data/`.

**UACNJ's slots.** `recordings/msk144-uacnj/*.wav` and their `.json`
sidecars: 20 slots of meteor-scatter's archive at UACNJ, 6 m, 30 s each. Cases
`real-uacnj-YYMMDD_HHMMSS`. Each runs at `-p 30`, as UACNJ's station runs
meteor-scatter (`tr_period_sec = 30` in its configuration; meteor-scatter's own
code default is 15), so the reference reads the whole slot. They are a
sample, not the hour the design (section 7) asks for.

## What the reference does here

- 97 cases, 395 runs (the float cases, ours, add 2 cases and 3 runs, 398 in all): 338 from the `corpus` step (66 cases), 29 from `focus` and 28 from `extra`.
  `release` and `o2` agree on all 268 runs of the plain variants.
- On the SNR sweep, depth 1 decodes down to -2 dB and depth 3 down to -4 dB. At -6 dB and below nothing decodes.
- `-F 200` reaches all four +-200 Hz offsets at depth 3; at depth 1 three of the four fail. `defaults` (`-F 20`) loses every offset but 0.
- MSK40 frames decode only under `shmsg` or `swl`, and not all of them: `msk40-rpt03` decodes under neither, `msk40-r03` only under `swl`. The port must reproduce that.
- The clipped file still decodes (two lines, as the unclipped sample).
- The stereo, 11025 Hz and 8-bit files end with status 2 and a line on stderr; the header-only and tiny files end 0.
