#!/usr/bin/env bash
# Run our jt9 (NEW: src/jt9_main.cpp) and WSJT-X's own jt9 (STOCK:
# scripts/fst4w-oracle.sh or scripts/ft8-oracle.sh builds it, at 567ad29ce)
# side by side, each in fresh directories of its own, and compare everything
# they print and write: stdout, stderr, the exit status, and the files left
# in the run's directories, subdirectories too (fst4_decodes.dat,
# decoded.txt, fst4w_calls.txt, fort.52; not timer.out or jt9_wisdom.dat,
# the clock's and FFTW's, which NEW does not write). FST4W's cases and
# FT8's (jt9 -8: the oracle's jt9 cases, every file-mode fixture case,
# every FT8 recording and generated period, odd WAVs and the command line).
#
#   scripts/jt9-compare.sh NEW_JT9 STOCK_JT9 [--docker IMAGE]
#
#   --docker IMAGE                STOCK runs in IMAGE (the oracle's jt9 needs
#                                 its libraries: wsjtx-oracle), the repository
#                                 and work dir mounted at their own paths, its
#                                 stdout and stderr files there as NEW's are
#                                 (gfortran orders stderr by what it is
#                                 connected to: jt9_cli.h, stderr_as_file).
#                                 Where the host's glibc is not the build's
#                                 (arecibo), run the whole script in the image
#                                 instead, so that NEW has its own libm too.
#   JT9_COMPARE_WORK=/tmp/...     where to run (default: a temp dir, kept)
#   JT9_COMPARE_JOBS=n            runs at once (default: nproc)
#   JT9_COMPARE_ONLY=regex        only the cases whose names match
#   JT9_COMPARE_IGNORE="fort.52"  files not compared (arm64: fort.52, the
#                                 spread's debug spectrum, differs in its last
#                                 printed digit, FFTW's NEON rounding)
#
# NEW runs with "--spread exact" first (WSJT-X's spread, so every byte is
# its). Allowed apart:
#   usage   (-h, no file, an unknown option): each prints its own option
#           list; both must say " Usage:" and end 0, stderr the same
#   version (-v): each its own line
#   crash   (an FST4W period fst4_decode has no sizes for): STOCK dies of
#           SIGFPE with a backtrace, NEW prints gfortran's line alone: the
#           status and that line compared, and everything written before it
#   ioerror (an empty -c, -x, -G or -g: Fortran's list read ends): STOCK
#           says "At line N of file DIR/jt9.f90", the runtime error and a
#           backtrace, NEW the two lines without DIR: compared but for DIR,
#           and the status
# Prints PASS/FAIL/SKIP a case, a diff under each FAIL; exits 1 on any FAIL.
set -uo pipefail

usage() { echo "usage: $0 NEW_JT9 STOCK_JT9 [--docker IMAGE]" >&2; exit 2; }
(( $# >= 2 )) || usage
new=$(realpath "$1") stock=$(realpath "$2"); shift 2
image=
while (( $# )); do
    case $1 in
        --docker) (( $# >= 2 )) || usage; image=$2; shift ;;
        *) usage ;;
    esac
    shift
done
repo=$(cd "$(dirname "$0")/.." && pwd)
work=$(realpath -m "${JT9_COMPARE_WORK:-$(mktemp -d "${TMPDIR:-/tmp}/jt9-compare.XXXXXX")}")
jobs_max=${JT9_COMPARE_JOBS:-$(nproc)}
only=${JT9_COMPARE_ONLY:-}
ignore=" ${JT9_COMPARE_IGNORE:-} "
rm -rf "${work:?}/new" "${work:?}/stock" "${work:?}/results"
mkdir -p "$work/new" "$work/stock" "$work/results" "$work/inputs"
echo "jt9-compare: work dir $work" >&2

# ---------------------------------------------------------------- inputs
# Odd WAVs for the file checks: 11025 Hz, stereo, 8-bit, a header alone, not
# RIFF, a short period, a tiny one; and for FT8's: a LIST chunk before the
# data, an odd data size, data before fmt, data said longer than the file, a
# short fmt chunk, 30 s of a tone (more than jt9 -8 reads), an FT8 period cut
# to 13 s.
JT9_COMPARE_REPO=$repo python3 -I - "$work/inputs" <<'EOF'
import struct, sys, math
out = sys.argv[1]
def wav(name, rate=12000, ch=1, bits=16, n=12000 * 120, data=None):
    if data is None:
        data = b"".join(struct.pack("<h", int(3000 * math.sin(i * 0.3))) for i in range(n)) if bits == 16 else bytes(n)
    hdr = b"RIFF" + struct.pack("<I", 36 + len(data)) + b"WAVEfmt " + struct.pack("<IHHIIHH", 16, 1, ch, rate, rate * ch * bits // 8, ch * bits // 8, bits)
    open(f"{out}/{name}", "wb").write(hdr + b"data" + struct.pack("<I", len(data)) + data)
wav("rate11025-261012_1200.wav", rate=11025, n=11025 * 2)
wav("stereo-261012_1200.wav", ch=2, n=12000 * 4)
wav("bits8-261012_1200.wav", bits=8, n=12000 * 2)
wav("header-261012_1200.wav", data=b"")
wav("short-261012_1200.wav", n=12000 * 20)
wav("tiny-261012_1200.wav", n=1000)
open(f"{out}/text-261012_1200.wav", "w").write("not a wav\n")
def raw(name, chunks):
    body = b"WAVE" + b"".join(cid + struct.pack("<I", len(d) if n is None else n) + d + (b"\0" if len(d) % 2 else b"")
                              for cid, d, n in chunks)
    open(f"{out}/{name}", "wb").write(b"RIFF" + struct.pack("<I", len(body)) + body)
fmt = struct.pack("<HHIIHH", 1, 1, 12000, 24000, 2, 16)
tone = b"".join(struct.pack("<h", int(2000 * math.sin(i * 2 * math.pi * 1500 / 12000))) for i in range(12000 * 30))
raw("list-261012_193500.wav", [(b"fmt ", fmt, None), (b"LIST", b"INFOISFT\x05\0\0\0abcd\0", None), (b"data", tone[:360000], None)])
raw("odd-261012_193500.wav", [(b"fmt ", fmt, None), (b"data", tone[:360001], None)])
raw("data-first-261012_193500.wav", [(b"data", tone[:2000], None), (b"fmt ", fmt, None)])
raw("long-data-261012_193500.wav", [(b"fmt ", fmt, None), (b"data", tone[:2000], 400000)])
raw("short-fmt-261012_193500.wav", [(b"fmt ", fmt[:14], None), (b"data", tone[:2000], None)])
wav("tone30-261012_193500.wav", n=0, data=tone)
import os
src = os.environ.get("JT9_COMPARE_REPO", ".") + "/recordings/ft8-live/40m-20261012-1934/261012_193515.wav"
if os.path.exists(src):
    b = open(src, "rb").read()
    i = b.index(b"data")
    wav("cut13-261012_193515.wav", data=b[i + 8:i + 8 + 2 * 13 * 12000])
open(f"{out}/notes.txt", "w").write("not a wav\n")
EOF

# ---------------------------------------------------------------- cases
# A case: name|kind|seed|args. kind: same, usage or crash (above). seed: an
# fst4w_calls.txt put in the run's dir first (repository-relative) or -.
# args split on blanks; recordings/, tests/ and @I/ (the odd inputs) made
# absolute; '' an empty argument, @BLANK one blank. Each runs in its own fresh dir, the cwd.
cases=()
add() { cases+=("$1"); }

# 1. the oracle's cases, as scripts/fst4w-oracle.sh ran them
casesdir=$repo/tests/data/fst4w/cases
while IFS= read -r a; do
    d=$(dirname "$a"); name=${d#"$casesdir"/}
    seed=-; [[ -f $d/fst4w_calls.in ]] && seed=${d#"$repo"/}/fst4w_calls.in
    bare=; [[ -f $d/flags.txt ]] && grep -q bare "$d/flags.txt" && bare=bare
    add "oracle/$name|same${bare:+-bare}|$seed|-a . $(cat "$a")"
done < <(find "$casesdir" -name args.txt | sort)

S=recordings/fst4w-sim
mix=$S/f2-mix-261012_1200.wav
sil=$S/f2-silence-261012_1252.wav
h1=$S/f2-hash-261012_1300.wav h2=$S/f2-hash-261012_1302.wav h3=$S/f2-hash-261012_1304.wav
k="-p 120 --fst4w -f 1500 -F 100"

# 2. the command line
add "cli/no-args|usage|-|"
add "cli/h|usage|-|-h"
add "cli/help-long|usage|-|--help"
add "cli/no-file|usage|-|$k -d 3"
add "cli/unknown|usage|-|-z $mix"
add "cli/unknown-long|usage|-|--nothing $mix"
add "cli/attached|usage|-|-p120 --fst4w $mix"
add "cli/cluster|usage|-|-Wq -p 120 $mix"
add "cli/missing-p|usage|-|--fst4w -p"
add "cli/missing-p-dash|usage|-|-p -d 3 $mix"
add "cli/missing-long|usage|-|--fst4w --tr-period"
add "cli/dash-alone|usage|-|- $mix"
add "cli/v|version|-|-v"
add "cli/v-after|version|-|-p 120 -v $mix"
add "cli/long-value-quirk|same|-|--tr-period 120 $mix $mix"
add "cli/long-fst4w|same|-|-p 120 --fst4w -d 1 $mix"
add "cli/W|same|-|-p 120 -W -d 1 $mix"
add "cli/dashdash|same|-|$k -d 1 -- $mix"
add "cli/option-after-file|same|-|$k $mix -d 3"
add "cli/p-garbage|same|-|-p x --fst4w $mix"
add "cli/p-exp|same|-|-p 1.2e2 --fst4w -d 1 $mix"
add "cli/p-plus|same|-|-p +120 --fst4w -d 1 $mix"
add "cli/p-frac|same|-|-p 120.5 --fst4w -d 1 $mix"
add "cli/p-zero|same|-|-p 0 --fst4w $mix"
add "cli/p-high|same|-|-p 1800.5 --fst4w $mix"
add "cli/p-1e|same|-|-p 1e --fst4w $mix"
add "cli/d-garbage|same|-|-d 3.0 $k $mix"
add "cli/f-garbage|same|-|-f + $k $mix"
add "cli/F-huge|same|-|$k -F 99999999999 $mix"
add "cli/X-garbage|same|-|$k -X 1x $mix"
add "cli/p100|crash|-|-p 100 --fst4w -d 1 $mix"
add "cli/p45-two|crash|-|-p 45 --fst4w -d 1 $sil $mix"
add "cli/p100-silent|same|-|-p 100 --fst4w $sil"
add "cli/nonexistent|same|-|$k nothere-261012_1200.wav"
add "cli/nonexistent-second|same|-|$k -d 1 $mix nothere-261012_1200.wav"
add "cli/empty-arg|same|-|$k ''"
for f in rate11025 stereo bits8 header short tiny text; do add "cli/wav-$f|same|-|$k -d 1 @I/$f-261012_1200.wav"; done
add "cli/not-wav-name|same|-|$k -d 1 @I/notes.txt"
add "cli/c-slash|same|-|$k -d 3 -c K1ABC/P $h2"
add "cli/c-calls|same|-|$k -d 3 -c W1XY -x K9AN $h1"
add "cli/c-b|same|-|$k -d 3 -c b -x b $h2"
add "cli/x-b|same|-|$k -d 3 -x b $mix"
add "cli/c-short|same|-|$k -d 3 -c K1 $mix"
add "cli/q|same|-|$k -d 1 -q $mix"
add "cli/Y|same|-|-p 120 -Y -d 3 $h1 $h2"
add "cli/t-dir|same|-|$k -d 1 -t . $mix"
for x in 0 256 512 1024 1792; do add "cli/X$x|same|-|$k -d 1 -X $x $mix"; done
for F in 10 200; do add "cli/F$F|same|-|-p 120 --fst4w -f 1500 -F $F -d 1 $mix"; done
add "cli/f1450|same|-|-p 120 --fst4w -f 1450 -F 50 -d 1 $mix"
add "cli/no-f-F|same|-|-p 120 --fst4w -d 1 $mix"
for d in 0 2 7; do add "cli/d$d|same|-|$k -d $d $mix"; done
add "cli/d-1|usage|-|$k -d -1 $mix"  # a value starting with '-': getopt says it is missing
add "cli/p60-on-120|same|-|-p 60 --fst4w -d 3 $mix"
add "cli/no-p-on-120|same|-|--fst4w -d 3 $mix"
add "cli/p300-on-120|same|-|-p 300 --fst4w -d 1 $mix"

# 3. every recording with the usual options (-f 1500 -F 100 -d 3) and a plain run
for f in "$repo"/$S/*.wav "$repo"/recordings/fst4w-real/*.wav "$repo"/recordings/40m-wspr-20261010-1914/*.wav; do
    [[ -f $f ]] || continue
    r=${f#"$repo"/} b=$(basename "$f" .wav)
    case $b in
        p15-*) p=15 ;; p30-*) p=30 ;; p60-*) p=60 ;; f5-*) p=300 ;; f15-*) p=900 ;; f30-*|201230_*) p=1800 ;; *) p=120 ;;
    esac
    add "matrix/$b/skim|same|-|-p $p --fst4w -f 1500 -F 100 -d 3 $r"
    add "matrix/$b/plain|same|-|-p $p --fst4w $r"
done

# 4. FT8 (jt9 -8). The FT8 oracle's jt9 cases and every file-mode fixture
# case (scripts/ft8-filemode-oracle.sh: its prep.txt as the seed steps).
for a in $(find "$repo/tests/data/ft8/cases" -path '*/jt9-*' -name args.txt | sort); do
    d=$(dirname "$a"); name=${d#"$repo/tests/data/ft8/cases/"}
    add "ft8-oracle/$name|same-bare|-|$(cut -d' ' -f2- "$a")"
done
for a in $(find "$repo/tests/data/ft8/filemode/cases" -name args.txt | sort); do
    d=$(dirname "$a"); name=${d#"$repo/tests/data/ft8/filemode/cases/"}
    prep=$(cat "$d/prep.txt"); [[ $prep == - ]] || prep="steps:$prep"
    add "ft8-filemode/$name|same-bare|$prep|$(cut -d' ' -f2- "$a")"
done
# Every recorded period alone at depth 2, both sequences whole at depths 1
# and 3 (blank calls, the skimmer's band) and with jt9's defaults; every
# generated and real period at depth 3 and plain.
F=recordings/ft8-live
for f in "$repo"/$F/40m-20261012-1934/*.wav "$repo"/$F/20m-20261012-1934/*.wav; do
    [[ -f $f ]] || continue
    r=${f#"$repo"/} b=$(basename "$(dirname "$f")" | cut -c1-3)-$(basename "$f" .wav)
    add "ft8-live/$b/d2|same-bare|-|-8 -d 2 $r"
done
for band in 40m 20m; do
    seq=$(cd "$repo" && ls $F/$band-20261012-1934/*.wav 2>/dev/null | tr '\n' ' ')
    [[ -n $seq ]] || continue
    for d in 1 3; do add "ft8-live/$band-sequence/d$d|same-bare|-|-8 -L 100 -H 3200 -c b -x b -d $d $seq"; done
    add "ft8-live/$band-sequence/defaults|same-bare|-|-8 $seq"
done
for f in "$repo"/recordings/ft8-sim/*.wav "$repo"/recordings/ft8-real/*.wav; do
    [[ -f $f ]] || continue
    r=${f#"$repo"/} b=$(basename "$f" .wav)
    add "ft8-matrix/$b/d3|same-bare|-|-8 -d 3 $r"
    add "ft8-matrix/$b/plain|same-bare|-|-8 $r"
done
# Odd inputs and several files.
k8="-8 -d 1"
for f in rate11025 stereo bits8 header short tiny text list odd data-first long-data short-fmt tone30 cut13; do
    for w in "$work"/inputs/$f-*.wav; do [[ -f $w ]] && add "ft8-cli/wav-$f|same-bare|-|$k8 @I/$(basename "$w")"; done
done
L4=recordings/ft8-live/40m-20261012-1934
# Without the recorded 40 m periods (a copy of the decoder without them), the
# command-line cases run on three generated busy periods under their names.
if [[ ! -f $repo/$L4/261012_193500.wav ]]; then
    mkdir -p "$work/inputs/l4"
    for p in 193500:140000 193515:140015 193530:140030; do
        cp "$repo/recordings/ft8-sim/busy-261013_${p#*:}.wav" "$work/inputs/l4/261012_${p%:*}.wav" 2> /dev/null
    done
    L4=@I/l4
fi
add "ft8-cli/not-wav-name|same-bare|-|$k8 @I/notes.txt"
add "ft8-cli/nonexistent|same-bare|-|$k8 nothere-261012_193500.wav"
add "ft8-cli/nonexistent-second|same-bare|-|-8 -d 2 $L4/261012_193500.wav nothere-261012_193515.wav"
add "ft8-cli/bad-rate-second|same-bare|-|-8 -d 2 $L4/261012_193500.wav @I/rate11025-261012_1200.wav"
add "ft8-cli/empty-arg|same-bare|-|$k8 ''"
add "ft8-cli/several-odd|same-bare|-|-8 -d 2 @I/header-261012_1200.wav @I/tiny-261012_1200.wav @I/cut13-261012_193515.wav $L4/261012_193530.wav @I/short-261012_1200.wav"
add "ft8-cli/same-file-thrice|same-bare|-|-8 -d 3 $L4/261012_193515.wav $L4/261012_193515.wav $L4/261012_193515.wav"
# The command line with -8.
for o in -c -x -G -g; do
    add "ft8-cli/empty$o|ioerror-bare|-|-8 -d 2 $o '' $L4/261012_193500.wav"
    add "ft8-cli/blank$o|ioerror-bare|-|-8 -d 2 $o @BLANK $L4/261012_193500.wav"
done
add "cli/empty-c|ioerror|-|-p 120 --fst4w -c '' $mix"
add "ft8-cli/p10-test12k|same-bare|-|-8 -p 10 -d 3 -c K1ABC -x W9XYZ recordings/ft8-real/261012_125100.wav"
add "ft8-cli/long-ft8|same-bare|-|--ft8 -d 2 $L4/261012_193500.wav"
add "ft8-cli/long-options|same-bare|-|-8 --lowest 500 x --highest 2500 y --QSOprog 2 z --my-grid FN20 w --his-grid EN37 v -d 2 $L4/261012_193500.wav"
add "ft8-cli/W-then-8|same-bare|-|-W -p 15 -8 -d 2 $L4/261012_193500.wav"
add "ft8-cli/8-then-W|same-bare|-|-8 -W -p 120 -d 1 @I/short-261012_1200.wav"
for o in "-Q 6" "-Q 1.5" "-L x" "-H 1e3" "-p 1801" "-p 0" "-d 3.0" "-f +" "-X 1x"; do
    add "ft8-cli/bad${o// /}|same-bare|-|-8 $o $L4/261012_193500.wav"
done
for o in "-d" "-Q" "-L" "-g"; do add "ft8-cli/missing$o|usage|-|-8 $o"; done
add "ft8-cli/f-negative|usage|-|-8 -f -50 $L4/261012_193500.wav"
n=0
for o in "-d 0" "-d 7" "-F 0" "-f 0" "-f 4000" "-L 4000" "-H 0" "-p 0.5" "-p 1800" "-p 14.399" "-X 8" "-X 7" "-X 64" "-q" \
         "-c K1ABC/P" "-c 'K1ABC/P'" "-x ,W9XYZ" "-c b -x W9XYZ" "-g EN" "-g 'EN37XX99'" "-c 'ABCDEFGHIJKLMNOP'" "-G ,FN20" \
         "-t ." "-a . -a ."; do
    n=$((n + 1)); add "ft8-cli/opt$n${o//[^A-Za-z0-9.-]/}|same-bare|-|-8 -d 2 $o $L4/261012_193500.wav"
done

# ---------------------------------------------------------------- running
run_one() {  # side spec
    local side=$1 spec=$2 name kind seed args
    IFS='|' read -r name kind seed args <<< "$spec"
    local dir=$work/$side/$name
    mkdir -p "$dir/cwd"
    # seed: an fst4w_calls.txt, or steps:mkdir:DIR;seed:PATH>NAME;...
    if [[ $seed == steps:* ]]; then
        local step steps
        IFS=';' read -r -a steps <<< "${seed#steps:}"
        for step in "${steps[@]}"; do
            case $step in
                mkdir:*) mkdir -p "$dir/cwd/${step#mkdir:}" ;;
                seed:*) step=${step#seed:}; cp "$repo/${step%%>*}" "$dir/cwd/${step#*>}" ;;
            esac
        done
    elif [[ $seed != - ]]; then
        cp "$repo/$seed" "$dir/cwd/fst4w_calls.txt"
    fi
    [[ $kind == *bare* ]] || touch "$dir/cwd/plotspec" "$dir/cwd/decdata"
    local argv=() tok
    [[ $side == new ]] && argv+=(--spread exact)
    for tok in $args; do
        case $tok in
            recordings/*|tests/*) argv+=("$repo/$tok") ;;
            @I/*) argv+=("$work/inputs/${tok#@I/}") ;;
            "''") argv+=("") ;;
            @BLANK) argv+=(" ") ;;
            *) argv+=("$tok") ;;
        esac
    done
    if [[ $side == new ]]; then
        (cd "$dir/cwd" && "$new" "${argv[@]}" > ../stdout 2> ../stderr < /dev/null)
    elif [[ -n $image ]]; then
        docker run --rm --network none -u "$(id -u):$(id -g)" -v "$repo:$repo:ro" -v "$work:$work" \
            -v "$(dirname "$stock"):$(dirname "$stock"):ro" -w "$dir/cwd" "$image" \
            sh -c 'exec "$@" > ../stdout 2> ../stderr < /dev/null' sh "$stock" "${argv[@]}"
    else
        (cd "$dir/cwd" && "$stock" "${argv[@]}" > ../stdout 2> ../stderr < /dev/null)
    fi
    echo $? > "$dir/status"
    # (and a crash's core file, where the system writes one)
    rm -f "$dir/cwd/plotspec" "$dir/cwd/decdata" "$dir"/cwd/core "$dir"/cwd/core.*
    find "$dir/cwd" \( -name timer.out -o -name jt9_wisdom.dat \) -delete
}

compare() {  # spec -> results/<name>
    local spec=$1 name kind seed args
    IFS='|' read -r name kind seed args <<< "$spec"
    local n=$work/new/$name s=$work/stock/$name out=$work/results/${name//\//:}
    local bad=0 report=
    cp "$n/stderr" "$n/stderr.cmp"
    if [[ $kind == usage* || $kind == version ]]; then
        cmp -s "$n/stderr.cmp" "$s/stderr" || { bad=1; report+=$(diff "$s/stderr" "$n/stderr.cmp")$'\n'; }
        if [[ $kind == version ]]; then
            grep -q '^jt9 ' "$n/stdout" && grep -q '^jt9 ' "$s/stdout" || { bad=1; report+="no version"$'\n'; }
        else
            grep -q '^ Usage: ' "$n/stdout" && grep -q '^ Usage: ' "$s/stdout" || { bad=1; report+="no usage"$'\n'; }
        fi
        cmp -s "$n/status" "$s/status" || { bad=1; report+="status $(cat "$s/status") / $(cat "$n/status")"$'\n'; }
    else
        cmp -s "$n/status" "$s/status" || { bad=1; report+="status: jt9 $(cat "$s/status"), jt9 $(cat "$n/status")"$'\n'; }
        cmp -s "$n/stdout" "$s/stdout" || { bad=1; report+="stdout:"$'\n'$(diff "$s/stdout" "$n/stdout" | head -20)$'\n'; }
        if [[ $kind == crash* ]]; then
            head -2 "$s/stderr" | cmp -s - "$n/stderr.cmp" || { bad=1; report+="stderr:"$'\n'$(diff <(head -2 "$s/stderr") "$n/stderr.cmp")$'\n'; }
        elif [[ $kind == ioerror* ]]; then
            head -2 "$s/stderr" | sed '1s|of file .*/jt9\.f90$|of file jt9.f90|' | cmp -s - "$n/stderr.cmp" ||
                { bad=1; report+="stderr:"$'\n'$(diff <(head -2 "$s/stderr") "$n/stderr.cmp")$'\n'; }
        else
            cmp -s "$n/stderr.cmp" "$s/stderr" || { bad=1; report+="stderr:"$'\n'$(diff "$s/stderr" "$n/stderr.cmp" | head -20)$'\n'; }
        fi
        local f
        for f in $( (cd "$s/cwd" && find . -type f; cd "$n/cwd" && find . -type f) | sort -u); do
            [[ $ignore == *" ${f#./} "* ]] && continue
            if ! cmp -s "$s/cwd/$f" "$n/cwd/$f"; then
                bad=1
                report+="$f:"$'\n'$(diff "$s/cwd/$f" "$n/cwd/$f" 2>&1 | head -10)$'\n'
            fi
        done
    fi
    if (( bad )); then { echo "FAIL $name"; echo "$report"; } > "$out"; else echo "PASS $name" > "$out"; fi
}

selected=()
for spec in "${cases[@]}"; do
    name=${spec%%|*}
    [[ -n $only && ! $name =~ $only ]] && continue
    # inputs present?
    skip=
    for tok in ${spec##*|}; do [[ ($tok == recordings/* || $tok == tests/*) && ! -e $repo/$tok ]] && skip=$tok; done
    if [[ -n $skip ]]; then echo "SKIP ${name} ($skip missing)" > "$work/results/${name//\//:}"; continue; fi
    selected+=("$spec")
done
running=0
for spec in "${selected[@]}"; do
    for side in new stock; do
        if (( running >= jobs_max )); then wait -n; running=$((running - 1)); fi
        run_one "$side" "$spec" &
        running=$((running + 1))
    done
done
wait
for spec in "${selected[@]}"; do compare "$spec"; done

cat "$work"/results/* | grep -E '^(PASS|FAIL|SKIP)' | sort | awk '{print $1}' | uniq -c >&2
fails=0
for r in "$work"/results/*; do
    if grep -q '^FAIL' "$r"; then cat "$r"; fails=1; fi
done
grep -h '^SKIP' "$work"/results/* 2>/dev/null
exit $fails
