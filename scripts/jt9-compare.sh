#!/usr/bin/env bash
# Run our jt9 (NEW: src/jt9_main.cpp) and WSJT-X's own jt9 (STOCK:
# scripts/fst4w-oracle.sh builds it, at 567ad29ce) side by side, each in
# fresh directories of its own, and compare everything they print and write:
# stdout, stderr, the exit status, and the files left in the run's
# directories (fst4_decodes.dat, decoded.txt, fst4w_calls.txt, fort.52; not
# timer.out or jt9_wisdom.dat, the clock's and FFTW's, which NEW does not
# write).
#
#   scripts/jt9-compare.sh NEW_JT9 STOCK_JT9 [--docker IMAGE]
#
#   --docker IMAGE                STOCK runs in IMAGE (the oracle's jt9 needs
#                                 its libraries: wsjtx-oracle), the repository
#                                 and work dir mounted at their own paths
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
# RIFF, a short period, a tiny one.
python3 -I - "$work/inputs" <<'EOF'
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
open(f"{out}/notes.txt", "w").write("not a wav\n")
EOF

# ---------------------------------------------------------------- cases
# A case: name|kind|seed|args. kind: same, usage or crash (above). seed: an
# fst4w_calls.txt put in the run's dir first (repository-relative) or -.
# args split on blanks; recordings/, tests/ and @I/ (the odd inputs) made
# absolute; '' an empty argument. Each runs in its own fresh dir, the cwd.
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

# ---------------------------------------------------------------- running
run_one() {  # side spec
    local side=$1 spec=$2 name kind seed args
    IFS='|' read -r name kind seed args <<< "$spec"
    local dir=$work/$side/$name
    mkdir -p "$dir/cwd"
    [[ $seed != - ]] && cp "$repo/$seed" "$dir/cwd/fst4w_calls.txt"
    [[ $kind == *bare* ]] || touch "$dir/cwd/plotspec" "$dir/cwd/decdata"
    local argv=() tok
    [[ $side == new ]] && argv+=(--spread exact)
    for tok in $args; do
        case $tok in
            recordings/*|tests/*) argv+=("$repo/$tok") ;;
            @I/*) argv+=("$work/inputs/${tok#@I/}") ;;
            "''") argv+=("") ;;
            *) argv+=("$tok") ;;
        esac
    done
    if [[ $side == new ]]; then
        (cd "$dir/cwd" && "$new" "${argv[@]}" > ../stdout 2> ../stderr < /dev/null)
    elif [[ -n $image ]]; then
        docker run --rm -i -u "$(id -u):$(id -g)" -v "$repo:$repo:ro" -v "$work:$work" \
            -v "$(dirname "$stock"):$(dirname "$stock"):ro" -w "$dir/cwd" "$image" \
            "$stock" "${argv[@]}" > "$dir/stdout" 2> "$dir/stderr" < /dev/null
    else
        (cd "$dir/cwd" && "$stock" "${argv[@]}" > ../stdout 2> ../stderr < /dev/null)
    fi
    echo $? > "$dir/status"
    # (and a crash's core file, where the system writes one)
    rm -f "$dir/cwd/timer.out" "$dir/cwd/jt9_wisdom.dat" "$dir/cwd/plotspec" "$dir/cwd/decdata" "$dir"/cwd/core "$dir"/cwd/core.*
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
