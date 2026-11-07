#!/usr/bin/env bash
# CPU time of a jt9 on the MSK144 corpus: every run of tests/data/msk144/cases
# once, one file a process (so no --feed), in a fresh directory; the whole sweep
# timed by the shell's `time`, REPEATS times (default 5).
#
#   scripts/msk144-cpu.sh JT9 [REPEATS]
#
# The runs timed: those without --msk-shmsg or --msk-swl (the only ones 3.0.2's
# jt9 can be given), less the 24 runs on which 3.0.2's jt9 stops early or
# differently (all of freq-edge, 7 of which differ: bounds-check abort;
# odd-bits8, odd-rate11025, odd-stereo: it reads what 3.2.0-rc1 refuses), and
# less the cases the final fix wave added (odd-trunc, odd-stream, real-uacnj-*):
# the README's table is the corpus of 276 processes as it stood at Task 8. The
# same work for the port, 3.2.0-rc1 and 3.0.2.
#
# Prints, per repeat, user, system and real seconds, and the load average before
# and after. Set OMP_NUM_THREADS=1 for a single thread; MSK144_ODD is the odd
# WAVs' directory (scripts/msk144-ab-302.sh). Run it where JT9 runs: for the
# reference, inside the wsjtx-oracle image (scripts/msk144-oracle.sh exec).
set -u
repo=$(cd "$(dirname "$0")/.." && pwd)
jt9=$1; reps=${2:-5}
odd=${MSK144_ODD:-$HOME/.cache/msk144-oracle/inputs}
unset OMP_DYNAMIC OMP_THREAD_LIMIT
work=$(mktemp -d "${TMPDIR:-/tmp}/msk144-cpu.XXXXXX"); mkdir "$work/cwd"
runs=$work/runs; : > "$runs"
skip='^(freq-edge/|odd-bits8/|odd-rate11025/|odd-stereo/|odd-trunc/|odd-stream/|real-uacnj-)'
for a in $(find "$repo/tests/data/msk144/cases" -path '*/jt9-*' -name args.txt | sort); do
    d=$(dirname "$a"); v=${d##*/jt9-}; c=$(basename "$(dirname "$d")")
    [[ $v == feed ]] && continue
    grep -q -e '--msk-shmsg' -e '--msk-swl' "$a" && continue
    [[ $c/$v =~ $skip ]] && continue
    read -r -a toks <<< "$(cut -d' ' -f2- "$a")"
    opts=() files=()
    for t in "${toks[@]}"; do
        case $t in recordings/*) files+=("$repo/$t") ;; @I/*) files+=("$odd/${t#@I/}") ;; *) opts+=("$t") ;; esac
    done
    for f in "${files[@]}"; do echo "${opts[*]}|$f" >> "$runs"; done
done
n=$(wc -l < "$runs")
echo "msk144-cpu: $jt9: $n process runs a sweep, $reps repeats, OMP_NUM_THREADS=${OMP_NUM_THREADS:-unset}"
sweep() {
    cd "$work/cwd" || exit 1
    while IFS='|' read -r o f; do
        # shellcheck disable=SC2086
        "$jt9" $o "$f" > /dev/null 2>&1 < /dev/null
    done < "$runs"
}
TIMEFORMAT='%U %S %R'
for r in $(seq "$reps"); do
    la0=$(cut -d' ' -f1-3 /proc/loadavg)
    t=$( { time sweep; } 2>&1 )
    la1=$(cut -d' ' -f1-3 /proc/loadavg)
    echo "repeat $r: user sys real = $t  loadavg before [$la0] after [$la1]"
done
rm -rf "$work"
