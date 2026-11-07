#!/usr/bin/env bash
# The A/B of Task 8: WSJT-X 3.0.2's jt9 against the stored 3.2.0-rc1 outputs
# (tests/data/msk144/cases, which the port matches), run by run.
#
#   scripts/msk144-ab-302.sh [JT9]     JT9 defaults to /opt/git/sigmond/wsjtx/build-sigmond/jt9
#   MSK144_ODD=DIR                     the odd WAVs (@I/ in args.txt), default ~/.cache/msk144-oracle/inputs,
#                                      where scripts/msk144-oracle.sh corpus puts them
#
# Each run of the corpus goes to JT9 in a fresh directory, one file a process
# (as the references were made; 3.0.2 has no --feed, so a --feed run is left
# to its jt9-files twin). The shmsg and swl variants are left out: 3.0.2's
# decode_msk144 has bshmsg and bswl fixed false as well, and we did not patch
# or rebuild it. OMP_DYNAMIC and OMP_THREAD_LIMIT are unset. Compared: stdout
# (the repository's directory taken out of the paths), the status, stderr and
# the files left in the directory. Prints SAME or DIFF per run, then a count.
set -u
repo=$(cd "$(dirname "$0")/.." && pwd)
jt9=${1:-/opt/git/sigmond/wsjtx/build-sigmond/jt9}
odd=${MSK144_ODD:-$HOME/.cache/msk144-oracle/inputs}
work=$(mktemp -d "${TMPDIR:-/tmp}/msk144-ab.XXXXXX")
unset OMP_DYNAMIC OMP_THREAD_LIMIT
same=0 diff=0 left=0
for a in $(find "$repo/tests/data/msk144/cases" -path '*/jt9-*' -name args.txt | sort); do
    d=$(dirname "$a"); v=${d##*/jt9-}; c=$(basename "$(dirname "$d")")
    [[ $v == shmsg || $v == swl || $v == feed ]] && { left=$((left + 1)); continue; }
    grep -q -e '--msk-shmsg' -e '--msk-swl' "$a" && { left=$((left + 1)); continue; }
    read -r -a toks <<< "$(cut -d' ' -f2- "$a")"
    opts=() files=()
    for t in "${toks[@]}"; do
        case $t in recordings/*) files+=("$repo/$t") ;; @I/*) files+=("$odd/${t#@I/}") ;; *) opts+=("$t") ;; esac
    done
    r=$work/$c-$v; mkdir -p "$r/cwd"; : > "$r/stdout"; : > "$r/stderr"; st=0
    for f in "${files[@]}"; do
        # the reference's one process a file: its status is the last, its stderr the sum
        (cd "$r/cwd" && "$jt9" "${opts[@]}" "$f" >> ../stdout 2>> ../stderr < /dev/null); st=$?
    done
    echo $st > "$r/status"
    find "$r/cwd" \( -name timer.out -o -name jt9_wisdom.dat \) -delete
    bad=
    cmp -s <(sed "s|$repo/||g;s|$odd/|@I/|g" "$r/stdout") "$d/stdout.txt" || bad+=" stdout"
    cmp -s "$r/status" "$d/status.txt" || bad+=" status($(cat "$d/status.txt")->$st)"
    if [[ -f $d/stderr.txt ]]; then cmp -s <(sed "s|$repo/||g;s|$odd/|@I/|g" "$r/stderr") "$d/stderr.txt" || bad+=" stderr"
    else [[ -s $r/stderr ]] && bad+=" stderr"; fi
    for f in $( (cd "$r/cwd" && find . -type f; [[ -f $d/decoded.txt ]] && echo ./decoded.txt) | sort -u); do
        cmp -s "$r/cwd/$f" "$d/$f" 2> /dev/null || bad+=" $f"
    done
    if [[ -z $bad ]]; then echo "SAME $c/$v"; same=$((same + 1))
    else
        echo "DIFF $c/$v:$bad"; diff=$((diff + 1))
        diff "$d/stdout.txt" <(sed "s|$repo/||g;s|$odd/|@I/|g" "$r/stdout") | sed 's/^/    /' | head -12
        [[ -s $r/stderr ]] && { echo "    stderr:"; sed "s|$repo/||g;s|$odd/|@I/|g" "$r/stderr" | head -4 | sed 's/^/    /'; }
    fi
done
echo "msk144-ab-302: $same same, $diff different, $left left out (shmsg, swl, feed, any --msk-* flag); runs in $work"
