#!/usr/bin/env bash
# Run a wsprd built on the in-process decoder (NEW) and WSJT-X's own wsprd
# (STOCK, scripts/wsprd-oracle.sh --build-only builds it) side by side, each
# in fresh directories of its own, and compare everything they print and
# write, byte for byte: stdout, stderr, exit status, and after every run the
# whole directory (ALL_WSPR.TXT, wspr_spots.txt, hashtable.txt, the -c file,
# whatever lands in an -a dir).
#
#   scripts/wsprd-compare.sh NEW_WSPRD STOCK_WSPRD [--spread-oracle N_WSPRD]
#
#   WSPRD_COMPARE_WORK=/tmp/...   where to run (default: a temp dir, kept)
#   WSPRD_COMPARE_JOBS=n          runs at once (default: nproc)
#   WSPRD_COMPARE_ONLY=regex      only the cases whose names match
#
# The cases:
#   oracle/...  every case in tests/data/wsprd/cases, replayed as
#               scripts/wsprd-oracle.sh ran it (-a . first, the data dir the
#               cwd, hashtable.in, files.txt's files and read-only ones);
#   cli/...     the command line: no arguments, bad and missing options,
#               clustered flags and attached values, operands anywhere, --,
#               files that are not there or not .wav/.c2, -a dirs, odd -C -e
#               -o -z -f values, -m, -t, and -0 with framed PCM on stdin;
#   matrix/...  every recording (tests/data/wsprd/inputs, recordings/wsprd-sim,
#               the 40 m periods, recordings/wspr2-150426_0918-14.0956.wav)
#               under each option alone and two common sets, and the five
#               40 m periods in one data dir, run after run;
#   spread/...  with --spread-oracle only: every recording with
#               "-n -c -C 500 -o 4 -d -f MHz" and "-n -f MHz", NEW against
#               N_WSPRD, another wsprd taking -n and printing the spread
#               column: ALL_WSPR.TXT's lines compared without
#               their last column (" %6.3f", the spread, which NEW's lines
#               must end in; the values are shown side by side, not
#               compared), everything else byte for byte.
# A NEW whose usage has -n's line (wsprd's does): STOCK's stdout
# with that line put after its usage's -z line is what NEW must print (each
# run's stdout only; every other byte still STOCK's). A NEW without it (STOCK
# against itself): no allowance at all.
# Cases on a missing input are skipped. Both binaries run as argv[0]
# "wsprd" (getopt's messages name it). Not compared to the byte, being the
# clock's or FFTW's: wspr_timer.out (each number, nan or inf, with the
# blanks before it, made " N": its lines, labels and fields still are) and wspr_wisdom.dat (only whether it
# is there). Files over 64 KiB (the -c file) are compared by sha256.
# Prints PASS/FAIL/SKIP a case, the diff under each FAIL; exits 1 on any
# FAIL. Everything is left in the work dir: new/, stock/ (and spread/), one
# directory a case, run.N/ under it holding run N's outputs and the data
# dir as that run left it (files/), cwd.live/ the dir itself.
set -uo pipefail

usage() { echo "usage: $0 NEW_WSPRD STOCK_WSPRD [--spread-oracle N_WSPRD]" >&2; exit 2; }
(( $# >= 2 )) || usage
new=$(realpath "$1") stock=$(realpath "$2"); shift 2
spread= with_n=0
while (( $# )); do
    case $1 in
        --spread-oracle) (( $# >= 2 )) || usage; spread=$(realpath "$2"); shift ;;
        *) usage ;;
    esac
    shift
done
for b in "$new" "$stock" ${spread:+"$spread"}; do
    [[ -x $b ]] || { echo "wsprd-compare: $b is not executable" >&2; exit 2; }
done
# NEW's usage (no arguments) with -n's line: its stdout to expect STOCK's plus that line
new_usage=$("$new" 2>/dev/null) || true  # wsprd exits 1 after its usage
grep -qx -- '       -n no drift, does not attempt to follow a drifting frequency' <<<"$new_usage" && with_n=1
repo=$(cd "$(dirname "$0")/.." && pwd)
work=${WSPRD_COMPARE_WORK:-$(mktemp -d "${TMPDIR:-/tmp}/wsprd-compare.XXXXXX")}
jobs_max=${WSPRD_COMPARE_JOBS:-$(nproc)}
only=${WSPRD_COMPARE_ONLY:-}
mkdir -p "$work"
if [[ -d $work/new || -d $work/stock || -d $work/spread || -d $work/streams ]]; then
    chmod -R u+w "$work/new" "$work/stock" "$work/spread" "$work/streams" 2>/dev/null
    rm -rf "${work:?}/new" "${work:?}/stock" "${work:?}/spread" "${work:?}/streams" "${work:?}/results" "${work:?}/expected"
fi
mkdir -p "$work/new" "$work/stock" "$work/streams" "$work/results"
echo "wsprd-compare: work dir $work" >&2

# ---------------------------------------------------------------- cases
# A case: name|seed|files|runs. seed: a hashtable.txt to start from
# (repository-relative) or -. files, ';'-separated, as the oracle's
# files.txt ("SRC NAME [ro]", SRC repository-relative or - for an empty
# file), plus "dir NAME" for a directory. runs, ';'-separated, each one
# wsprd's argv, split on blanks: a token "<FILE" is stdin (else /dev/null),
# '' an empty argument, recordings/... and tests/... made absolute, @W the
# streams dir. The oracle's cases get "-a ." in front, as it ran them; the
# others run as written, the case dir their cwd.
cases=()
add() { cases+=("$1"); }

# 1. tests/data/wsprd/cases
casesdir=$repo/tests/data/wsprd/cases
while IFS= read -r a; do
    d=$(dirname "$a"); name=${d#"$casesdir"/}
    seed=-; [[ -f $d/hashtable.in ]] && seed=${d#"$repo"/}/hashtable.in
    files=; [[ -f $d/files.txt ]] && files=$(paste -sd';' "$d/files.txt")
    runs=$(sed 's/^/-a . /' "$a" | paste -sd';')
    add "oracle/$name|$seed|$files|$runs"
done < <(find "$casesdir" -name args.txt | sort)

# 2. the command line
T=tests/data/wsprd/inputs
Rd=recordings/wsprd-sim
R40=recordings/40m-wspr-20261010-1914
st=$T/sim-types-261010_1200.c2
sw=$Rd/sim-wav2-261010_1218.wav
ct=$T/cut-tiny-261010_1218.wav
add "cli/no-args|-||"
add "cli/h|-||-h"
add "cli/question|-||-?"
add "cli/unknown|-||-x $st"
add "cli/unknown-in-cluster|-||-cxq $st"
add "cli/long-option|-||--help"
add "cli/stream-long|-||--stream"
for o in a C e f o z; do add "cli/missing-$o|-||-$o"; done
add "cli/missing-C-after-file|-||$st -C"
add "cli/missing-o-after-flags|-||-c -d -o"
add "cli/cluster-wd|-||-wd -f 14.0956 $st"
add "cli/cluster-qs|-||-qs $st"
add "cli/cluster-C500|-||-C500 -c $st"
add "cli/cluster-o4|-||-o4 -C500 $st"
add "cli/cluster-cdC500|-||-cdC500 $st"
add "cli/cluster-cdC-500|-||-cdC 500 $st"
add "cli/cluster-HJBqv|-||-HJBqv $st"
add "cli/cluster-a|-|dir data|-cadata $st"
add "cli/operand-first|-||$st -c -C 500 -o 2"
add "cli/operand-middle|-||-c $st -C 500 -d"
add "cli/dashdash|-||-c -- $st"
add "cli/dashdash-then-option|-||-- -c $st"
add "cli/dashdash-alone|-||-c --"
add "cli/two-operands|-||$st $T/sim-sweep-261010_1202.c2"
add "cli/two-operands-bad-first|-||nothere-261010_1200.c2 $st"
add "cli/nonexistent-c2|-||nothere-261010_1200.c2"
add "cli/nonexistent-wav|-||nothere-261010_1200.wav"
add "cli/nonexistent-wav-dir|-||nodir/nothere-261010_1200.wav"
add "cli/suffix-txt|-|$st sim-261010_1200.txt|sim-261010_1200.txt"
add "cli/suffix-none|-|$st sim-261010_1200|sim-261010_1200"
add "cli/suffix-inside|-|$st x.c2.txt|x.c2.txt"
add "cli/empty-operand|-||-c ''"
add "cli/dash-operand|-||-"
add "cli/a-nonexistent|-||-a nodir $st"
add "cli/a-nonexistent-no-file|-||-a nodir"
add "cli/a-relative|-|dir data|-a data -c $st"
add "cli/a-relative-slash|-|dir data|-a data/ -c $st"
add "cli/a-twice|-|dir data|-a nodir -a data $st"
add "cli/a-empty|-||-a '' $st"
add "cli/a-file|-|- plain|-a plain $st"
add "cli/a-readonly|-|dir rodir;- rodir/x;rodir rodir ro|-a rodir $st"
add "cli/e-garbage|-||-e x -f 14.0956 $st"
add "cli/e-negative|-||-e -3.5 -d $st"
add "cli/e-wav|-||-e 7 -f 7.0386 $sw"
add "cli/z-garbage|-||-z x $T/sim-sweep-261010_1202.c2"
add "cli/z-negative|-||-z -0.5 $T/sim-sweep-261010_1202.c2"
add "cli/z-trailing|-||-z 0.6abc $T/sim-sweep-261010_1202.c2"
add "cli/C-garbage|-||-C x $st"
# strtoul wraps a minus sign; the cast to unsigned keeps the low 32 bits:
# 500 (a plain -5 is 4294967291 cycles a bit, hours on noise)
add "cli/C-negative|-||-C -4294966796 $st"
add "cli/C-zero|-||-C 0 -d $st"
add "cli/C-hex|-||-C 0x100 $st"
add "cli/o-garbage|-||-o x $st"
add "cli/o-negative|-||-o -1 $st"
add "cli/o-six|-||-o 6 -C 500 $T/sim-sweep-261010_1202.c2"
add "cli/f-garbage-wav|-||-f x $sw"
add "cli/f-negative-wav|-||-f -7.0386 $sw"
add "cli/f-c2-ignored|-||-f 99 $st"
add "cli/no-f-wav|-||$sw"
add "cli/m-wav2|-||-m -f 7.0386 $sw"
add "cli/m-c2|-||-m $st"
add "cli/t-flag|-||-t $st"
add "cli/t-cluster|-||-tc $st"
add "cli/v-H|-||-v -H $st"
add "cli/H-seeded|tests/data/wsprd/hashtable-odd.txt||-H $st"
add "cli/seeded-odd|tests/data/wsprd/hashtable-odd.txt||-c $st"
add "cli/cut-tiny-default|-||-f 7.0386 $ct"
add "cli/empty-wav|-|- empty-261010_1200.wav|empty-261010_1200.wav"
add "cli/empty-c2|-|- empty-261010_1200.c2|empty-261010_1200.c2"
add "cli/timer-empty|-|- wspr_timer.out|$st"
add "cli/timer-garbage|-|tests/data/wsprd/hashtable-odd.txt wspr_timer.out|$st"
add "cli/rerun-same-dir|-||-c $st;-c $st;-H $st"
# -0: framed PCM on stdin (the streams are made below).
add "cli/stream-wav2|-||-0 <@W/wav2.bin"
add "cli/stream-wav2-options|-||-0 -c -C 500 -o 3 -d <@W/wav2.bin"
add "cli/stream-wav2-e|-||-e 3 -0 <@W/wav2.bin"
add "cli/stream-wav2-operand|-||-0 $st <@W/wav2.bin"
add "cli/stream-wav2-nohalt|-||-0 <@W/wav2-nohalt.bin"
add "cli/stream-wav15|-||-0 <@W/wav15.bin"
add "cli/stream-empty|-||-0"
add "cli/stream-bad-magic|-||-0 <@W/bad-magic.bin"
add "cli/stream-bad-rate|-||-0 <@W/bad-rate.bin"
add "cli/stream-no-configure|-||-0 <@W/no-configure.bin"
add "cli/stream-incomplete|-||-0 <@W/incomplete.bin"
add "cli/stream-short-pcm|-||-0 <@W/short-pcm.bin"
add "cli/stream-odd-frame|-||-0 <@W/odd-frame.bin"
add "cli/stream-bad-type|-||-0 <@W/bad-type.bin"
add "cli/stream-no-t|-||-0 <@W/no-t.bin"
add "cli/stream-short-header|-||-0 <@W/short-header.bin"
add "cli/stream-short-body|-||-0 <@W/short-body.bin"
add "cli/stream-too-large|-||-0 <@W/too-large.bin"
add "cli/stream-unknown-frame|-||-0 <@W/unknown-frame.bin"
add "cli/stream-a-nonexistent|-||-a nodir -0 <@W/wav2.bin"

# 3. the option matrix
dial() {
    case $1 in
        *150426_0918*|*sweep-*|*types*|*crowd*|*strong*) echo 14.0956 ;;
        *sweep2*) echo 10.1387 ;;
        *drift*|*wav2*|*cut-*|*40m-wspr*) echo 7.0386 ;;
        *wide*) echo 3.5686 ;;
        *edge*) echo 1.8366 ;;
        *noise*) echo 28.1246 ;;
        *wspr15*|*wav15*) echo 0.4742 ;;
    esac
}
recordings=()
for f in "$repo"/$T/*.c2 "$repo"/$T/*.wav "$repo"/$Rd/*.c2 "$repo"/$Rd/*.wav \
         "$repo"/$R40/*.wav "$repo"/recordings/wspr2-150426_0918-14.0956.wav; do
    [[ -f $f ]] && recordings+=("${f#"$repo"/}")
done
optsets=("-B" "-c" "-d" "-H" "-J" "-q" "-s" "-v" "-w"
         "-o 0" "-o 1" "-o 2" "-o 3" "-o 4" "-o 5" "-C 500" "-z 0.6" "-e 2")
recname() {  # recording -> a short case name
    local s=${1#recordings/}; s=${s#tests/data/wsprd/}; s=${s%.*}; echo "${s//\//:}"
}
for r in "${recordings[@]}"; do
    m=; [[ $r == *wav15* ]] && m="-m "
    f=$(dial "$r")
    for o in "${optsets[@]}"; do
        add "matrix/$(recname "$r")/${o// /}|-||${m}$o -f $f $r"
    done
    add "matrix/$(recname "$r")/common-skim|-||${m}-c -C 500 -o 4 -d -f $f $r"
    add "matrix/$(recname "$r")/common-plain|-||${m}-o 0 -q -s -H -f $f $r"
done
join() { local IFS=';'; echo "$*"; }
if [[ -f $repo/$R40/261010_1914.wav ]]; then
    for o in "-c -C 500 -o 4 -d" "" "-o 0 -q -s" "-w -v" "-B -J -o 1"; do
        seq=()
        for t in 1914 1916 1918 1920 1922; do seq+=("$o -f 7.0386 $R40/261010_$t.wav"); done
        n=${o// /}; add "matrix/40m-sequence/${n:-default}|-||$(join "${seq[@]}")"
    done
fi

# 4. -n against N_WSPRD
spread_cases=()
if [[ -n $spread ]]; then
    for r in "${recordings[@]}"; do
        m=; [[ $r == *wav15* ]] && m="-m "
        f=$(dial "$r")
        spread_cases+=("spread/$(recname "$r")/skim|-||${m}-n -c -C 500 -o 4 -d -f $f $r")
        spread_cases+=("spread/$(recname "$r")/plain|-||${m}-n -f $f $r")
    done
fi

# ---------------------------------------------------------------- streams
# wsprd_stream.c's framing: "WSJT", fmt, channels, kHz (LE u16), then
# frames [type][len LE u32][body]: 1 audio (int16 LE), 2 control JSON.
python3 -I - "$work/streams" "$repo/$Rd/sim-wav2-261010_1218.wav" "$repo/$Rd/sim-wav15-261010_1230.wav" <<'EOF'
import os, struct, sys
out, wav2, wav15 = sys.argv[1:4]
def hdr(magic=b"WSJT", fmt=0, ch=1, khz=12): return magic + bytes([fmt, ch]) + struct.pack("<H", khz)
def frame(t, body): return bytes([t]) + struct.pack("<I", len(body)) + body
def ctl(s): return frame(2, s.encode())
def pcm(path):
    if not os.path.exists(path): return None
    with open(path, "rb") as f: return f.read()[44:]
def chunks(b, n=48000 * 2):
    return b"".join(frame(1, b[i:i + n]) for i in range(0, len(b), n))
cfg2 = ctl('{"t":"configure","date":"261010","time":"1218","dialfreq":7.0386,"wspr_type":2,"mycall":"K1JT","mygrid":"FN20"}')
cfg15 = ctl('{"t":"configure","date":"261010","time":"1230","dialfreq":0.4742,"wspr_type":15}')
halt = ctl('{"t":"halt"}')
w2, w15 = pcm(wav2), pcm(wav15)
files = {
    "bad-magic.bin": b"WSJX" + bytes([0, 1, 12, 0]),
    "bad-rate.bin": hdr(khz=48) + cfg2 + halt,
    "no-configure.bin": hdr() + frame(1, b"\0\0" * 100) + halt,
    "incomplete.bin": hdr() + ctl('{"t":"configure","date":"261010"}') + halt,
    "short-pcm.bin": hdr() + cfg2 + frame(1, b"\1\0" * 1000) + halt,
    "odd-frame.bin": hdr() + cfg2 + frame(1, b"\0\0\0"),
    "bad-type.bin": hdr() + ctl('{"t":"configure","date":"261010","time":"1218","dialfreq":7.0386,"wspr_type":3}'),
    "no-t.bin": hdr() + ctl('{"x":1}') + cfg2 + frame(1, b"\0\0" * 10) + halt,
    "short-header.bin": b"WSJ",
    "short-body.bin": hdr() + cfg2 + bytes([1]) + struct.pack("<I", 1000) + b"\0" * 10,
    "too-large.bin": hdr() + bytes([2]) + struct.pack("<I", 5000),
}
if w2 is not None:
    files["wav2.bin"] = hdr() + cfg2 + chunks(w2) + halt
    files["wav2-nohalt.bin"] = hdr() + cfg2 + chunks(w2)
    files["unknown-frame.bin"] = hdr() + frame(9, b"abc") + cfg2 + chunks(w2) + halt
if w15 is not None:
    files["wav15.bin"] = hdr() + cfg15 + chunks(w15) + halt
for k, v in files.items():
    with open(os.path.join(out, k), "wb") as f: f.write(v)
EOF

# ---------------------------------------------------------------- running
# snapshot DIR OUT: DIR's files into OUT as compared (see the top).
snapshot() {
    local dir=$1 out=$2 f rel
    mkdir -p "$out"
    (cd "$dir" && find . -mindepth 1 -type d) | while IFS= read -r rel; do mkdir -p "$out/$rel"; done
    (cd "$dir" && find . -mindepth 1 ! -type d) | while IFS= read -r rel; do
        f=$dir/$rel
        case $rel in
            */wspr_wisdom.dat) echo "(FFTW wisdom, not compared)" > "$out/$rel" ;;
            */wspr_timer.out) sed -E 's/ *-?([0-9]+\.[0-9]+|nan|inf)/ N/g' "$f" > "$out/$rel" 2>&1 ;;
            *) if [[ ! -r $f ]]; then echo "(unreadable)" > "$out/$rel"
               elif (( $(stat -c %s "$f") > 65536 )); then
                   echo "sha256 $(sha256sum < "$f" | cut -d' ' -f1) size $(stat -c %s "$f")" > "$out/$rel"
               else cp "$f" "$out/$rel"; fi ;;
        esac
    done
}

# run_case BIN SPEC ROOT: the case in ROOT/NAME/cwd.live (its cwd), run
# N's stdout, stderr, status and snapshot in ROOT/NAME/run.N. No core
# files; the shell's "Segmentation fault" notice goes nowhere.
run_case() {
    local -; set -f
    local bin=$1 spec=$2 root=$3
    local name=${spec%%|*} rest=${spec#*|}
    local seed=${rest%%|*}; rest=${rest#*|}
    local files=${rest%%|*} runs=${rest#*|}
    local cdir=$root/$name dir=$root/$name/cwd.live
    mkdir -p "$dir"
    [[ $seed != - ]] && cp "$repo/$seed" "$dir/hashtable.txt"
    local f fsrc fname fro ros=()
    if [[ -n $files ]]; then
        IFS=';' read -r -a fl <<< "$files"
        for f in "${fl[@]}"; do
            read -r fsrc fname fro <<< "$f"
            if [[ $fsrc == dir ]]; then mkdir -p "$dir/$fname"
            elif [[ $fsrc == - ]]; then : > "$dir/$fname"
            elif [[ -d $dir/$fsrc ]]; then :   # a directory made above, ro below
            else cp "$repo/$fsrc" "$dir/$fname"; fi
            [[ $fro == ro ]] && ros+=("$dir/$fname")
        done
    fi
    for f in ${ros[@]+"${ros[@]}"}; do chmod a-w "$f"; done
    local i=0 line tok argv input lines
    IFS=';' read -r -a lines <<< "$runs"
    (( ${#lines[@]} )) || lines=("")
    for line in "${lines[@]}"; do
        i=$((i+1))
        argv=(); input=/dev/null
        for tok in $line; do
            case $tok in
                "<"*) input=${tok#<}; input=${input//@W/$work/streams} ;;
                "''") argv+=("") ;;
                recordings/*|tests/*) argv+=("$repo/$tok") ;;
                *) argv+=("$tok") ;;
            esac
        done
        local rd=$cdir/run.$i
        mkdir -p "$rd"
        printf '%s\n' "$line" > "$rd/args"
        { (ulimit -c 0; cd "$dir" && exec -a wsprd "$bin" ${argv[@]+"${argv[@]}"} < "$input" > "$rd/stdout" 2> "$rd/stderr")
          echo $? > "$rd/status"; } 2>/dev/null
        snapshot "$dir" "$rd/files"
    done
    chmod -R u+w "$dir" 2>/dev/null
}

# spread_fix ALL: ALL_WSPR.TXT's lines without their last column, the
# spread (" %6.3f"), to ALL.cut, the values to ALL.spread; a line that
# does not end in one to ALL.bad.
spread_fix() {
    awk -v cut="$1.cut" -v sp="$1.spread" -v bad="$1.bad" '
        BEGIN { printf "" > cut; printf "" > sp; printf "" > bad }
        {
            if (match($0, / +-?[0-9]+\.[0-9][0-9][0-9]$/)) {
                v = substr($0, RSTART) + 0; s = sprintf(" %6.3f", v)
                if (substr($0, length($0) - length(s) + 1) == s) {
                    print substr($0, 1, length($0) - length(s)) > cut
                    print substr(s, 2) > sp
                    next
                }
            }
            print NR ": " $0 > bad; print $0 > cut
        }' "$1"
}

compare_case() {  # NAME ROOT_A ROOT_B LABEL_A LABEL_B [spread]
    local name=$1 a=$2/$1 b=$3/$1 la=$4 lb=$5 mode=${6:-}
    local res=$work/results/${name//\//_}
    local d
    if [[ $mode == spread ]]; then
        local fails=0 r
        for r in "$a"/run.* "$b"/run.*; do
            [[ -f $r/files/ALL_WSPR.TXT ]] && spread_fix "$r/files/ALL_WSPR.TXT"
        done
        for r in "$a"/run.*; do
            if [[ -s $r/files/ALL_WSPR.TXT.bad ]]; then
                echo "$la $(basename "$r") ALL_WSPR.TXT lines without a \" %6.3f\" spread:"
                cat "$r/files/ALL_WSPR.TXT.bad"; fails=1
            fi
        done
        d=$(diff -r -x cwd.live -x ALL_WSPR.TXT -x ALL_WSPR.TXT.spread -x ALL_WSPR.TXT.bad \
                "$a" "$b" 2>&1 | head -c 20000) || fails=1
        [[ -n $d ]] && fails=1
        {
            [[ -n $d ]] && echo "$d"
            for r in "$a"/run.*; do
                local rb=$b/$(basename "$r")
                [[ -f $r/files/ALL_WSPR.TXT.spread ]] || continue
                echo "  spread ($la | $lb), $(basename "$r"):"
                paste -d'|' "$r/files/ALL_WSPR.TXT.spread" "$rb/files/ALL_WSPR.TXT.spread" 2>/dev/null \
                    | sed 's/^/    /'
            done
        } > "$res.txt"
        (( fails )) && echo FAIL > "$res.status" || echo PASS > "$res.status"
        return
    fi
    if (( with_n )) && [[ $3 == "$work/stock" ]]; then
        # NEW with -n: STOCK's runs, their stdout with the -n usage line after
        # the -z line, in expected/ (what NEW must match)
        local e=$work/expected/$1 r
        rm -rf "$e"; mkdir -p "$e"
        for r in "$b"/run.*; do
            [[ -d $r ]] || continue
            cp -r "$r" "$e/"
            sed '/^       -z x (x is fano metric table bias, default is 0\.45)$/a\
       -n no drift, does not attempt to follow a drifting frequency' \
                "$r/stdout" > "$e/$(basename "$r")/stdout"
        done
        b=$e
    fi
    if d=$(diff -r -x cwd.live "$a" "$b" 2>&1); then echo PASS > "$res.status"; : > "$res.txt"
    else echo FAIL > "$res.status"; printf '%s\n' "$d" | head -c 20000 > "$res.txt"; fi
}

skip_missing() {  # SPEC -> 0 when an input under recordings/ is missing
    local -; set -f
    local tok
    for tok in ${1//[|;]/ }; do
        [[ $tok == recordings/* && ! -e $repo/$tok ]] && return 0
        [[ $tok == "<@W/"* && ! -e $work/streams/${tok#<@W/} ]] && return 0
    done
    return 1
}

running=0
spawn() { if (( running >= jobs_max )); then wait -n; running=$((running-1)); fi; "$@" & running=$((running+1)); }
for spec in "${cases[@]}"; do
    name=${spec%%|*}
    [[ -n $only && ! $name =~ $only ]] && continue
    skip_missing "$spec" && continue
    spawn run_case "$new" "$spec" "$work/new"
    spawn run_case "$stock" "$spec" "$work/stock"
done
if [[ -n $spread ]]; then
    mkdir -p "$work/spread/new" "$work/spread/other"
    for spec in "${spread_cases[@]}"; do
        name=${spec%%|*}
        [[ -n $only && ! $name =~ $only ]] && continue
        skip_missing "$spec" && continue
        spawn run_case "$new" "$spec" "$work/spread/new"
        spawn run_case "$spread" "$spec" "$work/spread/other"
    done
fi
wait

pass=0 fail=0 skip=0
report() {  # NAME
    local res=$work/results/${1//\//_}
    if [[ $(cat "$res.status") == PASS ]]; then
        echo "PASS $1"; pass=$((pass+1))
        [[ $1 == spread/* ]] && cat "$res.txt"
    else
        echo "FAIL $1"; sed 's/^/    /' "$res.txt"; fail=$((fail+1))
    fi
}
for spec in "${cases[@]}"; do
    name=${spec%%|*}
    [[ -n $only && ! $name =~ $only ]] && continue
    if skip_missing "$spec"; then echo "SKIP $name (input missing)"; skip=$((skip+1)); continue; fi
    compare_case "$name" "$work/new" "$work/stock" new stock
    report "$name"
done
for spec in ${spread_cases[@]+"${spread_cases[@]}"}; do
    name=${spec%%|*}
    [[ -n $only && ! $name =~ $only ]] && continue
    if skip_missing "$spec"; then echo "SKIP $name (input missing)"; skip=$((skip+1)); continue; fi
    compare_case "$name" "$work/spread/new" "$work/spread/other" new other spread
    report "$name"
done
echo "wsprd-compare: $pass passed, $fail failed, $skip skipped ($work)"
(( fail == 0 ))
