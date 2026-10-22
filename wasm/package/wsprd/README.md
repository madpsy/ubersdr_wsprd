# ubersdr-wsprd-wasm

WSJT-X's WSPR decoder, `wsprd` (WSJT-X 567ad29ce's), in WebAssembly, for
browsers, Web Workers and Node 18+. WSPR-2 and WSPR-15.

Its decodes are WSJT-X's `wsprd`'s, byte for byte: on every one of the
test cases WSJT-X's `wsprd` was run on, every decode and every line it
prints and writes is the same. Two things differ: where WSJT-X's `wsprd`
crashes (a read-only data directory, `-J` with a large `-z`), what it had
printed is not lost; and the `.c2` file `-c` writes from a `.wav` differs
in its last bits (FFTW has no WebAssembly SIMD).

## Use

```js
import createWspr from "ubersdr-wsprd-wasm";

const wspr = await createWspr();
const dec = new wspr.WsprDecoder();

// One period's 12 kHz mono audio from its start (114 s are read):
// Int16Array, or Float32Array -1..1.
const r = dec.decode(samples, { dialFreq: 14.0956, date: "261022", time: "1200" });
for (const d of r.decodes) console.log(d.snr, d.dt, d.freq, d.drift, d.callsign, d.grid, d.power);

dec.delete();                                // when finished with it
```

Keep one `WsprDecoder` per dial and give it every period in turn. It keeps
wsprd's hash table (`hashtable.txt`: the calls behind `<...>` in type 2 and
3 messages) from period to period. `dec.hashtable()` returns it as wsprd
writes it; `dec.loadHashtable(text)` reads one in.

Options (all optional; the defaults are wsprd's as WSJT-X runs it, `-C 500
-o 4`):

| option | default | wsprd |
|---|---|---|
| `date`, `time` | `000000`, `0000` | YYMMDD and HHMM, from a .wav's name |
| `dialFreq` | 0 | `-f`, MHz |
| `mode` | 2 | `-m` (15: WSPR-15, 912 s read) |
| `cycles` | 500 | `-C` |
| `osdDepth` | 4 | `-o` (-1: none) |
| `deep` | false | `-d` |
| `quick` | false | `-q` |
| `subtraction` | true | false: `-s` |
| `blockDemod` | true | false: `-B` |
| `stackDecoder` | false | `-J` |
| `wide` | false | `-w` |
| `bias` | 0.45 | `-z` |
| `dialError` | 0 | `-e`, Hz |
| `noDrift` | false | `-n` |
| `hashtable` | true | false: `-H` |

Each decode has `snr`, `dt` (s), `freq` (MHz), `drift`, `message`,
`callsign`, `grid`, `power` (dBm), `sync`, `spread` (Hz), `osd`, `pass`,
and wsprd's lines: `line` (stdout), `allWsprLine` (ALL_WSPR.TXT),
`spotsLine` (wspr_spots.txt).

A WSPR-2 period takes under a second to decode (0.66 s under Node on
x86-64, `-C 500 -o 4`; `deep` and a higher `osdDepth` take longer). Run it
in a Web Worker in a browser.

## wsprd on Node

`wsprd.cjs` is `wsprd`'s command line on the real file system, for Node:

    npx wsprd-wasm -C 500 -o 4 -f 14.0956 150426_0918.wav

It writes `ALL_WSPR.TXT`, `wspr_spots.txt`, `hashtable.txt` and the rest as
`wsprd` does.

## Licence

GPL-3.0-or-later, as WSJT-X. `wsprd.c`: Joe Taylor, K1JT, and Steven
Franke, K9AN; `fano.c`: Phil Karn, KA9Q; `osdwspr.f90`: Steven Franke and
Joe Taylor.
