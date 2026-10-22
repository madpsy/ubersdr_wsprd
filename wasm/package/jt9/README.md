# ubersdr-jt9-wasm

FT8, FT4 and FT2 decoders in WebAssembly, decoding as WSJT-X's `jt9` does
(WSJT-X 567ad29ce for FT8 and FT4; WSJT-X improved 3.2.0's logic for FT2).
Runs in browsers, Web Workers and Node 18+.

Against WSJT-X's own outputs, 299 of 305 test cases are identical, byte
for byte. In the other six a few decodes differ: weak a-priori (a7)
decodes at -25 dB, and in one 43-period run three decodes of one crowded
period. WebAssembly has no fused multiply-add and FFTW has no WebAssembly
SIMD, so a few results round differently from x86-64's in their last bits.

## Use

```js
import createFt from "ubersdr-jt9-wasm";

const ft = await createFt();
const dec = new ft.Decoder("ft8");          // "ft8", "ft4" or "ft2"

// One period's 12 kHz mono audio from its start: FT8 15 s (180,000
// samples), FT4 7.5 s, FT2 3.75 s. Float32Array -1..1, or Int16Array.
const r = dec.decode(samples, 193500, { depth: 2 });   // utc as HHMMSS
for (const d of r.decodes) console.log(d.snr, d.dt, d.freq, d.message);

dec.delete();                                // when finished with it
```

Keep one `Decoder` per dial and feed it every period in turn. It carries
WSJT-X's state from period to period (the hash tables for `<...>` calls
and the a7 decodes), as one running `jt9` does.

Options (all optional; the defaults are the WSJT-X GUI's):

| option | default | |
|---|---|---|
| `depth` | 2 | 1-3, as `jt9 -d` |
| `freqLow`, `freqHigh` | 100, 3200 | search range, Hz |
| `rxFreq`, `txFreq` | 1500 | Hz |
| `apWidth` | 50 | a-priori search width, Hz |
| `myCall`, `hisCall`, `hisGrid` | blank | with blank calls only the CQ a-priori passes run |
| `qsoProgress` | 0 | |
| `contest` | 0 | WSJT-X's special operating activity (low 3 bits) |
| `apDecode` | true | |
| `apCqOnly` | false | |
| `fastSubtract` | false | quicker; the same decodes, not WSJT-X's to the bit |

Each decode has `snr`, `dt` (s), `freq` (Hz), `message`, `lowConfidence`
(jt9's `?`), `ap` (a-priori type, 0 for none), `quality`, `sync` and `line`
(the line `jt9` prints). `r.text` is the period's whole `jt9` output.

A decode takes time: about 1-2 s for a busy FT8 period at depth 2. Run it
in a Web Worker in a browser.

In a browser bundler, `ft.wasm` is fetched next to `ft.mjs`. To serve it
from somewhere else, pass `locateFile`:
`createFt({ locateFile: (f) => "/static/" + f })`.

## jt9 on Node

`jt9.cjs` is the `jt9` command line (FT8 `-8`, FT4 `-5`, FT2 `--ft2`,
FST4W `-W`) on the real file system, for Node:

    npx jt9-wasm -8 -d 2 261012_193500.wav

It writes `decoded.txt` and `jt9_wisdom.dat` as `jt9` does. `-s` (shared
memory) is not available.

## Licence

GPL-3.0-or-later, as WSJT-X. Joe Taylor, K1JT, Steven Franke, K9AN, and
the WSJT-X developers.
