// ft.mjs (the JavaScript API) given what jt9's file mode reads of each WAV
// and the parameters it decodes with (ft8_jt9.cpp's, ft4_jt9.cpp's and
// ft2_jt9.cpp's params_of, jt9's defaults), one Decoder for all the files as
// one jt9 run: its text must be the native jt9's stdout (less its
// " EOF on input file" lines). wasm/check.sh runs it.
//
//   node wasm/check-api.mjs -8|-5|--ft2 DEPTH NATIVE_STDOUT WAV...
import { readFileSync } from "node:fs";
import { basename } from "node:path";
import createFt from "./dist/ft.mjs";

const [flag, depthArg, nativeOut, ...wavs] = process.argv.slice(2);
const mode = { "-8": "ft8", "-5": "ft4", "--ft2": "ft2" }[flag];
// The samples jt9's file mode decodes: FT8 50 blocks of 3456, FT4 21, FT2 10.
const take = { ft8: 50 * 3456, ft4: 21 * 3456, ft2: 10 * 3456 }[mode];

// jt9's WAV reader's samples, 16-bit or 32-bit float.
function readWav(path) {
  const b = readFileSync(path);
  let o = 12, format = 1, bits = 16;
  while (o + 8 <= b.length) {
    const id = b.toString("ascii", o, o + 4);
    const n = b.readUInt32LE(o + 4);
    if (id === "fmt ") {
      format = b.readUInt16LE(o + 8);
      bits = b.readUInt16LE(o + 22);
    }
    if (id === "data") {
      const end = Math.min(b.length, o + 8 + n);
      const data = b.buffer.slice(b.byteOffset + o + 8, b.byteOffset + end);
      if (format === 3 && bits === 32) return new Float32Array(data.slice(0, data.byteLength & ~3));
      return new Int16Array(data.slice(0, data.byteLength & ~1));
    }
    o += 8 + n + (n & 1);
  }
  throw new Error(`${path}: no data chunk`);
}

// jt9's: the name's last six digits, or _HHMM.
function nutc(path) {
  const stem = basename(path).replace(/\.(wav|WAV)$/, "");
  let m = stem.match(/(\d{6})$/);
  if (m) return Number(m[1]);
  m = stem.match(/_(\d{4})$/);
  return m ? Number(m[1]) : 0;
}

const options = { depth: Number(depthArg), freqLow: 200, freqHigh: 4000, rxFreq: 1500, myCall: "K1ABC", hisCall: "W9XYZ" };
if (mode === "ft8") Object.assign(options, { txFreq: 0, apWidth: 75, hisGrid: "EN37" });

const ft = await createFt();
const dec = new ft.Decoder(mode);
let text = "";
for (const w of wavs) text += dec.decode(readWav(w).subarray(0, take), nutc(w), options).text;
dec.delete();

const want = readFileSync(nativeOut, "utf8").split("\n").filter((l) => !l.startsWith(" EOF on input file")).join("\n");
if (text !== want) {
  const a = want.split("\n"), b = text.split("\n");
  const i = a.findIndex((l, k) => l !== b[k]);
  console.log(`  ft.mjs differs at line ${i + 1}:\n    native: ${a[i]}\n    ft.mjs: ${b[i]}`);
  process.exit(1);
}
