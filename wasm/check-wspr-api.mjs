// wspr.mjs (the JavaScript API) given each WAV's samples and wsprd's options
// as the command line gives them, one WsprDecoder for all the files as one
// data directory's runs: its text must be the native wsprd's stdout for
// the same runs. wasm/check.sh runs it.
//
//   node wasm/check-wspr-api.mjs NATIVE_STDOUT "WSPRD OPTIONS" WAV...
import { readFileSync } from "node:fs";
import { basename } from "node:path";
import createWspr from "./dist/wsprd/wspr.mjs";

const [nativeOut, optArgs, ...wavs] = process.argv.slice(2);

// wsprd's options as WsprDecoder's.
function options(words) {
  const o = { cycles: 10000, osdDepth: -1 };  // wsprd's own defaults
  for (let i = 0; i < words.length; i++) {
    const w = words[i], v = () => Number(words[++i]);
    if (w === "-C") o.cycles = v();
    else if (w === "-o") o.osdDepth = v();
    else if (w === "-f") o.dialFreq = v();
    else if (w === "-e") o.dialError = v();
    else if (w === "-z") o.bias = v();
    else if (w === "-m") o.mode = v();
    else if (w === "-d") o.deep = true;
    else if (w === "-q") o.quick = true;
    else if (w === "-s") o.subtraction = false;
    else if (w === "-B") o.blockDemod = false;
    else if (w === "-J") o.stackDecoder = true;
    else if (w === "-w") o.wide = true;
    else if (w === "-n") o.noDrift = true;
    else if (w === "-H") o.hashtable = false;
    else throw new Error(`option ${w} not mapped`);
  }
  return o;
}

// wsprd's samples of a .wav: 16-bit, or 32-bit float.
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
      const data = b.buffer.slice(b.byteOffset + o + 8, b.byteOffset + Math.min(b.length, o + 8 + n));
      if (format === 3 && bits === 32) return new Float32Array(data.slice(0, data.byteLength & ~3));
      return new Int16Array(data.slice(0, data.byteLength & ~1));
    }
    o += 8 + n + (n & 1);
  }
  throw new Error(`${path}: no data chunk`);
}

const base = options(optArgs.split(/\s+/).filter(Boolean));
const wspr = await createWspr();
const dec = new wspr.WsprDecoder();
let text = "";
for (const w of wavs) {
  // wsprd's: the 6 characters 11 before ".wav" and the 4 before it.
  const stem = basename(w).replace(/\.wav$/, "");
  text += dec.decode(readWav(w), { ...base, date: stem.slice(-11, -5), time: stem.slice(-4) }).text;
}
dec.delete();

const want = readFileSync(nativeOut, "utf8");
if (text !== want) {
  const a = want.split("\n"), b = text.split("\n");
  const i = a.findIndex((l, k) => l !== b[k]);
  console.log(`  wspr.mjs differs at line ${i + 1}:\n    native:   ${a[i]}\n    wspr.mjs: ${b[i]}`);
  process.exit(1);
}
