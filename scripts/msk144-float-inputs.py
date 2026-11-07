#!/usr/bin/env python3
"""The two float WAVs of the MSK144 corpus's float-* cases, made from WSJT-X's
sample recordings/msk144-samples/181211_120500.wav (K1JT WA4CQG EM72 at 8.7 s):

  level-181211_120500.wav   every sample s as s/32768 * 2**-20: an rms far below
                            one 16-bit count, as radiod's F32 can sit
  nan-181211_120500.wav     the samples s/32768, with NaN, +inf and -inf written
                            over some of them: samples 10-19 (NaN), 5000 (+inf)
                            and 7000 (-inf), in the first 0.6 s (the first block,
                            where a NaN would set the noise level for the whole
                            file), then 2 s (50 NaN), 3 s (+inf), 4 s (-inf) and
                            13 s (1000 NaN), all away from the decode at 8.7 s

(both in recordings/msk144-float/).  Both are 12 kHz mono WAVE_FORMAT_IEEE_FLOAT under a 44-byte header (as
scripts/float-input-check.sh writes them). Run from the repository's root:

    python3 -I scripts/msk144-float-inputs.py
"""
import math, os, struct, wave

SRC = "recordings/msk144-samples/181211_120500.wav"
OUT = "recordings/msk144-float"


def write(path, rate, xs):
    data = struct.pack(f"<{len(xs)}f", *xs)
    with open(path, "wb") as f:
        f.write(b"RIFF" + struct.pack("<I", 36 + len(data)) + b"WAVE" + b"fmt " +
                struct.pack("<IHHIIHH", 16, 3, 1, rate, rate * 4, 4, 32) + b"data" + struct.pack("<I", len(data)) + data)


with wave.open(SRC, "rb") as w:
    assert w.getnchannels() == 1 and w.getsampwidth() == 2
    rate, raw = w.getframerate(), w.readframes(w.getnframes())
n = len(raw) // 2
ints = struct.unpack(f"<{n}h", raw)
os.makedirs(OUT, exist_ok=True)
write(f"{OUT}/level-181211_120500.wav", rate, [s / 32768.0 * 2.0 ** -20 for s in ints])
x = [s / 32768.0 for s in ints]
for i in range(10, 20): x[i] = math.nan
x[5000] = math.inf
x[7000] = -math.inf
for i in range(2 * rate, 2 * rate + 50): x[i] = math.nan
x[3 * rate] = math.inf
x[4 * rate] = -math.inf
for i in range(13 * rate, 13 * rate + 1000): x[i] = math.nan
write(f"{OUT}/nan-181211_120500.wav", rate, x)
