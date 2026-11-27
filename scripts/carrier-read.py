#!/usr/bin/env python3
"""Read a carrier IQ .wav (wsprd -I, jt9 --carrier-dir; src/carrier_iq.h): its JSON and its frequency over time.

  carrier-read.py FILE.wav [--window S] [--step S]

Prints the JSON, then one line a window: the window's centre (s from the input file's start), the carrier's
frequency (audio Hz: ref_hz + the baseband frequency, the window's periodogram peak) and its
SNR in the window. With --mean only the whole transmission's frequency (its periodogram peak).
"""
import json, struct, sys
import numpy as np


def read(path):
    b = open(path, "rb").read()
    assert b[:4] == b"RIFF" and b[8:12] == b"WAVE", "not a WAV"
    i, meta, rate, iq = 12, {}, None, None
    while i + 8 <= len(b):
        cid, n = b[i:i + 4], struct.unpack("<I", b[i + 4:i + 8])[0]
        body = b[i + 8:i + 8 + n]
        if cid == b"fmt ":
            fmt, ch, rate = struct.unpack("<HHI", body[:8])
            assert fmt == 3 and ch == 2, "not 2-channel float"
        elif cid == b"LIST" and body[:4] == b"INFO":
            j = 4
            while j + 8 <= len(body):
                sid, m = body[j:j + 4], struct.unpack("<I", body[j + 4:j + 8])[0]
                if sid == b"ICMT":
                    meta = json.loads(body[j + 8:j + 8 + m].split(b"\0")[0])
                j += 8 + m + (m & 1)
        elif cid == b"data":
            a = np.frombuffer(body, dtype="<f4")
            iq = a[0::2].astype(np.float64) + 1j * a[1::2].astype(np.float64)
        i += 8 + n + (n & 1)
    return meta, rate, iq


def freq_fit(z, rate):
    """Frequency (Hz) of the carrier in z: the peak of its periodogram (zero-padded 64 times, then a parabola through
    the top three points), the maximum-likelihood estimate of a single tone's frequency in noise."""
    n = 1 << int(np.ceil(np.log2(len(z) * 64)))
    S = np.abs(np.fft.fft(z, n)) ** 2
    k = int(np.argmax(S))
    a, b, c = S[k - 1], S[k], S[(k + 1) % n]
    d = 0.5 * (a - c) / (a - 2 * b + c) if (a - 2 * b + c) != 0 else 0.0
    f = (k + d) * rate / n
    return f - rate if f >= rate / 2 else f


def snr_db(z, rate):
    """Carrier power over the noise in 2500 Hz (as WSPR reports): the carrier's power within 0.5 Hz of the strongest
    bin, the noise's density the mean of the bins 3 Hz or more from it."""
    n = len(z)
    S = np.abs(np.fft.fft(z * np.hanning(n))) ** 2
    f = np.fft.fftfreq(n, 1 / rate)
    k = np.argmax(S)
    near = np.abs(f - f[k]) < 0.5
    far = np.abs(f - f[k]) > 3.0
    noise_per_bin = np.mean(S[far])
    sig = S[near].sum() - noise_per_bin * near.sum()
    bin_hz = rate / n
    return 10 * np.log10(max(sig, 1e-30) / (noise_per_bin / bin_hz * 2500))


def main():
    args = sys.argv[1:]
    win = float(args[args.index("--window") + 1]) if "--window" in args else 10.0
    step = float(args[args.index("--step") + 1]) if "--step" in args else 5.0
    meta, rate, iq = read(args[0])
    print(json.dumps(meta))
    ref, t0 = meta["ref_hz"], meta["t0_s"]
    if "--mean" in args:
        print(f"mean {ref + freq_fit(iq, rate):.6f} Hz  snr {snr_db(iq, rate):.1f} dB  n={len(iq)}")
        return
    nw, ns = int(win * rate), int(step * rate)
    for s in range(0, len(iq) - nw + 1, ns):
        z = iq[s:s + nw]
        print(f"{t0 + (s + nw / 2) / rate:8.2f} {ref + freq_fit(z, rate):.6f} {snr_db(z, rate):6.1f}")


if __name__ == "__main__":
    main()
