// fine-freq-wspr-gen: a synthetic WSPR-2 WAV at an exact frequency, double precision.
// usage: fine-freq-wspr-gen F_CENTRE_AUDIO_HZ SNR_DB SEED OUT.wav "CALL GRID DBM"
// The centre (tone 1.5) is at F; tones 12000/8192 Hz apart, phase continuous,
// starting 1 s into a 120 s file (DT 0). Noise: unit variance at 12 kHz, the
// signal's power set for SNR in 2500 Hz; all times 100 into 16 bits.
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <random>
#include <vector>
#include "wsprd_internal.h"
int main(int argc, char** argv) {
    if (argc < 6) { fprintf(stderr, "usage\n"); return 2; }
    double f = atof(argv[1]), snr = atof(argv[2]); unsigned seed = strtoul(argv[3], 0, 10);
    std::vector<char> hashtab(32768 * 13, 0), loctab(32768 * 5, 0);
    char msg[64]; snprintf(msg, sizeof msg, "%s", argv[5]);
    unsigned char sym[162];
    if (!wsprd::detail::get_wspr_channel_symbols(msg, hashtab.data(), loctab.data(), sym)) { fprintf(stderr, "bad message\n"); return 3; }
    const double fs = 12000, baud = 12000.0 / 8192, twopi = 2 * M_PI;
    const long nmax = 120 * 12000, i0 = 12000;
    std::mt19937_64 rng(seed);
    std::vector<double> x(nmax, 0.0);
    const double A = std::sqrt(2.0 * 2500.0 / 6000.0 * std::pow(10.0, snr / 10.0));
    double phi = std::uniform_real_distribution<double>(0, twopi)(rng);
    for (int i = 0; i < 162; i++)
        for (int m = 0; m < 8192; m++) {
            long k = i0 + (long)i * 8192 + m;
            if (k < nmax) x[k] = A * std::cos(phi);
            phi += twopi * (f + (sym[i] - 1.5) * baud) / fs;
            if (phi > twopi) phi -= twopi;
        }
    std::normal_distribution<double> g(0.0, 1.0);
    std::vector<int16_t> s(nmax);
    for (long i = 0; i < nmax; i++) { double v = 100.0 * (x[i] + g(rng)); v = v > 32767 ? 32767 : v < -32768 ? -32768 : v; s[i] = (int16_t)std::lrint(v); }
    FILE* o = fopen(argv[4], "wb"); uint32_t db = nmax * 2, rs = 36 + db, sr = 12000, br = 24000, fl = 16; uint16_t fmt = 1, ch = 1, ba = 2, bps = 16;
    fwrite("RIFF", 1, 4, o); fwrite(&rs, 4, 1, o); fwrite("WAVEfmt ", 1, 8, o); fwrite(&fl, 4, 1, o); fwrite(&fmt, 2, 1, o); fwrite(&ch, 2, 1, o);
    fwrite(&sr, 4, 1, o); fwrite(&br, 4, 1, o); fwrite(&ba, 2, 1, o); fwrite(&bps, 2, 1, o); fwrite("data", 1, 4, o); fwrite(&db, 4, 1, o);
    fwrite(s.data(), 2, nmax, o); fclose(o);
    return 0;
}
