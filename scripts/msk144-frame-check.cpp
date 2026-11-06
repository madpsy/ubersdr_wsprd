// The port's side of scripts/msk144-piece-check.sh frame: reads the oracle's
// output (scripts/msk144-frame-oracle.f90) on stdin, takes its PK, C, RUN, QS,
// EN, GM, FM, AT and AV lines as inputs, and prints every line as the
// reference printed it, the results from src/msk144_decode.cpp instead (the
// trace lines through the port's hooks). The two outputs must be identical
// byte for byte.
//
// The module packjt77's state is shared by the reference's calls: the driver
// makes the same pack77 and unpack77 calls on one Packjt77 in the order the
// oracle's lines come (PK, then RUN's unpack77, then QS's genmsk_128_90).
//
// License: GNU GPL v3 (as WSJT-X)
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "ft8_maths.h"
#include "msk144_decode.h"
#include "msk144_ldpc.h"

namespace {

using msk144::cfloat;

// the bits of a float, a NaN as the default one: BLIND SPOT: neither side compares the sign and payload of a NaN
// (the oracle's copy and the library are different objects, and x86 returns the first operand's NaN), so a port
// that gets a NaN's sign or payload wrong passes (the atan2 sweeps' hashw uses raw bits)
uint32_t bits_of(float x) {
    if (x != x) return 0x7FC00000u;
    uint32_t b;
    std::memcpy(&b, &x, 4);
    return b;
}
float float_of(const std::string& hex) {
    const uint32_t b = static_cast<uint32_t>(std::strtoul(hex.c_str(), nullptr, 16));
    float x;
    std::memcpy(&x, &b, 4);
    return x;
}

// prr and prc of the oracle: 16 values a line
void prr(const char* tag, int k, const float* x, int n) {
    for (int i = 1; i <= n; i += 16) {
        std::printf("%s %d %d", tag, k, i);
        for (int j = i; j <= std::min(i + 15, n); ++j) std::printf(" %08X", bits_of(x[j - 1]));
        std::printf("\n");
    }
}
void prc(const char* tag, int k, const cfloat* x, int n) {
    std::vector<float> v(2 * static_cast<size_t>(std::max(n, 0)));
    for (int j = 0; j < n; ++j) {
        v[2 * j] = x[j].real();
        v[2 * j + 1] = x[j].imag();
    }
    prr(tag, k, v.data(), 2 * n);
}

int kc = 0;

// the traces, as the oracle's hooks print
void t_consts(const float* pp, const cfloat* cb) {
    prr("KP", kc, pp, 12);
    prc("KB", kc, cb, 42);
}
void t_phase(cfloat cca, cfloat ccb, float phase0) {
    std::printf("KH %d %08X %08X %08X %08X %08X\n", kc, bits_of(cca.real()), bits_of(cca.imag()), bits_of(ccb.real()),
                bits_of(ccb.imag()), bits_of(phase0));
}
void t_corrected(cfloat cfac, const cfloat* c) {
    std::printf("KF %d %08X %08X\n", kc, bits_of(cfac.real()), bits_of(cfac.imag()));
    prc("KC", kc, c, 864);
}
void t_raw(const float* softbits) { prr("KS", kc, softbits, 144); }
void t_sync(int n1, int n2) { std::printf("KY %d %d %d\n", kc, n1, n2); }
void t_norm(float sav, float s2av, float ssig) {
    std::printf("KN %d %08X %08X %08X\n", kc, bits_of(sav), bits_of(s2av), bits_of(ssig));
}
void t_llr(const float* llr) { prr("KL", kc, llr, 128); }
void t_bp(const signed char* decoded77, const signed char* cw, int nharderror, int niterations) {
    std::string m, c;
    for (int i = 0; i < 128; ++i) c += static_cast<char>('0' + cw[i]);
    if (nharderror >= 0) {
        for (int i = 0; i < 77; ++i) m += static_cast<char>('0' + decoded77[i]);
    } else {
        m = "-";  // decoded77 is written only on success
    }
    std::printf("KD %d %d %d %s %s\n", kc, nharderror, niterations, m.c_str(), c.c_str());
}
void t_tag(const char* c77, int n3, int i3) { std::printf("KT %d %.77s %d %d\n", kc, c77, n3, i3); }

// one oracle line's values from its fourth field on: TAG k i v...
void store(const std::string& line, std::vector<float>& into) {
    std::istringstream is(line);
    std::string tag, hex;
    int k, i;
    is >> tag >> k >> i;
    size_t at = static_cast<size_t>(i - 1);
    while (is >> hex) {
        if (into.size() <= at) into.resize(at + 1);
        into[at++] = float_of(hex);
    }
}

// the oracle's second LCG, for the atan2 sweeps
uint64_t state2 = 777;
double nextu2() {
    state2 = (state2 * 1664525ULL + 1013904223ULL) % 4294967296ULL;
    return static_cast<double>(state2) / 4294967296.0;
}
void hashw(uint64_t& h1, uint64_t& h2, float z) {
    uint32_t raw;  // the atan2 sweeps hash every bit, a NaN's included
    std::memcpy(&raw, &z, 4);
    const uint64_t w = raw;
    h1 = (h1 * 1000003ULL + w) % 2147483647ULL;
    h2 = (h2 * 999983ULL + w) % 2147483629ULL;
}

// a 37-character field: the characters after the first '|' at or after `from`
std::string field37(const std::string& line, size_t from, size_t* after = nullptr) {
    const size_t bar = line.find('|', from);
    std::string f = line.substr(bar + 1, 37);
    if (after) *after = bar + 1 + 37 + 1;  // past the closing '|'
    return f;
}

}  // namespace

int main() {
    msk144::Packjt77 pk;
    std::map<int, std::vector<float>> frames;  // case -> cframe, re/im pairs
    std::vector<float> softbits(144, -3.0f);   // the last decoded frame's, as the oracle gives signal quality
    int nframes = 0;
    msk144::decode_trace.consts = t_consts;
    msk144::decode_trace.phase = t_phase;
    msk144::decode_trace.corrected = t_corrected;
    msk144::decode_trace.raw = t_raw;
    msk144::decode_trace.sync = t_sync;
    msk144::decode_trace.norm = t_norm;
    msk144::decode_trace.llr = t_llr;
    msk144::decode_trace.bp = t_bp;
    msk144::decode_trace.tag = t_tag;

    std::string line;
    while (std::getline(std::cin, line)) {
        std::istringstream is(line);
        std::string tag;
        is >> tag;
        if (tag == "PK" || tag == "PN") {
            int it, i3, n3;
            is >> it;
            size_t after;
            const std::string text = field37(line, 0, &after);
            std::istringstream is2(line.substr(after));
            is2 >> i3 >> n3;  // the oracle's, not used: the port makes its own
            const bool nosave = tag == "PN";  // record_tx_hashes false: pack77 (these texts pack, so as the fallback does)
            const msk144::Packjt77::Packed p = nosave ? pk.pack77(text, false) : pk.pack77_legacy_truncating_fallback(text);
            std::printf("%s %d |%s| %d %d %s\n", tag.c_str(), it, text.c_str(), p.i3, p.n3, p.c77.c_str());
        } else if (tag == "MB") {
            std::puts(line.c_str());  // an input, given back
        } else if (tag == "C") {
            std::puts(line.c_str());
            int k;
            is >> k;
            store(line, frames[k]);
        } else if (tag == "RUN") {
            std::puts(line.c_str());
            is >> kc;
            nframes = std::max(nframes, kc);
            const std::vector<float>& f = frames.at(kc);
            std::vector<cfloat> c(864);
            for (int i = 0; i < 864; ++i) c[i] = cfloat(f[2 * i], f[2 * i + 1]);
            std::fill(softbits.begin(), softbits.end(), -3.0f);
            char msg[37];
            std::memset(msg, ' ', 37);
            std::memcpy(msg, "unset", 5);
            int ns = -1;
            msk144::msk144decodeframe(pk, c.data(), softbits.data(), msg, ns);
            std::printf("R %d %d |%.37s|\n", kc, ns, msg);
            prr("S", kc, softbits.data(), 144);
            prc("O", kc, c.data(), 864);
        } else if (tag == "QS") {
            std::puts(line.c_str());
            int k;
            is >> k;
            size_t after;
            const std::string msg = field37(line, 0, &after);
            std::istringstream is2(line.substr(after));
            std::string a, b, d;
            is2 >> a >> b >> d;
            const float snr = float_of(a), freq = float_of(b), t0 = float_of(d);
            const double pcoef[5] = {0, 0, 0, 0, 0};
            int nbit = -99;
            float eye = -7.0f;
            std::vector<cfloat> cframe(864);
            msk144::msk144signalquality(pk, cframe.data(), snr, freq, t0, softbits.data(), msg.c_str(), "K1ABC       ", false, "/tmp",
                                        nbit, eye, pcoef);
            std::printf("Q %d %d %08X\n", k, nbit, bits_of(eye));
        } else if (tag == "EN") {
            int it, mask;
            std::string m;
            is >> it >> m >> mask;
            signed char m77[77], cw[128];
            for (int i = 0; i < 77; ++i) m77[i] = static_cast<signed char>(m[i] - '0');
            msk144::encode_128_90(m77, cw);
            std::string c;
            for (int i = 0; i < 128; ++i) c += static_cast<char>('0' + cw[i]);
            std::printf("EN %d %s %d %s\n", it, m.c_str(), mask, c.c_str());
        } else if (tag == "GM") {
            int it, ichk, itype;
            is >> it;
            size_t after;
            const std::string msg = field37(line, 0, &after);
            char msgsent[37];
            int tones[144];
            std::memset(msgsent, '?', 37);
            for (int& t : tones) t = -77;
            itype = -99;
            msk144::genmsk_128_90(pk, msg.c_str(), 0, msgsent, tones, itype);
            std::istringstream is2(line.substr(after));
            is2 >> ichk;
            std::string ts;
            for (int t : tones) ts += static_cast<char>('0' + std::max(0, std::min(9, t)));
            // (the message may hold a NUL: written as 37 characters)
            std::printf("GM %d |", it);
            std::fwrite(msg.data(), 1, 37, stdout);
            std::printf("| %d %d |", ichk, itype);
            std::fwrite(msgsent, 1, 37, stdout);
            std::printf("| %s\n", ts.c_str());
        } else if (tag == "FM") {
            int it, n;
            is >> it >> n;
            const size_t bar = line.find('|');
            std::string in = line.substr(bar + 1, static_cast<size_t>(n));
            std::string buf = in;
            int iz;
            msk144::fmtmsg(buf, iz);
            if (n > 0) {
                std::printf("FM %d %d |%s| %d |%s|\n", it, n, in.c_str(), iz, buf.c_str());
            } else {
                std::printf("FM %d %d || %d ||\n", it, n, iz);
            }
        } else if (tag == "AT") {
            int which, n;
            is >> which >> n;
            uint64_t h1 = 1, h2 = 1;
            int cnt = 0;
            if (which == 1) {
                static const double sc[14] = {1.0e-7, 1.0e-6, 1.0e-5, 1.0e-4, 1.0e-3, 1.0e-2, 1.0e-1,
                                              1.0e0,  1.0e1,  1.0e2,  1.0e3,  1.0e4,  1.0e5,  1.0e6};
                state2 = 777;
                for (int idx = 1; idx <= 1500000; ++idx) {
                    double u1 = nextu2(), u2 = nextu2();
                    const float y = static_cast<float>((u1 - 0.5) * 2.0 * sc[static_cast<int>(u2 * 14.0)]);
                    u1 = nextu2();
                    u2 = nextu2();
                    const float x = static_cast<float>((u1 - 0.5) * 2.0 * sc[static_cast<int>(u2 * 14.0)]);
                    hashw(h1, h2, ft8::atan2f_glibc(y, x));
                    ++cnt;
                }
            } else if (which == 2) {
                for (int64_t ib = -2147483647LL - 1; ib <= 2147483647LL - 251; ib += 251) {
                    const uint32_t u = static_cast<uint32_t>(static_cast<int32_t>(ib));
                    float y;
                    std::memcpy(&y, &u, 4);
                    hashw(h1, h2, ft8::atan2f_glibc(y, 1.0f));
                    ++cnt;
                }
            } else {
                state2 = 777;
                // the oracle's third sweep follows the first on the same generator: replay the first's draws
                static const double sc[14] = {1.0e-7, 1.0e-6, 1.0e-5, 1.0e-4, 1.0e-3, 1.0e-2, 1.0e-1,
                                              1.0e0,  1.0e1,  1.0e2,  1.0e3,  1.0e4,  1.0e5,  1.0e6};
                (void)sc;
                for (int draws = 0; draws < 1500000 * 4; ++draws) nextu2();
                for (int idx = 1; idx <= 1000000; ++idx) {
                    float z = static_cast<float>(nextu2() * 12.0);
                    if (idx > 500000) z = static_cast<float>(nextu2() * 0.75);
                    float y = z * 1.5f;
                    float x = 1.5f;
                    for (int ii = 1; ii <= 4; ++ii) {
                        if (ii == 2 || ii == 4) y = -y;
                        if (ii >= 3) x = -1.5f;
                        hashw(h1, h2, ft8::atan2f_glibc(y, x));
                        ++cnt;
                    }
                }
            }
            std::printf("AT %d %d %llu %llu\n", which, cnt, static_cast<unsigned long long>(h1), static_cast<unsigned long long>(h2));
        } else if (tag == "AV") {
            int kk;
            std::string yh, xh, zh;
            is >> kk >> yh >> xh >> zh;
            const float y = float_of(yh), x = float_of(xh);
            std::printf("AV %d %08X %08X %08X\n", kk, bits_of(y), bits_of(x), bits_of(ft8::atan2f_glibc(y, x)));
        } else if (tag == "KP" || tag == "KB" || tag == "KH" || tag == "KF" || tag == "KC" || tag == "KS" || tag == "KY" ||
                   tag == "KN" || tag == "KL" || tag == "KD" || tag == "KT" || tag == "R" || tag == "S" || tag == "O" ||
                   tag == "Q" || tag == "END") {
            continue;
        } else {
            std::fprintf(stderr, "msk144-frame-check: unknown line: %s\n", line.c_str());
            return 1;
        }
    }
    std::printf("END frame %d\n", nframes);
    return 0;
}
