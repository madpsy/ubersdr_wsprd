// The port's side of scripts/msk144-piece-check.sh frontend: reads the
// oracle's output (scripts/msk144-frontend-oracle.f90) on stdin, takes its A,
// D, B, T and S lines as inputs (and its O arrays as the inputs of tweak1 and
// the frequency search), and prints every line as the reference printed it,
// the results from src/msk_frontend.cpp instead. The two outputs must be
// identical byte for byte.
//
// License: GNU GPL v3 (as WSJT-X)
#include <array>
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

#include "msk_frontend.h"

namespace {

using msk144::cfloat;

uint32_t bits_of(float x) {
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
double double_of(const std::string& hex) {
    const uint64_t b = std::strtoull(hex.c_str(), nullptr, 16);
    double x;
    std::memcpy(&x, &b, 8);
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
    std::vector<float> v(2 * static_cast<size_t>(n));
    for (int j = 0; j < n; ++j) {
        v[2 * j] = x[j].real();
        v[2 * j + 1] = x[j].imag();
    }
    prr(tag, k, v.data(), 2 * n);
}

int kcase = 0;
// analytic_trace of the oracle: H, K and F before the filter, G after it, W after the DC term's half
void trace(int stage, const float* h, const cfloat* corr, const cfloat* c, int nfft) {
    const int nh = nfft / 2;
    if (stage == 0) {
        prr("H", kcase, h, nh + 1);
        prc("K", kcase, corr, nh + 1);
        prc("F", kcase, c, nfft);
    } else if (stage == 1) {
        prc("G", kcase, c, nh + 1);
    } else {
        prc("W", kcase, c, 2);
    }
}

// a line's values, from its third field on (TAG k i v...): stored at i-1...
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

}  // namespace

int main() {
    constexpr int NSPM = 864;
    std::string line;
    int kmax = 0;  // the last case number: the oracle's END count
    std::map<int, std::vector<float>> din, oraclo;  // case -> the input d, the oracle's O (re, im pairs)
    std::vector<float> cbv;
    struct Case {
        int npts = 0, nfft = 0;
        bool beq = false;
        std::array<double, 5> pc{};
    } cur;
    int curk = -1;
    int oskip = -1;  // the case whose H..O lines the port has printed
    int tskip = -1, sskip = -1;
    msk144::analytic_trace = trace;
    while (std::getline(std::cin, line)) {
        std::istringstream is(line);
        std::string tag;
        is >> tag;
        if (tag == "B") {
            std::puts(line.c_str());
            store(line, cbv);
        } else if (tag == "A") {
            std::puts(line.c_str());
            int k, beq;
            std::string h[5];
            is >> k >> cur.npts >> cur.nfft >> beq >> h[0] >> h[1] >> h[2] >> h[3] >> h[4];
            cur.beq = beq != 0;
            for (int i = 0; i < 5; ++i) cur.pc[i] = double_of(h[i]);
            curk = k;
            kmax = std::max(kmax, k);
        } else if (tag == "D") {
            std::puts(line.c_str());
            store(line, din[curk]);
        } else if (tag == "H" || tag == "K" || tag == "F" || tag == "G" || tag == "W" || tag == "O") {
            int k;
            is >> k;
            if (tag == "O") store(line, oraclo[k]);
            if (oskip == k) continue;
            // the first of the case's lines: run analytic, which prints H, K, F and G, then O
            oskip = k;
            kcase = k;
            std::vector<float> d = din[k];
            d.resize(static_cast<size_t>(std::max(cur.nfft, cur.npts)), 0.0f);
            std::vector<cfloat> c(static_cast<size_t>(cur.nfft));
            std::array<double, 5> pc = cur.pc;
            msk144::analytic(d.data(), cur.npts, cur.nfft, c.data(), pc, cur.beq);
            prc("O", k, c.data(), cur.nfft);
        } else if (tag == "T") {
            std::puts(line.c_str());
            int k, src, jz;
            std::string f0hex;
            is >> k >> src >> jz >> f0hex;
            kmax = std::max(kmax, k);
            const std::vector<float>& o = oraclo[src];
            std::vector<cfloat> ca(static_cast<size_t>(jz)), cb(static_cast<size_t>(jz));
            for (int i = 0; i < jz; ++i) ca[i] = cfloat(o[2 * i], o[2 * i + 1]);
            msk144::tweak1(ca.data(), jz, float_of(f0hex), cb.data());
            prc("U", k, cb.data(), jz);
            tskip = k;
        } else if (tag == "U") {
            continue;
        } else if (tag == "S") {
            std::puts(line.c_str());
            int k, src, nfr, mask[8], if1, if2;
            std::string fchex, delfhex;
            is >> k >> src >> nfr;
            kmax = std::max(kmax, k);
            for (int i = 0; i < 8; ++i) is >> mask[i];
            is >> fchex >> if1 >> if2 >> delfhex;
            const float fc = float_of(fchex), delf = float_of(delfhex);
            const std::vector<float>& o = oraclo[src];
            const int n = NSPM * nfr;
            std::vector<cfloat> cdat(static_cast<size_t>(n)), cdat2(static_cast<size_t>(n), cfloat(0.0f, 0.0f)), cd1(static_cast<size_t>(n));
            for (int i = 0; i < n; ++i) cdat[i] = cfloat(o[2 * i], o[2 * i + 1]);
            std::vector<cfloat> cb(42);
            for (int i = 0; i < 42; ++i) cb[i] = cfloat(cbv[2 * i], cbv[2 * i + 1]);
            std::vector<cfloat> cs(NSPM);
            std::vector<float> xccs(NSPM);
            float xmax, bestf;
            msk144::msk144_freq_search(cdat.data(), fc, if1, if2, delf, nfr, mask, cb.data(), cdat2.data(), xmax, bestf, cs.data(),
                                       xccs.data());
            std::printf("X %d %08X %08X\n", k, bits_of(xmax), bits_of(bestf));
            prc("Y", k, cs.data(), NSPM);
            prr("Z", k, xccs.data(), NSPM);
            prc("V", k, cdat2.data(), n);
            for (int ifr = if1; ifr <= if2; ++ifr) {
                std::vector<cfloat> cs1(NSPM);
                float x1, b1;
                msk144::msk144_freq_search(cdat.data(), fc, ifr, ifr, delf, nfr, mask, cb.data(), cd1.data(), x1, b1, cs1.data(),
                                           xccs.data());
                std::printf("Q %d %d %08X %08X\n", k, ifr, bits_of(x1), bits_of(b1));
            }
            sskip = k;
        } else if (tag == "X" || tag == "Y" || tag == "Z" || tag == "V" || tag == "Q" || tag == "END") {
            continue;
        } else {
            std::fprintf(stderr, "msk144-frontend-check: unknown line: %s\n", line.c_str());
            return 1;
        }
    }
    (void)tskip;
    (void)sskip;
    std::printf("END frontend %d\n", kmax);
    return 0;
}
