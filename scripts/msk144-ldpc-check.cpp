// The port's side of scripts/msk144-piece-check.sh ldpc: reads the oracle's
// output (scripts/msk144-ldpc-oracle.f90) on stdin, takes its L, C, P and T
// lines as inputs, and prints every line as the reference printed it, the
// results from src/msk144_ldpc.cpp (and ft8::tanhf_glibc) instead; the oracle's
// I, J and R lines are not read, the port's trace hook and decoder make them.
// The two outputs must be identical byte for byte.
//
// License: GNU GPL v3 (as WSJT-X)
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <string>

#include "ft8_ldpc.h"
#include "msk144_ldpc.h"

namespace {

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

// msk144_trace_print in the oracle: the same sums in the same order.
void trace(int iter, const float zn[129], const float tov[4][129]) {
    float s1 = 0.0f;
    for (int i = 1; i <= 128; ++i) s1 = s1 + zn[i];
    float s2 = 0.0f;
    for (int j = 1; j <= 128; ++j)
        for (int kk = 1; kk <= 3; ++kk) s2 = s2 + tov[kk][j];
    std::printf("I %d %08X %08X %08X %08X %08X %08X %08X\n", iter, bits_of(s1), bits_of(s2), bits_of(zn[1]),
                bits_of(zn[65]), bits_of(zn[128]), bits_of(tov[1][1]), bits_of(tov[2][64]));
    std::printf("J %08X\n", bits_of(tov[3][128]));
}

}  // namespace

int main() {
    std::string line;
    int nrun = 0;  // decodes run: the oracle's END count
    float llr[128] = {};
    signed char m77[77] = {};
    int maxiter = 10, want_trace = 0;
    while (std::getline(std::cin, line)) {
        if (line.rfind("L ", 0) == 0) {
            // L MAXITER TRACE HEX*128
            char* end;
            maxiter = static_cast<int>(std::strtol(line.c_str() + 2, &end, 10));
            want_trace = static_cast<int>(std::strtol(end, &end, 10));
            const std::string hex(end + 1);
            for (int i = 0; i < 128; ++i) llr[i] = float_of(hex.substr(8 * i, 8));
            std::puts(line.c_str());
        } else if (line.rfind("R ", 0) == 0) {
            int nharderror = -99, iter = -99;
            signed char cw[128] = {};
            msk144::bpdecode128_90_trace = want_trace ? trace : nullptr;
            msk144::bpdecode128_90(llr, maxiter, m77, cw, nharderror, iter);
            msk144::bpdecode128_90_trace = nullptr;
            signed char word[90];
            std::memcpy(word, cw, 90);
            const int nbadcrc = msk144::chkcrc13a(word) ? 0 : 1;
            std::string m, c;
            for (int i = 0; i < 77; ++i) m += static_cast<char>('0' + m77[i]);
            for (int i = 0; i < 128; ++i) c += static_cast<char>('0' + cw[i]);
            std::printf("R %d %d %d %s %s %d\n", maxiter, iter, nharderror, m.c_str(), c.c_str(), nbadcrc);
            ++nrun;
        } else if (line.rfind("C ", 0) == 0) {
            signed char word[90];
            for (int i = 0; i < 90; ++i) word[i] = static_cast<signed char>(line[2 + i] - '0');
            std::printf("C %s %d\n", line.substr(2, 90).c_str(), msk144::chkcrc13a(word) ? 0 : 1);
        } else if (line.rfind("P ", 0) == 0) {
            const float x = float_of(line.substr(2, 8));
            std::printf("P %08X %08X\n", bits_of(x), bits_of(msk144::platanh(x)));
        } else if (line.rfind("T ", 0) == 0) {
            const float x = float_of(line.substr(2, 8));
            std::printf("T %08X %08X\n", bits_of(x), bits_of(ft8::tanhf_glibc(x)));
        }
    }
    std::printf("END ldpc %d\n", nrun);
    return 0;
}
