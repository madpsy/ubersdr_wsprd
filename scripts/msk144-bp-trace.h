// The port's side of scripts/msk144-bp-trace.f90: the hooks of
// src/msk144_ldpc.h printing bpdecode128_90's iterations as the oracle's
// bpdecode128_90_trace does (the I, J, BZ, BV, BO and BH lines), for
// scripts/msk144-ldpc-check.cpp and scripts/msk144-frame-check.cpp.
//
// License: GNU GPL v3 (as WSJT-X)
#pragma once

#include <cstdint>
#include <cstdio>
#include <cstring>

#include "msk144_ldpc.h"

namespace bp_trace {

inline uint32_t bits(float x) {
    uint32_t b;
    std::memcpy(&b, &x, 4);
    return b;
}

// msk144_trace_print: the same sums in the same order, then every value.
inline void iter_lines(int iter, const float zn[129], const float tov[4][129]) {
    float s1 = 0.0f;
    for (int i = 1; i <= 128; ++i) s1 = s1 + zn[i];
    float s2 = 0.0f;
    for (int j = 1; j <= 128; ++j)
        for (int kk = 1; kk <= 3; ++kk) s2 = s2 + tov[kk][j];
    std::printf("I %d %08X %08X %08X %08X %08X %08X %08X\n", iter, bits(s1), bits(s2), bits(zn[1]), bits(zn[65]),
                bits(zn[128]), bits(tov[1][1]), bits(tov[2][64]));
    std::printf("J %08X\n", bits(tov[3][128]));
    std::printf("BZ %d ", iter);
    for (int i = 1; i <= 128; ++i) std::printf("%08X", bits(zn[i]));
    std::printf("\nBV %d ", iter);
    for (int j = 1; j <= 128; ++j)
        for (int kk = 1; kk <= 3; ++kk) std::printf("%08X", bits(tov[kk][j]));
    std::printf("\n");
}

// msk144_trace_msg: toc(1:nrw(j),j) and tanhtoc(1:nrw(j),j), j = 1..38.
inline void msg_lines(int iter, const float toc[12][39], const float tanhtoc[12][39]) {
    static const int nrw[38] = {10, 10, 11, 10, 11, 11, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10,
                                10, 10, 10, 10, 10, 10, 10, 11, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10};
    int n = 0;
    for (int j = 0; j < 38; ++j) n += nrw[j];
    for (int w = 0; w < 2; ++w) {
        const float(*a)[39] = w == 0 ? toc : tanhtoc;
        std::printf("%s %d %d ", w == 0 ? "BO" : "BH", iter, n);
        for (int j = 1; j <= 38; ++j)
            for (int i = 1; i <= nrw[j - 1]; ++i) std::printf("%08X", bits(a[i][j]));
        std::printf("\n");
    }
}

inline void install() {
    msk144::bpdecode128_90_trace = iter_lines;
    msk144::bpdecode128_90_msg_trace = msg_lines;
}
inline void remove() {
    msk144::bpdecode128_90_trace = nullptr;
    msk144::bpdecode128_90_msg_trace = nullptr;
}

}  // namespace bp_trace
