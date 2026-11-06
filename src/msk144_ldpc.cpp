// MSK144's (128,90) LDPC decoding and CRC-13 (msk144_ldpc.h): WSJT-X's
// bpdecode128_90.f90, chkcrc13a.f90 (lib/ft8/), crc13.cpp and
// lib/platanh.f90 at 567ad29ce, line for line. Arrays keep the Fortran's
// 1-based indices (element 0 unused) where that keeps the two readable side
// by side.
//
// Steven Franke, K9AN, and Joe Taylor, K1JT, WSJT-X (GNU GPL v3), and so is
// this.
#include "msk144_ldpc.h"

#include <cmath>
#include <cstring>

#include "ft8_ldpc.h"  // ft8::tanhf_glibc

namespace msk144 {

Bpdecode128_90Trace bpdecode128_90_trace = nullptr;

namespace {

#include "msk144_generator.inc"
#include "msk144_ldpc_tables.inc"

constexpr int N = 128, K = 90, M = N - K;

// crc13(data, length) (lib/crc13.cpp): boost::augmented_crc<13, 0x15D7>,
// initial remainder 0, bytes most significant bit first.
int crc13(const unsigned char* data, int length) {
    unsigned rem = 0;
    for (int i = 0; i < length; ++i)
        for (int b = 7; b >= 0; --b) {
            const unsigned top = (rem >> 12) & 1u;
            rem = ((rem << 1) | ((data[i] >> b) & 1u)) & 0x1FFFu;
            if (top) rem ^= 0x15D7u;
        }
    return static_cast<short>(rem);
}

}  // namespace

// platanh(x, y) (lib/platanh.f90): a piecewise-linear atanh.
float platanh(float x) {
    int isign = +1;
    float z = x;
    if (x < 0) {
        isign = -1;
        z = std::fabs(x);
    }
    if (z <= 0.664f) return x / 0.83f;
    if (z <= 0.9217f) return static_cast<float>(isign) * (z - 0.4064f) / 0.322f;
    if (z <= 0.9951f) return static_cast<float>(isign) * (z - 0.8378f) / 0.0524f;
    if (z <= 0.9998f) return static_cast<float>(isign) * (z - 0.9914f) / 0.0012f;
    return static_cast<float>(isign) * 7.0f;
}

// chkcrc13a (lib/ft8/chkcrc13a.f90): true when nbadcrc comes out 0.
bool chkcrc13a(const signed char decoded[90]) {
    // read(cbits,'(12b8)') i1Dec8BitBytes: bytes 1..11 are bits 1..88 (byte 12,
    // bits 89-90 of a short record, is zeroed below as byte 11 is).
    unsigned char bytes[12];
    for (int i = 0; i < 11; ++i) {
        unsigned v = 0;
        for (int b = 0; b < 8; ++b) v = v * 2 + static_cast<unsigned>(decoded[i * 8 + b]);
        bytes[i] = static_cast<unsigned char>(v);
    }
    bytes[11] = 0;
    // read(cbits,'(77x,b13)') ncrc13: the received CRC, bits 78-90
    int ncrc13 = 0;
    for (int b = 77; b < 90; ++b) ncrc13 = ncrc13 * 2 + decoded[b];
    // i1Dec8BitBytes(10)=iand(i1Dec8BitBytes(10),128+64+32+16+8); (11:12)=0
    bytes[9] &= static_cast<unsigned char>(128 + 64 + 32 + 16 + 8);
    bytes[10] = 0;
    bytes[11] = 0;
    const int icrc13 = crc13(bytes, 12);  // CRC13 computed from 77 msg bits
    int nbadcrc = 1;
    if (ncrc13 == icrc13) nbadcrc = 0;
    return nbadcrc == 0;
}

// bpdecode128_90 (lib/bpdecode128_90.f90). The Fortran's Nm(i,j), Mn(kk,j),
// nrw(j) are the lambdas below; tov, toc, tanhtoc, zn, llr, cw keep 1-based
// indices.
int bpdecode128_90(const float llr_in[128], int maxiterations, signed char message77[77], signed char cw_out[128],
                   int& nharderror, int& niterations) {
    signed char cw[N + 1];
    signed char decoded[K + 1];
    int synd[M + 1];
    float tov[4][N + 1];
    float toc[12][M + 1];
    float tanhtoc[12][M + 1];
    float zn[N + 1];
    float llr[N + 1];
    float Tmn;
    const int ncw = kNcw;
    auto Nm = [](int i, int j) { return static_cast<int>(kNm[j - 1][i - 1]); };
    auto Mn = [](int kk, int j) { return static_cast<int>(kMn[j - 1][kk - 1]); };
    auto nrw = [](int j) { return static_cast<int>(kNrw[j - 1]); };
    for (int i = 1; i <= N; ++i) llr[i] = llr_in[i - 1];
    // msk144decodeframe passes apmask = 0: every bit takes the llr + tov branch.

    std::memset(decoded, 0, sizeof decoded);
    std::memset(toc, 0, sizeof toc);
    std::memset(tov, 0, sizeof tov);
    std::memset(tanhtoc, 0, sizeof tanhtoc);
    // initialize messages to checks
    for (int j = 1; j <= M; ++j)
        for (int i = 1; i <= nrw(j); ++i) toc[i][j] = llr[Nm(i, j)];

    int ncnt = 0;
    int nclast = 0;

    int iter;
    for (iter = 0; iter <= maxiterations; ++iter) {
        // Update bit log likelihood ratios (tov=0 in iteration 0).
        for (int i = 1; i <= N; ++i) {
            float s = tov[1][i];  // sum(tov(1:ncw,i)), in order
            for (int kk = 2; kk <= ncw; ++kk) s += tov[kk][i];
            zn[i] = llr[i] + s;
        }

        if (bpdecode128_90_trace) bpdecode128_90_trace(iter, zn, tov);  // the piece check's hook, not in the Fortran

        // Check to see if we have a codeword (check before we do any iteration).
        for (int i = 1; i <= N; ++i) cw[i] = zn[i] > 0.f ? 1 : 0;
        int ncheck = 0;
        for (int i = 1; i <= M; ++i) {
            int s = 0;
            for (int k = 1; k <= nrw(i); ++k) s += cw[Nm(k, i)];
            synd[i] = s;
            if (synd[i] % 2 != 0) ++ncheck;
        }
        if (ncheck == 0) {  // we have a codeword
            for (int i = 1; i <= K; ++i) decoded[i] = cw[i];
            if (chkcrc13a(decoded + 1)) {
                for (int i = 1; i <= 77; ++i) message77[i - 1] = decoded[i];
                nharderror = 0;
                for (int i = 1; i <= N; ++i)
                    if (static_cast<float>(2 * cw[i] - 1) * llr[i] < 0.0f) ++nharderror;
                std::memcpy(cw_out, cw + 1, N);
                niterations = iter;
                return nharderror;
            }
        }

        if (iter > 0) {  // this code block implements an early stopping criterion
            const int nd = ncheck - nclast;
            if (nd < 0) {  // # of unsatisfied parity checks decreased
                ncnt = 0;  // reset counter
            } else {
                ++ncnt;
            }
            if (ncnt >= 3 && iter >= 5 && ncheck > 10) {
                nharderror = -1;
                std::memcpy(cw_out, cw + 1, N);
                niterations = iter;
                return nharderror;
            }
        }
        nclast = ncheck;

        // Send messages from bits to check nodes
        for (int j = 1; j <= M; ++j)
            for (int i = 1; i <= nrw(j); ++i) {
                const int ibj = Nm(i, j);
                toc[i][j] = zn[ibj];
                for (int kk = 1; kk <= ncw; ++kk)  // subtract off what the bit had received from the check
                    if (Mn(kk, ibj) == j) toc[i][j] = toc[i][j] - tov[kk][ibj];
            }

        // send messages from check nodes to variable nodes
        for (int i = 1; i <= M; ++i)
            for (int r = 1; r <= 11; ++r) tanhtoc[r][i] = ft8::tanhf_glibc(-toc[r][i] / 2);

        for (int j = 1; j <= N; ++j)
            for (int i = 1; i <= ncw; ++i) {
                const int ichk = Mn(i, j);  // Mn(:,j) are the checks that include bit j
                // product(tanhtoc(1:nrw(ichk),ichk), mask=Nm(1:nrw(ichk),ichk).ne.j)
                Tmn = 1.0f;
                for (int r = 1; r <= nrw(ichk); ++r)
                    if (Nm(r, ichk) != j) Tmn = Tmn * tanhtoc[r][ichk];
                const float y = platanh(-Tmn);
                tov[i][j] = 2 * y;
            }
    }
    nharderror = -1;
    std::memcpy(cw_out, cw + 1, N);
    niterations = iter;
    return nharderror;
}

// encode_128_90 (lib/encode_128_90.f90).
void encode_128_90(const signed char message77[77], signed char codeword[128]) {
    // fill the generator matrix: gen(i, icol) for the 38 rows, 90 columns (once, thread-safe)
    struct Gen {
        signed char g[M][K];
        Gen() {
            std::memset(g, 0, sizeof g);
            for (int i = 0; i < M; ++i) {
                for (int j = 1; j <= 23; ++j) {
                    const char h = kGen[i][j - 1];
                    const int istr = h <= '9' ? h - '0' : (h | 0x20) - 'a' + 10;  // read(g(i)(j:j),"(Z1)") istr
                    int ibmax = 4;
                    if (j == 23) ibmax = 2;
                    for (int jj = 1; jj <= ibmax; ++jj) {
                        const int icol = (j - 1) * 4 + jj;
                        if ((istr >> (4 - jj)) & 1) g[i][icol - 1] = 1;  // btest(istr, 4-jj)
                    }
                }
            }
        }
    };
    static const Gen gen_matrix;
    const auto& gen = gen_matrix.g;
    // Add 13 bit CRC to form 90-bit message+CRC13: tmpchar = the 77 bits and '000', read as 10b8
    unsigned char bytes[12] = {};
    for (int i = 0; i < 10; ++i) {
        unsigned v = 0;
        for (int b = 0; b < 8; ++b) {
            const int idx = i * 8 + b;  // 0-based bit of tmpchar(1:80)
            v = v * 2 + (idx < 77 ? static_cast<unsigned>(message77[idx]) : 0u);
        }
        bytes[i] = static_cast<unsigned char>(v);
    }
    const int ncrc13 = crc13(bytes, 12);
    signed char message[K];
    for (int i = 0; i < 77; ++i) message[i] = message77[i];
    for (int b = 0; b < 13; ++b) message[77 + b] = static_cast<signed char>((ncrc13 >> (12 - b)) & 1);  // write(tmpchar(78:90),'(b13)')
    signed char pchecks[M];
    for (int i = 0; i < M; ++i) {
        int nsum = 0;
        for (int j = 0; j < K; ++j) nsum = nsum + message[j] * gen[i][j];
        pchecks[i] = static_cast<signed char>(nsum % 2);
    }
    for (int i = 0; i < K; ++i) codeword[i] = message[i];
    for (int i = 0; i < M; ++i) codeword[K + i] = pchecks[i];
}

}  // namespace msk144
