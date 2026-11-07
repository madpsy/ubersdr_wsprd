// MSK144's (128,90) LDPC decoding and CRC-13 (msk144_ldpc.h): WSJT-X's
// bpdecode128_90.f90, chkcrc13a.f90 (lib/ft8/), crc13.cpp and
// lib/platanh.f90 at 567ad29ce, line for line (bpdecode128_90_reference);
// bpdecode128_90 is that decoder rearranged for speed, every float it
// computes the same (see there). Arrays keep the Fortran's 1-based indices
// (element 0 unused) where that keeps the two readable side by side.
//
// Steven Franke, K9AN, and Joe Taylor, K1JT, WSJT-X (GNU GPL v3), and so is
// this.
#include "msk144_ldpc.h"

#include <cmath>
#include <cstdint>
#include <cstring>

#include "ft8_ldpc.h"  // ft8::tanhf_glibc, ft8::bp_tanhf

namespace msk144 {

Bpdecode128_90Trace bpdecode128_90_trace = nullptr;
Bpdecode128_90MsgTrace bpdecode128_90_msg_trace = nullptr;

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

// bpdecode128_90 (lib/bpdecode128_90.f90), line for line. This is the routine
// the fast bpdecode128_90 below was made from, and the one it falls back on if
// its tables do not check out. The Fortran's
// Nm(i,j), Mn(kk,j), nrw(j) are the lambdas below; tov, toc, tanhtoc, zn,
// llr, cw keep 1-based indices.
int bpdecode128_90_reference(const float llr_in[128], int maxiterations, signed char message77[77],
                             signed char cw_out[128], int& nharderror, int& niterations) {
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
            float s = 0.f;  // sum(tov(1:ncw,i)): gfortran starts from +0, so -0 terms give +0
            for (int kk = 1; kk <= ncw; ++kk) s += tov[kk][i];
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
        if (bpdecode128_90_msg_trace && iter < maxiterations) bpdecode128_90_msg_trace(iter, toc, tanhtoc);

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

namespace {

// bpdecode128_90's code as edges (BpTables128): the 384 (check, bit) pairs,
// laid out by check, row r (the Fortran's Nm row, 0-based) times kW plus the
// check's column. Columns: the checks of 11 bits first, then those of 10, each
// group in check order; 38 checks padded to 40, a multiple of 8.
constexpr int kW = 40, kRows = 11;
constexpr int kNt = 10 * kW + 8;  // edges taken tanh of: rows 0-9, and row 10's first 8 columns

struct BpTables128 {
    bool ok = false;
    int n11 = 0;                  // checks of 11 bits: columns 0..n11-1 (at most 8: column group 0)
    int col_check[kW];            // the check (0-based) in a column; -1 for the padding
    int ebit[kRows * kW];         // the bit (0-based) on an edge; 0 for the unused slots
    static_assert(kNcw == 3, "slot[][3] and the three-term sum in bpdecode128_90 assume ncw == 3");
    int slot[N][3];               // the edge of bit i and its kk-th check Mn(kk+1, i+1)
    BpTables128() {
        int col = 0;
        for (int pass = 11; pass >= 10; --pass)
            for (int j = 0; j < M; ++j)
                if (kNrw[j] == pass) col_check[col++] = j;
        if (col != M) return;  // a check of another size: the reference
        for (int j = 0; j < M; ++j) n11 += kNrw[j] == 11;
        if (n11 > 8) return;
        for (int c = M; c < kW; ++c) col_check[c] = -1;
        std::memset(ebit, 0, sizeof ebit);
        int seen[kRows * kW] = {};
        for (int c = 0; c < M; ++c)
            for (int r = 0; r < kNrw[col_check[c]]; ++r) {
                const int b = kNm[col_check[c]][r];
                if (b < 1 || b > N) return;
                ebit[r * kW + c] = b - 1;
            }
        // Each bit's checks are distinct and each holds the bit once, so the
        // Fortran's toc subtracts exactly one tov and its product leaves out
        // exactly one factor: the edge's own.
        for (int i = 0; i < N; ++i)
            for (int kk = 0; kk < kNcw; ++kk) {
                const int j = kMn[i][kk] - 1;
                if (j < 0 || j >= M) return;
                for (int k2 = 0; k2 < kk; ++k2)
                    if (kMn[i][k2] - 1 == j) return;
                int c = 0;
                while (col_check[c] != j) ++c;
                int found = 0;
                for (int r = 0; r < kNrw[j]; ++r)
                    if (kNm[j][r] == i + 1) {
                        slot[i][kk] = r * kW + c;
                        ++found;
                    }
                if (found != 1) return;
                ++seen[slot[i][kk]];
            }
        int nedge = 0;
        for (int e = 0; e < kRows * kW; ++e) {
            if (seen[e] > 1) return;
            nedge += seen[e];
        }
        int nrws = 0;
        for (int j = 0; j < M; ++j) nrws += kNrw[j];
        ok = kNcw == 3 && nedge == N * kNcw && nedge == nrws;  // every edge some bit's, once
    }
};

const BpTables128& bp_tables128() {
    static const BpTables128 t;
    return t;
}
}  // namespace

bool bpdecode128_90_fast_tables_ok() { return bp_tables128().ok; }

namespace {

// Eight float lanes (GCC and Clang vector extensions: AVX on x86-64-v3, NEON
// pairs on arm64, scalar code elsewhere). Each lane's operations are IEEE
// operations on that lane alone, each rounded on its own (-ffp-contract=off:
// no multiply-add is fused), so a lane's result is the scalar code's.
typedef float v8f __attribute__((vector_size(32)));
typedef int32_t v8i __attribute__((vector_size(32)));
inline v8f splat(float x) { return v8f{x, x, x, x, x, x, x, x}; }
inline v8f vsel(v8i m, v8f a, v8f b) { return (v8f)((m & (v8i)a) | (~m & (v8i)b)); }  // m ? a : b
inline v8f vload(const float* p) {
    v8f v;
    std::memcpy(&v, p, sizeof v);
    return v;
}
inline void vstore(float* p, v8f v) { std::memcpy(p, &v, sizeof v); }

// 2*platanh(-Tmn) in eight lanes, platanh's values to the bit. The scalar
// takes the first of its tests that holds; the thresholds rise, so the first
// that holds is the least, and the selects below, made from the last test
// back, take the same branch. A branch's division is its own operands' (x and
// 0.83, or isign*(z - a) and b, a and b that branch's constants), so dividing
// once, after selecting them, gives what the branch's own division gives: an
// IEEE division depends on nothing but its operands. A lane no test holds for
// (z > 0.9998, or z a NaN: every comparison false) takes isign*7, as the
// scalar (isign +1 for a NaN, as `x < 0` is false there too).
inline v8f platanh2_8(v8f tmn) {
    const v8f x = -tmn;
    const v8i neg = x < splat(0.0f);
    const v8f z = vsel(neg, -x, x);  // fabs(x) where x < 0
    const v8f isign = vsel(neg, splat(-1.0f), splat(1.0f));
    const v8i c1 = z <= splat(0.664f), c2 = z <= splat(0.9217f), c3 = z <= splat(0.9951f),
              c4 = z <= splat(0.9998f);
    v8f a = vsel(c3, splat(0.8378f), splat(0.9914f));
    a = vsel(c2, splat(0.4064f), a);
    v8f b = vsel(c3, splat(0.0524f), splat(0.0012f));
    b = vsel(c2, splat(0.322f), b);
    b = vsel(c1, splat(0.83f), b);
    const v8f num = vsel(c1, x, isign * (z - a));
    const v8f y = vsel(c4, num / b, isign * splat(7.0f));
    return splat(2.0f) * y;
}

// The tov of the eight checks in columns c..c+7: for each slot k, Tmn = 1.0
// times tanhtoc of the check's other slots in row order, then 2*platanh(-Tmn)
// (the Fortran's product(tanhtoc(1:nrw), mask=Nm.ne.j): its mask leaves out
// exactly slot k, BpTables128). nr rows: 11 in column group 0, where a ten-bit
// check's row 10 holds 1.0, its last factor (Tmn*1.0 is Tmn, every float and
// NaN); 10 elsewhere. A ten-bit check's slot 10 is computed and never read.
void check_tov8(const float* tanhtoc, int c, int nr, float* tov) {
    v8f t[kRows];
    for (int r = 0; r < nr; ++r) t[r] = vload(tanhtoc + r * kW + c);
    for (int k = 0; k < nr; ++k) {
        v8f tmn = splat(1.0f);
        for (int r = 0; r < nr; ++r)
            if (r != k) tmn = tmn * t[r];
        vstore(tov + k * kW + c, platanh2_8(tmn));
    }
}

}  // namespace

// platanh2_8 for the side-by-side test: out[l] = 2*platanh(-tmn[l]), l < 8.
void platanh2_n8(const float tmn[8], float out[8]) { vstore(out, platanh2_8(vload(tmn))); }

// bpdecode128_90 as bpdecode128_90_reference above, every float it computes
// the same, with its messages kept by edge (BpTables128):
// - zn(i): llr(i) + (tov(1,i) + tov(2,i) + tov(3,i)), the same sum in the same
//   order, each tov(kk,i) read from bit i's kk-th edge;
// - toc on an edge: zn of its bit less the edge's own tov (the one term the
//   Fortran's kk loop subtracts), then tanh(-toc/2) on 408 edge slots by
//   ft8::bp_tanhf, eight at a time (glibc 2.39's tanhf to the bit, as
//   tanhf_glibc), not on the eleven rows of all 38 checks: rows past a check's
//   nrw are never read;
// - each tov from the same product in the same order (check_tov8), eight
//   checks at a time;
// - the last iteration's messages (iter == maxiterations) are not computed:
//   the loop ends after them and nothing reads them (cw, the answer, comes
//   from that iteration's zn, before them).
// Early exits, iteration counts and outputs are the reference's. One
// exception, which no output shows: a NaN llr makes NaN messages, and their
// sign and payload may differ from the reference's (the compiler may compute
// -toc/2 as toc*-0.5, which keeps a NaN's sign). Nothing reads a NaN but
// `zn > 0` and platanh's comparisons, all false for any NaN, so every tov
// (always finite) and every output stay the same.
int bpdecode128_90(const float llr[128], int maxiterations, signed char message77[77], signed char cw_out[128],
                   int& nharderror, int& niterations) {
    const BpTables128& bt = bp_tables128();
    if (!bt.ok) return bpdecode128_90_reference(llr, maxiterations, message77, cw_out, nharderror, niterations);
    signed char cw[N];
    float zn[N];
    alignas(32) float tov[kRows * kW], toc[kRows * kW], arg[kRows * kW], tanhtoc[kRows * kW];
    std::memset(tov, 0, sizeof tov);
    // Row 10 of a ten-bit check in column group 0: 1.0 (set again after each
    // tanh); past group 0 that row is never read.
    for (int c = 0; c < kW; ++c) tanhtoc[10 * kW + c] = 1.0f;

    int ncnt = 0;
    int nclast = 0;
    int iter;
    for (iter = 0; iter <= maxiterations; ++iter) {
        // Update bit log likelihood ratios (tov=0 in iteration 0).
        for (int i = 0; i < N; ++i) {
            float s = 0.f;  // sum(tov(1:ncw,i)): gfortran starts from +0, so -0 terms give +0
            s += tov[bt.slot[i][0]];
            s += tov[bt.slot[i][1]];
            s += tov[bt.slot[i][2]];
            zn[i] = llr[i] + s;
        }

        if (bpdecode128_90_trace) {  // the piece check's hook, not in the Fortran: its 1-based arrays
            float zn1[N + 1] = {}, tov1[4][N + 1] = {};
            for (int i = 0; i < N; ++i) {
                zn1[i + 1] = zn[i];
                for (int kk = 0; kk < 3; ++kk) tov1[kk + 1][i + 1] = tov[bt.slot[i][kk]];
            }
            bpdecode128_90_trace(iter, zn1, tov1);
        }

        // Check to see if we have a codeword (check before we do any iteration).
        for (int i = 0; i < N; ++i) cw[i] = zn[i] > 0.f ? 1 : 0;
        int ncheck = 0;
        for (int i = 0; i < M; ++i) {
            int s = 0;
            for (int k = 0; k < kNrw[i]; ++k) s += cw[kNm[i][k] - 1];
            if (s % 2 != 0) ++ncheck;
        }
        if (ncheck == 0) {  // we have a codeword
            if (chkcrc13a(cw)) {  // decoded = cw(1:K)
                for (int i = 0; i < 77; ++i) message77[i] = cw[i];
                nharderror = 0;
                for (int i = 0; i < N; ++i)
                    if (static_cast<float>(2 * cw[i] - 1) * llr[i] < 0.0f) ++nharderror;
                std::memcpy(cw_out, cw, N);
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
                std::memcpy(cw_out, cw, N);
                niterations = iter;
                return nharderror;
            }
        }
        nclast = ncheck;
        if (iter == maxiterations) continue;  // its messages are never read

        // Send messages from bits to check nodes (the padding and unused slots
        // take finite rubbish, or fall back to the scalar tanh, and are never read)
        for (int e = 0; e < kNt; ++e) {
            toc[e] = zn[bt.ebit[e]] - tov[e];
            arg[e] = -toc[e] / 2;
        }
        // send messages from check nodes to variable nodes
        ft8::bp_tanhf(arg, tanhtoc, kNt);
        for (int c = bt.n11; c < 8; ++c) tanhtoc[10 * kW + c] = 1.0f;

        if (bpdecode128_90_msg_trace) {  // the piece check's hook: toc(1:nrw(j),j), tanhtoc(1:nrw(j),j)
            float toc1[12][M + 1] = {}, tanh1[12][M + 1] = {};
            for (int c = 0; c < M; ++c) {
                const int j = bt.col_check[c];
                for (int r = 0; r < kNrw[j]; ++r) {
                    toc1[r + 1][j + 1] = toc[r * kW + c];
                    tanh1[r + 1][j + 1] = tanhtoc[r * kW + c];
                }
            }
            bpdecode128_90_msg_trace(iter, toc1, tanh1);
        }

        for (int c = 0; c < kW; c += 8) check_tov8(tanhtoc, c, c < bt.n11 ? kRows : kRows - 1, tov);
    }
    nharderror = -1;
    std::memcpy(cw_out, cw, N);
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
