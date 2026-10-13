// FT8's (174,91) LDPC decoding (ft8_ldpc.h): WSJT-X's lib/ft8
// decode174_91.f90, osd174_91.f90 (with mrbencode91 and nextpat91),
// encode174_91.f90, encode174_91_nocrc.f90, get_crc14.f90, lib/crc14.cpp and
// lib/platanh.f90 at 567ad29ce, line for line. Arrays keep the Fortran's
// 1-based indices (element 0 unused) where that keeps the two readable side
// by side.
//
// Steven Franke, K9AN, and Joe Taylor, K1JT, WSJT-X (GNU GPL v3), and so is
// this.
#include "ft8_ldpc.h"

#include <cmath>
#include <cstring>
#include <stdexcept>

#if defined(__AVX2__) && defined(__x86_64__)
#include <immintrin.h>
#elif defined(__aarch64__) && defined(__ARM_NEON)
#include <arm_neon.h>
#endif

#include "wsprd_osd.h"  // indexx

namespace ft8 {

namespace {

#include "ft8_ldpc_tables.inc"

constexpr int N = 174, K = 91, M = N - K;

// encode174_91's (and _nocrc's) generator matrix, from the hex rows (filled
// once, as the Fortran's SAVEd gen).
struct Generator {
    int8_t gen[M][K];
    Generator() {
        std::memset(gen, 0, sizeof gen);
        for (int i = 0; i < M; ++i) {
            for (int j = 1; j <= 23; ++j) {
                const char c = kGen[i][j - 1];
                const int istr = (c >= '0' && c <= '9') ? c - '0' : (c >= 'a' && c <= 'f') ? c - 'a' + 10 : c - 'A' + 10;
                const int ibmax = j == 23 ? 3 : 4;
                for (int jj = 1; jj <= ibmax; ++jj) {
                    const int icol = (j - 1) * 4 + jj;
                    if ((istr >> (4 - jj)) & 1) gen[i][icol - 1] = 1;
                }
            }
        }
    }
};

const Generator& generator() {
    static const Generator g;
    return g;
}

// The parity checks of message (91 bits) by gen: pchecks(i) = mod(nsum, 2).
void parity(const int8_t* message, int8_t* codeword) {
    const auto& g = generator().gen;
    int8_t pchecks[M];
    for (int i = 0; i < M; ++i) {
        int nsum = 0;
        for (int j = 0; j < K; ++j) nsum += message[j] * g[i][j];
        pchecks[i] = static_cast<int8_t>(nsum % 2);
    }
    std::memcpy(codeword, message, K);
    std::memcpy(codeword + K, pchecks, M);
}

// A 174-bit word: bit j-1 of w is position j (w[3] is padding).
struct alignas(32) OsdWord {
    uint64_t w[4];
};

// osd174_91's SAVEd gen(k,N) for k = 91, row by row: each unit message
// encoded without its CRC (filled on the first call in the Fortran, with the
// first caller's k; every WSJT-X caller passes 91). For i <= 77 the CRC14 of
// the unit message is put in bits 78-91 and then message91(78:k) zeroed:
// with k = 91 all of them.
struct OsdGenerator {
    int8_t gen[K][N];
    OsdGenerator() {
        const int k = K;
        int8_t message91[91], m96[96], cw[N];
        for (int i = 1; i <= k; ++i) {
            std::memset(message91, 0, sizeof message91);
            message91[i - 1] = 1;
            if (i <= 77) {
                std::memset(m96, 0, sizeof m96);
                std::memcpy(m96, message91, 91);
                const int ncrc14 = get_crc14(m96, 96);
                for (int b = 0; b < 14; ++b) message91[77 + b] = static_cast<int8_t>((ncrc14 >> (13 - b)) & 1);
                for (int b = 78; b <= k; ++b) message91[b - 1] = 0;
            }
            encode174_91_nocrc(message91, cw);
            std::memcpy(gen[i - 1], cw, N);
        }
        // The same rows bit-packed, for osd174_91.
        std::memset(rows, 0, sizeof rows);
        for (int r = 0; r < K; ++r)
            for (int j = 0; j < N; ++j) rows[r].w[j >> 6] |= static_cast<uint64_t>(gen[r][j]) << (j & 63);
    }
    OsdWord rows[K];
};

const OsdGenerator& osd_generator() {
    static const OsdGenerator g;
    return g;
}

float bitsf(uint32_t u) {
    float f;
    std::memcpy(&f, &u, 4);
    return f;
}
uint32_t fbits(float f) {
    uint32_t u;
    std::memcpy(&u, &f, 4);
    return u;
}

// glibc's (2.39) expm1f, as tanhf calls it (sysdeps/ieee754/flt-32
// s_expm1f.c, fdlibm's), line for line; finite |x| < 88.72 only (tanhf's
// arguments: no overflow branch).
float expm1f_glibc(float x) {
    const float one = 1.0f, tiny = 1.0e-30f;
    const float ln2_hi = bitsf(0x3f317180), ln2_lo = bitsf(0x3717f7d1), invln2 = bitsf(0x3fb8aa3b);
    const float Q1 = bitsf(0xbd088889), Q2 = bitsf(0x3ad00d01), Q3 = bitsf(0xb8a670cd), Q4 = bitsf(0x36867e54),
                Q5 = bitsf(0xb457edbb);
    float hi, lo, c = 0.0f, t, e, hxs, hfx, r1, y;
    int32_t k;
    uint32_t hx = fbits(x);
    const uint32_t xsb = hx & 0x80000000u;  // sign bit of x
    hx &= 0x7fffffffu;                      // high word of |x|
    if (hx >= 0x4195b844u && xsb != 0) return tiny - one;  // x < -27*ln2, return -1.0
    // argument reduction
    if (hx > 0x3eb17218u) {      // if  |x| > 0.5 ln2
        if (hx < 0x3F851592u) {  // and |x| < 1.5 ln2
            if (xsb == 0) {
                hi = x - ln2_hi;
                lo = ln2_lo;
                k = 1;
            } else {
                hi = x + ln2_hi;
                lo = -ln2_lo;
                k = -1;
            }
        } else {
            k = static_cast<int32_t>(invln2 * x + ((xsb == 0) ? 0.5f : -0.5f));
            t = static_cast<float>(k);
            hi = x - t * ln2_hi;  // t*ln2_hi is exact here
            lo = t * ln2_lo;
        }
        x = hi - lo;
        c = (hi - x) - lo;
    } else if (hx < 0x33000000u) {  // when |x|<2**-25, return x
        return x;
    } else {
        k = 0;
    }
    // x is now in primary range
    hfx = 0.5f * x;
    hxs = x * hfx;
    r1 = one + hxs * (Q1 + hxs * (Q2 + hxs * (Q3 + hxs * (Q4 + hxs * Q5))));
    t = 3.0f - r1 * hfx;
    e = hxs * ((r1 - t) / (6.0f - x * t));
    if (k == 0) return x - (x * e - hxs);  // c is 0
    e = (x * (e - c) - c);
    e -= hxs;
    if (k == -1) return 0.5f * (x - e) - 0.5f;
    if (k == 1) {
        if (x < -0.25f) return -2.0f * (e - (x + 0.5f));
        return one + 2.0f * (x - e);
    }
    if (k <= -2 || k > 56) {  // suffice to return exp(x)-1
        y = one - (e - x);
        y = bitsf(fbits(y) + (static_cast<uint32_t>(k) << 23));  // add k to y's exponent
        return y - one;
    }
    if (k < 23) {
        t = bitsf(0x3f800000u - (0x1000000u >> k));  // t=1-2^-k
        y = t - (e - x);
    } else {
        t = bitsf(static_cast<uint32_t>(0x7f - k) << 23);  // 2^-k
        y = x - (e + t);
        y += one;
    }
    return bitsf(fbits(y) + (static_cast<uint32_t>(k) << 23));  // add k to y's exponent
}

}  // namespace

float tanhf_glibc(float x) {
    const float one = 1.0f, two = 2.0f, tiny = 1.0e-30f;
    float t, z;
    const uint32_t jx = fbits(x);
    const uint32_t ix = jx & 0x7fffffffu;
    if (ix >= 0x7f800000u) {  // x is INF or NaN
        if (!(jx >> 31)) return one / x + one;  // tanh(+-inf)=+-1
        return one / x - one;                   // tanh(NaN) = NaN
    }
    if (ix < 0x41b00000u) {               // |x|<22
        if (ix == 0) return x;            // x == +-0
        if (ix < 0x24000000u) return x * (one + x);  // |x|<2**-55: tanh(small) = small
        if (ix >= 0x3f800000u) {          // |x|>=1
            t = expm1f_glibc(two * std::fabs(x));
            z = one - two / (t + two);
        } else {
            t = expm1f_glibc(-two * std::fabs(x));
            z = -t / (t + two);
        }
    } else {  // |x| >= 22, return +-1
        z = one - tiny;
    }
    return !(jx >> 31) ? z : -z;
}

namespace {

// platanh(x, y): a piecewise-linear atanh.
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

// mrbencode91(me, codeword, g2, N, K): g2 is N x k, here column c of g2 is
// row c of rows (g2(:,c) = genmrb(c,:)). 1-based me and codeword.
void mrbencode91(const int8_t* me, int8_t* codeword, const int8_t* rows, int k) {
    std::memset(codeword + 1, 0, N);
    for (int i = 1; i <= k; ++i)
        if (me[i] == 1) {
            const int8_t* g = rows + static_cast<size_t>(i - 1) * N;
            for (int r = 1; r <= N; ++r) codeword[r] ^= g[r - 1];
        }
}

// nextpat91(mi, k, iorder, iflag): the next test pattern of weight iorder.
void nextpat91(int8_t* mi, int k, int iorder, int& iflag) {
    int ind = -1;
    for (int i = 1; i <= k - 1; ++i)
        if (mi[i] == 0 && mi[i + 1] == 1) ind = i;
    if (ind < 0) {  // no more patterns of this order
        iflag = ind;
        return;
    }
    int8_t ms[K + 1];
    std::memset(ms, 0, sizeof ms);
    for (int i = 1; i <= ind - 1; ++i) ms[i] = mi[i];
    ms[ind] = 1;
    ms[ind + 1] = 0;
    if (ind + 1 < k) {
        int s = 0;
        for (int i = 1; i <= k; ++i) s += ms[i];
        const int nz = iorder - s;
        for (int i = k - nz + 1; i <= k; ++i) ms[i] = 1;
    }
    std::memcpy(mi + 1, ms + 1, static_cast<size_t>(k));
    for (int i = 1; i <= k; ++i)  // iflag will point to the lowest-index 1 in mi
        if (mi[i] == 1) {
            iflag = i;
            break;
        }
}

// The CRC of a 91-bit word laid out as the decoders check it: m96(1:77) =
// cw(1:77), m96(78:82) = 0, m96(83:96) = cw(78:91).
int crc_of_m96(const int8_t* cw) {
    int8_t m96[96];
    std::memset(m96, 0, sizeof m96);
    std::memcpy(m96, cw, 77);
    std::memcpy(m96 + 82, cw + 77, 14);
    return get_crc14(m96, 96);
}

// decode174_91's message passing check by check, eight checks a vector: the
// checks in columns c (the seven-bit ones first, then the six-bit ones, each
// in Nm's order; padding columns past M), Nm's slots r = 0..6 in rows, so
// toc(r+1, order(c)+1) is at r*kW + c. Row 6 (a seventh bit) is used only
// in the first n7 columns, the rest of it holding tanh 1.0: a six-bit
// check's Tmn is then 1.0 times its five other terms times 1.0, the
// reference's product in its order (Tmn*1.0 is Tmn, NaN included). For each
// edge, its bit (ebit), and for each bit's Mn slot kk, the edge that check
// sends it (slot: tov(kk,j)). ok is false if the code were ever not this
// shape (then the reference runs).
struct BpTables {
    static constexpr int kW = (M + 7) / 8 * 8;  // columns, eight a vector
    static constexpr int kSize = 7 * kW;
    int n7 = 0;                // seven-bit checks: columns 0..n7-1
    int16_t ebit[kSize];       // the edge's bit, 0-based (0 for padding)
    int16_t slot[N][3];        // tov(kk+1, j+1)'s edge
    bool ok = true;
    BpTables() {
        int order[kW], c = 0;
        for (int w : {7, 6})
            for (int i = 0; i < M; ++i)
                if (kNrw[i] == w) order[c++] = i;
        if (c != M || kNcw != 3) ok = false;
        for (int i = 0; i < M; ++i)
            if (kNrw[i] == 7) ++n7;
        for (; c < kW; ++c) order[c] = -1;
        int nslot[N] = {};
        std::memset(ebit, 0, sizeof ebit);
        for (c = 0; c < kW && ok; ++c) {
            const int i = order[c];
            if (i < 0) continue;
            for (int r = 0; r < kNrw[i]; ++r) {
                const int b = kNm[i][r] - 1;
                ebit[r * kW + c] = static_cast<int16_t>(b);
                int nmatch = 0;
                for (int kk = 0; kk < kNcw; ++kk)
                    if (kMn[b][kk] == i + 1) {
                        slot[b][kk] = static_cast<int16_t>(r * kW + c);
                        ++nmatch;
                    }
                if (nmatch != 1) ok = false;
                ++nslot[b];
            }
        }
        for (int b = 0; b < N && ok; ++b)
            if (nslot[b] != 3) ok = false;
    }
};

const BpTables& bp_tables() {
    static const BpTables t;
    return t;
}

#if defined(__AVX2__) && defined(__x86_64__)
#define FT8_TANHF_AVX2 1
// tanhf_glibc and the expm1f it calls, eight lanes, each lane the scalar's
// operations in its order. Proven equal to tanhf_glibc, and that to glibc
// 2.39's tanhf on x86-64, for every float (tests/test_ft8_ldpc.cpp,
// [ft8-tanhf-all] and [ft8-tanhf-libm]).
inline __m256i ival(int v) { return _mm256_set1_epi32(v); }
inline __m256 fval(uint32_t bits) { return _mm256_castsi256_ps(_mm256_set1_epi32(static_cast<int>(bits))); }
inline __m256 sel(__m256 mask, __m256 a, __m256 b) { return _mm256_blendv_ps(b, a, mask); }  // mask ? a : b
inline __m256 msk(__m256i m) { return _mm256_castsi256_ps(m); }

// expm1f(x) for finite |x| < 88.72 (0x42b17218: no overflow branch).
__m256 expm1f8(__m256 x) {
    const __m256 one = _mm256_set1_ps(1.0f), half = _mm256_set1_ps(0.5f);
    const __m256 ln2_hi = fval(0x3f317180), ln2_lo = fval(0x3717f7d1), invln2 = fval(0x3fb8aa3b);
    const __m256i hx0 = _mm256_castps_si256(x);
    const __m256i hx = _mm256_and_si256(hx0, ival(0x7fffffff));
    const __m256 neg = msk(_mm256_srai_epi32(hx0, 31));  // xsb != 0
    // argument reduction: |x| in (0.5 ln2, 1.5 ln2): k = +-1
    const __m256 hi1 = sel(neg, _mm256_add_ps(x, ln2_hi), _mm256_sub_ps(x, ln2_hi));
    const __m256 lo1 = sel(neg, _mm256_xor_ps(ln2_lo, fval(0x80000000)), ln2_lo);
    const __m256i k1 = _mm256_castps_si256(sel(neg, msk(ival(-1)), msk(ival(1))));
    // further: k = invln2*x + (+-0.5), truncated
    const __m256i k2 = _mm256_cvttps_epi32(
        _mm256_add_ps(_mm256_mul_ps(invln2, x), sel(neg, _mm256_set1_ps(-0.5f), half)));
    const __m256 t2 = _mm256_cvtepi32_ps(k2);
    const __m256 hi2 = _mm256_sub_ps(x, _mm256_mul_ps(t2, ln2_hi));
    const __m256 lo2 = _mm256_mul_ps(t2, ln2_lo);
    const __m256 near = msk(_mm256_cmpgt_epi32(ival(0x3F851592), hx));
    const __m256 hi = sel(near, hi1, hi2), lo = sel(near, lo1, lo2);
    const __m256i kr = _mm256_castps_si256(sel(near, msk(k1), msk(k2)));
    const __m256 red = msk(_mm256_cmpgt_epi32(hx, ival(0x3eb17218)));
    const __m256 xr_red = _mm256_sub_ps(hi, lo);
    const __m256 c = _mm256_sub_ps(_mm256_sub_ps(hi, xr_red), lo);
    const __m256 xr = sel(red, xr_red, x);
    const __m256i k = _mm256_and_si256(kr, _mm256_castps_si256(red));  // 0 unless reduced
    // x is now in primary range
    const __m256 hfx = _mm256_mul_ps(half, xr);
    const __m256 hxs = _mm256_mul_ps(xr, hfx);
    __m256 p = _mm256_mul_ps(hxs, fval(0xb457edbb));             // hxs*Q5
    p = _mm256_mul_ps(hxs, _mm256_add_ps(fval(0x36867e54), p));  // hxs*(Q4+..)
    p = _mm256_mul_ps(hxs, _mm256_add_ps(fval(0xb8a670cd), p));  // hxs*(Q3+..)
    p = _mm256_mul_ps(hxs, _mm256_add_ps(fval(0x3ad00d01), p));  // hxs*(Q2+..)
    p = _mm256_mul_ps(hxs, _mm256_add_ps(fval(0xbd088889), p));  // hxs*(Q1+..)
    const __m256 r1 = _mm256_add_ps(one, p);
    const __m256 t = _mm256_sub_ps(_mm256_set1_ps(3.0f), _mm256_mul_ps(r1, hfx));
    __m256 e = _mm256_mul_ps(
        hxs, _mm256_div_ps(_mm256_sub_ps(r1, t), _mm256_sub_ps(_mm256_set1_ps(6.0f), _mm256_mul_ps(xr, t))));
    const __m256 res0 = _mm256_sub_ps(xr, _mm256_sub_ps(_mm256_mul_ps(xr, e), hxs));  // k == 0
    e = _mm256_sub_ps(_mm256_mul_ps(xr, _mm256_sub_ps(e, c)), c);
    e = _mm256_sub_ps(e, hxs);
    const __m256 xme = _mm256_sub_ps(xr, e), emx = _mm256_sub_ps(e, xr);
    const __m256 resm1 = _mm256_sub_ps(_mm256_mul_ps(half, xme), half);  // k == -1
    const __m256 resp1 = sel(_mm256_cmp_ps(xr, _mm256_set1_ps(-0.25f), _CMP_LT_OQ),  // k == 1
                             _mm256_mul_ps(_mm256_set1_ps(-2.0f), _mm256_sub_ps(e, _mm256_add_ps(xr, half))),
                             _mm256_add_ps(one, _mm256_mul_ps(_mm256_set1_ps(2.0f), xme)));
    const __m256i kexp = _mm256_slli_epi32(k, 23);
    // k <= -2 || k > 56
    const __m256 yb = _mm256_sub_ps(one, emx);
    const __m256 resbig = _mm256_sub_ps(_mm256_castsi256_ps(_mm256_add_epi32(_mm256_castps_si256(yb), kexp)), one);
    // 2 <= k < 23: t = 1-2^-k
    const __m256 tm = _mm256_castsi256_ps(_mm256_sub_epi32(ival(0x3f800000), _mm256_srlv_epi32(ival(0x1000000), k)));
    const __m256 ym = _mm256_sub_ps(tm, emx);
    const __m256 resmid = _mm256_castsi256_ps(_mm256_add_epi32(_mm256_castps_si256(ym), kexp));
    // 23 <= k <= 56: t = 2^-k
    const __m256 th = _mm256_castsi256_ps(_mm256_slli_epi32(_mm256_sub_epi32(ival(0x7f), k), 23));
    const __m256 yh = _mm256_add_ps(_mm256_sub_ps(xr, _mm256_add_ps(e, th)), one);
    const __m256 reshi = _mm256_castsi256_ps(_mm256_add_epi32(_mm256_castps_si256(yh), kexp));
    __m256 res = sel(msk(_mm256_cmpgt_epi32(k, ival(22))), reshi, resmid);
    res = sel(msk(_mm256_or_si256(_mm256_cmpgt_epi32(ival(-1), k), _mm256_cmpgt_epi32(k, ival(56)))), resbig, res);
    res = sel(msk(_mm256_cmpeq_epi32(k, ival(1))), resp1, res);
    res = sel(msk(_mm256_cmpeq_epi32(k, ival(-1))), resm1, res);
    res = sel(msk(_mm256_cmpeq_epi32(k, ival(0))), res0, res);
    // |x| < 2^-25: x; x < -27 ln2: tiny-one (-1)
    res = sel(msk(_mm256_cmpgt_epi32(ival(0x33000000), hx)), x, res);
    res = sel(_mm256_and_ps(neg, msk(_mm256_cmpgt_epi32(hx, ival(0x4195b843)))), _mm256_set1_ps(-1.0f), res);
    return res;
}

// tanhf(x) for eight finite x; false (nothing written) if any is not.
bool tanhf8(const float* in, float* out) {
    const __m256 x = _mm256_loadu_ps(in);
    const __m256i jx = _mm256_castps_si256(x);
    const __m256i ix = _mm256_and_si256(jx, ival(0x7fffffff));
    if (!_mm256_testz_si256(_mm256_cmpgt_epi32(ix, ival(0x7f7fffff)), _mm256_set1_epi32(-1))) return false;
    const __m256 one = _mm256_set1_ps(1.0f), two = _mm256_set1_ps(2.0f);
    const __m256 ax = _mm256_castsi256_ps(ix);
    const __m256 ge1 = msk(_mm256_cmpgt_epi32(ix, ival(0x3f7fffff)));  // |x| >= 1
    // expm1f(2|x|) or expm1f(-2|x|); z = 1 - 2/(t+2) or -t/(t+2)
    const __m256 tw = _mm256_add_ps(ax, ax);
    const __m256 t = expm1f8(sel(ge1, tw, _mm256_xor_ps(tw, fval(0x80000000))));
    const __m256 tp2 = _mm256_add_ps(t, two);
    __m256 z = sel(ge1, _mm256_sub_ps(one, _mm256_div_ps(two, tp2)),
                   _mm256_div_ps(_mm256_xor_ps(t, fval(0x80000000)), tp2));
    z = sel(msk(_mm256_cmpgt_epi32(ix, ival(0x41afffff))), _mm256_sub_ps(one, fval(0x0da24260)), z);  // |x| >= 22: 1-tiny
    z = _mm256_xor_ps(z, _mm256_castsi256_ps(_mm256_and_si256(jx, ival(static_cast<int>(0x80000000u)))));  // jx < 0: -z
    // |x| < 2^-55: x*(1+x); x = +-0: x
    z = sel(msk(_mm256_cmpgt_epi32(ival(0x24000000), ix)), _mm256_mul_ps(x, _mm256_add_ps(one, x)), z);
    z = sel(msk(_mm256_cmpeq_epi32(ix, ival(0))), x, z);
    _mm256_storeu_ps(out, z);
    return true;
}

// 2*platanh(-Tmn) for eight Tmn, as the scalar: platanh's branches as
// selects, in their order (an unordered z takes none: 7*isign).
inline __m256 platanh2_8(__m256 tmn) {
    const __m256 sign = fval(0x80000000);
    const __m256 x = _mm256_xor_ps(tmn, sign);  // -Tmn
    const __m256 neg = _mm256_cmp_ps(x, _mm256_setzero_ps(), _CMP_LT_OQ);
    const __m256 z = sel(neg, _mm256_andnot_ps(sign, x), x);
    const __m256 isign = sel(neg, _mm256_set1_ps(-1.0f), _mm256_set1_ps(1.0f));
    auto piece = [&](float a, float b) { return _mm256_div_ps(_mm256_mul_ps(isign, _mm256_sub_ps(z, _mm256_set1_ps(a))), _mm256_set1_ps(b)); };
    __m256 y = _mm256_mul_ps(isign, _mm256_set1_ps(7.0f));
    y = sel(_mm256_cmp_ps(z, _mm256_set1_ps(0.9998f), _CMP_LE_OQ), piece(0.9914f, 0.0012f), y);
    y = sel(_mm256_cmp_ps(z, _mm256_set1_ps(0.9951f), _CMP_LE_OQ), piece(0.8378f, 0.0524f), y);
    y = sel(_mm256_cmp_ps(z, _mm256_set1_ps(0.9217f), _CMP_LE_OQ), piece(0.4064f, 0.322f), y);
    y = sel(_mm256_cmp_ps(z, _mm256_set1_ps(0.664f), _CMP_LE_OQ), _mm256_div_ps(x, _mm256_set1_ps(0.83f)), y);
    return _mm256_mul_ps(_mm256_set1_ps(2.0f), y);
}

// tov for the eight checks in columns c..c+7: for each slot k, Tmn = 1.0
// times tanhtoc of the check's other slots in r order (row 6's 1.0 last
// for a six-bit check), then 2*platanh(-Tmn). Slot 6 only if seven.
void check_tov8(const float* tanhtoc, int c, bool seven, float* tov) {
    constexpr int kW = BpTables::kW;
    __m256 t[7];
    for (int r = 0; r < 7; ++r) t[r] = _mm256_load_ps(tanhtoc + r * kW + c);
    const int nk = seven ? 7 : 6;
    for (int k = 0; k < nk; ++k) {
        __m256 tmn = _mm256_set1_ps(1.0f);
        for (int r = 0; r < 7; ++r)
            if (r != k) tmn = _mm256_mul_ps(tmn, t[r]);
        _mm256_store_ps(tov + k * kW + c, platanh2_8(tmn));
    }
}
#elif defined(__aarch64__) && defined(__ARM_NEON)
#define FT8_TANHF_NEON 1
// The AVX2 code above on NEON, four lanes a register (eight as two): each
// line the same operations on the same values in the same order, multiply
// then add (never vfmaq/vmlaq: -ffp-contract=off keeps vmulq and vaddq
// apart), selects as bit selects of full masks.
inline int32x4_t ival(int v) { return vdupq_n_s32(v); }
inline float32x4_t fval(uint32_t bits) { return vreinterpretq_f32_u32(vdupq_n_u32(bits)); }
inline float32x4_t sel(uint32x4_t mask, float32x4_t a, float32x4_t b) { return vbslq_f32(mask, a, b); }  // mask ? a : b
inline uint32x4_t gt(int32x4_t a, int32x4_t b) { return vcgtq_s32(a, b); }
inline uint32x4_t eq(int32x4_t a, int32x4_t b) { return vceqq_s32(a, b); }
inline float32x4_t xor_bits(float32x4_t a, float32x4_t b) {
    return vreinterpretq_f32_u32(veorq_u32(vreinterpretq_u32_f32(a), vreinterpretq_u32_f32(b)));
}
inline int32x4_t bits_of(float32x4_t x) { return vreinterpretq_s32_f32(x); }
inline float32x4_t float_of(int32x4_t x) { return vreinterpretq_f32_s32(x); }

// expm1f(x) for finite |x| < 88.72 (0x42b17218: no overflow branch).
float32x4_t expm1f4(float32x4_t x) {
    const float32x4_t one = vdupq_n_f32(1.0f), half = vdupq_n_f32(0.5f);
    const float32x4_t ln2_hi = fval(0x3f317180), ln2_lo = fval(0x3717f7d1), invln2 = fval(0x3fb8aa3b);
    const int32x4_t hx0 = bits_of(x);
    const int32x4_t hx = vandq_s32(hx0, ival(0x7fffffff));
    const uint32x4_t neg = vreinterpretq_u32_s32(vshrq_n_s32(hx0, 31));  // xsb != 0
    // argument reduction: |x| in (0.5 ln2, 1.5 ln2): k = +-1
    const float32x4_t hi1 = sel(neg, vaddq_f32(x, ln2_hi), vsubq_f32(x, ln2_hi));
    const float32x4_t lo1 = sel(neg, xor_bits(ln2_lo, fval(0x80000000)), ln2_lo);
    const int32x4_t k1 = vbslq_s32(neg, ival(-1), ival(1));
    // further: k = invln2*x + (+-0.5), truncated
    const int32x4_t k2 = vcvtq_s32_f32(vaddq_f32(vmulq_f32(invln2, x), sel(neg, vdupq_n_f32(-0.5f), half)));
    const float32x4_t t2 = vcvtq_f32_s32(k2);
    const float32x4_t hi2 = vsubq_f32(x, vmulq_f32(t2, ln2_hi));
    const float32x4_t lo2 = vmulq_f32(t2, ln2_lo);
    const uint32x4_t near = gt(ival(0x3F851592), hx);
    const float32x4_t hi = sel(near, hi1, hi2), lo = sel(near, lo1, lo2);
    const int32x4_t kr = vbslq_s32(near, k1, k2);
    const uint32x4_t red = gt(hx, ival(0x3eb17218));
    const float32x4_t xr_red = vsubq_f32(hi, lo);
    const float32x4_t c = vsubq_f32(vsubq_f32(hi, xr_red), lo);
    const float32x4_t xr = sel(red, xr_red, x);
    const int32x4_t k = vandq_s32(kr, vreinterpretq_s32_u32(red));  // 0 unless reduced
    // x is now in primary range
    const float32x4_t hfx = vmulq_f32(half, xr);
    const float32x4_t hxs = vmulq_f32(xr, hfx);
    float32x4_t p = vmulq_f32(hxs, fval(0xb457edbb));         // hxs*Q5
    p = vmulq_f32(hxs, vaddq_f32(fval(0x36867e54), p));      // hxs*(Q4+..)
    p = vmulq_f32(hxs, vaddq_f32(fval(0xb8a670cd), p));      // hxs*(Q3+..)
    p = vmulq_f32(hxs, vaddq_f32(fval(0x3ad00d01), p));      // hxs*(Q2+..)
    p = vmulq_f32(hxs, vaddq_f32(fval(0xbd088889), p));      // hxs*(Q1+..)
    const float32x4_t r1 = vaddq_f32(one, p);
    const float32x4_t t = vsubq_f32(vdupq_n_f32(3.0f), vmulq_f32(r1, hfx));
    float32x4_t e = vmulq_f32(hxs, vdivq_f32(vsubq_f32(r1, t), vsubq_f32(vdupq_n_f32(6.0f), vmulq_f32(xr, t))));
    const float32x4_t res0 = vsubq_f32(xr, vsubq_f32(vmulq_f32(xr, e), hxs));  // k == 0
    e = vsubq_f32(vmulq_f32(xr, vsubq_f32(e, c)), c);
    e = vsubq_f32(e, hxs);
    const float32x4_t xme = vsubq_f32(xr, e), emx = vsubq_f32(e, xr);
    const float32x4_t resm1 = vsubq_f32(vmulq_f32(half, xme), half);  // k == -1
    const float32x4_t resp1 = sel(vcltq_f32(xr, vdupq_n_f32(-0.25f)),  // k == 1
                                  vmulq_f32(vdupq_n_f32(-2.0f), vsubq_f32(e, vaddq_f32(xr, half))),
                                  vaddq_f32(one, vmulq_f32(vdupq_n_f32(2.0f), xme)));
    const int32x4_t kexp = vshlq_n_s32(k, 23);
    // k <= -2 || k > 56
    const float32x4_t yb = vsubq_f32(one, emx);
    const float32x4_t resbig = vsubq_f32(float_of(vaddq_s32(bits_of(yb), kexp)), one);
    // 2 <= k < 23: t = 1-2^-k (0x1000000 >> k; a lane outside 0..31 is never taken)
    const int32x4_t sh = vreinterpretq_s32_u32(vshlq_u32(vdupq_n_u32(0x1000000), vnegq_s32(k)));
    const float32x4_t tm = float_of(vsubq_s32(ival(0x3f800000), sh));
    const float32x4_t ym = vsubq_f32(tm, emx);
    const float32x4_t resmid = float_of(vaddq_s32(bits_of(ym), kexp));
    // 23 <= k <= 56: t = 2^-k
    const float32x4_t th = float_of(vshlq_n_s32(vsubq_s32(ival(0x7f), k), 23));
    const float32x4_t yh = vaddq_f32(vsubq_f32(xr, vaddq_f32(e, th)), one);
    const float32x4_t reshi = float_of(vaddq_s32(bits_of(yh), kexp));
    float32x4_t res = sel(gt(k, ival(22)), reshi, resmid);
    res = sel(vorrq_u32(gt(ival(-1), k), gt(k, ival(56))), resbig, res);
    res = sel(eq(k, ival(1)), resp1, res);
    res = sel(eq(k, ival(-1)), resm1, res);
    res = sel(eq(k, ival(0)), res0, res);
    // |x| < 2^-25: x; x < -27 ln2: tiny-one (-1)
    res = sel(gt(ival(0x33000000), hx), x, res);
    res = sel(vandq_u32(neg, gt(hx, ival(0x4195b843))), vdupq_n_f32(-1.0f), res);
    return res;
}

// tanhf(x) for four finite x.
float32x4_t tanhf4(float32x4_t x) {
    const int32x4_t jx = bits_of(x);
    const int32x4_t ix = vandq_s32(jx, ival(0x7fffffff));
    const float32x4_t one = vdupq_n_f32(1.0f), two = vdupq_n_f32(2.0f);
    const float32x4_t ax = float_of(ix);
    const uint32x4_t ge1 = gt(ix, ival(0x3f7fffff));  // |x| >= 1
    // expm1f(2|x|) or expm1f(-2|x|); z = 1 - 2/(t+2) or -t/(t+2)
    const float32x4_t tw = vaddq_f32(ax, ax);
    const float32x4_t t = expm1f4(sel(ge1, tw, xor_bits(tw, fval(0x80000000))));
    const float32x4_t tp2 = vaddq_f32(t, two);
    float32x4_t z = sel(ge1, vsubq_f32(one, vdivq_f32(two, tp2)), vdivq_f32(xor_bits(t, fval(0x80000000)), tp2));
    z = sel(gt(ix, ival(0x41afffff)), vsubq_f32(one, fval(0x0da24260)), z);  // |x| >= 22: 1-tiny
    z = xor_bits(z, float_of(vandq_s32(jx, ival(static_cast<int>(0x80000000u)))));  // jx < 0: -z
    // |x| < 2^-55: x*(1+x); x = +-0: x
    z = sel(gt(ival(0x24000000), ix), vmulq_f32(x, vaddq_f32(one, x)), z);
    z = sel(eq(ix, ival(0)), x, z);
    return z;
}

// tanhf(x) for eight finite x; false (nothing written) if any is not.
bool tanhf8(const float* in, float* out) {
    const float32x4_t lo = vld1q_f32(in), hi = vld1q_f32(in + 4);
    const uint32x4_t bad = vorrq_u32(gt(vandq_s32(bits_of(lo), ival(0x7fffffff)), ival(0x7f7fffff)),
                                     gt(vandq_s32(bits_of(hi), ival(0x7fffffff)), ival(0x7f7fffff)));
    if (vmaxvq_u32(bad) != 0) return false;
    vst1q_f32(out, tanhf4(lo));
    vst1q_f32(out + 4, tanhf4(hi));
    return true;
}

// 2*platanh(-Tmn) for four Tmn, as the scalar: platanh's branches as
// selects, in their order (an unordered z takes none: 7*isign).
inline float32x4_t platanh2_4(float32x4_t tmn) {
    const float32x4_t sign = fval(0x80000000);
    const float32x4_t x = xor_bits(tmn, sign);  // -Tmn
    const uint32x4_t neg = vcltq_f32(x, vdupq_n_f32(0.0f));
    const float32x4_t z =
        sel(neg, vreinterpretq_f32_u32(vbicq_u32(vreinterpretq_u32_f32(x), vreinterpretq_u32_f32(sign))), x);
    const float32x4_t isign = sel(neg, vdupq_n_f32(-1.0f), vdupq_n_f32(1.0f));
    auto piece = [&](float a, float b) { return vdivq_f32(vmulq_f32(isign, vsubq_f32(z, vdupq_n_f32(a))), vdupq_n_f32(b)); };
    float32x4_t y = vmulq_f32(isign, vdupq_n_f32(7.0f));
    y = sel(vcleq_f32(z, vdupq_n_f32(0.9998f)), piece(0.9914f, 0.0012f), y);
    y = sel(vcleq_f32(z, vdupq_n_f32(0.9951f)), piece(0.8378f, 0.0524f), y);
    y = sel(vcleq_f32(z, vdupq_n_f32(0.9217f)), piece(0.4064f, 0.322f), y);
    y = sel(vcleq_f32(z, vdupq_n_f32(0.664f)), vdivq_f32(x, vdupq_n_f32(0.83f)), y);
    return vmulq_f32(vdupq_n_f32(2.0f), y);
}

// check_tov8 as the AVX2 one, its eight checks as two halves of four.
void check_tov8(const float* tanhtoc, int c, bool seven, float* tov) {
    constexpr int kW = BpTables::kW;
    const int nk = seven ? 7 : 6;
    for (int h = 0; h < 8; h += 4) {
        float32x4_t t[7];
        for (int r = 0; r < 7; ++r) t[r] = vld1q_f32(tanhtoc + r * kW + c + h);
        for (int k = 0; k < nk; ++k) {
            float32x4_t tmn = vdupq_n_f32(1.0f);
            for (int r = 0; r < 7; ++r)
                if (r != k) tmn = vmulq_f32(tmn, t[r]);
            vst1q_f32(tov + k * kW + c + h, platanh2_4(tmn));
        }
    }
}
#endif
#if defined(FT8_TANHF_AVX2) || defined(FT8_TANHF_NEON)
#define FT8_TANHF_VEC 1
#endif

}  // namespace

bool bp_tanhf_vector() {
#ifdef FT8_TANHF_VEC
    return true;
#else
    return false;
#endif
}

void bp_tanhf(const float* x, float* y, int n) {
    int i = 0;
#ifdef FT8_TANHF_VEC
    for (; i + 8 <= n; i += 8)
        if (!tanhf8(x + i, y + i))
            for (int l = i; l < i + 8; ++l) y[l] = tanhf_glibc(x[l]);
#endif
    for (; i < n; ++i) y[i] = tanhf_glibc(x[i]);
}

int get_crc14(const int8_t* mc, int len) {
    // polynomial for 14-bit CRC 0x6757
    static const int8_t p[15] = {1, 1, 0, 0, 1, 1, 1, 0, 1, 0, 1, 0, 1, 1, 1};
    int8_t r[15];
    std::memcpy(r, mc, 15);
    for (int i = 0; i <= len - 15; ++i) {
        r[14] = mc[i + 14];
        const int8_t r1 = r[0];
        for (int j = 0; j < 15; ++j) r[j] = static_cast<int8_t>((r[j] + r1 * p[j]) % 2);
        const int8_t first = r[0];  // cshift(r, 1)
        for (int j = 0; j < 14; ++j) r[j] = r[j + 1];
        r[14] = first;
    }
    int ncrc = 0;
    for (int j = 0; j < 14; ++j) ncrc = 2 * ncrc + r[j];
    return ncrc;
}

int crc14(const unsigned char* data, int length) {
    // boost::augmented_crc<14, 0x2757>: the bits shifted through the
    // register most significant first, the polynomial added each time a 1
    // leaves it; no initial or final XOR, no reflection.
    unsigned rem = 0;
    for (int i = 0; i < length; ++i)
        for (int b = 7; b >= 0; --b) {
            const unsigned top = (rem >> 13) & 1u;
            rem = ((rem << 1) | ((data[i] >> b) & 1u)) & 0x3fffu;
            if (top) rem ^= 0x2757u;
        }
    return static_cast<int>(rem);
}

void encode174_91(const int8_t message77[77], int8_t codeword[174]) {
    // Add 14-bit CRC to form 91-bit message+CRC14: the 77 bits and three
    // zeros as ten bytes, two zero bytes after them for the CRC.
    unsigned char i1MsgBytes[12] = {};
    for (int i = 0; i < 80; ++i)
        if (i < 77 && message77[i] == 1) i1MsgBytes[i / 8] = static_cast<unsigned char>(i1MsgBytes[i / 8] | (0x80u >> (i % 8)));
    const int ncrc14 = crc14(i1MsgBytes, 12);
    int8_t message[K];
    std::memcpy(message, message77, 77);
    for (int b = 0; b < 14; ++b) message[77 + b] = static_cast<int8_t>((ncrc14 >> (13 - b)) & 1);
    parity(message, codeword);
}

void encode174_91_nocrc(const int8_t message[91], int8_t codeword[174]) { parity(message, codeword); }

// decode174_91 as decode174_91_reference below, its message passing check
// by check (BpTables): tov kept by edge, toc's subtraction by the edge's own
// tov, tanh only on the 522 edges (a six-bit check's seventh slot, tanh(-0),
// is never read) and eight at a time (bp_tanhf: glibc 2.39's tanhf to the bit),
// each Tmn the product over the same terms in the same order and platanh,
// eight checks at a time (check_tov8). Every float is what the reference
// computes.
void decode174_91(const float llr[174], int Keff, int maxosd, int norder, const int8_t apmask[174],
                  int8_t message91[91], int8_t cw[174], int& ntype, int& nharderror, float& dmin) {
    const BpTables& bt = bp_tables();
    if (!bt.ok) {
        decode174_91_reference(llr, Keff, maxosd, norder, apmask, message91, cw, ntype, nharderror, dmin);
        return;
    }
    constexpr int kW = BpTables::kW;
    const int maxiterations = 30;
    float zsave[3][N];  // the Fortran's local zsave(N,3)
    int nosd = 0;
    if (maxosd > 3) maxosd = 3;
    if (maxosd == 0) {  // osd with channel llrs
        nosd = 1;
        std::memcpy(zsave[0], llr, sizeof(float) * N);
    } else if (maxosd > 0) {
        nosd = maxosd;
    } else {  // just bp
        nosd = 0;
    }

    float zn[N], zsum[N];
    alignas(32) float tov[BpTables::kSize], arg[BpTables::kSize], tanhtoc[BpTables::kSize];
    std::memset(tov, 0, sizeof tov);
    for (int c = 0; c < kW; ++c) tanhtoc[6 * kW + c] = 1.0f;  // row 6 past n7: 1.0
    const int nt = 6 * kW + bt.n7;  // edges (and padding) to take tanh of

    int ncnt = 0, nclast = 0;
    for (int i = 0; i < N; ++i) zsum[i] = 0.0f;
    for (int iter = 0; iter <= maxiterations; ++iter) {
        // Update bit log likelihood ratios (tov=0 in iteration 0).
        for (int i = 0; i < N; ++i) {
            if (apmask[i] != 1) {
                float s = 0.0f;
                for (int kk = 0; kk < kNcw; ++kk) s += tov[bt.slot[i][kk]];
                zn[i] = llr[i] + s;
            } else {
                zn[i] = llr[i];
            }
        }
        for (int i = 0; i < N; ++i) zsum[i] = zsum[i] + zn[i];
        if (iter > 0 && iter <= maxosd) std::memcpy(zsave[iter - 1], zsum, sizeof zsum);

        // Check to see if we have a codeword (check before we do any iteration).
        for (int i = 0; i < N; ++i) cw[i] = zn[i] > 0.0f ? 1 : 0;
        int ncheck = 0;
        for (int i = 0; i < M; ++i) {
            int synd = 0;
            for (int r = 0; r < kNrw[i]; ++r) synd += cw[kNm[i][r] - 1];
            if (synd % 2 != 0) ++ncheck;
        }
        if (ncheck == 0) {  // we have a codeword - if crc is good, return it
            const int nbadcrc = crc_of_m96(cw);
            int nh = 0;  // count((2*cw-1)*llr .lt. 0.0)
            for (int i = 0; i < N; ++i)
                if (static_cast<float>(2 * cw[i] - 1) * llr[i] < 0.0f) ++nh;
            nharderror = nh;
            if (nbadcrc == 0) {
                std::memcpy(message91, cw, 91);
                float d = 0.0f;
                for (int i = 0; i < N; ++i) {
                    const int8_t hdec = llr[i] >= 0 ? 1 : 0;
                    const int8_t nxor = static_cast<int8_t>(hdec ^ cw[i]);
                    d += static_cast<float>(nxor) * std::fabs(llr[i]);
                }
                dmin = d;
                ntype = 1;
                return;
            }
        }

        if (iter > 0) {  // an early stopping criterion
            const int nd = ncheck - nclast;
            if (nd < 0) ncnt = 0;  // # of unsatisfied parity checks decreased
            else ++ncnt;
            if (ncnt >= 5 && iter >= 10 && ncheck > 15) {
                nharderror = -1;
                break;
            }
        }
        nclast = ncheck;

        // Send messages from bits to check nodes (toc = zn less what the
        // bit had received from the check), and from check nodes to
        // variable nodes: tanh(-toc/2) on each edge (the padding columns'
        // finite rubbish is never read).
        for (int e = 0; e < nt; ++e) arg[e] = -(zn[bt.ebit[e]] - tov[e]) / 2;
        bp_tanhf(arg, tanhtoc, nt);

#ifdef FT8_TANHF_VEC
        for (int c = 0; c < kW; c += 8) check_tov8(tanhtoc, c, c < bt.n7, tov);
#else
        for (int c = 0; c < kW; ++c)
            for (int k = 0; k < (c < bt.n7 ? 7 : 6); ++k) {
                float Tmn = 1.0f;
                for (int r = 0; r < 7; ++r)
                    if (r != k) Tmn = Tmn * tanhtoc[r * kW + c];
                tov[k * kW + c] = 2 * platanh(-Tmn);
            }
#endif
    }  // bp iterations

    for (int i = 1; i <= nosd; ++i) {
        std::memcpy(zn, zsave[i - 1], sizeof zn);
        float dminosd = 0.0f;
        osd174_91(zn, Keff, apmask, norder, message91, cw, nharderror, dminosd);
        if (nharderror > 0) {
            float d = 0.0f;
            for (int r = 0; r < N; ++r) {
                const int8_t hdec = llr[r] >= 0 ? 1 : 0;
                const int8_t nxor = static_cast<int8_t>(hdec ^ cw[r]);
                d += static_cast<float>(nxor) * std::fabs(llr[r]);
            }
            dmin = d;  // against the channel llrs; nharderror stays the OSD's
            ntype = 2;
            return;
        }
    }

    ntype = 0;
    nharderror = -1;
}

// osd174_91 on bit-packed words, giving exactly what the line-for-line
// osd174_91_reference below gives:
// - genmrb's rows stay in the generator's own column order (column j at bit
//   j-1 of a row's words) and its column reorders are kept in indices only:
//   the reference's G(r,i) is bit indices(i)-1 of row r. The elimination
//   picks the same pivots and XORs the same rows (XOR is exact in any
//   order), and the codewords are kept in that order too, so the final
//   re-order into [message bits][parity bits] is a plain copy;
// - only order 1 is reached (ndeep 1 and 2), where misub holds a single 1
//   at iflag and nextpat91 moves it down by one: a test pattern's message is
//   m0 ^ e(iflag) [^ e(n1)] and mrbencode91's codeword c0 XOR those rows;
//   the pattern's a-priori test is apmaskr(iflag) (masked: the whole group
//   is skipped) and apmaskr(n1);
// - every float sum is the reference's sum of nxor(i)*absrx(i) over the
//   reordered i in index order, with its zero terms left out: absrx >= +0
//   and finite, so each term left out is +0 and s + (+0) == s for the
//   s >= +0 seen here (d1 is the one term absrx(iflag)). Non-finite llrs go
//   to the reference.
void osd174_91(const float llr[174], int k, const int8_t apmask[174], int ndeep, int8_t message91[91],
               int8_t cw_out[174], int& nhardmin, float& dmin) {
    if (k != K) throw std::invalid_argument("osd174_91: k must be 91");
    // ndeep 3 and more (and below 0, the Fortran's else branch) reach
    // boxit91/fetchit91 (npre2 = 1): never from WSJT-X's (174,91) callers.
    if (ndeep < 0 || ndeep >= 3) throw std::invalid_argument("osd174_91: only ndeep 0-2 carried");
    for (int i = 0; i < N; ++i)
        if (!std::isfinite(llr[i])) {
            osd174_91_reference(llr, k, apmask, ndeep, message91, cw_out, nhardmin, dmin);
            return;
        }

    using Word = OsdWord;
    auto bit = [](const Word& w, int j) -> int { return static_cast<int>((w.w[j >> 6] >> (j & 63)) & 1); };
    auto popcount = [](const Word& w) {
        return __builtin_popcountll(w.w[0]) + __builtin_popcountll(w.w[1]) + __builtin_popcountll(w.w[2]);
    };

    // 1-based working arrays, as the Fortran's; absrx, apmaskr and hdec are
    // the reordered ones (rx is never read after its reorder: dropped).
    int8_t hdec[N + 1], apmaskr[N + 1];
    int indices[N + 1], indx[N + 1];
    float absrx0[N + 1], absrx[N + 1];
    Word rows[K + 1];  // rows[r] is genmrb(r,:), g2(:,r), in generator column order
    std::memcpy(rows + 1, osd_generator().rows, sizeof(Word) * K);

    // Use magnitude of received symbols as a measure of reliability.
    for (int i = 1; i <= N; ++i) absrx0[i] = std::fabs(llr[i - 1]);
    wsprd::indexx(absrx0 + 1, N, indx + 1);

    // Re-order the columns of the generator matrix in order of decreasing reliability.
    for (int i = 1; i <= N; ++i) indices[i] = indx[N + 1 - i];

    // Do gaussian elimination to create a generator matrix with the most reliable
    // received bits in positions 1:k in order of decreasing reliability (more or less).
    for (int id = 1; id <= k; ++id) {  // diagonal element indices
        for (int icol = id; icol <= k + 20; ++icol) {  // The 20 is ad hoc - beware
            if (bit(rows[id], indices[icol] - 1)) {
                if (icol != id) {  // reorder column
                    const int itmp = indices[id];
                    indices[id] = indices[icol];
                    indices[icol] = itmp;
                }
                // Every row ii /= id with G(ii,id) = 1 gets row id XORed in;
                // row id itself is put back after.
                const Word piv = rows[id];
                const int j = indices[id] - 1;
                const int jw = j >> 6, jb = j & 63;
                for (int ii = 1; ii <= k; ++ii) {
                    const uint64_t mask = uint64_t{0} - ((rows[ii].w[jw] >> jb) & 1);
                    for (int w = 0; w < 4; ++w) rows[ii].w[w] ^= piv.w[w] & mask;
                }
                rows[id] = piv;
                break;
            }
        }
    }

    // The hard decisions for the k MRB bits define the order 0 message, m0.
    Word hd{{0, 0, 0, 0}};  // the hard decisions, generator column order
    for (int j = 0; j < N; ++j) hd.w[j >> 6] |= static_cast<uint64_t>(llr[j] >= 0 ? 1 : 0) << (j & 63);
    for (int i = 1; i <= N; ++i) {
        hdec[i] = static_cast<int8_t>(llr[indices[i] - 1] >= 0 ? 1 : 0);
        absrx[i] = absrx0[indices[i]];
        apmaskr[i] = apmask[indices[i] - 1];
    }
    // The parity positions k+1:N and their first nt = 40, as masks.
    Word par{{0, 0, 0, 0}}, par40{{0, 0, 0, 0}};
    for (int i = k + 1; i <= N; ++i) {
        const int j = indices[i] - 1;
        par.w[j >> 6] |= uint64_t{1} << (j & 63);
        if (i <= k + 40) par40.w[j >> 6] |= uint64_t{1} << (j & 63);
    }
    // s plus absrx(i) for each reordered i in lo:N where e's bit is set, in
    // order.
    auto sum_set = [&](float s, const Word& e, int lo) {
        for (int i = lo; i <= N; ++i)
            if (bit(e, indices[i] - 1)) s += absrx[i];
        return s;
    };

    Word c0{{0, 0, 0, 0}};  // mrbencode91(m0): m0 = hdec(1:k)
    for (int i = 1; i <= k; ++i)
        if (hdec[i] == 1)
            for (int w = 0; w < 3; ++w) c0.w[w] ^= rows[i].w[w];
    Word nxor;
    for (int w = 0; w < 4; ++w) nxor.w[w] = c0.w[w] ^ hd.w[w];
    nhardmin = popcount(nxor);
    dmin = sum_set(0.0f, nxor, 1);
    Word cw = c0;

    if (ndeep != 0) {  // else goto 998 (order 0)
        int npre1 = 0, ntheta = 12;  // nord 1 and nt 40 for both
        if (ndeep == 1) {
            npre1 = 0;
            ntheta = 12;
        } else if (ndeep == 2) {
            npre1 = 1;
            ntheta = 10;
        }
        Word r40[K + 1];  // each row's bits at the first nt parity positions
        for (int r = 1; r <= k; ++r)
            for (int w = 0; w < 4; ++w) r40[r].w[w] = rows[r].w[w] & par40.w[w];

        for (int iflag = k; iflag >= 1; --iflag) {  // misub = e(iflag), iorder 1
            if ((apmaskr[iflag] & 1) == 1) continue;  // every n1 of the group masked
            const int iend = npre1 == 0 ? iflag : 1;
            // n1 == iflag: me = m0 ^ e(iflag), e2sub its parity errors.
            Word csub, e2sub;
            for (int w = 0; w < 4; ++w) {
                csub.w[w] = c0.w[w] ^ rows[iflag].w[w];
                e2sub.w[w] = (csub.w[w] ^ hd.w[w]) & par.w[w];
            }
            const float d1 = absrx[iflag];
            const uint64_t s0 = e2sub.w[0] & par40.w[0], s1 = e2sub.w[1] & par40.w[1], s2 = e2sub.w[2] & par40.w[2];
            const int nd1kpt = 1 + __builtin_popcountll(s0) + __builtin_popcountll(s1) + __builtin_popcountll(s2);
            if (nd1kpt <= ntheta) {
                const float dd = d1 + sum_set(0.0f, e2sub, k + 1);
                if (dd < dmin) {
                    dmin = dd;
                    cw = csub;
                    Word nx;
                    for (int w = 0; w < 4; ++w) nx.w[w] = csub.w[w] ^ hd.w[w];
                    nhardmin = popcount(nx);
                }
            }
            for (int n1 = iflag - 1; n1 >= iend; --n1) {
                if ((apmaskr[n1] & 1) == 1) continue;
                const Word& g = r40[n1];  // g2(k+1:k+nt,n1)
                const int nd1kpt2 = 2 + __builtin_popcountll(s0 ^ g.w[0]) + __builtin_popcountll(s1 ^ g.w[1]) +
                                    __builtin_popcountll(s2 ^ g.w[2]);
                if (nd1kpt2 <= ntheta) {
                    Word ce, nx, e2;
                    for (int w = 0; w < 4; ++w) {
                        ce.w[w] = csub.w[w] ^ rows[n1].w[w];
                        nx.w[w] = ce.w[w] ^ hd.w[w];
                        e2.w[w] = nx.w[w] & par.w[w];
                    }
                    const float s = sum_set(0.0f, e2, k + 1);
                    const float dd = d1 + static_cast<float>(bit(nx, indices[n1] - 1)) * absrx[n1] + s;
                    if (dd < dmin) {
                        dmin = dd;
                        cw = ce;
                        nhardmin = popcount(nx);
                    }
                }
            }
        }
    }

    // 998: Re-order the codeword to [message bits][parity bits] format: cw
    // is in the generator's column order already.
    for (int j = 0; j < N; ++j) cw_out[j] = static_cast<int8_t>(bit(cw, j));
    std::memcpy(message91, cw_out, 91);
    if (crc_of_m96(cw_out) != 0) nhardmin = -nhardmin;
}

void osd174_91_reference(const float llr[174], int k, const int8_t apmask[174], int ndeep, int8_t message91[91],
                         int8_t cw_out[174], int& nhardmin, float& dmin) {
    if (k != K) throw std::invalid_argument("osd174_91: k must be 91");
    // ndeep 3 and more (and below 0, the Fortran's else branch) reach
    // boxit91/fetchit91 (npre2 = 1): never from WSJT-X's (174,91) callers.
    if (ndeep < 0 || ndeep >= 3) throw std::invalid_argument("osd174_91: only ndeep 0-2 carried");
    const auto& gen = osd_generator().gen;

    // 1-based working arrays, as the Fortran's.
    int8_t hdec[N + 1], cw[N + 1], ce[N + 1], c0[N + 1], apmaskr[N + 1], t8[N + 1];
    int8_t m0[K + 1], me[K + 1], mi[K + 1], misub[K + 1], temp[K + 1];
    int8_t e2sub[N - K + 1], e2[N - K + 1];
    int indices[N + 1], nxor[N + 1], indx[N + 1], ti[N + 1];
    float rx[N + 1], absrx[N + 1], tf[N + 1];
    int8_t genmrb[K][N];  // genmrb(r,c) = genmrb[r-1][c-1]; row r is g2(:,r)
    auto G = [&](int r, int c) -> int8_t& { return genmrb[r - 1][c - 1]; };

    for (int i = 1; i <= N; ++i) rx[i] = llr[i - 1];
    for (int i = 1; i <= N; ++i) apmaskr[i] = apmask[i - 1];

    // Hard decisions on the received word.
    for (int i = 1; i <= N; ++i) hdec[i] = rx[i] >= 0 ? 1 : 0;

    // Use magnitude of received symbols as a measure of reliability.
    for (int i = 1; i <= N; ++i) absrx[i] = std::fabs(rx[i]);
    wsprd::indexx(absrx + 1, N, indx + 1);

    // Re-order the columns of the generator matrix in order of decreasing reliability.
    for (int i = 1; i <= N; ++i) {
        for (int r = 1; r <= k; ++r) G(r, i) = gen[r - 1][indx[N + 1 - i] - 1];
        indices[i] = indx[N + 1 - i];
    }

    // Do gaussian elimination to create a generator matrix with the most reliable
    // received bits in positions 1:k in order of decreasing reliability (more or less).
    for (int id = 1; id <= k; ++id) {  // diagonal element indices
        for (int icol = id; icol <= k + 20; ++icol) {  // The 20 is ad hoc - beware
            if (G(id, icol) == 1) {
                if (icol != id) {  // reorder column
                    for (int r = 1; r <= k; ++r) temp[r] = G(r, id);
                    for (int r = 1; r <= k; ++r) G(r, id) = G(r, icol);
                    for (int r = 1; r <= k; ++r) G(r, icol) = temp[r];
                    const int itmp = indices[id];
                    indices[id] = indices[icol];
                    indices[icol] = itmp;
                }
                for (int ii = 1; ii <= k; ++ii)
                    if (ii != id && G(ii, id) == 1)
                        for (int c = 1; c <= N; ++c) G(ii, c) ^= G(id, c);
                break;
            }
        }
    }
    const int8_t* g2 = &genmrb[0][0];  // g2 = transpose(genmrb)

    // The hard decisions for the k MRB bits define the order 0 message, m0.
    for (int i = 1; i <= N; ++i) t8[i] = hdec[indices[i]];
    std::memcpy(hdec + 1, t8 + 1, N);  // hard decisions from received symbols
    std::memcpy(m0 + 1, hdec + 1, static_cast<size_t>(k));  // zero'th order message
    for (int i = 1; i <= N; ++i) tf[i] = absrx[indices[i]];
    std::memcpy(absrx + 1, tf + 1, sizeof(float) * N);
    for (int i = 1; i <= N; ++i) tf[i] = rx[indices[i]];
    std::memcpy(rx + 1, tf + 1, sizeof(float) * N);
    for (int i = 1; i <= N; ++i) t8[i] = apmaskr[indices[i]];
    std::memcpy(apmaskr + 1, t8 + 1, N);

    mrbencode91(m0, c0, g2, k);
    for (int i = 1; i <= N; ++i) nxor[i] = c0[i] ^ hdec[i];
    nhardmin = 0;
    for (int i = 1; i <= N; ++i) nhardmin += nxor[i];
    dmin = 0.0f;
    for (int i = 1; i <= N; ++i) dmin += static_cast<float>(nxor[i]) * absrx[i];
    std::memcpy(cw + 1, c0 + 1, N);

    if (ndeep != 0) {  // else goto 998 (order 0)
        int nord = 1, npre1 = 0, nt = 40, ntheta = 12;
        if (ndeep == 1) {
            nord = 1;
            npre1 = 0;
            nt = 40;
            ntheta = 12;
        } else if (ndeep == 2) {
            nord = 1;
            npre1 = 1;
            nt = 40;
            ntheta = 10;
        }

        const int np = N - k;
        for (int iorder = 1; iorder <= nord; ++iorder) {
            for (int i = 1; i <= k - iorder; ++i) misub[i] = 0;
            for (int i = k - iorder + 1; i <= k; ++i) misub[i] = 1;
            int iflag = k - iorder + 1;
            while (iflag >= 0) {
                const int iend = (iorder == nord && npre1 == 0) ? iflag : 1;
                float d1 = 0.0f;
                for (int n1 = iflag; n1 >= iend; --n1) {
                    std::memcpy(mi + 1, misub + 1, static_cast<size_t>(k));
                    mi[n1] = 1;
                    bool masked = false;
                    for (int i = 1; i <= k; ++i)
                        if ((apmaskr[i] & mi[i]) == 1) {
                            masked = true;
                            break;
                        }
                    if (masked) continue;
                    for (int i = 1; i <= k; ++i) me[i] = static_cast<int8_t>(m0[i] ^ mi[i]);
                    int nd1kpt;
                    if (n1 == iflag) {
                        mrbencode91(me, ce, g2, k);
                        for (int i = 1; i <= np; ++i) e2sub[i] = static_cast<int8_t>(ce[k + i] ^ hdec[k + i]);
                        std::memcpy(e2 + 1, e2sub + 1, static_cast<size_t>(np));
                        nd1kpt = 1;
                        for (int i = 1; i <= nt; ++i) nd1kpt += e2sub[i];
                        d1 = 0.0f;
                        for (int i = 1; i <= k; ++i) d1 += static_cast<float>(me[i] ^ hdec[i]) * absrx[i];
                    } else {
                        const int8_t* g = g2 + static_cast<size_t>(n1 - 1) * N;  // g2(k+1:N,n1)
                        for (int i = 1; i <= np; ++i) e2[i] = static_cast<int8_t>(e2sub[i] ^ g[k + i - 1]);
                        nd1kpt = 2;
                        for (int i = 1; i <= nt; ++i) nd1kpt += e2[i];
                    }
                    if (nd1kpt <= ntheta) {
                        mrbencode91(me, ce, g2, k);
                        for (int i = 1; i <= N; ++i) nxor[i] = ce[i] ^ hdec[i];
                        float dd;
                        if (n1 == iflag) {
                            float s = 0.0f;
                            for (int i = 1; i <= np; ++i) s += static_cast<float>(e2sub[i]) * absrx[k + i];
                            dd = d1 + s;
                        } else {
                            float s = 0.0f;
                            for (int i = 1; i <= np; ++i) s += static_cast<float>(e2[i]) * absrx[k + i];
                            dd = d1 + static_cast<float>(ce[n1] ^ hdec[n1]) * absrx[n1] + s;
                        }
                        if (dd < dmin) {
                            dmin = dd;
                            std::memcpy(cw + 1, ce + 1, N);
                            nhardmin = 0;
                            for (int i = 1; i <= N; ++i) nhardmin += nxor[i];
                        }
                    }
                }
                // Get the next test error pattern, iflag will go negative
                // when the last pattern with weight iorder has been generated.
                nextpat91(misub, k, iorder, iflag);
            }
        }
    }

    // 998: Re-order the codeword to [message bits][parity bits] format.
    for (int i = 1; i <= N; ++i) ti[i] = cw[i];
    for (int i = 1; i <= N; ++i) cw[indices[i]] = static_cast<int8_t>(ti[i]);
    std::memcpy(cw_out, cw + 1, N);
    std::memcpy(message91, cw + 1, 91);
    if (crc_of_m96(cw + 1) != 0) nhardmin = -nhardmin;
}

// The line-for-line decode174_91.f90 decode174_91 was optimised from, kept
// for the side-by-side test and as its fallback.
void decode174_91_reference(const float llr[174], int Keff, int maxosd, int norder, const int8_t apmask[174],
                            int8_t message91[91], int8_t cw[174], int& ntype, int& nharderror, float& dmin) {
    const int maxiterations = 30;
    float zsave[3][N];  // the Fortran's local zsave(N,3)
    int nosd = 0;
    if (maxosd > 3) maxosd = 3;
    if (maxosd == 0) {  // osd with channel llrs
        nosd = 1;
        std::memcpy(zsave[0], llr, sizeof(float) * N);
    } else if (maxosd > 0) {
        nosd = maxosd;
    } else {  // just bp
        nosd = 0;
    }

    float toc[M][7], tov[N][3], tanhtoc[M][7], zn[N], zsum[N];
    std::memset(toc, 0, sizeof toc);
    std::memset(tov, 0, sizeof tov);
    std::memset(tanhtoc, 0, sizeof tanhtoc);
    // initialize messages to checks
    for (int j = 0; j < M; ++j)
        for (int i = 0; i < kNrw[j]; ++i) toc[j][i] = llr[kNm[j][i] - 1];

    int ncnt = 0, nclast = 0;
    for (int i = 0; i < N; ++i) zsum[i] = 0.0f;
    for (int iter = 0; iter <= maxiterations; ++iter) {
        // Update bit log likelihood ratios (tov=0 in iteration 0).
        for (int i = 0; i < N; ++i) {
            if (apmask[i] != 1) {
                float s = 0.0f;
                for (int kk = 0; kk < kNcw; ++kk) s += tov[i][kk];
                zn[i] = llr[i] + s;
            } else {
                zn[i] = llr[i];
            }
        }
        for (int i = 0; i < N; ++i) zsum[i] = zsum[i] + zn[i];
        if (iter > 0 && iter <= maxosd) std::memcpy(zsave[iter - 1], zsum, sizeof zsum);

        // Check to see if we have a codeword (check before we do any iteration).
        for (int i = 0; i < N; ++i) cw[i] = zn[i] > 0.0f ? 1 : 0;
        int ncheck = 0;
        for (int i = 0; i < M; ++i) {
            int synd = 0;
            for (int r = 0; r < kNrw[i]; ++r) synd += cw[kNm[i][r] - 1];
            if (synd % 2 != 0) ++ncheck;
        }
        if (ncheck == 0) {  // we have a codeword - if crc is good, return it
            const int nbadcrc = crc_of_m96(cw);
            int nh = 0;  // count((2*cw-1)*llr .lt. 0.0)
            for (int i = 0; i < N; ++i)
                if (static_cast<float>(2 * cw[i] - 1) * llr[i] < 0.0f) ++nh;
            nharderror = nh;
            if (nbadcrc == 0) {
                std::memcpy(message91, cw, 91);
                float d = 0.0f;
                for (int i = 0; i < N; ++i) {
                    const int8_t hdec = llr[i] >= 0 ? 1 : 0;
                    const int8_t nxor = static_cast<int8_t>(hdec ^ cw[i]);
                    d += static_cast<float>(nxor) * std::fabs(llr[i]);
                }
                dmin = d;
                ntype = 1;
                return;
            }
        }

        if (iter > 0) {  // an early stopping criterion
            const int nd = ncheck - nclast;
            if (nd < 0) ncnt = 0;  // # of unsatisfied parity checks decreased
            else ++ncnt;
            if (ncnt >= 5 && iter >= 10 && ncheck > 15) {
                nharderror = -1;
                break;
            }
        }
        nclast = ncheck;

        // Send messages from bits to check nodes
        for (int j = 0; j < M; ++j) {
            for (int i = 0; i < kNrw[j]; ++i) {
                const int ibj = kNm[j][i];
                toc[j][i] = zn[ibj - 1];
                for (int kk = 0; kk < kNcw; ++kk)  // subtract off what the bit had received from the check
                    if (kMn[ibj - 1][kk] == j + 1) toc[j][i] = toc[j][i] - tov[ibj - 1][kk];
            }
        }

        // send messages from check nodes to variable nodes (all seven
        // slots, as the Fortran: a six-bit check's seventh is tanh(-0))
        for (int i = 0; i < M; ++i)
            for (int r = 0; r < 7; ++r) tanhtoc[i][r] = tanhf_glibc(-toc[i][r] / 2);

        for (int j = 0; j < N; ++j) {
            for (int i = 0; i < kNcw; ++i) {
                const int ichk = kMn[j][i];  // Mn(:,j) are the checks that include bit j
                float Tmn = 1.0f;
                for (int r = 0; r < kNrw[ichk - 1]; ++r)
                    if (kNm[ichk - 1][r] != j + 1) Tmn = Tmn * tanhtoc[ichk - 1][r];
                const float y = platanh(-Tmn);
                tov[j][i] = 2 * y;
            }
        }
    }  // bp iterations

    for (int i = 1; i <= nosd; ++i) {
        std::memcpy(zn, zsave[i - 1], sizeof zn);
        float dminosd = 0.0f;
        osd174_91(zn, Keff, apmask, norder, message91, cw, nharderror, dminosd);
        if (nharderror > 0) {
            float d = 0.0f;
            for (int r = 0; r < N; ++r) {
                const int8_t hdec = llr[r] >= 0 ? 1 : 0;
                const int8_t nxor = static_cast<int8_t>(hdec ^ cw[r]);
                d += static_cast<float>(nxor) * std::fabs(llr[r]);
            }
            dmin = d;  // against the channel llrs; nharderror stays the OSD's
            ntype = 2;
            return;
        }
    }

    ntype = 0;
    nharderror = -1;
}

}  // namespace ft8
