// osdwspr.f90 and indexx.f90 of WSJT-X (lib/wsprd, lib): see wsprd_osd.h.
// osdwspr.f90: Steven Franke, K9AN, and Joe Taylor, K1JT, WSJT-X (GPLv3),
// and so is this.
//
// The Fortran's control flow is followed statement for statement so the two
// read side by side, with indices 1-based as there. Its 0/1 integer*1 vectors
// are bits (Bits below; messages a uint64_t, bit i-1 for element i), so an
// encode is a few XORs of g2's rows and a count a popcount. What had to
// change, and nothing else:
//   - the SAVEd generator matrix is a function-local static built once;
//   - genmrb is kept as its rows, g2's (g2[k] is g2(:,k)), as the Gaussian
//     elimination only XORs rows and swaps columns;
//   - mrbencode(me) is c0 ^ mrbencode(mi), me = m0 ^ mi (XOR is linear);
//   - boxit/fetchit's COMMON /boxes/ and fetchit's SAVEd lastpat/inext are a
//     caller's OsdTables, reset as the reference resets them (each call that
//     uses them starts with a reset) where anything reads them; fp covers
//     the 2^16 patterns ntau=16 can make, not 0:525000;
//   - the counters nothing reads (nrejected, ncount2, ntotal2, nd1Kptbest,
//     and ntotal outside boxit's index) are left out;
//   - indexx's STOP 'NSTACK too small' (unreachable for n=162) aborts.
// The float sums run in the Fortran's order, single precision, left to
// right (gfortran reorders none of them without -ffast-math); every product
// in them is a 0/1 times a magnitude, so no FMA contraction can change them,
// and leaving out the zero terms changes nothing either (add_set).
// For speed, with the same results: without an apmask the order 1..nord
// search visits only the n1 whose nd1Kpt, a popcount of g2's parity words,
// passes ntheta (with AVX2 or NEON, four n1 at a time, all found before the
// first visit); the generator's columns are reordered from their bits; the
// second pre-processing rule looks fp's slots up as bits, a pattern and its
// ntau one-flip neighbours at once, and calls fetchit only on those that are
// set; and its distance stops adding once it reaches dmin.
#include "wsprd_osd.h"

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <vector>

#if defined(__AVX2__) && (defined(__x86_64__) || defined(_M_X64))
#include <immintrin.h>
#define WSPRD_OSD_AVX2 1
#elif defined(__aarch64__) && defined(__ARM_NEON)
#include <arm_neon.h>
#define WSPRD_OSD_NEON 1
#endif

namespace wsprd {

namespace {

constexpr int N = 162, K = 50, L = 32;
using i1 = signed char;  // integer*1

const i1 gg[64] = {1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 0, 0, 1, 1, 0, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0,
                   0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1};

struct Gen {
    i1 g[N + 1][K + 1] = {};  // g[j][i] = gen(i,j)
    Gen() {
        for (int j = 1; j <= 2 * L; j++) g[j][1] = gg[j - 1];
        for (int i = 2; i <= K; i++)  // gen(i,:)=cshift(gen(i-1,:),-2)
            for (int j = 1; j <= N; j++) g[j][i] = g[(j - 1 - 2 + N) % N + 1][i - 1];
    }
};

const Gen& gen() {
    static const Gen g;
    return g;
}

// gen(:,j) as bits (bit k-1 for gen(k,j)), for the column reorder.
struct GenCols {
    uint64_t c[N + 1] = {};
    GenCols() {
        const Gen& G = gen();
        for (int j = 1; j <= N; j++)
            for (int k = 1; k <= K; k++)
                if (G.g[j][k] == 1) c[j] |= uint64_t{1} << (k - 1);
    }
};

const GenCols& gen_cols() {
    static const GenCols g;
    return g;
}

// A codeword (or a row of g2, or any 162 0/1 values in the reordered
// positions 1..N) as bits: w[0] holds positions 1..K (bit i-1 for position
// i), w[1] and w[2] the parity positions K+1..N (bit q%64 of w[1+q/64] for
// position K+1+q). Positions run upwards through w[0], w[1], w[2], so walking
// the set bits in that order visits them as the Fortran's loops do.
struct Bits {
    uint64_t w[3];
    Bits operator^(const Bits& o) const { return {{w[0] ^ o.w[0], w[1] ^ o.w[1], w[2] ^ o.w[2]}}; }
    Bits& operator^=(const Bits& o) { return *this = *this ^ o; }
    bool get(int p) const { return p <= K ? (w[0] >> (p - 1)) & 1 : (w[1 + (p - K - 1) / 64] >> ((p - K - 1) % 64)) & 1; }
    void flip(int p) {
        if (p <= K)
            w[0] ^= uint64_t{1} << (p - 1);
        else
            w[1 + (p - K - 1) / 64] ^= uint64_t{1} << ((p - K - 1) % 64);
    }
    int count() const { return __builtin_popcountll(w[0]) + __builtin_popcountll(w[1]) + __builtin_popcountll(w[2]); }
};

// s plus a[p] for each set bit p of x, in ascending p; a is 1-based and base
// is the position of x's bit 0. The Fortran adds (0|1)*a(i) for every i: a
// zero term leaves the sum as it was (the sums start at +0 and a is finite
// and not negative, contracted or not: fma(0,a,s) = s, fma(1,a,s) = s+a), so
// adding only the ones gives the same bits.
inline float add_set(float s, uint64_t x, const float* a, int base) {
    while (x) {
        s += a[base + __builtin_ctzll(x)];
        x &= x - 1;
    }
    return s;
}

inline float add_all(float s, const Bits& b, const float* a) {
    s = add_set(s, b.w[0], a, 1);
    s = add_set(s, b.w[1], a, K + 1);
    return add_set(s, b.w[2], a, K + 65);
}

// add_all(0, b, a) < dmin, the sum in dd when it is: as the terms are not
// negative and rounding is monotonic, a partial sum >= dmin makes the whole
// one >= dmin, so the rest is not added (a NaN never compares >= and is added
// through to fail the last test, as it would have).
inline bool add_all_below(const Bits& b, const float* a, float dmin, float& dd) {
    float s = add_set(0.f, b.w[0], a, 1);
    if (s >= dmin) return false;
    s = add_set(s, b.w[1], a, K + 1);
    if (s >= dmin) return false;
    dd = add_set(s, b.w[2], a, K + 65);
    return dd < dmin;
}

inline float add_parity(float s, const Bits& b, const float* a) {
    s = add_set(s, b.w[1], a, K + 1);
    return add_set(s, b.w[2], a, K + 65);
}

// codeword for message bits mi given c0, the codeword of m0: the XOR of
// g2's rows for m0^mi's 1s, which (XOR being linear) is c0 with the rows for
// mi's 1s added.
inline Bits mrbencode(const Bits& c0, uint64_t mi, const Bits* g2) {
    Bits c = c0;
    while (mi) {
        c ^= g2[__builtin_ctzll(mi) + 1];
        mi &= mi - 1;
    }
    return c;
}

// nextpat on mi's bits (bit i-1 for mi(i)): the Fortran's next pattern of
// iorder ones among k, by moving the highest 1 that has a 0 below it up one
// place and packing the 1s above it at the top.
void nextpat(uint64_t& mi, int k, int iorder, int& iflag) {
    uint64_t t = ~mi & (mi >> 1) & ((uint64_t{1} << (k - 1)) - 1);  // bit i-1: mi(i)=0, mi(i+1)=1
    if (!t) {  // no more patterns of this order
        iflag = -1;
        return;
    }
    int ind = 64 - __builtin_clzll(t);
    uint64_t ms = (mi & ((uint64_t{1} << (ind - 1)) - 1)) | (uint64_t{1} << (ind - 1));
    if (ind + 1 < k) {
        int nz = iorder - __builtin_popcountll(ms);
        if (nz > 0) ms |= ((uint64_t{1} << nz) - 1) << (k - nz);
    }
    mi = ms;
    iflag = __builtin_ctzll(mi) + 1;  // iflag will point to the lowest-index 1 in mi
}

// COMMON /boxes/ and fetchit's SAVEd state. A pattern is e2(1:ntau) as the
// bits of w[1] (e2(i) in bit i-1), not the Fortran's ipat (e2(i) in bit
// ntau-i): a one to one relabelling of the keys, so the same lists come back.
struct Boxes {
    int indexes[4000 + 1][2 + 1];
    std::vector<int> fp = std::vector<int>(1 << 16, -1);
    std::vector<int> used;  // the fp slots set since the last reset
    uint64_t set[(1 << 16) / 64] = {};  // a bit per fp slot that is not -1
    int np[4000 + 1];
    int lastpat = -1, inext = -1;

    // The reset is the reference's (fp, np, indexes all -1) for what anything
    // reads: fp's slots set last time, and np up to npindex as it grows.
    void boxit(bool& reset, int ipat, int npindex, int i1v, int i2v) {
        if (reset) {
            for (int k : used) {
                fp[k] = -1;
                set[k >> 6] = 0;
            }
            used.clear();
            reset = false;
        }
        np[npindex] = -1;
        indexes[npindex][1] = i1v;
        indexes[npindex][2] = i2v;
        int ip = fp[ipat];  // see what's currently stored in fp(ipat)
        if (ip == -1) {
            fp[ipat] = npindex;
            set[ipat >> 6] |= uint64_t{1} << (ipat & 63);
            used.push_back(ipat);
        } else {
            while (np[ip] != -1) ip = np[ip];
            np[ip] = npindex;
        }
    }

    void fetch_reset(bool& reset) {
        if (reset) {
            lastpat = -1;
            reset = false;
        }
    }

    void fetchit(bool& reset, int ipat, int& i1v, int& i2v) {
        fetch_reset(reset);
        int index = fp[ipat];
        if (lastpat != ipat && index > 0) {  // return first set of indices
            i1v = indexes[index][1];
            i2v = indexes[index][2];
            inext = np[index];
        } else if (lastpat == ipat && inext > 0) {
            i1v = indexes[inext][1];
            i2v = indexes[inext][2];
            inext = np[inext];
        } else {
            i1v = -1;
            i2v = -1;
            inext = -1;
        }
        lastpat = ipat;
    }

    // fetchit(reset, ipat, ...) for an ipat whose fp slot is -1, once reset
    // is false: the first test fails (index -1), and so does the second, as
    // inext > 0 only after a hit on lastpat. Most fetches are these, and
    // miss() is their effect. ball() says which of ipat and the ntau (6..16)
    // patterns one bit flip away have a slot set: bit 0 for ipat, bit i for
    // ipat with bit i-1 flipped. Flips of bits 0..5 stay in ipat's word of
    // set, the others move to another word at the same bit.
    uint32_t ball(int ipat, int ntau) const {
        int w = ipat >> 6, b = ipat & 63;
        uint64_t x = set[w];
        uint32_t m = (x >> b) & 1;
        for (int i = 0; i < 6; i++) m |= static_cast<uint32_t>((x >> (b ^ 1 << i)) & 1) << (i + 1);
        for (int i = 6; i < ntau; i++) m |= static_cast<uint32_t>((set[w ^ (1 << (i - 6))] >> b) & 1) << (i + 1);
        return m;
    }
    void miss(int ipat) {
        inext = -1;
        lastpat = ipat;
    }
};

#ifdef WSPRD_OSD_AVX2
// Bit n1-1 set for each n1 in 1..n (n <= 49) whose popcount(s1 ^ p1[n1]) +
// popcount(s2 ^ p2[n1]) is at most lim, four n1 to a vector (a nibble
// lookup per byte, summed per 64-bit lane); p1 and p2 hold n1's words at
// [n1-1] and reach at least 4 entries past n.
inline uint64_t near_rows(uint64_t s1, uint64_t s2, const uint64_t* p1, const uint64_t* p2, int n, int lim) {
    const __m256i lut = _mm256_setr_epi8(0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4, 0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3,
                                         2, 3, 3, 4);
    const __m256i lo = _mm256_set1_epi8(0x0f), zero = _mm256_setzero_si256();
    const __m256i v1 = _mm256_set1_epi64x(static_cast<long long>(s1)), v2 = _mm256_set1_epi64x(static_cast<long long>(s2));
    const __m256i top = _mm256_set1_epi64x(lim + 1);
    uint64_t m = 0;
    for (int b = 0; b < n; b += 4) {
        __m256i x = _mm256_xor_si256(v1, _mm256_loadu_si256(reinterpret_cast<const __m256i*>(p1 + b)));
        __m256i y = _mm256_xor_si256(v2, _mm256_loadu_si256(reinterpret_cast<const __m256i*>(p2 + b)));  // < 4: low nibble only
        __m256i c = _mm256_add_epi8(_mm256_shuffle_epi8(lut, _mm256_and_si256(x, lo)),
                                    _mm256_shuffle_epi8(lut, _mm256_and_si256(_mm256_srli_epi16(x, 4), lo)));
        c = _mm256_sad_epu8(_mm256_add_epi8(c, _mm256_shuffle_epi8(lut, y)), zero);
        m |= static_cast<uint64_t>(_mm256_movemask_pd(_mm256_castsi256_pd(_mm256_cmpgt_epi64(top, c)))) << b;
    }
    return m & ((uint64_t{1} << n) - 1);
}
#elif defined(WSPRD_OSD_NEON)
// As the AVX2 near_rows: a byte count per byte (cnt), summed per 64-bit
// lane, two lanes to a vector; four n1's counts narrowed to 32 bits are
// compared at once and their bits gathered (1, 2, 4, 8) by an add across.
inline uint64_t near_rows(uint64_t s1, uint64_t s2, const uint64_t* p1, const uint64_t* p2, int n, int lim) {
    const uint64x2_t v1 = vdupq_n_u64(s1), v2 = vdupq_n_u64(s2);
    const int32x4_t top = vdupq_n_s32(lim + 1);
    static const uint32_t kBit[4] = {1, 2, 4, 8};
    const uint32x4_t bit = vld1q_u32(kBit);
    uint64_t m = 0;
    for (int b = 0; b < n; b += 4) {
        uint64x2_t c[2];
        for (int h = 0; h < 2; h++) {
            uint8x16_t x = vreinterpretq_u8_u64(veorq_u64(v1, vld1q_u64(p1 + b + 2 * h)));
            uint8x16_t y = vreinterpretq_u8_u64(veorq_u64(v2, vld1q_u64(p2 + b + 2 * h)));
            c[h] = vpaddlq_u32(vpaddlq_u16(vpaddlq_u8(vaddq_u8(vcntq_u8(x), vcntq_u8(y)))));
        }
        int32x4_t c4 = vreinterpretq_s32_u32(vcombine_u32(vmovn_u64(c[0]), vmovn_u64(c[1])));
        m |= static_cast<uint64_t>(vaddvq_u32(vandq_u32(vcltq_s32(c4, top), bit))) << b;
    }
    return m & ((uint64_t{1} << n) - 1);
}
#endif

}  // namespace

void indexx(const float* arr0, int n, int* indx0) {
    constexpr int M = 7, NSTACK = 50;
    const float* arr = arr0 - 1;  // 1-based
    int* indx = indx0 - 1;
    int i, indxt, ir, itemp, j, jstack, k, l, istack[NSTACK + 1];
    float a;

    for (j = 1; j <= n; j++) indx[j] = j;

    jstack = 0;
    l = 1;
    ir = n;
    for (;;) {
        if (ir - l < M) {
            for (j = l + 1; j <= ir; j++) {
                indxt = indx[j];
                a = arr[indxt];
                for (i = j - 1; i >= 1; i--) {
                    if (arr[indx[i]] <= a) goto l2;
                    indx[i + 1] = indx[i];
                }
                i = 0;
            l2:
                indx[i + 1] = indxt;
            }
            if (jstack == 0) return;

            ir = istack[jstack];
            l = istack[jstack - 1];
            jstack = jstack - 2;
        } else {
            k = (l + ir) / 2;
            itemp = indx[k];
            indx[k] = indx[l + 1];
            indx[l + 1] = itemp;

            if (arr[indx[l + 1]] > arr[indx[ir]]) {
                itemp = indx[l + 1];
                indx[l + 1] = indx[ir];
                indx[ir] = itemp;
            }

            if (arr[indx[l]] > arr[indx[ir]]) {
                itemp = indx[l];
                indx[l] = indx[ir];
                indx[ir] = itemp;
            }

            if (arr[indx[l + 1]] > arr[indx[l]]) {
                itemp = indx[l + 1];
                indx[l + 1] = indx[l];
                indx[l] = itemp;
            }

            i = l + 1;
            j = ir;
            indxt = indx[l];
            a = arr[indxt];
            for (;;) {
                do i = i + 1;
                while (arr[indx[i]] < a);
                do j = j - 1;
                while (arr[indx[j]] > a);
                if (j < i) break;
                itemp = indx[i];
                indx[i] = indx[j];
                indx[j] = itemp;
            }

            indx[l] = indx[j];
            indx[j] = indxt;
            jstack = jstack + 2;
            if (jstack > NSTACK) std::abort();  // 'NSTACK too small in indexx': never for n=162
            if (ir - i + 1 >= j - l) {
                istack[jstack] = ir;
                istack[jstack - 1] = i;
                ir = j - 1;
            } else {
                istack[jstack] = j - 1;
                istack[jstack - 1] = l;
                l = i;
            }
        }
    }
}

struct OsdBoxes : Boxes {};
void OsdBoxesFree::operator()(OsdBoxes* b) const { delete b; }

void osdwspr(const float ss[162], const unsigned char apmask[162], int ndeep, unsigned char cw_out[162], int* nhardmin_out,
             float* dmin_out, OsdTables* tables) {
    const GenCols& GC = gen_cols();
    Bits g2[K + 1] = {};  // g2[k] = g2(:,k), the rows of genmrb
    Bits hdec = {}, cw, c0, ce;
    uint64_t apmaskr = 0, m0, mi, misub;
    int indices[N + 1];
    int indx[N + 1], nhardmin;
    float rx[N + 1], absrx[N + 1], dmin;
    bool reset;
    int nord, npre1, npre2, nt, ntheta, ntau = 0;
    int iflag, iend, in1, in2, nd1Kpt;
    float d1, dd;

    for (int i = 1; i <= N; i++) rx[i] = ss[i - 1] / 127.0f;

    // Use magnitude of received symbols as a measure of reliability.
    for (int i = 1; i <= N; i++) absrx[i] = std::fabs(rx[i]);
    indexx(absrx + 1, N, indx + 1);

    // Re-order the columns of the generator matrix in order of decreasing reliability.
    for (int i = 1; i <= N; i++) {
        indices[i] = indx[N + 1 - i];
        for (uint64_t c = GC.c[indices[i]]; c; c &= c - 1) g2[__builtin_ctzll(c) + 1].flip(i);
    }

    // Do gaussian elimination to create a generator matrix with the most reliable
    // received bits in positions 1:K in order of decreasing reliability (more or less).
    for (int id = 1; id <= K; id++) {         // diagonal element indices
        for (int icol = id; icol <= K + 20; icol++) {  // The 20 is ad hoc - beware
            if (g2[id].get(icol)) {
                if (icol != id) {  // reorder column
                    for (int k = 1; k <= K; k++)
                        if (g2[k].get(id) != g2[k].get(icol)) {
                            g2[k].flip(id);
                            g2[k].flip(icol);
                        }
                    int itmp = indices[id];
                    indices[id] = indices[icol];
                    indices[icol] = itmp;
                }
                for (int ii = 1; ii <= K; ii++)
                    if (ii != id && g2[ii].get(id)) g2[ii] ^= g2[id];
                break;
            }
        }
    }

    // The hard decisions for the K MRB bits define the order 0 message, m0.
    // Encode m0 using the modified generator matrix to find the "order 0" codeword.
    // Flip various combinations of bits in m0 and re-encode to generate a list of
    // codewords. Return the member of the list that has the smallest Euclidean
    // distance to the received word.
    {
        float tf[N + 1];
        for (int i = 1; i <= N; i++) {
            if (rx[indices[i]] >= 0) hdec.flip(i);  // hard decisions from received symbols
            tf[i] = absrx[indices[i]];
        }
        std::memcpy(absrx, tf, sizeof absrx);
        for (int i = 1; i <= K; i++)  // apmask past K is never read
            if (apmask[indices[i] - 1] & 1) apmaskr |= uint64_t{1} << (i - 1);
    }
    m0 = hdec.w[0];  // zero'th order message

    c0 = mrbencode(Bits{}, m0, g2);

    {
        Bits nxor = c0 ^ hdec;
        nhardmin = nxor.count();
        dmin = add_all(0.f, nxor, absrx);
    }
    cw = c0;

    if (ndeep <= 0) goto l998;  // norder=0
    if (ndeep > 5) ndeep = 5;
    if (ndeep == 1) {
        nord = 1;
        npre1 = 0;
        npre2 = 0;
        nt = 66;
        ntheta = 16;
    } else if (ndeep == 2) {
        nord = 1;
        npre1 = 1;
        npre2 = 0;
        nt = 66;
        ntheta = 22;
        ntau = 16;
    } else if (ndeep == 3) {
        nord = 1;
        npre1 = 1;
        npre2 = 1;
        nt = 66;
        ntheta = 22;
        ntau = 16;
    } else if (ndeep == 4) {
        nord = 2;
        npre1 = 1;
        npre2 = 1;
        nt = 66;
        ntheta = 22;
        ntau = 16;
    } else if (ndeep == 5) {
        nord = 3;
        npre1 = 1;
        npre2 = 1;
        nt = 66;
        ntheta = 22;
        ntau = 16;
    } else {
        goto l998;
    }

    {
        // sum(e2(1:nt)) for nt=66: parity positions K+1..K+66.
        const uint64_t ntmask = (uint64_t{1} << (nt - 64)) - 1;
        auto ones_nt = [ntmask](const Bits& e) { return __builtin_popcountll(e.w[1]) + __builtin_popcountll(e.w[2] & ntmask); };
        Bits e2sub = {}, e2;
        uint64_t p1[K + 4] = {}, p2[K + 4] = {};  // g2's rows' parity words on 1:nt, row k at [k-1]
        for (int k = 1; k <= K; k++) {
            p1[k - 1] = g2[k].w[1];
            p2[k - 1] = g2[k].w[2] & ntmask;
        }
        for (int iorder = 1; iorder <= nord; iorder++) {
            misub = ((uint64_t{1} << iorder) - 1) << (K - iorder);
            iflag = K - iorder + 1;
            while (iflag >= 0) {
                if (iorder == nord && npre1 == 0)
                    iend = iflag;
                else
                    iend = 1;
                d1 = 0.f;
                auto visit = [&](int n1) {
                    mi = misub | uint64_t{1} << (n1 - 1);
                    if (apmaskr & mi) return;
                    if (n1 == iflag) {
                        ce = mrbencode(c0, mi, g2);
                        e2sub = ce ^ hdec;  // only its parity words are e2sub
                        e2 = e2sub;
                        nd1Kpt = ones_nt(e2sub) + 1;
                        d1 = add_set(0.f, (m0 ^ mi) ^ hdec.w[0], absrx, 1);
                    } else {
                        e2 = e2sub ^ g2[n1];
                        nd1Kpt = ones_nt(e2) + 2;
                    }
                    if (nd1Kpt <= ntheta) {
                        ce = mrbencode(c0, mi, g2);
                        Bits nxor = ce ^ hdec;
                        if (n1 == iflag) {
                            dd = d1 + add_parity(0.f, e2sub, absrx);
                        } else {
                            float s = add_parity(0.f, e2, absrx);
                            float t = d1;
                            if (nxor.get(n1)) t += absrx[n1];
                            dd = t + s;
                        }
                        if (dd < dmin) {
                            dmin = dd;
                            cw = ce;
                            nhardmin = nxor.count();
                        }
                    }
                };
                if (apmaskr) {
                    for (int n1 = iflag; n1 >= iend; n1--) visit(n1);
                } else {
                    // Without an apmask, an n1 below iflag does something
                    // only when its nd1Kpt (e2sub ^ g2[n1] on 1:nt, plus 2)
                    // is at most ntheta: find those from the parity words.
                    // (Those visits leave e2sub as it was, so all are found first.)
                    visit(iflag);
                    const uint64_t s1 = e2sub.w[1], s2 = e2sub.w[2] & ntmask;
#if defined(WSPRD_OSD_AVX2) || defined(WSPRD_OSD_NEON)
                    if (iend < iflag) {  // iend is 1 here
                        uint64_t m = near_rows(s1, s2, p1, p2, iflag - 1, ntheta - 2);
                        while (m) {
                            int n1 = 64 - __builtin_clzll(m);
                            m ^= uint64_t{1} << (n1 - 1);
                            visit(n1);
                        }
                    }
#else
                    for (int n1 = iflag - 1; n1 >= iend; n1--)
                        if (__builtin_popcountll(s1 ^ p1[n1 - 1]) + __builtin_popcountll(s2 ^ p2[n1 - 1]) + 2 <= ntheta) visit(n1);
#endif
                }
                // Get the next test error pattern, iflag will go negative
                // when the last pattern with weight iorder has been generated.
                nextpat(misub, K, iorder, iflag);
            }
        }
    }

    if (npre2 == 1) {
        OsdTables own;
        if (!tables) tables = &own;
        if (!*tables) tables->reset(new OsdBoxes);
        Boxes& boxes = **tables;
        const uint64_t taumask = (uint64_t{1} << ntau) - 1;
        reset = true;
        int ntotal = 0;
        for (int i1v = K; i1v >= 1; i1v--) {
            for (int i2v = i1v - 1; i2v >= 1; i2v--) {
                ntotal = ntotal + 1;
                boxes.boxit(reset, static_cast<int>((g2[i1v].w[1] ^ g2[i2v].w[1]) & taumask), ntotal, i1v, i2v);
            }
        }

        reset = true;
        // Now run through again and do the second pre-processing rule
        misub = ((uint64_t{1} << nord) - 1) << (K - nord);
        iflag = K - nord + 1;
        boxes.fetch_reset(reset);  // as the first fetchit would, so the tests below may skip it
        while (iflag >= 0) {
            ce = mrbencode(c0, misub, g2);
            uint64_t e2sub = (ce.w[1] ^ hdec.w[1]) & taumask;
            auto r2pat_of = [e2sub](int i2) { return static_cast<int>(i2 > 0 ? e2sub ^ uint64_t{1} << (i2 - 1) : e2sub); };
            // Only the i2 whose r2pat has a slot set are fetched; a run of
            // the others leaves fetchit's state as its last miss would.
            int last = -1;  // the last i2 fetched
            for (uint32_t ball = boxes.ball(static_cast<int>(e2sub), ntau); ball; ball &= ball - 1) {
                int i2 = __builtin_ctz(ball);
                if (i2 > last + 1) boxes.miss(r2pat_of(i2 - 1));
                last = i2;
                int r2pat = r2pat_of(i2);
                for (;;) {  // 778
                    boxes.fetchit(reset, r2pat, in1, in2);
                    if (!(in1 > 0 && in2 > 0)) break;
                    mi = misub | uint64_t{1} << (in1 - 1) | uint64_t{1} << (in2 - 1);
                    if (__builtin_popcountll(mi) < nord + npre1 + npre2 || (apmaskr & mi)) break;  // cycle: on to the next i2
                    ce = mrbencode(c0, mi, g2);
                    Bits nxor = ce ^ hdec;
                    if (add_all_below(nxor, absrx, dmin, dd)) {
                        dmin = dd;
                        cw = ce;
                        nhardmin = nxor.count();
                    }
                }
            }
            if (last < ntau) boxes.miss(r2pat_of(ntau));
            nextpat(misub, K, nord, iflag);
        }
    }

l998:
    // Re-order the codeword to as-received order.
    for (int i = 1; i <= N; i++) cw_out[indices[i] - 1] = static_cast<unsigned char>(cw.get(i));
    *nhardmin_out = nhardmin;
    *dmin_out = dmin;
}

}  // namespace wsprd
