// FT4's front end (ft4_frontend.h): lib/ft4/get_ft4_bitmetrics.f90 at
// 567ad29ce, line for line, and the part of lib/ft4_decode.f90 that reads
// its result (the hard-sync count nsync_qual and the llra, llrb, llrc
// mapping with scalefac 2.83).
//
// 567ad29ce's get_ft4_bitmetrics has three metric columns (coherent
// sequences of 1, 2 and 4 symbols), no imetric argument, no fourth
// (bmetd-like) or fifth (the >= cherry-pick) column: those are WSJT-X
// improved's additions and are not here.
//
// As its Release object does it (from its disassembly): each symbol's
// 32-point forward c2c FFT as four2a plans it (fst4w::fft_c2c_inplace on a
// 16-byte or better aligned buffer, as the reference's csymb on its stack:
// the trace's '@ bm-csymb' line), bins 0..3 kept unscaled (no /1e3 as
// FT8's ft8b has); abs() as cabsf (cabs_, hypotf_glibc) at every site; the
// nsym 2 and 4 sums of complex terms left to right; gfortran's maxloc
// (first maximum, NaNs skipped) and masked maxval (ft8::fortran's
// maxloc, maxval_bit); then FT8's normalizebmet on each column (the same
// object at 567ad29ce).
//
// Joe Taylor, K1JT, Steven Franke, K9AN, and the WSJT-X developers (GNU GPL
// v3), and so is this.
#include "ft4_frontend.h"

#include <cmath>
#include <cstdint>

#if defined(__AVX2__) && (defined(__x86_64__) || defined(_M_X64))
#include <immintrin.h>
#define FT4_METRICS_AVX2 1
#elif defined(__aarch64__) && defined(__ARM_NEON)
#include <arm_neon.h>
#define FT4_METRICS_NEON 1
#endif

#include "fst4w_frontend.h"  // fft_c2c_inplace, FftwVec
#include "ft8_fortran.h"     // the Fortran helpers, FT8's
#include "ft8_frontend.h"    // normalizebmet
#include "ft8_maths.h"       // glibc 2.39's maths, carried

namespace ft4 {

using namespace ft8::fortran;

// This group's state. The reference SAVEs only `first` and the `one`
// table (pure functions of the index: maxval_bit tests the bit itself);
// csymb is a stack local there, kept here so that it is FFTW-allocated
// (aligned as the plan the cache holds), and al the alignment the trace's
// '@ bm-csymb' line last reported (as the instrumented reference writes
// it when loc(csymb) mod 16 changes).
struct Frontend::Metrics {
    fst4w::FftwVec<cfloat> csymb = fst4w::FftwVec<cfloat>(NSS);
    int al = -1;
};

Frontend::Metrics* Frontend::make_metrics() { return new Metrics(); }
void Frontend::del_metrics(Metrics* p) { delete p; }
size_t Frontend::metrics_bytes(const Metrics* p) { return p ? p->csymb.capacity() * sizeof(cfloat) : 0; }

namespace {

// For each bit b of the index, the maximum of s2(0:255) over the indices
// with that bit set (mx1[b]) and clear (mx0[b]), as maxval_bit gives them
// when no s2 is a NaN (then any order finds the same float: the values
// are magnitudes, never -0). false, nothing written, when one is a NaN.
bool masked_maxima256(const float* s2, float* mx1, float* mx0) {
#if defined(FT4_METRICS_AVX2)
    // Lanes are the index's bits 0-2; vectors (32 of 8) bits 3-7.
    __m256 nan = _mm256_setzero_ps();
    __m256 all = _mm256_set1_ps(-INFINITY);
    __m256 v1[5], v0[5];
    for (int b = 0; b < 5; ++b) v1[b] = v0[b] = all;
    for (int g = 0; g < 32; ++g) {
        const __m256 v = _mm256_load_ps(s2 + 8 * g);
        nan = _mm256_or_ps(nan, _mm256_cmp_ps(v, v, _CMP_UNORD_Q));
        all = _mm256_max_ps(all, v);
        for (int b = 0; b < 5; ++b) {
            if ((g >> b) & 1)
                v1[b] = _mm256_max_ps(v1[b], v);
            else
                v0[b] = _mm256_max_ps(v0[b], v);
        }
    }
    if (_mm256_movemask_ps(nan) != 0) return false;
    alignas(32) float a[8], t[8];
    _mm256_store_ps(a, all);
    for (int b = 0; b < 3; ++b) {
        float m1 = -INFINITY, m0 = -INFINITY;
        for (int l = 0; l < 8; ++l) {
            if ((l >> b) & 1)
                m1 = std::max(m1, a[l]);
            else
                m0 = std::max(m0, a[l]);
        }
        mx1[b] = m1;
        mx0[b] = m0;
    }
    auto hmax = [&](__m256 v) {
        _mm256_store_ps(t, v);
        float m = t[0];
        for (int l = 1; l < 8; ++l) m = std::max(m, t[l]);
        return m;
    };
    for (int b = 0; b < 5; ++b) {
        mx1[3 + b] = hmax(v1[b]);
        mx0[3 + b] = hmax(v0[b]);
    }
    return true;
#elif defined(FT4_METRICS_NEON)
    // Lanes are the index's bits 0-1; vectors (64 of 4) bits 2-7.
    uint32x4_t ok = vdupq_n_u32(0xffffffffu);
    float32x4_t all = vdupq_n_f32(-INFINITY);
    float32x4_t v1[6], v0[6];
    for (int b = 0; b < 6; ++b) v1[b] = v0[b] = all;
    for (int g = 0; g < 64; ++g) {
        const float32x4_t v = vld1q_f32(s2 + 4 * g);
        ok = vandq_u32(ok, vceqq_f32(v, v));
        all = vmaxq_f32(all, v);
        for (int b = 0; b < 6; ++b) {
            if ((g >> b) & 1)
                v1[b] = vmaxq_f32(v1[b], v);
            else
                v0[b] = vmaxq_f32(v0[b], v);
        }
    }
    if (vminvq_u32(ok) == 0) return false;
    alignas(16) float a[4];
    vst1q_f32(a, all);
    for (int b = 0; b < 2; ++b) {
        float m1 = -INFINITY, m0 = -INFINITY;
        for (int l = 0; l < 4; ++l) {
            if ((l >> b) & 1)
                m1 = std::max(m1, a[l]);
            else
                m0 = std::max(m0, a[l]);
        }
        mx1[b] = m1;
        mx0[b] = m0;
    }
    for (int b = 0; b < 6; ++b) {
        mx1[2 + b] = vmaxvq_f32(v1[b]);
        mx0[2 + b] = vmaxvq_f32(v0[b]);
    }
    return true;
#else
    for (int i = 0; i < 256; ++i)
        if (std::isnan(s2[i])) return false;
    for (int b = 0; b < 8; ++b) {
        float m1 = -INFINITY, m0 = -INFINITY;
        for (int i = 0; i < 256; ++i) {
            if ((i >> b) & 1)
                m1 = std::max(m1, s2[i]);
            else
                m0 = std::max(m0, s2[i]);
        }
        mx1[b] = m1;
        mx0[b] = m0;
    }
    return true;
#endif
}

}  // namespace

void Frontend::get_ft4_bitmetrics(const cfloat* cd, float* bitmetrics, bool& badsync) {
    static const int icos4a[4] = {0, 1, 3, 2};
    static const int icos4b[4] = {1, 0, 2, 3};
    static const int icos4c[4] = {2, 3, 1, 0};
    static const int icos4d[4] = {3, 2, 0, 1};
    static const int graymap[4] = {0, 1, 3, 2};
    constexpr int N2 = 2 * NN;  // 206 rows a column

    // cs(0:3,NN) and s4(0:3,NN): (q,k) at [(k-1)*4 + q].
    cfloat cs[4 * NN];
    float s4[4 * NN];
    alignas(32) float s2[256];
    auto bmet = [&](int i, int j) -> float& { return bitmetrics[(j - 1) * N2 + i - 1]; };

    Metrics& m = *metrics_;
    cfloat* csymb = m.csymb.data();
    for (int k = 1; k <= NN; ++k) {
        const int i1 = (k - 1) * NSS;
        for (int q = 0; q < NSS; ++q) csymb[q] = cd[i1 + q];
        if (trace_) {
            const int al = static_cast<int>(reinterpret_cast<uintptr_t>(csymb) % 16);
            if (al != m.al) {
                m.al = al;
                trace_("@ bm-csymb " + I(al));
            }
        }
        fst4w::fft_c2c_inplace(csymb, NSS, -1);
        for (int q = 0; q < 4; ++q) cs[(k - 1) * 4 + q] = csymb[q];
    }
    {
        // s4 = abs(cs): cabsf's bits, four at a time.
        alignas(32) float xr[4 * NN], xi[4 * NN];
        for (int j = 0; j < 4 * NN; ++j) {
            xr[j] = cs[j].real();
            xi[j] = cs[j].imag();
        }
        ft8::hypotf_n_glibc(xr, xi, s4, 4 * NN);
    }

    // Sync quality check
    int is1 = 0, is2 = 0, is3 = 0, is4 = 0;
    badsync = false;
    auto ip = [&](int k) { return maxloc(1, 4, [&](int i) { return s4[(k - 1) * 4 + i - 1]; }); };
    for (int k = 1; k <= 4; ++k) {
        if (icos4a[k - 1] == ip(k) - 1) ++is1;
        if (icos4b[k - 1] == ip(k + 33) - 1) ++is2;
        if (icos4c[k - 1] == ip(k + 66) - 1) ++is3;
        if (icos4d[k - 1] == ip(k + 99) - 1) ++is4;
    }
    const int nsync = is1 + is2 + is3 + is4;  // Number of correct hard sync symbols, 0-16
    if (trace_) trace_("M " + X(h_c(cs, 4 * NN)) + " " + X(h_r(s4, 4 * NN)) + " " + I(nsync));
    if (nsync < 8) {
        badsync = true;
        return;
    }

    auto c = [&](int q, int ks) { return cs[(ks - 1) * 4 + q]; };
    for (int nseq = 1; nseq <= 3; ++nseq) {  // Try coherent sequences of 1, 2, and 4 symbols
        const int nsym = nseq == 1 ? 1 : nseq == 2 ? 2 : 4;
        const int nt = 1 << (2 * nsym);
        for (int ks = 1; ks <= NN - nsym + 1; ks += nsym) {  // 87+16=103 symbols.
            const int ipt = 1 + (ks - 1) * 2;
            const int ibmax = nsym == 1 ? 1 : nsym == 2 ? 3 : 7;
            if (nsym == 4) {
                // The 256 sums made with their shared left parts
                // once (((c1+c2)+c3)+c4, the same additions in the same
                // order), the magnitudes by hypotf_n_glibc (cabsf's bits),
                // and with no NaN among them the masked maxima in any
                // order (the values are >= +0 or +inf: a maximum is the
                // same float whichever order finds it).
                cfloat p12[16], p123[64];
                for (int a = 0; a < 4; ++a)
                    for (int b = 0; b < 4; ++b) p12[a * 4 + b] = cadd(c(graymap[a], ks), c(graymap[b], ks + 1));
                for (int ab = 0; ab < 16; ++ab)
                    for (int q = 0; q < 4; ++q) p123[ab * 4 + q] = cadd(p12[ab], c(graymap[q], ks + 2));
                alignas(32) float xr[256], xi[256];
                for (int abc = 0; abc < 64; ++abc)
                    for (int q = 0; q < 4; ++q) {
                        const cfloat v = cadd(p123[abc], c(graymap[q], ks + 3));
                        xr[abc * 4 + q] = v.real();
                        xi[abc * 4 + q] = v.imag();
                    }
                ft8::hypotf_n_glibc(xr, xi, s2, 256);
                float mx1[8], mx0[8];
                if (masked_maxima256(s2, mx1, mx0)) {
                    for (int ib = 0; ib <= ibmax; ++ib) {
                        const float bm = mx1[ibmax - ib] - mx0[ibmax - ib];
                        if (ipt + ib > N2) continue;
                        bmet(ipt + ib, nseq) = bm;
                    }
                    continue;
                }
            } else if (nsym == 2) {
                // The 16 magnitudes four at a time (cabsf's bits).
                alignas(32) float xr[16], xi[16];
                for (int i = 0; i < 16; ++i) {
                    const cfloat v = cadd(c(graymap[(i & 15) / 4], ks), c(graymap[i & 3], ks + 1));
                    xr[i] = v.real();
                    xi[i] = v.imag();
                }
                ft8::hypotf_n_glibc(xr, xi, s2, 16);
            } else {
                for (int i = 0; i < nt; ++i) {
                    const int i1 = i / 64;
                    const int i2 = (i & 63) / 16;
                    const int i3 = (i & 15) / 4;
                    const int i4 = i & 3;
                    if (nsym == 1) {
                        s2[i] = s4[(ks - 1) * 4 + graymap[i4]];  // cabs_(c(graymap(i4), ks)), made above
                    } else if (nsym == 2) {
                        s2[i] = cabs_(cadd(c(graymap[i3], ks), c(graymap[i4], ks + 1)));
                    } else {
                        s2[i] = cabs_(cadd(cadd(cadd(c(graymap[i1], ks), c(graymap[i2], ks + 1)),
                                                c(graymap[i3], ks + 2)),
                                           c(graymap[i4], ks + 3)));
                    }
                }
            }
            for (int ib = 0; ib <= ibmax; ++ib) {
                const float bm = maxval_bit(s2, nt, ibmax - ib, true) - maxval_bit(s2, nt, ibmax - ib, false);
                if (ipt + ib > N2) continue;
                bmet(ipt + ib, nseq) = bm;
            }
        }
    }

    bmet(205, 2) = bmet(205, 1);
    bmet(206, 2) = bmet(206, 1);
    for (int i = 201; i <= 204; ++i) bmet(i, 3) = bmet(i, 2);
    bmet(205, 3) = bmet(205, 1);
    bmet(206, 3) = bmet(206, 1);

    if (trace_)
        trace_("m " + X(h_r(bitmetrics, N2)) + " " + X(h_r(bitmetrics + N2, N2)) + " " +
               X(h_r(bitmetrics + 2 * N2, N2)));
    ft8::normalizebmet(bitmetrics, N2);
    ft8::normalizebmet(bitmetrics + N2, N2);
    ft8::normalizebmet(bitmetrics + 2 * N2, N2);
    if (trace_)
        trace_("B " + X(h_r(bitmetrics, N2)) + " " + X(h_r(bitmetrics + N2, N2)) + " " +
               X(h_r(bitmetrics + 2 * N2, N2)));
}

// get_ft4_bitmetrics in its plain form, line for line with the
// reference (every magnitude by cabs_, maxval_bit for every bit), kept to
// check get_ft4_bitmetrics against (tests/test_ft4_frontend.cpp,
// [ft4-fe-reference]).
void Frontend::get_ft4_bitmetrics_reference(const cfloat* cd, float* bitmetrics, bool& badsync) {
    static const int icos4a[4] = {0, 1, 3, 2};
    static const int icos4b[4] = {1, 0, 2, 3};
    static const int icos4c[4] = {2, 3, 1, 0};
    static const int icos4d[4] = {3, 2, 0, 1};
    static const int graymap[4] = {0, 1, 3, 2};
    constexpr int N2 = 2 * NN;  // 206 rows a column

    // cs(0:3,NN) and s4(0:3,NN): (q,k) at [(k-1)*4 + q].
    cfloat cs[4 * NN];
    float s4[4 * NN];
    float s2[256];
    auto bmet = [&](int i, int j) -> float& { return bitmetrics[(j - 1) * N2 + i - 1]; };

    Metrics& m = *metrics_;
    cfloat* csymb = m.csymb.data();
    for (int k = 1; k <= NN; ++k) {
        const int i1 = (k - 1) * NSS;
        for (int q = 0; q < NSS; ++q) csymb[q] = cd[i1 + q];
        if (trace_) {
            const int al = static_cast<int>(reinterpret_cast<uintptr_t>(csymb) % 16);
            if (al != m.al) {
                m.al = al;
                trace_("@ bm-csymb " + I(al));
            }
        }
        fst4w::fft_c2c_inplace(csymb, NSS, -1);
        for (int q = 0; q < 4; ++q) {
            cs[(k - 1) * 4 + q] = csymb[q];
            s4[(k - 1) * 4 + q] = cabs_(csymb[q]);
        }
    }

    // Sync quality check
    int is1 = 0, is2 = 0, is3 = 0, is4 = 0;
    badsync = false;
    auto ip = [&](int k) { return maxloc(1, 4, [&](int i) { return s4[(k - 1) * 4 + i - 1]; }); };
    for (int k = 1; k <= 4; ++k) {
        if (icos4a[k - 1] == ip(k) - 1) ++is1;
        if (icos4b[k - 1] == ip(k + 33) - 1) ++is2;
        if (icos4c[k - 1] == ip(k + 66) - 1) ++is3;
        if (icos4d[k - 1] == ip(k + 99) - 1) ++is4;
    }
    const int nsync = is1 + is2 + is3 + is4;  // Number of correct hard sync symbols, 0-16
    if (trace_) trace_("M " + X(h_c(cs, 4 * NN)) + " " + X(h_r(s4, 4 * NN)) + " " + I(nsync));
    if (nsync < 8) {
        badsync = true;
        return;
    }

    auto c = [&](int q, int ks) { return cs[(ks - 1) * 4 + q]; };
    for (int nseq = 1; nseq <= 3; ++nseq) {  // Try coherent sequences of 1, 2, and 4 symbols
        const int nsym = nseq == 1 ? 1 : nseq == 2 ? 2 : 4;
        const int nt = 1 << (2 * nsym);
        for (int ks = 1; ks <= NN - nsym + 1; ks += nsym) {  // 87+16=103 symbols.
            for (int i = 0; i < nt; ++i) {
                const int i1 = i / 64;
                const int i2 = (i & 63) / 16;
                const int i3 = (i & 15) / 4;
                const int i4 = i & 3;
                if (nsym == 1) {
                    s2[i] = cabs_(c(graymap[i4], ks));
                } else if (nsym == 2) {
                    s2[i] = cabs_(cadd(c(graymap[i3], ks), c(graymap[i4], ks + 1)));
                } else {
                    s2[i] = cabs_(cadd(cadd(cadd(c(graymap[i1], ks), c(graymap[i2], ks + 1)), c(graymap[i3], ks + 2)),
                                       c(graymap[i4], ks + 3)));
                }
            }
            const int ipt = 1 + (ks - 1) * 2;
            const int ibmax = nsym == 1 ? 1 : nsym == 2 ? 3 : 7;
            for (int ib = 0; ib <= ibmax; ++ib) {
                const float bm = maxval_bit(s2, nt, ibmax - ib, true) - maxval_bit(s2, nt, ibmax - ib, false);
                if (ipt + ib > N2) continue;
                bmet(ipt + ib, nseq) = bm;
            }
        }
    }

    bmet(205, 2) = bmet(205, 1);
    bmet(206, 2) = bmet(206, 1);
    for (int i = 201; i <= 204; ++i) bmet(i, 3) = bmet(i, 2);
    bmet(205, 3) = bmet(205, 1);
    bmet(206, 3) = bmet(206, 1);

    if (trace_)
        trace_("m " + X(h_r(bitmetrics, N2)) + " " + X(h_r(bitmetrics + N2, N2)) + " " +
               X(h_r(bitmetrics + 2 * N2, N2)));
    ft8::normalizebmet(bitmetrics, N2);
    ft8::normalizebmet(bitmetrics + N2, N2);
    ft8::normalizebmet(bitmetrics + 2 * N2, N2);
    if (trace_)
        trace_("B " + X(h_r(bitmetrics, N2)) + " " + X(h_r(bitmetrics + N2, N2)) + " " +
               X(h_r(bitmetrics + 2 * N2, N2)));
}

// ft4_decode.f90:297-303.
int ft4_nsync_qual(const float* bitmetrics) {
    static const int8_t p1[8] = {0, 0, 0, 1, 1, 0, 1, 1}, p2[8] = {0, 1, 0, 0, 1, 1, 1, 0},
                        p3[8] = {1, 1, 1, 0, 0, 1, 0, 0}, p4[8] = {1, 0, 1, 1, 0, 0, 0, 1};
    // hbits=0; where(bitmetrics(:,1).ge.0) hbits=1 (a NaN gives 0)
    auto hb = [&](int i) { return bitmetrics[i - 1] >= 0.0f ? 1 : 0; };
    int ns1 = 0, ns2 = 0, ns3 = 0, ns4 = 0;
    for (int k = 0; k < 8; ++k) {
        ns1 += hb(1 + k) == p1[k];
        ns2 += hb(67 + k) == p2[k];
        ns3 += hb(133 + k) == p3[k];
        ns4 += hb(199 + k) == p4[k];
    }
    return ns1 + ns2 + ns3 + ns4;
}

// ft4_decode.f90:306-318: llrX(1:58)=bitmetrics(9:66,X),
// (59:116)=(75:132,X), (117:174)=(141:198,X) (the data bits, the four
// Costas arrays' 8 rows each skipped), then llrX=scalefac*llrX with the
// real*4 scalefac 2.83.
void ft4_llrs(const float* bitmetrics, float* llra, float* llrb, float* llrc) {
    const float scalefac = 2.83f;
    float* llr[3] = {llra, llrb, llrc};
    for (int j = 0; j < 3; ++j) {
        const float* b = bitmetrics + j * 2 * NN;
        float* l = llr[j];
        for (int i = 0; i < 58; ++i) {
            l[i] = scalefac * b[8 + i];
            l[58 + i] = scalefac * b[74 + i];
            l[116 + i] = scalefac * b[140 + i];
        }
    }
}

}  // namespace ft4
