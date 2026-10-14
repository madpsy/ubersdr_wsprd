// FT4's front end (ft4_frontend.h): lib/ft4/ft4_downsample.f90 and
// sync4d.f90 at 567ad29ce, line for line, as their Release objects
// compute them (from their disassembly, site by site):
//  - ft4_downsample's window: two raised-cosine ramps of iwt = 63 points,
//    each sent through libmvec's cosf four at a time on its first 60
//    points (cosf4_glibc), its last three folded by GCC (their bits
//    below); then cshift by iws = 126. c1 = c1*window/NFFT2 as gfortran's
//    complex-by-real product and Smith's complex-by-real quotient (the
//    zero imaginary part's terms kept: they decide zeros' signs). The big
//    r2c FFT (SAVEd x, NMAX = 72,576) only when newdata, which is not
//    cleared here. c1(NFFT2/2) written twice, in the loop's order.
//  - sync4d: the four Costas tables by sincosf and fmod (fprem: exact), the
//    phases running on across the four tones; stride-2 sums of 64 complex
//    products in index order, the partial head (i1 < 0) and tail (past
//    NP); p(z) = powf(re^2 + im^2, 0.5) of z*fac by scalar powf (not
//    sqrtf), the four added left to right.
//
// Joe Taylor, K1JT, Steven Franke, K9AN, and the WSJT-X developers (GNU GPL
// v3), and so is this.
#include "ft4_frontend.h"

#include <algorithm>
#include <cstdint>
#include <cstring>

#if defined(__AVX2__) && (defined(__x86_64__) || defined(_M_X64))
#include <immintrin.h>
#define FT4_SYNC_AVX2 1
#elif defined(__aarch64__) && defined(__ARM_NEON)
#include <arm_neon.h>
#define FT4_SYNC_NEON 1
#endif

#include "fst4w_frontend.h"  // FftwVec, the cached four2a plans
#include "ft8_fortran.h"     // the Fortran helpers, FT8's
#include "ft8_maths.h"       // glibc 2.39's maths and x86-64 libmvec, carried

namespace ft4 {

using namespace ft8::fortran;
using ft8::cosf4_glibc;
using ft8::powf_glibc;

namespace {

constexpr int NFFT2 = NMAX / NDOWN;  // 4,032 (= NP, NDMAX)
constexpr int NP = NMAX / NDOWN;
// sync4d_many's zero padding of its copies of cd0, as elements of cd0
// before 0 and past NP-1 (a stride-4 copy has a quarter as many floats).
constexpr int PADX = 512;

// ft4_downsample's, folded (df = 12000.0/NMAX; iwt, iwf, iws from
// 0.5*baud/df, 4*baud/df, baud/df, baud = 12000.0/NSPS).
const float DF = bitsf(0x3E294FEAu);
constexpr int IWT = 63, IWF = 504, IWS = 126;

// sync4d's Costas arrays.
const int icos4a[4] = {0, 1, 3, 2};
const int icos4b[4] = {1, 0, 2, 3};
const int icos4c[4] = {2, 3, 1, 0};
const int icos4d[4] = {3, 2, 0, 1};

// a*x, x real: gfortran's complex-by-real product, (x, 0) kept as such.
inline cfloat cmulr(cfloat a, float x) {
    return cfloat(a.real() * x - a.imag() * 0.0f, a.real() * 0.0f + a.imag() * x);
}

// a/d, d a real constant (not a power of two): Smith's rule with ratio
// 0/d = +0 and denominator d + 0*0 = d, folded.
inline cfloat cdivr(cfloat a, float d) {
    return cfloat((a.imag() * 0.0f + a.real()) / d, (a.imag() - a.real() * 0.0f) / d);
}

// sum(cd0(i:i+2*(n-1):2)*conjg(cs(0:n-1))), in index order from zero.
inline cfloat zsum(const cfloat* cd0, int i, const cfloat* cs, int n) {
    float re = 0.0f, im = 0.0f;
    for (int m = 0; m < n; ++m) {
        const cfloat p = cmul(cd0[i + 2 * m], cfloat(cs[m].real(), -cs[m].imag()));
        re = re + p.real();
        im = im + p.imag();
    }
    return cfloat(re, im);
}

int align16(const void* p) { return static_cast<int>(reinterpret_cast<uintptr_t>(p) % 16u); }

}  // namespace

// This group's SAVEd state.
struct Frontend::Sync {
    // ft4_downsample's
    bool ds_first = true;
    float window[NFFT2] = {};
    fst4w::FftwVec<float> x = fst4w::FftwVec<float>(NMAX + 2, 0.0f);  // x(NMAX), cx(0:NMAX/2) equivalenced
    fst4w::FftwVec<cfloat> c1 = fst4w::FftwVec<cfloat>(NFFT2);         // a local there, 16-byte aligned
    // sync4d's
    bool s4d_first = true;
    float twopi = 0.0f, fac = 0.0f;
    cfloat csynca[2 * NSS] = {}, csyncb[2 * NSS] = {}, csyncc[2 * NSS] = {}, csyncd[2 * NSS] = {};
    cfloat csync2[2 * NSS] = {};  // a local there; with itwk != 1 what it last held
    // The trace's '@' lines (loc() mod 16 of ds-x, ds-c1; -1 none yet).
    int al_x = -1, al_c1 = -1;
    // sync4d_many's copies of cd0 (sync4d_prepare): real and imaginary
    // parts apart (r1, i1: element k at [PADX + k]) and, for a stride of
    // 4, the four phases apart (dr[p], di[p]: element 4j+p at
    // [PADX/4 + j]); PADX elements of zeros before and after, so that the
    // lanes of a batch whose sums sync4d makes otherwise (partial, or
    // past the ends) still read inside, their results unused.
    const cfloat* prepared = nullptr;
    alignas(32) float r1[NP + 2 * PADX] = {}, i1[NP + 2 * PADX] = {};
    alignas(32) float dr[4][(NP + 2 * PADX) / 4] = {}, di[4][(NP + 2 * PADX) / 4] = {};
};

Frontend::Sync* Frontend::make_sync() { return new Sync(); }
void Frontend::del_sync(Sync* p) { delete p; }
size_t Frontend::sync_bytes(const Sync* p) {
    return sizeof(Sync) + p->x.size() * sizeof(float) + p->c1.size() * sizeof(cfloat);
}

void Frontend::ft4_downsample(const float* dd, bool newdata, float f0, cfloat* c) {
    Sync& s = *sync_;
    float* window = s.window;
    if (s.ds_first) {
        // window(0:iwt-1) = 0.5*(1+cos(pi*(/(i,i=iwt-1,0,-1)/)/iwt)):
        // elements 0..59 by libmvec, 60..62 folded.
        for (int k = 0; k < 60; k += 4) {
            float a[4], cs[4];
            for (int q = 0; q < 4; ++q) a[q] = (PI * static_cast<float>(IWT - 1 - (k + q))) / static_cast<float>(IWT);
            cosf4_glibc(a, cs);
            for (int q = 0; q < 4; ++q) window[k + q] = (1.0f + cs[q]) * 0.5f;
        }
        window[60] = bitsf(0x3F7F5D2Bu);
        window[61] = bitsf(0x3F7FD744u);
        window[62] = bitsf(0x3F800000u);
        for (int k = IWT; k <= IWT + IWF - 1; ++k) window[k] = 1.0f;
        // window(iwt+iwf:2*iwt+iwf-1) = 0.5*(1+cos(pi*(/(i,i=0,iwt-1)/)/iwt))
        for (int k = 0; k < 60; k += 4) {
            float a[4], cs[4];
            for (int q = 0; q < 4; ++q) a[q] = (PI * static_cast<float>(k + q)) / static_cast<float>(IWT);
            cosf4_glibc(a, cs);
            for (int q = 0; q < 4; ++q) window[IWT + IWF + k + q] = (1.0f + cs[q]) * 0.5f;
        }
        window[IWT + IWF + 60] = bitsf(0x3BB6FF00u);
        window[IWT + IWF + 61] = bitsf(0x3B22D500u);
        window[IWT + IWF + 62] = bitsf(0x3A22EE00u);
        for (int k = 2 * IWT + IWF; k < NFFT2; ++k) window[k] = 0.0f;
        // window = cshift(window, iws)
        float tmp[NFFT2];
        for (int k = 0; k < NFFT2; ++k) tmp[k] = window[(k + IWS) % NFFT2];
        std::memcpy(window, tmp, sizeof tmp);
        s.ds_first = false;
        if (trace_) trace_("D " + X(h_r(window, NFFT2)) + " " + I(IWT) + " " + I(IWF) + " " + I(IWS));
    }

    float* x = s.x.data();
    if (newdata) {
        std::memcpy(x, dd, sizeof(float) * NMAX);
        if (trace_ && align16(x) != s.al_x) {
            s.al_x = align16(x);
            trace_("@ ds-x " + I(s.al_x));
        }
        fst4w::fft_r2c_inplace(x, NMAX);  // r2c FFT to freq domain
        if (trace_) trace_("X " + X(h_r(x, NMAX + 2)));
    }
    const cfloat* cx = reinterpret_cast<const cfloat*>(x);
    const int i0 = nint(f0 / DF);
    cfloat* c1 = s.c1.data();
    std::fill(c1, c1 + NFFT2, cfloat(0.0f, 0.0f));
    if (i0 >= 0 && i0 <= NMAX / 2) c1[0] = cx[i0];
    for (int i = 1; i <= NFFT2 / 2; ++i) {
        if (i0 + i >= 0 && i0 + i <= NMAX / 2) c1[i] = cx[i0 + i];
        if (i0 - i >= 0 && i0 - i <= NMAX / 2) c1[NFFT2 - i] = cx[i0 - i];
    }
    // c1 = c1*window/NFFT2
    for (int k = 0; k < NFFT2; ++k) c1[k] = cdivr(cmulr(c1[k], window[k]), static_cast<float>(NFFT2));
    uint32_t h1 = 0;
    if (trace_) {
        h1 = h_c(c1, NFFT2);
        if (align16(c1) != s.al_c1) {
            s.al_c1 = align16(c1);
            trace_("@ ds-c1 " + I(s.al_c1));
        }
    }
    fst4w::fft_c2c_inplace(c1, NFFT2, 1);  // c2c FFT back to time domain
    std::memcpy(c, c1, sizeof(cfloat) * NFFT2);
    if (trace_) trace_("x " + XF(f0) + " " + I(i0) + " " + X(h1) + " " + X(h_c(c, NFFT2)));
}

// sync4d's first call: its four Costas tables, fac and twopi (SAVEd).
void Frontend::sync4d_init() {
    Sync& s = *sync_;
    if (s.s4d_first) {
        s.twopi = TWOPI;
        const float twopi = s.twopi;
        int k = 0;
        float phia = 0.0f, phib = 0.0f, phic = 0.0f, phid = 0.0f;
        for (int i = 0; i <= 3; ++i) {
            const float dphia = ((2.0f * twopi) * static_cast<float>(icos4a[i])) / static_cast<float>(NSS);
            const float dphib = ((2.0f * twopi) * static_cast<float>(icos4b[i])) / static_cast<float>(NSS);
            const float dphic = ((2.0f * twopi) * static_cast<float>(icos4c[i])) / static_cast<float>(NSS);
            const float dphid = ((2.0f * twopi) * static_cast<float>(icos4d[i])) / static_cast<float>(NSS);
            for (int j = 1; j <= NSS / 2; ++j) {
                s.csynca[k] = cis(phia);
                s.csyncb[k] = cis(phib);
                s.csyncc[k] = cis(phic);
                s.csyncd[k] = cis(phid);
                phia = std::fmod(phia + dphia, twopi);
                phib = std::fmod(phib + dphib, twopi);
                phic = std::fmod(phic + dphic, twopi);
                phid = std::fmod(phid + dphid, twopi);
                ++k;
            }
        }
        s.s4d_first = false;
        s.fac = 1.0f / (2.0f * static_cast<float>(NSS));
        if (trace_)
            trace_("T " + X(h_c(s.csynca, 2 * NSS)) + " " + X(h_c(s.csyncb, 2 * NSS)) + " " +
                   X(h_c(s.csyncc, 2 * NSS)) + " " + X(h_c(s.csyncd, 2 * NSS)) + " " + XF(s.fac) + " " +
                   XF(s.twopi));
    }
}

float Frontend::sync4d(const cfloat* cd0, int i0, const cfloat* ctwk, int itwk) {
    Sync& s = *sync_;
    sync4d_init();
    const float fac = s.fac;
    cfloat* csync2 = s.csync2;
    auto tweak = [&](const cfloat* cs) {
        if (itwk == 1)
            for (int j = 0; j < 2 * NSS; ++j) csync2[j] = cmul(ctwk[j], cs[j]);  // Tweak the frequency
    };

    const int i1 = i0;  // four Costas arrays
    const int i2 = i0 + 33 * NSS;
    const int i3 = i0 + 66 * NSS;
    const int i4 = i0 + 99 * NSS;
    const cfloat zero(0.0f, 0.0f);
    cfloat z1 = zero, z2 = zero, z3 = zero, z4 = zero;

    tweak(s.csynca);
    if (i1 >= 0 && i1 + 4 * NSS - 1 <= NP - 1) {
        z1 = zsum(cd0, i1, csync2, 2 * NSS);
    } else if (i1 < 0) {
        const int npts = (i1 + 4 * NSS - 1) / 2;
        if (npts <= 16)
            z1 = zero;
        else  // cd0(0:i1+4*NSS-1:2) with csync2(2*NSS-npts:), npts+1 terms
            z1 = zsum(cd0, 0, csync2 + (2 * NSS - npts - 1), npts + 1);
    }

    tweak(s.csyncb);
    if (i2 >= 0 && i2 + 4 * NSS - 1 <= NP - 1) z2 = zsum(cd0, i2, csync2, 2 * NSS);

    tweak(s.csyncc);
    if (i3 >= 0 && i3 + 4 * NSS - 1 <= NP - 1) z3 = zsum(cd0, i3, csync2, 2 * NSS);

    tweak(s.csyncd);
    if (i4 >= 0 && i4 + 4 * NSS - 1 <= NP - 1) {
        z4 = zsum(cd0, i4, csync2, 2 * NSS);
    } else if (i4 + 4 * NSS - 1 > NP - 1) {
        const int npts = (NP - 1 - i4 + 1) / 2;
        if (npts <= 16)
            z4 = zero;
        else  // cd0(i4:i4+2*npts-1:2) with csync2(1:npts)
            z4 = zsum(cd0, i4, csync2, npts);
    }

    // p(z) = (real(z*fac)**2 + aimag(z*fac)**2)**0.5, a powf
    auto p = [&](cfloat z) {
        const cfloat w = cmulr(z, fac);
        return powf_glibc(w.real() * w.real() + w.imag() * w.imag(), 0.5f);
    };
    return ((p(z1) + p(z2)) + p(z3)) + p(z4);
}

void Frontend::sync4d_prepare(const cfloat* cd0) {
    Sync& s = *sync_;
    for (int k = 0; k < NP; ++k) {
        s.r1[PADX + k] = cd0[k].real();
        s.i1[PADX + k] = cd0[k].imag();
        s.dr[k & 3][PADX / 4 + (k >> 2)] = cd0[k].real();
        s.di[k & 3][PADX / 4 + (k >> 2)] = cd0[k].imag();
    }
    s.prepared = cd0;
}

// sync4d(cd0, istart0 + k*step, ctwk, 1) for k = 0..n-1, into sync[k]:
// the same operations on the same numbers, sixteen starts at a time. Each
// lane is one call's arithmetic: its four sums of 64 products in index
// order from zero (separate accumulators per lane, no reordering, no fused
// multiply-add), the tweaked tables made by the same products
// (cmul(ctwk(j), csync(j))) and the conjugate taken as there. A lane whose
// sum for an array is partial or zero (i1 < 0, i4 past the end) gets it
// from sync4d's own branches; p(z)'s powf by powf_half_glibc (the same
// bits) and the four added left to right. csync2 is left as sync4d's
// last call leaves it. Needs sync4d_prepare(cd0) first, after cd0's last
// change.
void Frontend::sync4d_many(const cfloat* cd0, int istart0, int n, int step, const cfloat* ctwk, float* sync) {
    Sync& s = *sync_;
    if (n <= 0) return;
    sync4d_init();
    const float fac = s.fac;
    const cfloat* tabs[4] = {s.csynca, s.csyncb, s.csyncc, s.csyncd};
    cfloat tw[4][2 * NSS];  // csync2 for each array, as the tweak makes it
    alignas(32) float tre[4][2 * NSS], tim[4][2 * NSS];  // conjg(csync2): real, imaginary
    for (int t = 0; t < 4; ++t)
        for (int j = 0; j < 2 * NSS; ++j) {
            tw[t][j] = cmul(ctwk[j], tabs[t][j]);  // Tweak the frequency
            tre[t][j] = tw[t][j].real();
            tim[t][j] = -tw[t][j].imag();
        }
    std::memcpy(s.csync2, tw[3], sizeof s.csync2);
    const bool vec = s.prepared == cd0 && (step == 1 || step == 4);
    const int off[4] = {0, 33 * NSS, 66 * NSS, 99 * NSS};
    const cfloat zero(0.0f, 0.0f);
    // A table's sum at start i as sync4d makes it: whole, the head
    // (i1 < 0, first table) or tail (i4 past NP, last table), or zero.
    auto zone = [&](int t, int i) -> cfloat {
        const int it = i + off[t];
        if (it >= 0 && it + 4 * NSS - 1 <= NP - 1) return zsum(cd0, it, tw[t], 2 * NSS);
        if (t == 0 && it < 0) {
            const int npts = (it + 4 * NSS - 1) / 2;
            if (npts <= 16) return zero;
            return zsum(cd0, 0, tw[t] + (2 * NSS - npts - 1), npts + 1);
        }
        if (t == 3 && it + 4 * NSS - 1 > NP - 1) {
            const int npts = (NP - 1 - it + 1) / 2;
            if (npts <= 16) return zero;
            return zsum(cd0, it, tw[t], npts);
        }
        return zero;
    };
    constexpr int L = 16;
    for (int b0 = 0; b0 < n; b0 += L) {
        const int nl = std::min(L, n - b0);
        const int ib = istart0 + b0 * step;  // lane 0's start
        alignas(32) float zr[4][L], zi[4][L];
        bool whole[4][L];  // the lane's sum is the whole 64 terms (sync4d's first branch)
        for (int t = 0; t < 4; ++t) {
            bool any = false;
            for (int l = 0; l < L; ++l) {
                const int it = ib + l * step + off[t];
                whole[t][l] = vec && l < nl && it >= 0 && it + 4 * NSS - 1 <= NP - 1;
                any = any || whole[t][l];
            }
            // Every lane's reads (element x + l*step + 2m, m <= 63) inside
            // the padded copies.
            const int x = ib + off[t];
            const bool inside = x >= -PADX && x + (L - 1) * step + 4 * NSS - 2 < NP + PADX;
            if (!any || !inside) {
                for (int l = 0; l < L; ++l) whole[t][l] = false;
                continue;
            }
            // Lane l's element at m: cd0(x + l*step + 2m); even m from one
            // copy, odd m from the other (stride 4: the phase alternates),
            // each accumulator's terms in m's order.
            const float *re0, *im0, *re1, *im1;
            int inc;
            if (step == 1) {
                re0 = s.r1 + PADX + x, im0 = s.i1 + PADX + x;
                re1 = re0 + 2, im1 = im0 + 2;
                inc = 4;
            } else {
                re0 = s.dr[x & 3] + PADX / 4 + (x >> 2), im0 = s.di[x & 3] + PADX / 4 + (x >> 2);
                re1 = s.dr[(x + 2) & 3] + PADX / 4 + ((x + 2) >> 2);
                im1 = s.di[(x + 2) & 3] + PADX / 4 + ((x + 2) >> 2);
                inc = 1;
            }
            const float* tr = tre[t];
            const float* ti = tim[t];
#if defined(FT4_SYNC_AVX2)
            constexpr int h = 8;  // lanes 8-15: eight floats on in either copy
            __m256 ar0 = _mm256_setzero_ps(), ai0 = _mm256_setzero_ps();
            __m256 ar1 = _mm256_setzero_ps(), ai1 = _mm256_setzero_ps();
            // cmul(a, b): (ar*br - ai*bi, ar*bi + ai*br), b = conjg(csync2(m))
            auto term = [&](const float* r, const float* im, int m) {
                const __m256 br = _mm256_set1_ps(tr[m]), bi = _mm256_set1_ps(ti[m]);
                const __m256 xr0 = _mm256_loadu_ps(r), xi0 = _mm256_loadu_ps(im);
                const __m256 xr1 = _mm256_loadu_ps(r + h), xi1 = _mm256_loadu_ps(im + h);
                ar0 = _mm256_add_ps(ar0, _mm256_sub_ps(_mm256_mul_ps(xr0, br), _mm256_mul_ps(xi0, bi)));
                ai0 = _mm256_add_ps(ai0, _mm256_add_ps(_mm256_mul_ps(xr0, bi), _mm256_mul_ps(xi0, br)));
                ar1 = _mm256_add_ps(ar1, _mm256_sub_ps(_mm256_mul_ps(xr1, br), _mm256_mul_ps(xi1, bi)));
                ai1 = _mm256_add_ps(ai1, _mm256_add_ps(_mm256_mul_ps(xr1, bi), _mm256_mul_ps(xi1, br)));
            };
            for (int m2 = 0; m2 < NSS; ++m2) {
                const int o = m2 * inc;
                term(re0 + o, im0 + o, 2 * m2);
                term(re1 + o, im1 + o, 2 * m2 + 1);
            }
            _mm256_store_ps(zr[t], ar0);
            _mm256_store_ps(zi[t], ai0);
            _mm256_store_ps(zr[t] + 8, ar1);
            _mm256_store_ps(zi[t] + 8, ai1);
#elif defined(FT4_SYNC_NEON)
            // As the AVX2 kernel, four registers of four lanes (lanes 4q..4q+3
            // at r + 4q), each lane's terms in m's order, no fused
            // multiply-add (-ffp-contract=off).
            float32x4_t ar[4], ai[4];
            for (int q = 0; q < 4; ++q) ar[q] = ai[q] = vdupq_n_f32(0.0f);
            auto term = [&](const float* r, const float* im, int m) {
                const float32x4_t br = vdupq_n_f32(tr[m]), bi = vdupq_n_f32(ti[m]);
                for (int q = 0; q < 4; ++q) {
                    const float32x4_t xr = vld1q_f32(r + 4 * q), xi = vld1q_f32(im + 4 * q);
                    ar[q] = vaddq_f32(ar[q], vsubq_f32(vmulq_f32(xr, br), vmulq_f32(xi, bi)));
                    ai[q] = vaddq_f32(ai[q], vaddq_f32(vmulq_f32(xr, bi), vmulq_f32(xi, br)));
                }
            };
            for (int m2 = 0; m2 < NSS; ++m2) {
                const int o = m2 * inc;
                term(re0 + o, im0 + o, 2 * m2);
                term(re1 + o, im1 + o, 2 * m2 + 1);
            }
            for (int q = 0; q < 4; ++q) {
                vst1q_f32(zr[t] + 4 * q, ar[q]);
                vst1q_f32(zi[t] + 4 * q, ai[q]);
            }
#else
            for (int l = 0; l < L; ++l) zr[t][l] = zi[t][l] = 0.0f;
            for (int m2 = 0; m2 < NSS; ++m2)
                for (int odd = 0; odd < 2; ++odd) {
                    const int m = 2 * m2 + odd;
                    const float* r = (odd ? re1 : re0) + m2 * inc;
                    const float* im = (odd ? im1 : im0) + m2 * inc;
                    for (int l = 0; l < L; ++l) {
                        const float pr = r[l] * tr[m] - im[l] * ti[m];
                        const float pi = r[l] * ti[m] + im[l] * tr[m];
                        zr[t][l] = zr[t][l] + pr;
                        zi[t][l] = zi[t][l] + pi;
                    }
                }
#endif
        }
        // p(z) = (real(z*fac)**2 + aimag(z*fac)**2)**0.5, powf's bits
        // four at a time (powf_half_glibc), the four added left to right.
        float a[L * 4], pz[L * 4];
        for (int l = 0; l < nl; ++l) {
            const int i = ib + l * step;
            for (int t = 0; t < 4; ++t) {
                const cfloat z = whole[t][l] ? cfloat(zr[t][l], zi[t][l]) : zone(t, i);
                const cfloat w = cmulr(z, fac);
                a[4 * l + t] = w.real() * w.real() + w.imag() * w.imag();
            }
        }
        ft8::powf_half_glibc(a, pz, 4 * nl);
        for (int l = 0; l < nl; ++l)
            sync[b0 + l] = ((pz[4 * l] + pz[4 * l + 1]) + pz[4 * l + 2]) + pz[4 * l + 3];
    }
}

}  // namespace ft4
