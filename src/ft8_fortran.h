// Fortran's small helpers as FT8's decoder has them, shared with FT4's:
// moved here verbatim from ft8_frontend.cpp
// (bitsf, fbits, PI, TWOPI, nint, maxloc, maxval_bit, maxval, the trace
// hashes, shell, pctile, determ, polyfit, cis, cabs_, cmul, cadd, cscale)
// and ft8.cpp (fix, trim, findex, db, twkfreq1, isgrid4), unchanged.
//
// They do float maths that must round as the reference's (gfortran's
// code, glibc 2.39's maths through ft8_maths.h), so they keep the internal
// linkage they had in their files' anonymous namespaces: every includer
// gets its own copy, compiled with its own flags, and no COMDAT copy is
// merged at link time with one compiled otherwise. Every includer must be
// compiled as ft8_frontend.cpp is (FT8_FP_OPTIONS: -ffp-contract=off, and
// -mno-fma with GCC < 13; CMakeLists.txt). [[maybe_unused]] only silences
// -Wunused-function in includers that do not call a helper.
//
// pctile.f90, shell.f90, polyfit.f90, determ.f90, twkfreq1.f90, db.f90 and
// ft8_a7.f90's isgrid4: Joe Taylor, K1JT, Steven Franke, K9AN, and the
// WSJT-X developers (GNU GPL v3), and so is this.
#pragma once

#include <cmath>
#include <complex>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <limits>
#include <string>
#include <vector>

#include "ft8_maths.h"  // glibc 2.39's maths, carried

namespace ft8 {
namespace fortran {
namespace {

using cfloat_ = std::complex<float>;

[[maybe_unused]] float bitsf(uint32_t u) {
    float f;
    std::memcpy(&f, &u, 4);
    return f;
}
[[maybe_unused]] uint32_t fbits(float f) {
    uint32_t u;
    std::memcpy(&u, &f, 4);
    return u;
}

// Constants as gfortran folds them.
[[maybe_unused]] const float PI = bitsf(0x40490FDBu);     // 4.0*atan(1.0)
[[maybe_unused]] const float TWOPI = bitsf(0x40C90FDBu);  // 8.0*atan(1.0)

// Fortran's nint of a default real.
inline int nint(float x) { return static_cast<int>(std::lround(x)); }

// gfortran's maxloc over at(lo..hi): the first maximum; NaNs never win
// unless all are NaN (then the first).
template <typename F>
int maxloc(int lo, int hi, F at) {
    bool found = false;
    int best = lo;
    float m = 0.0f;
    for (int i = lo; i <= hi; ++i) {
        const float v = at(i);
        if (!found) {
            if (!std::isnan(v)) {
                found = true;
                best = i;
                m = v;
            }
        } else if (v > m) {
            best = i;
            m = v;
        }
    }
    return best;
}

// gfortran's maxval of s2(0:nt-1) under one(0:nt-1, bit) (set) or its
// .not. (!set): the first value in the mask that is not a NaN, then larger
// ones; NaN when all in the mask are NaN, -HUGE when none is in it.
inline float maxval_bit(const float* s2, int nt, int bit, bool set) {
    int i = 0;
    bool any = false;
    for (; i < nt; ++i) {
        if ((((i >> bit) & 1) != 0) != set) continue;
        any = true;
        if (s2[i] >= -std::numeric_limits<float>::infinity()) break;
    }
    if (i == nt) return any ? std::numeric_limits<float>::quiet_NaN() : -std::numeric_limits<float>::max();
    float m = s2[i];
    for (++i; i < nt; ++i)
        if ((((i >> bit) & 1) != 0) == set && s2[i] > m) m = s2[i];
    return m;
}

// gfortran's maxval of a(0:n-1).
inline float maxval(const float* a, int n) {
    int i = 0;
    for (; i < n; ++i)
        if (a[i] >= -std::numeric_limits<float>::infinity()) break;
    if (i == n) return n > 0 ? std::numeric_limits<float>::quiet_NaN() : -std::numeric_limits<float>::max();
    float m = a[i];
    for (++i; i < n; ++i)
        if (a[i] > m) m = a[i];
    return m;
}

// ---------------------------------------------------------------- trace
// The instrumented references' hashes (scripts/ft8-frontend-patch.py's
// h32: h = (h xor w) * 16777619 over 32-bit words from 2166136261) and
// their line fields.
[[maybe_unused]] uint32_t hw(uint32_t h, uint32_t w) { return (h ^ w) * 16777619u; }
[[maybe_unused]] uint32_t h_r(const float* a, size_t n) {
    uint32_t h = 2166136261u;
    for (size_t i = 0; i < n; ++i) h = hw(h, fbits(a[i]));
    return h;
}
[[maybe_unused]] uint32_t h_c(const cfloat_* a, size_t n) {
    uint32_t h = 2166136261u;
    for (size_t i = 0; i < n; ++i) {
        h = hw(h, fbits(a[i].real()));
        h = hw(h, fbits(a[i].imag()));
    }
    return h;
}
[[maybe_unused]] uint32_t h_i(const int* a, size_t n) {
    uint32_t h = 2166136261u;
    for (size_t i = 0; i < n; ++i) h = hw(h, static_cast<uint32_t>(a[i]));
    return h;
}
[[maybe_unused]] uint32_t h_i1(const int8_t* a, size_t n) {
    uint32_t h = 2166136261u;
    for (size_t i = 0; i < n; ++i) h = hw(h, static_cast<uint32_t>(static_cast<int>(a[i])));
    return h;
}
[[maybe_unused]] std::string X(uint32_t u) {
    char b[16];
    std::snprintf(b, sizeof b, "%08X", u);
    return b;
}
[[maybe_unused]] std::string XF(float f) { return X(fbits(f)); }
[[maybe_unused]] std::string I(long v) { return std::to_string(v); }
[[maybe_unused]] const char* L(bool b) { return b ? "T" : "F"; }

// ---------------------------------------------------------------- baseline
// shell(n, a): a 1-based.
[[maybe_unused]] void shell(int n, float* a) {
    int inc = 1;
    do inc = 3 * inc + 1;
    while (inc <= n);
    do {
        inc = inc / 3;
        for (int i = inc + 1; i <= n; ++i) {
            const float v = a[i];
            int j = i;
            while (a[j - inc] > v) {
                a[j] = a[j - inc];
                j = j - inc;
                if (j <= inc) break;
            }
            a[j] = v;
        }
    } while (inc > 1);
}

// pctile(x, npts, npct, xpct), x 1-based.
[[maybe_unused]] float pctile(const float* x, int npts, int npct) {
    if (npts < 0 || npct < 0 || npct > 100) return 1.0f;
    std::vector<float> tmp(static_cast<size_t>(npts) + 1);
    for (int i = 1; i <= npts; ++i) tmp[i] = x[i];
    shell(npts, tmp.data());
    int j = nint(static_cast<float>(npts) * 0.01f * static_cast<float>(npct));
    if (j < 1) j = 1;
    if (j > npts) j = npts;
    return tmp[j];
}

// determ(array, norder): array(10,10) column-major, modified as there.
[[maybe_unused]] double determ(double* array, int norder) {
    auto A = [&](int i, int j) -> double& { return array[(j - 1) * 10 + (i - 1)]; };
    double det = 1.0;
    for (int k = 1; k <= norder; ++k) {
        if (A(k, k) == 0) {
            int j;
            for (j = k; j <= norder; ++j)
                if (A(k, j) != 0) break;
            if (j > norder) return 0.0;
            for (int i = k; i <= norder; ++i) {
                const double s8 = A(i, j);
                A(i, j) = A(i, k);
                A(i, k) = s8;
            }
            det = -1. * det;
        }
        det = det * A(k, k);
        if (k < norder) {
            const int k1 = k + 1;
            for (int i = k1; i <= norder; ++i)
                for (int j = k1; j <= norder; ++j) A(i, j) = A(i, j) - A(i, k) * A(k, j) / A(k, k);
        }
    }
    return det;
}

// polyfit(x, y, sigmay, npts, nterms, mode=0, a, chisqr): x, y, a 1-based;
// chisqr not computed (nothing reads it).
[[maybe_unused]] void polyfit(const double* x, const double* y, int npts, int nterms, double* a) {
    const int nmax = 2 * nterms - 1;
    double sumx[11] = {}, sumy[11] = {}, array[100];
    for (int i = 1; i <= npts; ++i) {
        const double xi = x[i], yi = y[i];
        const double weight = 1;
        double xterm = weight;
        for (int n = 1; n <= nmax; ++n) {
            sumx[n] = sumx[n] + xterm;
            xterm = xterm * xi;
        }
        double yterm = weight * yi;
        for (int n = 1; n <= nterms; ++n) {
            sumy[n] = sumy[n] + yterm;
            yterm = yterm * xi;
        }
    }
    auto A = [&](int i, int j) -> double& { return array[(j - 1) * 10 + (i - 1)]; };
    for (int j = 1; j <= nterms; ++j)
        for (int k = 1; k <= nterms; ++k) A(j, k) = sumx[j + k - 1];
    const double delta = determ(array, nterms);
    if (delta == 0) {
        for (int l = 1; l <= nterms; ++l) a[l] = 0.;
        return;
    }
    for (int l = 1; l <= nterms; ++l) {
        for (int j = 1; j <= nterms; ++j) {
            for (int k = 1; k <= nterms; ++k) A(j, k) = sumx[j + k - 1];
            A(j, l) = sumy[j];
        }
        a[l] = determ(array, nterms) / delta;
    }
}

// ---------------------------------------------------------------- complex
// sincosf as gfortran calls it for cmplx(cos(x), sin(x)) (glibc 2.39's, carried).
inline cfloat_ cis(float phi) {
    float s, c;
    sincosf_glibc(phi, &s, &c);
    return cfloat_(c, s);
}

// cabsf (glibc: hypotf; 2.39's, carried).
inline float cabs_(cfloat_ z) { return hypotf_glibc(z.real(), z.imag()); }

// Fortran's complex product (no NaN recovery).
inline cfloat_ cmul(cfloat_ a, cfloat_ b) {
    return cfloat_(a.real() * b.real() - a.imag() * b.imag(), a.real() * b.imag() + a.imag() * b.real());
}
inline cfloat_ cadd(cfloat_ a, cfloat_ b) { return cfloat_(a.real() + b.real(), a.imag() + b.imag()); }
inline cfloat_ cscale(cfloat_ a, float r) { return cfloat_(a.real() * r, a.imag() * r); }

// ---------------------------------------------------------------- strings
// Fortran's character*n assignment: cut or blank padded.
[[maybe_unused]] std::string fix(std::string s, size_t n) {
    s.resize(n, ' ');
    return s;
}

// trim(s): trailing blanks only (NULs stay).
[[maybe_unused]] std::string trim(const std::string& s) {
    size_t n = s.size();
    while (n > 0 && s[n - 1] == ' ') --n;
    return s.substr(0, n);
}

// index(s, sub): 1-based, 0 when absent, 1 for an empty sub.
[[maybe_unused]] int findex(const std::string& s, const std::string& sub) {
    const size_t p = s.find(sub);
    return p == std::string::npos ? 0 : static_cast<int>(p) + 1;
}

// isgrid4 (ft8_a7_save's statement function) on g4, its character*4.
[[maybe_unused]] bool isgrid4(const std::string& w) {
    const std::string g4 = fix(w, 4);
    return trim(g4).size() == 4 && g4[0] >= 'A' && g4[0] <= 'R' && g4[1] >= 'A' && g4[1] <= 'R' && g4[2] >= '0' &&
           g4[2] <= '9' && g4[3] >= '0' && g4[3] <= '9';
}

// ---------------------------------------------------------------- maths
// db.f90
[[maybe_unused]] float db(float x) {
    float d = -99.0f;
    if (x > 1.259e-10f) d = 10.0f * log10f_glibc(x);
    return d;
}

// twkfreq1.f90 (ca and cb may be the same array, as ft8_a8d passes them).
[[maybe_unused]] void twkfreq1(const cfloat_* ca, int npts, float fsample, const float a[5], cfloat_* cb) {
    const float twopi = 6.283185307f;
    cfloat_ w(1.0f, 0.0f);
    const float x0 = 0.5f * static_cast<float>(npts + 1);
    const float s = 2.0f / static_cast<float>(npts);
    for (int i = 1; i <= npts; ++i) {
        const float x = s * (static_cast<float>(i) - x0);
        const float p2 = 1.5f * x * x - 0.5f;
        const float p3 = 2.5f * (x * x * x) - 1.5f * x;
        const float p4 = 4.375f * ((x * x) * (x * x)) - 3.75f * (x * x) + 0.375f;
        const float dphi = (a[0] + x * a[1] + p2 * a[2] + p3 * a[3] + p4 * a[4]) * (twopi / fsample);
        float sn, cs;
        sincosf_glibc(dphi, &sn, &cs);
        w = cmul(w, cfloat_(cs, sn));
        cb[i - 1] = cmul(w, ca[i - 1]);
    }
}

}  // namespace
}  // namespace fortran
}  // namespace ft8
