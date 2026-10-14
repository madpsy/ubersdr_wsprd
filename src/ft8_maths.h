// The maths FT8's front end takes from glibc where its results differ
// between platforms: carried here, as glibc 2.39 computes them on x86-64
// (the build WSJT-X 567ad29ce's fixtures were made with), so that every
// platform gets the reference's bits whatever its own libm and libmvec.
// Measured on every float: glibc 2.41 (on a Raspberry Pi 5) has
// CORE-MATH's correctly rounded erff and log10f (they differ from 2.39's
// in the last bit or two, 0.2-34% of the floats in the ranges
// tried), and AArch64's libmvec cosf, sinf and log10f are other algorithms
// (up to 2 ulp, 15% of FT8's table entries). mingw's libm differs again
// (its powf, cosf and sinf), so every maths function the FT8
// decode calls is carried: also glibc 2.39's scalar sinf, cosf, sincosf,
// expf, logf, powf and hypotf. tanhf is ft8_ldpc.h's tanhf_glibc. Left to
// the platform only what IEEE 754 defines exactly: sqrtf, fabsf, fmodf,
// lroundf.
//
// Each is checked equal to glibc 2.39's own on every float on x86-64
// (tests/test_ft8_maths.cpp, [ft8-maths-libm]).
#pragma once

namespace ft8 {

// erff and log10f: sysdeps/ieee754/flt-32 s_erff.c and e_log10f.c
// (fdlibm's), line for line, with the platform's expf and logf.
float erff_glibc(float x);
float log10f_glibc(float x);

// sinf, cosf, sincosf and expf: sysdeps/ieee754/flt-32 s_sinf.c, s_cosf.c,
// s_sincosf.c (with s_sincosf.h, sincosf_poly.h, s_sincosf_data.c) and
// e_expf.c (e_exp2f_data.c) as x86-64 glibc runs them on a CPU with FMA
// (sysdeps/x86_64/fpu/multiarch's *-fma.c builds, -mfma -mavx2): every
// multiply-add GCC fused there is a std::fma here, the rest rounded as
// written. The SSE2 builds give other bits on some floats in [16, 120)
// (sinf, cosf, sincosf) and [32, 64) and near -60 (expf).
float sinf_glibc(float x);
float cosf_glibc(float x);
void sincosf_glibc(float x, float* s, float* c);
// sincosf_glibc(x[i], &s[i], &c[i]), i = 0..n-1, the same bits (FT4's
// gen_ft4wave): four at a time on x86-64 with AVX2 and FMA and on
// aarch64 with NEON ([ft8-maths]: "sincosf four at a time").
void sincosf_n_glibc(const float* x, float* s, float* c, int n);
float expf_glibc(float x);

// logf and powf: e_logf.c (e_logf_data.c) and e_powf.c
// (e_powf_log2_data.c, e_exp2f_data.c), without fused operations: their
// FMA and SSE2 builds agree on every float (logf), on 10^y for every
// float y (powf: FT8's only base, and FT4's ft4_baseline) and on x^0.5 for
// every float x (FT4's sync4d).
float logf_glibc(float x);
float powf_glibc(float x, float y);
// out[i] = powf_glibc(x[i], 0.5f), i = 0..n-1, the same bits (FT4's
// sync4d): four at a time on x86-64 with AVX2 and on aarch64 with
// NEON, powf_glibc's operations in double lanes for positive
// normal x, the rest by
// powf_glibc ([ft8-maths]: "powf(x, 0.5) four at a time").
void powf_half_glibc(const float* x, float* out, int n);

// hypotf: e_hypotf.c, sqrt((double)x*x + (double)y*y) rounded to float
// (the products are exact in double, so fused or not is the same).
float hypotf_glibc(float x, float y);
// out[i] = hypotf_glibc(x[i], y[i]), i = 0..n-1, the same bits (FT4's
// bit metrics): four at a time on x86-64 with AVX2 when finite
// (aarch64: one at a time; NEON gained nothing there)
// ([ft8-maths]: "hypotf four at a time").
void hypotf_n_glibc(const float* x, const float* y, float* out, int n);
// The same from interleaved pairs: out[i] = hypotf_glibc(xy[2i], xy[2i+1])
// (complex magnitudes, cabsf's bits; [ft8-maths]: "hypotf four at a time,
// interleaved").
void hypotf_c_glibc(const float* xy, float* out, int n);

// out[l] = f(in[l]), l = 0..3: x86-64 libmvec's _ZGVbN4v_cosf, _ZGVbN4v_sinf
// and _ZGVbN4v_log10f as glibc 2.39 runs them on a CPU with SSE4.1
// (sysdeps/x86_64/fpu/multiarch svml_s_cosf4_core_sse4.S,
// svml_s_sinf4_core_sse4.S, svml_s_log10f4_core_sse4.S), each lane's
// operations in their order, lanes out of their range by the scalar cosf,
// sinf and log10f_glibc as there.
void cosf4_glibc(const float* in, float* out);
void sinf4_glibc(const float* in, float* out);
void log10f4_glibc(const float* in, float* out);

// out[l] = x[l]^y[l], l = 0..3: x86-64 libmvec's _ZGVbN4vv_powf as glibc
// 2.39 runs it on a CPU with SSE4.1 (svml_s_powf4_core_sse4.S), lanes out
// of its range by the scalar powf_glibc as there; its rcpps as AMD Zen 2
// and 3 give it (another x86's may round other mantissas differently: see
// ft8_maths.cpp). FT4's ft4_baseline calls it with x = 10 in every lane
// (checked against glibc's on every float y there).
void powf4_glibc(const float* x, const float* y, float* out);

}  // namespace ft8
