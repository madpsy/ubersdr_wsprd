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
float expf_glibc(float x);

// logf and powf: e_logf.c (e_logf_data.c) and e_powf.c
// (e_powf_log2_data.c, e_exp2f_data.c), without fused operations: their
// FMA and SSE2 builds agree on every float (logf) and on 10^y for every
// float y (powf: FT8's only base).
float logf_glibc(float x);
float powf_glibc(float x, float y);

// hypotf: e_hypotf.c, sqrt((double)x*x + (double)y*y) rounded to float
// (the products are exact in double, so fused or not is the same).
float hypotf_glibc(float x, float y);

// out[l] = f(in[l]), l = 0..3: x86-64 libmvec's _ZGVbN4v_cosf, _ZGVbN4v_sinf
// and _ZGVbN4v_log10f as glibc 2.39 runs them on a CPU with SSE4.1
// (sysdeps/x86_64/fpu/multiarch svml_s_cosf4_core_sse4.S,
// svml_s_sinf4_core_sse4.S, svml_s_log10f4_core_sse4.S), each lane's
// operations in their order, lanes out of their range by the scalar cosf,
// sinf and log10f_glibc as there.
void cosf4_glibc(const float* in, float* out);
void sinf4_glibc(const float* in, float* out);
void log10f4_glibc(const float* in, float* out);

}  // namespace ft8
