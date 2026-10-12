// The vector operations the FST4W decoder uses (fst4w_frontend.cpp,
// fst4w_ldpc.cpp, fst4w.cpp): eight float lanes as one AVX2 register on
// x86-64 (built for x86-64-v3), two NEON registers on aarch64, or a plain
// array elsewhere, every lane computed as the scalar code would compute it
// (multiply then add, never fused, unless a function says otherwise); and
// glibc's vector maths (libmvec) four lanes at a time where the reference's
// gfortran -O3 build calls it, under each platform's name.
//
// libmvec: on x86-64 the reference calls _ZGVbN4v_log10f, _ZGVbN4vv_powf,
// _ZGVbN4v_cosf and _ZGVbN4v_sinf (its disassembly); on aarch64 gfortran's
// math-vector-fortran.h declares cosf, sinf and log10f as SIMD but not powf,
// so there powf stays scalar; without glibc (mingw) everything is scalar.
#pragma once

#include <cmath>
#include <cstdint>

#if defined(__AVX2__) && (defined(__x86_64__) || defined(_M_X64))
#define FST4W_AVX2 1
#include <immintrin.h>
#elif defined(__aarch64__) && defined(__ARM_NEON)
#define FST4W_NEON 1
#include <arm_neon.h>
#endif

#if defined(__GLIBC__) && (defined(__x86_64__) || defined(__aarch64__))
#define FST4W_LIBMVEC 1
#if defined(__x86_64__)
#include <immintrin.h>
extern "C" __m128 _ZGVbN4v_log10f(__m128);
extern "C" __m128 _ZGVbN4vv_powf(__m128, __m128);
extern "C" __m128 _ZGVbN4v_cosf(__m128);
extern "C" __m128 _ZGVbN4v_sinf(__m128);
#else
#include <arm_neon.h>
extern "C" float32x4_t _ZGVnN4v_log10f(float32x4_t);
extern "C" float32x4_t _ZGVnN4v_cosf(float32x4_t);
extern "C" float32x4_t _ZGVnN4v_sinf(float32x4_t);
#endif
#endif

namespace fst4w::simd {

// ---------------------------------------------------------------- libmvec
// out[l] = f(in[l]), l = 0..3, as the platform's reference computes it.
inline void log10f4(const float* in, float* out) {
#if defined(FST4W_LIBMVEC) && defined(__x86_64__)
    _mm_storeu_ps(out, _ZGVbN4v_log10f(_mm_loadu_ps(in)));
#elif defined(FST4W_LIBMVEC)
    vst1q_f32(out, _ZGVnN4v_log10f(vld1q_f32(in)));
#else
    for (int l = 0; l < 4; ++l) out[l] = std::log10(in[l]);
#endif
}
inline void powf4(float base, const float* exp, float* out) {
#if defined(FST4W_LIBMVEC) && defined(__x86_64__)
    _mm_storeu_ps(out, _ZGVbN4vv_powf(_mm_set1_ps(base), _mm_loadu_ps(exp)));
#else
    for (int l = 0; l < 4; ++l) out[l] = std::pow(base, exp[l]);
#endif
}
inline void cosf4(const float* in, float* out) {
#if defined(FST4W_LIBMVEC) && defined(__x86_64__)
    _mm_storeu_ps(out, _ZGVbN4v_cosf(_mm_loadu_ps(in)));
#elif defined(FST4W_LIBMVEC)
    vst1q_f32(out, _ZGVnN4v_cosf(vld1q_f32(in)));
#else
    for (int l = 0; l < 4; ++l) out[l] = std::cos(in[l]);
#endif
}
inline void sinf4(const float* in, float* out) {
#if defined(FST4W_LIBMVEC) && defined(__x86_64__)
    _mm_storeu_ps(out, _ZGVbN4v_sinf(_mm_loadu_ps(in)));
#elif defined(FST4W_LIBMVEC)
    vst1q_f32(out, _ZGVnN4v_sinf(vld1q_f32(in)));
#else
    for (int l = 0; l < 4; ++l) out[l] = std::sin(in[l]);
#endif
}

// ---------------------------------------------------------------- eight lanes
#if defined(FST4W_AVX2)
typedef __m256 V8;
inline V8 zero() { return _mm256_setzero_ps(); }
inline V8 set1(float x) { return _mm256_set1_ps(x); }
inline V8 load(const float* p) { return _mm256_loadu_ps(p); }
inline void store(float* p, V8 x) { _mm256_storeu_ps(p, x); }
inline V8 add(V8 a, V8 b) { return _mm256_add_ps(a, b); }
inline V8 sub(V8 a, V8 b) { return _mm256_sub_ps(a, b); }
inline V8 mul(V8 a, V8 b) { return _mm256_mul_ps(a, b); }
inline V8 div(V8 a, V8 b) { return _mm256_div_ps(a, b); }
inline V8 max(V8 a, V8 b) { return _mm256_max_ps(a, b); }
// a*b + c, fused where the platform can (not exact: only for sums whose
// bits do not matter).
inline V8 fma_loose(V8 a, V8 b, V8 c) { return _mm256_fmadd_ps(a, b, c); }
// x where the mask (eight 0 / ~0 words) is set, else 0.
inline V8 and_mask(V8 x, const uint32_t* m) { return _mm256_and_ps(x, _mm256_load_ps(reinterpret_cast<const float*>(m))); }
// acc + p[l] for the lanes l whose bit is set in byte (lut[byte]: its eight
// lane masks, as and_mask takes them).
inline V8 add_masked_byte(V8 acc, const float* p, unsigned byte, const uint32_t* lut) {
    return add(acc, and_mask(load(p), lut + 8 * byte));
}
// The lanes' sum, in no particular order.
inline float hsum(V8 x) {
    __m128 v = _mm_add_ps(_mm256_castps256_ps128(x), _mm256_extractf128_ps(x, 1));
    v = _mm_add_ps(v, _mm_movehl_ps(v, v));
    v = _mm_add_ss(v, _mm_shuffle_ps(v, v, 1));
    return _mm_cvtss_f32(v);
}
// cabsf(x + iy) as glibc's hypotf: the double sum of the exact double
// squares, its square root rounded to float.
inline V8 abs_exact(V8 x, V8 y) {
    const __m256d xl = _mm256_cvtps_pd(_mm256_castps256_ps128(x)), xh = _mm256_cvtps_pd(_mm256_extractf128_ps(x, 1));
    const __m256d yl = _mm256_cvtps_pd(_mm256_castps256_ps128(y)), yh = _mm256_cvtps_pd(_mm256_extractf128_ps(y, 1));
    const __m256d sl = _mm256_sqrt_pd(_mm256_add_pd(_mm256_mul_pd(xl, xl), _mm256_mul_pd(yl, yl)));
    const __m256d sh = _mm256_sqrt_pd(_mm256_add_pd(_mm256_mul_pd(xh, xh), _mm256_mul_pd(yh, yh)));
    return _mm256_set_m128(_mm256_cvtpd_ps(sh), _mm256_cvtpd_ps(sl));
}
#elif defined(FST4W_NEON)
struct V8 {
    float32x4_t lo, hi;
};
inline V8 zero() { return {vdupq_n_f32(0.0f), vdupq_n_f32(0.0f)}; }
inline V8 set1(float x) { return {vdupq_n_f32(x), vdupq_n_f32(x)}; }
inline V8 load(const float* p) { return {vld1q_f32(p), vld1q_f32(p + 4)}; }
inline void store(float* p, V8 x) {
    vst1q_f32(p, x.lo);
    vst1q_f32(p + 4, x.hi);
}
inline V8 add(V8 a, V8 b) { return {vaddq_f32(a.lo, b.lo), vaddq_f32(a.hi, b.hi)}; }
inline V8 sub(V8 a, V8 b) { return {vsubq_f32(a.lo, b.lo), vsubq_f32(a.hi, b.hi)}; }
inline V8 mul(V8 a, V8 b) { return {vmulq_f32(a.lo, b.lo), vmulq_f32(a.hi, b.hi)}; }
inline V8 div(V8 a, V8 b) { return {vdivq_f32(a.lo, b.lo), vdivq_f32(a.hi, b.hi)}; }
inline V8 max(V8 a, V8 b) { return {vmaxq_f32(a.lo, b.lo), vmaxq_f32(a.hi, b.hi)}; }
inline V8 fma_loose(V8 a, V8 b, V8 c) { return {vfmaq_f32(c.lo, a.lo, b.lo), vfmaq_f32(c.hi, a.hi, b.hi)}; }
inline V8 and_mask(V8 x, const uint32_t* m) {
    return {vreinterpretq_f32_u32(vandq_u32(vreinterpretq_u32_f32(x.lo), vld1q_u32(m))),
            vreinterpretq_f32_u32(vandq_u32(vreinterpretq_u32_f32(x.hi), vld1q_u32(m + 4)))};
}
inline float hsum(V8 x) { return vaddvq_f32(vaddq_f32(x.lo, x.hi)); }
// The byte's lane masks made in registers (vtst against the lanes' bits):
// no table to load, which costs NEON more than AVX2.
inline V8 add_masked_byte(V8 acc, const float* p, unsigned byte, const uint32_t*) {
    static const uint32_t lo_bits[4] = {1, 2, 4, 8}, hi_bits[4] = {16, 32, 64, 128};
    const uint32x4_t b = vdupq_n_u32(byte);
    const uint32x4_t ml = vtstq_u32(b, vld1q_u32(lo_bits)), mh = vtstq_u32(b, vld1q_u32(hi_bits));
    return {vaddq_f32(acc.lo, vreinterpretq_f32_u32(vandq_u32(vreinterpretq_u32_f32(vld1q_f32(p)), ml))),
            vaddq_f32(acc.hi, vreinterpretq_f32_u32(vandq_u32(vreinterpretq_u32_f32(vld1q_f32(p + 4)), mh)))};
}
inline float32x4_t abs4(float32x4_t x, float32x4_t y) {
    const float64x2_t xl = vcvt_f64_f32(vget_low_f32(x)), xh = vcvt_high_f64_f32(x);
    const float64x2_t yl = vcvt_f64_f32(vget_low_f32(y)), yh = vcvt_high_f64_f32(y);
    const float64x2_t sl = vsqrtq_f64(vaddq_f64(vmulq_f64(xl, xl), vmulq_f64(yl, yl)));
    const float64x2_t sh = vsqrtq_f64(vaddq_f64(vmulq_f64(xh, xh), vmulq_f64(yh, yh)));
    return vcvt_high_f32_f64(vcvt_f32_f64(sl), sh);
}
inline V8 abs_exact(V8 x, V8 y) { return {abs4(x.lo, y.lo), abs4(x.hi, y.hi)}; }
#else
struct V8 {
    float v[8];
};
inline V8 zero() { return {{0, 0, 0, 0, 0, 0, 0, 0}}; }
inline V8 set1(float x) { return {{x, x, x, x, x, x, x, x}}; }
inline V8 load(const float* p) {
    V8 r;
    for (int l = 0; l < 8; ++l) r.v[l] = p[l];
    return r;
}
inline void store(float* p, V8 x) {
    for (int l = 0; l < 8; ++l) p[l] = x.v[l];
}
#define FST4W_V8_OP(name, expr)                         \
    inline V8 name(V8 a, V8 b) {                        \
        V8 r;                                           \
        for (int l = 0; l < 8; ++l) r.v[l] = (expr);    \
        return r;                                       \
    }
FST4W_V8_OP(add, a.v[l] + b.v[l])
FST4W_V8_OP(sub, a.v[l] - b.v[l])
FST4W_V8_OP(mul, a.v[l] * b.v[l])
FST4W_V8_OP(div, a.v[l] / b.v[l])
FST4W_V8_OP(max, a.v[l] > b.v[l] ? a.v[l] : b.v[l])
#undef FST4W_V8_OP
inline V8 fma_loose(V8 a, V8 b, V8 c) {
    V8 r;
    for (int l = 0; l < 8; ++l) r.v[l] = a.v[l] * b.v[l] + c.v[l];
    return r;
}
inline V8 and_mask(V8 x, const uint32_t* m) {
    V8 r;
    for (int l = 0; l < 8; ++l) r.v[l] = m[l] ? x.v[l] : 0.0f;
    return r;
}
inline V8 add_masked_byte(V8 acc, const float* p, unsigned byte, const uint32_t*) {
    for (int l = 0; l < 8; ++l)
        if ((byte >> l) & 1) acc.v[l] += p[l];
    return acc;
}
inline float hsum(V8 x) {
    float s = 0;
    for (int l = 0; l < 8; ++l) s += x.v[l];
    return s;
}
inline V8 abs_exact(V8 x, V8 y) {
    V8 r;
    for (int l = 0; l < 8; ++l)
        r.v[l] = static_cast<float>(std::sqrt(static_cast<double>(x.v[l]) * x.v[l] + static_cast<double>(y.v[l]) * y.v[l]));
    return r;
}
#endif

}  // namespace fst4w::simd
