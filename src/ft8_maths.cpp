// glibc 2.39's erff, log10f, sinf, cosf, sincosf, expf, logf, powf,
// hypotf and x86-64 libmvec cosf, sinf, log10f and powf (ft8_maths.h), carried. Built without contraction (CMakeLists.txt), as
// glibc's x86-64 build of them has none: every multiply and add rounded as
// written.
//
// s_erff.c and e_log10f.c: Sun Microsystems (fdlibm; freely
// distributable). s_sinf.c, s_cosf.c, s_sincosf.c, e_expf.c, e_logf.c,
// e_powf.c and their data: ARM Ltd. and the Free Software Foundation (GNU
// LGPL v2.1 or later); e_hypotf.c: the Free Software Foundation (GNU LGPL
// v2.1 or later). svml_s_cosf4_core_sse4.S, svml_s_sinf4_core_sse4.S,
// svml_s_log10f4_core_sse4.S, svml_s_powf4_core_sse4.S and
// svml_s_powf_data.S: the GNU C Library, Free Software Foundation
// (GNU LGPL v2.1 or later).
#include "ft8_maths.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>

#if defined(__AVX2__) && (defined(__x86_64__) || defined(_M_X64))
#include <immintrin.h>
#define FT8_MATHS_AVX2 1
#elif defined(__aarch64__) && defined(__ARM_NEON)
// NEON beside each four-wide AVX2 path, every lane the same
// operations in the same order (double lanes two to a register, four as
// two), no fused multiply-add but the explicit std::fma the scalar code has
// (vfmaq_f64 there); table reads lane by lane where AVX2 gathers.
#include <arm_neon.h>
#define FT8_MATHS_NEON 1
#endif

namespace ft8 {

namespace {

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

// svml_s_trig_data.S
const float sAbsRangeReductionVal = bitsf(0x461c4000);  // threshold for out-of-range values
const float sPI1 = bitsf(0x40490000), sPI2 = bitsf(0x3a7da000), sPI3 = bitsf(0x34222000),
            sPI4 = bitsf(0x2cb4611a);  // PI in four parts
const float sA3 = bitsf(0xbe2aaaa6), sA5 = bitsf(0x3c08876a), sA7 = bitsf(0xb94fb7ff),
            sA9 = bitsf(0x362edef8);  // polynomial, without FMA
const float sInvPI = bitsf(0x3ea2f983), sRShifter = bitsf(0x4b400000), sHalfPI = bitsf(0x3fc90fdb),
            sOneHalf = bitsf(0x3f000000);

// The polynomial and sign of both: R = X + X * X^2 * (A3 + x^2 * (A5 + ...,
// the sign S set on X after X^2.
float trig_poly(float r, uint32_t s) {
    const float x2 = r * r;
    const float rs = bitsf(fbits(r) ^ s);
    float p = sA9 * x2;
    p = p + sA7;
    p = p * x2;
    p = p + sA5;
    p = p * x2;
    p = p + sA3;
    float q = x2 * p;
    q = q * rs;
    return rs + q;
}

// ---- glibc 2.39 sysdeps/ieee754/flt-32: sinf, cosf, sincosf, expf, logf,
// powf, hypotf (math_config.h with TOINT_INTRINSICS 0, WANT_ROUNDING and
// WANT_ERRNO 1, errno not set here) ----

uint64_t asuint64(double d) {
    uint64_t u;
    std::memcpy(&u, &d, 8);
    return u;
}
double asdouble(uint64_t u) {
    double d;
    std::memcpy(&d, &u, 8);
    return d;
}

// math_errf.c (xflowf's (sign ? -y : y) * y; invalid's (x - x) / (x - x)).
float xflowf(uint32_t sign, float y) { return (sign ? -y : y) * y; }
float math_oflowf(uint32_t sign) { return xflowf(sign, 0x1p97f); }
float math_uflowf(uint32_t sign) { return xflowf(sign, 0x1p-95f); }
float math_may_uflowf(uint32_t sign) { return xflowf(sign, 0x1.4p-75f); }
float math_divzerof(uint32_t sign) {
    volatile float y = 0;
    return static_cast<float>(sign ? -1 : 1) / y;
}
float math_invalidf(float x) {
    volatile float d = x - x;
    return d / d;
}

// s_sincosf_data.c, sincosf_poly.h's sincos_t: the cosine polynomial is
// negated in the second entry.
struct sincos_t {
    double sign[4];          // Sign of sine in quadrants 0..3.
    double hpi_inv;          // 2 / PI * 2^24 (TOINT_INTRINSICS 0).
    double hpi;              // PI / 2.
    double c0, c1, c2, c3, c4;  // Cosine polynomial.
    double s1, s2, s3;          // Sine polynomial.
};
const sincos_t sincosf_table[2] = {
    {{1.0, -1.0, -1.0, 1.0}, 0x1.45F306DC9C883p+23, 0x1.921FB54442D18p0, 0x1p0, -0x1.ffffffd0c621cp-2,
     0x1.55553e1068f19p-5, -0x1.6c087e89a359dp-10, 0x1.99343027bf8c3p-16, -0x1.555545995a603p-3,
     0x1.1107605230bc4p-7, -0x1.994eb3774cf24p-13},
    {{1.0, -1.0, -1.0, 1.0}, 0x1.45F306DC9C883p+23, 0x1.921FB54442D18p0, -0x1p0, 0x1.ffffffd0c621cp-2,
     -0x1.55553e1068f19p-5, 0x1.6c087e89a359dp-10, -0x1.99343027bf8c3p-16, -0x1.555545995a603p-3,
     0x1.1107605230bc4p-7, -0x1.994eb3774cf24p-13}};
// Table with 4/PI to 192 bit precision.
const uint32_t inv_pio4[24] = {0xa2,       0xa2f9,     0xa2f983,   0xa2f9836e, 0xf9836e4e, 0x836e4e44,
                               0x6e4e4415, 0x4e441529, 0x441529fc, 0x1529fc27, 0x29fc2757, 0xfc2757d1,
                               0x2757d1f5, 0x57d1f534, 0xd1f534dd, 0xf534ddc0, 0x34ddc0db, 0xddc0db62,
                               0xc0db6295, 0xdb629599, 0x6295993c, 0x95993c43, 0x993c4390, 0x3c439041};
const double pi63 = 0x1.921FB54442D18p-62;  // 2PI * 2^-64
const float pio4 = 0x1.921FB6p-1f;          // PI / 4

// Top 12 bits of the float representation with the sign bit cleared.
uint32_t abstop12(float x) { return (fbits(x) >> 20) & 0x7ff; }

// The sine and cosine of X and X2 (X squared) by the polynomial P into
// SINP and COSP, swapped if N is odd. Fused as the FMA build: each
// "a + b * c" whose product has no other use.
void sincosf_poly(double x, double x2, const sincos_t* p, int n, float* sinp, float* cosp) {
    const double x4 = x2 * x2;
    const double x3 = x2 * x;
    const double c2 = std::fma(x2, p->c4, p->c3);
    const double s1 = std::fma(x2, p->s3, p->s2);
    // Swap sin/cos result based on quadrant.
    float* tmp = (n & 1 ? cosp : sinp);
    cosp = (n & 1 ? sinp : cosp);
    sinp = tmp;
    const double c1 = std::fma(x2, p->c1, p->c0);
    const double x5 = x3 * x2;
    const double x6 = x4 * x2;
    const double s = std::fma(x3, p->s1, x);
    const double c = std::fma(x4, p->c2, c1);
    *sinp = static_cast<float>(std::fma(x5, s1, s));
    *cosp = static_cast<float>(std::fma(x6, c2, c));
}

// The sine of X and X2 by P; N odd: the cosine polynomial.
float sinf_poly(double x, double x2, const sincos_t* p, int n) {
    if ((n & 1) == 0) {
        const double x3 = x * x2;
        const double s1 = std::fma(x2, p->s3, p->s2);
        const double x7 = x3 * x2;
        const double s = std::fma(x3, p->s1, x);
        return static_cast<float>(std::fma(x7, s1, s));
    }
    const double x4 = x2 * x2;
    const double c2 = std::fma(x2, p->c4, p->c3);
    const double c1 = std::fma(x2, p->c1, p->c0);
    const double x6 = x4 * x2;
    const double c = std::fma(x4, p->c2, c1);
    return static_cast<float>(std::fma(x6, c2, c));
}

// Fast range reduction using single multiply-subtract (fused): the modulo
// of X between -PI/4 and PI/4, the quadrant in NP; for |X| <= 120.
double reduce_fast(double x, const sincos_t* p, int* np) {
    const double r = x * p->hpi_inv;
    const int n = (static_cast<int32_t>(r) + 0x800000) >> 24;
    *np = n;
    return std::fma(-static_cast<double>(n), p->hpi, x);
}

// Reduce XI (a reinterpreted float >= 2.0f, sign ignored) to a multiple of
// PI/2 by integer arithmetic on 4/PI to 192 bits.
double reduce_large(uint32_t xi, int* np) {
    const uint32_t* arr = &inv_pio4[(xi >> 26) & 15];
    const int shift = (xi >> 23) & 7;
    uint64_t n, res0, res1, res2;
    xi = (xi & 0xffffff) | 0x800000;
    xi <<= shift;
    res0 = static_cast<uint32_t>(xi * arr[0]);
    res1 = static_cast<uint64_t>(xi) * arr[4];
    res2 = static_cast<uint64_t>(xi) * arr[8];
    res0 = (res2 >> 32) | (res0 << 32);
    res0 += res1;
    n = (res0 + (1ULL << 61)) >> 62;
    res0 -= n << 62;
    const double x = static_cast<double>(static_cast<int64_t>(res0));
    *np = static_cast<int>(n);
    return x * pi63;
}

// e_exp2f_data.c
constexpr int EXP2F_TABLE_BITS = 5;
constexpr int EXP2F_N = 1 << EXP2F_TABLE_BITS;
const uint64_t exp2f_tab[EXP2F_N] = {
    0x3ff0000000000000, 0x3fefd9b0d3158574, 0x3fefb5586cf9890f, 0x3fef9301d0125b51, 0x3fef72b83c7d517b,
    0x3fef54873168b9aa, 0x3fef387a6e756238, 0x3fef1e9df51fdee1, 0x3fef06fe0a31b715, 0x3feef1a7373aa9cb,
    0x3feedea64c123422, 0x3feece086061892d, 0x3feebfdad5362a27, 0x3feeb42b569d4f82, 0x3feeab07dd485429,
    0x3feea47eb03a5585, 0x3feea09e667f3bcd, 0x3fee9f75e8ec5f74, 0x3feea11473eb0187, 0x3feea589994cce13,
    0x3feeace5422aa0db, 0x3feeb737b0cdc5e5, 0x3feec49182a3f090, 0x3feed503b23e255d, 0x3feee89f995ad3ad,
    0x3feeff76f2fb5e47, 0x3fef199bdd85529c, 0x3fef3720dcef9069, 0x3fef5818dcfba487, 0x3fef7c97337b9b5f,
    0x3fefa4afa2a490da, 0x3fefd0765b6e4540};
const double exp2f_shift_scaled = 0x1.8p+52 / EXP2F_N;
const double exp2f_poly[3] = {0x1.c6af84b912394p-5, 0x1.ebfce50fac4f3p-3, 0x1.62e42ff0c52d6p-1};
const double exp2f_shift = 0x1.8p+52;
const double exp2f_invln2_scaled = 0x1.71547652b82fep+0 * EXP2F_N;
const double exp2f_poly_scaled[3] = {0x1.c6af84b912394p-5 / EXP2F_N / EXP2F_N / EXP2F_N,
                                     0x1.ebfce50fac4f3p-3 / EXP2F_N / EXP2F_N, 0x1.62e42ff0c52d6p-1 / EXP2F_N};

uint32_t top12(float x) { return fbits(x) >> 20; }

// e_logf_data.c
constexpr int LOGF_TABLE_BITS = 4;
struct invc_logc {
    double invc, logc;
};
const invc_logc logf_tab[1 << LOGF_TABLE_BITS] = {
    {0x1.661ec79f8f3bep+0, -0x1.57bf7808caadep-2}, {0x1.571ed4aaf883dp+0, -0x1.2bef0a7c06ddbp-2},
    {0x1.49539f0f010bp+0, -0x1.01eae7f513a67p-2},  {0x1.3c995b0b80385p+0, -0x1.b31d8a68224e9p-3},
    {0x1.30d190c8864a5p+0, -0x1.6574f0ac07758p-3}, {0x1.25e227b0b8eap+0, -0x1.1aa2bc79c81p-3},
    {0x1.1bb4a4a1a343fp+0, -0x1.a4e76ce8c0e5ep-4}, {0x1.12358f08ae5bap+0, -0x1.1973c5a611cccp-4},
    {0x1.0953f419900a7p+0, -0x1.252f438e10c1ep-5}, {0x1p+0, 0x0p+0},
    {0x1.e608cfd9a47acp-1, 0x1.aa5aa5df25984p-5},  {0x1.ca4b31f026aap-1, 0x1.c5e53aa362eb4p-4},
    {0x1.b2036576afce6p-1, 0x1.526e57720db08p-3},  {0x1.9c2d163a1aa2dp-1, 0x1.bc2860d22477p-3},
    {0x1.886e6037841edp-1, 0x1.1058bc8a07ee1p-2},  {0x1.767dcf5534862p-1, 0x1.4043057b6ee09p-2}};
const double logf_ln2 = 0x1.62e42fefa39efp-1;
const double logf_poly[3] = {-0x1.00ea348b88334p-2, 0x1.5575b0be00b6ap-2, -0x1.ffffef20a4123p-2};

// e_powf_log2_data.c (POWF_SCALE 1)
constexpr int POWF_LOG2_TABLE_BITS = 4;
const invc_logc powf_log2_tab[1 << POWF_LOG2_TABLE_BITS] = {
    {0x1.661ec79f8f3bep+0, -0x1.efec65b963019p-2}, {0x1.571ed4aaf883dp+0, -0x1.b0b6832d4fca4p-2},
    {0x1.49539f0f010bp+0, -0x1.7418b0a1fb77bp-2},  {0x1.3c995b0b80385p+0, -0x1.39de91a6dcf7bp-2},
    {0x1.30d190c8864a5p+0, -0x1.01d9bf3f2b631p-2}, {0x1.25e227b0b8eap+0, -0x1.97c1d1b3b7afp-3},
    {0x1.1bb4a4a1a343fp+0, -0x1.2f9e393af3c9fp-3}, {0x1.12358f08ae5bap+0, -0x1.960cbbf788d5cp-4},
    {0x1.0953f419900a7p+0, -0x1.a6f9db6475fcep-5}, {0x1p+0, 0x0p+0},
    {0x1.e608cfd9a47acp-1, 0x1.338ca9f24f53dp-4},  {0x1.ca4b31f026aap-1, 0x1.476a9543891bap-3},
    {0x1.b2036576afce6p-1, 0x1.e840b4ac4e4d2p-3},  {0x1.9c2d163a1aa2dp-1, 0x1.40645f0c6651cp-2},
    {0x1.886e6037841edp-1, 0x1.88e9c2c1b9ff8p-2},  {0x1.767dcf5534862p-1, 0x1.ce0a44eb17bccp-2}};
const double powf_log2_poly[5] = {0x1.27616c9496e0bp-2, -0x1.71969a075c67ap-2, 0x1.ec70a6ca7baddp-2,
                                  -0x1.7154748bef6c8p-1, 0x1.71547652ab82bp0};

// e_powf.c's log2_inline: subnormal input normalized so ix has a negative
// biased exponent.
double powf_log2_inline(uint32_t ix) {
    constexpr uint32_t OFF = 0x3f330000;
    constexpr int N = 1 << POWF_LOG2_TABLE_BITS;
    // x = 2^k z; where z is in range [OFF,2*OFF] and exact.
    const uint32_t tmp = ix - OFF;
    const int i = static_cast<int>((tmp >> (23 - POWF_LOG2_TABLE_BITS)) % N);
    const uint32_t top = tmp & 0xff800000;
    const uint32_t iz = ix - top;
    const int k = static_cast<int32_t>(top) >> 23;  // arithmetic shift
    const double invc = powf_log2_tab[i].invc;
    const double logc = powf_log2_tab[i].logc;
    const double z = static_cast<double>(bitsf(iz));
    // log2(x) = log1p(z/c-1)/ln2 + log2(c) + k
    const double r = z * invc - 1;
    const double y0 = logc + static_cast<double>(k);
    const double* A = powf_log2_poly;
    const double r2 = r * r;
    double y = A[0] * r + A[1];
    const double p = A[2] * r + A[3];
    const double r4 = r2 * r2;
    double q = A[4] * r + y0;
    q = p * r2 + q;
    y = y * r4 + q;
    return y;
}

// e_powf.c's exp2_inline: xd unscaled in [-1021,1023], sign_bias the
// result's sign.
double powf_exp2_inline(double xd, uint32_t sign_bias) {
    const double* C = exp2f_poly;
    // x = k/N + r with r in [-1/(2N), 1/(2N)]
    double kd = static_cast<double>(xd + exp2f_shift_scaled);  // Rounding to double precision is required.
    const uint64_t ki = asuint64(kd);
    kd -= exp2f_shift_scaled;  // k/N
    const double r = xd - kd;
    // exp2(x) = 2^(k/N) * 2^r ~= s * (C0*r^3 + C1*r^2 + C2*r + 1)
    uint64_t t = exp2f_tab[ki % EXP2F_N];
    const uint64_t ski = ki + sign_bias;
    t += ski << (52 - EXP2F_TABLE_BITS);
    const double s = asdouble(t);
    const double z = C[0] * r + C[1];
    const double r2 = r * r;
    double y = C[2] * r + 1;
    y = z * r2 + y;
    y = y * s;
    return y;
}

// 0 if not int, 1 if odd int, 2 if even int (iy: a non-zero finite float).
int powf_checkint(uint32_t iy) {
    const int e = iy >> 23 & 0xff;
    if (e < 0x7f) return 0;
    if (e > 0x7f + 23) return 2;
    if (iy & ((1u << (0x7f + 23 - e)) - 1)) return 0;
    if (iy & (1u << (0x7f + 23 - e))) return 1;
    return 2;
}

bool zeroinfnan(uint32_t ix) { return 2 * ix - 1 >= 2u * 0x7f800000 - 1; }

bool issignalingf(float x) {
    const uint32_t ix = fbits(x);
    return (ix & 0x7fc00000u) == 0x7f800000u && (ix & 0x003fffffu) != 0;
}

}  // namespace

float sinf_glibc(float y) {
    double x = y;
    double s;
    int n;
    const sincos_t* p = &sincosf_table[0];
    if (abstop12(y) < abstop12(pio4)) {
        s = x * x;
        if (abstop12(y) < abstop12(0x1p-12f)) return y;
        return sinf_poly(x, s, p, 0);
    } else if (abstop12(y) < abstop12(120.0f)) {
        x = reduce_fast(x, p, &n);
        // Setup the signs for sin and cos.
        s = p->sign[n & 3];
        if (n & 2) p = &sincosf_table[1];
        return sinf_poly(x * s, x * x, p, n);
    } else if (abstop12(y) < abstop12(INFINITY)) {
        const uint32_t xi = fbits(y);
        const int sign = static_cast<int>(xi >> 31);
        x = reduce_large(xi, &n);
        // Setup signs for sin and cos - include original sign.
        s = p->sign[(n + sign) & 3];
        if ((n + sign) & 2) p = &sincosf_table[1];
        return sinf_poly(x * s, x * x, p, n);
    }
    return math_invalidf(y);
}

float cosf_glibc(float y) {
    double x = y;
    double s;
    int n;
    const sincos_t* p = &sincosf_table[0];
    if (abstop12(y) < abstop12(pio4)) {
        const double x2 = x * x;
        if (abstop12(y) < abstop12(0x1p-12f)) return 1.0f;
        return sinf_poly(x, x2, p, 1);
    } else if (abstop12(y) < abstop12(120.0f)) {
        x = reduce_fast(x, p, &n);
        s = p->sign[n & 3];
        if (n & 2) p = &sincosf_table[1];
        return sinf_poly(x * s, x * x, p, n ^ 1);
    } else if (abstop12(y) < abstop12(INFINITY)) {
        const uint32_t xi = fbits(y);
        const int sign = static_cast<int>(xi >> 31);
        x = reduce_large(xi, &n);
        s = p->sign[(n + sign) & 3];
        if ((n + sign) & 2) p = &sincosf_table[1];
        return sinf_poly(x * s, x * x, p, n ^ 1);
    }
    return math_invalidf(y);
}

void sincosf_glibc(float y, float* sinp, float* cosp) {
    double x = y;
    double s;
    int n;
    const sincos_t* p = &sincosf_table[0];
    if (abstop12(y) < abstop12(pio4)) {
        const double x2 = x * x;
        if (abstop12(y) < abstop12(0x1p-12f)) {
            *sinp = y;
            *cosp = 1.0f;
            return;
        }
        sincosf_poly(x, x2, p, 0, sinp, cosp);
    } else if (abstop12(y) < abstop12(120.0f)) {
        x = reduce_fast(x, p, &n);
        s = p->sign[n & 3];
        if (n & 2) p = &sincosf_table[1];
        sincosf_poly(x * s, x * x, p, n, sinp, cosp);
    } else if (abstop12(y) < abstop12(INFINITY)) {
        const uint32_t xi = fbits(y);
        const int sign = static_cast<int>(xi >> 31);
        x = reduce_large(xi, &n);
        s = p->sign[(n + sign) & 3];
        if ((n + sign) & 2) p = &sincosf_table[1];
        sincosf_poly(x * s, x * x, p, n, sinp, cosp);
    } else {
        // Return NaN if Inf or NaN for both sin and cos.
        *sinp = *cosp = y - y;
    }
}

void sincosf_n_glibc(const float* x, float* sinp, float* cosp, int n) {
    int i = 0;
#if defined(FT8_MATHS_AVX2) && defined(__FMA__)
    // Four at a time in double lanes when every |x| is below 120 (by the
    // top 12 bits, as sincosf_glibc tests it): its reduce_fast path with
    // the same fused operations, which is also its |x| < pi/4 path
    // exactly (n = 0 there: x - 0*hpi = x, sign 1, the first table), and
    // the tiny lanes (|x| < 2^-12) put in as (x, 1). Other groups by
    // sincosf_glibc.
    const sincos_t* p0 = &sincosf_table[0];
    const __m256d sgnbit = _mm256_set1_pd(-0.0);
    const __m128i a120 = _mm_set1_epi32(static_cast<int>(abstop12(120.0f)));
    const __m128i atiny = _mm_set1_epi32(static_cast<int>(abstop12(0x1p-12f)));
    for (; i + 4 <= n; i += 4) {
        const __m128 y = _mm_loadu_ps(x + i);
        const __m128i a12 = _mm_and_si128(_mm_srli_epi32(_mm_castps_si128(y), 20), _mm_set1_epi32(0x7ff));
        if (_mm_movemask_ps(_mm_castsi128_ps(_mm_cmplt_epi32(a12, a120))) != 0xF) {
            for (int l = 0; l < 4; ++l) sincosf_glibc(x[i + l], sinp + i + l, cosp + i + l);
            continue;
        }
        const __m128 tiny = _mm_castsi128_ps(_mm_cmplt_epi32(a12, atiny));
        const __m256d xd = _mm256_cvtps_pd(y);
        // reduce_fast
        const __m256d r = _mm256_mul_pd(xd, _mm256_set1_pd(p0->hpi_inv));
        const __m128i nq = _mm_srai_epi32(_mm_add_epi32(_mm256_cvttpd_epi32(r), _mm_set1_epi32(0x800000)), 24);
        const __m256d nd = _mm256_cvtepi32_pd(nq);
        const __m256d xr = _mm256_fmadd_pd(_mm256_xor_pd(nd, sgnbit), _mm256_set1_pd(p0->hpi), xd);
        // sign[n & 3] = -1 for n & 3 = 1, 2; the second table (cosine
        // negated) for n & 2; the swap for n & 1.
        const __m128i n1 = _mm_add_epi32(nq, _mm_set1_epi32(1));
        const __m256d neg = _mm256_castsi256_pd(_mm256_cvtepi32_epi64(
            _mm_cmpeq_epi32(_mm_and_si128(n1, _mm_set1_epi32(2)), _mm_set1_epi32(2))));
        const __m256d t1 = _mm256_castsi256_pd(_mm256_cvtepi32_epi64(
            _mm_cmpeq_epi32(_mm_and_si128(nq, _mm_set1_epi32(2)), _mm_set1_epi32(2))));
        const __m256d odd = _mm256_castsi256_pd(_mm256_cvtepi32_epi64(
            _mm_cmpeq_epi32(_mm_and_si128(nq, _mm_set1_epi32(1)), _mm_set1_epi32(1))));
        const __m256d sg = _mm256_blendv_pd(_mm256_set1_pd(1.0), _mm256_set1_pd(-1.0), neg);
        auto cc = [&](double v) {  // a cosine coefficient of the lane's table
            return _mm256_blendv_pd(_mm256_set1_pd(v), _mm256_set1_pd(-v), t1);
        };
        const __m256d xs = _mm256_mul_pd(xr, sg);  // x * s
        const __m256d x2 = _mm256_mul_pd(xr, xr);  // x * x
        // sincosf_poly(xs, x2, p, n)
        const __m256d x4 = _mm256_mul_pd(x2, x2);
        const __m256d x3 = _mm256_mul_pd(x2, xs);
        const __m256d c2 = _mm256_fmadd_pd(x2, cc(p0->c4), cc(p0->c3));
        const __m256d s1 = _mm256_fmadd_pd(x2, _mm256_set1_pd(p0->s3), _mm256_set1_pd(p0->s2));
        const __m256d c1 = _mm256_fmadd_pd(x2, cc(p0->c1), cc(p0->c0));
        const __m256d x5 = _mm256_mul_pd(x3, x2);
        const __m256d x6 = _mm256_mul_pd(x4, x2);
        const __m256d sv = _mm256_fmadd_pd(x3, _mm256_set1_pd(p0->s1), xs);
        const __m256d cv = _mm256_fmadd_pd(x4, cc(p0->c2), c1);
        const __m256d so = _mm256_fmadd_pd(x5, s1, sv);
        const __m256d co = _mm256_fmadd_pd(x6, c2, cv);
        __m128 fs = _mm256_cvtpd_ps(_mm256_blendv_pd(so, co, odd));
        __m128 fc = _mm256_cvtpd_ps(_mm256_blendv_pd(co, so, odd));
        fs = _mm_blendv_ps(fs, y, tiny);
        fc = _mm_blendv_ps(fc, _mm_set1_ps(1.0f), tiny);
        _mm_storeu_ps(sinp + i, fs);
        _mm_storeu_ps(cosp + i, fc);
    }
#elif defined(FT8_MATHS_NEON)
    // As the AVX2 path: two double lanes a register, the four as lo
    // (lanes 0, 1) and hi (2, 3); fma(a, b, c) = vfmaq_f64(c, a, b).
    const sincos_t* p0 = &sincosf_table[0];
    const uint32x4_t a120 = vdupq_n_u32(abstop12(120.0f));
    const uint32x4_t atiny = vdupq_n_u32(abstop12(0x1p-12f));
    for (; i + 4 <= n; i += 4) {
        const float32x4_t y = vld1q_f32(x + i);
        const uint32x4_t a12 = vandq_u32(vshrq_n_u32(vreinterpretq_u32_f32(y), 20), vdupq_n_u32(0x7ff));
        if (vminvq_u32(vcltq_u32(a12, a120)) == 0) {
            for (int l = 0; l < 4; ++l) sincosf_glibc(x[i + l], sinp + i + l, cosp + i + l);
            continue;
        }
        const uint32x4_t tiny = vcltq_u32(a12, atiny);
        const float64x2_t xd[2] = {vcvt_f64_f32(vget_low_f32(y)), vcvt_high_f64_f32(y)};
        // reduce_fast: n = ((int32_t)r + 0x800000) >> 24 (|r| < 2^31 here)
        float64x2_t r[2];
        for (int h = 0; h < 2; ++h) r[h] = vmulq_f64(xd[h], vdupq_n_f64(p0->hpi_inv));
        const int32x4_t ri = vcombine_s32(vmovn_s64(vcvtq_s64_f64(r[0])), vmovn_s64(vcvtq_s64_f64(r[1])));
        const int32x4_t nq = vshrq_n_s32(vaddq_s32(ri, vdupq_n_s32(0x800000)), 24);
        const uint32x4_t neg4 = vceqq_s32(vandq_s32(vaddq_s32(nq, vdupq_n_s32(1)), vdupq_n_s32(2)), vdupq_n_s32(2));
        const uint32x4_t t14 = vceqq_s32(vandq_s32(nq, vdupq_n_s32(2)), vdupq_n_s32(2));
        const uint32x4_t odd4 = vceqq_s32(vandq_s32(nq, vdupq_n_s32(1)), vdupq_n_s32(1));
        float32x2_t fs[2], fc[2];
        for (int h = 0; h < 2; ++h) {
            auto wide = [&](uint32x4_t m) {  // lanes 2h, 2h+1 as 64-bit masks
                const int32x2_t m2 = h == 0 ? vget_low_s32(vreinterpretq_s32_u32(m)) : vget_high_s32(vreinterpretq_s32_u32(m));
                return vreinterpretq_u64_s64(vmovl_s32(m2));
            };
            const int32x2_t n2 = h == 0 ? vget_low_s32(nq) : vget_high_s32(nq);
            const float64x2_t nd = vcvtq_f64_s64(vmovl_s32(n2));
            const float64x2_t xr = vfmaq_f64(xd[h], vnegq_f64(nd), vdupq_n_f64(p0->hpi));
            const uint64x2_t neg = wide(neg4), t1 = wide(t14), odd = wide(odd4);
            const float64x2_t sg = vbslq_f64(neg, vdupq_n_f64(-1.0), vdupq_n_f64(1.0));
            auto cc = [&](double v) { return vbslq_f64(t1, vdupq_n_f64(-v), vdupq_n_f64(v)); };
            const float64x2_t xs = vmulq_f64(xr, sg);  // x * s
            const float64x2_t x2 = vmulq_f64(xr, xr);  // x * x
            // sincosf_poly(xs, x2, p, n)
            const float64x2_t x4 = vmulq_f64(x2, x2);
            const float64x2_t x3 = vmulq_f64(x2, xs);
            const float64x2_t c2 = vfmaq_f64(cc(p0->c3), x2, cc(p0->c4));
            const float64x2_t s1 = vfmaq_f64(vdupq_n_f64(p0->s2), x2, vdupq_n_f64(p0->s3));
            const float64x2_t c1 = vfmaq_f64(cc(p0->c0), x2, cc(p0->c1));
            const float64x2_t x5 = vmulq_f64(x3, x2);
            const float64x2_t x6 = vmulq_f64(x4, x2);
            const float64x2_t sv = vfmaq_f64(xs, x3, vdupq_n_f64(p0->s1));
            const float64x2_t cv = vfmaq_f64(c1, x4, cc(p0->c2));
            const float64x2_t so = vfmaq_f64(sv, x5, s1);
            const float64x2_t co = vfmaq_f64(cv, x6, c2);
            fs[h] = vcvt_f32_f64(vbslq_f64(odd, co, so));
            fc[h] = vcvt_f32_f64(vbslq_f64(odd, so, co));
        }
        float32x4_t s4 = vcombine_f32(fs[0], fs[1]), c4 = vcombine_f32(fc[0], fc[1]);
        s4 = vbslq_f32(tiny, y, s4);
        c4 = vbslq_f32(tiny, vdupq_n_f32(1.0f), c4);
        vst1q_f32(sinp + i, s4);
        vst1q_f32(cosp + i, c4);
    }
#endif
    for (; i < n; ++i) sincosf_glibc(x[i], sinp + i, cosp + i);
}

float expf_glibc(float x) {
    const double xd = static_cast<double>(x);
    const uint32_t abstop = top12(x) & 0x7ff;
    if (abstop >= top12(88.0f)) {
        // |x| >= 88 or x is nan.
        if (fbits(x) == fbits(-INFINITY)) return 0.0f;
        if (abstop >= top12(INFINITY)) return x + x;
        if (x > 0x1.62e42ep6f) return math_oflowf(0);   // x > log(0x1p128) ~= 88.72
        if (x < -0x1.9fe368p6f) return math_uflowf(0);  // x < log(0x1p-150) ~= -103.97
        if (x < -0x1.9d1d9ep6f) return math_may_uflowf(0);  // x < log(0x1p-149) ~= -103.28
    }
    // x*N/Ln2 = k + r with r in [-1/2, 1/2] and int k. The FMA build fuses
    // z = InvLn2N * xd into both its uses.
    double kd = std::fma(exp2f_invln2_scaled, xd, exp2f_shift);
    const uint64_t ki = asuint64(kd);
    kd -= exp2f_shift;
    const double r = std::fma(exp2f_invln2_scaled, xd, -kd);
    // exp(x) = 2^(k/N) * 2^(r/N) ~= s * (C0*r^3 + C1*r^2 + C2*r + 1)
    const double* C = exp2f_poly_scaled;
    uint64_t t = exp2f_tab[ki % EXP2F_N];
    t += ki << (52 - EXP2F_TABLE_BITS);
    const double s = asdouble(t);
    const double z = std::fma(C[0], r, C[1]);
    const double r2 = r * r;
    double y = std::fma(C[2], r, 1.0);
    y = std::fma(z, r2, y);
    y = y * s;
    return static_cast<float>(y);
}

float logf_glibc(float x) {
    constexpr uint32_t OFF = 0x3f330000;
    constexpr int N = 1 << LOGF_TABLE_BITS;
    uint32_t ix = fbits(x);
    // Fix sign of zero with downward rounding when x==1.
    if (ix == 0x3f800000) return 0;
    if (ix - 0x00800000 >= 0x7f800000 - 0x00800000) {
        // x < 0x1p-126 or inf or nan.
        if (ix * 2 == 0) return math_divzerof(1);
        if (ix == 0x7f800000) return x;  // log(inf) == inf.
        if ((ix & 0x80000000) || ix * 2 >= 0xff000000) return math_invalidf(x);
        // x is subnormal, normalize it.
        ix = fbits(x * 0x1p23f);
        ix -= 23 << 23;
    }
    // x = 2^k z; where z is in range [OFF,2*OFF] and exact.
    const uint32_t tmp = ix - OFF;
    const int i = static_cast<int>((tmp >> (23 - LOGF_TABLE_BITS)) % N);
    const int k = static_cast<int32_t>(tmp) >> 23;  // arithmetic shift
    const uint32_t iz = ix - (tmp & 0x1ffu << 23);
    const double invc = logf_tab[i].invc;
    const double logc = logf_tab[i].logc;
    const double z = static_cast<double>(bitsf(iz));
    // log(x) = log1p(z/c-1) + log(c) + k*Ln2
    const double r = z * invc - 1;
    const double y0 = logc + static_cast<double>(k) * logf_ln2;
    const double* A = logf_poly;
    const double r2 = r * r;
    double y = A[1] * r + A[2];
    y = A[0] * r2 + y;
    y = y * r2 + (y0 + r);
    return static_cast<float>(y);
}

float powf_glibc(float x, float y) {
    constexpr uint32_t SIGN_BIAS = 1u << (EXP2F_TABLE_BITS + 11);
    uint32_t sign_bias = 0;
    uint32_t ix = fbits(x);
    const uint32_t iy = fbits(y);
    if (ix - 0x00800000 >= 0x7f800000 - 0x00800000 || zeroinfnan(iy)) {
        // Either (x < 0x1p-126 or inf or nan) or (y is 0 or inf or nan).
        if (zeroinfnan(iy)) {
            if (2 * iy == 0) return issignalingf(x) ? x + y : 1.0f;
            if (ix == 0x3f800000) return issignalingf(y) ? x + y : 1.0f;
            if (2 * ix > 2u * 0x7f800000 || 2 * iy > 2u * 0x7f800000) return x + y;
            if (2 * ix == 2 * 0x3f800000) return 1.0f;
            if ((2 * ix < 2 * 0x3f800000) == !(iy & 0x80000000)) return 0.0f;  // |x|<1 && y==inf or |x|>1 && y==-inf.
            return y * y;
        }
        if (zeroinfnan(ix)) {
            float x2 = x * x;
            if (ix & 0x80000000 && powf_checkint(iy) == 1) {
                x2 = -x2;
                sign_bias = 1;
            }
            if (2 * ix == 0 && iy & 0x80000000) return math_divzerof(sign_bias);
            return iy & 0x80000000 ? 1 / x2 : x2;
        }
        // x and y are non-zero finite.
        if (ix & 0x80000000) {
            // Finite x < 0.
            const int yint = powf_checkint(iy);
            if (yint == 0) return math_invalidf(x);
            if (yint == 1) sign_bias = SIGN_BIAS;
            ix &= 0x7fffffff;
        }
        if (ix < 0x00800000) {
            // Normalize subnormal x so exponent becomes negative.
            ix = fbits(x * 0x1p23f);
            ix &= 0x7fffffff;
            ix -= 23 << 23;
        }
    }
    const double logx = powf_log2_inline(ix);
    const double ylogx = y * logx;  // Note: cannot overflow, y is single prec.
    if ((asuint64(ylogx) >> 47 & 0xffff) >= asuint64(126.0) >> 47) {
        // |y*log(x)| >= 126.
        if (ylogx > 0x1.fffffffd1d571p+6) return math_oflowf(sign_bias);  // |x^y| > 0x1.ffffffp127.
        if (ylogx > 0x1.fffffffa3aae2p+6) {
            // |x^y| > 0x1.fffffep127, check if we round away from 0.
            volatile float tiny = 0x1p-25f;
            if ((!sign_bias && 1.0f + tiny != 1.0f) || (sign_bias && -1.0f - tiny != -1.0f))
                return math_oflowf(sign_bias);
        }
        if (ylogx <= -150.0) return math_uflowf(sign_bias);
        if (ylogx < -149.0) return math_may_uflowf(sign_bias);
    }
    return static_cast<float>(powf_exp2_inline(ylogx, sign_bias));
}

void powf_half_glibc(const float* x, float* out, int n) {
    int i = 0;
#if defined(FT8_MATHS_AVX2)
    // Four at a time in double lanes: powf_glibc's own path for a positive
    // normal x and y = 0.5 (sign_bias 0; |0.5*log2(x)| < 64, so no
    // overflow test applies), each operation as there, none fused; +0
    // gives +0 as there. Groups with any other x go through powf_glibc.
    const double* A = powf_log2_poly;
    const double* C = exp2f_poly;
    for (; i + 4 <= n; i += 4) {
        const __m128i ix = _mm_loadu_si128(reinterpret_cast<const __m128i*>(x + i));
        const __m128i tn = _mm_sub_epi32(ix, _mm_set1_epi32(0x00800000));
        const __m128i ok = _mm_and_si128(_mm_cmpgt_epi32(tn, _mm_set1_epi32(-1)),
                                         _mm_cmplt_epi32(tn, _mm_set1_epi32(0x7f800000 - 0x00800000)));
        // x = +0: powf gives x*x = +0 (y positive), put in below.
        const __m128i zero = _mm_cmpeq_epi32(ix, _mm_setzero_si128());
        if (_mm_movemask_ps(_mm_castsi128_ps(_mm_or_si128(ok, zero))) != 0xF) {
            for (int l = 0; l < 4; ++l) out[i + l] = powf_glibc(x[i + l], 0.5f);
            continue;
        }
        // powf_log2_inline(ix)
        const __m128i tmp = _mm_sub_epi32(ix, _mm_set1_epi32(0x3f330000));
        const __m128i idx = _mm_and_si128(_mm_srli_epi32(tmp, 23 - POWF_LOG2_TABLE_BITS), _mm_set1_epi32(15));
        const __m128i top = _mm_and_si128(tmp, _mm_set1_epi32(static_cast<int>(0xff800000u)));
        const __m128i iz = _mm_sub_epi32(ix, top);
        const __m128i k = _mm_srai_epi32(top, 23);
        const __m128i i2 = _mm_slli_epi32(idx, 1);  // {invc, logc}: two doubles an entry
        const __m256d invc = _mm256_i32gather_pd(&powf_log2_tab[0].invc, i2, 8);
        const __m256d logc = _mm256_i32gather_pd(&powf_log2_tab[0].logc, i2, 8);
        const __m256d z = _mm256_cvtps_pd(_mm_castsi128_ps(iz));
        const __m256d r = _mm256_sub_pd(_mm256_mul_pd(z, invc), _mm256_set1_pd(1.0));
        const __m256d y0 = _mm256_add_pd(logc, _mm256_cvtepi32_pd(k));
        const __m256d r2 = _mm256_mul_pd(r, r);
        __m256d y = _mm256_add_pd(_mm256_mul_pd(_mm256_set1_pd(A[0]), r), _mm256_set1_pd(A[1]));
        const __m256d p = _mm256_add_pd(_mm256_mul_pd(_mm256_set1_pd(A[2]), r), _mm256_set1_pd(A[3]));
        const __m256d r4 = _mm256_mul_pd(r2, r2);
        __m256d q = _mm256_add_pd(_mm256_mul_pd(_mm256_set1_pd(A[4]), r), y0);
        q = _mm256_add_pd(_mm256_mul_pd(p, r2), q);
        y = _mm256_add_pd(_mm256_mul_pd(y, r4), q);
        const __m256d xd = _mm256_mul_pd(_mm256_set1_pd(0.5), y);  // ylogx = y * logx
        // powf_exp2_inline(ylogx, 0)
        __m256d kd = _mm256_add_pd(xd, _mm256_set1_pd(exp2f_shift_scaled));
        const __m256i ki = _mm256_castpd_si256(kd);
        kd = _mm256_sub_pd(kd, _mm256_set1_pd(exp2f_shift_scaled));
        const __m256d rr = _mm256_sub_pd(xd, kd);
        __m256i t = _mm256_i64gather_epi64(reinterpret_cast<const long long*>(exp2f_tab),
                                           _mm256_and_si256(ki, _mm256_set1_epi64x(EXP2F_N - 1)), 8);
        t = _mm256_add_epi64(t, _mm256_slli_epi64(ki, 52 - EXP2F_TABLE_BITS));
        const __m256d sc = _mm256_castsi256_pd(t);
        const __m256d zz = _mm256_add_pd(_mm256_mul_pd(_mm256_set1_pd(C[0]), rr), _mm256_set1_pd(C[1]));
        const __m256d rr2 = _mm256_mul_pd(rr, rr);
        __m256d yy = _mm256_add_pd(_mm256_mul_pd(_mm256_set1_pd(C[2]), rr), _mm256_set1_pd(1.0));
        yy = _mm256_add_pd(_mm256_mul_pd(zz, rr2), yy);
        yy = _mm256_mul_pd(yy, sc);
        _mm_storeu_ps(out + i, _mm_andnot_ps(_mm_castsi128_ps(zero), _mm256_cvtpd_ps(yy)));
    }
#elif defined(FT8_MATHS_NEON)
    // As the AVX2 path, two double lanes a register; the table entries
    // read lane by lane.
    const double* A = powf_log2_poly;
    const double* C = exp2f_poly;
    for (; i + 4 <= n; i += 4) {
        const uint32x4_t ix = vreinterpretq_u32_f32(vld1q_f32(x + i));
        // positive normal: ix - 0x00800000 < 0x7f800000 - 0x00800000, unsigned
        const uint32x4_t ok = vcltq_u32(vsubq_u32(ix, vdupq_n_u32(0x00800000)), vdupq_n_u32(0x7f800000 - 0x00800000));
        // x = +0: powf gives x*x = +0 (y positive), put in below.
        const uint32x4_t zero = vceqq_u32(ix, vdupq_n_u32(0));
        if (vminvq_u32(vorrq_u32(ok, zero)) == 0) {
            for (int l = 0; l < 4; ++l) out[i + l] = powf_glibc(x[i + l], 0.5f);
            continue;
        }
        // powf_log2_inline(ix)
        const uint32x4_t tmp = vsubq_u32(ix, vdupq_n_u32(0x3f330000));
        alignas(16) uint32_t idx[4];
        vst1q_u32(idx, vandq_u32(vshrq_n_u32(tmp, 23 - POWF_LOG2_TABLE_BITS), vdupq_n_u32(15)));
        const uint32x4_t top = vandq_u32(tmp, vdupq_n_u32(0xff800000u));
        const float32x4_t izf = vreinterpretq_f32_u32(vsubq_u32(ix, top));
        const int32x4_t k = vshrq_n_s32(vreinterpretq_s32_u32(top), 23);
        float32x2_t res[2];
        for (int h = 0; h < 2; ++h) {
            const int a = idx[2 * h], b = idx[2 * h + 1];
            const double ic[2] = {powf_log2_tab[a].invc, powf_log2_tab[b].invc};
            const double lc[2] = {powf_log2_tab[a].logc, powf_log2_tab[b].logc};
            const float64x2_t invc = vld1q_f64(ic), logc = vld1q_f64(lc);
            const float64x2_t z = h == 0 ? vcvt_f64_f32(vget_low_f32(izf)) : vcvt_high_f64_f32(izf);
            const float64x2_t kd0 = vcvtq_f64_s64(vmovl_s32(h == 0 ? vget_low_s32(k) : vget_high_s32(k)));
            const float64x2_t r = vsubq_f64(vmulq_f64(z, invc), vdupq_n_f64(1.0));
            const float64x2_t y0 = vaddq_f64(logc, kd0);
            const float64x2_t r2 = vmulq_f64(r, r);
            float64x2_t y = vaddq_f64(vmulq_f64(vdupq_n_f64(A[0]), r), vdupq_n_f64(A[1]));
            const float64x2_t p = vaddq_f64(vmulq_f64(vdupq_n_f64(A[2]), r), vdupq_n_f64(A[3]));
            const float64x2_t r4 = vmulq_f64(r2, r2);
            float64x2_t q = vaddq_f64(vmulq_f64(vdupq_n_f64(A[4]), r), y0);
            q = vaddq_f64(vmulq_f64(p, r2), q);
            y = vaddq_f64(vmulq_f64(y, r4), q);
            const float64x2_t xd = vmulq_f64(vdupq_n_f64(0.5), y);  // ylogx = y * logx
            // powf_exp2_inline(ylogx, 0)
            float64x2_t kd = vaddq_f64(xd, vdupq_n_f64(exp2f_shift_scaled));
            const uint64x2_t ki = vreinterpretq_u64_f64(kd);
            kd = vsubq_f64(kd, vdupq_n_f64(exp2f_shift_scaled));
            const float64x2_t rr = vsubq_f64(xd, kd);
            const uint64_t tb[2] = {exp2f_tab[vgetq_lane_u64(ki, 0) % EXP2F_N], exp2f_tab[vgetq_lane_u64(ki, 1) % EXP2F_N]};
            const uint64x2_t t = vaddq_u64(vld1q_u64(tb), vshlq_n_u64(ki, 52 - EXP2F_TABLE_BITS));
            const float64x2_t sc = vreinterpretq_f64_u64(t);
            const float64x2_t zz = vaddq_f64(vmulq_f64(vdupq_n_f64(C[0]), rr), vdupq_n_f64(C[1]));
            const float64x2_t rr2 = vmulq_f64(rr, rr);
            float64x2_t yy = vaddq_f64(vmulq_f64(vdupq_n_f64(C[2]), rr), vdupq_n_f64(1.0));
            yy = vaddq_f64(vmulq_f64(zz, rr2), yy);
            yy = vmulq_f64(yy, sc);
            res[h] = vcvt_f32_f64(yy);
        }
        const uint32x4_t bits = vreinterpretq_u32_f32(vcombine_f32(res[0], res[1]));
        vst1q_f32(out + i, vreinterpretq_f32_u32(vbicq_u32(bits, zero)));
    }
#endif
    for (; i < n; ++i) out[i] = powf_glibc(x[i], 0.5f);
}

void hypotf_n_glibc(const float* x, const float* y, float* out, int n) {
    // aarch64: hypotf_glibc a value (a NEON copy gained nothing
    // measurable on a Pi 5, under 0.2% of the FT4 decode).
    int i = 0;
#if defined(FT8_MATHS_AVX2)
    // Four at a time in double lanes when all eight inputs are finite:
    // the exact squares added once and the square root, both correctly
    // rounded, as hypotf_glibc's finite path; other groups by it.
    const __m128 inf = _mm_set1_ps(INFINITY);
    const __m128 absm = _mm_castsi128_ps(_mm_set1_epi32(0x7fffffff));
    for (; i + 4 <= n; i += 4) {
        const __m128 vx = _mm_loadu_ps(x + i), vy = _mm_loadu_ps(y + i);
        const __m128 fin = _mm_and_ps(_mm_cmplt_ps(_mm_and_ps(vx, absm), inf), _mm_cmplt_ps(_mm_and_ps(vy, absm), inf));
        if (_mm_movemask_ps(fin) != 0xF) {
            for (int l = 0; l < 4; ++l) out[i + l] = hypotf_glibc(x[i + l], y[i + l]);
            continue;
        }
        const __m256d dx = _mm256_cvtps_pd(vx), dy = _mm256_cvtps_pd(vy);
        const __m256d ss = _mm256_add_pd(_mm256_mul_pd(dx, dx), _mm256_mul_pd(dy, dy));
        _mm_storeu_ps(out + i, _mm256_cvtpd_ps(_mm256_sqrt_pd(ss)));
    }
#endif
    for (; i < n; ++i) out[i] = hypotf_glibc(x[i], y[i]);
}

void hypotf_c_glibc(const float* xy, float* out, int n) {
    int i = 0;
#if defined(FT8_MATHS_AVX2)
    // As hypotf_n_glibc, from interleaved pairs: each pair's two exact
    // squares added (horizontally, x*x + y*y) and the square root.
    const __m256 inf = _mm256_set1_ps(INFINITY);
    const __m256 absm = _mm256_castsi256_ps(_mm256_set1_epi32(0x7fffffff));
    for (; i + 4 <= n; i += 4) {
        const __m256 v = _mm256_loadu_ps(xy + 2 * i);
        if (_mm256_movemask_ps(_mm256_cmp_ps(_mm256_and_ps(v, absm), inf, _CMP_LT_OQ)) != 0xFF) {
            for (int l = 0; l < 4; ++l) out[i + l] = hypotf_glibc(xy[2 * (i + l)], xy[2 * (i + l) + 1]);
            continue;
        }
        const __m256d lo = _mm256_cvtps_pd(_mm256_castps256_ps128(v));    // x0 y0 x1 y1
        const __m256d hi = _mm256_cvtps_pd(_mm256_extractf128_ps(v, 1));  // x2 y2 x3 y3
        const __m256d ss = _mm256_hadd_pd(_mm256_mul_pd(lo, lo), _mm256_mul_pd(hi, hi));  // 0 2 1 3
        const __m256d r = _mm256_sqrt_pd(_mm256_permute4x64_pd(ss, 0xD8));
        _mm_storeu_ps(out + i, _mm256_cvtpd_ps(r));
    }
#endif
    for (; i < n; ++i) out[i] = hypotf_glibc(xy[2 * i], xy[2 * i + 1]);
}

float hypotf_glibc(float x, float y) {
    if (!std::isfinite(x) || !std::isfinite(y)) {
        if ((std::isinf(x) || std::isinf(y)) && !issignalingf(x) && !issignalingf(y)) return INFINITY;
        return x + y;
    }
    return static_cast<float>(std::sqrt(static_cast<double>(x) * static_cast<double>(x) +
                                        static_cast<double>(y) * static_cast<double>(y)));
}


float erff_glibc(float x) {
    const float tiny = 1e-30f, one = 1.0f, erx = bitsf(0x3f58560b);
    // Coefficients for approximation to  erf on [0,0.84375]
    const float efx = bitsf(0x3e0375d4), pp0 = bitsf(0x3e0375d4), pp1 = bitsf(0xbea66beb), pp2 = bitsf(0xbce9528f),
                pp3 = bitsf(0xbbbd1489), pp4 = bitsf(0xb7c756b1), qq1 = bitsf(0x3ecbbbce), qq2 = bitsf(0x3d852a63),
                qq3 = bitsf(0x3ba68116), qq4 = bitsf(0x390aee49), qq5 = bitsf(0xb684e21a);
    // Coefficients for approximation to  erf  in [0.84375,1.25]
    const float pa0 = bitsf(0xbb1acdc6), pa1 = bitsf(0x3ed46805), pa2 = bitsf(0xbebe9208), pa3 = bitsf(0x3ea2fe54),
                pa4 = bitsf(0xbde31cc2), pa5 = bitsf(0x3d1151b3), pa6 = bitsf(0xbb0df9c0), qa1 = bitsf(0x3dd9f331),
                qa2 = bitsf(0x3f0a5785), qa3 = bitsf(0x3d931ae7), qa4 = bitsf(0x3e013307), qa5 = bitsf(0x3c5f6e13),
                qa6 = bitsf(0x3c445aa3);
    // Coefficients for approximation to  erfc in [1.25,1/0.35]
    const float ra0 = bitsf(0xbc21a093), ra1 = bitsf(0xbf31a0b7), ra2 = bitsf(0xc128f022), ra3 = bitsf(0xc2798057),
                ra4 = bitsf(0xc322658c), ra5 = bitsf(0xc3389ae7), ra6 = bitsf(0xc2a2932b), ra7 = bitsf(0xc11d077e),
                sa1 = bitsf(0x419d35ce), sa2 = bitsf(0x4309a863), sa3 = bitsf(0x43d9486f), sa4 = bitsf(0x442158c9),
                sa5 = bitsf(0x43d6810b), sa6 = bitsf(0x42d9451f), sa7 = bitsf(0x40d23f7c), sa8 = bitsf(0xbd777f97);
    // Coefficients for approximation to  erfc in [1/.35,28]
    const float rb0 = bitsf(0xbc21a092), rb1 = bitsf(0xbf4c9dd4), rb2 = bitsf(0xc18e104b), rb3 = bitsf(0xc320a2ea),
                rb4 = bitsf(0xc41f6441), rb5 = bitsf(0xc480230b), rb6 = bitsf(0xc3f1c275), sb1 = bitsf(0x41f2b459),
                sb2 = bitsf(0x43a2e571), sb3 = bitsf(0x44c01759), sb4 = bitsf(0x4547fdbb), sb5 = bitsf(0x451f90ce),
                sb6 = bitsf(0x43ed43a7), sb7 = bitsf(0xc1b38712);
    const int32_t hx = static_cast<int32_t>(fbits(x));
    int32_t ix = hx & 0x7fffffff;
    if (ix >= 0x7f800000) {  // erf(nan)=nan
        const int32_t i = static_cast<int32_t>((static_cast<uint32_t>(hx) >> 31) << 1);
        return static_cast<float>(1 - i) + one / x;  // erf(+-inf)=+-1
    }
    if (ix < 0x3f580000) {      // |x|<0.84375
        if (ix < 0x31800000) {  // |x|<2**-28
            if (ix < 0x04000000) return 0.0625f * (16.0f * x + (16.0f * efx) * x);  // avoid spurious underflow
            return x + efx * x;
        }
        const float z = x * x;
        const float r = pp0 + z * (pp1 + z * (pp2 + z * (pp3 + z * pp4)));
        const float s = one + z * (qq1 + z * (qq2 + z * (qq3 + z * (qq4 + z * qq5))));
        const float y = r / s;
        return x + x * y;
    }
    if (ix < 0x3fa00000) {  // 0.84375 <= |x| < 1.25
        const float s = std::fabs(x) - one;
        const float P = pa0 + s * (pa1 + s * (pa2 + s * (pa3 + s * (pa4 + s * (pa5 + s * pa6)))));
        const float Q = one + s * (qa1 + s * (qa2 + s * (qa3 + s * (qa4 + s * (qa5 + s * qa6)))));
        if (hx >= 0) return erx + P / Q;
        return -erx - P / Q;
    }
    if (ix >= 0x40c00000) {  // inf>|x|>=6
        if (hx >= 0) return one - tiny;
        return tiny - one;
    }
    x = std::fabs(x);
    const float s = one / (x * x);
    float R, S;
    if (ix < 0x4036DB6E) {  // |x| < 1/0.35
        R = ra0 + s * (ra1 + s * (ra2 + s * (ra3 + s * (ra4 + s * (ra5 + s * (ra6 + s * ra7))))));
        S = one + s * (sa1 + s * (sa2 + s * (sa3 + s * (sa4 + s * (sa5 + s * (sa6 + s * (sa7 + s * sa8)))))));
    } else {  // |x| >= 1/0.35
        R = rb0 + s * (rb1 + s * (rb2 + s * (rb3 + s * (rb4 + s * (rb5 + s * rb6)))));
        S = one + s * (sb1 + s * (sb2 + s * (sb3 + s * (sb4 + s * (sb5 + s * (sb6 + s * sb7))))));
    }
    ix = static_cast<int32_t>(fbits(x));
    const float z = bitsf(static_cast<uint32_t>(ix) & 0xfffff000u);
    const float r = expf_glibc(-z * z - 0.5625f) * expf_glibc((z - x) * (z + x) + R / S);
    if (hx >= 0) return one - r / x;
    return r / x - one;
}

float log10f_glibc(float x) {
    const float two25 = bitsf(0x4c000000), ivln10 = bitsf(0x3ede5bd9), log10_2hi = bitsf(0x3e9a2080),
                log10_2lo = bitsf(0x355427db);
    int32_t hx = static_cast<int32_t>(fbits(x));
    int32_t k = 0;
    if (hx < 0x00800000) {  // x < 2**-126
        if ((hx & 0x7fffffff) == 0) return -two25 / std::fabs(x);  // log(+-0)=-inf
        if (hx < 0) return (x - x) / (x - x);                       // log(-#) = NaN
        k -= 25;
        x *= two25;  // subnormal number, scale up x
        hx = static_cast<int32_t>(fbits(x));
    }
    if (hx >= 0x7f800000) return x + x;
    k += (hx >> 23) - 127;
    const int32_t i = static_cast<int32_t>((static_cast<uint32_t>(k) & 0x80000000u) >> 31);
    hx = (hx & 0x007fffff) | ((0x7f - i) << 23);
    const float y = static_cast<float>(k + i);
    x = bitsf(static_cast<uint32_t>(hx));
    const float z = y * log10_2lo + ivln10 * logf_glibc(x);
    return z + y * log10_2hi;
}

void cosf4_glibc(const float* in, float* out) {
    for (int l = 0; l < 4; ++l) {
        const float x = in[l];
        if (!(std::fabs(x) <= sAbsRangeReductionVal)) {  // large and special arguments
            out[l] = cosf_glibc(x);
            continue;
        }
        // Add Pi/2 to X (Cos to Sin), the octant Y by 1/Pi and the Right
        // Shifter; its last bit is the sign S; Y less the Right Shifter and
        // 0.5; X less Y*PI in four parts.
        float y = sHalfPI + x;
        y = y * sInvPI;
        y = y + sRShifter;
        const uint32_t s = fbits(y) << 31;
        y = y - sRShifter;
        y = y - sOneHalf;
        float r = x - sPI1 * y;
        r = r - sPI2 * y;
        r = r - sPI3 * y;
        r = r - y * sPI4;
        out[l] = trig_poly(r, s);
    }
}

void sinf4_glibc(const float* in, float* out) {
    for (int l = 0; l < 4; ++l) {
        const float x = in[l];
        const uint32_t sign = fbits(x) & 0x80000000u;
        const float ax = bitsf(fbits(x) & 0x7fffffffu);
        if (!(ax <= sAbsRangeReductionVal)) {  // large and special values
            out[l] = sinf_glibc(x);
            continue;
        }
        // The octant Y of |X| by 1/Pi and the Right Shifter; its last bit is
        // the sign S; Y less the Right Shifter; |X| less Y*PI in four parts;
        // the source's sign on the result.
        float y = sInvPI * ax;
        y = y + sRShifter;
        const uint32_t s = fbits(y) << 31;
        y = y - sRShifter;
        float r = ax - sPI1 * y;
        r = r - sPI2 * y;
        r = r - sPI3 * y;
        r = r - y * sPI4;
        out[l] = bitsf(fbits(trig_poly(r, s)) ^ sign);
    }
}

void log10f4_glibc(const float* in, float* out) {
    // __svml_slog10_data_internal
    const float MinNorm = bitsf(0x00800000), MaxNorm = bitsf(0x7f7fffff), L2H = bitsf(0x3e9a2100),
                L2L = bitsf(0xb64AF600), One = bitsf(0x3f800000);
    const uint32_t iBrkValue = 0x3f2aaaab, iOffExpoMask = 0x007fffff;
    const float c9 = bitsf(0x3d8063B4), c8 = bitsf(0xbd890073), c7 = bitsf(0x3d775317), c6 = bitsf(0xbd91FB27),
                c5 = bitsf(0x3dB20B96), c4 = bitsf(0xbdDE6E20), c3 = bitsf(0x3e143CE5), c2 = bitsf(0xbe5E5BC5),
                c1 = bitsf(0x3eDE5BD9);
    for (int l = 0; l < 4; ++l) {
        const float x = in[l];
        if (x < MinNorm || !(x <= MaxNorm)) {  // special inputs
            out[l] = log10f_glibc(x);
            continue;
        }
        // reduction: compute r, n
        const uint32_t t = fbits(x) - iBrkValue;
        const float m = bitsf((t & iOffExpoMask) + iBrkValue);
        const float n = static_cast<float>(static_cast<int32_t>(t) >> 23);
        const float R = m - One;
        const float R2 = R * R;
        float p5 = c9 * R;
        float p6 = c7 * R;
        p5 = p5 + c8;
        float p7 = c5 * R;
        p6 = p6 + c6;
        float p9 = c3 * R;
        p5 = p5 * R2;
        p7 = p7 + c4;
        p9 = p9 + c2;
        p6 = p6 + p5;
        p6 = p6 * R2;
        p7 = p7 + p6;
        float res = L2L * n;
        const float q = R2 * p7;
        p9 = p9 + q;
        const float nh = n * L2H;
        p9 = p9 * R;
        p9 = p9 + c1;
        res = res + R * p9;
        out[l] = res + nh;
    }
}


namespace {

// svml_s_powf_data.S: _Log_LA_table (ln of the reciprocal, for each
// rounded reciprocal 512..1024 of the mantissa at 2^-10) and __dbT
// (2^(j/2048), j = 0..2048).
const uint64_t powf4_log_tab[513] = {
    0x8000000000000000, 0xbf5ff802a9ab10e6, 0xbf6ff00aa2b10bc0, 0xbf77ee11ebd82e94, 0xbf7fe02a6b106789,
    0xbf83e7295d25a7d9, 0xbf87dc475f810a77, 0xbf8bcf712c74384c, 0xbf8fc0a8b0fc03e4, 0xbf91d7f7eb9eebe7,
    0xbf93cea44346a575, 0xbf95c45a51b8d389, 0xbf97b91b07d5b11b, 0xbf99ace7551cc514, 0xbf9b9fc027af9198,
    0xbf9d91a66c543cc4, 0xbf9f829b0e783300, 0xbfa0b94f7c196176, 0xbfa1b0d98923d980, 0xbfa2a7ec2214e873,
    0xbfa39e87b9febd60, 0xbfa494acc34d911c, 0xbfa58a5bafc8e4d5, 0xbfa67f94f094bd98, 0xbfa77458f632dcfc,
    0xbfa868a83083f6cf, 0xbfa95c830ec8e3eb, 0xbfaa4fe9ffa3d235, 0xbfab42dd711971bf, 0xbfac355dd0921f2d,
    0xbfad276b8adb0b52, 0xbfae19070c276016, 0xbfaf0a30c01162a6, 0xbfaffae9119b9303, 0xbfb075983598e471,
    0xbfb0ed839b5526fe, 0xbfb16536eea37ae1, 0xbfb1dcb263db1944, 0xbfb253f62f0a1417, 0xbfb2cb0283f5de1f,
    0xbfb341d7961bd1d1, 0xbfb3b87598b1b6ee, 0xbfb42edcbea646f0, 0xbfb4a50d3aa1b040, 0xbfb51b073f06183f,
    0xbfb590cafdf01c28, 0xbfb60658a93750c4, 0xbfb67bb0726ec0fc, 0xbfb6f0d28ae56b4c, 0xbfb765bf23a6be13,
    0xbfb7da766d7b12cd, 0xbfb84ef898e8282a, 0xbfb8c345d6319b21, 0xbfb9375e55595ede, 0xbfb9ab42462033ad,
    0xbfba1ef1d8061cd4, 0xbfba926d3a4ad563, 0xbfbb05b49bee43fe, 0xbfbb78c82bb0eda1, 0xbfbbeba818146765,
    0xbfbc5e548f5bc743, 0xbfbcd0cdbf8c13e1, 0xbfbd4313d66cb35d, 0xbfbdb5270187d927, 0xbfbe27076e2af2e6,
    0xbfbe98b549671467, 0xbfbf0a30c01162a6, 0xbfbf7b79fec37ddf, 0xbfbfec9131dbeabb, 0xbfc02ebb42bf3d4b,
    0xbfc0671512ca596e, 0xbfc09f561ee719c3, 0xbfc0d77e7cd08e59, 0xbfc10f8e422539b1, 0xbfc14785846742ac,
    0xbfc17f6458fca611, 0xbfc1b72ad52f67a0, 0xbfc1eed90e2dc2c3, 0xbfc2266f190a5acb, 0xbfc25ded0abc6ad2,
    0xbfc29552f81ff523, 0xbfc2cca0f5f5f251, 0xbfc303d718e47fd3, 0xbfc33af575770e4f, 0xbfc371fc201e8f74,
    0xbfc3a8eb2d31a376, 0xbfc3dfc2b0ecc62a, 0xbfc41682bf727bc0, 0xbfc44d2b6ccb7d1e, 0xbfc483bccce6e3dd,
    0xbfc4ba36f39a55e5, 0xbfc4f099f4a230b2, 0xbfc526e5e3a1b438, 0xbfc55d1ad4232d6f, 0xbfc59338d9982086,
    0xbfc5c940075972b9, 0xbfc5ff3070a793d4, 0xbfc6350a28aaa758, 0xbfc66acd4272ad51, 0xbfc6a079d0f7aad2,
    0xbfc6d60fe719d21d, 0xbfc70b8f97a1aa75, 0xbfc740f8f54037a5, 0xbfc7764c128f2127, 0xbfc7ab890210d909,
    0xbfc7e0afd630c274, 0xbfc815c0a14357eb, 0xbfc84abb75865139, 0xbfc87fa06520c911, 0xbfc8b46f8223625b,
    0xbfc8e928de886d41, 0xbfc91dcc8c340bde, 0xbfc9525a9cf456b4, 0xbfc986d3228180ca, 0xbfc9bb362e7dfb83,
    0xbfc9ef83d2769a34, 0xbfca23bc1fe2b563, 0xbfca57df28244dcd, 0xbfca8becfc882f19, 0xbfcabfe5ae46124c,
    0xbfcaf3c94e80bff3, 0xbfcb2797ee46320c, 0xbfcb5b519e8fb5a4, 0xbfcb8ef670420c3b, 0xbfcbc286742d8cd6,
    0xbfcbf601bb0e44e2, 0xbfcc2968558c18c1, 0xbfcc5cba543ae425, 0xbfcc8ff7c79a9a22, 0xbfccc320c0176502,
    0xbfccf6354e09c5dc, 0xbfcd293581b6b3e7, 0xbfcd5c216b4fbb91, 0xbfcd8ef91af31d5e, 0xbfcdc1bca0abec7d,
    0xbfcdf46c0c722d2f, 0xbfce27076e2af2e6, 0xbfce598ed5a87e2f, 0xbfce8c0252aa5a60, 0xbfcebe61f4dd7b0b,
    0xbfcef0adcbdc5936, 0xbfcf22e5e72f105d, 0xbfcf550a564b7b37, 0xbfcf871b28955045, 0xbfcfb9186d5e3e2b,
    0xbfcfeb0233e607cc, 0xbfd00e6c45ad501d, 0xbfd0274dc16c232f, 0xbfd0402594b4d041, 0xbfd058f3c703ebc6,
    0xbfd071b85fcd590d, 0xbfd08a73667c57af, 0xbfd0a324e27390e3, 0xbfd0bbccdb0d24bd, 0xbfd0d46b579ab74b,
    0xbfd0ed005f657da4, 0xbfd1058bf9ae4ad5, 0xbfd11e0e2dad9cb7, 0xbfd136870293a8b0, 0xbfd14ef67f88685a,
    0xbfd1675cababa60e, 0xbfd17fb98e15095d, 0xbfd1980d2dd4236f, 0xbfd1b05791f07b49, 0xbfd1c898c16999fb,
    0xbfd1e0d0c33716be, 0xbfd1f8ff9e48a2f3, 0xbfd211255986160c, 0xbfd22941fbcf7966, 0xbfd241558bfd1404,
    0xbfd2596010df763a, 0xbfd27161913f853d, 0xbfd2895a13de86a3, 0xbfd2a1499f762bc9, 0xbfd2b9303ab89d25,
    0xbfd2d10dec508583, 0xbfd2e8e2bae11d31, 0xbfd300aead06350c, 0xbfd31871c9544185, 0xbfd3302c16586588,
    0xbfd347dd9a987d55, 0xbfd35f865c93293e, 0xbfd3772662bfd85b, 0xbfd38ebdb38ed321, 0xbfd3a64c556945ea,
    0xbfd3bdd24eb14b6a, 0xbfd3d54fa5c1f710, 0xbfd3ecc460ef5f50, 0xbfd404308686a7e4, 0xbfd41b941cce0bee,
    0xbfd432ef2a04e814, 0xbfd44a41b463c47c, 0xbfd4618bc21c5ec2, 0xbfd478cd5959b3d9, 0xbfd49006804009d1,
    0xbfd4a7373cecf997, 0xbfd4be5f957778a1, 0xbfd4d57f8fefe27f, 0xbfd4ec973260026a, 0xbfd503a682cb1cb3,
    0xbfd51aad872df82d, 0xbfd531ac457ee77e, 0xbfd548a2c3add263, 0xbfd55f9107a43ee2, 0xbfd5767717455a6c,
    0xbfd58d54f86e02f2, 0xbfd5a42ab0f4cfe2, 0xbfd5baf846aa1b19, 0xbfd5d1bdbf5809ca, 0xbfd5e87b20c2954a,
    0xbfd5ff3070a793d4, 0xbfd615ddb4bec13c, 0xbfd62c82f2b9c795, 0x3fd61965cdb02c1f, 0x3fd602d08af091ec,
    0x3fd5ec433d5c35ae, 0x3fd5d5bddf595f30, 0x3fd5bf406b543db2, 0x3fd5a8cadbbedfa1, 0x3fd5925d2b112a59,
    0x3fd57bf753c8d1fb, 0x3fd565995069514c, 0x3fd54f431b7be1a9, 0x3fd538f4af8f72fe, 0x3fd522ae0738a3d8,
    0x3fd50c6f1d11b97c, 0x3fd4f637ebba9810, 0x3fd4e0086dd8baca, 0x3fd4c9e09e172c3c, 0x3fd4b3c077267e9a,
    0x3fd49da7f3bcc41f, 0x3fd487970e958770, 0x3fd4718dc271c41b, 0x3fd45b8c0a17df13, 0x3fd44591e0539f49,
    0x3fd42f9f3ff62642, 0x3fd419b423d5e8c7, 0x3fd403d086cea79c, 0x3fd3edf463c1683e, 0x3fd3d81fb5946dba,
    0x3fd3c25277333184, 0x3fd3ac8ca38e5c5f, 0x3fd396ce359bbf54, 0x3fd3811728564cb2, 0x3fd36b6776be1117,
    0x3fd355bf1bd82c8b, 0x3fd3401e12aecba1, 0x3fd32a84565120a8, 0x3fd314f1e1d35ce4, 0x3fd2ff66b04ea9d4,
    0x3fd2e9e2bce12286, 0x3fd2d46602adccee, 0x3fd2bef07cdc9354, 0x3fd2a982269a3dbf, 0x3fd2941afb186b7c,
    0x3fd27ebaf58d8c9d, 0x3fd269621134db92, 0x3fd25410494e56c7, 0x3fd23ec5991eba49, 0x3fd22981fbef797b,
    0x3fd214456d0eb8d4, 0x3fd1ff0fe7cf47a7, 0x3fd1e9e1678899f4, 0x3fd1d4b9e796c245, 0x3fd1bf99635a6b95,
    0x3fd1aa7fd638d33f, 0x3fd1956d3b9bc2fa, 0x3fd180618ef18adf, 0x3fd16b5ccbacfb73, 0x3fd1565eed455fc3,
    0x3fd14167ef367783, 0x3fd12c77cd00713b, 0x3fd1178e8227e47c, 0x3fd102ac0a35cc1c, 0x3fd0edd060b78081,
    0x3fd0d8fb813eb1ef, 0x3fd0c42d676162e3, 0x3fd0af660eb9e279, 0x3fd09aa572e6c6d4, 0x3fd085eb8f8ae797,
    0x3fd07138604d5862, 0x3fd05c8be0d9635a, 0x3fd047e60cde83b8, 0x3fd03346e0106062, 0x3fd01eae5626c691,
    0x3fd00a1c6adda473, 0x3fcfeb2233ea07cd, 0x3fcfc218be620a5e, 0x3fcf991c6cb3b379, 0x3fcf702d36777df0,
    0x3fcf474b134df229, 0x3fcf1e75fadf9bde, 0x3fcef5ade4dcffe6, 0x3fceccf2c8fe920a, 0x3fcea4449f04aaf5,
    0x3fce7ba35eb77e2a, 0x3fce530effe71012, 0x3fce2a877a6b2c12, 0x3fce020cc6235ab5, 0x3fcdd99edaf6d7e9,
    0x3fcdb13db0d48940, 0x3fcd88e93fb2f450, 0x3fcd60a17f903515, 0x3fcd38666871f465, 0x3fcd1037f2655e7b,
    0x3fcce816157f1988, 0x3fccc000c9db3c52, 0x3fcc97f8079d44ec, 0x3fcc6ffbc6f00f71, 0x3fcc480c0005ccd1,
    0x3fcc2028ab17f9b4, 0x3fcbf851c067555f, 0x3fcbd087383bd8ad, 0x3fcba8c90ae4ad19, 0x3fcb811730b823d2,
    0x3fcb5971a213acdb, 0x3fcb31d8575bce3d, 0x3fcb0a4b48fc1b46, 0x3fcae2ca6f672bd4, 0x3fcabb55c31693ad,
    0x3fca93ed3c8ad9e3, 0x3fca6c90d44b704e, 0x3fca454082e6ab05, 0x3fca1dfc40f1b7f1, 0x3fc9f6c407089664,
    0x3fc9cf97cdce0ec3, 0x3fc9a8778debaa38, 0x3fc981634011aa75, 0x3fc95a5adcf7017f, 0x3fc9335e5d594989,
    0x3fc90c6db9fcbcd9, 0x3fc8e588ebac2dbf, 0x3fc8beafeb38fe8c, 0x3fc897e2b17b19a5, 0x3fc871213750e994,
    0x3fc84a6b759f512f, 0x3fc823c16551a3c2, 0x3fc7fd22ff599d4f, 0x3fc7d6903caf5ad0, 0x3fc7b0091651528c,
    0x3fc7898d85444c73, 0x3fc7631d82935a86, 0x3fc73cb9074fd14d, 0x3fc716600c914054, 0x3fc6f0128b756abc,
    0x3fc6c9d07d203fc7, 0x3fc6a399dabbd383, 0x3fc67d6e9d785771, 0x3fc6574ebe8c133a, 0x3fc6313a37335d76,
    0x3fc60b3100b09476, 0x3fc5e533144c1719, 0x3fc5bf406b543db2, 0x3fc59958ff1d52f1, 0x3fc5737cc9018cdd,
    0x3fc54dabc26105d2, 0x3fc527e5e4a1b58d, 0x3fc5022b292f6a45, 0x3fc4dc7b897bc1c8, 0x3fc4b6d6fefe22a4,
    0x3fc4913d8333b561, 0x3fc46baf0f9f5db7, 0x3fc4462b9dc9b3dc, 0x3fc420b32740fdd4, 0x3fc3fb45a59928cc,
    0x3fc3d5e3126bc27f, 0x3fc3b08b6757f2a9, 0x3fc38b3e9e027479, 0x3fc365fcb0159016, 0x3fc340c59741142e,
    0x3fc31b994d3a4f85, 0x3fc2f677cbbc0a96, 0x3fc2d1610c86813a, 0x3fc2ac55095f5c59, 0x3fc28753bc11aba5,
    0x3fc2625d1e6ddf57, 0x3fc23d712a49c202, 0x3fc2188fd9807263, 0x3fc1f3b925f25d41, 0x3fc1ceed09853752,
    0x3fc1aa2b7e23f72a, 0x3fc185747dbecf34, 0x3fc160c8024b27b1, 0x3fc13c2605c398c3, 0x3fc1178e8227e47c,
    0x3fc0f301717cf0fb, 0x3fc0ce7ecdccc28d, 0x3fc0aa06912675d5, 0x3fc08598b59e3a07, 0x3fc06135354d4b18,
    0x3fc03cdc0a51ec0d, 0x3fc0188d2ecf6140, 0x3fbfe89139dbd566, 0x3fbfa01c9db57ce2, 0x3fbf57bc7d9005db,
    0x3fbf0f70cdd992e3, 0x3fbec739830a1120, 0x3fbe7f1691a32d3e, 0x3fbe3707ee30487b, 0x3fbdef0d8d466db9,
    0x3fbda727638446a2, 0x3fbd5f55659210e2, 0x3fbd179788219364, 0x3fbccfedbfee13a8, 0x3fbc885801bc4b23,
    0x3fbc40d6425a5cb1, 0x3fbbf968769fca11, 0x3fbbb20e936d6974, 0x3fbb6ac88dad5b1c, 0x3fbb23965a52ff00,
    0x3fbadc77ee5aea8c, 0x3fba956d3ecade63, 0x3fba4e7640b1bc38, 0x3fba0792e9277cac, 0x3fb9c0c32d4d2548,
    0x3fb97a07024cbe74, 0x3fb9335e5d594989, 0x3fb8ecc933aeb6e8, 0x3fb8a6477a91dc29, 0x3fb85fd927506a48,
    0x3fb8197e2f40e3f0, 0x3fb7d33687c293c9, 0x3fb78d02263d82d3, 0x3fb746e100226ed9, 0x3fb700d30aeac0e1,
    0x3fb6bad83c1883b6, 0x3fb674f089365a7a, 0x3fb62f1be7d77743, 0x3fb5e95a4d9791cb, 0x3fb5a3abb01ade25,
    0x3fb55e10050e0384, 0x3fb518874226130a, 0x3fb4d3115d207eac, 0x3fb48dae4bc31018, 0x3fb4485e03dbdfad,
    0x3fb403207b414b7f, 0x3fb3bdf5a7d1ee64, 0x3fb378dd7f749714, 0x3fb333d7f8183f4b, 0x3fb2eee507b40301,
    0x3fb2aa04a44717a5, 0x3fb26536c3d8c369, 0x3fb2207b5c78549e, 0x3fb1dbd2643d190b, 0x3fb1973bd1465567,
    0x3fb152b799bb3cc9, 0x3fb10e45b3cae831, 0x3fb0c9e615ac4e17, 0x3fb08598b59e3a07, 0x3fb0415d89e74444,
    0x3faffa6911ab9301, 0x3faf723b517fc523, 0x3faeea31c006b87c, 0x3fae624c4a0b5e1b, 0x3fadda8adc67ee4e,
    0x3fad52ed6405d86f, 0x3faccb73cdddb2cc, 0x3fac441e06f72a9e, 0x3fabbcebfc68f420, 0x3fab35dd9b58baad,
    0x3faaaef2d0fb10fc, 0x3faa282b8a936171, 0x3fa9a187b573de7c, 0x3fa91b073efd7314, 0x3fa894aa149fb343,
    0x3fa80e7023d8ccc4, 0x3fa788595a3577ba, 0x3fa70265a550e777, 0x3fa67c94f2d4bb58, 0x3fa5f6e73078efb8,
    0x3fa5715c4c03ceef, 0x3fa4ebf43349e26f, 0x3fa466aed42de3ea, 0x3fa3e18c1ca0ae92, 0x3fa35c8bfaa1306b,
    0x3fa2d7ae5c3c5bae, 0x3fa252f32f8d183f, 0x3fa1ce5a62bc353a, 0x3fa149e3e4005a8d, 0x3fa0c58fa19dfaaa,
    0x3fa0415d89e74444, 0x3f9f7a9b16782856, 0x3f9e72bf2813ce51, 0x3f9d6b2725979802, 0x3f9c63d2ec14aaf2,
    0x3f9b5cc258b718e6, 0x3f9a55f548c5c43f, 0x3f994f6b99a24475, 0x3f98492528c8cabf, 0x3f974321d3d006d3,
    0x3f963d6178690bd6, 0x3f9537e3f45f3565, 0x3f9432a925980cc1, 0x3f932db0ea132e22, 0x3f9228fb1fea2e28,
    0x3f912487a5507f70, 0x3f90205658935847, 0x3f8e38ce3033310c, 0x3f8c317384c75f06, 0x3f8a2a9c6c170462,
    0x3f882448a388a2aa, 0x3f861e77e8b53fc6, 0x3f841929f96832f0, 0x3f82145e939ef1e9, 0x3f8010157588de71,
    0x3f7c189cbb0e27fb, 0x3f78121214586b54, 0x3f740c8a747878e2, 0x3f70080559588b35, 0x3f680904828985c0,
    0x3f60040155d5889e, 0x3f50020055655889, 0x0000000000000000};
const uint64_t powf4_exp_tab[2049] = {
    0x3feffffffc27dd9e, 0x3ff00162f1a4047d, 0x3ff002c603f68252, 0x3ff00429350e12af, 0x3ff0058c84ed6032,
    0x3ff006eff39715b2, 0x3ff00853810dde41, 0x3ff009b72d54652f, 0x3ff00b1af86d5604, 0x3ff00c7ee25b5c86,
    0x3ff00de2eb2124b3, 0x3ff00f4712c15ac8, 0x3ff010ab593eab39, 0x3ff0120fbe9bc2ba, 0x3ff0137442db4e38,
    0x3ff014d8e5fffada, 0x3ff0163da80c7604, 0x3ff017a289036d56, 0x3ff0190788e78eab, 0x3ff01a6ca7bb8818,
    0x3ff01bd1e58207ef, 0x3ff01d37423dbcbc, 0x3ff01e9cbdf15549, 0x3ff02002589f8099, 0x3ff02168124aedec,
    0x3ff022cdeaf64cbc, 0x3ff02433e2a44cc1, 0x3ff02599f9579ded, 0x3ff027002f12f06d, 0x3ff0286683d8f4ac,
    0x3ff029ccf7ac5b4d, 0x3ff02b338a8fd532, 0x3ff02c9a3c861379, 0x3ff02e010d91c778, 0x3ff02f67fdb5a2c4,
    0x3ff030cf0cf4572d, 0x3ff032363b5096bc, 0x3ff0339d88cd13bc, 0x3ff03504f56c80ae, 0x3ff0366c81319053,
    0x3ff037d42c1ef5a2, 0x3ff0393bf63763d5, 0x3ff03aa3df7d8e5f, 0x3ff03c0be7f428eb, 0x3ff03d740f9de766,
    0x3ff03edc567d7df7, 0x3ff04044bc95a0fe, 0x3ff041ad41e9051d, 0x3ff04315e67a5f2a, 0x3ff0447eaa4c643e,
    0x3ff045e78d61c9ac, 0x3ff047508fbd4502, 0x3ff048b9b1618c0b, 0x3ff04a22f25154cd, 0x3ff04b8c528f558b,
    0x3ff04cf5d21e44c4, 0x3ff04e5f7100d935, 0x3ff04fc92f39c9d4, 0x3ff051330ccbcdd5, 0x3ff0529d09b99ca8,
    0x3ff054072605edfb, 0x3ff0557161b379b3, 0x3ff056dbbcc4f7f8, 0x3ff05846373d212a, 0x3ff059b0d11eade5,
    0x3ff05b1b8a6c5706, 0x3ff05c866328d5a2, 0x3ff05df15b56e30a, 0x3ff05f5c72f938cf, 0x3ff060c7aa1290bd,
    0x3ff0623300a5a4db, 0x3ff0639e76b52f6e, 0x3ff0650a0c43eaf6, 0x3ff06675c1549232, 0x3ff067e195e9e01a,
    0x3ff0694d8a068fe7, 0x3ff06ab99dad5d0c, 0x3ff06c25d0e10338, 0x3ff06d9223a43e58, 0x3ff06efe95f9ca95,
    0x3ff0706b27e46455, 0x3ff071d7d966c83a, 0x3ff07344aa83b324, 0x3ff074b19b3de22f, 0x3ff0761eab9812b4,
    0x3ff0778bdb950247, 0x3ff078f92b376ebc, 0x3ff07a669a821621, 0x3ff07bd42977b6c4, 0x3ff07d41d81b0f2b,
    0x3ff07eafa66ede1e, 0x3ff0801d9475e2a0, 0x3ff0818ba232dbee, 0x3ff082f9cfa88985, 0x3ff084681cd9ab21,
    0x3ff085d689c900b6, 0x3ff0874516794a79, 0x3ff088b3c2ed48d9, 0x3ff08a228f27bc86, 0x3ff08b917b2b6667,
    0x3ff08d0086fb07a6, 0x3ff08e6fb29961a8, 0x3ff08fdefe09360d, 0x3ff0914e694d46b6, 0x3ff092bdf46855c0,
    0x3ff0942d9f5d2582, 0x3ff0959d6a2e7893, 0x3ff0970d54df11c8, 0x3ff0987d5f71b432, 0x3ff099ed89e9231e,
    0x3ff09b5dd448221a, 0x3ff09cce3e9174ec, 0x3ff09e3ec8c7df9d, 0x3ff09faf72ee2670, 0x3ff0a1203d070de5,
    0x3ff0a29127155abd, 0x3ff0a402311bd1f0, 0x3ff0a5735b1d38bb, 0x3ff0a6e4a51c5493, 0x3ff0a8560f1beb2c,
    0x3ff0a9c7991ec278, 0x3ff0ab394327a0a7, 0x3ff0acab0d394c25, 0x3ff0ae1cf7568b9d, 0x3ff0af8f018225f7,
    0x3ff0b1012bbee259, 0x3ff0b273760f8825, 0x3ff0b3e5e076defc, 0x3ff0b5586af7aebc, 0x3ff0b6cb1594bf84,
    0x3ff0b83de050d9ab, 0x3ff0b9b0cb2ec5ca, 0x3ff0bb23d6314cb7, 0x3ff0bc97015b3783, 0x3ff0be0a4caf4f81,
    0x3ff0bf7db8305e3f, 0x3ff0c0f143e12d8a, 0x3ff0c264efc4876c, 0x3ff0c3d8bbdd362e, 0x3ff0c54ca82e0455,
    0x3ff0c6c0b4b9bca6, 0x3ff0c834e1832a24, 0x3ff0c9a92e8d180e, 0x3ff0cb1d9bda51e1, 0x3ff0cc92296da35b,
    0x3ff0ce06d749d876, 0x3ff0cf7ba571bd6a, 0x3ff0d0f093e81eab, 0x3ff0d265a2afc8f1, 0x3ff0d3dad1cb892b,
    0x3ff0d550213e2c8c, 0x3ff0d6c5910a8081, 0x3ff0d83b213352b8, 0x3ff0d9b0d1bb711b, 0x3ff0db26a2a5a9d4,
    0x3ff0dc9c93f4cb4a, 0x3ff0de12a5aba423, 0x3ff0df88d7cd0344, 0x3ff0e0ff2a5bb7cd, 0x3ff0e2759d5a9121,
    0x3ff0e3ec30cc5edd, 0x3ff0e562e4b3f0df, 0x3ff0e6d9b9141745, 0x3ff0e850adefa265, 0x3ff0e9c7c34962db,
    0x3ff0eb3ef924297d, 0x3ff0ecb64f82c75e, 0x3ff0ee2dc6680dd6, 0x3ff0efa55dd6ce75, 0x3ff0f11d15d1db0c,
    0x3ff0f294ee5c05ab, 0x3ff0f40ce77820a2, 0x3ff0f5850128fe7a, 0x3ff0f6fd3b717200, 0x3ff0f87596544e3f,
    0x3ff0f9ee11d4667f, 0x3ff0fb66adf48e46, 0x3ff0fcdf6ab7995c, 0x3ff0fe5848205bc4, 0x3ff0ffd14631a9c2,
    0x3ff1014a64ee57d9, 0x3ff102c3a4593ac9, 0x3ff1043d04752792, 0x3ff105b68544f373, 0x3ff1073026cb73e9,
    0x3ff108a9e90b7eb2, 0x3ff10a23cc07e9c6, 0x3ff10b9dcfc38b63, 0x3ff10d17f44139fe, 0x3ff10e923983cc53,
    0x3ff1100c9f8e1955, 0x3ff111872662f83e, 0x3ff11301ce054081, 0x3ff1147c9677c9d2, 0x3ff115f77fbd6c23,
    0x3ff1177289d8ffa9, 0x3ff118edb4cd5cd4, 0x3ff11a69009d5c54, 0x3ff11be46d4bd71a, 0x3ff11d5ffadba653,
    0x3ff11edba94fa36e, 0x3ff1205778aaa817, 0x3ff121d368ef8e3b, 0x3ff1234f7a213008, 0x3ff124cbac4267e5,
    0x3ff12647ff56107f, 0x3ff127c4735f04bd, 0x3ff1294108601fcb, 0x3ff12abdbe5c3d0f, 0x3ff12c3a95563832,
    0x3ff12db78d50ed19, 0x3ff12f34a64f37ed, 0x3ff130b1e053f513, 0x3ff1322f3b62012e, 0x3ff133acb77c3927,
    0x3ff1352a54a57a1f, 0x3ff136a812e0a17c, 0x3ff13825f2308ce0, 0x3ff139a3f2981a2e, 0x3ff13b22141a278a,
    0x3ff13ca056b99356, 0x3ff13e1eba793c33, 0x3ff13f9d3f5c0103, 0x3ff1411be564c0e7, 0x3ff1429aac965b40,
    0x3ff1441994f3afae, 0x3ff145989e7f9e13, 0x3ff14717c93d068e, 0x3ff14897152ec980, 0x3ff14a168257c787,
    0x3ff14b9610bae185, 0x3ff14d15c05af897, 0x3ff14e95913aee1e, 0x3ff15015835da3b8, 0x3ff1519596c5fb46,
    0x3ff15315cb76d6e5, 0x3ff15496217318f6, 0x3ff1561698bda417, 0x3ff1579731595b27, 0x3ff15917eb492145,
    0x3ff15a98c68fd9d1, 0x3ff15c19c330686b, 0x3ff15d9ae12db0ef, 0x3ff15f1c208a977f, 0x3ff1609d814a007b,
    0x3ff1621f036ed081, 0x3ff163a0a6fbec71, 0x3ff165226bf4396d, 0x3ff166a4525a9cd3, 0x3ff168265a31fc44,
    0x3ff169a8837d3da3, 0x3ff16b2ace3f4710, 0x3ff16cad3a7afeeb, 0x3ff16e2fc8334bd8, 0x3ff16fb2776b14b8,
    0x3ff17135482540ad, 0x3ff172b83a64b71a, 0x3ff1743b4e2c5fa4, 0x3ff175be837f222d, 0x3ff17741da5fe6da,
    0x3ff178c552d1960f, 0x3ff17a48ecd71873, 0x3ff17bcca87356e9, 0x3ff17d5085a93a9b, 0x3ff17ed4847bacec,
    0x3ff18058a4ed9787, 0x3ff181dce701e451, 0x3ff183614abb7d75, 0x3ff184e5d01d4d5b, 0x3ff1866a772a3ead,
    0x3ff187ef3fe53c58, 0x3ff189742a513185, 0x3ff18af9367109a1, 0x3ff18c7e6447b059, 0x3ff18e03b3d8119c,
    0x3ff18f8925251997, 0x3ff1910eb831b4ba, 0x3ff192946d00cfb6, 0x3ff1941a4395577c, 0x3ff195a03bf2393e,
    0x3ff19726561a626d, 0x3ff198ac9210c0c1, 0x3ff19a32efd8422c, 0x3ff19bb96f73d4e5, 0x3ff19d4010e66763,
    0x3ff19ec6d432e85c, 0x3ff1a04db95c46cc, 0x3ff1a1d4c06571ed, 0x3ff1a35be9515937, 0x3ff1a4e33422ec69,
    0x3ff1a66aa0dd1b81, 0x3ff1a7f22f82d6bc, 0x3ff1a979e0170e9a, 0x3ff1ab01b29cb3dd, 0x3ff1ac89a716b786,
    0x3ff1ae11bd880ada, 0x3ff1af99f5f39f5d, 0x3ff1b122505c66d5, 0x3ff1b2aaccc5534b, 0x3ff1b4336b315705,
    0x3ff1b5bc2ba3648e, 0x3ff1b7450e1e6eb3, 0x3ff1b8ce12a56880, 0x3ff1ba57393b4544, 0x3ff1bbe081e2f88e,
    0x3ff1bd69ec9f762f, 0x3ff1bef37973b23b, 0x3ff1c07d2862a105, 0x3ff1c206f96f3724, 0x3ff1c390ec9c696f,
    0x3ff1c51b01ed2cfe, 0x3ff1c6a53964772e, 0x3ff1c82f93053d99, 0x3ff1c9ba0ed2761e, 0x3ff1cb44accf16dc,
    0x3ff1cccf6cfe1634, 0x3ff1ce5a4f626acb, 0x3ff1cfe553ff0b83, 0x3ff1d1707ad6ef85, 0x3ff1d2fbc3ed0e37,
    0x3ff1d4872f445f44, 0x3ff1d612bcdfda99, 0x3ff1d79e6cc27863, 0x3ff1d92a3eef3111, 0x3ff1dab63368fd56,
    0x3ff1dc424a32d624, 0x3ff1ddce834fb4b0, 0x3ff1df5adec29273, 0x3ff1e0e75c8e6927, 0x3ff1e273fcb632c5,
    0x3ff1e400bf3ce98b, 0x3ff1e58da42587fa, 0x3ff1e71aab7308d1, 0x3ff1e8a7d5286717, 0x3ff1ea3521489e0e,
    0x3ff1ebc28fd6a942, 0x3ff1ed5020d5847a, 0x3ff1eeddd4482bc3, 0x3ff1f06baa319b6b, 0x3ff1f1f9a294d004,
    0x3ff1f387bd74c660, 0x3ff1f515fad47b95, 0x3ff1f6a45ab6ecfa, 0x3ff1f832dd1f1829, 0x3ff1f9c1820ffafe,
    0x3ff1fb50498c9397, 0x3ff1fcdf3397e057, 0x3ff1fe6e4034dfdf, 0x3ff1fffd6f669115, 0x3ff2018cc12ff324,
    0x3ff2031c35940574, 0x3ff204abcc95c7b4, 0x3ff2063b863839d2, 0x3ff207cb627e5c01, 0x3ff2095b616b2eb7,
    0x3ff20aeb8301b2aa, 0x3ff20c7bc744e8d5, 0x3ff20e0c2e37d274, 0x3ff20f9cb7dd7108, 0x3ff2112d6438c651,
    0x3ff212be334cd455, 0x3ff2144f251c9d5b, 0x3ff215e039ab23ee, 0x3ff2177170fb6adb, 0x3ff21902cb107531,
    0x3ff21a9447ed4643, 0x3ff21c25e794e1a7, 0x3ff21db7aa0a4b34, 0x3ff21f498f508707, 0x3ff220db976a997e,
    0x3ff2226dc25b8739, 0x3ff224001026551c, 0x3ff2259280ce084e, 0x3ff227251455a63b, 0x3ff228b7cac0348e,
    0x3ff22a4aa410b938, 0x3ff22bdda04a3a6b, 0x3ff22d70bf6fbea0, 0x3ff22f0401844c8d, 0x3ff23097668aeb31,
    0x3ff2322aee86a1ca, 0x3ff233be997a77db, 0x3ff235526769752b, 0x3ff236e65856a1c4, 0x3ff2387a6c4505ef,
    0x3ff23a0ea337aa40, 0x3ff23ba2fd319789, 0x3ff23d377a35d6de, 0x3ff23ecc1a47719b, 0x3ff24060dd69715c,
    0x3ff241f5c39ee002, 0x3ff2438acceac7b2, 0x3ff2451ff95032d1, 0x3ff246b548d22c0c, 0x3ff2484abb73be50,
    0x3ff249e05137f4cf, 0x3ff24b760a21daff, 0x3ff24d0be6347c98, 0x3ff24ea1e572e597, 0x3ff2503807e0223a,
    0x3ff251ce4d7f3f08, 0x3ff25364b65348c6, 0x3ff254fb425f4c80, 0x3ff25691f1a65784, 0x3ff25828c42b7763,
    0x3ff259bfb9f1b9f7, 0x3ff25b56d2fc2d55, 0x3ff25cee0f4ddfdd, 0x3ff25e856ee9e031, 0x3ff2601cf1d33d35,
    0x3ff261b4980d0613, 0x3ff2634c619a4a36, 0x3ff264e44e7e1952, 0x3ff2667c5ebb8358, 0x3ff2681492559883,
    0x3ff269ace94f694f, 0x3ff26b4563ac067d, 0x3ff26cde016e8111, 0x3ff26e76c299ea53, 0x3ff2700fa73153d1,
    0x3ff271a8af37cf5b, 0x3ff27341dab06f07, 0x3ff274db299e452c, 0x3ff276749c046468, 0x3ff2780e31e5df9c,
    0x3ff279a7eb45c9ef, 0x3ff27b41c82736c8, 0x3ff27cdbc88d39d8, 0x3ff27e75ec7ae70f, 0x3ff2801033f352a4,
    0x3ff281aa9ef99111, 0x3ff283452d90b716, 0x3ff284dfdfbbd9b7, 0x3ff2867ab57e0e3a, 0x3ff28815aeda6a2d,
    0x3ff289b0cbd4035f, 0x3ff28b4c0c6defe6, 0x3ff28ce770ab461b, 0x3ff28e82f88f1c9e, 0x3ff2901ea41c8a50,
    0x3ff291ba7356a657, 0x3ff2935666408820, 0x3ff294f27cdd475a, 0x3ff2968eb72ffbfc, 0x3ff2982b153bbe3d,
    0x3ff299c79703a69e, 0x3ff29b643c8acdde, 0x3ff29d0105d44d08, 0x3ff29e9df2e33d66, 0x3ff2a03b03bab88b,
    0x3ff2a1d8385dd84a, 0x3ff2a37590cfb6bf, 0x3ff2a5130d136e49, 0x3ff2a6b0ad2c198d, 0x3ff2a84e711cd374,
    0x3ff2a9ec58e8b729, 0x3ff2ab8a6492e024, 0x3ff2ad28941e6a18, 0x3ff2aec6e78e7104, 0x3ff2b0655ee61129,
    0x3ff2b203fa28670e, 0x3ff2b3a2b9588f7f, 0x3ff2b5419c79a78c, 0x3ff2b6e0a38ecc8b, 0x3ff2b87fce9b1c18,
    0x3ff2ba1f1da1b412, 0x3ff2bbbe90a5b29e, 0x3ff2bd5e27aa3627, 0x3ff2befde2b25d5c, 0x3ff2c09dc1c14733,
    0x3ff2c23dc4da12e4, 0x3ff2c3ddebffdff0, 0x3ff2c57e3735ce1a, 0x3ff2c71ea67efd6d, 0x3ff2c8bf39de8e38,
    0x3ff2ca5ff157a10f, 0x3ff2cc00cced56cd, 0x3ff2cda1cca2d08f, 0x3ff2cf42f07b2fbd, 0x3ff2d0e4387995fe,
    0x3ff2d285a4a12544, 0x3ff2d42734f4ffc3, 0x3ff2d5c8e97847f6, 0x3ff2d76ac22e209d, 0x3ff2d90cbf19acbd,
    0x3ff2daaee03e0fa3, 0x3ff2dc51259e6ce0, 0x3ff2ddf38f3de848, 0x3ff2df961d1fa5f9, 0x3ff2e138cf46ca57,
    0x3ff2e2dba5b67a06, 0x3ff2e47ea071d9f4, 0x3ff2e621bf7c0f57, 0x3ff2e7c502d83fa4, 0x3ff2e9686a89909e,
    0x3ff2eb0bf6932845, 0x3ff2ecafa6f82ce7, 0x3ff2ee537bbbc512, 0x3ff2eff774e1179d, 0x3ff2f19b926b4ba5,
    0x3ff2f33fd45d888a, 0x3ff2f4e43abaf5f5, 0x3ff2f688c586bbd5, 0x3ff2f82d74c4025c, 0x3ff2f9d24875f205,
    0x3ff2fb77409fb390, 0x3ff2fd1c5d447003, 0x3ff2fec19e6750ab, 0x3ff30067040b7f19, 0x3ff3020c8e342527,
    0x3ff303b23ce46cf4, 0x3ff30558101f80e3, 0x3ff306fe07e88ba0, 0x3ff308a42442b81d, 0x3ff30a4a65313191,
    0x3ff30bf0cab7237a, 0x3ff30d9754d7b99d, 0x3ff30f3e03962005, 0x3ff310e4d6f58302, 0x3ff3128bcef90f2b,
    0x3ff31432eba3f15f, 0x3ff315da2cf956c2, 0x3ff3178192fc6cbf, 0x3ff319291db06106, 0x3ff31ad0cd186190,
    0x3ff31c78a1379c9b, 0x3ff31e209a1140ab, 0x3ff31fc8b7a87c8c, 0x3ff32170fa007f51, 0x3ff32319611c7851,
    0x3ff324c1ecff972d, 0x3ff3266a9dad0bca, 0x3ff3281373280654, 0x3ff329bc6d73b741, 0x3ff32b658c934f48,
    0x3ff32d0ed089ff6c, 0x3ff32eb8395af8f4, 0x3ff33061c7096d6f, 0x3ff3320b79988eb2, 0x3ff333b5510b8eda,
    0x3ff3355f4d65a04a, 0x3ff337096ea9f5ab, 0x3ff338b3b4dbc1f0, 0x3ff33a5e1ffe384f, 0x3ff33c08b0148c49,
    0x3ff33db36521f1a3, 0x3ff33f5e3f299c69, 0x3ff341093e2ec0f2, 0x3ff342b4623493d7, 0x3ff3445fab3e49fa,
    0x3ff3460b194f1887, 0x3ff347b6ac6a34ec, 0x3ff349626492d4e2, 0x3ff34b0e41cc2e67, 0x3ff34cba441977c4,
    0x3ff34e666b7de784, 0x3ff35012b7fcb47d, 0x3ff351bf299915c9, 0x3ff3536bc05642cf, 0x3ff355187c377337,
    0x3ff356c55d3fdef4, 0x3ff358726372be40, 0x3ff35a1f8ed3499b, 0x3ff35bccdf64b9cf, 0x3ff35d7a552a47ea,
    0x3ff35f27f0272d45, 0x3ff360d5b05ea37f, 0x3ff3628395d3e47d, 0x3ff36431a08a2a6f, 0x3ff365dfd084afc9,
    0x3ff3678e25c6af48, 0x3ff3693ca05363f3, 0x3ff36aeb402e0915, 0x3ff36c9a0559da43, 0x3ff36e48efda1358,
    0x3ff36ff7ffb1f078, 0x3ff371a734e4ae11, 0x3ff373568f7588d3, 0x3ff375060f67bdb9, 0x3ff376b5b4be8a0a,
    0x3ff378657f7d2b4c, 0x3ff37a156fa6df54, 0x3ff37bc5853ee43d, 0x3ff37d75c0487869, 0x3ff37f2620c6da82,
    0x3ff380d6a6bd497c, 0x3ff38287522f048e, 0x3ff38438231f4b3f, 0x3ff385e919915d57, 0x3ff3879a35887ae9,
    0x3ff3894b7707e450, 0x3ff38afcde12da2f, 0x3ff38cae6aac9d71, 0x3ff38e601cd86f48, 0x3ff39011f4999132,
    0x3ff391c3f1f344f1, 0x3ff3937614e8cc90, 0x3ff395285d7d6a65, 0x3ff396dacbb4610c, 0x3ff3988d5f90f36a,
    0x3ff39a40191664ac, 0x3ff39bf2f847f847, 0x3ff39da5fd28f1f8, 0x3ff39f5927bc95c8, 0x3ff3a10c78062804,
    0x3ff3a2bfee08ed45, 0x3ff3a47389c82a68, 0x3ff3a6274b472498, 0x3ff3a7db32892144, 0x3ff3a98f3f916626,
    0x3ff3ab4372633941, 0x3ff3acf7cb01e0de, 0x3ff3aeac4970a392, 0x3ff3b060edb2c837, 0x3ff3b215b7cb95f3,
    0x3ff3b3caa7be5434, 0x3ff3b57fbd8e4aae, 0x3ff3b734f93ec163, 0x3ff3b8ea5ad30097, 0x3ff3ba9fe24e50dd,
    0x3ff3bc558fb3fb0d, 0x3ff3be0b6307484a, 0x3ff3bfc15c4b81ff, 0x3ff3c1777b83f1e0, 0x3ff3c32dc0b3e1ea,
    0x3ff3c4e42bde9c62, 0x3ff3c69abd076bd8, 0x3ff3c85174319b24, 0x3ff3ca0851607565, 0x3ff3cbbf54974607,
    0x3ff3cd767dd958bd, 0x3ff3cf2dcd29f984, 0x3ff3d0e5428c749e, 0x3ff3d29cde04169e, 0x3ff3d4549f942c57,
    0x3ff3d60c874002ee, 0x3ff3d7c4950ae7cb, 0x3ff3d97cc8f828a2, 0x3ff3db35230b136f, 0x3ff3dceda346f679,
    0x3ff3dea649af204f, 0x3ff3e05f1646dfca, 0x3ff3e2180911840c, 0x3ff3e3d122125c7f, 0x3ff3e58a614cb8db,
    0x3ff3e743c6c3e91b, 0x3ff3e8fd527b3d8a, 0x3ff3eab7047606b7, 0x3ff3ec70dcb7957e, 0x3ff3ee2adb433b04,
    0x3ff3efe5001c48b5, 0x3ff3f19f4b46104c, 0x3ff3f359bcc3e3c8, 0x3ff3f51454991573, 0x3ff3f6cf12c8f7e5,
    0x3ff3f889f756ddfa, 0x3ff3fa4502461adb, 0x3ff3fc00339a01f9, 0x3ff3fdbb8b55e710, 0x3ff3ff77097d1e26,
    0x3ff40132ae12fb8a, 0x3ff402ee791ad3d5, 0x3ff404aa6a97fbea, 0x3ff40666828dc8f6, 0x3ff40822c0ff9071,
    0x3ff409df25f0a81b, 0x3ff40b9bb16465fe, 0x3ff40d58635e2070, 0x3ff40f153be12e0f, 0x3ff410d23af0e5c5,
    0x3ff4128f60909ec4, 0x3ff4144cacc3b08a, 0x3ff4160a1f8d72dd, 0x3ff417c7b8f13dd0, 0x3ff4198578f269be,
    0x3ff41b435f944f4c, 0x3ff41d016cda476c, 0x3ff41ebfa0c7ab57, 0x3ff4207dfb5fd491, 0x3ff4223c7ca61cea,
    0x3ff423fb249dde7b, 0x3ff425b9f34a73a7, 0x3ff42778e8af371d, 0x3ff4293804cf83d5, 0x3ff42af747aeb511,
    0x3ff42cb6b1502661, 0x3ff42e7641b73399, 0x3ff43035f8e738de, 0x3ff431f5d6e3929c, 0x3ff433b5dbaf9d8b,
    0x3ff43576074eb6ac, 0x3ff4373659c43b4c, 0x3ff438f6d3138902, 0x3ff43ab7733ffdb1, 0x3ff43c783a4cf784,
    0x3ff43e39283dd4f5, 0x3ff43ffa3d15f4c3, 0x3ff441bb78d8b5fc, 0x3ff4437cdb8977f9, 0x3ff4453e652b9a59,
    0x3ff4470015c27d0e, 0x3ff448c1ed51804b, 0x3ff44a83ebdc0497, 0x3ff44c4611656abf, 0x3ff44e085df113da,
    0x3ff44fcad182614e, 0x3ff4518d6c1cb4c9, 0x3ff453502dc37045, 0x3ff455131679f608, 0x3ff456d62643a8a0,
    0x3ff458995d23eae9, 0x3ff45a5cbb1e2009, 0x3ff45c204035ab72, 0x3ff45de3ec6df0df, 0x3ff45fa7bfca5459,
    0x3ff4616bba4e3a34, 0x3ff4632fdbfd070c, 0x3ff464f424da1fcc, 0x3ff466b894e8e9a6, 0x3ff4687d2c2cca1e,
    0x3ff46a41eaa926fc, 0x3ff46c06d061665a, 0x3ff46dcbdd58ee98, 0x3ff46f9111932664, 0x3ff471566d1374b7,
    0x3ff4731befdd40d6, 0x3ff474e199f3f251, 0x3ff476a76b5af103, 0x3ff4786d6415a514, 0x3ff47a33842776f6,
    0x3ff47bf9cb93cf67, 0x3ff47dc03a5e1772, 0x3ff47f86d089b86d, 0x3ff4814d8e1a1bf8, 0x3ff483147312ac00,
    0x3ff484db7f76d2be, 0x3ff486a2b349fab7, 0x3ff4886a0e8f8ebb, 0x3ff48a31914af9e6, 0x3ff48bf93b7fa79f,
    0x3ff48dc10d31039b, 0x3ff48f89066279da, 0x3ff49151271776a6, 0x3ff493196f536698, 0x3ff494e1df19b693,
    0x3ff496aa766dd3c6, 0x3ff4987335532bad, 0x3ff49a3c1bcd2c0f, 0x3ff49c0529df4300, 0x3ff49dce5f8cdee0,
    0x3ff49f97bcd96e5a, 0x3ff4a16141c86066, 0x3ff4a32aee5d2449, 0x3ff4a4f4c29b2993, 0x3ff4a6bebe85e020,
    0x3ff4a888e220b819, 0x3ff4aa532d6f21f4, 0x3ff4ac1da0748e6f, 0x3ff4ade83b346e9c, 0x3ff4afb2fdb233d4,
    0x3ff4b17de7f14fbb, 0x3ff4b348f9f53446, 0x3ff4b51433c153b3, 0x3ff4b6df9559208f, 0x3ff4b8ab1ec00db1,
    0x3ff4ba76cff98e3e, 0x3ff4bc42a90915a7, 0x3ff4be0ea9f217aa, 0x3ff4bfdad2b80852, 0x3ff4c1a7235e5bf6,
    0x3ff4c3739be88737, 0x3ff4c5403c59ff09, 0x3ff4c70d04b638a6, 0x3ff4c8d9f500a999, 0x3ff4caa70d3cc7b9,
    0x3ff4cc744d6e0926, 0x3ff4ce41b597e454, 0x3ff4d00f45bdcffe, 0x3ff4d1dcfde3432d, 0x3ff4d3aade0bb539,
    0x3ff4d578e63a9dc2, 0x3ff4d747167374bd, 0x3ff4d9156eb9b264, 0x3ff4dae3ef10cf42, 0x3ff4dcb2977c442f,
    0x3ff4de8167ff8a4e, 0x3ff4e050609e1b11, 0x3ff4e21f815b7036, 0x3ff4e3eeca3b03c9, 0x3ff4e5be3b405022,
    0x3ff4e78dd46ecfe6, 0x3ff4e95d95c9fe0b, 0x3ff4eb2d7f5555ce, 0x3ff4ecfd911452bd, 0x3ff4eecdcb0a70b3,
    0x3ff4f09e2d3b2bd8, 0x3ff4f26eb7aa00a1, 0x3ff4f43f6a5a6bd2, 0x3ff4f610454fea79, 0x3ff4f7e1488df9f4,
    0x3ff4f9b2741817ee, 0x3ff4fb83c7f1c25e, 0x3ff4fd55441e778b, 0x3ff4ff26e8a1b608, 0x3ff500f8b57efcb6,
    0x3ff502caaab9cac1, 0x3ff5049cc8559fa7, 0x3ff5066f0e55fb31, 0x3ff508417cbe5d77, 0x3ff50a14139246db,
    0x3ff50be6d2d53812, 0x3ff50db9ba8ab21a, 0x3ff50f8ccab63642, 0x3ff51160035b4625, 0x3ff51333647d63ad,
    0x3ff51506ee201112, 0x3ff516daa046d0d6, 0x3ff518ae7af525ce, 0x3ff51a827e2e931a, 0x3ff51c56a9f69c2a,
    0x3ff51e2afe50c4b9, 0x3ff51fff7b4090d2, 0x3ff521d420c984ce, 0x3ff523a8eeef2553, 0x3ff5257de5b4f757,
    0x3ff52753051e801a, 0x3ff529284d2f4530, 0x3ff52afdbdeacc76, 0x3ff52cd357549c19, 0x3ff52ea919703a95,
    0x3ff5307f04412eb4, 0x3ff5325517caff8d, 0x3ff5342b54113485, 0x3ff53601b9175551, 0x3ff537d846e0e9f5,
    0x3ff539aefd717ac0, 0x3ff53b85dccc9053, 0x3ff53d5ce4f5b39a, 0x3ff53f3415f06dd2, 0x3ff5410b6fc04885,
    0x3ff542e2f268cd8c, 0x3ff544ba9ded870f, 0x3ff546927251ff84, 0x3ff5486a6f99c1b1, 0x3ff54a4295c858a6,
    0x3ff54c1ae4e14fc7, 0x3ff54df35ce832c3, 0x3ff54fcbfde08d9b, 0x3ff551a4c7cdec9c, 0x3ff5537dbab3dc60,
    0x3ff55556d695e9d6, 0x3ff557301b77a236, 0x3ff55909895c9309, 0x3ff55ae320484a28, 0x3ff55cbce03e55b7,
    0x3ff55e96c942442b, 0x3ff56070db57a44b, 0x3ff5624b16820529, 0x3ff564257ac4f625, 0x3ff56600082406f4,
    0x3ff567dabea2c792, 0x3ff569b59e44c851, 0x3ff56b90a70d99ce, 0x3ff56d6bd900ccf5, 0x3ff56f473421f304,
    0x3ff57122b8749d87, 0x3ff572fe65fc5e58, 0x3ff574da3cbcc7a0, 0x3ff576b63cb96bd8, 0x3ff5789265f5ddca,
    0x3ff57a6eb875b08c, 0x3ff57c4b343c7786, 0x3ff57e27d94dc66d, 0x3ff58004a7ad3148, 0x3ff581e19f5e4c6a,
    0x3ff583bec064ac79, 0x3ff5859c0ac3e669, 0x3ff587797e7f8f7c, 0x3ff589571b9b3d44, 0x3ff58b34e21a85a7,
    0x3ff58d12d200fed2, 0x3ff58ef0eb523f4a, 0x3ff590cf2e11ddde, 0x3ff592ad9a4371af, 0x3ff5948c2fea922e,
    0x3ff5966aef0ad71b, 0x3ff59849d7a7d883, 0x3ff59a28e9c52ec9, 0x3ff59c0825667299, 0x3ff59de78a8f3cf4,
    0x3ff59fc719432727, 0x3ff5a1a6d185cad3, 0x3ff5a386b35ac1e4, 0x3ff5a566bec5a699, 0x3ff5a746f3ca1381,
    0x3ff5a927526ba378, 0x3ff5ab07daadf1af, 0x3ff5ace88c9499a3, 0x3ff5aec968233721, 0x3ff5b0aa6d5d6649,
    0x3ff5b28b9c46c389, 0x3ff5b46cf4e2eb9d, 0x3ff5b64e77357b97, 0x3ff5b830234210d3, 0x3ff5ba11f90c4902,
    0x3ff5bbf3f897c221, 0x3ff5bdd621e81a81, 0x3ff5bfb87500f0c1, 0x3ff5c19af1e5e3d2, 0x3ff5c37d989a92f2,
    0x3ff5c56069229db4, 0x3ff5c7436381a3f7, 0x3ff5c92687bb45ed, 0x3ff5cb09d5d3241a, 0x3ff5cced4dccdf4c,
    0x3ff5ced0efac18aa, 0x3ff5d0b4bb7471a6, 0x3ff5d298b1298c02, 0x3ff5d47cd0cf09d4, 0x3ff5d6611a688d81,
    0x3ff5d8458df9b9bf, 0x3ff5da2a2b863193, 0x3ff5dc0ef3119855, 0x3ff5ddf3e49f91ad, 0x3ff5dfd90033c193,
    0x3ff5e1be45d1cc4f, 0x3ff5e3a3b57d567d, 0x3ff5e5894f3a0506, 0x3ff5e76f130b7d27, 0x3ff5e95500f5646d,
    0x3ff5eb3b18fb60b3, 0x3ff5ed215b21182a, 0x3ff5ef07c76a314e, 0x3ff5f0ee5dda52f4, 0x3ff5f2d51e752439,
    0x3ff5f4bc093e4c90, 0x3ff5f6a31e3973bf, 0x3ff5f88a5d6a41d9, 0x3ff5fa71c6d45f43, 0x3ff5fc595a7b74b4,
    0x3ff5fe4118632b33, 0x3ff60029008f2c1b, 0x3ff6021113032116, 0x3ff603f94fc2b41f, 0x3ff605e1b6d18f82,
    0x3ff607ca48335ddf, 0x3ff609b303ebca24, 0x3ff60b9be9fe7f93, 0x3ff60d84fa6f29bf, 0x3ff60f6e3541748a,
    0x3ff611579a790c29, 0x3ff613412a199d25, 0x3ff6152ae426d453, 0x3ff61714c8a45edf, 0x3ff618fed795ea44,
    0x3ff61ae910ff244e, 0x3ff61cd374e3bb1b, 0x3ff61ebe03475d1c, 0x3ff620a8bc2db914, 0x3ff622939f9a7e14,
    0x3ff6247ead915b83, 0x3ff62669e6160116, 0x3ff62855492c1ed7, 0x3ff62a40d6d76521, 0x3ff62c2c8f1b84a0,
    0x3ff62e1871fc2e52, 0x3ff630047f7d1386, 0x3ff631f0b7a1e5e1, 0x3ff633dd1a6e5753, 0x3ff635c9a7e61a25,
    0x3ff637b6600ce0ed, 0x3ff639a342e65e97, 0x3ff63b905076465f, 0x3ff63d7d88c04bd0, 0x3ff63f6aebc822cd,
    0x3ff6415879917f88, 0x3ff6434632201685, 0x3ff6453415779c9b, 0x3ff64722239bc6f3, 0x3ff649105c904b09,
    0x3ff64afec058dea9, 0x3ff64ced4ef937f3, 0x3ff64edc08750d5b, 0x3ff650caecd015a3, 0x3ff652b9fc0e07e3,
    0x3ff654a936329b85, 0x3ff656989b418844, 0x3ff658882b3e862e, 0x3ff65a77e62d4da4, 0x3ff65c67cc119759,
    0x3ff65e57dcef1c54, 0x3ff6604818c995ed, 0x3ff662387fa4bdce, 0x3ff6642911844df6, 0x3ff66619ce6c00b4,
    0x3ff6680ab65f90ac, 0x3ff669fbc962b8d3, 0x3ff66bed07793473, 0x3ff66dde70a6bf26, 0x3ff66fd004ef14db,
    0x3ff671c1c455f1d2, 0x3ff673b3aedf129f, 0x3ff675a5c48e342a, 0x3ff67798056713aa, 0x3ff6798a716d6eaf,
    0x3ff67b7d08a50316, 0x3ff67d6fcb118f12, 0x3ff67f62b8b6d12a, 0x3ff68155d1988835, 0x3ff6834915ba7361,
    0x3ff6853c8520522a, 0x3ff687301fcde464, 0x3ff68923e5c6ea33, 0x3ff68b17d70f2412, 0x3ff68d0bf3aa52c9,
    0x3ff68f003b9c3779, 0x3ff690f4aee89395, 0x3ff692e94d9328e0, 0x3ff694de179fb976, 0x3ff696d30d1207c0,
    0x3ff698c82dedd681, 0x3ff69abd7a36e8c9, 0x3ff69cb2f1f101ff, 0x3ff69ea8951fe5e0, 0x3ff6a09e63c75876,
    0x3ff6a2945deb1e23, 0x3ff6a48a838efb9d, 0x3ff6a680d4b6b5ec, 0x3ff6a8775166126a, 0x3ff6aa6df9a0d6c8,
    0x3ff6ac64cd6ac90a, 0x3ff6ae5bccc7af86, 0x3ff6b052f7bb50e6, 0x3ff6b24a4e497429, 0x3ff6b441d075e0a1,
    0x3ff6b6397e445df5, 0x3ff6b83157b8b41e, 0x3ff6ba295cd6ab6a, 0x3ff6bc218da20c7a, 0x3ff6be19ea1ea046,
    0x3ff6c01272503016, 0x3ff6c20b263a8587, 0x3ff6c40405e16a8b, 0x3ff6c5fd1148a969, 0x3ff6c7f648740cb9,
    0x3ff6c9efab675f6a, 0x3ff6cbe93a266cbe, 0x3ff6cde2f4b5004b, 0x3ff6cfdcdb16e5fb, 0x3ff6d1d6ed4fea0d,
    0x3ff6d3d12b63d914, 0x3ff6d5cb95567ff7, 0x3ff6d7c62b2babf1, 0x3ff6d9c0ece72a93, 0x3ff6dbbbda8cc9c0,
    0x3ff6ddb6f42057b0, 0x3ff6dfb239a5a2f3, 0x3ff6e1adab207a67, 0x3ff6e3a94894ad43, 0x3ff6e5a512060b13,
    0x3ff6e7a1077863b4, 0x3ff6e99d28ef875a, 0x3ff6eb99766f468d, 0x3ff6ed95effb722a, 0x3ff6ef929597db64,
    0x3ff6f18f674853bf, 0x3ff6f38c6510ad16, 0x3ff6f5898ef4b99a, 0x3ff6f786e4f84bcd, 0x3ff6f984671f368a,
    0x3ff6fb82156d4cfe, 0x3ff6fd7fefe662ac, 0x3ff6ff7df68e4b6c, 0x3ff7017c2968db6b, 0x3ff7037a8879e729,
    0x3ff7057913c5437d, 0x3ff70777cb4ec594, 0x3ff70976af1a42ec, 0x3ff70b75bf2b915c, 0x3ff70d74fb868710,
    0x3ff70f74642efa85, 0x3ff71173f928c291, 0x3ff71373ba77b65f, 0x3ff71573a81fad6d, 0x3ff71773c2247f90,
    0x3ff71974088a04f2, 0x3ff71b747b541612, 0x3ff71d751a868bc4, 0x3ff71f75e6253f32, 0x3ff72176de3409db,
    0x3ff7237802b6c593, 0x3ff7257953b14c84, 0x3ff7277ad127792f, 0x3ff7297c7b1d2667, 0x3ff72b7e51962f56,
    0x3ff72d8054966f7e, 0x3ff72f828421c2b3, 0x3ff73184e03c0520, 0x3ff7338768e91346, 0x3ff7358a1e2cc9fc,
    0x3ff7378d000b066d, 0x3ff739900e87a61c, 0x3ff73b9349a686df, 0x3ff73d96b16b86e5, 0x3ff73f9a45da84b1,
    0x3ff7419e06f75f1a, 0x3ff743a1f4c5f551, 0x3ff745a60f4a26da, 0x3ff747aa5687d38f, 0x3ff749aeca82dba1,
    0x3ff74bb36b3f1f98, 0x3ff74db838c0804e, 0x3ff74fbd330adef7, 0x3ff751c25a221d1c, 0x3ff753c7ae0a1c9b,
    0x3ff755cd2ec6bfaa, 0x3ff757d2dc5be8d3, 0x3ff759d8b6cd7af8, 0x3ff75bdebe1f5950, 0x3ff75de4f2556769,
    0x3ff75feb53738927, 0x3ff761f1e17da2c4, 0x3ff763f89c7798d0, 0x3ff765ff84655034, 0x3ff76806994aae2c,
    0x3ff76a0ddb2b984c, 0x3ff76c154a0bf47d, 0x3ff76e1ce5efa903, 0x3ff77024aeda9c72, 0x3ff7722ca4d0b5ba,
    0x3ff77434c7d5dc1c, 0x3ff7763d17edf738, 0x3ff77845951ceefb, 0x3ff77a4e3f66abb0, 0x3ff77c5716cf15f4,
    0x3ff77e601b5a16bd, 0x3ff780694d0b9758, 0x3ff78272abe78169, 0x3ff7847c37f1bee8, 0x3ff78685f12e3a27,
    0x3ff7888fd7a0ddcc, 0x3ff78a99eb4d94d8, 0x3ff78ca42c384a9f, 0x3ff78eae9a64eacc, 0x3ff790b935d76165,
    0x3ff792c3fe939ac3, 0x3ff794cef49d8396, 0x3ff796da17f908e9, 0x3ff798e568aa181a, 0x3ff79af0e6b49ee0,
    0x3ff79cfc921c8b49, 0x3ff79f086ae5cbba, 0x3ff7a11471144eef, 0x3ff7a320a4ac03fa, 0x3ff7a52d05b0da48,
    0x3ff7a7399426c19b, 0x3ff7a9465011aa0b, 0x3ff7ab5339758409, 0x3ff7ad605056405d, 0x3ff7af6d94b7d027,
    0x3ff7b17b069e24de, 0x3ff7b388a60d3050, 0x3ff7b5967308e4a3, 0x3ff7b7a46d953455, 0x3ff7b9b295b6123a,
    0x3ff7bbc0eb6f7180, 0x3ff7bdcf6ec545aa, 0x3ff7bfde1fbb8295, 0x3ff7c1ecfe561c73, 0x3ff7c3fc0a9907d2,
    0x3ff7c60b44883993, 0x3ff7c81aac27a6f1, 0x3ff7ca2a417b4580, 0x3ff7cc3a04870b28, 0x3ff7ce49f54eee2d,
    0x3ff7d05a13d6e52a, 0x3ff7d26a6022e710, 0x3ff7d47ada36eb2a, 0x3ff7d68b8216e919, 0x3ff7d89c57c6d8d7,
    0x3ff7daad5b4ab2b8, 0x3ff7dcbe8ca66f64, 0x3ff7decfebde07de, 0x3ff7e0e178f57582, 0x3ff7e2f333f0b201,
    0x3ff7e5051cd3b766, 0x3ff7e71733a28014, 0x3ff7e929786106c7, 0x3ff7eb3beb134693, 0x3ff7ed4e8bbd3ae3,
    0x3ff7ef615a62df7a, 0x3ff7f17457083077, 0x3ff7f38781b12a4e, 0x3ff7f59ada61c9cc, 0x3ff7f7ae611e0c17,
    0x3ff7f9c215e9eeae, 0x3ff7fbd5f8c96f66, 0x3ff7fdea09c08c72, 0x3ff7fffe48d34457, 0x3ff80212b60595f7,
    0x3ff80427515b808b, 0x3ff8063c1ad903a4, 0x3ff8085112821f2e, 0x3ff80a66385ad36d, 0x3ff80c7b8c6720fb,
    0x3ff80e910eab08ce, 0x3ff810a6bf2a8c34, 0x3ff812bc9de9acd3, 0x3ff814d2aaec6cab, 0x3ff816e8e636ce15,
    0x3ff818ff4fccd3c0, 0x3ff81b15e7b280b9, 0x3ff81d2cadebd863, 0x3ff81f43a27cde79, 0x3ff8215ac5699711,
    0x3ff8237216b60699, 0x3ff82589966631da, 0x3ff827a1447e1df3, 0x3ff829b92101d060, 0x3ff82bd12bf54ef1,
    0x3ff82de9655c9fd6, 0x3ff83001cd3bc993, 0x3ff8321a6396d307, 0x3ff834332871c36a, 0x3ff8364c1bd0a24e,
    0x3ff838653db7779f, 0x3ff83a7e8e2a4ba1, 0x3ff83c980d2d26f1, 0x3ff83eb1bac41287, 0x3ff840cb96f317b4,
    0x3ff842e5a1be4023, 0x3ff844ffdb2995d7, 0x3ff8471a4339232f, 0x3ff84934d9f0f2e1, 0x3ff84b4f9f550fff,
    0x3ff84d6a936985f4, 0x3ff84f85b6326082, 0x3ff851a107b3abca, 0x3ff853bc87f17443, 0x3ff855d836efc6bd,
    0x3ff857f414b2b067, 0x3ff85a10213e3ec4, 0x3ff85c2c5c967fb5, 0x3ff85e48c6bf8171, 0x3ff860655fbd528d,
    0x3ff86282279401f7, 0x3ff8649f1e479ef5, 0x3ff866bc43dc392a, 0x3ff868d99855e08f, 0x3ff86af71bb8a57c,
    0x3ff86d14ce08989e, 0x3ff86f32af49cb03, 0x3ff87150bf804e0b, 0x3ff8736efeb03378, 0x3ff8758d6cdd8d61,
    0x3ff877ac0a0c6e38, 0x3ff879cad640e8cc, 0x3ff87be9d17f1044, 0x3ff87e08fbcaf822, 0x3ff880285528b444,
    0x3ff88247dd9c58df, 0x3ff884679529fa86, 0x3ff886877bd5ae23, 0x3ff888a791a388ff, 0x3ff88ac7d697a0b9,
    0x3ff88ce84ab60b4f, 0x3ff88f08ee02df15, 0x3ff89129c08232be, 0x3ff8934ac2381d54, 0x3ff8956bf328b63f,
    0x3ff8978d53581541, 0x3ff899aee2ca5273, 0x3ff89bd0a183864e, 0x3ff89df28f87c9a5, 0x3ff8a014acdb35a2,
    0x3ff8a236f981e3cd, 0x3ff8a459757fee0b, 0x3ff8a67c20d96e96, 0x3ff8a89efb928009, 0x3ff8aac205af3d57,
    0x3ff8ace53f33c1ce, 0x3ff8af08a8242919, 0x3ff8b12c40848f3b, 0x3ff8b35008591095, 0x3ff8b573ffa5c9e2,
    0x3ff8b798266ed839, 0x3ff8b9bc7cb8590d, 0x3ff8bbe102866a27, 0x3ff8be05b7dd29b2, 0x3ff8c02a9cc0b632,
    0x3ff8c24fb1352e86, 0x3ff8c474f53eb1e8, 0x3ff8c69a68e15fed, 0x3ff8c8c00c215887, 0x3ff8cae5df02bc04,
    0x3ff8cd0be189ab0a, 0x3ff8cf3213ba46a0, 0x3ff8d1587598b023, 0x3ff8d37f07290950, 0x3ff8d5a5c86f743d,
    0x3ff8d7ccb970135d, 0x3ff8d9f3da2f097f, 0x3ff8dc1b2ab079ca, 0x3ff8de42aaf887c7, 0x3ff8e06a5b0b5758,
    0x3ff8e2923aed0cb7, 0x3ff8e4ba4aa1cc81, 0x3ff8e6e28a2dbba9, 0x3ff8e90af994ff81, 0x3ff8eb3398dbbdb7,
    0x3ff8ed5c68061c54, 0x3ff8ef85671841bc, 0x3ff8f1ae961654b0, 0x3ff8f3d7f5047c4f, 0x3ff8f60183e6e012,
    0x3ff8f82b42c1a7cf, 0x3ff8fa553198fbb8, 0x3ff8fc7f5071045a, 0x3ff8fea99f4deaa1, 0x3ff900d41e33d7d1,
    0x3ff902fecd26f58f, 0x3ff90529ac2b6dda, 0x3ff90754bb456b0e, 0x3ff9097ffa7917e2, 0x3ff90bab69ca9f6c,
    0x3ff90dd7093e2d1b, 0x3ff91002d8d7ecbd, 0x3ff9122ed89c0a7e, 0x3ff9145b088eb2e4, 0x3ff9168768b412d0,
    0x3ff918b3f9105783, 0x3ff91ae0b9a7ae9b, 0x3ff91d0daa7e4610, 0x3ff91f3acb984c37, 0x3ff921681cf9efc3,
    0x3ff923959ea75fc4, 0x3ff925c350a4cba7, 0x3ff927f132f66333, 0x3ff92a1f45a05690, 0x3ff92c4d88a6d63f,
    0x3ff92e7bfc0e1323, 0x3ff930aa9fda3e74, 0x3ff932d9740f89d1, 0x3ff9350878b2272d, 0x3ff93737adc648dd,
    0x3ff9396713502192, 0x3ff93b96a953e45b, 0x3ff93dc66fd5c4a2, 0x3ff93ff666d9f630, 0x3ff942268e64ad2b,
    0x3ff94456e67a1e16, 0x3ff946876f1e7dd2, 0x3ff948b82856019b, 0x3ff94ae91224df0d, 0x3ff94d1a2c8f4c1e,
    0x3ff94f4b77997f27, 0x3ff9517cf347aeda, 0x3ff953ae9f9e1246, 0x3ff955e07ca0e0dd, 0x3ff958128a545266,
    0x3ff95a44c8bc9f0e, 0x3ff95c7737ddff5a, 0x3ff95ea9d7bcac2f, 0x3ff960dca85cdecf, 0x3ff9630fa9c2d0da,
    0x3ff96542dbf2bc4e, 0x3ff967763ef0db86, 0x3ff969a9d2c1693a, 0x3ff96bdd9768a084, 0x3ff96e118ceabcd7,
    0x3ff97045b34bfa05, 0x3ff9727a0a90943f, 0x3ff974ae92bcc816, 0x3ff976e34bd4d273, 0x3ff9791835dcf0a3,
    0x3ff97b4d50d9604e, 0x3ff97d829cce5f7c, 0x3ff97fb819c02c8f, 0x3ff981edc7b3064d, 0x3ff98423a6ab2bd5,
    0x3ff98659b6acdca7, 0x3ff9888ff7bc58a2, 0x3ff98ac669dde001, 0x3ff98cfd0d15b35d, 0x3ff98f33e16813b0,
    0x3ff9916ae6d94251, 0x3ff993a21d6d80f4, 0x3ff995d9852911ae, 0x3ff998111e1036f2, 0x3ff99a48e827338e,
    0x3ff99c80e3724ab5, 0x3ff99eb90ff5bff1, 0x3ff9a0f16db5d730, 0x3ff9a329fcb6d4be, 0x3ff9a562bcfcfd42,
    0x3ff9a79bae8c95c8, 0x3ff9a9d4d169e3b4, 0x3ff9ac0e25992ccd, 0x3ff9ae47ab1eb739, 0x3ff9b08161fec979,
    0x3ff9b2bb4a3daa71, 0x3ff9b4f563dfa161, 0x3ff9b72faee8f5e9, 0x3ff9b96a2b5df009, 0x3ff9bba4d942d81f,
    0x3ff9bddfb89bf6e9, 0x3ff9c01ac96d9580, 0x3ff9c2560bbbfd60, 0x3ff9c4917f8b7866, 0x3ff9c6cd24e050c8,
    0x3ff9c908fbbed121, 0x3ff9cb45042b4467, 0x3ff9cd813e29f5f2, 0x3ff9cfbda9bf3179, 0x3ff9d1fa46ef430e,
    0x3ff9d43715be772a, 0x3ff9d67416311aa0, 0x3ff9d8b1484b7aa2, 0x3ff9daeeac11e4c5, 0x3ff9dd2c4188a6fb,
    0x3ff9df6a08b40f94, 0x3ff9e1a801986d45, 0x3ff9e3e62c3a0f1d, 0x3ff9e624889d448d, 0x3ff9e86316c65d65,
    0x3ff9eaa1d6b9a9d6, 0x3ff9ece0c87b7a6f, 0x3ff9ef1fec102020, 0x3ff9f15f417bec36, 0x3ff9f39ec8c33062,
    0x3ff9f5de81ea3eb2, 0x3ff9f81e6cf56995, 0x3ff9fa5e89e903d9, 0x3ff9fc9ed8c960ac, 0x3ff9fedf599ad39d,
    0x3ffa01200c61b09a, 0x3ffa0360f1224bf2, 0x3ffa05a207e0fa53, 0x3ffa07e350a210ca, 0x3ffa0a24cb69e4c7,
    0x3ffa0c66783ccc19, 0x3ffa0ea8571f1ced, 0x3ffa10ea68152dd4, 0x3ffa132cab2355bc, 0x3ffa156f204debf5,
    0x3ffa17b1c7994830, 0x3ffa19f4a109c27b, 0x3ffa1c37aca3b348, 0x3ffa1e7aea6b7367, 0x3ffa20be5a655c0a,
    0x3ffa2301fc95c6c4, 0x3ffa2545d1010d86, 0x3ffa2789d7ab8aa3, 0x3ffa29ce109998cf, 0x3ffa2c127bcf931c,
    0x3ffa2e571951d502, 0x3ffa309be924ba55, 0x3ffa32e0eb4c9f4a, 0x3ffa35261fcde079, 0x3ffa376b86acdad9,
    0x3ffa39b11fedebc2, 0x3ffa3bf6eb9570ef, 0x3ffa3e3ce9a7c878, 0x3ffa40831a2950d8, 0x3ffa42c97d1e68ec,
    0x3ffa4510128b6ff1, 0x3ffa4756da74c583, 0x3ffa499dd4dec9a2, 0x3ffa4be501cddcad, 0x3ffa4e2c61465f66,
    0x3ffa5073f34cb2f0, 0x3ffa52bbb7e538cc, 0x3ffa5503af1452e0, 0x3ffa574bd8de6371, 0x3ffa59943547cd25,
    0x3ffa5bdcc454f307, 0x3ffa5e25860a387d, 0x3ffa606e7a6c0154, 0x3ffa62b7a17eb1b8, 0x3ffa6500fb46ae37,
    0x3ffa674a87c85bbf, 0x3ffa699447081fa2, 0x3ffa6bde390a5f91, 0x3ffa6e285dd3819f, 0x3ffa7072b567ec43,
    0x3ffa72bd3fcc0653, 0x3ffa7507fd043708, 0x3ffa7752ed14e5fb, 0x3ffa799e10027b29, 0x3ffa7be965d15ef0,
    0x3ffa7e34ee85fa0f, 0x3ffa8080aa24b5a6, 0x3ffa82cc98b1fb3a, 0x3ffa8518ba3234b0, 0x3ffa87650ea9cc4d,
    0x3ffa89b1961d2cbb, 0x3ffa8bfe5090c106, 0x3ffa8e4b3e08f499, 0x3ffa90985e8a3344, 0x3ffa92e5b218e937,
    0x3ffa953338b98307, 0x3ffa9780f2706da6, 0x3ffa99cedf42166e, 0x3ffa9c1cff32eb19, 0x3ffa9e6b524759c1,
    0x3ffaa0b9d883d0e6, 0x3ffaa30891ecbf66, 0x3ffaa5577e869486, 0x3ffaa7a69e55bfea, 0x3ffaa9f5f15eb19b,
    0x3ffaac4577a5da02, 0x3ffaae95312fa9ec, 0x3ffab0e51e009287, 0x3ffab3353e1d0565, 0x3ffab5859189747c,
    0x3ffab7d6184a5220, 0x3ffaba26d264110c, 0x3ffabc77bfdb245d, 0x3ffabec8e0b3ff90, 0x3ffac11a34f31687,
    0x3ffac36bbc9cdd87, 0x3ffac5bd77b5c936, 0x3ffac80f66424e9f, 0x3ffaca618846e330, 0x3ffaccb3ddc7fcb7,
    0x3ffacf0666ca1167, 0x3ffad159235197d6, 0x3ffad3ac136306fc, 0x3ffad5ff3702d636, 0x3ffad8528e357d43,
    0x3ffadaa618ff7445, 0x3ffadcf9d76533bf, 0x3ffadf4dc96b349b, 0x3ffae1a1ef15f025, 0x3ffae3f64869e00c,
    0x3ffae64ad56b7e60, 0x3ffae89f961f4598, 0x3ffaeaf48a89b08d, 0x3ffaed49b2af3a7a, 0x3ffaef9f0e945eff,
    0x3ffaf1f49e3d9a1f, 0x3ffaf44a61af6840, 0x3ffaf6a058ee462d, 0x3ffaf8f683feb114, 0x3ffafb4ce2e52685,
    0x3ffafda375a62474, 0x3ffafffa3c46293a, 0x3ffb025136c9b394, 0x3ffb04a8653542a2, 0x3ffb06ffc78d55e6,
    0x3ffb09575dd66d48, 0x3ffb0baf28150913, 0x3ffb0e07264da9f8, 0x3ffb105f5884d106, 0x3ffb12b7bebeffb8,
    0x3ffb15105900b7e6, 0x3ffb1769274e7bcf, 0x3ffb19c229acce18, 0x3ffb1c1b602031c6, 0x3ffb1e74caad2a44,
    0x3ffb20ce69583b61, 0x3ffb23283c25e951, 0x3ffb2582431ab8ab, 0x3ffb27dc7e3b2e6b, 0x3ffb2a36ed8bcff1,
    0x3ffb2c9191112300, 0x3ffb2eec68cfadc2, 0x3ffb314774cbf6c3, 0x3ffb33a2b50a84f5, 0x3ffb35fe298fdfad,
    0x3ffb3859d2608ea7, 0x3ffb3ab5af811a00, 0x3ffb3d11c0f60a3b, 0x3ffb3f6e06c3e840, 0x3ffb41ca80ef3d5d,
    0x3ffb44272f7c9343, 0x3ffb468412707405, 0x3ffb48e129cf6a20, 0x3ffb4b3e759e0071, 0x3ffb4d9bf5e0c23e,
    0x3ffb4ff9aa9c3b30, 0x3ffb525793d4f751, 0x3ffb54b5b18f8319, 0x3ffb571403d06b5b, 0x3ffb59728a9c3d55,
    0x3ffb5bd145f786a7, 0x3ffb5e3035e6d559, 0x3ffb608f5a6eb7d6, 0x3ffb62eeb393bcee, 0x3ffb654e415a73d6,
    0x3ffb67ae03c76c2a, 0x3ffb6a0dfadf35e8, 0x3ffb6c6e26a66177, 0x3ffb6ece87217fa1, 0x3ffb712f1c552196,
    0x3ffb738fe645d8e9, 0x3ffb75f0e4f83795, 0x3ffb78521870cffb, 0x3ffb7ab380b434df, 0x3ffb7d151dc6f96c,
    0x3ffb7f76efadb132, 0x3ffb81d8f66cf026, 0x3ffb843b32094aa4, 0x3ffb869da287556c, 0x3ffb890047eba5a5,
    0x3ffb8b63223ad0da, 0x3ffb8dc631796cfe, 0x3ffb902975ac1068, 0x3ffb928ceed751d6, 0x3ffb94f09cffc869,
    0x3ffb9754802a0bab, 0x3ffb99b8985ab38a, 0x3ffb9c1ce596585d, 0x3ffb9e8167e192dc, 0x3ffba0e61f40fc29,
    0x3ffba34b0bb92dca, 0x3ffba5b02d4ec1ab, 0x3ffba81584065220, 0x3ffbaa7b0fe479e1, 0x3ffbace0d0edd40c,
    0x3ffbaf46c726fc27, 0x3ffbb1acf2948e1f, 0x3ffbb413533b2643, 0x3ffbb679e91f614c, 0x3ffbb8e0b445dc58,
    0x3ffbbb47b4b334eb, 0x3ffbbdaeea6c08f0, 0x3ffbc0165574f6bb, 0x3ffbc27df5d29d00, 0x3ffbc4e5cb899adf,
    0x3ffbc74dd69e8fdc, 0x3ffbc9b617161be5, 0x3ffbcc1e8cf4df48, 0x3ffbce87383f7ac1, 0x3ffbd0f018fa8f6d,
    0x3ffbd3592f2abed3, 0x3ffbd5c27ad4aae0, 0x3ffbd82bfbfcf5e7, 0x3ffbda95b2a842a2, 0x3ffbdcff9edb3432,
    0x3ffbdf69c09a6e20, 0x3ffbe1d417ea945a, 0x3ffbe43ea4d04b36, 0x3ffbe6a967503772, 0x3ffbe9145f6efe30,
    0x3ffbeb7f8d3144fc, 0x3ffbedeaf09bb1c7, 0x3ffbf05689b2eaec, 0x3ffbf2c2587b9729, 0x3ffbf52e5cfa5da6,
    0x3ffbf79a9733e5f3, 0x3ffbfa07072cd804, 0x3ffbfc73ace9dc39, 0x3ffbfee0886f9b53, 0x3ffc014d99c2be80,
    0x3ffc03bae0e7ef53, 0x3ffc06285de3d7c7, 0x3ffc089610bb223d, 0x3ffc0b03f9727980, 0x3ffc0d72180e88c1,
    0x3ffc0fe06c93fb98, 0x3ffc124ef7077e06, 0x3ffc14bdb76dbc74, 0x3ffc172cadcb63b0, 0x3ffc199bda2520f2,
    0x3ffc1c0b3c7fa1d9, 0x3ffc1e7ad4df946e, 0x3ffc20eaa349a71c, 0x3ffc235aa7c288be, 0x3ffc25cae24ee890,
    0x3ffc283b52f37637, 0x3ffc2aabf9b4e1c5, 0x3ffc2d1cd697dbaf, 0x3ffc2f8de9a114d2, 0x3ffc31ff32d53e76,
    0x3ffc3470b2390a49, 0x3ffc36e267d12a62, 0x3ffc395453a25140, 0x3ffc3bc675b131cb, 0x3ffc3e38ce027f50,
    0x3ffc40ab5c9aed89, 0x3ffc431e217f3095, 0x3ffc45911cb3fcfd, 0x3ffc48044e3e07b0, 0x3ffc4a77b6220609,
    0x3ffc4ceb5464adc8, 0x3ffc4f5f290ab517, 0x3ffc51d33418d28a, 0x3ffc54477593bd1c, 0x3ffc56bbed802c30,
    0x3ffc59309be2d792, 0x3ffc5ba580c07778, 0x3ffc5e1a9c1dc47f, 0x3ffc608fedff77ae, 0x3ffc6305766a4a74,
    0x3ffc657b3562f6a9, 0x3ffc67f12aee368d, 0x3ffc6a675710c4cc, 0x3ffc6cddb9cf5c77, 0x3ffc6f54532eb909,
    0x3ffc71cb23339668, 0x3ffc744229e2b0e1, 0x3ffc76b96740c52b, 0x3ffc7930db529065, 0x3ffc7ba8861cd01a,
    0x3ffc7e2067a44239, 0x3ffc80987feda51f, 0x3ffc8310cefdb791, 0x3ffc858954d938bc, 0x3ffc88021184e837,
    0x3ffc8a7b05058602, 0x3ffc8cf42f5fd289, 0x3ffc8f6d90988e9c, 0x3ffc91e728b47b79, 0x3ffc9460f7b85ac7,
    0x3ffc96dafda8ee95, 0x3ffc99553a8af95b, 0x3ffc9bcfae633dfe, 0x3ffc9e4a59367fca, 0x3ffca0c53b098273,
    0x3ffca34053e10a1b, 0x3ffca5bba3c1db4b, 0x3ffca8372ab0baf6, 0x3ffcaab2e8b26e78, 0x3ffcad2eddcbbb9a,
    0x3ffcafab0a01688c, 0x3ffcb2276d583be7, 0x3ffcb4a407d4fcb3, 0x3ffcb720d97c725c, 0x3ffcb99de25364bb,
    0x3ffcbc1b225e9c14, 0x3ffcbe9899a2e114, 0x3ffcc1164824fcd0, 0x3ffcc3942de9b8ca, 0x3ffcc6124af5deee,
    0x3ffcc8909f4e3990, 0x3ffccb0f2af79372, 0x3ffccd8dedf6b7bd, 0x3ffcd00ce8507204, 0x3ffcd28c1a098e48,
    0x3ffcd50b8326d8f2, 0x3ffcd78b23ad1ed5, 0x3ffcda0afba12d30, 0x3ffcdc8b0b07d1aa, 0x3ffcdf0b51e5da58,
    0x3ffce18bd04015b7, 0x3ffce40c861b52b1, 0x3ffce68d737c6096, 0x3ffce90e98680f28, 0x3ffceb8ff4e32e8c,
    0x3ffcee1188f28f58, 0x3ffcf093549b0289, 0x3ffcf31557e15988, 0x3ffcf59792ca6629, 0x3ffcf81a055afaab,
    0x3ffcfa9caf97e9b8, 0x3ffcfd1f91860666, 0x3ffcffa2ab2a2432, 0x3ffd0225fc891709, 0x3ffd04a985a7b341,
    0x3ffd072d468acd9b, 0x3ffd09b13f373b42, 0x3ffd0c356fb1d1ce, 0x3ffd0eb9d7ff6743, 0x3ffd113e7824d20f,
    0x3ffd13c35026e90b, 0x3ffd1648600a837b, 0x3ffd18cda7d4790f, 0x3ffd1b532789a1e4, 0x3ffd1dd8df2ed681,
    0x3ffd205ecec8efd8, 0x3ffd22e4f65cc746, 0x3ffd256b55ef3696, 0x3ffd27f1ed8517fd, 0x3ffd2a78bd23461a,
    0x3ffd2cffc4ce9bfe, 0x3ffd2f87048bf51e, 0x3ffd320e7c602d5e, 0x3ffd34962c50210e, 0x3ffd371e1460aced,
    0x3ffd39a63496ae1e, 0x3ffd3c2e8cf70237, 0x3ffd3eb71d868736, 0x3ffd413fe64a1b88, 0x3ffd43c8e7469e02,
    0x3ffd46522080edeb, 0x3ffd48db91fdeaf0, 0x3ffd4b653bc2752c, 0x3ffd4def1dd36d29, 0x3ffd50793835b3da,
    0x3ffd53038aee2a9f, 0x3ffd558e1601b344, 0x3ffd5818d9753003, 0x3ffd5aa3d54d8381, 0x3ffd5d2f098f90cf,
    0x3ffd5fba76403b6b, 0x3ffd62461b64673f, 0x3ffd64d1f900f8a4, 0x3ffd675e0f1ad45a, 0x3ffd69ea5db6df94,
    0x3ffd6c76e4d9ffed, 0x3ffd6f03a4891b6e, 0x3ffd71909cc9188f, 0x3ffd741dcd9ede30, 0x3ffd76ab370f53a1,
    0x3ffd7938d91f609f, 0x3ffd7bc6b3d3ed53, 0x3ffd7e54c731e251, 0x3ffd80e3133e289e, 0x3ffd837197fda9a8,
    0x3ffd860055754f4c, 0x3ffd888f4baa03d3, 0x3ffd8b1e7aa0b1f5, 0x3ffd8dade25e44d5, 0x3ffd903d82e7a803,
    0x3ffd92cd5c41c77f, 0x3ffd955d6e718fb2, 0x3ffd97edb97bed76, 0x3ffd9a7e3d65ce10, 0x3ffd9d0efa341f33,
    0x3ffd9f9fefebceff, 0x3ffda2311e91cc02, 0x3ffda4c2862b0536, 0x3ffda75426bc6a05, 0x3ffda9e6004aea45,
    0x3ffdac7812db7638, 0x3ffdaf0a5e72fe91, 0x3ffdb19ce316746e, 0x3ffdb42fa0cac95a, 0x3ffdb6c29794ef50,
    0x3ffdb955c779d8b8, 0x3ffdbbe9307e7867, 0x3ffdbe7cd2a7c1a1, 0x3ffdc110adfaa815, 0x3ffdc3a4c27c1fe4,
    0x3ffdc63910311d9a, 0x3ffdc8cd971e9631, 0x3ffdcb6257497f13, 0x3ffdcdf750b6ce17, 0x3ffdd08c836b797f,
    0x3ffdd321ef6c7800, 0x3ffdd5b794bec0bc, 0x3ffdd84d73674b3f, 0x3ffddae38b6b0f89, 0x3ffddd79dccf0603,
    0x3ffde01067982789, 0x3ffde2a72bcb6d61, 0x3ffde53e296dd143, 0x3ffde7d560844d54, 0x3ffdea6cd113dc26,
    0x3ffded047b2178bb, 0x3ffdef9c5eb21e83, 0x3ffdf2347bcac95e, 0x3ffdf4ccd2707596, 0x3ffdf76562a81feb,
    0x3ffdf9fe2c76c585, 0x3ffdfc972fe163fd, 0x3ffdff306cecf95b, 0x3ffe01c9e39e8418, 0x3ffe046393fb0315,
    0x3ffe06fd7e0775aa, 0x3ffe0997a1c8db99, 0x3ffe0c31ff443512, 0x3ffe0ecc967e82b9, 0x3ffe1167677cc59c,
    0x3ffe14027243ff3b, 0x3ffe169db6d93183, 0x3ffe193935415ed1, 0x3ffe1bd4ed8189f2, 0x3ffe1e70df9eb621,
    0x3ffe210d0b9de709, 0x3ffe23a9718420c3, 0x3ffe2646115667d9, 0x3ffe28e2eb19c142, 0x3ffe2b7ffed33266,
    0x3ffe2e1d4c87c11e, 0x3ffe30bad43c73ae, 0x3ffe335895f650cf, 0x3ffe35f691ba5fa4, 0x3ffe3894c78da7c2,
    0x3ffe3b333775312f, 0x3ffe3dd1e176045e, 0x3ffe4070c5952a35, 0x3ffe430fe3d7ac06, 0x3ffe45af3c429394,
    0x3ffe484ecedaeb14, 0x3ffe4aee9ba5bd26, 0x3ffe4d8ea2a814df, 0x3ffe502ee3e6fdc2, 0x3ffe52cf5f6783c0,
    0x3ffe5570152eb33c, 0x3ffe581105419909, 0x3ffe5ab22fa54269, 0x3ffe5d53945ebd0f, 0x3ffe5ff53373171e,
    0x3ffe62970ce75f28, 0x3ffe653920c0a430, 0x3ffe67db6f03f5ab, 0x3ffe6a7df7b6637a, 0x3ffe6d20badcfdf3,
    0x3ffe6fc3b87cd5d9, 0x3ffe7266f09afc62, 0x3ffe750a633c8332, 0x3ffe77ae10667c5d, 0x3ffe7a51f81dfa6b,
    0x3ffe7cf61a681052, 0x3ffe7f9a7749d178, 0x3ffe823f0ec851b6, 0x3ffe84e3e0e8a554, 0x3ffe8788edafe10a,
    0x3ffe8a2e35231a01, 0x3ffe8cd3b74765d6, 0x3ffe8f797421da93, 0x3ffe921f6bb78eb2, 0x3ffe94c59e0d9924,
    0x3ffe976c0b291144, 0x3ffe9a12b30f0ee0, 0x3ffe9cb995c4aa3b, 0x3ffe9f60b34efc02, 0x3ffea2080bb31d5a,
    0x3ffea4af9ef627d4, 0x3ffea7576d1d3575, 0x3ffea9ff762d60b2, 0x3ffeaca7ba2bc471, 0x3ffeaf50391d7c09,
    0x3ffeb1f8f307a346, 0x3ffeb4a1e7ef5660, 0x3ffeb74b17d9b203, 0x3ffeb9f482cbd34b, 0x3ffebc9e28cad7ca,
    0x3ffebf4809dbdd7c, 0x3ffec1f2260402d5, 0x3ffec49c7d4866b8, 0x3ffec7470fae2879, 0x3ffec9f1dd3a67df,
    0x3ffecc9ce5f24521, 0x3ffecf4829dae0eb, 0x3ffed1f3a8f95c56, 0x3ffed49f6352d8ef, 0x3ffed74b58ec78b7,
    0x3ffed9f789cb5e20, 0x3ffedca3f5f4ac0a, 0x3ffedf509d6d85cb, 0x3ffee1fd803b0f2a, 0x3ffee4aa9e626c5f,
    0x3ffee757f7e8c217, 0x3ffeea058cd3356e, 0x3ffeecb35d26ebf2, 0x3ffeef6168e90ba5, 0x3ffef20fb01ebafb,
    0x3ffef4be32cd20da, 0x3ffef76cf0f9649a, 0x3ffefa1beaa8ae04, 0x3ffefccb1fe02556, 0x3ffeff7a90a4f33f,
    0x3fff022a3cfc40e1, 0x3fff04da24eb37d0, 0x3fff078a48770213, 0x3fff0a3aa7a4ca23, 0x3fff0ceb4279baea,
    0x3fff0f9c18faffca, 0x3fff124d2b2dc491, 0x3fff14fe79173584, 0x3fff17b002bc7f5a, 0x3fff1a61c822cf3c,
    0x3fff1d13c94f52c7, 0x3fff1fc606473809, 0x3fff22787f0fad85, 0x3fff252b33ade22f, 0x3fff27de24270571,
    0x3fff2a9150804723, 0x3fff2d44b8bed796, 0x3fff2ff85ce7e78a, 0x3fff32ac3d00a832, 0x3fff3560590e4b38,
    0x3fff3814b11602b5, 0x3fff3ac9451d0138, 0x3fff3d7e152879c2, 0x3fff4033213d9fc8, 0x3fff42e86961a731,
    0x3fff459ded99c45a, 0x3fff4853adeb2c11, 0x3fff4b09aa5b1398, 0x3fff4dbfe2eeb0a6, 0x3fff507657ab3963,
    0x3fff532d0895e46e, 0x3fff55e3f5b3e8d8, 0x3fff589b1f0a7e23, 0x3fff5b52849edc4a, 0x3fff5e0a26763bb8,
    0x3fff60c20495d54d, 0x3fff637a1f02e25c, 0x3fff663275c29cab, 0x3fff68eb08da3e7a, 0x3fff6ba3d84f0275,
    0x3fff6e5ce42623c1, 0x3fff71162c64ddf3, 0x3fff73cfb1106d1b, 0x3fff7689722e0db5, 0x3fff79436fc2fcb6,
    0x3fff7bfda9d47787, 0x3fff7eb82067bc04, 0x3fff8172d382087c, 0x3fff842dc3289bb5, 0x3fff86e8ef60b4ea,
    0x3fff89a4582f93c7, 0x3fff8c5ffd9a786e, 0x3fff8f1bdfa6a377, 0x3fff91d7fe5955eb, 0x3fff949459b7d14b,
    0x3fff9750f1c7578c, 0x3fff9a0dc68d2b16, 0x3fff9ccad80e8ec8, 0x3fff9f882650c5f2, 0x3fffa245b159145c,
    0x3fffa503792cbe42, 0x3fffa7c17dd10856, 0x3fffaa7fbf4b37bd, 0x3fffad3e3da09211, 0x3fffaffcf8d65d61,
    0x3fffb2bbf0f1e031, 0x3fffb57b25f8617d, 0x3fffb83a97ef28b2, 0x3fffbafa46db7db4, 0x3fffbdba32c2a8db,
    0x3fffc07a5ba9f2f6, 0x3fffc33ac196a548, 0x3fffc5fb648e098a, 0x3fffc8bc449569e9, 0x3fffcb7d61b21108,
    0x3fffce3ebbe94a01, 0x3fffd10053406061, 0x3fffd3c227bca02c, 0x3fffd684396355da, 0x3fffd9468839ce5a,
    0x3fffdc0914455712, 0x3fffdecbdd8b3dd8, 0x3fffe18ee410d0ff, 0x3fffe45227db5f4b, 0x3fffe715a8f037f6,
    0x3fffe9d96754aab1, 0x3fffec9d630e07a4, 0x3fffef619c219f69, 0x3ffff2261294c314, 0x3ffff4eac66cc42c,
    0x3ffff7afb7aef4b0, 0x3ffffa74e660a715, 0x3ffffd3a52872e44, 0x3ffffffffc27dd9e};

// rcpps's reciprocal of X = 1.m * 2^-10, rounded to an integer by roundpd
// (nearest, even): on AMD Zen 2 and Zen 3 (the fixtures' machines; their
// rcpps agree on every float) it is a function of the top 12 bits of m,
// 1024 less the number of these breakpoints (m >> 11) at or below them.
// The exact reciprocal rounds otherwise on 628,533 of the 2^23 mantissas,
// so rcpps is carried as this table. FT4's only base, 10 (m = 0x200000),
// gives 819 (rcpps 819.0 there; 1/X = 819.2, so any x86's rcpps,
// within 1.5 * 2^-12, gives 819 too).
const uint16_t powf4_rcp_breaks[512] = {
    2, 6, 11, 14, 18, 22, 27, 30, 35, 38, 43, 47, 51, 54, 59, 63,
    68, 71, 76, 80, 84, 88, 93, 96, 101, 104, 109, 113, 118, 121, 126, 130,
    135, 139, 143, 147, 152, 156, 160, 164, 169, 173, 178, 182, 187, 190, 196, 199,
    204, 208, 213, 217, 222, 226, 231, 235, 239, 244, 248, 252, 257, 261, 267, 271,
    275, 280, 285, 288, 294, 298, 303, 307, 313, 317, 322, 326, 331, 335, 341, 345,
    350, 354, 360, 364, 369, 373, 379, 383, 388, 392, 398, 402, 407, 411, 417, 421,
    427, 431, 437, 441, 446, 451, 456, 460, 466, 470, 476, 480, 486, 490, 495, 501,
    506, 510, 516, 521, 527, 530, 537, 541, 546, 551, 557, 561, 567, 572, 578, 583,
    589, 592, 599, 603, 609, 614, 620, 624, 631, 635, 641, 646, 652, 656, 663, 667,
    674, 679, 685, 689, 695, 700, 706, 711, 717, 722, 729, 733, 740, 745, 751, 756,
    762, 767, 773, 778, 784, 790, 796, 800, 807, 812, 819, 824, 831, 835, 842, 847,
    854, 859, 865, 870, 877, 882, 889, 894, 901, 906, 912, 918, 924, 930, 937, 942,
    949, 954, 961, 966, 973, 979, 985, 991, 998, 1003, 1011, 1016, 1023, 1028, 1035, 1040,
    1048, 1054, 1061, 1066, 1073, 1079, 1086, 1092, 1099, 1104, 1112, 1118, 1125, 1130, 1138, 1144,
    1151, 1157, 1164, 1170, 1177, 1183, 1190, 1196, 1204, 1209, 1217, 1223, 1230, 1237, 1244, 1250,
    1257, 1263, 1272, 1278, 1285, 1292, 1299, 1306, 1313, 1319, 1327, 1333, 1341, 1347, 1356, 1362,
    1370, 1376, 1384, 1390, 1398, 1404, 1412, 1419, 1427, 1434, 1442, 1448, 1457, 1463, 1471, 1477,
    1485, 1492, 1500, 1507, 1516, 1523, 1531, 1537, 1546, 1552, 1561, 1568, 1577, 1584, 1592, 1599,
    1607, 1614, 1623, 1630, 1639, 1646, 1654, 1661, 1670, 1677, 1686, 1693, 1702, 1709, 1718, 1725,
    1734, 1741, 1751, 1758, 1767, 1774, 1783, 1790, 1800, 1807, 1817, 1824, 1833, 1840, 1850, 1857,
    1867, 1874, 1884, 1891, 1901, 1908, 1918, 1925, 1935, 1943, 1953, 1961, 1970, 1978, 1988, 1996,
    2005, 2013, 2024, 2031, 2041, 2049, 2060, 2068, 2078, 2085, 2095, 2103, 2114, 2123, 2133, 2141,
    2151, 2159, 2169, 2178, 2189, 2197, 2207, 2215, 2226, 2235, 2246, 2254, 2264, 2273, 2284, 2293,
    2303, 2312, 2323, 2331, 2342, 2351, 2362, 2372, 2383, 2392, 2403, 2411, 2422, 2431, 2443, 2452,
    2463, 2473, 2484, 2493, 2505, 2515, 2526, 2534, 2547, 2555, 2568, 2577, 2589, 2599, 2610, 2620,
    2632, 2641, 2653, 2663, 2676, 2685, 2697, 2706, 2719, 2729, 2741, 2751, 2764, 2773, 2786, 2796,
    2809, 2819, 2832, 2842, 2855, 2865, 2879, 2889, 2902, 2912, 2925, 2935, 2948, 2959, 2972, 2982,
    2996, 3007, 3020, 3031, 3044, 3055, 3068, 3079, 3093, 3104, 3119, 3128, 3143, 3154, 3168, 3179,
    3194, 3203, 3218, 3229, 3243, 3255, 3271, 3282, 3296, 3307, 3322, 3333, 3349, 3360, 3375, 3386,
    3402, 3413, 3429, 3440, 3455, 3468, 3482, 3495, 3511, 3523, 3538, 3549, 3565, 3578, 3593, 3606,
    3622, 3634, 3650, 3663, 3679, 3693, 3709, 3721, 3737, 3750, 3768, 3780, 3796, 3810, 3826, 3839,
    3856, 3869, 3887, 3899, 3917, 3931, 3947, 3961, 3979, 3993, 4009, 4023, 4041, 4055, 4073, 4087};

}  // namespace

void powf4_glibc(const float* x, const float* y, float* out) {
    // __svml_spow_data
    const uint64_t ExpMask = 0x000fffffffffffffull, Two10 = 0x3f50000000000000ull, LOWKBITS = 0x7ff;
    const double poly_coeff_3 = asdouble(0x3fd55555555543c5ull), poly_coeff_4 = asdouble(0xbfdffffffffff81full),
                 One = asdouble(0x3ff0000000000000ull), Threshold = asdouble(0x4086a00000000000ull),
                 Bias = asdouble(0x408ff80000000000ull), Bias1 = asdouble(0x408ff00000000000ull),
                 L2 = asdouble(0x3fe62e42fefa39efull), dbInvLn2 = asdouble(0x40a71547652b82feull),
                 dbShifter = asdouble(0x4338000000000000ull), dbHALF = asdouble(0x3fe0000000000000ull),
                 dbC1 = asdouble(0x3f362f3904051fa1ull);
    const int32_t NMINNORM = static_cast<int32_t>(0x80800000u), NMAXVAL = static_cast<int32_t>(0xfeffffffu),
                  INF = 0x7f800000, iDomainRange = 0x4059fe36;
    for (int l = 0; l < 4; ++l) {
        const uint32_t ix = fbits(x[l]), iy = fbits(y[l]);
        // log: the mantissa at 2^-10, its rounded reciprocal R, the
        // biased exponent; r = X*R - 1, ln|x| = (exponent - bias)*ln2 +
        // table(R) + (r + c4*r^2 + c3*r^3).
        const uint64_t xb = asuint64(static_cast<double>(x[l]));
        const double X = asdouble((xb & ExpMask) | Two10);
        const int32_t e = static_cast<int32_t>(static_cast<uint32_t>((xb >> 20) >> 32));
        const int key = static_cast<int>((ix & 0x7fffffu) >> 11);
        const double R = 1024.0 - static_cast<double>(std::upper_bound(powf4_rcp_breaks, powf4_rcp_breaks + 512, key) -
                                                     powf4_rcp_breaks);
        const double r = X * R - One;
        const double bias = Threshold < R ? Bias : Bias1;
        const double kl = (static_cast<double>(e) - bias) * L2;
        const double r2 = r * r;
        double p = poly_coeff_3 * r;
        p = p + poly_coeff_4;
        p = p * r2;
        const double s = r + p;
        const uint64_t ti = (asuint64(R) >> 40) - 0x408000u;
        const double t = asdouble(powf4_log_tab[ti / 8]) + s;
        const double ylnx = (kl + t) * static_cast<double>(y[l]);
        // exp: d = y*ln|x| * 2048/ln2 rounded down to N*2048 + j, the
        // fraction by 2^(1/2048)-1; 2^(j/2048)*(1 + fraction), N added to
        // the exponent.
        const double d = dbInvLn2 * ylnx;
        const double sh = (d - dbHALF) + dbShifter;
        const double fl = sh - dbShifter;
        const double fr = (d - fl) * dbC1;
        const uint64_t shb = asuint64(sh);
        const double T = asdouble(powf4_exp_tab[shb & LOWKBITS]);
        const double res = T + fr * T;
        const float v = static_cast<float>(asdouble(asuint64(res) + ((shb >> 11) << 52)));
        // Lanes out of range by the scalar powf: x not in [2^-126,
        // FLT_MAX), y infinite or NaN, |y*ln|x|| over the domain range
        // (by the high word of the double).
        const int32_t xm = static_cast<int32_t>(ix - static_cast<uint32_t>(NMINNORM));
        const int32_t ay = static_cast<int32_t>(iy & 0x7fffffffu);
        const int32_t hi = static_cast<int32_t>(static_cast<uint32_t>(asuint64(ylnx) >> 32) & 0x7fffffffu);
        const bool special = xm > NMAXVAL || xm == NMAXVAL || ay > INF || ay == INF || hi > iDomainRange;
        out[l] = special ? powf_glibc(x[l], y[l]) : v;
    }
}

}  // namespace ft8
