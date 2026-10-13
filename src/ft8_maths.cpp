// glibc 2.39's erff, log10f, sinf, cosf, sincosf, expf, logf, powf,
// hypotf and x86-64 libmvec cosf, sinf and log10f (ft8_maths.h), carried. Built without contraction (CMakeLists.txt), as
// glibc's x86-64 build of them has none: every multiply and add rounded as
// written.
//
// s_erff.c and e_log10f.c: Sun Microsystems (fdlibm; freely
// distributable). s_sinf.c, s_cosf.c, s_sincosf.c, e_expf.c, e_logf.c,
// e_powf.c and their data: ARM Ltd. and the Free Software Foundation (GNU
// LGPL v2.1 or later); e_hypotf.c: the Free Software Foundation (GNU LGPL
// v2.1 or later). svml_s_cosf4_core_sse4.S, svml_s_sinf4_core_sse4.S,
// svml_s_log10f4_core_sse4.S: the GNU C Library, Free Software Foundation
// (GNU LGPL v2.1 or later).
#include "ft8_maths.h"

#include <cmath>
#include <cstdint>
#include <cstring>

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

}  // namespace ft8
