// jt9 --fast-subtract's smoothing (ft_fast_subtract.h).
#include "ft_fast_subtract.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace ft8 {

using cfloat = std::complex<float>;

FastSubtract::FastSubtract(int nframe, int nfilt, int d) : nframe_(nframe), d_(d), m_(nframe / d) {
    if (nframe % d != 0 || d % 2 != 0) throw std::invalid_argument("FastSubtract: whole, even blocks");
    const int q = nfilt / 2 / d + 1;  // the kernel's reach, in blocks
    nc_ = 1;
    while (nc_ < m_ + 2 * q + 1) nc_ *= 2;  // the convolution's FFT: long enough that nothing wraps
    std::vector<double> win(static_cast<size_t>(nfilt + 1));
    double sumw = 0;
    for (int j = -nfilt / 2; j <= nfilt / 2; ++j) {
        const double a = std::cos(3.14159265358979323846 * j / nfilt);
        win[static_cast<size_t>(j + nfilt / 2)] = a * a;
        sumw += a * a;
    }
    // K[k]: the window at lag k*d, d of them per decimated sample, at k mod nc; its spectrum, scaled by 1/nc for
    // the unnormalized inverse
    kspec_.assign(static_cast<size_t>(nc_), cfloat(0.0f, 0.0f));
    for (int k = -q; k <= q; ++k) {
        const int j = k * d;
        if (j >= -nfilt / 2 && j <= nfilt / 2)
            kspec_[static_cast<size_t>((k + nc_) % nc_)] =
                cfloat(static_cast<float>(d * win[static_cast<size_t>(j + nfilt / 2)] / sumw / nc_), 0.0f);
    }
    fst4w::fft_c2c_inplace(kspec_.data(), nc_, -1);
    s0_.resize(static_cast<size_t>(m_ + 1));
    s1_.resize(static_cast<size_t>(m_ + 1));
    y_.assign(static_cast<size_t>(nc_), cfloat(0.0f, 0.0f));
}

size_t FastSubtract::bytes() const {
    return (s0_.capacity() + s1_.capacity() + kspec_.capacity() + y_.capacity()) * sizeof(cfloat);
}

void FastSubtract::subtract(const float* x, const cfloat* cref, int n0, int n1, const float* gain, float* out) {
    const int d = d_, m = m_;
    // The triangle's halves are linear in a sample's place k in a d-sample half-block (half-block n: the samples
    // n*d - d/2 + k, k = 0..d-1), so each half-block's sums of c and k*c give it, each sample read once:
    // b[j] = (S1[j] + S0[j]/2 + (d - 1/2)*S0[j+1] - S1[j+1])/d^2, centred on t_c = j*d + (d-1)/2.
    cfloat* s0 = s0_.data();
    cfloat* s1 = s1_.data();
    for (int n = 0; n <= m; ++n) {
        const int base = n * d - d / 2;
        const int lo = std::max(base, n0), hi = std::min(base + d, n1);
        float s0r = 0, s0i = 0, s1r = 0, s1i = 0;
        for (int i = lo; i < hi; ++i) {
            const float k = static_cast<float>(i - base);
            const float cr = x[i] * cref[i].real(), ci = -x[i] * cref[i].imag();
            s0r += cr;
            s0i += ci;
            s1r += k * cr;
            s1i += k * ci;
        }
        s0[n] = cfloat(s0r, s0i);
        s1[n] = cfloat(s1r, s1i);
    }
    cfloat* y = y_.data();
    const float inv_d2 = 1.0f / static_cast<float>(d * d), dh = static_cast<float>(d) - 0.5f;
    for (int j = 0; j < m; ++j) y[j] = (s1[j] + 0.5f * s0[j] + dh * s0[j + 1] - s1[j + 1]) * inv_d2;
    std::fill(y + m, y + nc_, cfloat(0.0f, 0.0f));
    // y at the block centres t_j = j*d + (d-1)/2 - 1 (the cshift's -1): b convolved with K
    fst4w::fft_c2c_inplace(y, nc_, -1);
    const cfloat* ks = kspec_.data();
    for (int k = 0; k < nc_; ++k) y[k] = y[k] * ks[k];
    fst4w::fft_c2c_inplace(y, nc_, 1);
    // Each sample: y linear between the centres, a segment at a time (the samples between two centres:
    // u = 2*(i - t_0) in [2*d*j, 2*d*(j + 1)), t_0 = (d-1)/2 - 1), f stepping 1/d a sample; before the first centre
    // y[0], past the last y[m-1]
    for (int i = n0; i < n1;) {
        const int u = 2 * i - (d - 3);
        int j = u >= 0 ? u / (2 * d) : -((-u + 2 * d - 1) / (2 * d));
        float f0 = static_cast<float>(u - 2 * d * j) / static_cast<float>(2 * d), df = 1.0f / static_cast<float>(d);
        int iend;
        if (j < 0) {
            j = 0, f0 = 0, df = 0;
            iend = std::min(n1, (d - 2) / 2);  // the first i with u >= 0
        } else if (j >= m - 1) {
            j = m - 2, f0 = 1, df = 0;
            iend = n1;
        } else {
            iend = std::min(n1, i + (2 * d * (j + 1) - u + 1) / 2);
        }
        const float y0r = y[j].real(), y0i = y[j].imag();
        const float dyr = y[j + 1].real() - y0r, dyi = y[j + 1].imag() - y0i;
        const int len = iend - i;
        if (gain) {
            for (int k = 0; k < len; ++k) {
                const float f = f0 + static_cast<float>(k) * df, g = gain[i + k];
                const float yr = (y0r + f * dyr) * g, yi = (y0i + f * dyi) * g;
                out[i + k] = x[i + k] - 2.0f * (yr * cref[i + k].real() - yi * cref[i + k].imag());
            }
        } else {
            for (int k = 0; k < len; ++k) {
                const float f = f0 + static_cast<float>(k) * df;
                const float yr = y0r + f * dyr, yi = y0i + f * dyi;
                out[i + k] = x[i + k] - 2.0f * (yr * cref[i + k].real() - yi * cref[i + k].imag());
            }
        }
        i = iend;
    }
}

}  // namespace ft8
