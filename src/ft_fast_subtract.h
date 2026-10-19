// Ours: jt9 --fast-subtract's smoothing, shared by FT8's sqf (ft8_frontend.cpp) and FT4's subtractft4
// (ft4_wave.cpp, FT2's too). Both reconstruct a decoded signal from the received one: c = x*conj(cref) over the frame
// (zero outside [n0, n1)), smoothed with the normalized cos^2 window w of nfilt+1 taps, y[i] = sum_j w[j] c[i+1-j]
// (the +1 is the references' cshift), times a gain at the frame's ends (FT8's end correction; FT4 has none), and
// subtracted: out[i] = x[i] - 2*Re(y[i]*cref[i]). The references do the smoothing with two FFTs of the whole period.
// y's band is a few Hz wide, so here c is taken down by d through a 2d-sample triangle (as two cascaded d-sample block
// means: aliases ~80 dB down, where one block mean's ~40 dB moved decodes on crowded bands), convolved with w
// sampled every d samples (through a small FFT), and taken back to each sample linearly between the block centres.
// Not the references' to the bit; the same decodes (scripts/ft-fast-subtract-check.sh).
#pragma once

#include <complex>
#include <cstddef>
#include <vector>

#include "fst4w_frontend.h"  // FftwVec

namespace ft8 {

class FastSubtract {
public:
    // nframe a multiple of d, d even; the window w[j] = cos^2(pi*j/nfilt), j = -nfilt/2..nfilt/2, normalized.
    FastSubtract(int nframe, int nfilt, int d);

    // out[i] = x[i] - 2*Re(gain[i]*y[i]*cref[i]) for i in [n0, n1) (gain nullptr: 1); out may be x. Nothing else
    // of out is written.
    void subtract(const float* x, const std::complex<float>* cref, int n0, int n1, const float* gain, float* out);

    size_t bytes() const;

private:
    int nframe_, d_, m_, nc_;
    std::vector<std::complex<float>> s0_, s1_;  // each half-block's sums of c and k*c
    fst4w::FftwVec<std::complex<float>> kspec_, y_;
};

}  // namespace ft8
