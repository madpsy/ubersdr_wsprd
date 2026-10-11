// A forward complex FFT of fixed size over FFTW, transforming in place in a
// buffer it owns. Planning is serialised because FFTW's planner is not
// thread-safe; executing a plan is.
#pragma once

#include <complex>
#include <cstddef>
#include <mutex>

namespace skimmer {

using C32 = std::complex<float>;

// The lock every FFTW plan of the process is made and destroyed under, for
// code making its own plans.
std::mutex& fft_planner_lock();

// The wisdom file (fft_wisdom.cpp): FFTW's wisdom loaded from it, once a
// process, before each plan; what the planner knows saved to it after a
// measured plan. A program with no file of its own defines both to do
// nothing.
void fft_load_wisdom();
void fft_save_wisdom();

class Fft {
public:
    // quick: without wisdom for n, FFTW's estimate rather than its
    // measured plan, made at once: for transforms run seldom, whose
    // measured plans take minutes on a slow CPU (19200 points on an
    // RK3588).
    explicit Fft(size_t n, bool quick = false);
    ~Fft();
    Fft(const Fft&) = delete;
    Fft& operator=(const Fft&) = delete;

    size_t size() const { return n_; }
    // Fill this, call forward(), read the spectrum back out of it. Unnormalised,
    // bin 0 at DC and the top half negative frequencies.
    C32* data() { return buf_; }
    const C32* data() const { return buf_; }
    void forward();

private:
    size_t n_;
    C32* buf_;
    void* plan_;
};

}  // namespace skimmer
