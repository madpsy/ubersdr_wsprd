#include "fft.h"

#include <fftw3.h>

#include <mutex>
#include <new>

namespace skimmer {

namespace {

std::mutex& planner_lock() {
    static std::mutex m;
    return m;
}

}  // namespace

std::mutex& fft_planner_lock() { return planner_lock(); }

Fft::Fft(size_t n, bool quick) : n_(n) {
    std::lock_guard<std::mutex> lock(planner_lock());
    buf_ = reinterpret_cast<C32*>(fftwf_malloc(sizeof(fftwf_complex) * n));
    if (!buf_) throw std::bad_alloc();
    auto* p = reinterpret_cast<fftwf_complex*>(buf_);
    const int size = static_cast<int>(n);
    // MEASURE: timed against the alternatives on this machine, 15-30% faster
    // than ESTIMATE's guess and found in under a second at these sizes, then
    // saved, so the next run starts at once. PATIENT and EXHAUSTIVE gain
    // another 10-15% only at 8192 points (iq384), for 14 s and 5 minutes of
    // planning; wisdom from either is used in preference whenever the file
    // has it (fft_load_wisdom(), fft_save_wisdom(): fft_wisdom.cpp).
    fft_load_wisdom();
    plan_ = fftwf_plan_dft_1d(size, p, p, FFTW_FORWARD, FFTW_MEASURE | FFTW_WISDOM_ONLY);
    if (!plan_ && quick) plan_ = fftwf_plan_dft_1d(size, p, p, FFTW_FORWARD, FFTW_ESTIMATE);
    if (!plan_) {
        plan_ = fftwf_plan_dft_1d(size, p, p, FFTW_FORWARD, FFTW_MEASURE);
        fft_save_wisdom();
    }
    // Planning without wisdom overwrites the buffer.
    for (size_t i = 0; i < n; ++i) buf_[i] = C32();
}

Fft::~Fft() {
    std::lock_guard<std::mutex> lock(planner_lock());
    fftwf_destroy_plan(static_cast<fftwf_plan>(plan_));
    fftwf_free(buf_);
}

void Fft::forward() { fftwf_execute(static_cast<fftwf_plan>(plan_)); }

}  // namespace skimmer
