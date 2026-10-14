// A forward complex FFT of fixed size over FFTW, transforming in place in a
// buffer it owns. Planning is serialised because FFTW's planner is not
// thread-safe; executing a plan is.
#pragma once

#include <complex>
#include <cstddef>
#include <mutex>
#include <string>

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

// The wisdom WSJT-X's decoders' own plans (four2a's, FFTW_ESTIMATE_PATIENT:
// fst4w_frontend.cpp, ft8_frontend.cpp) are made with, apart from the
// process's (set aside while they plan, so the skimmer's measured plans
// never change their output). Off, as in the skimmer: none, each plan the
// planner's own. On, as in the jt9 program: what WSJT-X's jt9 has, its
// jt9_wisdom.dat imported at the start and every plan's added, exported at
// the end. All under fft_planner_lock() (the functions take it).
void decoder_wisdom_on(const std::string& import_file);  // the file missing or bad: none to start with
bool decoder_wisdom_export(const std::string& file);      // as fftwf_export_wisdom_to_filename
void decoder_wisdom_off();
// For the planners, holding fft_planner_lock(), after forgetting the
// process's wisdom: before a plan, the decoders' wisdom imported (when on);
// after it, what the planner then knows kept as theirs (when on).
void decoder_wisdom_before_plan();
void decoder_wisdom_after_plan();

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
