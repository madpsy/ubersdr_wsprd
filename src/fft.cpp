#include "fft.h"

#include <fftw3.h>

#include <cstdlib>
#include <mutex>
#include <new>
#include <string>

namespace skimmer {

namespace {

std::mutex& planner_lock() {
    static std::mutex m;
    return m;
}

}  // namespace

std::mutex& fft_planner_lock() { return planner_lock(); }

namespace {

bool decoder_wisdom_enabled = false;
std::string decoder_wisdom;  // FFTW's text; "" none

// The process's wisdom set aside while fn runs on the decoders' alone.
template <typename Fn>
void with_decoder_wisdom(Fn fn) {
    char* const saved = fftwf_export_wisdom_to_string();
    fftwf_forget_wisdom();
    if (!decoder_wisdom.empty()) fftwf_import_wisdom_from_string(decoder_wisdom.c_str());
    fn();
    fftwf_forget_wisdom();
    if (saved) {
        fftwf_import_wisdom_from_string(saved);
        std::free(saved);
    }
}

std::string wisdom_now() {
    char* w = fftwf_export_wisdom_to_string();
    std::string s = w ? w : "";
    std::free(w);
    return s;
}

}  // namespace

namespace {
int plan_patience = 1;
}
void set_decoder_plan_patience(int npatience) { plan_patience = npatience; }
int decoder_plan_patience() { return plan_patience; }

void decoder_wisdom_on(const std::string& import_file) {
    std::lock_guard<std::mutex> lock(planner_lock());
    decoder_wisdom.clear();
    decoder_wisdom_enabled = true;
    with_decoder_wisdom([&] {
        if (fftwf_import_wisdom_from_filename(import_file.c_str())) decoder_wisdom = wisdom_now();
    });
}

bool decoder_wisdom_export(const std::string& file) {
    std::lock_guard<std::mutex> lock(planner_lock());
    bool ok = false;
    with_decoder_wisdom([&] { ok = fftwf_export_wisdom_to_filename(file.c_str()) != 0; });
    return ok;
}

void decoder_wisdom_off() {
    std::lock_guard<std::mutex> lock(planner_lock());
    decoder_wisdom_enabled = false;
    decoder_wisdom.clear();
}

void decoder_wisdom_before_plan() {
    if (decoder_wisdom_enabled && !decoder_wisdom.empty()) fftwf_import_wisdom_from_string(decoder_wisdom.c_str());
}

void decoder_wisdom_after_plan() {
    if (decoder_wisdom_enabled) decoder_wisdom = wisdom_now();
}

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
