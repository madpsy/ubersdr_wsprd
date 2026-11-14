// Ours: a neural soft demodulator as one more FT8 bit metric (jt9 -8 --nn MODEL).
//
// Off unless --nn is given; then ft8b tries the network's LLRs as extra passes when WSJT-X's own five metrics fail
// (before its a-priori passes, which then start from the network's LLRs), on the same downsampled signal (cd0) at
// the same refined frequency and start plus small offsets (--nn-level), through the same decode174_91, and drops
// implausible OSD and a-priori decodes (ft8_plausible.h). With --nn every decode can differ from WSJT-X's. Needs
// ONNX Runtime's library at run time, loaded on the first model (ft8_nn.cpp); nothing else does.
#pragma once

#include <complex>
#include <cstddef>
#include <memory>
#include <string>

namespace ft8 {

class NnMetric {
public:
    // The model (an .onnx file), or nullptr with the reason in err.
    static std::unique_ptr<NnMetric> open(const std::string& path, std::string& err);
    // The same from the model's bytes in memory (size of them at data; kept only while opening).
    static std::unique_ptr<NnMetric> open_memory(const void* data, size_t size, std::string& err);
    // ONNX Runtime's library to load (a program carrying its own; before the first open): this file, not the one
    // beside the program or on the library path. No effect once one is loaded.
    static void use_library(const std::string& path);
    ~NnMetric();
    NnMetric(const NnMetric&) = delete;
    NnMetric& operator=(const NnMetric&) = delete;

    // 174 LLRs (log P(1)/P(0), WSJT-X's codeword order) for the signal whose symbol 0 starts at cd0[ibest];
    // cd0 is ft8b's 200 S/s baseband (n samples, tone 0 at 0 Hz). False if the network could not run.
    // df: assume the signal sits df Hz above cd0's tone 0 (the input is shifted down by df first)
    // Thread-safe: one NnMetric can serve every decoder.
    bool llr(const std::complex<float>* cd0, int n, int ibest, float out[174], float df = 0.0f);

private:
    NnMetric() = default;
    static std::unique_ptr<NnMetric> make(const std::string& path, const void* data, size_t size, std::string& err);
    struct Impl;
    std::unique_ptr<Impl> p_;
};

}  // namespace ft8
