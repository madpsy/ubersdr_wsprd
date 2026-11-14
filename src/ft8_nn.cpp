// ft8_nn.h: the network through ONNX Runtime's C API; its input made exactly as in training (to the
// bit), so it sees exactly what it was trained on. ONNX Runtime's library is loaded on the first model opened
// (dlopen, LoadLibrary): only its headers are built in (vendor/onnxruntime), so nothing needs it otherwise.
#include "ft8_nn.h"

#ifndef __EMSCRIPTEN__  // the WebAssembly build: no ONNX Runtime, every model refused
#include "onnxruntime_c_api.h"

#include <algorithm>
#include <cmath>
#include <mutex>
#include <vector>
#ifdef _WIN32
#include <windows.h>
#else
#include <dlfcn.h>
#include <unistd.h>
#endif
#endif

namespace ft8 {

#ifndef __EMSCRIPTEN__
namespace {
#ifdef _WIN32
constexpr const char* kOrtLib = "onnxruntime.dll";
#elif defined(__APPLE__)
constexpr const char* kOrtLib = "libonnxruntime.dylib";
#else
constexpr const char* kOrtLib = "libonnxruntime.so.1";
#endif

std::mutex ort_mu;
const OrtApi* ort_loaded = nullptr;
std::string ort_why, ort_given;  // why it failed; NnMetric::use_library's path
bool ort_tried = false;

// ONNX Runtime's API, its library loaded once: use_library's if given, else the program's own directory (Windows'
// LoadLibrary looks there first anyway), then the system's library path. Null with the reason in err.
const OrtApi* ort_api(std::string& err) {
    std::lock_guard<std::mutex> lock(ort_mu);
    const OrtApi*& api = ort_loaded;
    std::string& why = ort_why;
    if (ort_tried) {
        err = why;
        return api;
    }
    ort_tried = true;
    void* fn = nullptr;
#ifdef _WIN32
    HMODULE h = !ort_given.empty() ? LoadLibraryExA(ort_given.c_str(), nullptr, LOAD_WITH_ALTERED_SEARCH_PATH)
                                   : LoadLibraryA(kOrtLib);
    if (h) fn = reinterpret_cast<void*>(GetProcAddress(h, "OrtGetApiBase"));
#else
    void* h = nullptr;
    if (!ort_given.empty()) {
        h = dlopen(ort_given.c_str(), RTLD_NOW | RTLD_LOCAL);
    } else {
        char exe[4096];
        const ssize_t n = readlink("/proc/self/exe", exe, sizeof exe - 1);
        if (n > 0) {
            std::string dir(exe, static_cast<size_t>(n));
            dir = dir.substr(0, dir.rfind('/') + 1);
            h = dlopen((dir + kOrtLib).c_str(), RTLD_NOW | RTLD_LOCAL);
        }
        if (!h) h = dlopen(kOrtLib, RTLD_NOW | RTLD_LOCAL);
    }
    if (h) fn = dlsym(h, "OrtGetApiBase");
#endif
    if (!fn && !ort_given.empty()) {
        why = "cannot load ONNX Runtime from " + ort_given;
    } else if (!fn) {
        why = std::string("ONNX Runtime not found: put ") + kOrtLib +
              " (its official release, 1.31.1 or newer) beside the program or on the library path";
    } else {
        const OrtApiBase* base = reinterpret_cast<const OrtApiBase*(ORT_API_CALL*)(void)>(fn)();
        api = base->GetApi(ORT_API_VERSION);
        if (!api) why = std::string("ONNX Runtime ") + base->GetVersionString() + " is too old (1.31.1 or newer)";
    }
    err = why;
    return api;
}

constexpr int kNBB = 3200;  // the baseband the features read (cd0's length)
constexpr int kNN = 79, kSps = 32, kNbinLo = -4, kNbins = 24, kNFEAT = 2 * kNbins;

// The network's input, to the bit: for each
// of the 79 symbols the 64-point DFT of its 32 samples (as float, zero outside cd0) at bins -4..19, scaled by
// 1/sqrt(32), then every value divided by the square root of the median bin power (the lower middle one).
struct Dft64 {
    std::complex<double> w[kNbins][kSps];
    Dft64() {
        for (int q = 0; q < kNbins; ++q)
            for (int n = 0; n < kSps; ++n) {
                const int b = ((kNbinLo + q) % 64 + 64) % 64;
                const double a = -2 * M_PI * b * n / 64.0;
                w[q][n] = std::complex<double>(std::cos(a), std::sin(a));
            }
    }
};

void nn_features(const std::complex<double>* y, int start, float* out) {
    static const Dft64 dft;  // thread-safe initialisation
    static const double isq = 0.17677669529663687;  // 1/sqrt(32)
    double re[kNN][kNbins], im[kNN][kNbins];
    float p[kNN * kNbins];
    for (int s = 0; s < kNN; ++s) {
        std::complex<double> seg[kSps];
        for (int n = 0; n < kSps; ++n) {
            const int idx = start + s * kSps + n;
            seg[n] = idx >= 0 && idx < kNBB ? std::complex<double>(static_cast<float>(y[idx].real()),
                                                                   static_cast<float>(y[idx].imag()))
                                            : std::complex<double>(0, 0);
        }
        for (int q = 0; q < kNbins; ++q) {
            std::complex<double> z = 0;
            for (int n = 0; n < kSps; ++n) z += seg[n] * dft.w[q][n];
            z *= isq;
            re[s][q] = z.real();
            im[s][q] = z.imag();
            const float a = static_cast<float>(std::abs(z));
            p[s * kNbins + q] = a * a;
        }
    }
    float* mid = p + (kNN * kNbins - 1) / 2;
    std::nth_element(p, mid, p + kNN * kNbins);
    const double scale = std::sqrt(std::max(static_cast<double>(*mid), 1e-6));
    for (int s = 0; s < kNN; ++s)
        for (int q = 0; q < kNbins; ++q) {
            out[s * kNFEAT + q] = static_cast<float>(re[s][q] / scale);
            out[s * kNFEAT + kNbins + q] = static_cast<float>(im[s][q] / scale);
        }
}
}  // namespace

struct NnMetric::Impl {
    const OrtApi* api = nullptr;
    OrtEnv* env = nullptr;
    OrtSession* sess = nullptr;
    OrtMemoryInfo* mem = nullptr;
    ~Impl() {
        if (!api) return;
        if (mem) api->ReleaseMemoryInfo(mem);
        if (sess) api->ReleaseSession(sess);
        if (env) api->ReleaseEnv(env);
    }
    bool check(OrtStatus* st, std::string* err) {
        if (!st) return true;
        if (err) *err = api->GetErrorMessage(st);
        api->ReleaseStatus(st);
        return false;
    }
};

namespace {
// A model's path as ONNX Runtime takes it: UTF-16 on Windows (ORTCHAR_T wchar_t), as given elsewhere.
std::basic_string<ORTCHAR_T> ort_path(const std::string& path) {
#ifdef _WIN32
    const int n = MultiByteToWideChar(CP_UTF8, 0, path.c_str(), -1, nullptr, 0);
    std::wstring w(n > 0 ? static_cast<size_t>(n) : 1, L'\0');
    if (n > 0) MultiByteToWideChar(CP_UTF8, 0, path.c_str(), -1, &w[0], n);
    w.resize(w.size() - 1);
    return w;
#else
    return path;
#endif
}
}  // namespace

void NnMetric::use_library(const std::string& path) {
    std::lock_guard<std::mutex> lock(ort_mu);
    if (ort_loaded) return;  // one library for the program's life
    ort_given = path;
    ort_tried = false;
}

// A session for the model at path, or (path empty) the size bytes at data.
std::unique_ptr<NnMetric> NnMetric::make(const std::string& path, const void* data, size_t size, std::string& err) {
    std::unique_ptr<NnMetric> m(new NnMetric);
    m->p_ = std::make_unique<Impl>();
    Impl& p = *m->p_;
    p.api = ort_api(err);
    if (!p.api) return nullptr;
    OrtSessionOptions* so = nullptr;
    bool ok = p.check(p.api->CreateEnv(ORT_LOGGING_LEVEL_WARNING, "ft8", &p.env), &err) &&
              p.check(p.api->CreateSessionOptions(&so), &err) &&
              p.check(p.api->SetIntraOpNumThreads(so, 1), &err) &&
              p.check(p.api->SetInterOpNumThreads(so, 1), &err) &&
              p.check(p.api->SetSessionGraphOptimizationLevel(so, ORT_ENABLE_ALL), &err);
    if (ok)
        ok = path.empty() ? p.check(p.api->CreateSessionFromArray(p.env, data, size, so, &p.sess), &err)
                          : p.check(p.api->CreateSession(p.env, ort_path(path).c_str(), so, &p.sess), &err);
    ok = ok && p.check(p.api->CreateCpuMemoryInfo(OrtArenaAllocator, OrtMemTypeDefault, &p.mem), &err);
    if (so) p.api->ReleaseSessionOptions(so);
    if (!ok) return nullptr;
    return m;
}

std::unique_ptr<NnMetric> NnMetric::open(const std::string& path, std::string& err) {
    if (path.empty()) {
        err = "no model";
        return nullptr;
    }
    return make(path, nullptr, 0, err);
}

std::unique_ptr<NnMetric> NnMetric::open_memory(const void* data, size_t size, std::string& err) {
    return make("", data, size, err);
}

NnMetric::~NnMetric() = default;

bool NnMetric::llr(const std::complex<float>* cd0, int n, int ibest, float out[174], float df) {
    // Safe from several threads at once (one model for every decoder): the buffers are the thread's own, nn_features
    // only reads its table, and ONNX Runtime runs one session concurrently.
    thread_local std::vector<std::complex<double>> y(kNBB);
    thread_local std::vector<float> feat(kNN * kNFEAT);
    Impl& p = *p_;
    const double w = -2.0 * 3.14159265358979323846 * df / 200.0;  // cd0: 200 S/s
    for (int i = 0; i < kNBB; ++i)
        y[i] = i < n ? std::complex<double>(cd0[i]) * (df != 0.0f ? std::polar(1.0, w * i) : std::complex<double>(1, 0))
                       : std::complex<double>(0, 0);
    nn_features(y.data(), ibest, feat.data());
    const int64_t shape[3] = {1, kNN, kNFEAT};
    const char* in_names[] = {"features"};
    const char* out_names[] = {"llr"};
    OrtValue *in = nullptr, *res = nullptr;
    bool ok = p.check(p.api->CreateTensorWithDataAsOrtValue(p.mem, feat.data(), sizeof(float) * feat.size(),
                                                            shape, 3, ONNX_TENSOR_ELEMENT_DATA_TYPE_FLOAT, &in),
                      nullptr) &&
              p.check(p.api->Run(p.sess, nullptr, in_names, &in, 1, out_names, 1, &res), nullptr);
    float* o = nullptr;
    if (ok) ok = p.check(p.api->GetTensorMutableData(res, reinterpret_cast<void**>(&o)), nullptr);
    if (ok)
        for (int i = 0; i < 174; ++i) out[i] = o[i];
    if (in) p.api->ReleaseValue(in);
    if (res) p.api->ReleaseValue(res);
    return ok;
}

#else

struct NnMetric::Impl {};

std::unique_ptr<NnMetric> NnMetric::open(const std::string&, std::string& err) {
    err = "this build has no ONNX Runtime (WebAssembly), so no model can be used";
    return nullptr;
}
std::unique_ptr<NnMetric> NnMetric::open_memory(const void*, size_t, std::string& err) { return open("", err); }
void NnMetric::use_library(const std::string&) {}
NnMetric::~NnMetric() = default;
bool NnMetric::llr(const std::complex<float>*, int, int, float*, float) { return false; }

#endif

}  // namespace ft8
