// FST4W's front end: what jt9's fst4_decode does from a period's samples to
// each candidate's soft bits (WSJT-X lib/fst4_decode.f90 at 567ad29ce:
// blanker, four2a's FFTs, get_candidates_fst4 with fst4_baseline, pctile,
// shell, polyfit and determ, fst4_downsample, fst4_sync_search and
// sync_fst4; lib/fst4/get_fst4_bitmetrics.f90 with ft8b.f90's
// normalizebmet), giving the reference's numbers to the bit on x86-64.
//
// What it rests on: no multiply-adds fused (-ffp-contract=off,
// CMakeLists.txt); float sums in the Fortran's order; glibc's maths
// functions where gfortran calls them (sincosf for cos/sin pairs, cabsf for
// a complex abs, lroundf for nint, log10f), and glibc's vector maths where
// the reference's -O3 vectorised a loop with it: fst4_baseline's dB
// conversion (_ZGVbN4v_log10f) and back (_ZGVbN4vv_powf), four lanes at a
// time from the loop's first element, the last 1-3 elements by the scalar
// function, as its disassembly shows; FFTW in place with
// FFTW_ESTIMATE_PATIENT, as four2a plans (its output does not depend on
// the buffers' alignment for these sizes: checked).
//
// fst4_decode.f90, get_fst4_bitmetrics.f90, fst4_baseline.f90: Steven
// Franke, K9AN, and Joe Taylor, K1JT, WSJT-X (GNU GPL v3), and so is this.
#pragma once

#include <complex>
#include <cstddef>
#include <cstdint>
#include <new>
#include <vector>

namespace fst4w {

using cfloat = std::complex<float>;

// Memory from fftwf_malloc, so every FFT buffer has the alignment the cached
// plans were made for (fft_*_inplace).
void* fftw_alloc_bytes(size_t n);
void fftw_free_bytes(void* p);
template <typename T>
struct FftwAllocator {
    using value_type = T;
    FftwAllocator() = default;
    template <typename U>
    FftwAllocator(const FftwAllocator<U>&) {}
    T* allocate(size_t n) {
        void* p = fftw_alloc_bytes(n * sizeof(T));
        if (!p) throw std::bad_alloc();
        return static_cast<T*>(p);
    }
    void deallocate(T* p, size_t) { fftw_free_bytes(p); }
    template <typename U>
    bool operator==(const FftwAllocator<U>&) const { return true; }
    template <typename U>
    bool operator!=(const FftwAllocator<U>&) const { return false; }
};
template <typename T>
using FftwVec = std::vector<T, FftwAllocator<T>>;

// fst4_decode's per-period numbers.
struct Period {
    int ntrperiod = 120;
    int nsps = 0, nmax = 0, ndown = 0, nfft1 = 0, nfft2 = 0, nss = 0;
    float fs = 12000.0f, fs2 = 0, dt = 0, dt2 = 0, tt = 0, baud = 0, sigbw = 0;
    int nspsec = 0;
    static Period of(int ntrperiod);  // 15, 30, 60, 120, 300, 900 or 1800
};

// blanker(iwave, nz, ndropmax, npct, c_bigfft): iwave's -32768s become
// -32767 (in the caller's samples, as the reference's), blanking applied,
// samples copied as interleaved pairs into c_bigfft (complex 0:nz/2; its
// last element left as it was).
void blanker(int16_t* iwave, int nz, int ndropmax, int npct, float* c_bigfft);
// The same on float samples at the 16-bit scale: |x| binned as the 16-bit
// samples' (truncated, 32768 and over in the top bin).
void blanker(float* iwave, int nz, int ndropmax, int npct, float* c_bigfft);

// four2a's in-place FFTs, planned as four2a plans them (FFTW_ESTIMATE_PATIENT,
// in place), with only the planner's own wisdom (the process's is set aside
// while planning). r2c: data[0..n) real in, n/2+1 complex out. Plans are
// made once per size and kind for the process and run on the caller's
// buffer (FFTW's new-array execute: give it FftwVec memory; a buffer
// aligned otherwise gets a plan of its own).
void fft_r2c_inplace(float* data, int n);
void fft_c2c_inplace(cfloat* data, int n, int isign);  // isign -1 forward, +1 backward

// get_candidates_fst4. candidates is 200 x 5 (Fortran's candidates(200,5):
// element (i,j) at [(j-1)*200 + i-1]). Returns false (ncand untouched) when
// the baseline has a value <= 0, as the reference returns early.
bool get_candidates_fst4(const float* c_bigfft, int nfft1, int nsps, int hmod, float fs, float fa, float fb, int nfa,
                         int nfb, float minsync, int& ncand, float* candidates);

// fst4_downsample: c1 (nfft1/ndown complex) from c_bigfft around f0.
void fst4_downsample(const float* c_bigfft, int nfft1, int ndown, float f0, float sigbw, cfloat* c1);

class Frontend {
public:
    // sync_fst4(cd0, i0, f0, hmod, ncoh, np, nss, ntr, fs, sync), its SAVEd
    // tables (common /sync240com/, f0save, nss0, ntr0) kept here.
    float sync_fst4(const cfloat* cd0, int i0, float f0, int hmod, int ncoh, int np, int nss, int ntr, float fs);

    // fst4_sync_search(c2, nfft2, hmod, fs2, nss, ntrperiod, nsyncoh,
    // emedelay, sbest, fcbest, isbest): fcbest and isbest are left as they
    // were when no sync beats the start, as there.
    void fst4_sync_search(const cfloat* c2, int nfft2, int hmod, float fs2, int nss, int ntrperiod, int nsyncoh,
                          float emedelay, float& sbest, float& fcbest, int& isbest);

    // get_fst4_bitmetrics(cd, nss, bitmetrics, s4, nsync_qual, badsync):
    // bitmetrics is 320 x 4 (element (i,j) at [(j-1)*320 + i-1]), s4 4 x 160
    // ((itone,k) at [(k-1)*4 + itone]). On badsync, what the reference left
    // filled is filled (nsync_qual only once the sequence correlations ran).
    void get_fst4_bitmetrics(const cfloat* cd, int nss, float* bitmetrics, float* s4, int& nsync_qual, bool& badsync);
    // The same, line for line (the fast one falls back to it when a value
    // is not finite; the tests compare the two).
    void get_fst4_bitmetrics_reference(const cfloat* cd, int nss, float* bitmetrics, float* s4, int& nsync_qual,
                                       bool& badsync);

private:
    // sync_fst4's table set-up (its saved state below) for f0.
    void sync_tables(float f0, int hmod, int nss, int ntr, float fs);
    // sync_fst4 at istart = first, first+step, ... (n of them), nsyncoh 8:
    // eight starts at a time, each lane its own sums in sync_fst4's order.
    void sync_lanes(const cfloat* cd0, int first, int step, int n, float f0, int hmod, int np, int nss, int ntr,
                    float fs, float* out);
    void split(const cfloat* cd0, int np, int step);
    std::vector<float> dre_, dim_;  // cd0 split by start residue (split, sync_lanes)
    int split_len_ = 0;

    // sync_fst4's saved state.
    std::vector<cfloat> csync1_ = std::vector<cfloat>(3200), csync2_ = std::vector<cfloat>(3200),
                        csynct1_ = std::vector<cfloat>(3200), csynct2_ = std::vector<cfloat>(3200);
    float twopi_ = 0, sdt_ = 0, fac_ = 0, f0save_ = -99.9f;
    int nss0_ = -1, ntr0_ = -1;
    // get_fst4_bitmetrics' saved state.
    std::vector<cfloat> ci_;  // ci(nss, 0:3): (j,itone) at [itone*nss + j-1]
    int bm_nss0_ = -1;
    bool bm_first_ = true;
    std::vector<float> s2_ = std::vector<float>(65536, 0.0f);
};

}  // namespace fst4w
