// MSK144's short-ping detector and sync: WSJT-X's msk144spd.f90,
// msk144sync.f90, msk_spectrum.f90 (msk_clip_spectrum_window, msk_peak_offset)
// and indexx.f90 at 567ad29ce (3.2.0-rc1), line for line, giving the
// reference's numbers to the bit on x86-64. It rests on what
// src/msk_frontend.h rests on (no multiply-adds fused, complex products
// written out as gfortran does, glibc's maths where gfortran calls it, FFTW as
// four2a plans it), and on what the reference's compiler did with these
// routines: nint of a double is lround, of a float lroundf; abs of a complex is
// cabsf; MAX and MAXLOC as gfortran -O3 writes them for NaNs (the comments in
// src/msk144_sync.cpp say where).
//
// msk144spd's frame decoder (the reference's msk144decodeframe, a later
// piece) comes in as a callback.
//
// Joe Taylor, K1JT, Steven Franke, K9AN, and the WSJT-X team (GNU GPL v3),
// and so is this.
#pragma once

#include <complex>
#include <functional>

#include "msk_frontend.h"

namespace msk144 {

// msk144decodeframe(c, softbits, msgreceived, nsuccess): c the NSPM = 864
// points of a frame, which the decoder corrects in place (the phase
// correction, msk144decodeframe.f90:59: msk144signalquality later reads the
// corrected frame, so ct comes back from msk144spd as the decoder left it),
// softbits 144 values out, msgreceived 37 characters (blank
// padded, as Fortran's character*37) out, nsuccess out (> 0: decoded). The
// short-ping detector hands each frame it cuts to this.
using DecodeSink = std::function<void(cfloat* c, float* softbits, char* msgreceived, int& nsuccess)>;

// msk144sync(cdat, nframes, ntol, delf, navmask, npeaks, fc, fest, npklocs,
// nsuccess, xmax, c): the Fortran's argument list. cdat is 864*nframes
// points, navmask nframes entries, npklocs npeaks entries (1-based positions
// in 1..864, as the Fortran gives them), c 864 points. Not threaded, though the
// reference's jt9 links the OpenMP library (libwsjt_fort_omp.a, libgomp), where
// msk144sync runs nthreads = min(4, OMP_GET_MAX_THREADS()) threads, each
// with a slice of the trial frequencies, merged by a strict >. With a full team
// the serial search gives the same bits (the piece check runs the oracle against
// the OpenMP library at four threads). Two divergences remain, neither seen in
// the decoder's output: (a) with nfreqs < nthreads and every trial's peak zero
// or not finite, fest differs (the reference's last thread has the slice
// holding if2 only), but nsuccess is 0 and the callers cycle; (b) with a short
// team (OMP_THREAD_LIMIT below the request, or omp_set_dynamic(.true.), which
// decoder.f90 sets after the first file's MSK pass) the reference leaves
// slices unsearched, so its answer depends on the load in multi-file runs.
void msk144sync(const cfloat* cdat, int nframes, int ntol, float delf, const int* navmask, int npeaks, float fc, float& fest,
                int* npklocs, int& nsuccess, float& xmax, cfloat* c);

// msk144spd(cbig, n, ntol, nsuccess, msgreceived, fc, fret, tret, navg, ct,
// softbits): the Fortran's argument list, the frame decoder as `decode`.
// msgreceived is 37 characters, ct 864 points, softbits 144 values; fret, tret,
// navg, ct and softbits are written only as the Fortran writes them.
void msk144spd(const cfloat* cbig, int n, int ntol, int& nsuccess, char* msgreceived, float fc, float& fret, float& tret,
               int& navg, cfloat* ct, float* softbits, const DecodeSink& decode);

// The pieces msk144spd is made of, public for the piece check.
// indexx(arr, n, indx): the quicksort of indices, indx(1:n) = 1-based positions.
void indexx(const float* arr, int n, int* indx);
// maxloc(a(1:n)) and maxloc(a, mask) with the mask true on lo..hi (1-based),
// gfortran's way with NaNs: the first non-NaN element is the start, later ones
// replace it only if greater; all NaN gives the first element in the mask; 0 if
// the mask is empty.
int maxloc_f(const float* a, int n);
int maxloc_masked_f(const float* a, int lo, int hi);
// msk_clip_spectrum_window(nfft, requested_lo, requested_hi, lo, hi, usable)
void msk_clip_spectrum_window(int nfft, double requested_lo, double requested_hi, int& lo, int& hi, bool& usable);
// msk_peak_offset(spectrum(1:n), peak_bin)
float msk_peak_offset(const cfloat* spectrum, int n, int peak_bin);

// Debug-only traces for the piece check: null in production (an untaken
// branch at each place). Each is called where scripts/msk144-piece-check.sh
// sync's instrumented copy of the reference calls its hook, with what that
// hook prints.
struct SyncTrace {
    void (*consts)(const float* pp /*12*/, const cfloat* cb /*42*/) = nullptr;                          // first call
    void (*peak)(float xmax, float fest, const cfloat* c, const float* xcc, int nfreqs) = nullptr;     // before the peak search
};
struct SpdTrace {
    void (*consts)(const float* rcw /*12*/) = nullptr;  // first call
    void (*windows)(int ihlo, int ihhi, int illo, int ilhi, int i2000, int i4000) = nullptr;
    void (*step)(int istp, int ihpk, int ilpk, float deltah, float deltal, float ah, float al, float trath, float tratl,
                 float ferr) = nullptr;
    // after the median: detmet, detmet2, detfer (1..nstep, as the steps made them), the indices (1-based), xmed
    void (*median)(int nstep, const int* indices, float xmed, const float* detmet, const float* detmet2, const float* detfer) = nullptr;
    void (*candidates)(int ndet, const int* nstart, const float* ferrs, const float* snrs) = nullptr;
    void (*sync)(int icand, int iav, int nsyncsuccess, float fest, float xmax, const int* npkloc) = nullptr;
};
extern SyncTrace sync_trace;
extern SpdTrace spd_trace;

}  // namespace msk144
