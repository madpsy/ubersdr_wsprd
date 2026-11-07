// MSK40, the short-message form of MSK144 (WSJT-X's bshmsg), and the
// decoder pieces it takes: WSJT-X's msk40spd.f90, msk40sync.f90,
// msk40_freq_search.f90, msk40decodeframe.f90, bpdecode40.f90,
// encode_msk40.f90, update_msk40_hasharray.f90 and hash.f90 (over lookup3.c's
// hashlittle) at 567ad29ce (3.2.0-rc1), line for line, giving the reference's
// numbers to the bit on x86-64. It rests on what src/msk144_sync.h and
// src/msk144_decode.h rest on (no multiply-adds fused, complex products
// written out as gfortran does, glibc's maths where gfortran calls it, FFTW
// as four2a plans it).
//
// An MSK40 message carries a 12-bit hash of "MYCALL HISCALL" and a 4-bit
// report: the receiver prints "<MYCALL HISCALL> RPT". With bswl (the SWL
// mode) it also tries the hashes of every pair of the calls it has heard
// (packjt77's recent_calls), and prints "<HHHH> RPT" for a hash it cannot
// resolve.
//
// What the reference saves from call to call, and the port keeps in an
// Msk40State: mskrtd's nhasharray (the hashes of the recent calls, zero
// until update_msk40_hasharray first runs), and msk40decodeframe's mycall0,
// hiscall0 and ihash. It is part of mskrtd's RtdState (src/msk144.h), so it
// restarts with every file, --feed included, as a fresh jt9 process would.
//
// msk40sync is not threaded here, though the reference's jt9 links the OpenMP
// library (msk40sync runs nthreads = min(8, OMP_GET_MAX_THREADS()) threads,
// each with a slice of the trial frequencies, merged by a strict >). The
// slices cover every trial frequency, and the merge takes the first of the
// largest, as the serial loop does: the same bits. The piece check runs the
// oracle against the OpenMP library at 8 and at 4 threads to show it.
//
// Joe Taylor, K1JT, Steven Franke, K9AN, and the WSJT-X team (GNU GPL v3),
// and so is this. lookup3.c: Bob Jenkins, 2006, public domain.
#pragma once

#include <string>

#include "msk_frontend.h"
#include "packjt77.h"

namespace msk144 {

using packjt77::Packjt77;

// What the reference's MSK40 routines save (see above). nhasharray(i,j) is
// nhasharray[i-1][j-1].
struct Msk40State {
    int nhasharray[Packjt77::MAXRECENT][Packjt77::MAXRECENT] = {};  // mskrtd: saved, no initialiser: zero
    std::string mycall0 = "dummy       ";                           // msk40decodeframe: data mycall0/'dummy'/
    std::string hiscall0 = "dummy       ";                          //   and hiscall0/'dummy'/ (character*12)
    int ihash = 0;                                                   //   saved, no initialiser: zero
};

// hash(string, len, ihash): the 15 low bits of nhash(string, len, 146).
int hash_f90(const char* string, int len);

// encode_msk40(message, codeword): the 16 message bits (0/1) to the 32-bit
// codeword (0/1).
void encode_msk40(const signed char message[16], signed char codeword[32]);

// bpdecode40(llr, maxiterations, decoded, niterations): decoded (16 bits,
// 0/1) is written only when niterations >= 0 (a codeword was found after
// that many iterations), as the Fortran's.
void bpdecode40(const float llr[32], int maxiterations, signed char decoded[16], int& niterations);

// update_msk40_hasharray(nhasharray): from the calls pk holds (recent_calls).
void update_msk40_hasharray(Msk40State& ms, const Packjt77& pk);

// msk40decodeframe(c, mycall, hiscall, xsnr, bswl, nhasharray, msgreceived,
// nsuccess): c the NSPM = 240 points of a frame, which the decoder corrects
// in place; mycall and hiscall character*12 (cut or blank padded to 12);
// msgreceived 37 characters (blank padded), nsuccess out: 0 not decoded, 1 a
// message for mycall and hiscall, 2 (bswl) the hash of two recent calls, 3
// (bswl) a hash no recent call explains.
void msk40decodeframe(Msk40State& ms, const Packjt77& pk, cfloat* c, const std::string& mycall12,
                      const std::string& hiscall12, float xsnr, bool bswl, char* msgreceived, int& nsuccess);

// msk40_freq_search(cdat, fc, if1, if2, delf, nframes, navmask, cb, cdat2,
// xmax, bestf, cs, xccs), with NSPM = 240: cdat and cdat2 are 240*nframes
// points, xccs 240 values.
void msk40_freq_search(const cfloat* cdat, float fc, int if1, int if2, float delf, int nframes, const int* navmask,
                       const cfloat* cb, cfloat* cdat2, float& xmax, float& bestf, cfloat* cs, float* xccs);

// msk40sync(cdat, nframes, ntol, delf, navmask, npeaks, fc, fest, npklocs,
// nsuccess, c): cdat is 240*nframes points, npklocs npeaks entries (1-based
// positions in 1..240), c 240 points.
void msk40sync(const cfloat* cdat, int nframes, int ntol, float delf, const int* navmask, int npeaks, float fc,
               float& fest, int* npklocs, int& nsuccess, cfloat* c);

// msk40spd(cbig, n, ntol, mycall, hiscall, bswl, nhasharray, nsuccess,
// msgreceived, fc, fret, tret, navg): msgreceived is 37 characters; fret, tret
// and navg are written only on a decode, as the Fortran writes them.
void msk40spd(Msk40State& ms, const Packjt77& pk, const cfloat* cbig, int n, int ntol, const std::string& mycall12,
              const std::string& hiscall12, bool bswl, int& nsuccess, char* msgreceived, float fc, float& fret,
              float& tret, int& navg);

// Debug-only traces for the piece check: null in production (an untaken
// branch at each place). Each is called where scripts/msk144-piece-check.sh
// msk40's instrumented copy of the reference calls its hook, with what that
// hook prints.
struct Msk40Trace {
    // bpdecode40: at the top of each iteration, once zn is updated; zn[1..32], tov[1..3][1..32] (1-based)
    void (*bp)(int iter, const float zn[33], const float tov[4][33]) = nullptr;
    // msk40decodeframe
    void (*df_consts)(const float* pp /*12*/, const cfloat* cb /*42*/) = nullptr;  // first call
    void (*df_phase)(cfloat cca, float phase0) = nullptr;
    void (*df_corrected)(cfloat cfac, const cfloat* c /*240*/) = nullptr;  // after c=c*conjg(cfac)
    void (*df_raw)(const float* softbits /*40*/) = nullptr;                 // after the matched filter
    void (*df_sync)(int nbadsync1) = nullptr;
    void (*df_norm)(float sav, float s2av, float ssig) = nullptr;
    void (*df_llr)(float sigma, const float* llr /*32*/) = nullptr;
    void (*df_bp)(const signed char* decoded /*16*/, int niterations) = nullptr;  // decoded as left
    void (*df_msg)(int nhammd, float cord, int imsg, int nrxrpt, int nrxhash, int ihash) = nullptr;
    // msk40sync: pp and cb at its first call; xmax, fest, c and xcc before its peak search
    void (*sync_consts)(const float* pp /*12*/, const cfloat* cb /*42*/) = nullptr;
    void (*sync_peak)(float xmax, float fest, const cfloat* c, const float* xcc, int nfreqs) = nullptr;
    // msk40spd
    void (*spd_consts)(const float* rcw /*12*/) = nullptr;  // first call
    void (*spd_windows)(int ihlo, int ihhi, int illo, int ilhi, int i2000, int i4000) = nullptr;
    void (*spd_step)(int istp, int ihpk, int ilpk, float deltah, float deltal, float ah, float al, float trath,
                     float tratl, float ferr) = nullptr;
    void (*spd_median)(int nstep, const int* indices, float xmed, const float* detmet, const float* detmet2,
                       const float* detfer) = nullptr;
    void (*spd_candidates)(int ndet, const int* nstart, const float* ferrs, const float* snrs) = nullptr;
    void (*spd_sync)(int icand, int iav, int nsyncsuccess, float fest, const int* npkloc) = nullptr;
};
extern Msk40Trace msk40_trace;

}  // namespace msk144
