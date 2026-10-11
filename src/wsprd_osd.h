// WSPR ordered-statistics decoding: wsprd's fallback when the Fano decoder
// gives up on a candidate's 162 soft symbols. WSJT-X's lib/wsprd/osdwspr.f90
// (osdwspr, mrbencode, nextpat, boxit, fetchit) and the lib/indexx.f90 sort
// it ranks the symbols with, as C++ that gives the reference's results bit
// for bit: the same codeword, nhardmin and dmin from the same inputs.
//
// osdwspr.f90: Steven Franke, K9AN, and Joe Taylor, K1JT, WSJT-X (GPLv3),
// and so is this.
#pragma once

#include <memory>

namespace wsprd {

// boxit and fetchit's tables (the Fortran's COMMON /boxes/, and fetchit's
// SAVEd lastpat and inext), about 330 KB, made on first use. Every call
// that reads them resets them first, as the reference does, so which a call
// is given changes nothing: a caller decoding many candidates keeps one
// (wsprd's Decoder does), rather than each call making its own.
struct OsdBoxes;
struct OsdBoxesFree {
    void operator()(OsdBoxes* b) const;
};
using OsdTables = std::unique_ptr<OsdBoxes, OsdBoxesFree>;

// osdwspr(ss, apmask, ndeep, cw, nhardmin, dmin) as WSJT-X's Fortran: ss the
// 162 soft symbols (wsprd passes symbol - 128, so -128..127), apmask/cw 0/1
// bytes. ndeep 0 (or less) is order 0 only, 1..5 the deeper searches, more
// than 5 is 5 (the reference also writes 5 back to its caller's ndepth;
// wsprd never notices, as 5 and above decode alike). Thread safe, given
// tables of its own (nullptr: tables made for the call). Not a thread_local
// (the reference's are global): mingw-w64's run a thread's thread_local
// destructors after its storage is freed, which crashes as a thread ends.
void osdwspr(const float ss[162], const unsigned char apmask[162], int ndeep, unsigned char cw[162], int* nhardmin, float* dmin,
             OsdTables* tables = nullptr);

// indexx(arr, n, indx) of lib/indexx.f90: indx (1-based, as the Fortran's)
// lists arr in ascending order, ties broken as its Numerical Recipes
// quicksort does. Exposed for the tests.
void indexx(const float* arr, int n, int* indx);

}  // namespace wsprd
