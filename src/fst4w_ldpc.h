// FST4W's forward error correction: WSJT-X's (240,74) LDPC code with its
// 24-bit CRC (lib/fst4 at 567ad29ce: decode240_74.f90, fastosd240_74.f90,
// encode240_74.f90, get_crc24.f90, ldpc_240_74_*.f90, and lib/platanh.f90),
// giving the reference's results to the bit: the same message, codeword,
// hard errors, dmin and decoder type from the same soft bits, and the same
// cw and message74 left behind when decoding fails (fst4_decode unpacks a
// failed Keff=50 decode's message74 all the same).
// decode240_74's hard-error count is the sum of an integer*1 array, so
// 128 or more wrap negative (240 reads -16), as the reference's do.
//
// What it rests on: no multiply-adds fused (-ffp-contract=off,
// CMakeLists.txt), as the reference's Release build has none on x86-64;
// float sums in the Fortran's order (no reassociation); glibc's tanhf, as
// gfortran calls for tanh of a default real (scalar here: the reference's
// loop is not vectorised); indexx is lib/indexx.f90's (wsprd_osd.h).
//
// decode240_74.f90, fastosd240_74.f90, encode240_74.f90, get_crc24.f90:
// Steven Franke, K9AN, and Joe Taylor, K1JT, WSJT-X (GNU GPL v3), and so is
// this.
#pragma once

#include <cstdint>
#include <vector>

namespace fst4w {

// get_crc24(mc, len, ncrc): the CRC register after mc(1:len) (bits 0/1);
// zero when mc ends in its own CRC.
int get_crc24(const int8_t* mc, int len);

// encode240_74(message, codeword): 74 message bits (50 + CRC24) to 240.
void encode240_74(const int8_t message[74], int8_t codeword[240]);

class Ldpc {
public:
    // decode240_74(llr, Keff, maxosd, norder, apmask, message74, cw, ntype,
    // nharderror, dmin). message74, cw, ntype, nharderror and dmin are
    // written as the Fortran writes its arguments, and only then: dmin is
    // left alone when decoding fails, as there.
    void decode240_74(const float llr[240], int Keff, int maxosd, int norder, const int8_t apmask[240],
                      int8_t message74[74], int8_t cw[240], int& ntype, int& nharderror, float& dmin);

    // fastosd240_74(llr, k, apmask, ndeep, message74, cw, nhardmin, dmin),
    // on 240-bit words: the same results as the line-for-line version below,
    // which it falls back to when a soft bit is not finite.
    void fastosd240_74(const float llr[240], int k, const int8_t apmask[240], int ndeep, int8_t message74[74],
                       int8_t cw[240], int& nhardmin, float& dmin);

    // The Fortran line for line (byte arrays), for the fallback and the tests.
    void fastosd240_74_reference(const float llr[240], int k, const int8_t apmask[240], int ndeep, int8_t message74[74],
                                 int8_t cw[240], int& nhardmin, float& dmin);

private:
    // fastosd's SAVEd generator matrix (k x 240, row by row) and its k.
    std::vector<int8_t> gen_;
    int ksave_ = 64;
    bool first_ = true;
    // The fast version's generators, one per k (50..74), as 240-bit rows.
    std::vector<std::vector<uint64_t>> genbits_ = std::vector<std::vector<uint64_t>>(75);
};

}  // namespace fst4w
