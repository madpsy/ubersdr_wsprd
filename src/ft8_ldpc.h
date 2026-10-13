// FT8's forward error correction: WSJT-X's (174,91) LDPC code with its
// 14-bit CRC (lib/ft8 at 567ad29ce: decode174_91.f90, osd174_91.f90,
// encode174_91.f90, encode174_91_nocrc.f90, get_crc14.f90,
// ldpc_174_91_c_*.f90; lib/crc14.cpp and lib/platanh.f90), giving the
// reference's results to the bit: the same message, codeword, hard errors,
// dmin and decoder type from the same soft bits, and the same cw and
// message91 left behind when decoding fails.
// Unlike decode240_74, decode174_91's hard-error count is a default-integer
// COUNT: no integer*1 wrap.
//
// What it rests on: no multiply-adds fused (-ffp-contract=off,
// CMakeLists.txt), as the reference's Release build has none on x86-64;
// float sums in the Fortran's order (no reassociation); glibc's tanhf, as
// gfortran calls for tanh of a default real (scalar: the reference's loop
// is unrolled, not vectorised, and calls tanhf seven times a check), here
// a copy of glibc 2.39's (tanhf_glibc);
// indexx is lib/indexx.f90's (wsprd_osd.h).
//
// decode174_91.f90, osd174_91.f90, encode174_91*.f90, get_crc14.f90,
// crc14.cpp: Steven Franke, K9AN, and Joe Taylor, K1JT, WSJT-X (GNU GPL
// v3), and so is this.
#pragma once

#include <cstdint>

namespace ft8 {

// get_crc14(mc, len, ncrc): the CRC register after mc(1:len) (bits 0/1);
// zero when mc ends in its own CRC. The decoders call it on m96: message
// bits 1-77, five zeros, the 14 CRC bits.
int get_crc14(const int8_t* mc, int len);

// crc14(data, length) of lib/crc14.cpp: boost::augmented_crc<14, 0x2757>
// over length bytes, the last 14 bits being the zeros the CRC goes in.
int crc14(const unsigned char* data, int length);

// encode174_91(message77, codeword): 77 message bits, their CRC14 (by
// crc14() over 12 bytes, as the Fortran), to 174.
void encode174_91(const int8_t message77[77], int8_t codeword[174]);

// encode174_91_nocrc(message, codeword): 91 bits (CRC included) to 174.
void encode174_91_nocrc(const int8_t message[91], int8_t codeword[174]);

// decode174_91(llr, Keff, maxosd, norder, apmask, message91, cw, ntype,
// nharderror, dmin): belief propagation, then maxosd OSD runs on its saved
// sums (maxosd 0: one on llr; below 0: none). message91, cw, ntype,
// nharderror and dmin are written as the Fortran writes its arguments, and
// only then: dmin is left alone when decoding fails, message91 when belief
// propagation alone fails. Keff must be 91 (as every caller in WSJT-X
// passes); norder 0-2 (3 and more reach boxit91, not carried).
// Thread safe: no state but a generator matrix made once.
void decode174_91(const float llr[174], int Keff, int maxosd, int norder, const int8_t apmask[174],
                  int8_t message91[91], int8_t cw[174], int& ntype, int& nharderror, float& dmin);

// The line-for-line translation decode174_91 was optimised from, kept for
// the side-by-side test and as its fallback.
void decode174_91_reference(const float llr[174], int Keff, int maxosd, int norder, const int8_t apmask[174],
                            int8_t message91[91], int8_t cw[174], int& ntype, int& nharderror, float& dmin);

// tanhf as the reference has it: glibc 2.39's (sysdeps/ieee754/flt-32
// s_tanhf.c and the s_expm1f.c it calls, fdlibm's, no FMA in x86-64's
// build), line for line, so that every platform has the reference's tanh
// whatever its libm (glibc 2.41's tanhf, the Pi's, is CORE-MATH's correctly
// rounded one and differs in the last bit now and then; mingw's differs).
// Equal to glibc 2.39's tanhf on every float ([ft8-tanhf-libm], x86-64).
float tanhf_glibc(float x);

// y(i) = tanhf_glibc(x(i)), i < n: eight at a time on x86-64 with AVX2 and
// on aarch64 with NEON (bp_tanhf_vector() true), one by one elsewhere and
// for non-finite x. decode174_91's tanh(-toc/2).
void bp_tanhf(const float* x, float* y, int n);
bool bp_tanhf_vector();

// osd174_91(llr, k, apmask, ndeep, message91, cw, nhardmin, dmin), with the
// reference's results bit for bit (on bit-packed words): k 91 and ndeep 0-2
// as above; nhardmin negative when the CRC fails.
void osd174_91(const float llr[174], int k, const int8_t apmask[174], int ndeep, int8_t message91[91],
               int8_t cw[174], int& nhardmin, float& dmin);
// The line-for-line translation osd174_91 was optimised from, kept for the
// side-by-side test and as its fallback for non-finite llrs.
void osd174_91_reference(const float llr[174], int k, const int8_t apmask[174], int ndeep, int8_t message91[91],
                         int8_t cw[174], int& nhardmin, float& dmin);

}  // namespace ft8
