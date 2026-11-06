// MSK144's forward error correction: WSJT-X's (128,90) LDPC belief-
// propagation decoder and its 13-bit CRC check (lib at 567ad29ce:
// bpdecode128_90.f90, ldpc_128_90_reordered_parity.f90, lib/ft8/chkcrc13a.f90,
// lib/crc13.cpp, lib/platanh.f90), giving the reference's results to the bit:
// the same message, codeword, hard-error count and iteration count from the
// same soft bits, and the same message77 left behind when decoding fails.
//
// What it rests on: no multiply-adds fused (-ffp-contract=off, CMakeLists.txt);
// float sums and products in the Fortran's order; glibc's tanhf (tanhf_glibc,
// ft8_ldpc.h), which gfortran calls for tanh of a default real.
//
// bpdecode128_90.f90, chkcrc13a.f90, crc13.cpp: Steven Franke, K9AN, and Joe
// Taylor, K1JT, WSJT-X (GNU GPL v3), and so is this.
#pragma once

namespace msk144 {

// platanh(x, y) (lib/platanh.f90): a piecewise-linear atanh, exported for
// scripts/msk144-piece-check.sh's sweep against the Fortran's.
float platanh(float x);

// A debug-only trace of bpdecode128_90, for the piece check: null in
// production (one untaken branch per iteration). When set, it is called at
// the top of each iteration, once zn is updated, with iter, zn(1:128) and
// tov(1:3,1:128), both 1-based in their last index (zn[0] and tov[*][0]
// unused), as the Fortran's arrays are.
using Bpdecode128_90Trace = void (*)(int iter, const float zn[129], const float tov[4][129]);
extern Bpdecode128_90Trace bpdecode128_90_trace;

// bpdecode128_90(llr, apmask, maxiterations, message77, cw, nharderror, iter),
// with apmask all zero, as msk144decodeframe passes it. message77 is written
// only on success (it keeps its value otherwise, as the Fortran's does); cw is
// the last hard decision; nharderror is -1 on failure; niterations is `iter`
// as the Fortran leaves it (maxiterations+1 when the loop runs out). Bits are
// 0/1. The return value equals nharderror.
int bpdecode128_90(const float llr[128], int maxiterations, signed char message77[77], signed char cw[128],
                   int& nharderror, int& niterations);

// chkcrc13a(decoded, nbadcrc) with nbadcrc == 0: true when the 13 bits after
// the 77 of decoded(1:90) are the CRC of the 77 (bits 0/1).
bool chkcrc13a(const signed char decoded[90]);

// encode_128_90(message77, codeword) (lib/encode_128_90.f90, with
// lib/ldpc_128_90_generator.f90): the 77 bits (0/1), the 13-bit CRC of them
// and the 38 parity checks, the 128-bit codeword. Added with the frame
// decoder (msk144_decode.h), whose signal-quality step regenerates the
// transmitted bits; checked there (scripts/msk144-frame-oracle.f90).
void encode_128_90(const signed char message77[77], signed char codeword[128]);

}  // namespace msk144
