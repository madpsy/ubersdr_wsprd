// jt9's file mode for FT2, around ft2::Decoder (ft2.h): what
// `jt9 --ft2 [options] FILE.wav ...` prints and writes, byte for byte with
// the FT2 reference (WSJT-X improved 3.2.0's FT2-own code on
// 567ad29ce's FT4, tests/data/ft2/README.md). improved's
// jt9.f90 runs FT2 as mode 52: each WAV read in jt9's blocks of 3456
// samples, npts = 21*1728, so 10 blocks, 34,560 samples (2.88 s), whatever
// -p says; " EOF on input file" when it ends first; one multimode_decoder
// call on them (the FT2 branch stretches the first 39,744, zeros past what
// was read); ft2_decoded's lines (" * "), decoded.txt in the temp dir
// (" FT2"), <DecodeFinished>. Several files in one run share one decoder,
// as one jt9 process does (ft4_decode's saved state, packjt77's hashes).
//
// FT2 is --ft2 only: improved's getopt string has no '2', so "-2" is an
// unrecognised option (the usage, status 0), as there (jt9_cli.cpp).
// The options are FT8's (ft8::FileMode, read by jt9_cli.h); FT2 reads
// -d, -L, -H, -f, -Q, -X, -c, -x, -q, -t (and -p, -F, -g, -G only as jt9
// checks them), mode 5's parameters (jt9_params_init) but the mode.
//
// jt9.f90, jt9_params_init.f90, decoder.f90, decoder_callbacks.f90: Joe
// Taylor, K1JT, and the WSJT-X developers (GNU GPL v3); jt9 --ft2 from
// WSJT-X improved 3.2.0 (Uwe Risse, DG2YCB, and the same); and so is this.
#pragma once

#include <ostream>
#include <string>
#include <vector>

#include "ft8_jt9.h"

namespace ft2 {

// Each file decoded in turn by one Decoder, as jt9 --ft2 does (cwd and the
// return as ft8::run_files).
int run_files(const ft8::FileMode& fm, const std::vector<std::string>& files, const std::string& cwd,
              std::ostream& out, std::ostream& err);

// Ours, --feed: as ft8::run_feed, for FT2.
int run_feed(const ft8::FileMode& fm, std::istream& in, const std::string& cwd, std::ostream& out, std::ostream& err);

}  // namespace ft2
