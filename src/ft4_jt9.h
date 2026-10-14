// jt9's file mode for FT4, around ft4::Decoder (ft4.h): what
// `jt9 -5 [options] FILE.wav ...` prints and writes, byte for byte (WSJT-X
// lib/jt9.f90's WAV loop, decoder.f90's FT4 branch, decoder_callbacks.f90's
// ft4_decoded, at 567ad29ce): each WAV read as jt9 reads it for FT4 (21
// blocks of 3456 samples, 72,576, whatever -p says; " EOF on input file"
// when it ends first), one multimode_decoder call on them, the decode
// lines, decoded.txt in the temp dir as gfortran leaves it, and
// <DecodeFinished>. Several files in one run share one decoder, as one jt9
// process does (ft4_decode's saved state, packjt77's hashes).
//
// The options are FT8's (ft8::FileMode, read by jt9_cli.h); FT4 reads
// -d, -L, -H, -f, -Q, -X, -c, -x, -q, -t (and -p, -F, -g, -G only as jt9
// checks them).
//
// jt9.f90, jt9_params_init.f90, decoder.f90, decoder_callbacks.f90: Joe
// Taylor, K1JT, and the WSJT-X developers (GNU GPL v3), and so is this.
#pragma once

#include <ostream>
#include <string>
#include <vector>

#include "ft8_jt9.h"

namespace ft4 {

// Each file decoded in turn by one Decoder, as jt9 -5 does (cwd and the
// return as ft8::run_files).
int run_files(const ft8::FileMode& fm, const std::vector<std::string>& files, const std::string& cwd,
              std::ostream& out, std::ostream& err);

}  // namespace ft4
