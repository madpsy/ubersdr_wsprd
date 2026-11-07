// jt9's file mode for MSK144: what `jt9 -k [options] FILE.wav ...` (or
// --msk144) prints and writes, byte for byte (WSJT-X lib/jt9.f90's WAV loop,
// decode_msk144.f90, decoder.f90's one multimode_decoder call after it, at
// 567ad29ce): each WAV read as jt9 reads it for MSK144 (blocks of 3456
// samples, TRperiod*12000 of them at most and never past the block that
// reaches nhsym 181; " EOF on input file" when it ends first), mskrtd run on
// every 7168 samples at half-block steps, the decode lines (a80, on stdout:
// MSK144 writes none to decoded.txt), the empty decoded.txt that
// multimode_decoder opens once the first 180,000 samples have an rms of 0.5,
// and <DecodeFinished> (always 0 0 0: nothing of MSK144 is counted).
//
// WSJT-X's jt9 runs one process per file, so mskrtd's saved state and
// packjt77's tables never cross files. Several files in one run here (and
// --feed) each start afresh, as that many processes would.
//
// The options are FT8's (ft8::FileMode, read by jt9_cli.h); MSK144 reads
// -p, -f, -F, -d, -c, -x, -t, -q (and -g, -G, -L, -H, -Q, -X only as jt9
// checks them), and ours alone: --msk-shmsg (MSK40 short messages, WSJT-X's
// bshmsg) and --msk-swl (its bswl), both off by default, as they are in jt9.
//
// jt9.f90, jt9_params_init.f90, decode_msk144.f90, mskrtd.f90: Joe Taylor,
// K1JT, and the WSJT-X developers (GNU GPL v3), and so is this.
#pragma once

#include <istream>
#include <ostream>
#include <string>
#include <vector>

#include "ft8_jt9.h"

namespace msk144 {

// Each file decoded in turn, as one jt9 process per file would (cwd and the
// return as ft8::run_files).
int run_files(const ft8::FileMode& fm, const std::vector<std::string>& files, const std::string& cwd,
              std::ostream& out, std::ostream& err);

// Ours, --feed: as ft8::run_feed, for MSK144. "<FileDone> PATH N" gives N
// the decode lines the file printed.
int run_feed(const ft8::FileMode& fm, std::istream& in, const std::string& cwd, std::ostream& out, std::ostream& err);

}  // namespace msk144
