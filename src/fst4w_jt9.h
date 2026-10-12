// jt9's file mode for FST4W, around fst4w::Decoder (fst4w.h): what
// `jt9 -p P --fst4w ... FILE.wav ...` prints and writes, byte for byte
// (WSJT-X lib/jt9.f90, decoder.f90, decoder_callbacks.f90 and
// fst4_decode.f90's output at 567ad29ce): each WAV read as jt9 reads it,
// the decoder's RMS check, the decode lines on stdout, decoded.txt,
// fst4_decodes.dat, fort.52 (dopspread's spectra) and fst4w_calls.txt, and
// <DecodeFinished>. Several files in one run share one decoder, as one jt9
// process does (its hash table and calls list).
//
// The command line is jt9's (jt9_cli.h): run_files() takes it read.
//
// jt9.f90, decoder.f90, decoder_callbacks.f90: Joe Taylor, K1JT, and the
// WSJT-X developers (GNU GPL v3), and so is this.
#pragma once

#include <cstddef>
#include <cstdint>
#include <ostream>
#include <string>
#include <vector>

#include "fst4w.h"

namespace fst4w {

// jt9's file mode for FST4W, its command line read (jt9_cli.h does that):
// each file decoded in turn by one Decoder (one hash table and calls list
// for the run), as jt9 -p P --fst4w (or -Y) ... FILE... does. cwd is the
// directory jt9 runs in (plotspec is looked for there, fort.52 written
// there, and relative paths start there). Returns jt9's exit status.
struct FileMode {
    std::string program = "jt9";  // in messages: "jt9: cannot read WAV file ..."
    std::string data_dir = ".", temp_dir = ".";  // -a, -t
    double trperiod = 60.0;                      // -p (FST4W: 15, 30, 60, 120, 300, 900 or 1800)
    bool hash22 = false;                         // -Y: --fst4w printing hash22 values
    int ndepth = 1, nrxfreq = 1500, ntol = 20;   // -d, -f, -F
    int nexp_decode = 0;                         // -X
    bool nexp_decode_set = false;
    bool quiet = false;                          // -q
    std::string mycall = "K1ABC", hiscall = "W9XYZ";  // -c, -x as read ("b": none)
    Options::Spread spread = Options::Spread::Fast;   // --spread
};
int run_files(const FileMode& fm, const std::vector<std::string>& files, const std::string& cwd, std::ostream& out,
              std::ostream& err);

// A period's samples as jt9's file mode reads a WAV holding them: a buffer
// of trperiod*12000 (what the decoder reads of jt9's id2), filled in
// blocks of 3456 (so a part block at the end is left out: FST4W-120 reads
// 1,437,696 samples) and zero after; true when it passes
// multimode_decoder_core's check (the first 15 s not silent: RMS 0.5 or
// more), else jt9 does not decode it.
bool jt9_samples(const int16_t* samples, size_t n, int trperiod, std::vector<int16_t>& id2);

// A decode's fst4_decodes.dat line (with decdata), no newline.
std::string decodes_dat_line(const Result& r);

// fst4w_calls.txt as jt9 writes it ('(a20)' each: right-justified).
std::string calls_file(const std::vector<std::string>& calls);

// Fortran's Fw.d and Iw editing as gfortran writes them (asterisks when the
// value does not fit; Fw.0 ends in a point). For the line formats.
std::string fortran_f(float x, int w, int d);
std::string fortran_i(long v, int w, int m = 0);

}  // namespace fst4w
