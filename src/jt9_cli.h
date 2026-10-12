// jt9: our own build of WSJT-X's jt9 program, named as WSJT-X names it,
// run on WAV files: a mode at a time, each decoding as WSJT-X's jt9 does (at
// 567ad29ce), its command line, files and output the same, so it can stand
// in for WSJT-X's jt9 file for file. So far FST4W (--fst4w, -W, or -Y:
// fst4w_jt9.h, fst4w.h), with the options FST4W reads; other modes and
// their options are not taken yet (unknown to it, as WSJT-X's jt9 treats an
// option it does not know).
//
//   jt9 [OPTIONS] file1 [file2 ...]
//
// The command line is read as WSJT-X's jt9 reads it (lib/options.f90's
// getopt and jt9.f90's checks), its quirks too: an option and its value are
// separate words ("-p 120", never "-p120"), no value may start with '-', the
// first word that is not an option ends them (the rest are files), a long
// option taking a value takes the word after its value instead (getopt
// checks it twice), an unknown option is said and the usage printed (status
// 0), and -c and -x are read as Fortran reads a list ("-c PJ4/K1ABC" is
// PJ4). FST4W on a period fst4_decode has no sizes for (not 15, 30, 60, 120,
// 300, 900 or 1800 s) ends as WSJT-X's does, of SIGFPE: gfortran's line on
// stderr, status 136 (no backtrace).
//
// Ours alone: --spread fast|exact|off (FST4W's dopspread: fast by default,
// exact for WSJT-X's byte for byte; a long option alone, taking its one
// value); -v's line. Still apart: -c and -x only standard calls (WSJT-X's
// takes any: fst4w_unpack.h); no timer.out or jt9_wisdom.dat written (its
// timings and FFTW wisdom).
//
// jt9.f90, options.f90, jt9_input_validation.f90: Joe Taylor, K1JT, and the
// WSJT-X developers (GNU GPL v3), and so is this.
#pragma once

#include <ostream>
#include <string>
#include <vector>

namespace jt9_cli {

// argv without the program name; cwd as fst4w::run_files() takes it.
// version: what -v prints after "jt9 ". Returns the exit status.
int run(const std::vector<std::string>& args, const std::string& cwd, std::ostream& out, std::ostream& err,
        const std::string& version = "dev");

}  // namespace jt9_cli
