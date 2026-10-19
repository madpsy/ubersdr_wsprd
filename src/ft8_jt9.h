// jt9's file mode for FT8, around ft8::Decoder (ft8.h): what
// `jt9 -8 [options] FILE.wav ...` prints and writes, byte for byte (WSJT-X
// lib/jt9.f90's WAV loop and FT8 branch, decoder.f90, decoder_callbacks.f90's
// ft8_decoded, at 567ad29ce): each WAV read as jt9 reads it (blocks of
// 3456 samples, npts = TRperiod*12000 of them at most, never past 179,712;
// " EOF on input file" when it ends first), the three multimode_decoder
// calls of jt9's FT8 file mode (nzhsym 41, 47 and 50, the last on 172,800
// samples), the decode lines, decoded.txt in the temp dir as gfortran leaves
// it, and <DecodeFinished>. Several files in one run share one decoder, as
// one jt9 process does (the a7 tables, packjt77's hashes, ft8_decode's
// saved state).
//
// The command line is jt9's (jt9_cli.h): run_files() takes it read.
//
// jt9.f90, jt9_params_init.f90, decoder.f90, decoder_callbacks.f90: Joe
// Taylor, K1JT, and the WSJT-X developers (GNU GPL v3), and so is this.
#pragma once

#include <istream>
#include <ostream>
#include <string>
#include <vector>

namespace ft8 {

// jt9's options FT8 reads, as jt9.f90 holds them (its defaults).
struct FileMode {
    std::string program = "jt9";                  // in messages
    std::string data_dir = ".", temp_dir = ".";   // -a, -t (decoded.txt goes to the temp dir)
    double trperiod = 60.0;                       // -p
    int ndepth = 1;                               // -d
    int nrxfreq = 1500;                           // -f (nfqso)
    int flow = 200, fhigh = 4000;                 // -L, -H (nfa, nfb)
    int nqso_progress = 0;                        // -Q (0-5)
    int nexp_decode = 0;                          // -X (ncontest = its low 3 bits)
    bool quiet = false;                           // -q: no <DecodeFinished>
    bool fast_subtract = false;                   // ours, --fast-subtract (ft8.h's Params; FT4's, FT2's too)
    std::string mycall = "K1ABC", hiscall = "W9XYZ", hisgrid = "EN37";  // -c, -x, -g as read ("b": none)
};

// Each file decoded in turn by one Decoder, as jt9 -8 does. cwd is the
// directory jt9 runs in (relative paths and dirs start there). Returns
// jt9's exit status (2 with its message and "STOP 2" for a file it cannot
// read or of another rate).
int run_files(const FileMode& fm, const std::vector<std::string>& files, const std::string& cwd, std::ostream& out,
              std::ostream& err);

// Ours, --feed: the files named on in, one a line ("PATH" or "PATH HHMMSS"),
// each decoded as it arrives, as run_files decodes its files (one Decoder for
// all: the hashes and a7 tables carried from period to period, as one jt9 -s
// carries them for WSJT-X), its lines and <DecodeFinished> printed, then
// "<FileDone> PATH N" (N its decodes; 0 for a file that cannot be read, said
// on err, the feed going on) and out flushed. The period's UTC: HHMMSS when
// given, else jt9's from the name (..._HHMMSS.wav), else ka9q-radio
// pcmrecord's (YYYYMMDDTHHMMSSZ_...), else 0. Ends at the end of in (0).
int run_feed(const FileMode& fm, std::istream& in, const std::string& cwd, std::ostream& out, std::ostream& err);

// --feed's: a line's path and the UTC it gives (feed_nutc's rules), and
// pcmrecord's name's HHMMSS (-1 when the name is not one).
bool feed_line(const std::string& line, std::string& path, int& nutc);
int pcmrecord_nutc(const std::string& path);

// Shared with FT4's file mode (ft4_jt9.h). A path as jt9 opens it from cwd
// (absolute ones as they are), dir/name, and decoded.txt as decoder.f90
// leaves it after a period: absent until a call passes the RMS gate (exists),
// then created, its lines those of the last call that wrote any.
std::string jt9_join(const std::string& dir, const std::string& name);
std::string jt9_at(const std::string& cwd, const std::string& path);
void jt9_write_decoded(bool exists, const std::string& text, const std::string& path);

}  // namespace ft8
