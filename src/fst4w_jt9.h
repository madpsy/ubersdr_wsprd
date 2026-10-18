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
#include <cstring>
#include <fstream>
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

// lib/readwav.f90: a RIFF WAV's header read and checked as it reads it
// (its checks in its order, its messages), then its 16-bit samples from the
// data chunk on demand. jt9's WAV reader for every mode's file mode.
struct Wav {
    std::ifstream f;
    int rate = 0, channels = 0, bits = 0;
    bool is_float = false;  // 32-bit IEEE float samples (not readwav.f90's)
    long remaining = 0;  // bytes of data left
    // allow_float: a 32-bit float file too (FST4W's file mode alone).
    bool open(const std::string& path, std::string& why, bool allow_float = false) {
        f.open(path, std::ios::binary);
        if (!f) {
            why = "cannot open file";
            return false;
        }
        f.seekg(0, std::ios::end);
        const int64_t filesize = static_cast<int64_t>(f.tellg());
        f.seekg(0);
        auto fail = [&](const char* text) {
            why = text;
            f.close();
            return false;
        };
        if (filesize < 12) return fail("file is too short for a RIFF/WAVE header");
        char riff[12];
        if (!f.read(riff, 12)) return fail("cannot read RIFF/WAVE header");
        if (std::memcmp(riff, "RIFF", 4) != 0 || std::memcmp(riff + 8, "WAVE", 4) != 0)
            return fail("not a RIFF/WAVE file");
        int64_t filepos = 13;  // Fortran's, 1-based
        bool format_found = false;
        for (;;) {
            if (filepos + 8 > filesize + 1) return fail("missing data chunk");
            char id[4];
            int32_t size;
            f.seekg(filepos - 1);
            if (!f.read(id, 4) || !f.read(reinterpret_cast<char*>(&size), 4)) return fail("cannot read WAV chunk header");
            filepos += 8;
            if (size < 0) return fail("WAV chunk has an invalid size");
            if (filepos + size > filesize + 1) return fail("WAV chunk extends beyond end of file");
            const int64_t nextpos = filepos + size + (size % 2);
            if (std::memcmp(id, "fmt ", 4) == 0) {
                if (size < 16) return fail("fmt chunk is too short");
                char b[40] = {};
                if (!f.read(b, size < 40 ? size : 40)) return fail("cannot read fmt chunk");
                int16_t fmt_tag, ch, bps;
                int32_t sr;
                std::memcpy(&fmt_tag, &b[0], 2);
                std::memcpy(&ch, &b[2], 2);
                std::memcpy(&sr, &b[4], 4);
                std::memcpy(&bps, &b[14], 2);
                rate = sr;
                channels = ch;
                bits = bps;
                uint16_t tag = static_cast<uint16_t>(fmt_tag);
                if (tag == 0xFFFE && size >= 40) std::memcpy(&tag, &b[24], 2);  // the subformat GUID's first two bytes
                is_float = allow_float && tag == 3 && bps == 32 && ch == 1;
                if (!is_float && (fmt_tag != 1 || ch != 1 || bps != 16)) return fail("unsupported WAV format; expected mono 16-bit PCM");
                if (sr != 11025 && sr != 12000) return fail("unsupported WAV sample rate");
                format_found = true;
            } else if (std::memcmp(id, "data", 4) == 0) {
                if (!format_found) return fail("data chunk appears before fmt chunk");
                remaining = size;
                f.seekg(filepos - 1);
                return true;
            }
            if (nextpos > filesize + 1) return fail("WAV chunk extends beyond end of file");
            filepos = nextpos;
        }
    }
    // Up to n float samples, scaled to the 16-bit samples' full scale
    // (x 32768); how many were read.
    long read(float* out, long n) {
        const long want = std::min(n, remaining / 4);
        f.read(reinterpret_cast<char*>(out), want * 4);
        const long got = static_cast<long>(f.gcount()) / 4;
        for (long i = 0; i < got; ++i) out[i] *= 32768.0f;
        remaining -= got * 4;
        return got;
    }
    // Up to n samples; how many were read.
    long read(int16_t* out, long n) {
        const long want = std::min(n, remaining / 2);
        f.read(reinterpret_cast<char*>(out), want * 2);
        const long got = static_cast<long>(f.gcount()) / 2;
        remaining -= got * 2;
        return got;
    }
};


// jt9's parse_wav_filename_nutc: the last 6 characters of the name before
// .wav/.WAV if all digits, else '_' and the last 4, else 0.
int jt9_nutc(const std::string& filename);

// Fortran's Fw.d and Iw editing as gfortran writes them (asterisks when the
// value does not fit; Fw.0 ends in a point). For the line formats.
std::string fortran_f(float x, int w, int d);
std::string fortran_i(long v, int w, int m = 0);

}  // namespace fst4w
