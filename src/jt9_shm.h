// jt9 -s KEY: our jt9 as a back end over shared memory, as WSJT-X 2.7.0's
// jt9 runs under WSJT-X and under jt9-decode (lib/jt9a.f90 at v2.7.0, the
// classic protocol: ipc[1] 1 to decode, 999 to quit, ipc[2] 1 to go on;
// not 567ad29ce's generation-bound one). FT8, FT4 and FT2, by the same
// decoders as jt9's file mode (ft8.h, ft4.h, ft2.h), their parameters read
// from the memory's params as decoder.f90 reads them; the memory is Qt 5's
// QSharedMemory (qt_shm.h: System V on Linux, named objects on Windows),
// laid out as jt9_commons.h has it. What it prints is what jt9 prints: each
// decode's line and <DecodeFinished>, and decoded.txt in the temp dir.
//
// A request this jt9 cannot decode as jt9 would (another mode, FT8's
// multithreaded decoder, its Fox or SuperFox mode, an EME delay) is said
// once on err and answered with no decodes, so the caller is never left
// waiting.
#pragma once

#include <ostream>
#include <string>

namespace jt9_shm {

struct Options {
    std::string key;       // -s
    std::string temp_dir;  // -t: decoded.txt's
    std::string program = "jt9";
};

// Attaches to the memory and decodes each request until told to quit (0).
// The memory not there: said, and abort() as jt9a's "call abort". Memory
// smaller than jt9_commons.h's layout: said, status 2.
int run(const Options& o, const std::string& cwd, std::ostream& out, std::ostream& err);

}  // namespace jt9_shm
