// wsprd: WSJT-X's wsprd command line, files and output, through the
// in-process decoder (wsprd.h's run()), a drop-in for the reference:
//
//   wsprd [-a dir] [-BcdHJmnqstvw0] [-C n] [-e x] [-f MHz] [-o n] [-z x] file.wav|file.c2
//
// It also takes -n (Io::with_n): no
// coarse drift search, and each ALL_WSPR.TXT line ending in the decode's
// spread; its usage has -n's line. Without -n it is wsprd exactly. And
// "wsprd --version" (alone) prints the release it was built as, where
// wsprd's getopt would refuse "--version".
//
// WSPRD_TRACE=1 prints each stage to stderr as the oracle's instrumented
// wsprd does (tests/data/wsprd/README.md).
#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

#include "wsprd.h"

#ifndef WSPRD_VERSION
#define WSPRD_VERSION "dev"
#endif

int main(int argc, char** argv) {
    if (argc == 2 && std::string(argv[1]) == "--version") {
        std::printf("wsprd %s (WSJT-X 567ad29ce's decoder)\n", WSPRD_VERSION);
        return 0;
    }
    std::vector<std::string> args(argv, argv + argc);
    wsprd::Io io;
    io.out = [](const std::string& s) { std::fwrite(s.data(), 1, s.size(), stdout); };
    io.err = [](const std::string& s) { std::fwrite(s.data(), 1, s.size(), stderr); };
    io.flush = [] { std::fflush(stdout); };
    io.in = stdin;
    // as wsprd's SIGSEGV: what stdio holds is lost, as it buffers it here
    io.crash = [] { std::_Exit(139); };
    io.with_n = true;
    return wsprd::run(args, io);
}
