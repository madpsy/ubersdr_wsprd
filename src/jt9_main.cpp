// jt9: WSJT-X's jt9 run on WAV files, a mode at a time (jt9_cli.h; FST4W so
// far), its command line, output and files as jt9's.
#include <iostream>
#include <string>
#include <vector>

#include "jt9_cli.h"

#ifndef JT9_VERSION
#define JT9_VERSION "dev"
#endif

int main(int argc, char** argv) {
    std::ios::sync_with_stdio(false);
    const std::vector<std::string> args(argv + 1, argv + argc);
    const int status = jt9_cli::run(args, ".", std::cout, std::cerr, JT9_VERSION);
    std::cout.flush();
    return status;
}
