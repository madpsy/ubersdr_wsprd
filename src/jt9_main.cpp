// jt9: WSJT-X's jt9 run on WAV files, a mode at a time (jt9_cli.h; FST4W,
// FT8, FT4 and FT2 so far), its command line, output and files as jt9's;
// or with -s as a back end over shared memory (jt9_shm.h).
#include <sys/stat.h>

#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "jt9_cli.h"

#ifndef JT9_VERSION
#define JT9_VERSION "dev"
#endif

int main(int argc, char** argv) {
    std::ios::sync_with_stdio(false);
    const std::vector<std::string> args(argv + 1, argv + argc);
    std::ostringstream err;
    const int status = jt9_cli::run(args, ".", std::cout, err, JT9_VERSION, &std::cerr);
    std::cout.flush();
    // gfortran's order on stderr: a regular file gets "STOP 2" first.
    struct stat st;
    const bool file = fstat(2, &st) == 0 && S_ISREG(st.st_mode);
    std::cerr << (file ? jt9_cli::stderr_as_file(err.str()) : err.str());
    std::cerr.flush();
    return status;
}
