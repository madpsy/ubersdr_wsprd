// jt9 -s in WebAssembly: there is no SysV shared memory to attach to.
#include "jt9_shm.h"

namespace jt9_shm {

int run(const Options& o, const std::string&, std::ostream&, std::ostream& err) {
    err << o.program << ": -s (shared memory) is not available in this build\n";
    return 1;
}

}  // namespace jt9_shm
