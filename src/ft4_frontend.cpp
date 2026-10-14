// FT4's front end (ft4_frontend.h): the Frontend's construction and size;
// each routine group's state and code are in its own file
// (ft4_candidates.cpp, ft4_sync.cpp, ft4_metrics.cpp, ft4_wave.cpp).
//
// Joe Taylor, K1JT, Steven Franke, K9AN, and the WSJT-X developers (GNU GPL
// v3), and so is this.
#include "ft4_frontend.h"

namespace ft4 {

Frontend::Frontend()
    : cand_(make_cand()), sync_(make_sync()), metrics_(make_metrics()), wave_(make_wave()) {}

Frontend::~Frontend() {
    del_wave(wave_);
    del_metrics(metrics_);
    del_sync(sync_);
    del_cand(cand_);
}

size_t Frontend::heap_bytes() const {
    return cand_bytes(cand_) + sync_bytes(sync_) + metrics_bytes(metrics_) + wave_bytes(wave_);
}

}  // namespace ft4
