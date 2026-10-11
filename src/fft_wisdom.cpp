// fft.h's wisdom file: wsprd keeps none. Its one plan of fft.h's (each
// decode's spread) is quick, FFTW's estimate; the decoder's own plans use
// wsprd's wspr_wisdom.dat, as WSJT-X's wsprd does.
#include "fft.h"

namespace skimmer {

void fft_load_wisdom() {}
void fft_save_wisdom() {}

}  // namespace skimmer
