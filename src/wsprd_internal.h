// What wsprd.cpp, wsprd_unpack.cpp, wsprd_fano.cpp and wsprd_files.cpp share
// of WSJT-X's wsprd (lib/wsprd at 567ad29ce): its message packing and
// unpacking, hash, Fano and stack decoders, under their own names. Not for
// use outside the wsprd files: wsprd.h is the interface.
#pragma once

#include <cstddef>
#include <cstdint>
#include <ctime>

namespace wsprd::detail {

// Thrown where wsprd would crash mid-decode; Decoder sets status 139.
struct Crash {};

// wsprd's clock() for its timers. wsprd is a process of one thread, so its
// clock() is that thread's CPU time: here the calling thread's, in clock()'s
// ticks (as glibc's clock() makes them), so other threads in the process do
// not count. On Windows clock() is the wall time, which threads do not add
// to: clock() itself.
clock_t cpu_clock();

// wsprd_utils.h
enum {
    WSPRD_CALLSIGN_SIZE = 13,
    WSPRD_GRID4_SIZE = 5,
    WSPRD_GRID6_SIZE = 7,
    WSPRD_MESSAGE_SIZE = 23,
    WSPRD_HASH_COUNT = 32768
};

// wsprd_unpack.cpp: nhash.c (over lookup3.c's wsjt_hashlittle),
// wsprd_utils.c and wsprsim_utils.c.
uint32_t nhash(const void* key, size_t length, uint32_t initval);
void unpack50(signed char* dat, int32_t* n1, int32_t* n2);
int unpackcall(int32_t ncall, char* call);
int unpackgrid(int32_t ngrid, char* grid);
int unpackpfx(int32_t nprefix, char* call);
void deinterleave(unsigned char* sym);
int floatcomp(const void* elem1, const void* elem2);
int unpk_(signed char* message, char* hashtab, char* loctab, char* call_loc_pow, char* callsign);
int wsprd_load_hash_line(char const* line, char* hashtab, char* loctab);
int get_wspr_channel_symbols(char* message, char* hashtab, char* loctab, unsigned char* symbols);

// wsprd_fano.cpp: fano.c, jelinek.c, tab.c and metric_tables.c.
struct snode {
    uint64_t encstate;     // Encoder state
    int gamma;             // Cumulative metric to this node
    unsigned int depth;    // depth of this node
    unsigned int jpointer;
};
int fano(unsigned int* metric, unsigned int* cycles, unsigned int* maxnp, unsigned char* data, unsigned char* symbols,
         unsigned int nbits, int mettab[2][256], int delta, unsigned int maxcycles);
int encode(unsigned char* symbols, unsigned char* data, unsigned int nbytes);
int jelinek(unsigned int* metric, unsigned int* cycles, unsigned char* data, unsigned char* symbols, unsigned int nbits,
            unsigned int stacksize, struct snode* stack, int mettab[2][256], unsigned int maxcycles);
extern const float metric_tables[5][256];

}  // namespace wsprd::detail
