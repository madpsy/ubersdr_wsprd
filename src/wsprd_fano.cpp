// wsprd's sequential decoders for the K=32 r=1/2 convolutional code, its
// encoder and Fano metric tables: WSJT-X's lib/wsprd/fano.c, jelinek.c,
// tab.c and metric_tables.c at 567ad29ce, for wsprd.cpp (see wsprd.h).
//
// fano.c, fano.h: Copyright 1994, Phil Karn, KA9Q, minor modifications by
// Joe Taylor, K1JT. jelinek.c: Steve Franke, K9AN, for WSJT-X (July 2015).
// WSJT-X is GNU GPL v3, and so is this.
//
// The reference's code, kept as it was. What had to change, and nothing
// else: C++ casts on malloc, fano()'s nodes calloc'd (see there),
// jelinek()'s buckets a growing array (see there; wsprd leaks them),
// namespace wsprd::detail, the metric tables a const array at file scope
// (wsprd.c includes metric_tables.c inside main), and the NASA_STANDARD/MJ
// codes (fano.c builds LL) left out. unsigned long
// and long stay as they were: 32 bits on Windows, 64 elsewhere, as in the
// reference built for each. For speed, with the same results: ENCODE takes
// the compiler's parity rather than Partab's (on arm64 a fold by shifts),
// and fano() keeps *maxnp in a local (the int stores to its nodes made the
// compiler reload it), as the deepest node, updated only where np moves
// forward, keeps the current node's gamma, tm[i] and encstate in locals as
// well, and takes a first threshold step by a select. A node's two
// metrics are sorted once for each of the four syms its 0-branch can have,
// before the search, rather than on every visit (the bit that says the
// 1-branch is the better in one word, so encstate does not wait on a load
// that waits on it). And fano() does not ENCODE: a node keeps its
// encstate's parities (par), which a step forward shifts and a branch
// flips, so the next lsym is two bits of a shift, not two 32-bit parities
// on the encstate chain.
#include "wsprd_internal.h"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <vector>

namespace wsprd::detail {

namespace {

// tab.c: 8-bit parity lookup table.
const unsigned char Partab[] = {
 0, 1, 1, 0, 1, 0, 0, 1,
 1, 0, 0, 1, 0, 1, 1, 0,
 1, 0, 0, 1, 0, 1, 1, 0,
 0, 1, 1, 0, 1, 0, 0, 1,
 1, 0, 0, 1, 0, 1, 1, 0,
 0, 1, 1, 0, 1, 0, 0, 1,
 0, 1, 1, 0, 1, 0, 0, 1,
 1, 0, 0, 1, 0, 1, 1, 0,
 1, 0, 0, 1, 0, 1, 1, 0,
 0, 1, 1, 0, 1, 0, 0, 1,
 0, 1, 1, 0, 1, 0, 0, 1,
 1, 0, 0, 1, 0, 1, 1, 0,
 0, 1, 1, 0, 1, 0, 0, 1,
 1, 0, 0, 1, 0, 1, 1, 0,
 1, 0, 0, 1, 0, 1, 1, 0,
 0, 1, 1, 0, 1, 0, 0, 1,
 1, 0, 0, 1, 0, 1, 1, 0,
 0, 1, 1, 0, 1, 0, 0, 1,
 0, 1, 1, 0, 1, 0, 0, 1,
 1, 0, 0, 1, 0, 1, 1, 0,
 0, 1, 1, 0, 1, 0, 0, 1,
 1, 0, 0, 1, 0, 1, 1, 0,
 1, 0, 0, 1, 0, 1, 1, 0,
 0, 1, 1, 0, 1, 0, 0, 1,
 0, 1, 1, 0, 1, 0, 0, 1,
 1, 0, 0, 1, 0, 1, 1, 0,
 1, 0, 0, 1, 0, 1, 1, 0,
 0, 1, 1, 0, 1, 0, 0, 1,
 1, 0, 0, 1, 0, 1, 1, 0,
 0, 1, 1, 0, 1, 0, 0, 1,
 0, 1, 1, 0, 1, 0, 0, 1,
 1, 0, 0, 1, 0, 1, 1, 0,
};

// Layland-Lushbaugh code
// Nonsystematic, non-quick look-in, dmin=?, dfree=?
#define	POLY1	0xf2d05351
#define	POLY2	0xe4613c47

/* Convolutional encoder macro. Takes the encoder state, generates
 * a rate 1/2 symbol pair and stores it in 'sym'. The symbol generated from
 * POLY1 goes into the 2-bit of sym, and the symbol generated from POLY2
 * goes into the 1-bit.
 */
#if defined(__GNUC__) && defined(__aarch64__)
// The same two parities by shifts and XORs, side by side in one word
// (POLY1's 32 bits above POLY2's): arm64 has no scalar popcount, and the
// compiler's parity goes through a vector register and back.
#define	ENCODE(sym,encstate){\
unsigned long _w = (unsigned long)((encstate) & POLY1) << 32 | ((encstate) & POLY2);\
_w ^= _w >> 16; _w ^= _w >> 8; _w ^= _w >> 4; _w ^= _w >> 2; _w ^= _w >> 1;\
(sym) = (_w >> 31 & 2) | (_w & 1);\
}
#elif defined(__GNUC__)
// The parity of the polynomial's 32 bits, which is what Partab gives of
// them folded to a byte.
#define	ENCODE(sym,encstate){\
(sym) = __builtin_parityl((encstate) & POLY1) << 1;\
(sym) |= __builtin_parityl((encstate) & POLY2);\
}
#else
#define	ENCODE(sym,encstate){\
unsigned long _tmp;\
\
_tmp = (encstate) & POLY1;\
_tmp ^= _tmp >> 16;\
(sym) = Partab[(_tmp ^ (_tmp >> 8)) & 0xff] << 1;\
_tmp = (encstate) & POLY2;\
_tmp ^= _tmp >> 16;\
(sym) |= Partab[(_tmp ^ (_tmp >> 8)) & 0xff];\
}
#endif

// fano()'s par for an encstate: bit 2k+1 the parity of encstate & (POLY1>>k),
// bit 2k that of encstate & (POLY2>>k), k=0..31, so bits 1,0 are ENCODE's
// sym. encstate<<1 makes it par>>2 (k=32 masks nothing), and encstate^1
// flips bit k of each where the polynomial has it: par^PARSTEP, PARSTEP
// being POLY1's bit k at bit 2k+1 and POLY2's at bit 2k.
constexpr uint64_t interleave(uint64_t p1, uint64_t p2) {
  uint64_t r = 0;
  for (int k = 0; k < 32; k++) r |= (p1 >> k & 1) << (2*k+1) | (p2 >> k & 1) << (2*k);
  return r;
}
constexpr uint64_t PARSTEP = interleave(POLY1, POLY2);

struct node {
  unsigned long encstate;	// Encoder state of next node
  uint64_t par;			// encstate's parities (see PARSTEP)
  long gamma;		        // Cumulative metric to this node
  int best[4], worst[4];	// tm[0], tm[1] by the 0-branch's sym
  unsigned int up;		// bit lsym: whether the 1-branch is the better
  int tm[2];		        // Sorted metrics for current hypotheses
  int i;			// Current branch being tested
};

}  // namespace

// ------------------------------------------------------------------- fano.c

/* Convolutionally encode a packet. The input data bytes are read
 * high bit first and the encoded packet is written into 'symbols',
 * one symbol per byte. The first symbol is generated from POLY1,
 * the second from POLY2.
 *
 * Storing only one symbol per byte uses more space, but it is faster
 * and easier than trying to pack them more compactly.
 */
int encode(
	   unsigned char *symbols,	// Output buffer, 2*nbytes*8
	   unsigned char *data,		// Input buffer, nbytes
	   unsigned int nbytes)		// Number of bytes in data
{
  unsigned long encstate;
  int sym;
  int i;

  encstate = 0;
  while(nbytes-- != 0) {
    for(i=7;i>=0;i--) {
      encstate = (encstate << 1) | ((*data >> i) & 1);
      ENCODE(sym,encstate);
      *symbols++ = sym >> 1;
      *symbols++ = sym & 1;
    }
    data++;
  }
  return 0;
}

/* Decode packet with the Fano algorithm.
 * Return 0 on success, -1 on timeout
 */
int fano(
	 unsigned int  *metric,	   // Final path metric (returned value)
	 unsigned int  *cycles,	   // Cycle count (returned value)
	 unsigned int  *maxnp,     // Progress before timeout (returned value)
	 unsigned char *data,	   // Decoded output data
	 unsigned char *symbols,   // Raw deinterleaved input symbols
	 unsigned int nbits,	   // Number of output bits
	 int mettab[2][256],	   // Metric table, [sent sym][rx symbol]
	 int delta,		   // Threshold adjust parameter
	 unsigned int maxcycles)   // Decoding timeout in cycles per bit
{
  struct node *nodes;		   // First node
  struct node *np;	           // Current node
  struct node *lastnode;	   // Last node
  struct node *tail;		   // First node of tail
  int t;			   // Threshold
  int  m0,m1;
  int ngamma;
  unsigned int lsym;
  unsigned int i;

  // calloc, not malloc: with fewer cycles than bits (-C 0) the reference
  // returns nodes it never wrote, whatever its heap held; here zeros.
  if((nodes = (struct node *)calloc(nbits+1,sizeof(struct node))) == NULL) {
    printf("malloc failed\n");
    return 0;
  }
  lastnode = &nodes[nbits-1];
  tail = &nodes[nbits-31];
  *maxnp = 0;
  struct node *deep = nodes;  // the deepest np a cycle starts at (*maxnp)

/* Compute all possible branch metrics for each symbol pair
 * This is the only place we actually look at the raw input symbols
 */
  int metrics[4];
  for(np=nodes;np <= lastnode;np++) {
    metrics[0] = mettab[0][symbols[0]] + mettab[0][symbols[1]];
    metrics[1] = mettab[0][symbols[0]] + mettab[1][symbols[1]];
    metrics[2] = mettab[1][symbols[0]] + mettab[0][symbols[1]];
    metrics[3] = mettab[1][symbols[0]] + mettab[1][symbols[1]];
    symbols += 2;
    // A node's sorted metrics for each 0-branch sym lsym, as the loop
    // below sorts them: off the tail, if(m0 > m1) the 0-branch is the
    // better, else the 1-branch (the low bit set); on it, the 0-branch
    // only (tm[1] is never read there: i stays 0 on the tail).
    np->up = 0;
    for(lsym=0;lsym<4;lsym++) {
      m0 = metrics[lsym];
      m1 = metrics[3^lsym];
      if(np >= tail) {
        np->best[lsym] = m0;
        np->worst[lsym] = 0;
      } else {
        np->best[lsym] = m0 > m1 ? m0 : m1;
        np->worst[lsym] = m0 > m1 ? m1 : m0;
        np->up |= (unsigned int)(m0 <= m1) << lsym;
      }
    }
    if(np == nodes) {
      // The root: both branches, wherever the tail starts
      np->encstate = 0;
      ENCODE(lsym,np->encstate);	// 0-branch (LSB is 0)
      m0 = metrics[lsym];

/* Now do the 1-branch. To save another ENCODE call here and
 * inside the loop, we assume that both polynomials are odd,
 * providing complementary pairs of branch symbols.

 * This code should be modified if a systematic code were used.
 */

      m1 = metrics[3^lsym];
      if(m0 > m1) {
        np->tm[0] = m0;                         // 0-branch has better metric
        np->tm[1] = m1;
      } else {
        np->tm[0] = m1;                         // 1-branch is better
        np->tm[1] = m0;
        np->encstate++;	                        // Set low bit
      }
      np->par = np->encstate ? PARSTEP : 0;
    }
  }
  np = nodes;
  np->i = 0;	                                // Start with best branch
  maxcycles *= nbits;
  np->gamma = t = 0;

  // Start the Fano decoder. The current node's gamma, tm[i] and encstate
  // are kept in gam, tmc and es too (the nodes as the reference has them),
  // so a step forward does not wait on what the one before stored.
  long gam = 0;
  int tmc = np->tm[0];
  unsigned long es = np->encstate;
  uint64_t par = np->par;
  for(i=1;i <= maxcycles;i++) {
    // deep: the deepest np a cycle starts at, taken below where np
    // deepens (when another cycle follows)
#ifdef	debug
    printf("k=%ld, g=%ld, t=%d, m[%d]=%d, maxnp=%d, encstate=%lx\n",
	   np-nodes,np->gamma,t,np->i,np->tm[np->i],(int)(deep-nodes),np->encstate);
#endif
// Look forward */
    ngamma = gam + tmc;
    if(ngamma >= t) {
      if(gam < t + delta) {                     // Node is acceptable
	/* First time we've visited this node;
	 * Tighten threshold.
	 *
	 * This loop could be replaced with
	 *   t += delta * ((ngamma - t)/delta);
	 * but the multiply and divide are slower.
	 */
	t += ngamma >= t + delta ? delta : 0;     // mostly once at most: a select
	while(ngamma >= t + delta) t += delta;
      }
      np[1].gamma = gam = ngamma;               // Move forward
      np[1].encstate = es = es << 1;
      par >>= 2;
      if( ++np == (lastnode+1) ) {
	break;	                                // Done!
      }
      if(i < maxcycles && np > deep) deep=np;

      /* Compute and sort metrics, starting with the 
       * zero branch: from the node's table
       */
      lsym = par & 3;                           // ENCODE(lsym,es), see PARSTEP
      np->tm[0] = tmc = np->best[lsym];
      np->tm[1] = np->worst[lsym];
      const unsigned int up = (np->up >> lsym) & 1;
      np->encstate = es += up;
      np->par = par ^= -(uint64_t)up & PARSTEP;
      np->i = 0;	                        // Start with best branch
      continue;
    }
    // Threshold violated, can't go forward
    for(;;) {                                   // Look backward
      if(np == nodes || np[-1].gamma < t) {
	/* Can't back up either.
	 * Relax threshold and and look
	 * forward again to better branch.
	 */
	t -= delta;
	if(np->i != 0) {
	  np->i = 0;
	  np->encstate ^= 1;
	  np->par ^= PARSTEP;
	}
	tmc = np->tm[0];
	break;
      }
      // Back up
      if(--np < tail && np->i != 1) {
	np->i++;                          // Search next best branch (i is 0 or 1)
	np->encstate ^= 1;
	np->par ^= PARSTEP;
	tmc = np->tm[1];
	break;
      }                                   // else keep looking back
    }
    gam = np->gamma;
    es = np->encstate;
    par = np->par;
  }
  *maxnp = (int)(deep-nodes);
  *metric =  np->gamma;	                  // Return the final path metric  

  // Copy decoded data to user's buffer
  nbits >>= 3;
  np = &nodes[7];
  while(nbits-- != 0) {
    *data++ = np->encstate;
    np += 8;
  }
  *cycles = i+1;

  free(nodes);
  if(i >= maxcycles) return -1;	          // Decoder timed out
  return 0;		                  // Successful completion
}

// ---------------------------------------------------------------- jelinek.c

namespace {

// jelinek()'s buckets, malloc'd (and never freed) there. Its metric can
// take a bucket index past them (above 999 when the metric rises past
// about 25600, or wrapped below 0 past -6400: -J with a large -z), and
// wsprd then reads and writes whatever its heap holds there: other
// decodes, or a crash, as the heap lies. Here the array grows, zeroed, to
// take such an index; one wrapped or too far out to be in any heap
// (2^20 on) is where wsprd crashes for certain, and throws Crash.
class Buckets {
public:
    explicit Buckets(unsigned int n) : v_(n, 0) {}
    unsigned int& operator[](unsigned int i) {
        if (i >= (1u << 20)) throw Crash();
        if (i >= v_.size()) v_.resize((size_t)i + 1, 0);
        return v_[i];
    }

private:
    std::vector<unsigned int> v_;
};

}  // namespace

//Decoder - returns 0 on success, -1 on timeout
int jelinek(
            unsigned int *metric,	/* Final path metric (returned value) */
            unsigned int *cycles,	/* Cycle count (returned value) */
            unsigned char *data,	/* Decoded output data */
            unsigned char *symbols,	/* Raw deinterleaved input symbols */
            unsigned int nbits,	/* Number of output bits */
            unsigned int stacksize,
            struct snode *stack,
            int mettab[2][256],	/* Metric table, [sent sym][rx symbol] */
            unsigned int maxcycles)/* Decoding timeout in cycles per bit */
{
    
    // Compute branch metrics for each symbol pair
    // The sequential decoding algorithm only uses the metrics, not the
    // symbol values.
    unsigned int i;
    long int metrics[81][4];
    for(i=0; i<nbits; i++){
        metrics[i][0] = mettab[0][symbols[0]] + mettab[0][symbols[1]];
        metrics[i][1] = mettab[0][symbols[0]] + mettab[1][symbols[1]];
        metrics[i][2] = mettab[1][symbols[0]] + mettab[0][symbols[1]];
        metrics[i][3] = mettab[1][symbols[0]] + mettab[1][symbols[1]];
        symbols += 2;
    }
    
    // zero the stack
    memset(stack,0,stacksize*sizeof(struct snode));
    
    // initialize the loop variables
    unsigned int lsym, ntail=31;
    uint64_t encstate=0;
    unsigned int nbuckets=1000;
    unsigned int low_bucket=nbuckets-1; //will be set on first run-through
    unsigned int high_bucket=0;
    Buckets buckets(nbuckets);
    unsigned int bucket;
    unsigned int ptr=1;
    unsigned int stackptr=1; //pointer values of 0 are reserved (they mean that a bucket is empty)
    unsigned int depth=0, nbits_minus_ntail=nbits-ntail;
    unsigned int stacksize_minus_1=stacksize-1;
    long int totmet0, totmet1, gamma=0;
    
    unsigned int ncycles=maxcycles*nbits;
    /********************* Start the stack decoder *****************/
    for (i=1; i <= ncycles; i++) {
#ifdef DEBUG
        printf("***stackptr=%ld, depth=%d, gamma=%d, encstate=%lx, bucket %d, low_bucket %d, high_bucket %d\n",
               stackptr, depth, gamma, encstate, bucket, low_bucket, high_bucket);
#endif
        // no need to store more than 7 bytes (56 bits) for encoder state because
        // only 50 bits are not 0's.
        if( depth < 56 ) {
            encstate=encstate<<1;
            ENCODE(lsym,encstate); // get channel symbols associated with the 0 branch
        } else {
            ENCODE(lsym,encstate<<(depth-55));
        }

        // lsym are the 0-branch channel symbols and 3^lsym are the 1-branch
        // channel symbols (due to a special property of our generator polynomials)
        totmet0 = gamma+metrics[depth][lsym];   // total metric for 0-branch daughter node
        totmet1 = gamma+metrics[depth][3^lsym]; // total metric for 1-branch daughter node
        depth++; //the depth of the daughter nodes

        bucket=(totmet0>>5)+200; //fast, but not particularly safe - totmet can be negative
        if( bucket > high_bucket ) high_bucket=bucket;
        if( bucket < low_bucket ) low_bucket=bucket;
       
        // place the 0 node on the stack, overwriting the parent (current) node
        stack[ptr].encstate=encstate;
        stack[ptr].gamma=totmet0;
        stack[ptr].depth=depth;
        stack[ptr].jpointer=buckets[bucket];
        buckets[bucket]=ptr;
        
        // if in the tail, only need to evaluate the "0" branch.
        // Otherwise, enter this "if" and place the 1 node on the stack,
        if( depth <= nbits_minus_ntail ) {
            if( stackptr < stacksize_minus_1 ) {
                stackptr++;
                ptr=stackptr;
            } else { // stack full
                while( buckets[low_bucket] == 0 ) { //write latest to where the top of the lowest bucket points
                    low_bucket++;
                }
                ptr=buckets[low_bucket];
                buckets[low_bucket]=stack[ptr].jpointer; //make bucket point to next older entry
            }

            bucket=(totmet1>>5)+200; //this may not be safe on all compilers
            if( bucket > high_bucket ) high_bucket=bucket;
            if( bucket < low_bucket ) low_bucket=bucket;
            
            stack[ptr].encstate=encstate+1;
            stack[ptr].gamma=totmet1;
            stack[ptr].depth=depth;
            stack[ptr].jpointer=buckets[bucket];
            buckets[bucket]=ptr;
        }

    // pick off the latest entry from the high bucket
        while( buckets[high_bucket] == 0 ) {
            high_bucket--;
        }
        ptr=buckets[high_bucket];
        buckets[high_bucket]=stack[ptr].jpointer;
        depth=stack[ptr].depth;
        gamma=stack[ptr].gamma;
        encstate=stack[ptr].encstate;

        // we are done if the top entry on the stack is at depth nbits
        if (depth == nbits) {
            break;
        }
    }
    
    *cycles = i+1;
    *metric =  gamma;	/* Return final path metric */

    //    printf("cycles %d stackptr=%d, depth=%d, gamma=%d, encstate=%lx\n",
    //           *cycles, stackptr, depth, *metric, encstate);
    
    for (i=0; i<7; i++) {
        data[i]=(encstate>>(48-i*8))&(0x00000000000000ff);
    }
    for (i=7; i<11; i++) {
        data[i]=0;
    }

    if(*cycles/nbits >= maxcycles) //timed out
    {
        return -1;
    }
    return 0;		//success
}

// ----------------------------------------------------------- metric_tables.c

/*******************************************************************************
* 4 metric tables calculated via simulation for 2-FSK with Es/No=0,3,6,9 dB
* tables were calculated for constant rms noise level of 50. The symbol vector
* should be normalized to have rms amplitude equal to "symbol_scale".
********************************************************************************/
//float symbol_scale[5]={42.6, 53.3, 72.7, 100.2, 125.4};
const float metric_tables[5][256]={
    {0.9782,      0.9695,      0.9689,      0.9669,      0.9666,      0.9653,      0.9638,      0.9618,      0.9599,      0.9601, 
     0.9592,      0.9570,      0.9556,      0.9540,      0.9525,      0.9527,      0.9486,      0.9477,      0.9450,      0.9436, 
     0.9424,      0.9400,      0.9381,      0.9360,      0.9340,      0.9316,      0.9301,      0.9272,      0.9254,      0.9224, 
     0.9196,      0.9171,      0.9154,      0.9123,      0.9076,      0.9061,      0.9030,      0.9000,      0.8965,      0.8934, 
     0.8903,      0.8874,      0.8834,      0.8792,      0.8760,      0.8726,      0.8685,      0.8639,      0.8599,      0.8550, 
     0.8504,      0.8459,      0.8422,      0.8364,      0.8320,      0.8262,      0.8215,      0.8159,      0.8111,      0.8052, 
     0.7996,      0.7932,      0.7878,      0.7812,      0.7745,      0.7685,      0.7616,      0.7550,      0.7479,      0.7405, 
     0.7336,      0.7255,      0.7184,      0.7102,      0.7016,      0.6946,      0.6860,      0.6769,      0.6687,      0.6598, 
     0.6503,      0.6416,      0.6325,      0.6219,      0.6122,      0.6016,      0.5920,      0.5818,      0.5711,      0.5606, 
     0.5487,      0.5374,      0.5266,      0.5142,      0.5020,      0.4908,      0.4784,      0.4663,      0.4532,      0.4405, 
     0.4271,      0.4144,      0.4006,      0.3865,      0.3731,      0.3594,      0.3455,      0.3304,      0.3158,      0.3009, 
     0.2858,      0.2708,      0.2560,      0.2399,      0.2233,      0.2074,      0.1919,      0.1756,      0.1590,      0.1427, 
     0.1251,      0.1074,      0.0905,      0.0722,      0.0550,      0.0381,      0.0183,      0.0000,     -0.0185,     -0.0391, 
    -0.0571,     -0.0760,     -0.0966,     -0.1160,     -0.1370,     -0.1584,     -0.1787,     -0.1999,     -0.2214,     -0.2423, 
    -0.2643,     -0.2879,     -0.3114,     -0.3336,     -0.3568,     -0.3806,     -0.4050,     -0.4293,     -0.4552,     -0.4798, 
    -0.5046,     -0.5296,     -0.5564,     -0.5836,     -0.6093,     -0.6372,     -0.6645,     -0.6933,     -0.7208,     -0.7495, 
    -0.7763,     -0.8065,     -0.8378,     -0.8660,     -0.8964,     -0.9293,     -0.9592,     -0.9907,     -1.0214,     -1.0509, 
    -1.0850,     -1.1168,     -1.1528,     -1.1847,     -1.2157,     -1.2511,     -1.2850,     -1.3174,     -1.3540,     -1.3900, 
    -1.4201,     -1.4580,     -1.4956,     -1.5292,     -1.5683,     -1.6030,     -1.6411,     -1.6789,     -1.7147,     -1.7539, 
    -1.7887,     -1.8289,     -1.8699,     -1.9043,     -1.9469,     -1.9849,     -2.0267,     -2.0610,     -2.1028,     -2.1391, 
    -2.1855,     -2.2215,     -2.2712,     -2.3033,     -2.3440,     -2.3870,     -2.4342,     -2.4738,     -2.5209,     -2.5646, 
    -2.6016,     -2.6385,     -2.6868,     -2.7356,     -2.7723,     -2.8111,     -2.8524,     -2.9009,     -2.9428,     -2.9879, 
    -3.0103,     -3.0832,     -3.1340,     -3.1628,     -3.2049,     -3.2557,     -3.3101,     -3.3453,     -3.4025,     -3.4317, 
    -3.4828,     -3.5270,     -3.5745,     -3.6181,     -3.6765,     -3.7044,     -3.7410,     -3.8118,     -3.8368,     -3.9549, 
    -3.9488,     -3.9941,     -4.0428,     -4.0892,     -4.1648,     -4.1965,     -4.1892,     -4.2565,     -4.3356,     -4.3948, 
    -4.4481,     -4.4607,     -4.5533,     -4.5809,     -4.5927,     -5.1047}, 
    {0.9978,      0.9962,      0.9961,      0.9959,      0.9958,      0.9954,      0.9949,      0.9950,      0.9947,      0.9942, 
     0.9940,      0.9939,      0.9933,      0.9931,      0.9928,      0.9924,      0.9921,      0.9916,      0.9911,      0.9909, 
     0.9903,      0.9900,      0.9892,      0.9887,      0.9883,      0.9877,      0.9869,      0.9863,      0.9857,      0.9848, 
     0.9842,      0.9835,      0.9825,      0.9817,      0.9808,      0.9799,      0.9791,      0.9777,      0.9767,      0.9757, 
     0.9744,      0.9729,      0.9716,      0.9704,      0.9690,      0.9674,      0.9656,      0.9641,      0.9625,      0.9609, 
     0.9587,      0.9567,      0.9548,      0.9524,      0.9501,      0.9478,      0.9453,      0.9426,      0.9398,      0.9371, 
     0.9339,      0.9311,      0.9277,      0.9242,      0.9206,      0.9168,      0.9131,      0.9087,      0.9043,      0.8999, 
     0.8953,      0.8907,      0.8857,      0.8803,      0.8747,      0.8690,      0.8632,      0.8572,      0.8507,      0.8439, 
     0.8368,      0.8295,      0.8217,      0.8138,      0.8058,      0.7972,      0.7883,      0.7784,      0.7694,      0.7597, 
     0.7489,      0.7378,      0.7269,      0.7152,      0.7030,      0.6911,      0.6782,      0.6643,      0.6506,      0.6371, 
     0.6211,      0.6054,      0.5897,      0.5740,      0.5565,      0.5393,      0.5214,      0.5027,      0.4838,      0.4643, 
     0.4436,      0.4225,      0.4004,      0.3787,      0.3562,      0.3324,      0.3089,      0.2839,      0.2584,      0.2321, 
     0.2047,      0.1784,      0.1499,      0.1213,      0.0915,      0.0628,      0.0314,      0.0000,     -0.0321,     -0.0657, 
    -0.0977,     -0.1324,     -0.1673,     -0.2036,     -0.2387,     -0.2768,     -0.3150,     -0.3538,     -0.3936,     -0.4327, 
    -0.4739,     -0.5148,     -0.5561,     -0.6000,     -0.6438,     -0.6889,     -0.7331,     -0.7781,     -0.8247,     -0.8712, 
    -0.9177,     -0.9677,     -1.0142,     -1.0631,     -1.1143,     -1.1686,     -1.2169,     -1.2680,     -1.3223,     -1.3752, 
    -1.4261,     -1.4806,     -1.5356,     -1.5890,     -1.6462,     -1.7041,     -1.7591,     -1.8124,     -1.8735,     -1.9311, 
    -1.9891,     -2.0459,     -2.1048,     -2.1653,     -2.2248,     -2.2855,     -2.3466,     -2.4079,     -2.4668,     -2.5263, 
    -2.5876,     -2.6507,     -2.7142,     -2.7761,     -2.8366,     -2.8995,     -2.9620,     -3.0279,     -3.0973,     -3.1576, 
    -3.2238,     -3.2890,     -3.3554,     -3.4215,     -3.4805,     -3.5518,     -3.6133,     -3.6812,     -3.7473,     -3.8140, 
    -3.8781,     -3.9450,     -4.0184,     -4.0794,     -4.1478,     -4.2241,     -4.2853,     -4.3473,     -4.4062,     -4.4839, 
    -4.5539,     -4.6202,     -4.6794,     -4.7478,     -4.8309,     -4.9048,     -4.9669,     -5.0294,     -5.1194,     -5.1732, 
    -5.2378,     -5.3094,     -5.3742,     -5.4573,     -5.5190,     -5.5728,     -5.6637,     -5.7259,     -5.7843,     -5.8854, 
    -5.9553,     -6.0054,     -6.0656,     -6.1707,     -6.2241,     -6.3139,     -6.3393,     -6.4356,     -6.5153,     -6.5758, 
    -6.6506,     -6.7193,     -6.7542,     -6.8942,     -6.9219,     -6.9605,     -7.1013,     -7.1895,     -7.1549,     -7.2799, 
     -7.4119,     -7.4608,     -7.5256,     -7.5879,     -7.7598,     -8.4120}, 
    {0.9999,      0.9998,      0.9998,      0.9998,      0.9998,      0.9998,      0.9997,      0.9997,      0.9997,      0.9997, 
     0.9997,      0.9996,      0.9996,      0.9996,      0.9995,      0.9995,      0.9994,      0.9994,      0.9994,      0.9993, 
     0.9993,      0.9992,      0.9991,      0.9991,      0.9990,      0.9989,      0.9988,      0.9988,      0.9988,      0.9986, 
     0.9985,      0.9984,      0.9983,      0.9982,      0.9980,      0.9979,      0.9977,      0.9976,      0.9974,      0.9971, 
     0.9969,      0.9968,      0.9965,      0.9962,      0.9960,      0.9957,      0.9953,      0.9950,      0.9947,      0.9941, 
     0.9937,      0.9933,      0.9928,      0.9922,      0.9917,      0.9911,      0.9904,      0.9897,      0.9890,      0.9882, 
     0.9874,      0.9863,      0.9855,      0.9843,      0.9832,      0.9819,      0.9806,      0.9792,      0.9777,      0.9760, 
     0.9743,      0.9724,      0.9704,      0.9683,      0.9659,      0.9634,      0.9609,      0.9581,      0.9550,      0.9516, 
     0.9481,      0.9446,      0.9406,      0.9363,      0.9317,      0.9270,      0.9218,      0.9160,      0.9103,      0.9038, 
     0.8972,      0.8898,      0.8822,      0.8739,      0.8647,      0.8554,      0.8457,      0.8357,      0.8231,      0.8115, 
     0.7984,      0.7854,      0.7704,      0.7556,      0.7391,      0.7210,      0.7038,      0.6840,      0.6633,      0.6408, 
     0.6174,      0.5939,      0.5678,      0.5410,      0.5137,      0.4836,      0.4524,      0.4193,      0.3850,      0.3482, 
     0.3132,      0.2733,      0.2315,      0.1891,      0.1435,      0.0980,      0.0493,      0.0000,     -0.0510,     -0.1052, 
    -0.1593,     -0.2177,     -0.2759,     -0.3374,     -0.4005,     -0.4599,     -0.5266,     -0.5935,     -0.6626,     -0.7328, 
    -0.8051,     -0.8757,     -0.9498,     -1.0271,     -1.1019,     -1.1816,     -1.2642,     -1.3459,     -1.4295,     -1.5077, 
    -1.5958,     -1.6818,     -1.7647,     -1.8548,     -1.9387,     -2.0295,     -2.1152,     -2.2154,     -2.3011,     -2.3904, 
    -2.4820,     -2.5786,     -2.6730,     -2.7652,     -2.8616,     -2.9546,     -3.0526,     -3.1445,     -3.2445,     -3.3416, 
    -3.4357,     -3.5325,     -3.6324,     -3.7313,     -3.8225,     -3.9209,     -4.0248,     -4.1278,     -4.2261,     -4.3193, 
    -4.4220,     -4.5262,     -4.6214,     -4.7242,     -4.8234,     -4.9245,     -5.0298,     -5.1250,     -5.2232,     -5.3267, 
    -5.4332,     -5.5342,     -5.6431,     -5.7270,     -5.8401,     -5.9350,     -6.0407,     -6.1418,     -6.2363,     -6.3384, 
    -6.4536,     -6.5429,     -6.6582,     -6.7433,     -6.8438,     -6.9478,     -7.0789,     -7.1894,     -7.2714,     -7.3815, 
    -7.4810,     -7.5575,     -7.6852,     -7.8071,     -7.8580,     -7.9724,     -8.1000,     -8.2207,     -8.2867,     -8.4017, 
    -8.5287,     -8.6347,     -8.7082,     -8.8319,     -8.9448,     -9.0355,     -9.1885,     -9.2095,     -9.2863,     -9.4186, 
    -9.5064,     -9.6386,     -9.7207,     -9.8286,     -9.9453,    -10.0701,    -10.1735,    -10.3001,    -10.2858,    -10.5427, 
   -10.5982,    -10.7361,    -10.7042,    -10.9212,    -11.0097,    -11.0469,    -11.1155,    -11.2812,    -11.3472,    -11.4988, 
   -11.5327,    -11.6692,    -11.9376,    -11.8606,    -12.1372,    -13.2539}, 
    {1.0000,      1.0000,      1.0000,      1.0000,      1.0000,      1.0000,      1.0000,      1.0000,      1.0000,      1.0000, 
     1.0000,      1.0000,      1.0000,      1.0000,      1.0000,      1.0000,      1.0000,      1.0000,      1.0000,      1.0000, 
     1.0000,      1.0000,      1.0000,      1.0000,      1.0000,      1.0000,      1.0000,      1.0000,      1.0000,      1.0000, 
     0.9999,      0.9999,      0.9999,      0.9999,      0.9999,      0.9999,      0.9999,      0.9999,      0.9999,      0.9999, 
     0.9999,      0.9998,      0.9998,      0.9998,      0.9998,      0.9997,      0.9997,      0.9997,      0.9997,      0.9996, 
     0.9996,      0.9995,      0.9995,      0.9994,      0.9994,      0.9993,      0.9992,      0.9991,      0.9991,      0.9989, 
     0.9988,      0.9986,      0.9985,      0.9983,      0.9981,      0.9980,      0.9977,      0.9974,      0.9971,      0.9968, 
     0.9965,      0.9962,      0.9956,      0.9950,      0.9948,      0.9941,      0.9933,      0.9926,      0.9919,      0.9910, 
     0.9899,      0.9889,      0.9877,      0.9863,      0.9845,      0.9829,      0.9811,      0.9791,      0.9769,      0.9741, 
     0.9716,      0.9684,      0.9645,      0.9611,      0.9563,      0.9519,      0.9463,      0.9406,      0.9344,      0.9272, 
     0.9197,      0.9107,      0.9016,      0.8903,      0.8791,      0.8653,      0.8523,      0.8357,      0.8179,      0.7988, 
     0.7779,      0.7562,      0.7318,      0.7024,      0.6753,      0.6435,      0.6089,      0.5700,      0.5296,      0.4860, 
     0.4366,      0.3855,      0.3301,      0.2735,      0.2114,      0.1443,      0.0682,      0.0000,     -0.0715,     -0.1604, 
    -0.2478,     -0.3377,     -0.4287,     -0.5277,     -0.6291,     -0.7384,     -0.8457,     -0.9559,     -1.0742,     -1.1913, 
    -1.3110,     -1.4238,     -1.5594,     -1.6854,     -1.8093,     -1.9414,     -2.0763,     -2.2160,     -2.3611,     -2.4876, 
    -2.6374,     -2.7710,     -2.9225,     -3.0591,     -3.2077,     -3.3452,     -3.4916,     -3.6316,     -3.7735,     -3.9296, 
    -4.0682,     -4.2334,     -4.3607,     -4.5270,     -4.6807,     -4.8108,     -4.9753,     -5.1212,     -5.2631,     -5.4042, 
    -5.5510,     -5.7227,     -5.8794,     -6.0244,     -6.1677,     -6.3271,     -6.4862,     -6.6130,     -6.7449,     -6.9250, 
    -7.1232,     -7.1736,     -7.3628,     -7.5596,     -7.6906,     -7.8129,     -7.9817,     -8.1440,     -8.3016,     -8.4797, 
    -8.5734,     -8.7692,     -8.9198,     -9.0610,     -9.1746,     -9.3536,     -9.5939,     -9.6957,     -9.8475,     -9.9639, 
   -10.1730,    -10.2427,    -10.4573,    -10.5413,    -10.7303,    -10.9339,    -11.0215,    -11.2047,    -11.2894,    -11.4572, 
   -11.6256,    -11.7794,    -11.8801,    -12.1717,    -12.2354,    -12.3686,    -12.6195,    -12.6527,    -12.8247,    -12.9560, 
   -13.3265,    -13.1667,    -13.4274,    -13.6064,    -13.5515,    -13.9501,    -13.9926,    -14.4049,    -14.1653,    -14.4348, 
   -14.7983,    -14.7807,    -15.2349,    -15.3536,    -15.3026,    -15.2739,    -15.7170,    -16.2161,    -15.9185,    -15.9490, 
   -16.6258,    -16.5568,    -16.4318,    -16.7999,    -16.4101,    -17.6393,    -17.7643,    -17.2644,    -17.5973,    -17.0403, 
   -17.7039,    -18.0073,    -18.1840,    -18.3848,    -18.6286,    -20.7063}, 
    {1.0000,      1.0000,      1.0000,      1.0000,      1.0000,      1.0000,      1.0000,      1.0000,      1.0000,      1.0000, 
     1.0000,      1.0000,      1.0000,      1.0000,      1.0000,      1.0000,      1.0000,      1.0000,      1.0000,      1.0000, 
     1.0000,      1.0000,      1.0000,      1.0000,      1.0000,      1.0000,      1.0000,      1.0000,      1.0000,      1.0000, 
     0.9999,      0.9999,      0.9999,      0.9999,      0.9999,      0.9999,      0.9999,      0.9999,      0.9999,      0.9999, 
     0.9999,      0.9998,      0.9998,      0.9998,      0.9998,      0.9997,      0.9997,      0.9997,      0.9997,      0.9996, 
     0.9996,      0.9995,      0.9995,      0.9994,      0.9994,      0.9993,      0.9992,      0.9991,      0.9991,      0.9989, 
     0.9988,      0.9986,      0.9985,      0.9983,      0.9981,      0.9980,      0.9977,      0.9974,      0.9971,      0.9968, 
     0.9965,      0.9962,      0.9956,      0.9950,      0.9948,      0.9941,      0.9933,      0.9926,      0.9919,      0.9910, 
     0.9899,      0.9889,      0.9877,      0.9863,      0.9845,      0.9829,      0.9811,      0.9791,      0.9769,      0.9741, 
     0.9716,      0.9684,      0.9645,      0.9611,      0.9563,      0.9519,      0.9463,      0.9406,      0.9344,      0.9272, 
     0.9197,      0.9107,      0.9016,      0.8903,      0.8791,      0.8653,      0.8523,      0.8357,      0.8179,      0.7988, 
     0.7779,      0.7562,      0.7318,      0.7024,      0.6753,      0.6435,      0.6089,      0.5700,      0.5296,      0.4860, 
     0.4366,      0.3855,      0.3301,      0.2735,      0.2114,      0.1443,      0.0682,      0.0000,     -0.0715,     -0.1604, 
    -0.2478,     -0.3377,     -0.4287,     -0.5277,     -0.6291,     -0.7384,     -0.8457,     -0.9559,     -1.0742,     -1.1913, 
    -1.3110,     -1.4238,     -1.5594,     -1.6854,     -1.8093,     -1.9414,     -2.0763,     -2.2160,     -2.3611,     -2.4876, 
    -2.6374,     -2.7710,     -2.9225,     -3.0591,     -3.2077,     -3.3452,     -3.4916,     -3.6316,     -3.7735,     -3.9296, 
    -4.0682,     -4.2334,     -4.3607,     -4.5270,     -4.6807,     -4.8108,     -4.9753,     -5.1212,     -5.2631,     -5.4042, 
    -5.5510,     -5.7227,     -5.8794,     -6.0244,     -6.1677,     -6.3271,     -6.4862,     -6.6130,     -6.7449,     -6.9250, 
    -7.1232,     -7.1736,     -7.3628,     -7.5596,     -7.6906,     -7.8129,     -7.9817,     -8.1440,     -8.3016,     -8.4797, 
    -8.5734,     -8.7692,     -8.9198,     -9.0610,     -9.1746,     -9.3536,     -9.5939,     -9.6957,     -9.8475,     -9.9639, 
   -10.1730,    -10.2427,    -10.4573,    -11.7794,    -11.8801,    -12.1717,    -12.2354,    -12.3686,    -12.6195,    -12.6527,
   -12.8247,    -12.9560,    -13.3265,    -13.1667,    -13.4274,    -13.6064,    -13.5515,    -13.9501,    -13.9926,    -14.4049,
   -14.1653,    -14.4348,    -14.7983,    -14.7807,    -15.2349,    -15.3536,    -15.3026,    -15.2739,    -15.7170,    -16.2161,
   -15.9185,    -15.9490,    -16.6258,    -16.5568,    -16.4318,    -16.7999,    -16.4101,    -17.6393,    -17.7643,    -17.2644,
   -17.5973,    -17.0403,    -17.7039,    -18.0073,    -18.1840,    -18.3848,    -18.6286,    -20.7063,          1.43370769e-019,
     2.64031087e-006,  6.6908396e+031,          1.77537994e+028,          2.79322819e+020,             1.94326e-019,
     0.00019371575,   2.80722121e-041}};

}  // namespace wsprd::detail
