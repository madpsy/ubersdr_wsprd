// wsprd's message unpacking and packing, call hash and hash table lines:
// WSJT-X's lib/wsprd/wsprd_utils.c, wsprsim_utils.c, nhash.c and the
// lib/lookup3.c hash under it, at 567ad29ce, for wsprd.cpp (see wsprd.h).
//
// wsprd_utils.c, wsprsim_utils.c: Copyright 2001-2015, Joe Taylor, K1JT, and
// 2014-2015, Steven Franke, K9AN; GNU GPL v3 (or later), and so is this.
// lookup3.c: Bob Jenkins, 2006, public domain; WSJT-X's wrapper (C) 2008-2014
// Joseph Taylor, K1JT, GNU GPL v3+.
//
// The reference's code, kept as it was so the two read side by side. What
// had to change, and nothing else:
//   - C++: casts on malloc/calloc, namespace wsprd::detail;
//   - strtok() (process-wide state) is glibc's, its position per thread;
//   - wsprsim's printdata branch is gone (printdata is always 0 in wsprd);
//   - wsjt_hashlittle reads the key a byte at a time: its aligned 32- and
//     16-bit paths give the same hash on a little-endian machine, but read
//     char data through wider pointers.
#include "wsprd_internal.h"

#include <cctype>
#include <cerrno>
#include <cstdarg>
#include <cstdio>
#include <cstdlib>
#include <cstring>

// strtok keeps a pointer into its caller's buffer, as glibc's does; it is
// never followed once the caller returns (each use starts with a buffer).
#if defined(__GNUC__) && !defined(__clang__) && __GNUC__ >= 12
#pragma GCC diagnostic ignored "-Wdangling-pointer"
#endif

namespace wsprd::detail {

namespace {

// glibc's strtok, its saved position per thread.
thread_local char* strtok_save = nullptr;

char* strtok(char* s, const char* delim) {
    if (s == nullptr) s = strtok_save;
    if (*s == '\0') {
        strtok_save = s;
        return nullptr;
    }
    s += std::strspn(s, delim);
    if (*s == '\0') {
        strtok_save = s;
        return nullptr;
    }
    char* end = s + std::strcspn(s, delim);
    if (*end == '\0') {
        strtok_save = end;
        return s;
    }
    *end = '\0';
    strtok_save = end + 1;
    return s;
}

// ---------------------------------------------------------------- lookup3.c

#define rot(x,k) (((x)<<(k)) | ((x)>>(32-(k))))

#define mix(a,b,c) \
{ \
  a -= c;  a ^= rot(c, 4);  c += b; \
  b -= a;  b ^= rot(a, 6);  a += c; \
  c -= b;  c ^= rot(b, 8);  b += a; \
  a -= c;  a ^= rot(c,16);  c += b; \
  b -= a;  b ^= rot(a,19);  a += c; \
  c -= b;  c ^= rot(b, 4);  b += a; \
}

#define final(a,b,c) \
{ \
  c ^= b; c -= rot(b,14); \
  a ^= c; a -= rot(c,11); \
  b ^= a; b -= rot(a,25); \
  c ^= b; c -= rot(b,16); \
  a ^= c; a -= rot(c,4);  \
  b ^= a; b -= rot(a,14); \
  c ^= b; c -= rot(b,24); \
}

uint32_t wsjt_hashlittle(const void *key, size_t length, uint32_t initval)
{
  uint32_t a,b,c;                                          /* internal state */

  /* Set up the internal state */
  a = b = c = 0xdeadbeef + ((uint32_t)length) + initval;

  {                               /* need to read the key one byte at a time */
    const uint8_t *k = (const uint8_t *)key;

    /*--------------- all but the last block: affect some 32 bits of (a,b,c) */
    while (length > 12)
    {
      a += k[0];
      a += ((uint32_t)k[1])<<8;
      a += ((uint32_t)k[2])<<16;
      a += ((uint32_t)k[3])<<24;
      b += k[4];
      b += ((uint32_t)k[5])<<8;
      b += ((uint32_t)k[6])<<16;
      b += ((uint32_t)k[7])<<24;
      c += k[8];
      c += ((uint32_t)k[9])<<8;
      c += ((uint32_t)k[10])<<16;
      c += ((uint32_t)k[11])<<24;
      mix(a,b,c);
      length -= 12;
      k += 12;
    }

    /*-------------------------------- last block: affect all 32 bits of (c) */
    switch(length)                   /* all the case statements fall through */
    {
    case 12: c+=((uint32_t)k[11])<<24;
      /* fall through */
    case 11: c+=((uint32_t)k[10])<<16;
      /* fall through */
    case 10: c+=((uint32_t)k[9])<<8;
      /* fall through */
    case 9 : c+=k[8];
      /* fall through */
    case 8 : b+=((uint32_t)k[7])<<24;
      /* fall through */
    case 7 : b+=((uint32_t)k[6])<<16;
      /* fall through */
    case 6 : b+=((uint32_t)k[5])<<8;
      /* fall through */
    case 5 : b+=k[4];
      /* fall through */
    case 4 : a+=((uint32_t)k[3])<<24;
      /* fall through */
    case 3 : a+=((uint32_t)k[2])<<16;
      /* fall through */
    case 2 : a+=((uint32_t)k[1])<<8;
      /* fall through */
    case 1 : a+=k[0];
             break;
    case 0 : return c;
    }
  }

  final(a,b,c);
  return c;
}

#undef rot
#undef mix
#undef final

// ------------------------------------------------------------ wsprd_utils.c

int format_checked(char *dest, size_t dest_size, char const *format, ...)
{
    va_list args;
    va_start(args, format);
    int nwritten = vsnprintf(dest, dest_size, format, args);
    va_end(args);
    return nwritten >= 0 && (size_t)nwritten < dest_size;
}

int read_bounded_token(char const **cursor, char *dest,
                              size_t dest_size, int required)
{
    char const *p = *cursor;
    while( isspace((unsigned char)*p) ) p++;

    if( *p == '\0' ) {
        if( required ) return 0;
        dest[0] = '\0';
        *cursor = p;
        return 1;
    }

    char const *start = p;
    while( *p != '\0' && !isspace((unsigned char)*p) ) p++;

    size_t len = (size_t)(p - start);
    if( len >= dest_size ) return 0;

    memcpy(dest, start, len);
    dest[len] = '\0';
    *cursor = p;
    return 1;
}

int store_hash_call(int ihash, char const *callsign, char *hashtab)
{
    if( ihash < 0 || ihash >= WSPRD_HASH_COUNT ) return 0;
    if( strlen(callsign) >= WSPRD_CALLSIGN_SIZE ) return 0;
    if( !format_checked(hashtab + ihash * WSPRD_CALLSIGN_SIZE,
                        WSPRD_CALLSIGN_SIZE, "%s", callsign) ) return 0;
    return 1;
}

int store_hash_grid(int ihash, char const *grid, char *loctab)
{
    if( ihash < 0 || ihash >= WSPRD_HASH_COUNT ) return 0;
    if( strlen(grid) >= WSPRD_GRID4_SIZE ) return 0;
    if( !format_checked(loctab + ihash * WSPRD_GRID4_SIZE,
                        WSPRD_GRID4_SIZE, "%s", grid) ) return 0;
    return 1;
}

int store_hash_entry(int ihash, char const *callsign, char const *grid,
                            char *hashtab, char *loctab)
{
    if( !store_hash_call(ihash, callsign, hashtab) ) return 0;
    if( !store_hash_grid(ihash, grid, loctab) ) return 0;
    return 1;
}

}  // namespace

// nhash.c
uint32_t nhash(const void *key, size_t length, uint32_t initval)
{
  /* WSPR table indices use 15 bits; empty keys historically bypass the mask. */
  uint32_t hash = wsjt_hashlittle(key, length, initval);
  return length == 0 ? hash : hash & (uint32_t)0x7fff;
}

int wsprd_load_hash_line(char const *line, char *hashtab, char *loctab)
{
    char const *p = line;
    errno = 0;
    char *end = NULL;
    long ihash = strtol(p, &end, 10);
    if( p == end || errno == ERANGE ||
        ihash < 0 || ihash >= WSPRD_HASH_COUNT ) return 0;
    if( !isspace((unsigned char)*end) ) return 0;

    char callsign[WSPRD_CALLSIGN_SIZE];
    char grid[WSPRD_GRID4_SIZE];
    p = end;
    if( !read_bounded_token(&p, callsign, sizeof callsign, 1) ) return 0;
    if( !read_bounded_token(&p, grid, sizeof grid, 0) ) return 0;

    if( !store_hash_call((int)ihash, callsign, hashtab) ) return 0;
    /* Legacy loader semantics: missing grid preserves the existing locator. */
    if( grid[0] != '\0' &&
        !store_hash_grid((int)ihash, grid, loctab) ) return 0;
    return 1;
}

void unpack50( signed char *dat, int32_t *n1, int32_t *n2 )
{
    int32_t i,i4;

    i=dat[0];
    i4=i&255;
    *n1=i4<<20;

    i=dat[1];
    i4=i&255;
    *n1=*n1+(i4<<12);

    i=dat[2];
    i4=i&255;
    *n1=*n1+(i4<<4);

    i=dat[3];
    i4=i&255;
    *n1=*n1+((i4>>4)&15);
    *n2=(i4&15)<<18;

    i=dat[4];
    i4=i&255;
    *n2=*n2+(i4<<10);

    i=dat[5];
    i4=i&255;
    *n2=*n2+(i4<<2);

    i=dat[6];
    i4=i&255;
    *n2=*n2+((i4>>6)&3);
}

int unpackcall( int32_t ncall, char *call )
{
    char c[]={'0','1','2','3','4','5','6','7','8','9','A','B','C','D','E',
        'F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T',
        'U','V','W','X','Y','Z',' '};
    int32_t n;
    int i;
    char tmp[7];

    n=ncall;
    strcpy(call,"......");
    if (n < 262177560 ) {
        i=n%27+10;
        tmp[5]=c[i];
        n=n/27;
        i=n%27+10;
        tmp[4]=c[i];
        n=n/27;
        i=n%27+10;
        tmp[3]=c[i];
        n=n/27;
        i=n%10;
        tmp[2]=c[i];
        n=n/10;
        i=n%36;
        tmp[1]=c[i];
        n=n/36;
        i=n;
        tmp[0]=c[i];
        tmp[6]='\0';
        // remove leading whitespace
        for(i=0; i<5; i++) {
            if( tmp[i] != c[36] )
                break;
        }
        sprintf(call,"%-6s",&tmp[i]);
        // remove trailing whitespace
        for(i=0; i<6; i++) {
            if( call[i] == c[36] ) {
                call[i]='\0';
            }
        }
    } else {
        return 0;
    }
    return 1;
}

int unpackgrid( int32_t ngrid, char *grid)
{
    char c[]={'0','1','2','3','4','5','6','7','8','9','A','B','C','D','E',
        'F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T',
        'U','V','W','X','Y','Z',' '};
    int dlat, dlong;

    ngrid=ngrid>>7;
    if( ngrid < 32400 ) {
        dlat=(ngrid%180)-90;
        dlong=(ngrid/180)*2 - 180 + 2;
        if( dlong < -180 )
            dlong=dlong+360;
        if( dlong > 180 )
            dlong=dlong+360;
        int nlong = 60.0*(180.0-dlong)/5.0;
        int n1 = nlong/240;
        int n2 = (nlong - 240*n1)/24;
        grid[0] = c[10+n1];
        grid[2]=  c[n2];

        int nlat = 60.0*(dlat+90)/2.5;
        n1 = nlat/240;
        n2 = (nlat-240*n1)/24;
        grid[1]=c[10+n1];
        grid[3]=c[n2];
    } else {
        strcpy(grid,"XXXX");
        return 0;
    }
    return 1;
}

int unpackpfx( int32_t nprefix, char *call)
{
    char nc, pfx[4]={'\0'}, tmpcall[7];
    int i;
    int32_t n;

    if( !format_checked(tmpcall, sizeof tmpcall, "%s", call) ) return 0;
    if( nprefix < 60000 ) {
        // add a prefix of 1 to 3 characters
        n=nprefix;
        for (i=2; i>=0; i--) {
            nc=n%37;
            if( (nc >= 0) & (nc <= 9) ) {
                pfx[i]=nc+48;
            }
            else if( (nc >= 10) & (nc <= 35) ) {
                pfx[i]=nc+55;
            }
            else {
                pfx[i]=' ';
            }
            n=n/37;
        }

        char * p = strrchr(pfx,' ');
        if( !format_checked(call, WSPRD_CALLSIGN_SIZE, "%s/%s",
                            p ? p + 1 : pfx, tmpcall) ) return 0;

    } else {
        // add a suffix of 1 or 2 characters
        nc=nprefix-60000;
        if( (nc >= 0) & (nc <= 9) ) {
            pfx[0]=nc+48;
            if( !format_checked(call, WSPRD_CALLSIGN_SIZE, "%s/%.1s",
                                tmpcall, pfx) ) return 0;
        }
        else if( (nc >= 10) & (nc <= 35) ) {
            pfx[0]=nc+55;
            if( !format_checked(call, WSPRD_CALLSIGN_SIZE, "%s/%.1s",
                                tmpcall, pfx) ) return 0;
        }
        else if( (nc >= 36) & (nc <= 125) ) {
            pfx[0]=(nc-26)/10+48;
            pfx[1]=(nc-26)%10+48;
            if( !format_checked(call, WSPRD_CALLSIGN_SIZE, "%s/%.2s",
                                tmpcall, pfx) ) return 0;
        }
        else {
            return 0;
        }
    }
    return 1;
}

void deinterleave(unsigned char *sym)
{
    unsigned char tmp[162];
    unsigned char p, i, j;

    p=0;
    i=0;
    while (p<162) {
        j=((i * 0x80200802ULL) & 0x0884422110ULL) * 0x0101010101ULL >> 32;
        if (j < 162 ) {
            tmp[p]=sym[j];
            p=p+1;
        }
        i=i+1;
    }
    for (i=0; i<162; i++) {
        sym[i]=tmp[i];
    }
}

int floatcomp(const void* elem1, const void* elem2)
{
    if(*(const float*)elem1 < *(const float*)elem2)
        return -1;
    return *(const float*)elem1 > *(const float*)elem2;
}

int unpk_(signed char *message, char *hashtab, char *loctab, char *call_loc_pow, char *callsign)
{
    int n1,n2,n3,ndbm,ihash,nadd,noprint=0;
    char grid[WSPRD_GRID4_SIZE],grid6[WSPRD_GRID6_SIZE];

    unpack50(message,&n1,&n2);
    if( !unpackcall(n1,callsign) ) return 1;
    if( !unpackgrid(n2, grid) ) return 1;
    int ntype = (n2&127) - 64;
    callsign[WSPRD_CALLSIGN_SIZE - 1]=0;
    grid[4]=0;

    /*
     Based on the value of ntype, decide whether this is a Type 1, 2, or
     3 message.

     * Type 1: 6 digit call, grid, power - ntype is positive and is a member
     of the set {0,3,7,10,13,17,20...60}

     * Type 2: extended callsign, power - ntype is positive but not
     a member of the set of allowed powers

     * Type 3: hash, 6 digit grid, power - ntype is negative.
     */

    if( (ntype >= 0) && (ntype <= 62) ) {
        int nu=ntype%10;
        if( nu == 0 || nu == 3 || nu == 7 ) {
            ndbm=ntype;
            if( !format_checked(call_loc_pow, WSPRD_MESSAGE_SIZE, "%s %.4s %2d",
                                callsign, grid, ndbm) ) return 1;
            ihash=nhash(callsign,strlen(callsign),(uint32_t)146);
            if( !store_hash_entry(ihash, callsign, grid, hashtab, loctab) ) return 1;
        } else {
            nadd=nu;
            if( nu > 3 ) nadd=nu-3;
            if( nu > 7 ) nadd=nu-7;
            n3=n2/128+32768*(nadd-1);
            if( !unpackpfx(n3,callsign) ) return 1;
            ndbm=ntype-nadd;
            if( !format_checked(call_loc_pow, WSPRD_MESSAGE_SIZE, "%s %2d",
                                callsign, ndbm) ) return 1;
            int nu=ndbm%10;
            if( nu == 0 || nu == 3 || nu == 7 || nu == 10 ) { //make sure power is OK
                ihash=nhash(callsign,strlen(callsign),(uint32_t)146);
                if( !store_hash_call(ihash, callsign, hashtab) ) return 1;
            } else noprint=1;
        }
    } else if ( ntype < 0 ) {
        ndbm=-(ntype+1);
        memset(grid6,0,sizeof grid6);
//        size_t len=strlen(callsign);
        size_t len=6;
        if( !format_checked(grid6, sizeof grid6, "%.1s%.*s",
                            callsign+len-1, (int)(len-1), callsign) ) return 1;
        int nu=ndbm%10;
        if ((nu != 0 && nu != 3 && nu != 7 && nu != 10) ||
            !isalpha(grid6[0]) || !isalpha(grid6[1]) ||
            !isdigit(grid6[2]) || !isdigit(grid6[3])) {
               // not testing 4'th and 5'th chars because of this case: <PA0SKT/2> JO33 40
               // grid is only 4 chars even though this is a hashed callsign...
               //         isalpha(grid6[4]) && isalpha(grid6[5]) ) ) {
            noprint=1;
        }

        ihash=(n2-ntype-64)/128;
        if( strncmp(hashtab+ihash*13,"\0",1) != 0 ) {
            if( !format_checked(callsign, WSPRD_CALLSIGN_SIZE, "<%s>",
                                hashtab+ihash*13) ) return 1;
        } else {
            if( !format_checked(callsign, WSPRD_CALLSIGN_SIZE, "%5s",
                                "<...>") ) return 1;
        }

        if( !format_checked(call_loc_pow, WSPRD_MESSAGE_SIZE, "%s %s %2d",
                            callsign, grid6, ndbm) ) return 1;


        // I don't know what to do with these... They show up as "A000AA" grids.
        if( ntype == -64 ) noprint=1;
    }
    return noprint;
}

// ---------------------------------------------------------- wsprsim_utils.c

namespace {

char get_locator_character_code(char ch) {
    if( ch >=48 && ch <=57 ) { //0-9
        return ch-48;
    }
    if( ch == 32 ) {  //space
        return 36;
    }
    if( ch >= 65 && ch <= 82 ) { //A-Z
        return ch-65;
    }
    return -1;
}

char get_callsign_character_code(char ch) {
    if( ch >=48 && ch <=57 ) { //0-9
        return ch-48;
    }
    if( ch == 32 ) {  //space
        return 36;
    }
    if( ch >= 65 && ch <= 90 ) { //A-Z
        return ch-55;
    }
    return -1;
}

long unsigned int pack_grid4_power(char const *grid4, int power) {
    long unsigned int m;

    m=(179-10*grid4[0]-grid4[2])*180+10*grid4[1]+grid4[3];
    m=m*128+power+64;
    return m;
}

long unsigned int pack_call(char const *callsign) {
    unsigned int i;
    long unsigned int n;
    // 7, not the reference's 6: a 6-character call with a digit second and
    // not third ("K1ABCD") is shifted one along, its last character into
    // call6[6], past the reference's array. Only call6[0..5] are read.
    char call6[7];
    memset(call6,' ',sizeof(call6));
    // callsign is 6 characters in length. Exactly.
    size_t call_len = strlen(callsign);
    if( call_len > 6 ) {
        return 0;
    }
    if( isdigit(callsign[2]) ) {
        for (i=0; i<call_len; i++) {
            call6[i]=callsign[i];
        }
    } else if( isdigit(callsign[1]) ) {
        for (i=1; i<call_len+1; i++) {
            call6[i]=callsign[i-1];
        }
    }
    for (i=0; i<6; i++) {
        call6[i]=get_callsign_character_code(call6[i]);
    }
    n = call6[0];
    n = n*36+call6[1];
    n = n*10+call6[2];
    n = n*27+call6[3]-10;
    n = n*27+call6[4]-10;
    n = n*27+call6[5]-10;
    return n;
}

void pack_prefix(char *callsign, int32_t *n, int32_t *m, int32_t *nadd ) {
    size_t i;
    char * call6 = (char *)calloc(7,sizeof (char));
    size_t i1=strcspn(callsign,"/");

    if( callsign[i1+2] == 0 ) {
        //single char suffix
        for (i=0; i<i1; i++) {
            call6[i]=callsign[i];
        }
        call6[i] = '\0';
        *n=pack_call(call6);
        *nadd=1;
        int nc = callsign[i1+1];
        if( nc >= 48 && nc <= 57 ) {
            *m=nc-48;
        } else if ( nc >= 65 && nc <= 90 ) {
            *m=nc-65+10;
        } else {
            *m=38;
        }
        *m=60000-32768+*m;
    } else if( callsign[i1+3]==0 ) {
        //two char suffix
        for (i=0; i<i1; i++) {
            call6[i]=callsign[i];
        }
        *n=pack_call(call6);
        *nadd=1;
        *m=10*(callsign[i1+1]-48)+(callsign[i1+2]-48);
        *m=60000 + 26 + *m;
    } else {
        char const * pfx = strtok (callsign,"/");
        char const * call = strtok(NULL," ");
        *n = pack_call (call);
        size_t plen=strlen (pfx);
        if( plen ==1 ) {
            *m=36;
            *m=37*(*m)+36;
        } else if( plen == 2 ) {
            *m=36;
        } else {
            *m=0;
        }
        for (i=0; i<plen; i++) {
            int nc = callsign[i];
            if( nc >= 48 && nc <= 57 ) {
                nc=nc-48;
            } else if ( nc >= 65 && nc <= 90 ) {
                nc=nc-65+10;
            } else {
                nc=36;
            }
            *m=37*(*m)+nc;
        }
        *nadd=0;
        if( *m > 32768 ) {
            *m=*m-32768;
            *nadd=1;
        }
    }
    free (call6);
}

void interleave(unsigned char *sym)
{
    unsigned char tmp[162];
    unsigned char p, i, j;

    p=0;
    i=0;
    while (p<162) {
        j=((i * 0x80200802ULL) & 0x0884422110ULL) * 0x0101010101ULL >> 32;
        if (j < 162 ) {
            tmp[j]=sym[p];
            p=p+1;
        }
        i=i+1;
    }
    for (i=0; i<162; i++) {
        sym[i]=tmp[i];
    }
}

}  // namespace

int get_wspr_channel_symbols(char* rawmessage, char* hashtab, char* loctab, unsigned char* symbols) {
    int m=0, ntype=0;
    long unsigned int n=0;
    int i, j, ihash;
    unsigned char pr3[162]=
    {1,1,0,0,0,0,0,0,1,0,0,0,1,1,1,0,0,0,1,0,
        0,1,0,1,1,1,1,0,0,0,0,0,0,0,1,0,0,1,0,1,
        0,0,0,0,0,0,1,0,1,1,0,0,1,1,0,1,0,0,0,1,
        1,0,1,0,0,0,0,1,1,0,1,0,1,0,1,0,1,0,0,1,
        0,0,1,0,1,1,0,0,0,1,1,0,1,0,1,0,0,0,1,0,
        0,0,0,0,1,0,0,1,0,0,1,1,1,0,1,1,0,0,1,1,
        0,1,0,0,0,1,1,1,0,0,0,0,0,1,0,1,0,0,1,1,
        0,0,0,0,0,0,0,1,1,0,1,0,1,1,0,0,0,1,1,0,
        0,0};
    int nu[10]={0,-1,1,0,-1,2,1,0,-1,1};
    char *callsign, *grid, *powstr;
    char grid4[5], message[23];

    memset(message,0,sizeof(char)*23);
    i=0;
    while ( rawmessage[i] != 0 && i<23 ) {
        message[i]=rawmessage[i];
        i++;
    }

    size_t i1=strcspn(message," ");
    size_t i2=strcspn(message,"/");
    size_t i3=strcspn(message,"<");
    size_t i4=strcspn(message,">");
    size_t mlen=strlen(message);

    // Use the presence and/or absence of "<" and "/" to decide what
    // type of message. No sanity checks! Beware!

    if( i1 >= 3 && i1 < 7 && i2 == mlen && i3 == mlen ) {
        // Type 1 message: K9AN EN50 33
        //                 xxnxxxx xxnn nn
        callsign = strtok(message," ");
        grid = strtok(NULL," ");
        powstr = strtok(NULL," ");
        int power = atoi(powstr);
        n = pack_call(callsign);

        for (i=0; i<4; i++) {
            grid4[i]=get_locator_character_code(*(grid+i));
        }
        m = pack_grid4_power(grid4,power);

    } else if ( i3 == 0 && i4 < mlen ) {
        // Type 3:      <K1ABC> EN50WC 33
        //          <PJ4/K1ABC> FK52UD 37
        // send hash instead of callsign to make room for 6 char grid.
        // if 4-digit locator is specified, 2 spaces are added to the end.
        callsign=strtok(message,"<> ");
        grid=strtok(NULL," ");
        powstr=strtok(NULL," ");
        int power = atoi(powstr);
        if( power < 0 ) power=0;
        if( power > 60 ) power=60;
        power=power+nu[power%10];
        ntype=-(power+1);
        ihash=nhash(callsign,strlen(callsign),(uint32_t)146);
        m=128*ihash + ntype + 64;

        char grid6[7];
        memset(grid6,0,sizeof(char)*7);
        j=strlen(grid);
        for(i=0; i<j-1; i++) {
            grid6[i]=grid[i+1];
        }
        grid6[5]=grid[0];
        n = pack_call(grid6);
    } else if ( i2 < mlen ) {  // just looks for a right slash
        // Type 2: PJ4/K1ABC 37
        callsign = strtok (message," ");
        if( i2==0 || i2>strlen(callsign) ) return 0; //guards against pathological case
        powstr = strtok (NULL," ");
        int power = atoi (powstr);
        if( power < 0 ) power=0;
        if( power > 60 ) power=60;
        power=power+nu[power%10];
        int n1, ng, nadd;
        pack_prefix(callsign, &n1, &ng, &nadd);
        ntype=power + 1 + nadd;
        m=128*ng+ntype+64;
        n=n1;
    } else {
        return 0;
    }

    // pack 50 bits + 31 (0) tail bits into 11 bytes
    unsigned char it, data[11];
    memset(data,0,sizeof(char)*11);
    it=0xFF & (n>>20);
    data[0]=it;
    it=0xFF & (n>>12);
    data[1]=it;
    it=0xFF & (n>>4);
    data[2]=it;
    it= ((n&(0x0F))<<4) + ((m>>18)&(0x0F));
    data[3]=it;
    it=0xFF & (m>>10);
    data[4]=it;
    it=0xFF & (m>>2);
    data[5]=it;
    it=(m & 0x03)<<6 ;
    data[6]=it;
    data[7]=0;
    data[8]=0;
    data[9]=0;
    data[10]=0;

    // make sure that the 11-byte data vector is unpackable
    // unpack it with the routine that the decoder will use and display
    // the result. let the operator decide whether it worked.

    char *check_call_loc_pow, *check_callsign;
    check_call_loc_pow=(char *)malloc(sizeof(char)*23);
    check_callsign=(char *)malloc(sizeof(char)*13);
    signed char check_data[11];
    memcpy(check_data,data,sizeof(char)*11);

    unpk_(check_data,hashtab,loctab,check_call_loc_pow,check_callsign);
//    printf("Will decode as: %s\n",check_call_loc_pow);

    enum { nbytes=11 }; // The message with tail is packed into almost 11 bytes.
    unsigned char channelbits[nbytes*8*2]; /* 162 rounded up */
    memset(channelbits,0,sizeof(char)*nbytes*8*2);

    encode(channelbits,data,nbytes);

    interleave(channelbits);

    for (i=0; i<162; i++) {
        symbols[i]=2*channelbits[i]+pr3[i];
    }
    free(check_call_loc_pow);
    free(check_callsign);
    return 1;
}

}  // namespace wsprd::detail
