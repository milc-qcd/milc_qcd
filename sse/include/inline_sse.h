#ifndef _INLINE_SSE_H_
#define _INLINE_SSE_H_

/* These assembly inline versions work only with the gnu C compiler
   on Opteron processors */

#include "inline_headers.h"

/* Corrections by A. Alexandru */
typedef struct
{
   unsigned int c1,c2,c3,c4;
} sse_mask __attribute__ ((__aligned__ (16)));

static sse_mask _sse_sgn24 __attribute__ ((__aligned__ (16), __unused__)) ={0x00000000, 0x80000000, 0x00000000, 0x80000000};
static sse_mask _sse_sgn2 __attribute__  ((__aligned__ (16), __unused__)) ={0x00000000, 0x80000000, 0x00000000, 0x00000000};
static sse_mask _sse_sgn4 __attribute__  ((__aligned__ (16), __unused__)) ={0x00000000, 0x00000000, 0x00000000, 0x80000000};


#endif
