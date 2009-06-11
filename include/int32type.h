/* Some architectures have 64-bit ints.  Our binary file formats use
   32 bit integers.  So we define a type se we can get 32-bit integers
   on all platforms.  For 64-bit architectures, the code must be
   compiled with -DSHORT32.  */
/*
   7/26/01 Changed name of macro to make it more obvious CD
   4/17/98 Added an unsigned u_int32type for checksums CD
   2/26/98 Changed int32type to signed CD
*/

#ifndef _TYPE32_H
#define _TYPE32_H

#include "../include/config.h"

/* One and only one should be defined */
#if defined(SHORT_IS_32BIT) && defined(INT_IS_32BIT)
MAKE UP YOUR MIND!!  SEE config.h
#endif

#if !defined(SHORT_IS_32BIT) && !defined(INT_IS_32BIT)
MAKE UP YOUR MIND!!  SEE config.h
#endif

#ifdef SHORT_IS_32BIT
typedef short int32type;
typedef unsigned short u_int32type;

#else

typedef int int32type;
typedef unsigned int u_int32type;

#endif

typedef unsigned long long u_int64type;

#endif /* _TYPE32_H */
