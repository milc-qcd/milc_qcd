/******************** prefetch.c ***************************************/

#include "../include/config.h"
#ifdef HAVE_64_BYTE_CACHELINE

#include "prefetch64.c"

#else

#include "prefetch32.c"

#endif  /* HAVE_64_BYTE_CACHELINE */

/* prefetch.c */ 
