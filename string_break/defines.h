#ifndef _DEFINES_H
#define _DEFINES_H

/* Compiler macros common to all targets in this application */
#include "../include/su3.h"

#define SITERAND	/* Use site-based random number generators */
#define MAX_SRC 128	/* maximal number of gaussian random sources */

typedef struct { su3_vector s[MAX_SRC]; } su3_vector_src;
typedef struct { su3_vector_src n[2]; } dble_su3_vec_src;

#endif /* _DEFINES_H */
