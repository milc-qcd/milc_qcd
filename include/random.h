#ifndef _RANDOM_H
#define _RANDOM_H

#include "../include/precision.h"

/* random number structures */

typedef struct {
  /* We assume long is at least 32 bits */
  unsigned long r0,r1,r2,r3,r4,r5,r6;
  unsigned long multiplier,addend,ic_state;
  float scale;
} double_prn;

/* Generic random number generator returning a uniformly distributed
   random value on [0,1] */
Real myrand(double_prn *prn_pt);

#endif /* _RANDOM_H */
