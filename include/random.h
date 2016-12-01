#ifndef _RANDOM_H
#define _RANDOM_H

#include "../include/precision.h"

/* random number structures */

typedef struct {
  /* We assume long is at least 32 bits and long long at least 64 (checked at initialization) */
  unsigned long r0,r1,r2,r3,r4,r5,r6;
  unsigned long long multiplier,addend,ic_state;
  float scale;
  //Real gauss_gset; int gauss_toggle; //state information from gaussian_random_number
} double_prn;

/* Generic random number generator returning a uniformly distributed
   random value on [0,1] */
Real myrand(double_prn *prn_pt);

#endif /* _RANDOM_H */
