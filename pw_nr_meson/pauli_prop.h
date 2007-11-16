#ifndef _PAULI_PROP_H
#define _PAULI_PROP_H

#include "../include/su3.h"

typedef struct { su3_vector d[2]; } pauli_vector;
typedef struct { pauli_vector d[2]; } spin_pauli_vector;
typedef struct { spin_pauli_vector c[3]; } pauli_propagator;
typedef struct { 
  pauli_propagator up;
  pauli_propagator dn;
} block_pauli_propagator;

#endif /* _PAULI_PROP_H */
