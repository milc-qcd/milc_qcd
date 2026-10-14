This application is used to generate U(1) gauge configurations. In quantum electrodynamics, the momentum space distribution of the photon fields $A_\mu$ is Gaussian, so they can be generated independently and trivially (no need for a Markov chain) in momentum space. They are then converted to position space using a fast Fourier transform. Several gauge fixing choices are available including Lorenz, Coulomb, and Feynman.

Some key points:

* This application generates photon fields $A_\mu$, which can later be converted to U(1) link variables $U_\mu=e^{ieA_\mu}$

* Must be linked with FFTW library

A typical build command is:

```
MY_CC=cc \
MY_CXX=CC \
COMPILER="gnu" \
OPT="-O3 -Ofast -g" \
LDFLAGS="-g" \
PRECISION=2 \
WANTFFTW=true \
make -j 1 u1_g-Coulomb-QED_TL
```

## Calling the generator from another application

The field generation itself lives in `generic_u1/u1_gauge_field.c` and is
usable from any MILC application.  The caller owns the storage for the
`u1_A` field:

```c
#include "../include/generic_u1.h"

initialize_u1_gauge_field(U1_COULOMB, U1_QED_L);  /* once, after make_lattice() */

Real *A = create_u1_A_field();                    /* 4*sites_on_node Reals */
for(i = 0; i < nconf; i++){
  generate_u1_gauge_field(A);                     /* A[4*i+dir] */
  ...
}
destroy_u1_A_field(A);

finalize_u1_gauge_field();
```

`initialize_u1_gauge_field()` builds the momentum table and configures the
four-dimensional Fourier transform; neither depends on the random numbers,
so it is done only once no matter how many configurations are generated.

The gauge choice (`U1_FEYNMAN`, `U1_LORENZ`, `U1_COULOMB`) and the treatment
of the modes with vanishing spatial momentum (`U1_QED_NONE`, `U1_QED_L`,
`U1_QED_TL`) are run-time arguments, so no compilation macros are needed to
call this routine.  Add `u1_gauge_field.o` to the application's list of
`generic_u1` objects and link with FFTW.

The choice of algorithm is still the compile-time macro `GAUGE_METHOD`.  Left
undefined, it selects the general algorithm, which is the only one that
honors the two arguments above; `U1_MILC_ORIGINAL` and `U1_HATTON_REVISION`
select the two legacy algorithms, which have their gauge and scheme built in.

The `test` directory contains sample input files.

To run the tests, try:

```
MY_CC=cc \
MY_CXX=CC \
COMPILER="gnu" \
WANTFFTW=true \
make check
```

Check the output for errors and instances of "NOT OK".
