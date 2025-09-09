This application is used to generate U(1) gauge configurations. In quantum electrodynamics, the momentum space distribution of the photon fields $A_\mu$ is Gaussian, so they can be generated independently and trivially (no need for a Markov chain) in momentum space. They are then converted to position space using a fast Fourier transform. Several gauge fixing choices are available including Lorenz, Coulomb, and Feynman.

Some key points:

* This application generates photon fields $A_\mu$, which can later be converted to U(1) link variables $U_\mu=e^{ieA_\mu}$

* This application only supports scalar running (i.e. no MPI). Fortunately, it has a small memory footprint and is very fast compared to SU(3) gauge generation

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
