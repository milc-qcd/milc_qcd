                          MILC Version 7

This code was developed by the MILC collaboration for simulations of
SU3 lattice gauge theory on MIMD parallel machines.  This code is
publicly available for research purposes.  Publications of work
done using this code or derivatives of this code should acknowledge
this use.  Development of this code was supported in part by grants
from the US Department of Energy and National Science Foundation.
Use this code at your own risk.

Since this is our working code, it is continually in a state of
development.  We will informally support the code as best we can by
answering questions and fixing bugs.  We will be very grateful for
reports of problems and suggestions for improvements.  These may be
sent to 

   detar@physics.utah.edu.

      or 

   doug@klingon.physics.arizona.edu 


Target architectures:

Currently code is supposed to run on:
1. Any scalar machine ("vanilla" version)
2. Linux/Unix switched clusters via MPI
3. MPP GigE mesh architectures

Overview of the code:

Each "application", or major variant of the code, has its own
directory.  For example, the ks_imp_dyn applications directory,
contains code for simulating full QCD in the staggered fermion scheme.
The compilation requires code from the "libraries", "generic", and
"generic_ks" directories, as well as the "ks_dynamical" directory
itself.  All applications share the "libraries" directory containing
low-level stuff, and the "generic" directory containing high level
stuff that is more or less independent of the physics.  The various
staggered fermion applications share the "generic_ks" directory.
Examples of "generic" code are the random number routines, the lattice
layout routines, routines to evaluate the plaquette or Polyakov loop,
etc.  Of the shared code, only the libraries must be built separately
before building any application code.  Code in the various generic
directories is compiled automatically as needed, and should not be
compiled separately.

doc:  More detailed documentation of the code

libraries:	Low level routines and include files
	
	complex.h:	Definitions and macros for complex numbers
	su3.h		Definitions and macros for SU(3)
	complex.1.a 	Routines for complex numbers (single precision)
	complex.2.a 	Routines for complex numbers (double precision)
	su3.1.a		Routines for SU(3) operations (single precision)
	su3.2.a		Routines for SU(3) operations (double precision)

include:   Header files required by the code

generic: High level code for generic SU(3) simulation. The other
	directories, which are for real applications, should use 
	the routines in this directory where possible, otherwise 
	copy routines from this directory and modify them or write 
	new routines.

generic_ks:     High level code shared by staggered fermion applications.

generic_wilson:   High level code shared by Wilson fermion applications.


The remaining directories are "applications" directories.  Most of the
application code can be built in various ways, depending on the
requirements of the project.  For example, the ks_imp_dyn directory
contains code for simulating full QCD in the staggered fermion scheme.
This code can be built to use the "R", "phi", and hybrid Monte Carlo
algorithms, and may also include measurements of the hadron spectrum.
These variants are obtained by selecting the appropriate compilation
target, as indicated in the application Make_template file.

arb_overlap
       Computes eigenvalues and eigenvectors of the overlap operator.

clover_dynamical
	Simulations with improved dynamical Wilson fermions.  Variants
	include the "R", "phi" and hybrid Monte Carlo updating
	algorithms.  Measurements include plaquette, Polyakof loop,
	psi-bar-psi and fermion energy and pressure.  Optional
	measurements include hadron spectrum, screening spectrum,
	axial current quark mass, and Landau gauge quark propagators.

clover_invert2
        Calculation of meson and baryon propagators with improved
        Wilson fermions or naive Dirac fermions.

ext_src
	Utility for extracting an extended source from a propagator file.

file_utilities
       A variety of utilities for converting lattice file formats, running
       a checksum of a gauge file, and  comparing some binary files.

gauge_utilities
        A variety of utilities for manipulating the gauge field, including
        coordinate translations, gauge fixing, and boundary twists.

gluon_prop
        Calculation of the gluon propagator in a specific gauge and
        the nonperturbative renormalization of the vector and axial
        vector current.

hvy_qpot
       Measures static quark potential as a function of separation.
       Also a variety of Wilson loops.

ks_eigen
	Eigenvalues of the staggered Dirac operator.

ks_imp_dyn
	Simulations with dynamical Kogut-Susskind fermions.  Variants
	include the "R", "phi" and hybrid Monte Carlo updating
	algorithms.  Measurements include plaquette, Polyakov loop,
	psi-bar-psi, and fermions energy and pressure.  Optional
	measurements include hadron spectrum, screening spectrum,
	and some wave functions. (includes FFT routines in wave function)

ks_imp_rhmc
        Simulations as above, but using the rational hybrid Monte
        Carlo method.

ks_imp_utilities
      Test code for the staggered fermion force and staggered inverter.

ks_measure
	Scalar operators, mainly used for the equation of state ahd
	quark number suscptibilities.

ks_spectrum
	Calculation of meson and baryon spectra from a wide variety of
	sources and sinks.  

pure_gauge
	Simulation of the pure gauge theory with the plaquette gauge action.

schroed_cl_inv
        Schroedinger functional computations with improved Wilson
        fermions.

smooth_inst
	Topological charge.

symanzik_sl32
	Pure gauge theory with the symanzik, 1-loop improved action.

Each directory should contain a README file which gives some further
information.

How to build the code:

1. Start with the particular application on a particular machine and
   edit the Makefile according to instructions there.

2. Edit the libraries/Make_vanilla make file as needed.

3. Select a particular make target by consulting the Make_template file
   in the desired applications directory.

4. Sample input and output files are provided for most targets.
   It is a good idea to run the compiled code with the sample input
   file and compare the results with the sample output.

Good luck!

References:

For documentation, see the directory doc.

Please also see the README files in various subdirectories.




