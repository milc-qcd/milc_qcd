/****************** rcorr_includes.h ******************************/
/*
 *  Include files for the clover_invert application
 */

#include "../include/config.h"  /* Keep this first */
#include "../include/precision.h"
#include "params.h"
#include "lattice.h"
#include "../include/complex.h"
#include "../include/su3.h"
#include "../include/generic.h"
#include "../include/io_scidac.h"
#include "../include/comdefs.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include <qio.h>
#include <fftw3.h>

/* setup.c */
int setup(void);
int readin(int prompt);

/* accumulate_density.c */
void 
accumulate_current_density(char *filename, complex *qin[], 
			   Real charge, Real *mass, int nrand);
/* print_corr.c */
void
print_result(Real *q, int nrand);

/* rcorr.c */
Real *
rcorr(complex *qin[], int nrand);

/* symmetrize.c */
void
symmetrize(Real *q);

/* rcorr_includes.h */
