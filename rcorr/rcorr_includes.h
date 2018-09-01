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

/* linearlsq2 */
double linearlsq(double *m, double *sdm, double *b, double *sdb,
		 double x[], double y[], double sd[], int n);

/* print_corr.c */
void
print_result(Real *q[], Real *q2[], int nb, int rb[]);

/* rcorr.c */
void
rcorr(Real *qblock[], Real *qblock2[], 
      complex *qin_sloppy[], int nrand_sloppy, 
      complex *qin_diff[], int nrand_diff,
      int nblock, int rand_block[]);

/* symmetrize.c */
void
symmetrize(Real *q, Real *q2);

/* rcorr_includes.h */
