/****************** wi_dyn_includes.h ******************************/
/*
*  Include files for the wilson_invert application
*/

/* Include files */
#include "../include/config.h"  /* Keep this first */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../include/complex.h"
#include "../include/su3.h"
#include "jacobi.h"
#include "../include/macros.h"
#include "lattice.h"
#include "../include/comdefs.h"	/* definitions and variables for communications */
#include "../include/io_lat.h"
#include "../include/io_wprop.h"
#include "../include/generic.h"
#include "../include/generic_wilson.h"
#include "../include/generic_clover.h"
#include "../include/dirs.h"

/* prototypes for functions in high level code */
int setup_p();
int readin(int prompt);
EXTERN gauge_file *startlat_p;
EXTERN gauge_file *savelat_p;
/*  wi_inv_includes.h */

int bicgstab( /* Return value is number of iterations taken */
    field_offset src,   /* type wilson_vector (where source is to be created)*/
    field_offset dest,  /* type wilson_vector (answer and initial guess) */
    int MaxCG,          /* maximum number of iterations per restart */
    Real RsdCG,        /* desired residual - 
                           normalized as sqrt(r*r)/sqrt(src_e*src_e */
    Real *size_r,      /* resulting residual */
    int start_flag     /* 0: use a zero initial guess; 1: use dest */
    );

int congrad_t( /* Return value is number of iterations taken */
    field_offset src,   /* type wilson_vector (where source is to be created)*/
    field_offset dest,  /* type wilson_vector (answer and initial guess) */
    int MaxCG,          /* maximum number of iterations per restart */
    Real RsdCG,        /* desired residual -
                           normalized as sqrt(r*r)/sqrt(src_e*src_e */
    Real *size_r,      /* resulting residual */
    int start_flag     /* 0: use a zero initial guess; 1: use dest */
    );



void setup_offset();
void setup_links(int iflag);
void path_stuff(int iflag, int ioffset, int k, int npath[5],
 int dir[5][24][4], int sign[5][24][4], int length[5][24],
Real cc[5][24]);

void build_params(Real mass_0);
void monte_block_ape_b(int input);
#ifdef PAULI
void setup_pauli(),pauli();
#endif


int bicgstab_iter(
field_offset sol,
int niterr,
double rsqstop,
double *rsq_ptr,
 int plmin);

void delta0(field_offset src,field_offset dest,  int isign);

void path(int *dir,int *sign,int length);

void make_clov1();
void mult_ldu1(field_offset src,field_offset dest,
field_offset triang,field_offset diag, int parity);




#ifdef EIG
int f_measure2(int flag);
void boundary_flip(int sign );
#endif

void herm_delt(field_offset src,field_offset dest);

/* routines from ks_eigen... */

void dot_product(wilson_vector *vec1, wilson_vector *vec2, double_complex *dot, 
                 int parity) ;
void Matrix_Vec_mult(wilson_vector *src, wilson_vector *res, int parity) ;
void cleanup_Matrix() ;
int Rayleigh_min(wilson_vector *vec,wilson_vector **eigVec,Real Tolerance, 
                 Real RelTol,int Nvecs,int MaxIter,int Restart,int parity);
int Kalkreuter(wilson_vector **eigVec, double *eigVal, Real Tolerance, 
               Real RelTol, int Nvecs, int MaxIter, 
               int Restart, int iters, int parity) ;

void measure_chirality(wilson_vector *src, double *chirality, int parity) ;
void print_densities(wilson_vector *src, char *tag, int y,int z,int t,int parity);

void grsource();
void RotateBasis(wilson_vector **eigVec, Matrix *V, int parity) ;
