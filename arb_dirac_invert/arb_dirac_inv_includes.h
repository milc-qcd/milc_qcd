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

void rotate_field(field_offset src,field_offset dest, Real mass_0);

void w_meson_p(field_offset src1,field_offset src2,
Real prop[10][MAX_P][MAX_NT] );
void w_baryon_p(field_offset src1,field_offset src2,field_offset src3,
Real prop[4][MAX_P][MAX_NT] );
 


#ifdef CVC
void cvc(int source, int sink,
field_offset src1,field_offset src2,
Real cvc_prop[MAX_P][MAX_NT], Real hvc_prop[MAX_P][MAX_NT] );



void hax(field_offset src1,field_offset src2,Real hax_prop[MAX_P][MAX_NT]) ;
void cvc_pauli(field_offset src1,field_offset src2) ;
#endif

#ifdef EIG
int f_measure2(int flag);
void boundary_flip(int sign );
#endif
