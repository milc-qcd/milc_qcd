/****************** ks_spectrum_includes.h ******************************/
/*
 *  Include files for the clover_invert application
 */

#include "../include/config.h"  /* Keep this first */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "../include/complex.h"
#include "../include/su3.h"
#include "../include/macros.h"
#include "lattice.h"
#include "../include/comdefs.h"
#include "../include/io_lat.h"
#include "../include/io_ksprop.h"
#include "../include/io_wprop.h"
#include "../include/io_ks_eigen.h"
#include "../include/generic.h"
#include "../include/generic_wilson.h"
#include "../include/generic_clover.h"
#include "../include/dirs.h"
#include "../include/io_u1lat.h"
#include "../include/generic_u1.h"
#include "../include/openmp_defs.h"
#define MULTIMASS_SET 0
#define MULTISOURCE_SET 1

#ifdef PRTIME
#define STARTTIME dtime = -dclock();
#define ENDTIME(string) dtime += dclock(); node0_printf("Aggregate time to %s %e\n",(string),dtime);
#else
#define STARTTIME
#define ENDTIME(string)
#endif


/* prototypes for functions in high level code */
int setup(void);
int readin(int prompt);


#ifdef ANISOTROPY

#ifdef FREE_KS_ANI_TEST
void free_KS_ani_test( su3_vector **src,
                       quark_source *my_ksqs,
                       ks_param my_ksp[],
                       int num_prop,
                       su3_vector **dst ) ;
#endif

#ifdef FUNNYLINKS
void funnylinks( Real umu[] );
void noreunit( void );
#endif

#endif

/* functions from ../symanzik_sl32/symanzik_sl32_includes.h */
int update();
void update_h(Real eps);
void update_u(Real eps);
void relax(int NumStp);
void monte(int NumStp);
double d_action();
double imp_gauge_action();
double hmom_action();
void make_loop_table();
void gauge_field_copy(field_offset src, field_offset dest);



/*  ani_non_prod_tests_includes.h */
