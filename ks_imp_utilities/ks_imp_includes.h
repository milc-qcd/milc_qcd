/****************** ks_imp_includes.h ******************************/
/*
*  Include files for Kogut-Susskind dynamical improved action application
*/

/* Include files */
#include "../include/config.h"  /* Keep this first */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "../include/complex.h"
#include "../include/su3.h"
#include "lattice.h"
#include "../include/macros.h"
#include "../include/comdefs.h"
#include "../include/io_lat.h"
#include "../include/generic_ks.h"
#include "../include/generic.h"
#include "../include/dirs.h"
#include "../include/dirs.h"
#ifdef HAVE_QIO
#include "../include/io_scidac.h"
#include "../include/io_scidac_ks.h"
#endif
#ifdef HAVE_QDP
#include <qdp.h>
#endif

#ifdef FN
#define dslash dslash_fn
#endif
#ifdef EO
#define dslash dslash_eo
#endif

/* prototypes for functions in high level code */
int setup();
int readin(int prompt);

int ask_color_vector( int prompt, int *flag, char *filename );
int ask_color_matrix( int prompt, int *flag, char *filename );
void check_fermion_force( char *srcfile, int srcflag, field_offset src, 
			  char *ansfile, int ansflag, Real mass);
void check_ks_invert( char *srcfile, int srcflag, field_offset src, 
		      char *ansfile, int ansflag, field_offset ans, 
		      field_offset tmp, Real mass);
void check_invert2( field_offset src, field_offset dest, 
		    field_offset temp, Real mass,
		    Real tol, int parity);

char *create_QCDML();
void free_QCDML(char *qcdml);

