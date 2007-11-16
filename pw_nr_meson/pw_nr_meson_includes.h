/********************** pw_nr_meson_includes.h **************************
*									*
*  This header is included in all codes in this directory               *
*  MIMD version 7 							*
*									*
*/
#include "../include/config.h"  /* Keep this first */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "../include/complex.h"
#include "../include/su3.h"
#include "lattice.h"
#include "pauli_prop.h"
#include "../include/macros.h"
#include "../include/comdefs.h"
#include "../include/io_lat.h"
#include "../include/io_ksprop.h"
#include "../include/io_wprop.h"
#include "../include/generic.h"
#include "../include/generic_ks.h"
#include "../include/generic_wilson.h"
#include "../include/generic_clover.h"
#include "../include/dirs.h"
#include "../include/generic_quark_types.h"

#ifdef PRTIME
#define STARTTIME dtime = -dclock();
#define ENDTIME(string) dtime += dclock(); node0_printf("Time to %s %e\n",(string),dtime);
#else
#define STARTTIME
#define ENDTIME(string)
#endif


/* prototypes for functions in high level code */

int setup();

int readin(int prompt);

void init_pw_prop(void);

void all_pw_prop(int dir1, int dir2);

void assemble_pw_prop(void);

void av_pw_prop_out(void);

void print_pw_prop(int ksrc, int ksnk);

void rotate_w_quark(field_offset, field_offset, float);

void smear_quark(int ksnk, int dir1, int dir2);

void load_smearing(wilson_quark_source source_wqs[], 
		   wilson_quark_source sink_wqs[], int n);
void free_smearing(int n);

void make_pwave_source(wilson_quark_source source_wqs[], int ksrc, int dir);

int get_wprop_to_field(int startflag, char startfile[], 
		      int saveflag, char savefile[],
		      block_pauli_propagator *prop,
		      wilson_quark_source *wqs,
		      quark_invert_control *qic,
		      dirac_clover_param *dcp);

