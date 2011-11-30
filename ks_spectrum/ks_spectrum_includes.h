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
#include "../include/io_ksprop.h"
#include "../include/io_lat.h"
#include "../include/io_wprop.h"
#include "../include/generic.h"
#include "../include/generic_wilson.h"
#include "../include/generic_clover.h"
#include "../include/dirs.h"

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
int num_mes_report(void);

/* ksprop_info.c */
char *create_ks_XML(void);

/* ks_source_info.c */
char *create_kss_XML(char *filename, quark_source *ksqs);

/* make_prop.c */
void read_ksprop_to_ksp_field(int startflag, char startfile[], 
			      quark_source *my_ksqs, ks_prop_field *ksp);

int solve_ksprop(int num_prop, int startflag[], char startfile[][MAXFILENAME],
		 int saveflag[], char savefile[][MAXFILENAME],
		 ks_prop_field *ksprop[],
		 quark_source *my_ksqs,
		 quark_invert_control my_qic[],
		 ks_param my_ksp[],
		 Real bdry_phase[],
		 int r0[4],
		 int check);

void dump_ksprop_from_ksp_field(int saveflag, char savefile[], 
				ks_prop_field *ksp);
ks_prop_field *reread_ksprop_to_ksp_field(int saveflag, char savefile[], int nc);

/* setup.c */
int setup(void);
int readin(int prompt);

/* spectrum_ks.c */

int ask_corr_file( FILE *fp, int prompt, int *flag, char* filename);
void spectrum_ks(ks_prop_field *qp0, int naik_index0, ks_prop_field *qp1, int naik_index1, int pair);
void spectrum_ks_baryon(ks_prop_field *qp0, ks_prop_field *qp1, ks_prop_field *qp2, int triplet);

/*  ks_spectrum_includes.h */









