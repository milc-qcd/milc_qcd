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
#include "../include/io_ks_eigen.h"
#include "../include/generic.h"
#include "../include/generic_wilson.h"
#include "../include/generic_clover.h"
#include "../include/dirs.h"
#include "../include/io_u1lat.h"
#include "../include/generic_u1.h"

#ifdef GB_BARYON
#include "../include/gb_ops.h"
#ifdef BLIND
#include "../include/blind_data.h"
#endif
#endif

#ifdef PRTIME
#define STARTTIME dtime = -dclock();
#define ENDTIME(string) dtime += dclock(); node0_printf("Aggregate time to %s %e\n",(string),dtime);  fflush(stdout);
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

int solve_ksprop(enum set_type set_type, enum inv_type inv_type,
		 int num_prop, int startflag[], char startfile[][MAXFILENAME],
		 int saveflag[], char savefile[][MAXFILENAME],
		 ks_prop_field *ksprop[],
		 ks_prop_field *source[],
                 quark_source *my_ksqs[],
		 quark_invert_control my_qic[],
		 ks_param my_ksp[],
		 Real charge,
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
#ifdef GB_BARYON
void spectrum_ks_gb_baryon(ks_prop_field **qko0, ks_prop_field **qko1, ks_prop_field **qko2,
  su3_matrix *links, int triplet);
/* gb_baryon_snk.c */
void gb_baryon(ks_prop_field **qko0, ks_prop_field **qko1, ks_prop_field **qko2,
              su3_matrix *links, enum gb_baryon_op *src_op,
              enum gb_baryon_op *snk_op,
              int stIdx, short *dowall, short *docube, int num_d, int num_s, int *r0,
              int *mom, char *par, complex *momfld, int *flip_snk,
              int num_corr_gb, int *phase, Real *fact, complex **prop);


#endif

/*  ks_spectrum_includes.h */
