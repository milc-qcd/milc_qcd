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
#include "../include/io_u1lat.h"
#include "../include/generic_u1.h"

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

/* setup.c */
//int setup(void);
//int readin(int prompt);


/*  ks_spectrum_includes.h */
