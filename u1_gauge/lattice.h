/* ************************************************************	*/
/*			    LATTICE_U1G.H		       	*/
/* 								*/
/* Defines lattice parameters and declares lattice structure	*/
/* for pure U(1) gauge fields.					*/
/*								*/	
/* Header file for the main file control.c			*/
/*								*/
/* Last Updated on 05.01.07					*/
/*								*/
/* ************************************************************	*/
#ifndef _LATTICE_H
#define _LATTICE_H

#include "defines.h"
#include "../include/io_u1lat.h"
#include "../include/macros.h"
#include "../include/complex.h"
#include "../include/random.h"
#include "../include/su3.h"

#ifdef CONTROL
#define EXTERN
#else
#define EXTERN extern
#endif

/* Lattice structure for U(1) gauge fields */
typedef struct{

	int x,y,z,t;		/* coordinates of this site */
	char parity;		/* site is even or odd */
	int index;		/* index of the site */

        double_prn site_prn;	/* site random number info */
        int space1;

	complex u1tmp, loop, gftmp2;

	Real mom[4];		/* lattice momenta */

	      } site;
EXTERN site *lattice;

/* Global Lattice Parameters */
EXTERN int nx,ny,nz,nt; 	/* lattice dimensions */
EXTERN int volume;		/* nx*ny*nz*nt */
EXTERN int *latin;		/* lattice site index */

EXTERN Real echarge;            /* electron charge == e */
EXTERN Real *u1_A;              /* global u1 field (vector potential) */
EXTERN complex *u1gf;           /* temporary global u1 field (complex vector potential) */
EXTERN Real g_splaq,g_tplaq;	/* global plaquette measures */

EXTERN double g_ssplaq,g_stplaq;	/* global plaquette measures */
EXTERN double_complex linktrsum;	/* not used */
EXTERN  u_int32type nersc_checksum;	/* not used */

EXTERN int this_node;		/* node number of this node */
EXTERN int sites_on_node;	/* = volume */
EXTERN int number_of_nodes;	/* number of nodes in use */
EXTERN int even_sites_on_node;	/* = volume / 2 */
EXTERN int odd_sites_on_node;	/* = volume / 2 */

EXTERN int start_u1flag;	/* begin u(1) lattice: FRESH or RELOAD */
EXTERN char start_u1file[MAXFILENAME];
EXTERN gauge_file *start_u1lat_p;
EXTERN int save_u1flag;		/* end u(1) lattice: SAVE or FORGET */
EXTERN char save_u1file[MAXFILENAME];
EXTERN gauge_file *save_u1lat_p;

EXTERN int iseed;		/* random number seed */
EXTERN double_prn node_prn;	/* nodes' random number structure */
EXTERN int junk_id;		/* junk number, for checking */

/* pointers for accessing and gathering fields */
#define N_POINTERS 8
EXTERN char **gen_pt[N_POINTERS];

#endif /* _LATTICE_H */

/* ************************************************************	*/
