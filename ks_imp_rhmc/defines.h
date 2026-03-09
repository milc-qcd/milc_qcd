#ifndef _DEFINES_H
#define _DEFINES_H

/* Compiler macros common to all targets in this application */

#define SITERAND	/* Use site-based random number generators */
#define GAUGE_FIX_TOL 2e-6

#define MAX_DYN_MASSES 4
#define MAX_FPI_NMASSES 32
#define MAX_SPECTRUM_REQUEST 512
#define MAX_N_PSEUDO 10

#endif /* _DEFINES_H */

/*********************** QED U(1) *************************/
/* This can be moved to the Makefile.
 *   Added by Yuzhi Liu on 05/05/2016.
 */
//#define HAVE_U1
#define MAX_CHARGES 12
