#ifndef _DEFINES_H
#define _DEFINES_H

/* Compiler macros common to all targets in this application */

/* #define SITERAND */	/* Use site-based random number generators */
#define GAUGE_FIX_TOL 1e-7

#define HOPILU 91
#define BICGILU 92

/*** The number of q and k momentum values used in the code ****/
#define MAXMOM 45
/* The maximum number of B-meson momenta p */
#define MAXPMOM 3
/**** The maximum number of kappa values ***/
#define MAX_KAPPA 10
/**** The maximum number of stored zonked light propagators ****/
/* No greater than MAX_KAPPA */
#define MAX_ZONKED_LIGHT 3
/**** The maximum number of stored zonked heavy propagators ****/
/* No greater than MAX_KAPPA */
#define MAX_ZONKED_HEAVY 5

#endif /* _DEFINES_H */
