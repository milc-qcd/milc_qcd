#ifndef _DEFINES_H
#define _DEFINES_H

/* Compiler macros common to all targets in this application */

#ifdef RANDOM 
#define SITERAND 	/* Use site-based random number generators */
#endif
#define GAUGE_FIX_TOL 1.0e-7 /* For gauge fixing */

/* definitions to switch setup_links */
#define SIMPLE 0
#define MAX_P 5
#define MAX_NT 24
#define NLINK 40

#endif /* _DEFINES_H */

