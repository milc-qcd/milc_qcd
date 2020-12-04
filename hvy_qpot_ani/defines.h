#ifndef _DEFINES_H
#define _DEFINES_H

/* Compiler macros common to all targets in this application */
#ifdef FIX_GAUGE_FIX_STEPS
#define GAUGE_FIX_TOL 0
#define GAUGE_FIX_STEPS FIX_GAUGE_FIX_STEPS
#else
#define GAUGE_FIX_TOL 2e-6
#define GAUGE_FIX_STEPS 500
#endif

/* #define SITERAND */	/* Use site-based random number generators */
#endif /* _DEFINES_H */
