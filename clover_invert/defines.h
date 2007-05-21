#ifndef _DEFINES_H
#define _DEFINES_H

/* Compiler macros common to all targets in this application */

/* No random numbers in this application, so we save space by commenting out */

#define MAXHTMP 8

/*#define SITERAND*//* Use site-based random number generators */

#define GAUGE_FIX_TOL 1.0e-7 /* For gauge fixing */
#define MAX_SPECTRUM_REQUEST 512

#endif /* _DEFINES_H */
