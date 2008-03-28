#ifndef _DEFINES_H
#define _DEFINES_H

/* Compiler macros common to all targets in this application */

/* No random numbers in this application, so we save space by commenting out */

#define MAXHTMP 8

/*#define SITERAND*//* Use site-based random number generators */

#define GAUGE_FIX_TOL 1.0e-7 /* For gauge fixing */

/* defines for 3rd nearest neighbor (NAIK) stuff */
#define X3UP 8
#define Y3UP 9
#define Z3UP 10
#define T3UP 11
#define T3DOWN 12
#define Z3DOWN 13
#define Y3DOWN 14
#define X3DOWN 15
 
#define OPP_3_DIR(dir) (23-(dir))
#define DIR3(dir) ((dir)+8)

#endif /* _DEFINES_H */
