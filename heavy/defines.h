#ifndef _DEFINES_H
#define _DEFINES_H

/* Compiler macros common to all targets in this application */

/* No random number generator used in this application, so save space */
/* #define SITERAND */	/* Use site-based random number generators */
/* Uses hack to distinguish single and double precision files */
#define VERSION_NUMBER (59354 + 8*(sizeof(Real) - 4))

#define GAUGE_FIX_TOL 1.0e-7 /* For gauge fixing */

enum meson_io_options {SAVE_MESON_BINARY =11 , SAVE_MESON_ASCII  , FORGET_MESON }  ; 

#define NCHANNELS 4  /* number of meson channels */

#endif /* _DEFINES_H */

