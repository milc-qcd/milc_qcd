#ifndef _DEFINES_H
#define _DEFINES_H

/* Compiler macros common to all targets in this application */

/* #define SITERAND */	/* Use site-based random number generators */
#define GAUGE_FIX_TOL 1.0e-7 /* For gauge fixing */
#define VERSION_NUMBER 59354  /* For our own binary meson propagator file */

/** Parameters required for the variational code ****/

enum meson_io_options {SAVE_MESON_BINARY =11 , SAVE_MESON_ASCII  , FORGET_MESON }  ; 

#define NCHANNELS 4  /* number of meson channels */

/* some macros that are useful for the static variational code **/

/** memory location in the variational matrix ***/
#define VM_PT  t + nt*(aloop + nosmear*(bloop)) 

/*** memory location for the smeared meson operators ******/
/** This is represented in a fortran code as 
    complex avmeson(nc,nd,nc,nd,nt,nosmear)
**********/
#define MESON_WHERE colour + 3*(spin+4*(ic+3*(ispin+4*(t+nt*ismear))))

#define IF_MASTER  if(this_node==0)

enum  smear_type { LOCALsmear, WALLsmear , EXPON , twoSpoly , threeSpoly    } ;


#endif /* _DEFINES_H */


