#ifndef _DEFINES_H
#define _DEFINES_H

/* Compiler macros common to all targets in this application */

#define SITERAND	/* Use site-based random number generators */
#define GAUGE_FIX_TOL 2e-6 /* For gauge fixing */

enum boundary_choice { PERIODIC_EVERYWHERE = 10 , ANTI_PERIODIC_IN_TIME  }  ;
enum invert_quark_options { LOCAL_QUARK_SRC , WALL_QUARK_SRC   } ;


/** a useful debug macro **/
#define IF_VERBOSE_ON(ff)  if(this_node==0 && verbose_flag >= ff )

#endif /* _DEFINES_H */
