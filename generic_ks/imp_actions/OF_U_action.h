#ifndef _OF_U_ACTION_H
#define _OF_U_ACTION_H

#include "../include/dirs.h"
#include "../generic_ks/imp_actions/imp_action_types.h"
#define FERM_ACTION FN_TYPE

    /* The one-link plus staple(fat) action unitarized
	to lowest order*/
    /* Specify paths in orientation in which they appear in the
       forward part of the x component of dslash().  Rotations and
       reflections will be automatically included. Be careful
       about signs of coefficients.  See long comment at bottom
       of quark_stuff.c. */

#define QUARK_ACTION_DESCRIPTION "\"one-link + staple(fat) action (first order unitarized)\""
#define MAX_BASIC_PATHS 3
#define MAX_LENGTH 5
#define MAX_NUM 104	/* 8 + 48 + 48 */
#define ALPHA 0.5
#define FUDGE (1.0/(1+ALPHA*ALPHA))
#ifdef IMP_QUARK_ACTION_DEFINE_PATH_TABLES
    static int path_ind[MAX_BASIC_PATHS][MAX_LENGTH] = {
    { XUP, NODIR, NODIR, NODIR, NODIR },
    { YUP, XUP, YDOWN, NODIR, NODIR },
    { XUP, YUP, XDOWN, YDOWN, XUP }
    };
    static int path_length_in[MAX_BASIC_PATHS] = {1,3,5};
    static int quark_action_npaths = MAX_BASIC_PATHS;
    static Real path_coeff[MAX_BASIC_PATHS] = {
        1.0*FUDGE,	/* one link */
       -0.5*ALPHA*FUDGE,      /* staple(fat) */
        0.5*ALPHA*FUDGE       /* staple(fat) traversed backwards */
		/* minus sign for backwards, and one for enclosing
		   plaquette */
    };

#endif
#endif /* _OF_U_ACTION_H */
