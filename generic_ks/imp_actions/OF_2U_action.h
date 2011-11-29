#ifndef _OF_2U_ACTION_H
#define _OF_2U_ACTION_H

#include "../include/dirs.h"
#include "../generic_ks/imp_actions/imp_action_types.h"
#define FERM_ACTION FN_TYPE
    /* The one-link plus staple(fat) action unitarized
	to second order*/
    /* Specify paths in orientation in which they appear in the
       forward part of the x component of dslash().  Rotations and
       reflections will be automatically included. Be careful
       about signs of coefficients.  See long comment at bottom
       of quark_stuff.c. */

#define QUARK_ACTION_DESCRIPTION "\"one-link + staple(fat) action (second order unitarized)\""
#define MAX_BASIC_PATHS 13
#define MAX_LENGTH 9
#define MAX_NUM 10000	/* ??? */
#define ALPHA 0.20
#ifdef IMP_QUARK_ACTION_DEFINE_PATH_TABLES
    static int path_ind[MAX_BASIC_PATHS][MAX_LENGTH] = {
    { XUP, NODIR, NODIR, NODIR, NODIR, NODIR, NODIR, NODIR, NODIR },

    { YUP, XUP, YDOWN, NODIR, NODIR, NODIR, NODIR, NODIR, NODIR },
    { XUP, YUP, XDOWN, YDOWN, XUP, NODIR, NODIR, NODIR, NODIR },

    { YUP, XUP, YDOWN, XDOWN, YDOWN, XUP, YUP, NODIR, NODIR },
    { YUP, XUP, YDOWN, XDOWN, YUP, XUP, YDOWN, NODIR, NODIR },
    { YUP, XUP, YDOWN, XDOWN, ZUP, XUP, ZDOWN, NODIR, NODIR },

    /* renormalizes one link */
    { YUP, XUP, YDOWN, ZUP, XDOWN, ZDOWN, XUP, NODIR, NODIR },
    { YUP, XUP, YDOWN, YDOWN, XDOWN, YUP, XUP, NODIR, NODIR },

    /* renormalizes one link */
    { XUP, YUP, XDOWN, YDOWN, YDOWN, XUP, YUP, NODIR, NODIR },
    { XUP, YUP, XDOWN, YDOWN, ZUP, XUP, ZDOWN, NODIR, NODIR },

    { XUP, YUP, XDOWN, YDOWN, XUP, YUP, XDOWN, YDOWN, XUP },
    { XUP, YUP, XDOWN, YDOWN, XUP, YDOWN, XDOWN, YUP, XUP },
    { XUP, YUP, XDOWN, YDOWN, XUP, ZUP, XDOWN, ZDOWN, XUP }
    };
    static int path_length_in[MAX_BASIC_PATHS] = 
	{1, 3,5, 7,7,7, 7,7, 7,7, 9,9,9};
    static int quark_action_npaths = MAX_BASIC_PATHS;
    static Real path_coeff[MAX_BASIC_PATHS] = {
        1.0 - 2*6*0.125*ALPHA*ALPHA,	/* one link */
       -0.5*ALPHA,      /* staple(fat) */
        0.5*ALPHA,       /* staple(fat) traversed backwards */
		/* minus sign for backwards, and one for enclosing
		   plaquette */
	 0.125*ALPHA*ALPHA,
	 0.125*ALPHA*ALPHA,
	 0.125*ALPHA*ALPHA,

	/* renormalizes one link */
	-0.125*ALPHA*ALPHA,
	-0.125*ALPHA*ALPHA,

	/* renormalizes one link */
	-0.125*ALPHA*ALPHA,
	-0.125*ALPHA*ALPHA,

	 0.125*ALPHA*ALPHA,
	 0.125*ALPHA*ALPHA,
	 0.125*ALPHA*ALPHA
    };

#endif
#endif /* _OF_2U_ACTION_H */
