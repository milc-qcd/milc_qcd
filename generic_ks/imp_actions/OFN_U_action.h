#ifndef _OFN_U_ACTION_H
#define _OFN_U_ACTION_H

#include "../include/dirs.h"

    /* The one-link plus Naik plus staple(fat) action unitarized
	to lowest order*/
    /* Specify paths in orientation in which they appear in the
       forward part of the x component of dslash().  Rotations and
       reflections will be automatically included. Be careful
       about signs of coefficients.  See long comment at bottom
       of quark_stuff.c. */
#define MAX_BASIC_PATHS 4
#define MAX_LENGTH 5
#define MAX_NUM 112	/* 8 + 8 + 48 + 48 */


    static int path_ind[MAX_BASIC_PATHS][MAX_LENGTH] = {
    { XUP, NODIR, NODIR, NODIR, NODIR },
    { XUP, XUP, XUP, NODIR, NODIR },
    { YUP, XUP, YDOWN, NODIR, NODIR },
    { XUP, YUP, XDOWN, YDOWN, XUP }
    };
    static int path_length_in[MAX_BASIC_PATHS] = {1,3,3,5};
    static int quark_action_npaths = MAX_BASIC_PATHS;
    static Real path_coeff[MAX_BASIC_PATHS] = {
        (9.0/8.0),	/* one link */
       -(9.0/8.0)/27.0,     /* three link (Naik) */
       -(9.0/8.0)*0.25,      /* staple(fat) */
        (9.0/8.0)*0.25      /* staple(fat) traversed backwards */
		/* minus sign for backwards, and one for enclosing
		   plaquette */
    };
    static char quark_action_description[] =
      "\"one-link + Naik + staple(fat) action (first order unitarized)\"";
#endif /* _OFN_U_ACTION_H */
