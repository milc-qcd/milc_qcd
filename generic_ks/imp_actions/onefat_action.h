#ifndef _ONEFAT_ACTION_H
#define _ONEFAT_ACTION_H

#include "../include/dirs.h"
#include "../generic_ks/imp_actions/imp_action_types.h"
#define FERM_ACTION FN_TYPE

    /* The one-link plus staple(fat) action */
    /* Specify paths in orientation in which they appear in the
       forward part of the x component of dslash().  Rotations and
       reflections will be automatically included. Be careful
       about signs of coefficients.  See long comment at bottom
       of quark_stuff.c. */

#define QUARK_ACTION_DESCRIPTION "\"one-link + staple(fat) action\""
#define MAX_BASIC_PATHS 2
#define MAX_LENGTH 3
#define MAX_NUM 56
#define ALPHA 0.25
#ifdef IMP_QUARK_ACTION_DEFINE_PATH_TABLES
    static int path_ind[MAX_BASIC_PATHS][MAX_LENGTH] = {
    { XUP, NODIR, NODIR },
    { YUP, XUP, YDOWN }
    };
    static int path_length_in[MAX_BASIC_PATHS] = {1,3};
    static int quark_action_npaths = MAX_BASIC_PATHS;
    static Real path_coeff[MAX_BASIC_PATHS] = {
        1.0/(1.0+6.0*ALPHA),	/* one link */
       - ALPHA/(1.0+6.0*ALPHA)     /* staple(fat) */
    };

#endif
#endif /* _ONEFAT_ACTION_H */
