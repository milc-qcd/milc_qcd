#ifndef _P4FAT3_ACTION_H
#define _P4FAT3_ACTION_H

#include "../include/dirs.h"
#include "../generic_ks/imp_actions/imp_action_types.h"
#define FERM_ACTION EO_TYPE

    /* The Bielefeld "P4fat3" action */
    /* Specify paths in orientation in which they appear in the
       forward part of the x component of dslash().  Rotations and
       reflections will be automatically included. Be careful
       about signs of coefficients.  See long comment at bottom
       of quark_stuff.c. */
#ifdef FN
BOMB THE COMPILE p4fat3 is not a "FAT" action
#endif
#define QUARK_ACTION_DESCRIPTION "\"Bielefeld P4fat3, 1 staple + (2-1) paths\""

#define MAX_BASIC_PATHS 4
#define MAX_LENGTH 3
#define MAX_NUM 152
#ifdef IMP_QUARK_ACTION_DEFINE_PATH_TABLES
    static int path_ind[MAX_BASIC_PATHS][MAX_LENGTH] = {
    { XUP, NODIR, NODIR },	/* One Link */
    { YUP, XUP, YDOWN },	/* Staple */
    { XUP, YUP, YUP },
    { YUP, YUP, XUP }
    };
    static int path_length_in[MAX_BASIC_PATHS] = {1,3,3,3};
    static int quark_action_npaths = MAX_BASIC_PATHS;
    static Real path_coeff[MAX_BASIC_PATHS] = {
        (3.0/4.0)/2.2,		/* one link */
        -(3.0/4.0)*0.2/2.2,	/* staple (minus sign from staggered phases) */
        1.0/48.0,		/* three link "XXY" */
        1.0/48.0		/* three link "YYX" */
    };

#endif
#endif /* _P4FAT3_ACTION_H */
