#ifndef _P4_ACTION_H
#define _P4_ACTION_H

#include "../include/dirs.h"

    /* The Bielefeld "P4" action */
    /* Specify paths in orientation in which they appear in the
       forward part of the x component of dslash().  Rotations and
       reflections will be automatically included. Be careful
       about signs of coefficients.  See long comment at bottom
       of quark_stuff.c. */
#define MAX_BASIC_PATHS 4
#define MAX_LENGTH 3
#define MAX_NUM 152
#ifdef FN
BOMB THE COMPILE p4 is not a "FAT" action
#endif
#define QUARK_ACTION_DESCRIPTION "\"Bielefeld P4, 1+3 paths\""

#ifndef IMP_QUARK_ACTION_INFO_ONLY
    static int path_ind[MAX_BASIC_PATHS][MAX_LENGTH] = {
    { XUP, NODIR, NODIR },
    { XUP, YUP, YUP },
    { YUP, YUP, XUP },
    { YUP, XUP, YUP }
    };
    static int path_length_in[MAX_BASIC_PATHS] = {1,3,3,3};
    static int quark_action_npaths = MAX_BASIC_PATHS;
    static Real path_coeff[MAX_BASIC_PATHS] = {
        3.0/4.0,     /* one link */
       1.0/72.0,    /* three link "XXY" */
       1.0/72.0,    /* three link "YYX" */
       -1.0/72.0    /* three link "YXY", sign from crossing plaq. */
    };

#endif
#endif /* _P4_ACTION_H */
