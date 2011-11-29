#ifndef _ONELINK_ACTION_H
#define _ONELINK_ACTION_H

#include "../include/dirs.h"
#include "../generic_ks/imp_actions/imp_action_types.h"
#define FERM_ACTION FN_TYPE

    /* Include file for the conventional "one link" action */
    /* Specify paths in orientation in which they appear in the
       forward part of the x component of dslash().  Rotations and
       reflections will be automatically included. Be careful
       about signs of coefficients.  See long comment at bottom
       of quark_stuff.c. */
#define QUARK_ACTION_DESCRIPTION "\"Single link action\""

#define MAX_BASIC_PATHS 1
#define MAX_LENGTH 1
#define MAX_NUM 8
#ifdef IMP_QUARK_ACTION_DEFINE_PATH_TABLES
    static int path_ind[MAX_BASIC_PATHS][MAX_LENGTH] = {
    { XUP }
    };
    static int path_length_in[MAX_BASIC_PATHS] = {1};
    static int quark_action_npaths = MAX_BASIC_PATHS;
    static Real path_coeff[MAX_BASIC_PATHS] = {
        1.0	/* one link */
    };
#endif
#endif /* _ONELINK_ACTION_H */
