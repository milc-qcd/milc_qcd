#ifndef _NAIK_ACTION_H
#define _NAIK_ACTION_H

#include "../include/dirs.h"
#include "../generic_ks/imp_actions/imp_action_types.h"
#define FERM_ACTION FN_TYPE

    /* The Naik action */
    /* Specify paths in orientation in which they appear in the
       forward part of the x component of dslash().  Rotations and
       reflections will be automatically included. Be careful
       about signs of coefficients.  See long comment at bottom
       of quark_stuff.c. */
#define QUARK_ACTION_DESCRIPTION "\"The Naik action\""
#define MAX_BASIC_PATHS 2
#define MAX_LENGTH 3
#define MAX_NUM 16
#ifdef IMP_QUARK_ACTION_DEFINE_PATH_TABLES
    static int path_ind[MAX_BASIC_PATHS][MAX_LENGTH] = {
    { XUP, NODIR, NODIR },	/* One Link */
    { XUP, XUP  , XUP   }	/* Naik */
    };
    static int path_length_in[MAX_BASIC_PATHS] = {1,3};
    static int quark_action_npaths = MAX_BASIC_PATHS ;
    static Real path_coeff[MAX_BASIC_PATHS] = {
       ( 1.0 + (1.0/8.0) ),        /* one link */
	                           /*One link is 1 + 1/8 for Naik */
       (-1.0/24.0)	           /* Naik */
    };

#endif
#endif /* _NAIK_ACTION_H */
