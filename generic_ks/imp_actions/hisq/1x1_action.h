#ifndef _1X1_ACTION_H
#define _1X1_ACTION_H


#include "../include/umethod.h"
#include "../include/dirs.h"

    /* Specify paths in orientation in which they appear in the
       forward part of the x component of dslash().  Rotations and
       reflections will be automatically included. Be careful
       about signs of coefficients.  See long comment at bottom
       of quark_stuff.c. */

#define MAX_LENGTH 1	// Maximum length of path in any path table
#define MAX_BASIC_PATHS 1  // Max. no. of basic paths in any path table

// Smearing for first level
#define NUM_BASIC_PATHS_1 1
#define MAX_NUM_1 8
#define QUARK_ACTION_DESCRIPTION_1 "\"1 link\""
#ifndef IMP_QUARK_ACTION_INFO_ONLY
    static int path_ind_1[NUM_BASIC_PATHS_1][MAX_LENGTH] = {
    { XUP },	/* One Link */
    };
    static int quark_action_npaths_1 = NUM_BASIC_PATHS_1 ;
    static int path_length_in_1[NUM_BASIC_PATHS_1] = {1};
    static Real path_coeff_1[NUM_BASIC_PATHS_1] = {
       1.0,        /* one link */
    };
#endif

//#define UNITARIZATION_METHOD UNITARIZE_NONE
#define UNITARIZATION_METHOD UNITARIZE_ROOT

// Smearing for second level
#define NUM_BASIC_PATHS_2 1
#define MAX_NUM_2 8
#define QUARK_ACTION_DESCRIPTION_2 "\"1 link\""
#ifndef IMP_QUARK_ACTION_INFO_ONLY
    static int path_ind_2[NUM_BASIC_PATHS_2][MAX_LENGTH] = {
    { XUP },	/* One Link */
    };
    static int path_length_in_2[NUM_BASIC_PATHS_2] = {1};
    static int quark_action_npaths_2 = NUM_BASIC_PATHS_2 ;
    static Real path_coeff_2[NUM_BASIC_PATHS_2] = {
       1.0,        /* one link */
    };
#endif

//FOR THE MOMENT, we keep the old stuff around.  It will be superceded by
//two new sets of variables

#define NUM_BASIC_PATHS 1
#define MAX_NUM 8
#define QUARK_ACTION_DESCRIPTION "\"1 Link\""
#ifndef IMP_QUARK_ACTION_INFO_ONLY
    static int path_ind[NUM_BASIC_PATHS][MAX_LENGTH] = {
    { XUP },	/* One Link */
    };
    static int path_length_in[NUM_BASIC_PATHS] = {1};
    static int quark_action_npaths = NUM_BASIC_PATHS ;
    static Real path_coeff[NUM_BASIC_PATHS] = {
	1.0,
    };
#endif

#endif // _1X1_ACTION_H
