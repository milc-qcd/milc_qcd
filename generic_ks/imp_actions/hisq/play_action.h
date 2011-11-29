#ifndef _3X3_ACTION_H
#define _3X3_ACTION_H

#include "../include/umethod.h"
#include "../include/dirs.h"
#include "../generic_ks/imp_actions/imp_action_types.h"
#define FERM_ACTION HISQ

    /* Specify paths in orientation in which they appear in the
       forward part of the x component of dslash().  Rotations and
       reflections will be automatically included. Be careful
       about signs of coefficients.  See long comment at bottom
       of quark_stuff.c. */

#define MAX_LENGTH 5	// Maximum length of path in any path table
#define MAX_BASIC_PATHS 9  // Max. no. of basic paths in any path table

#define Clink_1		(1.0)	// single link coeff in first smearing
#define Cstaple_1	(-0.1)	// staple coeff in first smearing
#define Clink_2		(1.0)	// single link coeff in second smearing
#define Cstaple_2	(-0.1)	// staple coeff in second smearing

// Smearing for first level
#define QUARK_ACTION_DESCRIPTION_1 "\"1 link plus 3-staple\""
#define NUM_BASIC_PATHS_1 3
#define MAX_NUM_1 56
#ifdef IMP_QUARK_ACTION_DEFINE_PATH_TABLES
    static int path_ind_1[NUM_BASIC_PATHS_1][MAX_LENGTH] = {
    { XUP, NODIR, NODIR, NODIR, NODIR },	/* One Link */
    { YUP, XUP, YDOWN, NODIR, NODIR },	/* Staple */
    { XUP, XUP, XDOWN, NODIR, NODIR },	/* Fake link */
    };
    static int quark_action_npaths_1 = NUM_BASIC_PATHS_1 ;
    static int path_length_in_1[NUM_BASIC_PATHS_1] = {1,3,3};
    static Real path_coeff_1[NUM_BASIC_PATHS_1] = {
       Clink_1,        /* one link */
       Cstaple_1,     /* simple staple */
       Cstaple_1,     /* simple staple */
    };
#endif

// Unitarization algorithm
// Choices are UNITARIZE_NONE, UNITARIZE_APE, UNITARIZE_HISQ
#define UNITARIZATION_METHOD UNITARIZE_NONE

// Smearing for second level
#define QUARK_ACTION_DESCRIPTION_2 "\"1 link plus 3-staple\""

#ifdef IMP_QUARK_ACTION_DEFINE_PATH_TABLES
#define NUM_BASIC_PATHS_2 2
#define MAX_NUM_2 56
    static int path_ind_2[NUM_BASIC_PATHS_2][MAX_LENGTH] = {
    { XUP, NODIR, NODIR, NODIR, NODIR },	/* One Link */
    { YUP, XUP, YDOWN, NODIR, NODIR },	/* Staple */
    };
    static int path_length_in_2[NUM_BASIC_PATHS_2] = {1,3};
    static int quark_action_npaths_2 = NUM_BASIC_PATHS_2 ;
    static Real path_coeff_2[NUM_BASIC_PATHS_1] = {
       Clink_2,        /* one link */
       Cstaple_2,     /* simple staple */
    };
#endif

//FOR THE MOMENT, we keep the old stuff around.  It will be superceded by
//two new sets of variables

#define QUARK_ACTION_DESCRIPTION "\"3-staple smearing iterated twice\""

#ifdef IMP_QUARK_ACTION_DEFINE_PATH_TABLES
#define NUM_BASIC_PATHS 9
#define MAX_NUM (8+6*8+ 6*8+6*8+ 4*6*8+4*6*8+4*6*8+ 6*8+8) 
    static int path_ind[NUM_BASIC_PATHS][MAX_LENGTH] = {
    { XUP,   NODIR, NODIR, NODIR, NODIR },	/* One Link */
    { YUP,   XUP,   YDOWN, NODIR, NODIR },	/* Staple */

    { XDOWN, YUP,   XUP,   XUP,   YDOWN },
    { YUP,   XUP,   XUP,   YDOWN, XDOWN },

    { YUP,   ZUP,   YDOWN, XUP,   ZDOWN },
    { YUP,   XUP,   ZUP,   YDOWN, ZDOWN },
    { YUP,   ZUP,   XUP,   ZDOWN, YDOWN },

    { YUP,   YUP,   XUP,   YDOWN, YDOWN },

    { XUP,   XUP,   XUP,   NODIR, NODIR },	// Naik
    };
    static int path_length_in[NUM_BASIC_PATHS] = {1,3,5,5,5,5,5,5,3};
    static int quark_action_npaths = NUM_BASIC_PATHS ;
    static Real path_coeff[NUM_BASIC_PATHS] = {
	Clink_1*Clink_2 + 6*3*Cstaple_1*Cstaple_2,
	Clink_1*Cstaple_2 + Cstaple_1*Clink_2,

	Cstaple_1*Cstaple_2,
	Cstaple_1*Cstaple_2,

	Cstaple_1*Cstaple_2,
	Cstaple_1*Cstaple_2,
	Cstaple_1*Cstaple_2,

	Cstaple_1*Cstaple_2,

	0.0,	// Naik
    };
#endif

#endif // _3X3_ACTION_H
