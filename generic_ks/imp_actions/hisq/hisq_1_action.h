#ifndef _HISQ_1_ACTION_H
#define _HISQ_1_ACTION_H

#include "../include/umethod.h"
#include "../include/dirs.h"

    /* Specify paths in orientation in which they appear in the
       forward part of the x component of dslash().  Rotations and
       reflections will be automatically included. Be careful
       about signs of coefficients.  See long comment at bottom
       of quark_stuff.c. */

#define MAX_LENGTH 7	// Maximum length of path in any path table
#define MAX_BASIC_PATHS 688  // Max. no. of basic paths in any path table

// Smearing for first level
// This is Fat7
#define NUM_BASIC_PATHS_1 5
#define MAX_NUM_1 640
//#define ASQ_OPTIMIZED_FATTENING_1
//#define ASQ_OPTIMIZED_FORCE_1
#define QUARK_ACTION_DESCRIPTION_1 "\"Fat 7 (level 1)\""
#ifndef IMP_QUARK_ACTION_INFO_ONLY
    static int path_ind_1[NUM_BASIC_PATHS_1][MAX_LENGTH] = {
    { XUP, NODIR, NODIR, NODIR, NODIR, NODIR, NODIR },  /* One Link */
    { XUP, XUP, XUP, NODIR, NODIR, NODIR, NODIR },      /* Naik */
    { YUP, XUP, YDOWN, NODIR, NODIR, NODIR, NODIR },    /* Staple */
    { YUP, ZUP, XUP, ZDOWN, YDOWN, NODIR, NODIR },      /* 5-link for flavor sym. */
    { YUP, ZUP, TUP, XUP, TDOWN, ZDOWN, YDOWN}, /* 7-link for flavor sym. */
    };
    static int quark_action_npaths_1 = NUM_BASIC_PATHS_1 ;
    static int path_length_in_1[NUM_BASIC_PATHS_1] = {1,3,3,5,7};
    static Real path_coeff_1[NUM_BASIC_PATHS_1] = {
       (( 1.0/8.0)+(1.0/8.0)),        /* one link */
            /*One link is 1/8 as in fat7 + 1/8 for Naik */
       (-1.0/24.0),                 /* Naik */
       (-1.0/8.0)*0.5,              /* simple staple */
       ( 1.0/8.0)*0.25*0.5,         /* displace link in two directions */
       (-1.0/8.0)*0.125*(1.0/6.0),  /* displace link in three directions */
    };
#define INDEX_ONELINK 0
#define INDEX_NAIK 1
    static Real onelink_mass_renorm_fact = (1.0/8.0)*(-27.0/40.0);
    static Real naik_mass_renorm_fact = (-1.0/24.0)*(-27.0/40.0);
#endif

// Unitarization algorithm
// Choices are UNITARIZE_NONE, UNITARIZE_APE, UNITARIZE_HISQ
#define UNITARIZATION_METHOD UNITARIZE_APE

// Smearing for second level
#define NUM_BASIC_PATHS_2 6
#define MAX_NUM_2 688
//#define ASQ_OPTIMIZED_FATTENING_2
//#define ASQ_OPTIMIZED_FORCE_2
#define QUARK_ACTION_DESCRIPTION_2 "\"Fat7 + 2xLepage\""
#ifndef IMP_QUARK_ACTION_INFO_ONLY
    static int path_ind_2[NUM_BASIC_PATHS_2][MAX_LENGTH] = {
    { XUP, NODIR, NODIR, NODIR, NODIR, NODIR, NODIR },  /* One Link */
    { XUP, XUP, XUP, NODIR, NODIR, NODIR, NODIR },      /* Naik */
    { YUP, XUP, YDOWN, NODIR, NODIR, NODIR, NODIR },    /* Staple */
    { YUP, ZUP, XUP, ZDOWN, YDOWN, NODIR, NODIR },      /* 5-link for flavor sym. */
    { YUP, ZUP, TUP, XUP, TDOWN, ZDOWN, YDOWN}, /* 7-link for flavor sym. */
    { YUP, YUP, XUP, YDOWN, YDOWN, NODIR, NODIR },      /* 5-link compensation    */
    };
    static int path_length_in_2[NUM_BASIC_PATHS_2] = {1,3,3,5,7,5};
    static int quark_action_npaths_2 = NUM_BASIC_PATHS_2 ;
    static Real path_coeff_2[NUM_BASIC_PATHS_2] = {
       (( 1.0/8.0)+(2.0*6.0/16.0)+(1.0/8.0)),        /* one link */
            /*One link is 1/8 as in fat7 + 2*3/8 for Lepage + 1/8 for Naik */
       (-1.0/24.0),                 /* Naik */
       (-1.0/8.0)*0.5,              /* simple staple */
       ( 1.0/8.0)*0.25*0.5,         /* displace link in two directions */
       (-1.0/8.0)*0.125*(1.0/6.0),  /* displace link in three directions */
       (-2.0/16 ),                  /* Correct O(a^2) errors, 2X as much as Asqtad  */
    };
#endif

//FOR THE MOMENT, we keep the old stuff around.  It has been  superceded by
//two new sets of variables
// For this action, the old stuff is meaningless, so leave it set to trivial

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

#endif // _HISQ_1_ACTION_H