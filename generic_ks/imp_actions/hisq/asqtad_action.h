#ifndef _ASQTAD_ACTION_H
#define _ASQTAD_ACTION_H

#include "../include/umethod.h"
#include "../include/dirs.h"
#include "../generic_ks/imp_actions/imp_action_types.h"
#define FERM_ACTION HISQ

    /* The fat link action with seven link paths designed to zero
       couplings at momentum pi in any direction. The term introduced
       by Lepage to cancel the additional O(a^2) errors introduced
       by the fattening is added.  The Naik term corrects the dispersion
       relation  */
    /* Specify paths in orientation in which they appear in the
       forward part of the x component of dslash().  Rotations and
       reflections will be automatically included. Be careful
       about signs of coefficients.  See long comment at bottom
       of quark_stuff.c. */

//FOR THE MOMENT, we keep the old stuff around.  It will be superceded by
//two new sets of variables

//#define TADPOLE_IMPROVE	/* use tadpole improvement in quark action */
//#define ASQ_OPTIMIZED_FATTENING
//#define ASQ_OPTIMIZED_FORCE
#define ASQ_ACTION

#define QUARK_ACTION_DESCRIPTION "\"O(a^2): couplings(pi)=0, Naik term, No O(a^2) errors, tadpole weights\""

#define MAX_LENGTH 7	// Max length of path in any table
#define MAX_BASIC_PATHS_2 6	// Max. no. of basic paths in any table

#define NUM_BASIC_PATHS 6
#define MAX_NUM 688
#ifdef IMP_QUARK_ACTION_DEFINE_PATH_TABLES
    static int path_ind[NUM_BASIC_PATHS][MAX_LENGTH] = {
    { XUP, NODIR, NODIR, NODIR, NODIR, NODIR, NODIR },	/* One Link */
    { XUP, XUP, XUP, NODIR, NODIR, NODIR, NODIR },	/* Naik */
    { YUP, XUP, YDOWN, NODIR, NODIR, NODIR, NODIR },	/* Staple */
    { YUP, ZUP, XUP, ZDOWN, YDOWN, NODIR, NODIR },	/* 5-link for flavor sym. */
    { YUP, ZUP, TUP, XUP, TDOWN, ZDOWN, YDOWN},	/* 7-link for flavor sym. */
    { YUP, YUP, XUP, YDOWN, YDOWN, NODIR, NODIR },	/* 5-link compensation    */
    };
    static int path_length_in[NUM_BASIC_PATHS] = {1,3,3,5,7,5};
    static int quark_action_npaths = NUM_BASIC_PATHS ;
    static Real path_coeff[NUM_BASIC_PATHS] = {
       ( 1.0/8.0)+(6.0/16.0)+(1.0/8.0),        /* one link */
	    /*One link is 1/8 as in fat7 +3/8 for Lepage + 1/8 for Naik */
       (-1.0/24.0),	            /* Naik */
       (-1.0/8.0)*0.5,	            /* simple staple */
       ( 1.0/8.0)*0.25*0.5,         /* displace link in two directions */
       (-1.0/8.0)*0.125*(1.0/6.0),  /* displace link in three directions */
       (-1.0/16 ),                  /* Correct O(a^2) errors */
    };
#define INDEX_ONELINK 0
#define INDEX_NAIK 1
    static Real onelink_mass_renorm_fact = (1.0/8.0)*(-27.0/40.0);
    static Real naik_mass_renorm_fact = (-1.0/24.0)*(-27.0/40.0);
#endif


// Smearing for first level
#define QUARK_ACTION_DESCRIPTION_1 "\"No smearing\""

#ifdef IMP_QUARK_ACTION_DEFINE_PATH_TABLES
#define NUM_BASIC_PATHS_1 1
#define MAX_NUM_1 8
    static int path_ind_1[NUM_BASIC_PATHS_1][MAX_LENGTH] = {
    { XUP, NODIR, NODIR, NODIR, NODIR, NODIR, NODIR },	/* One Link */
    };
    static int path_length_in_1[NUM_BASIC_PATHS_1] = {1};
    static int quark_action_npaths_1 = NUM_BASIC_PATHS_1 ;
    static Real path_coeff_1[NUM_BASIC_PATHS_1] = {
	1.0,
    };
#endif

// Unitarization algorithm
// Choices are UNITARIZE_NONE, UNITARIZE_APE, UNITARIZE_HISQ
#define UNITARIZATION_METHOD UNITARIZE_NONE

// Smearing for second level
//#define TADPOLE_IMPROVE_2	/* use tadpole improvement in quark action */
#define HAS_NAIK_2	// second level smearing action includes 3-link term
//#define ASQ_OPTIMIZED_FATTENING_2
//#define ASQ_OPTIMIZED_FORCE_2
#define ASQ_ACTION_2
#define QUARK_ACTION_DESCRIPTION_2 "\"O(a^2): couplings(pi)=0, Naik term, No O(a^2) errors, tadpole weights\""

#ifdef IMP_QUARK_ACTION_DEFINE_PATH_TABLES
#define NUM_BASIC_PATHS_2 6
#define MAX_NUM_2 688
    static int path_ind_2[NUM_BASIC_PATHS_2][MAX_LENGTH] = {
    { XUP, NODIR, NODIR, NODIR, NODIR, NODIR, NODIR },	/* One Link */
    { XUP, XUP, XUP, NODIR, NODIR, NODIR, NODIR },	/* Naik */
    { YUP, XUP, YDOWN, NODIR, NODIR, NODIR, NODIR },	/* Staple */
    { YUP, ZUP, XUP, ZDOWN, YDOWN, NODIR, NODIR },	/* 5-link for flavor sym. */
    { YUP, ZUP, TUP, XUP, TDOWN, ZDOWN, YDOWN},	/* 7-link for flavor sym. */
    { YUP, YUP, XUP, YDOWN, YDOWN, NODIR, NODIR },	/* 5-link compensation    */
    };
    static int path_length_in_2[NUM_BASIC_PATHS_2] = {1,3,3,5,7,5};
    static int quark_action_npaths_2 = NUM_BASIC_PATHS_2 ;
    static Real path_coeff_2[NUM_BASIC_PATHS_2] = {
       ( 1.0/8.0)+(6.0/16.0)+(1.0/8.0),        /* one link */
	    /*One link is 1/8 as in fat7 +3/8 for Lepage + 1/8 for Naik */
       (-1.0/24.0),	            /* Naik */
       (-1.0/8.0)*0.5,	            /* simple staple */
       ( 1.0/8.0)*0.25*0.5,         /* displace link in two directions */
       (-1.0/8.0)*0.125*(1.0/6.0),  /* displace link in three directions */
       (-1.0/16 ),                  /* Correct O(a^2) errors */
    };
#endif


#endif /* _ASQTAD_ACTION_H */
