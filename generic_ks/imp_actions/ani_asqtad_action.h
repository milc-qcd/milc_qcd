#ifndef _ASQTAD_ACTION_H
#define _ASQTAD_ACTION_H

#include "../include/dirs.h"
#include "../generic_ks/imp_actions/imp_action_types.h"
#define FERM_ACTION FN_TYPE

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
#define TADPOLE_IMPROVE	/* use tadpole improvement in quark action */
#define ASQ_OPTIMIZED_FATTENING
#define ASQ_OPTIMIZED_FORCE
#define ASQ_ACTION
#define QUARK_ACTION_DESCRIPTION "\"O(a^2): couplings(pi)=0, Naik term, No O(a^2) errors, tadpole weights\""
#  ifndef ANISOTROPY
#define MAX_BASIC_PATHS 6
#define MAX_LENGTH 7
#define MAX_NUM 688
#ifdef IMP_QUARK_ACTION_DEFINE_PATH_TABLES
    static int path_ind[MAX_BASIC_PATHS][MAX_LENGTH] = {
    { XUP, NODIR, NODIR, NODIR, NODIR, NODIR, NODIR },	/* One Link */
    { XUP, XUP, XUP, NODIR, NODIR, NODIR, NODIR },	/* Naik */
    { YUP, XUP, YDOWN, NODIR, NODIR, NODIR, NODIR },	/* Staple */
    { YUP, ZUP, XUP, ZDOWN, YDOWN, NODIR, NODIR },	/* 5-link for flavor sym. */
    { YUP, ZUP, TUP, XUP, TDOWN, ZDOWN, YDOWN},	/* 7-link for flavor sym. */
    { YUP, YUP, XUP, YDOWN, YDOWN, NODIR, NODIR },	/* 5-link compensation    */
    };
    static int path_length_in[MAX_BASIC_PATHS] = {1,3,3,5,7,5};
    static int quark_action_npaths = MAX_BASIC_PATHS ;
    static Real path_coeff[MAX_BASIC_PATHS] = {
       ( 1.0/8.0)+(6.0/16.0)+(1.0/8.0),        /* one link */
	    /*One link is 1/8 as in fat7 +3/8 for Lepage + 1/8 for Naik */
       (-1.0/24.0),	            /* Naik */
       (-1.0/8.0)*0.5,	            /* simple staple */
       ( 1.0/8.0)*0.25*0.5,         /* displace link in two directions */
       (-1.0/8.0)*0.125*(1.0/6.0),  /* displace link in three directions */
       (-1.0/16 ),                  /* Correct O(a^2) errors */
    };
#endif
#  else
#    ifndef ABSORB_ANI_XIQ
#define MAX_BASIC_PATHS 11
#define MAX_LENGTH 7
#define MAX_NUM 688
enum ani_path_type { ANI_NK = 6, ANI_LP = 10 };
#define ISO_NUM 278 /* Number of path with less than 2 anisotropic links */
#ifdef IMP_QUARK_ACTION_DEFINE_PATH_TABLES
    static int path_ind[MAX_BASIC_PATHS][MAX_LENGTH] = {
    { XUP, NODIR, NODIR, NODIR, NODIR, NODIR, NODIR },	/* One Link */
    { XUP, XUP, XUP, NODIR, NODIR, NODIR, NODIR },	/* Naik */
    { YUP, XUP, YDOWN, NODIR, NODIR, NODIR, NODIR },	/* Staple */
    { YUP, ZUP, XUP, ZDOWN, YDOWN, NODIR, NODIR },	/* 5-link for flavor sym. */
    { YUP, ZUP, TUP, XUP, TDOWN, ZDOWN, YDOWN},	/* 7-link for flavor sym. */
    { YUP, YUP, XUP, YDOWN, YDOWN, NODIR, NODIR },	/* 5-link compensation    */
    { XUP, XUP, XUP, NODIR, NODIR, NODIR, NODIR },	/* Naik w three anisotropic links */
    { YUP, XUP, YDOWN, NODIR, NODIR, NODIR, NODIR },	/* Staple w two anisotropic links */
    { YUP, ZUP, XUP, ZDOWN, YDOWN, NODIR, NODIR },	/* 5-link for flavor sym. w two anisotropic links  */
    { YUP, ZUP, TUP, XUP, TDOWN, ZDOWN, YDOWN},	/* 7-link for flavor sym. w two anisotropic links*/
    { YUP, YUP, XUP, YDOWN, YDOWN, NODIR, NODIR },	/* 5-link compensation w four anisotropic links   */
    };
    static int path_length_in[MAX_BASIC_PATHS] = {1,3,3,5,7,5,3,3,5,7,5};
    static int quark_action_npaths = MAX_BASIC_PATHS ;
    static Real path_coeff[MAX_BASIC_PATHS] = {
       ( 1.0/8.0)+(6.0/16.0)+(1.0/8.0),        /* one link */
	    /*One link is 1/8 as in fat7 +3/8 for Lepage + 1/8 for Naik */
       (-1.0/24.0),	            /* Naik */
       (-1.0/8.0)*0.5,	            /* simple staple */
       ( 1.0/8.0)*0.25*0.5,         /* displace link in two directions */
       (-1.0/8.0)*0.125*(1.0/6.0),  /* displace link in three directions */
       (-1.0/16 ),                  /* Correct O(a^2) errors */
       (-1.0/24.0),	            /* Naik w three anisotropic links */
       (-1.0/8.0)*0.5,	            /* simple staple w two anisotropic links */
       ( 1.0/8.0)*0.25*0.5,         /* displace link in two directions w two anisotropic links */
       (-1.0/8.0)*0.125*(1.0/6.0),  /* displace link in three directions w two anisotropic links */
       (-1.0/16 ),                  /* Correct O(a^2) errors w four anisotropic links */
    };
#endif
#    else
#define MAX_BASIC_PATHS 15
#define MAX_LENGTH 7
#define MAX_NUM 688
enum ani_path_type { ANI0_1L, ANI0_NK, ANI0_3L, ANI0_5L, ANI0_LP, 
                     ANI1_1L, ANI1_3L, ANI1_5L, ANI1_7L, ANI1_LP,
                     ANI2_3L, ANI2_5L, ANI2_7L, ANI3_NK, ANI4_LP };
#define ANI0_1L_MAX   6 // + 3*2 
#define ANI0_NK_MAX  12 //       + 3*2
#define ANI0_3L_MAX  36 //             + 6*4
#define ANI0_5L_MAX  84 //                   +  6*8
#define ANI0_LP_MAX 108 //                                  + 6*4
#define ANI1_1L_MAX 110 // + 1*2
#define ANI1_3L_MAX 122 //             + 3*4
#define ANI1_5L_MAX 170 //                   +  6*8
#define ANI1_7L_MAX 266 //                          +  6*16
#define ANI1_LP_MAX 278 //                                  + 3*4
#define ANI2_3L_MAX 290 //             + 3*4
#define ANI2_5L_MAX 386 //                   + 12*8
#define ANI2_7L_MAX 674 //                          + 18*16
#define ANI3_NK_MAX 676 //       + 1*2
#define ANI4_LP_MAX 688 //                                  + 3*4
#ifdef IMP_QUARK_ACTION_DEFINE_PATH_TABLES
    static int path_ind[MAX_BASIC_PATHS][MAX_LENGTH] = {
    { XUP, NODIR, NODIR, NODIR, NODIR, NODIR, NODIR },	/* One Link iso */
    { XUP, XUP, XUP, NODIR, NODIR, NODIR, NODIR },	/* Naik iso-iso-iso */
    { YUP, XUP, YDOWN, NODIR, NODIR, NODIR, NODIR },	/* Staple iso-iso-iso */
    { YUP, ZUP, XUP, ZDOWN, YDOWN, NODIR, NODIR },	/* 5-link for flavor sym. iso-iso-iso-iso-iso */
    { YUP, YUP, XUP, YDOWN, YDOWN, NODIR, NODIR },	/* 5-link compensation iso-iso-iso-iso-iso   */
    { XUP, NODIR, NODIR, NODIR, NODIR, NODIR, NODIR },	/* One Link ani */
    { YUP, XUP, YDOWN, NODIR, NODIR, NODIR, NODIR },	/* Staple iso-ani-iso */
    { YUP, ZUP, XUP, ZDOWN, YDOWN, NODIR, NODIR },	/* 5-link for flavor sym. iso-iso-ani-iso-iso */
    { YUP, ZUP, TUP, XUP, TDOWN, ZDOWN, YDOWN},	/* 7-link for flavor sym. iso-iso-iso-ani-iso-iso-iso */
    { YUP, YUP, XUP, YDOWN, YDOWN, NODIR, NODIR },	/* 5-link compensation iso-iso-ani-iso-iso   */
    { YUP, XUP, YDOWN, NODIR, NODIR, NODIR, NODIR },	/* Staple ani-iso-ani */
    { YUP, ZUP, XUP, ZDOWN, YDOWN, NODIR, NODIR },	/* 5-link for flavor sym. iso-ani-iso-ani-iso or perm */
    { YUP, ZUP, TUP, XUP, TDOWN, ZDOWN, YDOWN},	/* 7-link for flavor sym. iso-iso-ani-iso-ani-iso-iso or perm */
    { XUP, XUP, XUP, NODIR, NODIR, NODIR, NODIR },	/* Naik ani-ani-ani */
    { YUP, YUP, XUP, YDOWN, YDOWN, NODIR, NODIR },	/* 5-link compensation ani-ani-iso-ani-ani   */
    };
    static int path_length_in[MAX_BASIC_PATHS] = {1,3,3,5,5, 1,3,5,7,5, 3,5,7, 3,5};
    static int path_u0rat_pow[MAX_BASIC_PATHS] = {0,0,0,0,0, 0,0,0,0,0, 2,2,2, 2, 4};
    static int quark_action_npaths = MAX_BASIC_PATHS ;
    static Real path_coeff[MAX_BASIC_PATHS] = {
       ( 1.0/8.0)+(6.0/16.0)+(1.0/8.0),        /* one link iso */
	    /*One link is 1/8 as in fat7 +3/8 for Lepage + 1/8 for Naik */
       (-1.0/24.0),	            /* Naik iso-iso-iso */
       (-1.0/8.0)*0.5,	            /* simple staple iso-iso-iso */
       ( 1.0/8.0)*0.25*0.5,         /* displace link in two directions iso-iso-iso-iso-iso */
       (-1.0/16 ),                  /* Correct O(a^2) errors iso-iso-iso-iso-iso */
       ( 1.0/8.0)+(6.0/16.0)+(1.0/8.0),        /* one link ani */
	    /*One link is 1/8 as in fat7 +3/8 for Lepage + 1/8 for Naik */
       (-1.0/8.0)*0.5,	            /* simple staple iso-ani-iso */
       ( 1.0/8.0)*0.25*0.5,         /* displace link in two directions iso-iso-ani-iso-iso */
       (-1.0/8.0)*0.125*(1.0/6.0),  /* displace link in three directions iso-iso-iso-ani-iso-iso-iso */
       (-1.0/16 ),                  /* Correct O(a^2) errors iso-iso-ani-iso-iso */
       (-1.0/8.0)*0.5,	            /* simple staple ani-iso-ani */
       ( 1.0/8.0)*0.25*0.5,         /* displace link in two directions iso-ani-iso-ani-iso or perm */
       (-1.0/8.0)*0.125*(1.0/6.0),  /* displace link in three directions iso-iso-ani-iso-ani-iso-iso or perm */
       (-1.0/24.0),	            /* Naik ani-ani-ani */
       (-1.0/16 ),                  /* Correct O(a^2) errors ani-ani-iso-ani-ani */
    };
#endif
#    endif
#  endif
#endif /* _ASQTAD_ACTION_H */
