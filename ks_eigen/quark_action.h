#ifndef _ASQTAD_ACTION_H
#define _ASQTAD_ACTION_H

#include "../include/dirs.h"

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
#define MAX_BASIC_PATHS 6
#define MAX_LENGTH 7
#define MAX_NUM 688
#define TADPOLE_IMPROVE	/* use tadpole improvement in quark action */
#define ASQ_OPTIMIZED_FATTENING
#define ASQ_OPTIMIZED_FORCE
#define ASQ_ACTION
#define QUARK_ACTION_DESCRIPTION "\"O(a^2): couplings(pi)=0, Naik term, No O(a^2) errors, tadpole weights\""
#ifndef IMP_QUARK_ACTION_INFO_ONLY
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
#endif /* _ASQTAD_ACTION_H */
