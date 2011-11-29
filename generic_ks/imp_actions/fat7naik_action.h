#ifndef _FAT7NAIK_ACTION_H
#define _FAT7NAIK_ACTION_H

#include "../include/dirs.h"
#include "../generic_ks/imp_actions/imp_action_types.h"
#define FERM_ACTION FN_TYPE

    /* The fat link action with seven link paths designed to zero
       couplings at momentum pi in any direction. 
       The Naik term corrects the dispersion relation  */
    /* Specify paths in orientation in which they appear in the
       forward part of the x component of dslash().  Rotations and
       reflections will be automatically included. Be careful
       about signs of coefficients.  See long comment at bottom
       of quark_stuff.c. */
#define QUARK_ACTION_DESCRIPTION "\"Fat7-Naik: couplings(pi)=0 plus Naik term\""
#define MAX_BASIC_PATHS 5
#define MAX_LENGTH 7
#define MAX_NUM 688
#ifdef IMP_QUARK_ACTION_DEFINE_PATH_TABLES
    static int path_ind[MAX_BASIC_PATHS][MAX_LENGTH] = {
    { XUP, NODIR, NODIR, NODIR, NODIR, NODIR, NODIR },	/* One Link */
    { XUP, XUP, XUP, NODIR, NODIR, NODIR, NODIR },	/* Naik */
    { YUP, XUP, YDOWN, NODIR, NODIR, NODIR, NODIR },	/* Staple */
    { YUP, ZUP, XUP, ZDOWN, YDOWN, NODIR, NODIR },	/* 5-link for flavor sym. */
    { YUP, ZUP, TUP, XUP, TDOWN, ZDOWN, YDOWN}	/* 7-link for flavor sym. */
    };
    static int path_length_in[MAX_BASIC_PATHS] = {1,3,3,5,7};
    static int quark_action_npaths = MAX_BASIC_PATHS ;
    static Real path_coeff[MAX_BASIC_PATHS] = {
       ( 1.0/8.0)+(1.0/8.0),        /* one link */
	    /*One link is 1/8 as in fat7 + 1/8 for Naik */
       (-1.0/24.0),	            /* Naik */
       (-1.0/8.0)*0.5,	            /* simple staple */
       ( 1.0/8.0)*0.25*0.5,         /* displace link in two directions */
       (-1.0/8.0)*0.125*(1.0/6.0)  /* displace link in three directions */
    };

#endif
#endif /* _FAT7NAIK_ACTION_H */
