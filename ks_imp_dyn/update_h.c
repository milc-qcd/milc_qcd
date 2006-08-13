/****** update_h.c  -- ******************/
/* updates momentum matrices for improved action */
/* D.T. & J.H., naik term    8/96
*  D.T., fat link fermion term 5/97
*  D.T. general quark action 1/98
*  D.T. two types of quarks 3/99
*  T.D. and A.H. improved gauge updating spliced in 5/97
*
* MIMD version 7 */

#include "ks_imp_includes.h"	/* definitions files and prototypes */

void update_h( Real eps ){
#ifdef FN
    free_longlinks();
    free_fatlinks();
#endif
    /* gauge field force */
    rephase(OFF);
    imp_gauge_force(eps,F_OFFSET(mom));
    rephase(ON);
    /* fermionic force */
    /* First compute M*xxx in temporary vector xxx_odd */
    /* See long comment at end of file */
	/* The diagonal term in M doesn't matter */
#ifdef ONEMASS
    eo_fermion_force_oneterm( eps, ((Real)nflavors)/4., F_OFFSET(xxx) );
#else
/**
    eo_fermion_force_oneterm( eps, ((Real)nflavors1)/4., F_OFFSET(xxx1) );
    eo_fermion_force_oneterm( eps, ((Real)nflavors2)/4., F_OFFSET(xxx2) );
**/
/**/
    eo_fermion_force_twoterms( eps, ((Real)nflavors1)/4., 
      ((Real)nflavors2)/4., F_OFFSET(xxx1), F_OFFSET(xxx2) );
/**/
#endif
} /* update_h */


