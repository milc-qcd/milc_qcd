/****** update_h.c  -- ******************/
/* updates momentum matrices for improved action */
/* D.T. & J.H., naik term    8/96
*  D.T., fat link fermion term 5/97
*  D.T. general quark action 1/98
*  D.T. two types of quarks 3/99
*  T.D. and A.H. improved gauge updating spliced in 5/97
*
* MIMD version 6 */

#include "ks_imp_includes.h"	/* definitions files and prototypes */

void update_h( Real eps ){
    /* gauge field force */
    rephase(OFF);
    imp_gauge_force(eps,F_OFFSET(mom));
    rephase(ON);
    /* fermionic force */
    /* First compute M*xxx in temporary vector xxx_odd */
    /* See long comment at end of file */
	/* The diagonal term in M doesn't matter */
    dslash( F_OFFSET(xxx1), F_OFFSET(xxx1), ODD );
    dslash( F_OFFSET(xxx2), F_OFFSET(xxx2), ODD );
/**
    eo_fermion_force( eps, nflavors1, F_OFFSET(xxx1) );
    eo_fermion_force( eps, nflavors2, F_OFFSET(xxx2) );
**/
/**/
    eo_fermion_force_3f( eps, nflavors1, F_OFFSET(xxx1),
      nflavors2, F_OFFSET(xxx2) );
/**/
} /* update_h */


