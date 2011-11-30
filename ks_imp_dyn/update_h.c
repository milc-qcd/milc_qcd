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
  int ff_prec = PRECISION;  /* Just use prevailing precision for now */
#ifdef FN
//    free_fn_links(&fn_links);
//    free_fn_links(&fn_links_dmdu0);
  invalidate_fermion_links(fn_links);
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
    eo_fermion_force_oneterm_site( eps, ((Real)nflavors)/4., F_OFFSET(xxx),
				   ff_prec, fn_links );
#else
/**
   eo_fermion_force_oneterm( eps, ((Real)nflavors1)/4., F_OFFSET(xxx1),
        ff_prec, fn_links );
    eo_fermion_force_oneterm( eps, ((Real)nflavors2)/4., F_OFFSET(xxx2),
        ff_prec, fn_links );
**/
/**/
    eo_fermion_force_twoterms_site( eps, ((Real)nflavors1)/4., 
				    ((Real)nflavors2)/4., F_OFFSET(xxx1), 
				    F_OFFSET(xxx2), ff_prec, fn_links );
    /**/
#endif
} /* update_h */


