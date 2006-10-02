/****** fermion_force_asqtad3_rhmc2.c  -- ******************/
/* MIMD version 7 */
/* D.T. 12/05 First try at RHMC version
 *
 */

// TEMP: define NEWFORCE, REVFORCE, or nothing
#define NEWFORCE

#include "ks_imp_includes.h"	/* definitions files and prototypes */

void eo_fermion_force_rhmc( int alg_flag, Real eps, int order, Real mass, Real *residues,
  Real *roots, su3_vector **multi_x, field_offset phi_off ){
    // alg_flag passes info for integration algorithms that work differently
    // at different time steps

    Real final_rsq;
    int i,j;
    site *s;
    int stoporder = 0;

    if( alg_flag==0 ){
	stoporder=order;
    } else if ( alg_flag < 0){
	stoporder=order+alg_flag;  // skip some force terms
    } else if ( alg_flag > 0 ){    // weight of some force terms X 3
	stoporder = order;
	for( i=0; i<alg_flag; i++)if(order-i>0)residues[order-i] *= 3.0;
    }

    // Compute ( M^\dagger M)^{-1} in xxx_even
    // Then compute M*xxx in temporary vector xxx_odd 
    /* See long comment at end of file */
	/* The diagonal term in M doesn't matter */
    ks_ratinv( phi_off, multi_x, mass, roots, order, niter,
	md_rsqmin, EVEN, &final_rsq );

    for(j=0;j<stoporder;j++){ dslash_field( multi_x[j], multi_x[j],  ODD ); }

    // NEWFORCE is force routine with sum over terms factored out (transport
    // 3x3 matrices around)
    // REVFORCE is same thing, with indices on CG solutions reversed to maybe
    // improve cache hits
    // OLDFORCE uses old multiterm force until number of terms left is less than
    // VECLENGTH, then uses two term or one term force routines to finish up
    j=1; // j counts orders done (order 0 is trivial constant term in rat. function)
#ifdef NEWFORCE
    fn_fermion_force_rhmc( eps, &(residues[j]), multi_x+j-1, stoporder ); //TRY_ONE
    j+= stoporder;
#elif defined REVFORCE
    fn_fermion_force_rhmc_reverse( eps, &(residues[j]), multi_x+j-1, stoporder ); //TRY_ONE
    j+= stoporder;
#elif defined OLDFORCE

    while( j<= stoporder-VECLENGTH+1 ){
      eo_fermion_force_multi( eps, &(residues[j]), multi_x+j-1, VECLENGTH ); //TRY_ONE
      j+=VECLENGTH;
    }
    for( ;j<=stoporder-1;j+=2){
        FOREVENSITES(i,s){
            s->xxx1 = multi_x[j-1][i] ;
            s->xxx2 = multi_x[j  ][i] ;
 }
        dslash_site( F_OFFSET(xxx1), F_OFFSET(xxx1), ODD );
        dslash_site( F_OFFSET(xxx2), F_OFFSET(xxx2), ODD );
        eo_fermion_force_twoterms( eps, residues[j], residues[j+1],
     F_OFFSET(xxx1), F_OFFSET(xxx2) );
    }
    for( ; j<=stoporder; j++ ){
        FOREVENSITES(i,s){ s->xxx1 = multi_x[j-1][i] ; }
        dslash_site( F_OFFSET(xxx1), F_OFFSET(xxx1), ODD );
        eo_fermion_force_oneterm( eps, residues[j], F_OFFSET(xxx1) );
    }

    if ( alg_flag > 0 ){    // weight of some force terms X 3
        for( i=0; i<alg_flag; i++)if(order-i>0)residues[order-i] /= 3.0;
    }
#endif
}

