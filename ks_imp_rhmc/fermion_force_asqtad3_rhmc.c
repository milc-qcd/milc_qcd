/****** fermion_force_asqtad3_rhmc2.c  -- ******************/
/* MIMD version 7 */
/* D.T. 12/05 First try at RHMC version
 * C.D. 10/06 Moved some code to generic_ks/fermion_force_fn_multi.c
 */


#include "ks_imp_includes.h"	/* definitions files and prototypes */

void eo_fermion_force_rhmc( int alg_flag, Real eps, params_ratfunc *rf, 
			    su3_vector **multi_x, field_offset phi_off ){
    // alg_flag passes info for integration algorithms that work differently
    // at different time steps

    Real final_rsq;
    int i,j;
    int stoporder = 0;
    int order = rf->order;
    Real *residues = rf->res;
    Real *roots = rf->pole;

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
    ks_ratinv( phi_off, multi_x, roots, order, niter,
	md_rsqmin, EVEN, &final_rsq );

    for(j=0;j<stoporder;j++){ dslash_field( multi_x[j], multi_x[j],  ODD ); }

    eo_fermion_force_multi( eps, &(residues[1]), multi_x, stoporder );

    if ( alg_flag > 0 ){    // weight of some force terms X 3
        for( i=0; i<alg_flag; i++)if(order-i>0)residues[order-i] /= 3.0;
    }
}

