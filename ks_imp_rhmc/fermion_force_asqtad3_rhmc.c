/* OBSOLETE. RENAMED eo_fermion_force_rhmc.c */
/****** fermion_force_asqtad3_rhmc2.c  -- ******************/
/* MIMD version 7 */
/* D.T. 12/05 First try at RHMC version
 * C.D. 10/06 Moved some code to generic_ks/fermion_force_fn_multi.c
 * D.T. 3/07 Removed alg_flag stuff
 */


#include "ks_imp_includes.h"	/* definitions files and prototypes */

void eo_fermion_force_rhmc( Real eps, params_ratfunc *rf, 
			    su3_vector **multi_x, field_offset phi_off,
			    Real my_rsqmin, int my_niter, int cg_prec,
			    int ff_prec, ferm_links_t *fn )
{
    // at different time steps

    Real final_rsq;
    int j;
    int order = rf->order;
    Real *residues = rf->res;
    Real *roots = rf->pole;

    // Compute ( M^\dagger M)^{-1} in xxx_even
    // Then compute M*xxx in temporary vector xxx_odd 
    /* See long comment at end of file */
	/* The diagonal term in M doesn't matter */
    load_ferm_links(&fn_links);
    ks_ratinv( phi_off, multi_x, roots, order, my_niter,
	       my_rsqmin, cg_prec, EVEN, &final_rsq, fn );

    for(j=0;j<order;j++){ dslash_field( multi_x[j], multi_x[j],  ODD, fn); }

    eo_fermion_force_multi( eps, &(residues[1]), multi_x, order, ff_prec, fn );
}

