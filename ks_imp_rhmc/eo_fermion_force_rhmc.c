/****** eo_fermion_force_rhmc.c  -- ******************/
/* MIMD version 7 */
/* D.T. 12/05 First try at RHMC version
 * C.D. 10/06 Moved some code to generic_ks/fermion_force_fn_multi.c
 * D.T. 3/07 Removed alg_flag stuff
 */


#include "ks_imp_includes.h"	/* definitions files and prototypes */

/* NOT SUPPORTED FOR HISQ! */

void eo_fermion_force_rhmc( Real eps, params_ratfunc *rf, 
			    su3_vector **multi_x, field_offset phi_off, 
			    Real my_rsqmin, int my_niter,
			    int cg_prec, int ff_prec, fermion_links_t *fl)
{
    // at different time steps

    Real final_rsq;
    int j;
    int order = rf->order;
    Real *residues = rf->res;
    Real *roots = rf->pole;
    imp_ferm_links_t **fn;

    // Compute ( M^\dagger M)^{-1} in xxx_even
    // Then compute M*xxx in temporary vector xxx_odd 
    /* See long comment at end of file */
	/* The diagonal term in M doesn't matter */
    //    load_ferm_links(fn);
    restore_fermion_links_from_site(fl, cg_prec);
    fn = get_fm_links(fl);

    /* Do the inversion for zero Naik term epsilon */
    ks_ratinv( phi_off, multi_x, roots, order, my_niter,
	       my_rsqmin, cg_prec, EVEN, &final_rsq, fn[0], 0, 0. );

    /* Do dslash for zero Naik term epsilon */
    for(j=0;j<order;j++){ dslash_field( multi_x[j], multi_x[j], ODD, fn[0] ); }

    eo_fermion_force_multi( eps, &(residues[1]), multi_x, order, ff_prec, fl );
}

