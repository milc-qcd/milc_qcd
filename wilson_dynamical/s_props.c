/**************** s_props.c ***************************************/
/* MIMD version 6 */
/* Wilson fermions */

/* Measure propagators for measuring m_quark by method of
   Iwasaki et al.

   Only works with LU preconditioned matrix !!!
*/
#ifndef LU
BOTCH: s_props needs LU
#endif

#include "wi_dyn_includes.h"

int s_props( ) {
register site *st;
register int i,spin,iters;
Real final_rsq,*piprop,*aprop;
register Real theta;
complex cc;
wilson_vector tvec;

    iters=0;
    piprop = (Real *)malloc( nz*sizeof(Real) );
    aprop  = (Real *)malloc( nz*sizeof(Real) );
    /* spin zero color zero for now */
    spin=0;

    FORALLSITES(i,st){
	clear_wvec( &(st->g_rand) );
    }

    /* point source at zero - an even site, so LU doesn't matter */
    /**
    if( this_node==node_number(0,0,0,0) )
	lattice[node_index(0,0,0,0)].g_rand.d[spin].c[0].real =1.0;
    **/
    /* Wall source, modulated with Matsubara phase */
    /* Even sites only for now */
    FOREVENSITES(i,st){
        if( st->z ==0){
	    theta =  PI*(Real)( st->t)/(Real)nt;
	    st->g_rand.d[spin].c[0].real = (Real)cos( (double)theta );
	}
    }

    /* Multiply by M_adjoint, result in chi */
    dslash_w_site( F_OFFSET(g_rand), F_OFFSET(g_rand), MINUS, ODD);
    dslash_w_site( F_OFFSET(g_rand), F_OFFSET(chi)   , MINUS, EVEN);
    FOREVENSITES(i,st){
        scalar_mult_add_wvec( &(st->g_rand), &(st->chi), -kappa*kappa,
            &(st->chi) );
    }

    /* Conjugate gradient to get M_inverse * source, result in psi */
    iters += congrad_w( niter,rsqprop,&final_rsq);

    /* Multiply by U_inverse to get solution on all sites */
    /* Since the source was zero on odd sites there is no diagonal term */
    dslash_w_site( F_OFFSET(psi), F_OFFSET(psi), PLUS, ODD);
    FORODDSITES(i,st){scalar_mult_wvec( &(st->psi), kappa, &(st->psi) );}
    /* Wingate found bug 10/4/95 */

    /* Sum pion propagator and gamma_3 propagator over z slices  */
    for(i=0;i<nz;i++) piprop[i] = aprop[i] = 0.0;
    FORALLSITES(i,st){
	piprop[st->z] += magsq_wvec( &(st->psi) );
	mult_by_gamma( &(st->psi), &tvec, ZUP );
	cc = wvec_dot( &(st->psi), &tvec );
	aprop[st->z] += cc.real;
    }
    for(i=0;i<nz;i++){
	g_floatsum( &(piprop[i]) );
	g_floatsum( &(aprop[i]) );
    }

    /* Dump results */
    if(this_node==0)for(i=0;i<nz;i++)
	printf("SPROPS %d %d %.7e %.7e\n",spin,i,(double)piprop[i],
	    (double)aprop[i]);

    free(piprop); free(aprop);
    return(iters);
}
