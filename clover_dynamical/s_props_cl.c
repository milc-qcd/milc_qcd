/**************** s_props_cl.c ***************************************/
/* MIMD version 7 */
/* Clover fermions */

/* Measure propagators for measuring m_quark by method of
   Iwasaki et al.

   Only works with LU preconditioned matrix !!!
*/
#ifndef LU
BOTCH: s_props needs LU
#endif

#include "cl_dyn_includes.h"

int s_props_cl( ) {
register site *st;
register int i,spin,iters;
Real *piprop,*aprop;
register Real theta;
complex cc;
wilson_vector tvec;


    iters=0;
    piprop = (Real *)malloc( nz*sizeof(Real) );
    aprop  = (Real *)malloc( nz*sizeof(Real) );

    /* spin zero color zero for now */
    spin=0;

    /* Make source directly in chi */
    FORALLSITES(i,st){
	clear_wvec( &(st->chi) );
    }

    /* point source at zero - an even site, so LU doesn't matter */
    /**
    if( this_node==node_number(0,0,0,0) )
	lattice[node_index(0,0,0,0)].chi.d[spin].c[0].real =1.0;
    **/
    /* Wall source, modulated with Matsubara phase */
    /* Even sites only for now */
    FOREVENSITES(i,st){
        if( st->z ==0){
	    theta = PI * (Real)(st->t) / (Real)nt;
	    st->chi.d[spin].c[0].real = (Real)cos( (double)theta );
	}
    }

    /* take M inverse, result in psi */
    
    /* Load inversion control structure */
    qic.start_flag = 0;   /* Use zero initial guess for psi */
    qic.prec = PRECISION;

    /* Load Dirac matrix parameters */
    
#ifdef BI
    iters += 
      wilson_invert_site(F_OFFSET(chi),F_OFFSET(psi)
			 bicgilu_cl_site,&qic,(void *)&dcp);
#else
    iters += 
      wilson_invert_site(F_OFFSET(chi),F_OFFSET(psi),
			 cgilu_cl_site,&qic,(void *)&dcp);
#endif


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
