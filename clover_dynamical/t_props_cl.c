/**************** t_props_cl.c ***************************************/
/* MIMD version 6 */
/* Clover fermions */

/* Measure propagators for measuring m_quark by method of
   Iwasaki et al.

   Only works with LU preconditioned matrix !!!
*/
#ifndef LU
BOTCH: t_props needs LU
#endif

#include "cl_dyn_includes.h"

int t_props_cl( ) {
register site *st;
register int i,spin,iters;
Real *piprop,*aprop;
complex cc;
wilson_vector tvec;


    iters=0;
    piprop = (Real *)malloc( nt*sizeof(Real) );
    aprop  = (Real *)malloc( nt*sizeof(Real) );

    /* spin zero color zero for now */
    spin=0;

    /* Make source directly in chi */
    FORALLSITES(i,st){
	clear_wvec( &(st->chi) );
    }

    /* Wall source */
    /* Even sites only for now */
    FOREVENSITES(i,st){
        if( st->t ==0) st->chi.d[spin].c[0].real = 1.0;
    }

    /* take M inverse, result in psi */
    
    /* Load inversion control structure */
    qic.start_flag = 0;   /* Use zero initial guess for psi */
    
    /* Load Dirac matrix parameters */
    dcp.work_f_mn = F_OFFSET(tempmat1);
    
#ifdef BI
    /* Load temporaries specific to inverter */
    qic.wv1 = F_OFFSET(tmp);
    qic.wv2 = F_OFFSET(mp);
    qic.wv3 = F_OFFSET(p);
    qic.wv4 = F_OFFSET(sss);
    
    iters += 
      wilson_invert(F_OFFSET(chi),F_OFFSET(psi),F_OFFSET(r),
		    bicgilu_cl,&qic,(void *)&dcp);
#else
    /* Load temporaries specific to inverter */
    qic.wv1 = F_OFFSET(tmp);
    qic.wv2 = F_OFFSET(mp);
    
    iters += 
      wilson_invert(F_OFFSET(chi),F_OFFSET(psi),F_OFFSET(r),
		    cgilu_cl,&qic,(void *)&dcp);
#endif

    /* Sum pion propagator and gamma_3 propagator over z slices  */
    for(i=0;i<nt;i++) piprop[i] = aprop[i] = 0.0;
    FORALLSITES(i,st){
	piprop[st->t] += magsq_wvec( &(st->psi) );
	mult_by_gamma( &(st->psi), &tvec, TUP );
	cc = wvec_dot( &(st->psi), &tvec );
	aprop[st->t] += cc.real;
    }
    for(i=0;i<nt;i++){
	g_floatsum( &(piprop[i]) );
	g_floatsum( &(aprop[i]) );
    }

    /* Dump results */
    if(this_node==0)for(i=0;i<nt;i++)
	printf("TPROPS %d %d %.7e %.7e\n",spin,i,(double)piprop[i],
	    (double)aprop[i]);

    free(piprop); free(aprop);
    free_clov();
    return(iters);
}
