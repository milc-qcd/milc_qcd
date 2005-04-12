/********** update.c ****************************************************/
/* MIMD version 6 */

/*
 Update lattice.
 Improved method for 1-4 flavors:
	update U by (epsilon/2)*(1-Nf/4)
	compute PHI
	update U to epsilon/2
	compute X
	update H, full step
	update U to next time needed

 This routine does not refresh the antihermitian momenta.
 This routine begins at "integral" time, with H and U evaluated
 at same time.
*/
#include "ks_imp_includes.h"	/* definitions files and prototypes */

int update()  {
int step, iters=0;
Real final_rsq;
#ifdef HMC_ALGORITHM
double startaction,endaction,d_action();
Real xrandom;
#endif

    /* refresh the momenta */
    ranmom();

    /* do "steps" microcanonical steps"  */
    for(step=1; step <= steps; step++){
 
#ifdef PHI_ALGORITHM
        /* generate a pseudofermion configuration only at start*/
	/* also clear xxx, since zero is our best guess for the solution
	   with a new random phi field. */
     	if(step==1){
	    clear_latvec( F_OFFSET(xxx1), EVENANDODD );
	    grsource_imp( F_OFFSET(phi1), mass1, EVEN);
	    clear_latvec( F_OFFSET(xxx2), EVENANDODD );
	    grsource_imp( F_OFFSET(phi2), mass2, EVEN);
	}

#ifdef HMC_ALGORITHM
        /* find action */
        /* do conjugate gradient to get (Madj M)inverse * phi */
        if(step==1){
            /* do conjugate gradient to get (Madj M)inverse * phi */
	    iters += ks_congrad( F_OFFSET(phi1), F_OFFSET(xxx1), mass1,
		 niter, rsqmin, EVEN, &final_rsq );
	    iters += ks_congrad( F_OFFSET(phi2), F_OFFSET(xxx2), mass2,
		niter, rsqmin, EVEN, &final_rsq );

     	    startaction=d_action();
            /* copy link field to old_link */
	    gauge_field_copy( F_OFFSET(link[0]), F_OFFSET(old_link[0]));
        }
#endif

	/* update U's to middle of interval */
     	update_u(0.5*epsilon);

#else /* "R" algorithm */
       	/* first update the U's to special time interval */
        /* and generate a pseudofermion configuration */
	/* probably makes most sense if nflavors1 >= nflavors2 */

       	update_u(epsilon*(0.5-nflavors1/8.0));
	clear_latvec( F_OFFSET(xxx1), EVENANDODD );
     	grsource_imp( F_OFFSET(phi1), mass1, EVEN);

       	update_u(epsilon*((nflavors1-nflavors2)/8.0));
	clear_latvec( F_OFFSET(xxx2), EVENANDODD );
     	grsource_imp( F_OFFSET(phi2), mass2, EVEN);

	/* update U's to middle of interval */
     	update_u(epsilon*nflavors2/8.0);
#endif

        /* do conjugate gradient to get (Madj M)inverse * phi */
#if 1
     	iters += ks_congrad( F_OFFSET(phi1), F_OFFSET(xxx1), mass1,
	    niter, rsqmin, EVEN, &final_rsq );
     	iters += ks_congrad( F_OFFSET(phi2), F_OFFSET(xxx2), mass2,
	    niter, rsqmin, EVEN, &final_rsq );
#else
	/* For future use */
	iters += ks_congrad_two_src( F_OFFSET(phi1), F_OFFSET(phi2),
				     F_OFFSET(xxx1), F_OFFSET(xxx2),
				     mass1, mass2, niter, rsqmin, 
				     EVEN, &final_rsq);
#endif
	/* now update H by full time interval */
    	update_h(epsilon);

    	/* update U's by half time step to get to even time */
    	update_u(epsilon*0.5);

        /* reunitarize the gauge field */
	rephase( OFF );
        reunitarize();
	rephase( ON );

    }	/* end loop over microcanonical steps */

#ifdef HMC_ALGORITHM
    /* find action */
    /* do conjugate gradient to get (Madj M)inverse * phi */
    iters += ks_congrad( F_OFFSET(phi1), F_OFFSET(xxx1), mass1,
	 niter, rsqmin, EVEN, &final_rsq );
    iters += ks_congrad( F_OFFSET(phi2), F_OFFSET(xxx2), mass2,
	 niter, rsqmin, EVEN, &final_rsq );
    endaction=d_action();
    /* decide whether to accept, if not, copy old link field back */
    /* careful - must generate only one random number for whole lattice */
    if(this_node==0)xrandom = myrand(&node_prn);
    broadcast_float(&xrandom);
    if( exp( (double)(startaction-endaction) ) < xrandom ){
	if(steps > 0)
	    gauge_field_copy( F_OFFSET(old_link[0]), F_OFFSET(link[0]) );
#ifdef FN
	valid_longlinks=valid_fatlinks=0;
#endif
	node0_printf("REJECT: delta S = %e\n", (double)(endaction-startaction));
    }
    else {
	node0_printf("ACCEPT: delta S = %e\n", (double)(endaction-startaction));
    }
#endif

    if(steps > 0)return (iters/steps);
    else return(-99);
}

