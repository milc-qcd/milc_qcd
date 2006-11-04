/******* d_congrad5_eo.c - conjugate gradient for SU3/fermions ****/
/* MIMD version 7 */
/* Kogut-Susskind fermions  -- improved */
/* for general even-odd dslash, using dslash_eo(), no longlinks
   or fatlinks */

/* Jim Hetrick, Kari Rummukainen, Doug Toussaint */

/* This version looks at the initial vector every "niter" passes */
/* The source vector is in "src", and the initial guess and answer
   in "dest".  "resid" is the residual vector, and "cg_p" and "ttt" are
   working vectors for the conjugate gradient.
   niter = maximum number of iterations.
   rsqmin = desired rsq, quit when we reach rsq <= rsqmin*source_norm.
	This is different than our old definition of the stopping
	criterion.  To convert an old stopping residual to the new
	one, multiply the old one by sqrt( (2/3)/(8+2*m) )
        This is because the source is obtained from
        a random vector with average squared magnitude 3 on each site.
        Then, on 1/2 the sites, we gather and sum the eight neighboring
        random vectors and add 2*m times the local vector.
            source = M_adjoint*R, on even sites
   reinitialize after niters iterations and try once more.
   parity=EVEN = do only even sites, parity=ODD = do odd sites,
   parity=EVENANDODD = do all sites
*/
#include "generic_ks_includes.h"	/* definitions files and prototypes */
/* debug */
/*#define CGTIME*/  /* Commented out - Defined at compilation time */

void cleanup_gathers(msg_tag *t1[16],msg_tag *t2[16]);

#include "../include/prefetch.h"
#define FETCH_UP 1
#define LOOPEND
#include "../include/loopend.h"

int ks_congrad( field_offset src, field_offset dest, Real mass,
    int niter, int nrestart, Real rsqmin, int parity, Real *final_rsq_ptr ){
  register int i;
  register site *s;
  int iteration;	/* counter for iterations */
  Real a,b;	/* Sugar's a,b */
  double rsq,oldrsq,pkp;	/* resid**2,last resid*2,pkp = cg_p.K.cg_p */
  Real msq_x4;	/* 4*mass*mass */
  double source_norm;	/* squared magnitude of source vector */
  double rsqstop;	/* stopping residual normalized by source norm */
  int l_parity=0;	/* parity we are currently doing */
  int l_otherparity=0;	/* the other parity */
  
#ifdef CGTIME
double dtimec;
#endif
double nflop;

nflop = 1187;
if(parity==EVENANDODD)nflop *=2;
	
	/* if we want both parities, we will do even first. */
	switch(parity){
	    case(EVEN): l_parity=EVEN; l_otherparity=ODD; break;
	    case(ODD):  l_parity=ODD; l_otherparity=EVEN; break;
	    case(EVENANDODD):  l_parity=EVEN; l_otherparity=ODD; break;
	}
	msq_x4 = 4.0*mass*mass;
        iteration = 0;

#ifdef CGTIME
 dtimec = -dclock(); 
#endif

	/* initialization process */
start:
	/**node0_printf("ks_congrad4: start, parity = %d\n",parity);**/
        /* ttt <-  (-1)*M_adjoint*M*dest
           resid,cg_p <- src + ttt
           rsq = |resid|^2
           source_norm = |src|^2
        */
/**if(this_node==0)if(iteration>1)printf("CONGRAD: restart rsq = %.10e\n",rsq);**/
        rsq = source_norm = 0.0;
	dslash_eo_site( dest, F_OFFSET(ttt), l_otherparity);
	dslash_eo_site(F_OFFSET(ttt),F_OFFSET(ttt),l_parity);
	/* ttt  <- ttt - msq_x4*src	(msq = mass squared) */
	FORSOMEPARITY(i,s,l_parity){
	  if( i < loopend-FETCH_UP ){
	    prefetch_VVVV( &((s+FETCH_UP)->ttt), 
			   (su3_vector *)F_PT(s+FETCH_UP,dest),
			   (su3_vector *)F_PT(s+FETCH_UP,src),
			   &((s+FETCH_UP)->resid));
	  }
	    scalar_mult_add_su3_vector( &(s->ttt), (su3_vector *)F_PT(s,dest),
		-msq_x4, &(s->ttt) );
	    add_su3_vector( (su3_vector *)F_PT(s,src), &(s->ttt), &(s->resid) );
		/* remember ttt contains -M_adjoint*M*src */
	    s->cg_p = s->resid;
	    source_norm += (double) magsq_su3vec( (su3_vector *)F_PT(s,src) );
            rsq += (double) magsq_su3vec( &(s->resid) );
	} END_LOOP
	g_doublesum( &source_norm );
        g_doublesum( &rsq );
/**if(this_node==0)printf("CONGRAD: start rsq = %.10e\n",rsq);**/
        iteration++ ;  /* iteration counts number of multiplications
                           by M_adjoint*M */
	total_iters++;
	rsqstop = rsqmin * source_norm;
	/**node0_printf("congrad: source_norm = %e\n", (double)source_norm);**/
        if( rsq <= rsqstop ){
    	    /* if parity==EVENANDODD, set up to do odd sites and go back */
            if(parity == EVENANDODD) {
		l_parity=ODD; l_otherparity=EVEN;
		parity=EVEN;	/* so we won't loop endlessly */
		iteration = 0;
		/**node0_printf("instant goto start\n"); **/
		goto start;
	    }
            *final_rsq_ptr=(Real)rsq;
	    /**node0_printf("instant return\n"); fflush(stdout);**/
             return (iteration);
        }
	/**pkp=0.0;
	if(mynode()==0){printf("iter=%d, rsq= %e, pkp=%e\n",
	iteration,(double)rsq,(double)pkp);fflush(stdout);}**/

    /* main loop - do until convergence or time to restart */
        /*
           oldrsq <- rsq
           ttt <- (-1)*M_adjoint*M*cg_p
           pkp <- (-1)*cg_p.M_adjoint*M.cg_p
           a <- -rsq/pkp
           dest <- dest + a*cg_p
           resid <- resid + a*ttt
           rsq <- |resid|^2
           b <- rsq/oldrsq
           cg_p <- resid + b*cg_p
        */
    do{
        oldrsq = rsq;
        pkp = 0.0;
	/* sum of neighbors */

	dslash_eo_site( F_OFFSET(cg_p), F_OFFSET(ttt), l_otherparity);
	dslash_eo_site( F_OFFSET(ttt), F_OFFSET(ttt), l_parity);

	/* finish computation of M_adjoint*m*p and p*M_adjoint*m*Kp */
	/* ttt  <- ttt - msq_x4*cg_p	(msq = mass squared) */
	/* pkp  <- cg_p.(ttt - msq*cg_p) */
	pkp = 0.0;
	FORSOMEPARITY(i,s,l_parity){
	  if( i < loopend-FETCH_UP ){
	    prefetch_VV( &((s+FETCH_UP)->ttt), &((s+FETCH_UP)->cg_p) );
	  }
	  scalar_mult_add_su3_vector( &(s->ttt), &(s->cg_p), -msq_x4,
				      &(s->ttt) );
	  pkp += (double)su3_rdot( &(s->cg_p), &(s->ttt) );
	} END_LOOP
	g_doublesum( &pkp );
	iteration++;
	total_iters++;

	a = (Real) (-rsq/pkp);

	/* dest <- dest - a*cg_p */
	/* resid <- resid - a*ttt */
	rsq=0.0;
	FORSOMEPARITY(i,s,l_parity){
	  if( i < loopend-FETCH_UP ){
	    prefetch_VVVV( 
			  (su3_vector *)F_PT((s+FETCH_UP),dest), 
			  &((s+FETCH_UP)->cg_p),
			  &((s+FETCH_UP)->resid), 
			  &((s+FETCH_UP)->ttt) );
	  }
	  scalar_mult_add_su3_vector( (su3_vector *)F_PT(s,dest), &(s->cg_p),
				      a, (su3_vector *)F_PT(s,dest) );
	  scalar_mult_add_su3_vector( &(s->resid), &(s->ttt), a, &(s->resid));
	  rsq += (double)magsq_su3vec( &(s->resid) );
	} END_LOOP
	g_doublesum(&rsq);
	/**if(mynode()==0){printf("iter=%d, rsq= %e, pkp=%e\n",
	   iteration,(double)rsq,(double)pkp);fflush(stdout);}**/

        if( rsq <= rsqstop ){
    	    /* if parity==EVENANDODD, set up to do odd sites and go back */
            if(parity == EVENANDODD) {
		l_parity=ODD; l_otherparity=EVEN;
		parity=EVEN;	/* so we won't loop endlessly */
		iteration = 0;
		/**node0_printf("normal goto start\n"); **/
		goto start;
	    }
            *final_rsq_ptr=(Real)rsq;
	    
	    /**node0_printf("normal return\n"); fflush(stdout);**/
#ifdef CGTIME
 dtimec += dclock();
if(this_node==0){printf("CONGRAD5: time = %e iters = %d mflops = %e\n",
dtimec,iteration,(double)(nflop*volume*iteration/(1.0e6*dtimec*numnodes())) );
fflush(stdout);}
#endif
             return (iteration);
        }

	b = (Real)rsq/oldrsq;
	/* cg_p  <- resid + b*cg_p */
	scalar_mult_add_latvec( F_OFFSET(resid), F_OFFSET(cg_p),
	    b, F_OFFSET(cg_p), l_parity);

    } while( iteration%niter != 0);

    if( iteration < nrestart*niter ){
	/**node0_printf("try again goto start\n");**/
	 goto start;
    }

    /* if parity==EVENANDODD, set up to do odd sites and go back */
    if(parity == EVENANDODD) {
	l_parity=ODD; l_otherparity=EVEN;
	parity=EVEN;	/* so we won't loop endlessly */
	iteration = 0;
	goto start;
    }

    *final_rsq_ptr=rsq;
    node0_printf(
        "ks_congrad: CG not converged after %d iterations, res. = %e wanted %e\n",
        iteration,rsq,rsqstop);
    fflush(stdout);
    return(iteration);
}
