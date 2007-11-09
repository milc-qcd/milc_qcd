/******* d_congrad5_fn_fewsums.c - conjugate gradient for SU3/fermions ****/
/* MIMD version 7 */
/* OBSOLETE! See d_congrad5_fn.c 5/7/07 C. DeTar */

/* TEST VERSION  4/18/03, TIMING FOR GLOBAL REDUCTIONS */
/* REDUCE NUMBER OF GLOBAL SUMS */
/* 4/18/03 D.T. */

/* Kogut-Susskind fermions -- this version for "fat plus Naik" quark
   actions.  

   This code combines d_congrad5_fn.c and d_congrad5_fn_tmp.c
   Calls dslash_fn_site or dslash_fn_field depending accordingly. */

/* Jim Hetrick, Kari Rummukainen, Doug Toussaint, Steven Gottlieb */
/* 10/02/01 C. DeTar Consolidated with tmp version */
/* 11/11/06 C. DeTar Merged with d_congrad5_fn.c using ifdef FEWSUMS */

/* This version looks at the initial vector every "niter" passes */
/* The source vector is in "src", and the initial guess and answer
   in "dest".  "resid" is the residual vector, and "cg_p" and "ttt" are
   working vectors for the conjugate gradient.
   niter = maximum number of iterations.
   nrestart = number of restarts
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
#include "../include/prefetch.h"
#define FETCH_UP 1
/**#define CG_DEBUG **/


void cleanup_gathers(msg_tag *t1[16],msg_tag *t2[16]); /* dslash_fn_field.c */

#define LOOPEND
#include "../include/loopend.h"

su3_vector *ttt,*cg_p;
su3_vector *resid;
su3_vector *t_dest;
static int first_congrad = 1;

#ifdef CGTIME
static const char *milc_prec[2] = {"F", "D"};
#endif

/* prec argument is ignored */
int ks_congrad( field_offset src, field_offset dest, Real mass,
		int niter, int nrestart, Real rsqmin, int prec,
		int parity, Real *final_rsq_ptr,
		ferm_links_t *fn){
  register int i;
  register site *s;
  int iteration;	/* counter for iterations */
  Real a,b;	/* Sugar's a,b */
#ifdef FEWSUMS
  double true_rsq; /* true_rsq = rsq from actual summation of resid */
  double c_tr,c_tt,tempsum[4];	/* Re<resid|ttt>, <ttt|ttt> */
#endif
  double rsq,oldrsq,pkp;	/* resid**2,last resid*2,pkp = cg_p.K.cg_p */
  Real msq_x4;	/* 4*mass*mass */
  double source_norm;	/* squared magnitude of source vector */
  double rsqstop;	/* stopping residual normalized by source norm */
  int l_parity = 0;	/* parity we are currently doing */
  int l_otherparity = 0; /* the other parity */
  msg_tag * tags1[16], *tags2[16];	/* tags for gathers to parity and opposite */
  int special_started;	/* 1 if dslash_fn_field_special has been called */

/* Timing */

double dtimec;
double reduce_time; /*TEST*/
double nflop;

 dtimec = -dclock(); 
 reduce_time = 0.0; /*TEST*/
 
 nflop = 1187;
 if(parity==EVENANDODD)nflop *=2;
	
	special_started=0;
	/* if we want both parities, we will do even first. */
	switch(parity){
	    case(EVEN): l_parity=EVEN; l_otherparity=ODD; break;
	    case(ODD):  l_parity=ODD; l_otherparity=EVEN; break;
	    case(EVENANDODD):  l_parity=EVEN; l_otherparity=ODD; break;
	}
	msq_x4 = 4.0*mass*mass;
        iteration = 0;

	/* now we can allocate temporary variables and copy then */
	/* PAD may be used to avoid cache trashing */
#define PAD 0

 	if(first_congrad) {
	  ttt = (su3_vector *) malloc((sites_on_node+PAD)*sizeof(su3_vector));
	  cg_p = (su3_vector *) malloc((sites_on_node+PAD)*sizeof(su3_vector));
	  resid = (su3_vector *) malloc((sites_on_node+PAD)*sizeof(su3_vector));
	  t_dest = (su3_vector *) malloc((sites_on_node+PAD)*sizeof(su3_vector));
	  first_congrad = 0;
 	}

 dtimec = -dclock(); 

	/* now we copy dest to temporaries */
  FORALLSITES(i,s) {
    t_dest[i] = *(su3_vector *)F_PT(s,dest);
  }

	/* initialization process */
start:
#ifdef CG_DEBUG
	node0_printf("ks_congrad: start, parity = %d\n",parity);
#endif
        /* ttt <-  (-1)*M_adjoint*M*dest
           resid,cg_p <- src + ttt
           rsq = |resid|^2
           source_norm = |src|^2
        */
	if(special_started==1) {	/* clean up gathers */
	    cleanup_gathers(tags1,tags2);
	    special_started=0;
	}
#ifdef CG_DEBUG
	if(this_node==0)if(iteration>1)printf("CONGRAD: restart rsq = %.10e\n",rsq);
#endif
        rsq = source_norm = 0.0;
	dslash_fn_field_special(t_dest, ttt,l_otherparity,tags2,1,fn);
	dslash_fn_field_special(ttt,ttt,l_parity,tags1,1,fn);
	cleanup_gathers(tags1,tags2);
	/* ttt  <- ttt - msq_x4*src	(msq = mass squared) */
	FORSOMEPARITY(i,s,l_parity){
	  if( i < loopend-FETCH_UP ){
	    prefetch_VVVV( &ttt[i+FETCH_UP], 
			   &t_dest[i+FETCH_UP],
			   (su3_vector *)F_PT(s+FETCH_UP,src),
			   &resid[i+FETCH_UP]);
	  }
	  scalar_mult_add_su3_vector( &ttt[i], &t_dest[i],
				      -msq_x4, &ttt[i] );
	    /* note that we go back to the site structure for src */
	  add_su3_vector( (su3_vector *)F_PT(s,src),
			  &ttt[i], &resid[i] );
	  /* remember ttt contains -M_adjoint*M*src */
	  cg_p[i] = resid[i];
	  /* note that we go back to the site structure for src */
	  source_norm += (double)magsq_su3vec( (su3_vector *)F_PT(s,src) );
	  rsq += (double)magsq_su3vec( &resid[i] );
	} END_LOOP
reduce_time -= dclock();
#ifdef FEWSUMS
	true_rsq = rsq; /* not yet summed over nodes */
	tempsum[0] = source_norm; tempsum[1] = rsq;
        g_vecdoublesum( tempsum, 2 );
	source_norm = tempsum[0]; rsq = tempsum[1];
#else
	g_doublesum( &source_norm );
        g_doublesum( &rsq );
#endif
reduce_time += dclock();
#ifdef CG_DEBUG
	if(this_node==0)printf("CONGRAD: start rsq = %.10e\n",rsq);
#endif
        iteration++ ;  /* iteration counts number of multiplications
                           by M_adjoint*M */
	total_iters++;
	rsqstop = rsqmin * source_norm;
#ifdef CG_DEBUG
	node0_printf("congrad: source_norm = %e\n", (double)source_norm);
#endif
        if( rsq <= rsqstop ){
    	    /* if parity==EVENANDODD, set up to do odd sites and go back */
            if(parity == EVENANDODD) {
		l_parity=ODD; l_otherparity=EVEN;
		parity=EVEN;	/* so we won't loop endlessly */
		iteration = 0;
#ifdef CG_DEUBG
		node0_printf("instant goto start\n");
#endif
		goto start;
	    }
            *final_rsq_ptr=(Real)rsq;
#ifdef CG_DEBUG
	    node0_printf("instant return\n"); fflush(stdout);
#endif
	    cleanup_dslash_temps();
	    free(ttt); free(cg_p); free(resid); free(t_dest); first_congrad = 1;
             return (iteration);
        }
#ifdef CG_DEBUG
	pkp=0.0;
	if(mynode()==0){printf("iter=%d, rsq= %e, pkp=%e\n",
	iteration,(double)rsq,(double)pkp);fflush(stdout);}
#endif

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
#ifdef FEWSUMS
        oldrsq = true_rsq;	/* not yet summed over nodes */
#else
        oldrsq = rsq;
#endif
        pkp = 0.0;
	/* sum of neighbors */

	if(special_started==0){
	    dslash_fn_field_special( cg_p, ttt, l_otherparity, tags2, 1, fn );
	    dslash_fn_field_special( ttt, ttt, l_parity, tags1, 1, fn);
	    special_started=1;
	}
	else {
	    dslash_fn_field_special( cg_p, ttt, l_otherparity, tags2, 0, fn );
	    dslash_fn_field_special( ttt, ttt, l_parity, tags1, 0, fn);
	}

	/* finish computation of M_adjoint*m*p and p*M_adjoint*m*Kp */
	/* ttt  <- ttt - msq_x4*cg_p	(msq = mass squared) */
	/* pkp  <- cg_p.(ttt - msq*cg_p) */
	pkp = 0.0;
#ifdef FEWSUMS
	c_tr=0.0; c_tt=0.0;
#endif
	FORSOMEPARITY(i,s,l_parity){
	  if( i < loopend-FETCH_UP ){
	    prefetch_VV( &ttt[i+FETCH_UP], &cg_p[i+FETCH_UP] );
	  }
	  scalar_mult_add_su3_vector( &ttt[i], &cg_p[i], -msq_x4,
				      &ttt[i] );
	  pkp += (double)su3_rdot( &cg_p[i], &ttt[i] );
#ifdef FEWSUMS
	  c_tr += (double)su3_rdot( &ttt[i], &resid[i] );
	  c_tt += (double)su3_rdot( &ttt[i], &ttt[i] );
#endif
	} END_LOOP
reduce_time -= dclock();
#ifdef FEWSUMS
	/* finally sum oldrsq over nodes, also other sums */
	tempsum[0] = pkp; tempsum[1] = c_tr; tempsum[2] = c_tt; tempsum[3] = oldrsq;
	g_vecdoublesum( tempsum, 4 );
	pkp = tempsum[0]; c_tr = tempsum[1]; c_tt = tempsum[2]; oldrsq = tempsum[3];
#else
	g_doublesum( &pkp );
#endif
reduce_time += dclock();
	iteration++;
	total_iters++;

	a = (Real) (-rsq/pkp);

	/* dest <- dest - a*cg_p */
	/* resid <- resid - a*ttt */
#ifdef FEWSUMS
	true_rsq=0.0;
#else
	rsq=0.0;
#endif
	FORSOMEPARITY(i,s,l_parity){
	  if( i < loopend-FETCH_UP ){
	    prefetch_VVVV( &t_dest[i+FETCH_UP], 
			   &cg_p[i+FETCH_UP], 
			   &resid[i+FETCH_UP], 
			   &ttt[i+FETCH_UP] );
	  }
	  scalar_mult_add_su3_vector( &t_dest[i], &cg_p[i], a, &t_dest[i] );
	  scalar_mult_add_su3_vector( &resid[i], &ttt[i], a, &resid[i]);
#ifdef FEWSUMS
	  true_rsq += (double)magsq_su3vec( &resid[i] );
#else
	  rsq += (double)magsq_su3vec( &resid[i] );
#endif
	} END_LOOP
#ifdef FEWSUMS
	/**printf("XXX:  node %d\t%e\t%e\t%e\n",this_node,oldrsq,c_tr,c_tt);**/
	rsq = oldrsq + 2.0*a*c_tr + a*a*c_tt; /*TEST - should equal true_rsq */
	/**c_tt = true_rsq;**/ /* TEMP for test */
	/**g_doublesum(&c_tt);**/ /* TEMP true value for rsq */
	/**node0_printf("RSQTEST: %e\t%e\t%e\n",rsq,c_tt,rsq-c_tt);**/
	/**if(mynode()==0){printf("iter=%d, rsq= %e, pkp=%e\n",
	   iteration,(double)rsq,(double)pkp);fflush(stdout);}**/
#else
	g_doublesum(&rsq);
#ifdef CG_DEBUG
	if(mynode()==0){printf("iter=%d, rsq= %e, pkp=%e\n",
	   iteration,(double)rsq,(double)pkp);fflush(stdout);}
#endif
#endif	

        if( rsq <= rsqstop ){
	  /* copy t_dest back to site structure */
          FORSOMEPARITY(i,s,l_parity){
                  *(su3_vector *)F_PT(s,dest) = t_dest[i];
          } END_LOOP
    	    /* if parity==EVENANDODD, set up to do odd sites and go back */
            if(parity == EVENANDODD) {
		l_parity=ODD; l_otherparity=EVEN;
		parity=EVEN;	/* so we won't loop endlessly */
		iteration = 0;
#ifdef CG_DEBUG
		node0_printf("normal goto start\n");
#endif
		goto start;
	    }
            *final_rsq_ptr=(Real)rsq;
	    if(special_started==1) {
	      cleanup_gathers(tags1,tags2);
	      special_started = 0;
	    }
	    
#ifdef CG_DEBUG
	    node0_printf("normal return\n"); fflush(stdout);
#endif
 dtimec += dclock();
#ifdef CGTIME
 if(this_node==0){
   printf("CONGRAD5: time = %e (fn_fewsums %s) masses = 1 iters = %d mflops = %e\n",
	  dtimec,milc_prec[PRECISION-1], iteration,
	  (double)(nflop*volume*iteration/(1.0e6*dtimec*numnodes())) );
   printf("TESTCONG: reduce_time = %e iters = %d time/iter = %e\n",
	  reduce_time,iteration,reduce_time/iteration );
//{ /* time stamp for NERSC performance studies */
//      time_t time_now;
//      char time_out[26];
//      time(&time_now);
//      ctime_r(&time_now,time_out);
//      printf("  Time stamp:   %s",time_out);
//}
fflush(stdout);}
#endif
	    cleanup_dslash_temps();
	    free(ttt); free(cg_p); free(resid); free(t_dest); first_congrad = 1;
             return (iteration);
        }

	b = (Real)rsq/oldrsq;
	/* cg_p  <- resid + b*cg_p */
        FORSOMEPARITY(i,s,l_parity){
           scalar_mult_add_su3_vector( &resid[i],
                                      &cg_p[i] , b , &cg_p[i]);
        } END_LOOP

    } while( iteration%niter != 0);

    if( iteration < nrestart*niter ){
#ifdef CG_DEBUG
	node0_printf("try again goto start\n");
#endif
	 goto start;
    }
    /* if we have gotten here, no convergence after several restarts: must
	copy t_dest back to site structure */
          FORSOMEPARITY(i,s,l_parity){
                  *(su3_vector *)F_PT(s,dest) = t_dest[i];
          } END_LOOP

    /* if parity==EVENANDODD, set up to do odd sites and go back */
    if(parity == EVENANDODD) {
	l_parity=ODD; l_otherparity=EVEN;
	parity=EVEN;	/* so we won't loop endlessly */
	iteration = 0;
	goto start;
    }

    *final_rsq_ptr=rsq;
    if(special_started==1){	/* clean up gathers */
      cleanup_gathers(tags1,tags2);
      special_started = 0;
    }
    node0_printf(
        "ks_congrad: CG not converged after %d iterations, res. = %e wanted %e\n",
        iteration,rsq,rsqstop);
    fflush(stdout);
    cleanup_dslash_temps();
    free(ttt); free(cg_p); free(resid); free(t_dest); first_congrad = 1;
    return(iteration);
}
