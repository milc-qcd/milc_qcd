/**** d_congrad5_fn_qop_milc_P.c - conjugate gradient for SU3/fermions ****/
/* MIMD version 7 */
/* Kogut-Susskind fermions -- this version for "fat plus Naik" quark
 * actions.  

 * This version implements the QOP API with the standard MILC algorithm
 * It is intended for comparing results with other QOP routines.
 * It is based only on d_congrad5_fn.c
 * At present it does not use the multicg inverter when it might.

 */

/* Jim Hetrick, Kari Rummukainen, Doug Toussaint, Steven Gottlieb */
/* 10/02/01 C. DeTar Consolidated with tmp version */
/* C.D. 12/05 Converted to the QOP API */

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

#if ( QOP_Precision == 1 )

#define CLEANUP_GATHERS_QOP_MILC cleanup_gathers_qop_milc_F
#define CLEANUP_DSLASH_QOP_MILC_TEMPS cleanup_dslash_qop_milc_temps_F
#define DSLASH_FN_QOP_MILC dslash_fn_qop_milc_F
#define DSLASH_FN_QOP_MILC_FIELD_SPECIAL dslash_fn_qop_milc_field_special_F
#define MYREAL float

#else

#define CLEANUP_GATHERS_QOP_MILC cleanup_gathers_qop_milc_D
#define CLEANUP_DSLASH_QOP_MILC_TEMPS cleanup_dslash_qop_milc_temps_D
#define DSLASH_FN_QOP_MILC dslash_fn_qop_milc_D
#define DSLASH_FN_QOP_MILC_FIELD_SPECIAL dslash_fn_qop_milc_field_special_D
#define MYREAL double

#endif

#include "generic_ks_includes.h"	/* definitions files and prototypes */
#include "../include/generic_qop.h"
#include "../include/generic_ks_qop.h"
#include "../include/qop_milc.h"
#include "../include/prefetch.h"
#define FETCH_UP 1
/**#define CG_DEBUG**/


void CLEANUP_GATHERS_QOP_MILC(msg_tag *t1[16],msg_tag *t2[16]);

#define LOOPEND
#include "../include/loopend.h"
#include <qop.h>

su3_vector *ttt,*cg_p;
su3_vector *resid;
su3_vector *t_dest;
static int first_congrad = 1;


/*  MILC imitation of Asqtad Level 3 inverter */

void QOP_asqtad_invert(QOP_info_t *info,
		       QOP_FermionLinksAsqtad *links,
		       QOP_invert_arg_t *inv_arg, 
		       QOP_resid_arg_t *res_arg, 
		       MYREAL mass,
		       QOP_ColorVector *dest_pt,
		       QOP_ColorVector *src_pt)
{
  register int i;
  register site *s;
  int iteration;	/* counter for iterations */
  MYREAL a,b;	/* Sugar's a,b */
  double rsq,oldrsq,pkp;	/* resid**2,last resid*2,pkp = cg_p.K.cg_p */
  MYREAL msq_x4;	/* 4*mass*mass */
  double source_norm;	/* squared magnitude of source vector */
  double rsqstop;	/* stopping residual normalized by source norm */
  int l_parity=0;	/* parity we are currently doing */
  int l_otherparity=0;	/* the other parity */
  msg_tag * tags1[16], *tags2[16];	/* tags for gathers to parity and opposite */
  int special_started;	/* 1 if dslash_fn_field_special has been called */
  QOP_evenodd_t qop_parity = inv_arg->evenodd;
  int parity = qop2milc_parity(qop_parity);   /* parity requested */
  MYREAL rsqmin         = res_arg->rsqmin;
  int niter           = inv_arg->max_iter;
  int max_restart     = inv_arg->restart;

  /* Convert specific QOP precision to prevailing MILC precision */
#if ( QOP_Precision == 1 )
  su3_vector *srcp = create_latvec_from_qop_milc_F(src_pt->v);
  su3_vector *solp = create_latvec_from_qop_milc_F(dest_pt->v);
  su3_matrix *fatlinks = create_links_from_qop_milc_F(links->fat->g);
  su3_matrix *longlinks = create_links_from_qop_milc_F(links->lng->g);
#else
  su3_vector *srcp = create_latvec_from_qop_milc_D(src_pt->v);
  su3_vector *solp = create_latvec_from_qop_milc_D(dest_pt->v);
  su3_matrix *fatlinks = create_links_from_qop_milc_D(links->fat->g);
  su3_matrix *longlinks = create_links_from_qop_milc_D(links->lng->g);
#endif

  MYREAL final_flop;
#ifdef CGTIME
  char *qop_prec[2] = { "F", "D" };
#endif
  
  /* Timing */
  
  double dtimec = -dclock();;
  double nflop = 1187;
  
#ifndef OLD_QOPQDP_NORM
  Real xnorm, xnorminv;
#endif

#ifndef OLD_QOPQDP_NORM
  /* qopqdp-0.9.0 changed the convention for the inverter.  For this
     version we simply change the normalization for the cases EVEN and
     ODD */
  if(parity==EVENANDODD){
    printf("QOP_asqtad_invert (MILC): EVENANDODD not supported for now\n");
    terminate(1);
  }

#define NORMFACT(a) 4.*(a)

  /* Change the normalization of the proposed solution vector */
  xnorm = NORMFACT(mass);
  xnorminv = 1./xnorm;
  FORSOMEPARITY(i,s,parity){
    scalar_mult_su3_vector(solp+i, xnorminv, solp+i);
  } END_LOOP
#endif

  /* Parity consistency is required */
  if(src_pt->evenodd != dest_pt->evenodd ||
     links->evenodd != QOP_EVENODD       ||
     qop_parity != src_pt->evenodd )
    {
      printf("QOP_asqtad_invert: Bad parity src %d dest %d links %d request %d\n",
	     src_pt->evenodd,dest_pt->evenodd,links->evenodd,qop_parity);
      terminate(1);
    }

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
	  t_dest = solp;
	  first_congrad = 0;
 	}

 dtimec = -dclock(); 

	/* initialization process */
start:
#ifdef CG_DEBUG
	node0_printf("QOP_asqtad_invert: start, parity = %d\n",parity);
#endif
        /* ttt <-  (-1)*M_adjoint*M*dest
           resid,cg_p <- src + ttt
           rsq = |resid|^2
           source_norm = |src|^2
        */
	if(special_started==1) {	/* clean up gathers */
	    CLEANUP_GATHERS_QOP_MILC(tags1,tags2);
	    special_started=0;
	}
#ifdef CG_DEBUG
	if(this_node==0)if(iteration>1)printf("CONGRAD: restart rsq = %.10e\n",rsq);
#endif
        rsq = source_norm = 0.0;
	DSLASH_FN_QOP_MILC_FIELD_SPECIAL(fatlinks, longlinks, 
					 t_dest, ttt,l_otherparity,tags2,1);
	DSLASH_FN_QOP_MILC_FIELD_SPECIAL(fatlinks, longlinks,
					 ttt,ttt,l_parity,tags1,1);
	CLEANUP_GATHERS_QOP_MILC(tags1,tags2);
	/* ttt  <- ttt - msq_x4*src	(msq = mass squared) */
	FORSOMEPARITY(i,s,l_parity){
	  if( i < loopend-FETCH_UP ){
	    prefetch_VVVV( ttt+i+FETCH_UP, 
			   t_dest+i+FETCH_UP,
			   srcp+i+FETCH_UP,
			   resid+i+FETCH_UP);
	  }
	  scalar_mult_add_su3_vector( &ttt[i], &t_dest[i],
				      -msq_x4, &ttt[i] );
	    /* note that we go back to the site structure for src */
	  add_su3_vector( srcp+i, ttt+i, resid+i );
	  /* remember ttt contains -M_adjoint*M*src */
	  cg_p[i] = resid[i];
	  /* note that we go back to the site structure for src */
	  source_norm += (double)magsq_su3vec( srcp + i );
	  rsq += (double)magsq_su3vec( resid + i );
	} END_LOOP
	g_doublesum( &source_norm );
        g_doublesum( &rsq );
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

#ifdef CG_DEBUG
	    node0_printf("instant return\n"); fflush(stdout);
#endif
	    CLEANUP_DSLASH_QOP_MILC_TEMPS();
	    free(ttt); free(cg_p); free(resid); first_congrad = 1;

	    /* Save diagnostics */
            res_arg->final_rsq=(MYREAL)rsq/source_norm;
	    res_arg->final_iter = iteration;
	    res_arg->final_restart = iteration/niter;
	    final_flop = (double)(nflop*volume*iteration)/(double)numnodes();
	    info->final_flop += final_flop;
	    dtimec += dclock();
	    info->final_sec  += dtimec;
#ifdef CGTIME
	    node0_printf("CONGRAD5: time = %e (fn_qop_milc) masses = 1 iters = %d mflops = %e\n",
			 dtimec,iteration,final_flop/(1.0e6*dtimec) );
	    fflush(stdout);
#endif
    /* Copy the solution back and free memory */

#ifndef OLD_QOPQDP_NORM
	    /* Change the normalization of the solution */
	    FORSOMEPARITY(i,s,parity){
	      scalar_mult_su3_vector(solp+i, xnorm, solp+i);
	    } END_LOOP
#endif

#if ( QOP_Precision == 1 )
	    copy_latvec_to_qop_milc_F(dest_pt->v, solp);
	    destroy_latvec_from_qop_milc_F(srcp);
	    destroy_latvec_from_qop_milc_F(solp);
	    destroy_links_from_qop_milc_F(fatlinks);
	    destroy_links_from_qop_milc_F(longlinks);
#else
	    copy_latvec_to_qop_milc_D(dest_pt->v, solp);
	    destroy_latvec_from_qop_milc_D(srcp);
	    destroy_latvec_from_qop_milc_D(solp);
	    destroy_links_from_qop_milc_D(fatlinks);
	    destroy_links_from_qop_milc_D(longlinks);
#endif
	    info->status = QOP_SUCCESS;
	    return;
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
        oldrsq = rsq;
        pkp = 0.0;
	/* sum of neighbors */

	if(special_started==0){
	    DSLASH_FN_QOP_MILC_FIELD_SPECIAL( fatlinks, longlinks,
				      cg_p, ttt, l_otherparity, tags2, 1 );
	    DSLASH_FN_QOP_MILC_FIELD_SPECIAL( fatlinks, longlinks,
				      ttt, ttt, l_parity, tags1, 1);
	    special_started=1;
	}
	else {
	    DSLASH_FN_QOP_MILC_FIELD_SPECIAL( fatlinks, longlinks,
				      cg_p, ttt, l_otherparity, tags2, 0 );
	    DSLASH_FN_QOP_MILC_FIELD_SPECIAL( fatlinks, longlinks,
				      ttt, ttt, l_parity, tags1, 0);
	}

	/* finish computation of M_adjoint*m*p and p*M_adjoint*m*Kp */
	/* ttt  <- ttt - msq_x4*cg_p	(msq = mass squared) */
	/* pkp  <- cg_p.(ttt - msq*cg_p) */
	pkp = 0.0;
	FORSOMEPARITY(i,s,l_parity){
	  if( i < loopend-FETCH_UP ){
	    prefetch_VV( &ttt[i+FETCH_UP], &cg_p[i+FETCH_UP] );
	  }
	  scalar_mult_add_su3_vector( &ttt[i], &cg_p[i], -msq_x4,
				      &ttt[i] );
	  pkp += (double)su3_rdot( &cg_p[i], &ttt[i] );
	} END_LOOP
	g_doublesum( &pkp );
	iteration++;
	total_iters++;

	a = (MYREAL) (-rsq/pkp);

	/* dest <- dest - a*cg_p */
	/* resid <- resid - a*ttt */
	rsq=0.0;
	FORSOMEPARITY(i,s,l_parity){
	  if( i < loopend-FETCH_UP ){
	    prefetch_VVVV( &t_dest[i+FETCH_UP], 
			   &cg_p[i+FETCH_UP], 
			   &resid[i+FETCH_UP], 
			   &ttt[i+FETCH_UP] );
	  }
	  scalar_mult_add_su3_vector( &t_dest[i], &cg_p[i], a, &t_dest[i] );
	  scalar_mult_add_su3_vector( &resid[i], &ttt[i], a, &resid[i]);
	  rsq += (double)magsq_su3vec( &resid[i] );
	} END_LOOP
	g_doublesum(&rsq);
#ifdef CG_DEBUG
	if(mynode()==0){printf("iter=%d, rsq= %e, pkp=%e\n",
	   iteration,(double)rsq,(double)pkp);fflush(stdout);}
#endif
	
        if( rsq <= rsqstop ){
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
	    if(special_started==1) {
	      CLEANUP_GATHERS_QOP_MILC(tags1,tags2);
	      special_started = 0;
	    }
	    
#ifdef CG_DEBUG
	    node0_printf("normal return\n"); fflush(stdout);
#endif
	    CLEANUP_DSLASH_QOP_MILC_TEMPS();
	    free(ttt); free(cg_p); free(resid); first_congrad = 1;

	    /* Save diagnostics */
            res_arg->final_rsq  = (MYREAL)rsq/source_norm;
	    res_arg->final_iter = iteration;
	    res_arg->final_restart = iteration/niter;
	    final_flop = (double)(nflop*volume*iteration)/(double)numnodes();
	    info->final_flop += final_flop;
	    dtimec += dclock();
	    info->final_sec  += dtimec;
#ifdef CGTIME
	    node0_printf("CONGRAD5: time = %e (fn_qop_milc) masses = 1 iters = %d mflops = %e\n",
			 dtimec,iteration,final_flop/(1.0e6*dtimec) );
	    fflush(stdout);
#endif

#ifndef OLD_QOPQDP_NORM
	    /* Change the normalization of the solution */
	    FORSOMEPARITY(i,s,parity){
	      scalar_mult_su3_vector(solp+i, xnorm, solp+i);
	    } END_LOOP
#endif
	    /* Copy the solution and free memory */
#if ( QOP_Precision == 1 )
	    copy_latvec_to_qop_milc_F(dest_pt->v, solp);
	    destroy_latvec_from_qop_milc_F(srcp);
	    destroy_latvec_from_qop_milc_F(solp);
	    destroy_links_from_qop_milc_F(fatlinks);
	    destroy_links_from_qop_milc_F(longlinks);
#else
	    copy_latvec_to_qop_milc_D(dest_pt->v, solp);
	    destroy_latvec_from_qop_milc_D(srcp);
	    destroy_latvec_from_qop_milc_D(solp);
	    destroy_links_from_qop_milc_D(fatlinks);
	    destroy_links_from_qop_milc_D(longlinks);
#endif
	    info->status = QOP_SUCCESS;
	    return;
        }

	b = (MYREAL)rsq/oldrsq;
	/* cg_p  <- resid + b*cg_p */
        FORSOMEPARITY(i,s,l_parity){
           scalar_mult_add_su3_vector( &resid[i],
                                      &cg_p[i] , b , &cg_p[i]);
        } END_LOOP

    } while( iteration%niter != 0);

    if( iteration < max_restart*niter ){
#ifdef CG_DEBUG
	node0_printf("try again goto start\n");
#endif
	 goto start;
    }
    /* if we have gotten here, no convergence after several restarts: */

    /* if parity==EVENANDODD, set up to do odd sites and go back */
    if(parity == EVENANDODD) {
	l_parity=ODD; l_otherparity=EVEN;
	parity=EVEN;	/* so we won't loop endlessly */
	iteration = 0;
	goto start;
    }

    if(special_started==1){	/* clean up gathers */
      CLEANUP_GATHERS_QOP_MILC(tags1,tags2);
      special_started = 0;
    }

    node0_printf(
        "QOP_asqtad_invert: CG not converged after %d iterations, res. = %e wanted %e\n",
        iteration,rsq,rsqstop);
    fflush(stdout);
    CLEANUP_DSLASH_QOP_MILC_TEMPS();
    free(ttt); free(cg_p); free(resid); first_congrad = 1;

#ifndef OLD_QOPQDP_NORM
    /* Change the normalization of the solution */
    FORSOMEPARITY(i,s,parity){
      scalar_mult_su3_vector(solp+i, xnorm, solp+i);
    } END_LOOP
#endif
    /* Copy the solution and free memory */
#if ( QOP_Precision == 1 )
    copy_latvec_to_qop_milc_F(dest_pt->v, solp);
    destroy_latvec_from_qop_milc_F(srcp);
    destroy_latvec_from_qop_milc_F(solp);
    destroy_links_from_qop_milc_F(fatlinks);
    destroy_links_from_qop_milc_F(longlinks);
#else
    copy_latvec_to_qop_milc_D(dest_pt->v, solp);
    destroy_latvec_from_qop_milc_D(srcp);
    destroy_latvec_from_qop_milc_D(solp);
    destroy_links_from_qop_milc_D(fatlinks);
    destroy_links_from_qop_milc_D(longlinks);
#endif

    /* Save diagnostics */

    res_arg->final_rsq  =(MYREAL)rsq/source_norm;
    res_arg->final_iter = iteration;
    res_arg->final_restart = iteration/niter;
    final_flop = (double)(nflop*volume*iteration)/(double)numnodes();
    info->final_flop += final_flop;
    dtimec += dclock();
    info->final_sec  += dtimec;
#ifdef CGTIME
    node0_printf("CONGRAD5: time = %e (qop_milc %s) masses = 1 iters = %d mflops = %e\n",
		 dtimec,qop_prec[QOP_Precision-1],
		 iteration,final_flop/(1.0e6*dtimec) );
    fflush(stdout);
#endif
    info->status = QOP_FAIL;
}

