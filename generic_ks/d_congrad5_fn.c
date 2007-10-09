/******* d_congrad5_fn.c - conjugate gradient for SU3/fermions ****/
/* MIMD version 7 */
/* Kogut-Susskind fermions -- this version for "fat plus Naik" quark
   actions.  

   Previously called d_congrad5_fn_fewsums.c and d_congrad5_fn_tmp.c
*/

/* Jim Hetrick, Kari Rummukainen, Doug Toussaint, Steven Gottlieb */
/* 10/02/01 C. DeTar Consolidated with tmp version */
/* 4/18/03  D.T.     Consolidate global sums (fewer reductions) */
/* 5/6/07   C. DeTar True residual stopping criterion now. */

/* This version looks at the initial vector every "niter" passes */
/* The source vector is in "src", and the initial guess and answer
   in "dest".  "resid" is the residual vector, and "cg_p" and "ttt" are
   working vectors for the conjugate gradient.
   niter = maximum number of iterations before restarting.
   max_restarts = max number of restarts
   rsqmin = desired rsq, quit when we reach rsq <= rsqmin*source_norm.

   reinitialize after niters iterations and try once more.
*/
#include "generic_ks_includes.h"
#include "../include/prefetch.h"
#define FETCH_UP 1

#ifdef CGTIME
static const char *prec_label[2] = {"F", "D"};
#endif

/*#define CG_DEBUG*/

#define LOOPEND
#include "../include/loopend.h"

/* The Fermilab relative residue */

static Real 
relative_residue(su3_vector *p, su3_vector *q, int parity)
{
  double residue, num, den;
  int i;
  site *s;
  
  residue = 0;
  FORSOMEPARITY(i,s,parity){
    num = (double)magsq_su3vec( &(p[i]) );
    den = (double)magsq_su3vec( &(q[i]) );
    residue += (den==0) ? 1.0 : (num/den);
  } END_LOOP

  g_doublesum(&residue);

  if(parity == EVENANDODD)
    return sqrt(residue/volume);
  else
    return sqrt(2*residue/volume);
}

static int
ks_congrad_parity( su3_vector *t_src, su3_vector *t_dest, 
		   quark_invert_control *qic, Real mass,
		   fn_links_t *fn, ks_action_paths *ap){
  register int i;
  register site *s;
  int iteration;	/* counter for iterations */
  Real a,b;           	/* Sugar's a,b */
#ifdef FEWSUMS
  double actual_rsq;      /* rsq from actual summation of resid */
  double c_tr,c_tt,tempsum[4];	/* Re<resid|ttt>, <ttt|ttt> */
#endif
  double rsq = 0,relrsq = 0; /* resid**2, rel resid*2 */
  double oldrsq,pkp;	/*last resid*2,pkp = cg_p.K.cg_p */
  Real msq_x4;	/* 4*mass*mass */
  double source_norm;	/* squared magnitude of source vector */
  int otherparity = 0; /* the other parity */
  msg_tag * tags1[16], *tags2[16];	/* tags for gathers to parity and opposite */
  int special_started = 0; /* 1 if dslash_fn_field_special has been called */
  int nrestart;  /* Restart counter */
  su3_vector *ttt, *cg_p, *resid;
  char myname[] = "ks_congrad_parity";

  /* Unpack structure */
  int niter        = qic->max;      /* maximum number of iters per restart */
  int max_restarts = qic->nrestart; /* maximum restarts */
  Real rsqmin      = qic->resid;    /* desired residual - 
			 normalized as sqrt(r*r)/sqrt(src_e*src_e) */
  Real relrsqmin   = qic->relresid; /* desired relative residual (FNAL)*/
  int parity       = qic->parity;   /* EVEN, ODD */

  int max_cg = max_restarts*niter; /* Maximum number of iterations */

  msq_x4 = 4.0*mass*mass;

  switch(parity){
  case(EVEN): otherparity=ODD; break;
  case(ODD):  otherparity=EVEN; break;
  }

  /* Allocate temporary variables */
  /* PAD may be used to avoid cache trashing */
#define PAD 0
  ttt = (su3_vector *) malloc((sites_on_node+PAD)*sizeof(su3_vector));
  cg_p = (su3_vector *) malloc((sites_on_node+PAD)*sizeof(su3_vector));
  resid = (su3_vector *) malloc((sites_on_node+PAD)*sizeof(su3_vector));

  if(ttt == NULL || cg_p == NULL || resid == NULL){
    printf("%s(%d): No room for temporaries\n",myname,this_node);
  }

  /* Source norm */
  source_norm = 0.0;
  FORSOMEPARITY(i,s,parity){
    source_norm += (double)magsq_su3vec( &t_src[i] );
  } END_LOOP
  g_doublesum( &source_norm );
#ifdef CG_DEBUG
  node0_printf("congrad: source_norm = %e\n", (double)source_norm);
#endif

  /* Start CG iterations */
  
  nrestart = 0;
  iteration = 0;
  qic->size_r = 0;
  qic->size_relr = 0;

  while(1) {
    /* Check for completion */
    if( ( iteration % niter == 0 ) || 
	( ( rsqmin    <= 0 || rsqmin    > qic->size_r   ) &&
	  ( relrsqmin <= 0 || relrsqmin > qic->size_relr) ) ) 
      {
	
	/* (re)initialization process */
	
	/* Compute true residual and relative residual */

	/* ttt <-  (-1)*M_adjoint*M*dest
	   resid,cg_p <- src + ttt
	   rsq = |resid|^2
	   source_norm = |src|^2
	*/
	if(special_started==1) {	/* clean up gathers */
	    cleanup_gathers(tags1,tags2);
	    special_started=0;
	}
	rsq = 0.0;
	dslash_fn_field(t_dest, ttt, otherparity, fn, ap);
	dslash_fn_field(ttt, ttt, parity, fn, ap);
	/* ttt  <- ttt - msq_x4*src	(msq = mass squared) */
	FORSOMEPARITYDOMAIN(i,s,parity){
	  if( i < loopend-FETCH_UP ){
	    prefetch_VVVV( &ttt[i+FETCH_UP], 
			   &t_dest[i+FETCH_UP],
			   &t_src[i+FETCH_UP],
			   &resid[i+FETCH_UP]);
	  }
	  scalar_mult_add_su3_vector( &ttt[i], &t_dest[i], -msq_x4, &ttt[i] );
	  add_su3_vector( &t_src[i], &ttt[i], &resid[i] );
	  /* remember ttt contains -M_adjoint*M*src */
	  cg_p[i] = resid[i];
	  rsq += (double)magsq_su3vec( &resid[i] );
	} END_LOOP
#ifdef FEWSUMS
	actual_rsq = rsq; /* not yet summed over nodes */
#endif
        g_doublesum( &rsq );

	if(relrsqmin > 0)
	  relrsq = relative_residue(resid, t_dest, parity);

	qic->final_rsq    = (Real)rsq/source_norm;
	qic->final_relrsq = (Real)relrsq;


	iteration++ ;  /* iteration counts number of multiplications
			  by M_adjoint*M */
	total_iters++;

#ifdef CG_DEBUG
	if(this_node==0)printf("CONGRAD: (re)start %d rsq = %.10e relrsq %.10e\n",
			       nrestart, qic->final_rsq, qic->final_relrsq);
#endif
	/* Quit when true residual and true relative residual are within
	   tolerance or when we exhaust iterations or restarts */
	
	if( iteration >= max_cg || 
	    nrestart  >= max_restarts ||
	    ( ( rsqmin    <= 0 || rsqmin    > qic->final_rsq   ) &&
	      ( relrsqmin <= 0 || relrsqmin > qic->final_relrsq) ) ) break;
	
	nrestart++;
      }

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

#ifdef FEWSUMS
    oldrsq = actual_rsq;	/* not yet summed over nodes */
#else
    oldrsq = rsq;
#endif
    /* sum of neighbors */
    
    if(special_started==0){
      dslash_fn_field_special( cg_p, ttt, otherparity, tags2, 1, fn, ap );
      dslash_fn_field_special( ttt, ttt, parity, tags1, 1, fn, ap);
      special_started=1;
    }
    else {
      dslash_fn_field_special( cg_p, ttt, otherparity, tags2, 0, fn, ap );
      dslash_fn_field_special( ttt, ttt, parity, tags1, 0, fn, ap);
    }
    
    /* finish computation of M_adjoint*m*p and p*M_adjoint*m*Kp */
    /* ttt  <- ttt - msq_x4*cg_p	(msq = mass squared) */
    /* pkp  <- cg_p.(ttt - msq*cg_p) */
    pkp = 0.0;
#ifdef FEWSUMS
    c_tr=0.0; c_tt=0.0;
#endif
    FORSOMEPARITYDOMAIN(i,s,parity){
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
#ifdef FEWSUMS
    /* finally sum oldrsq over nodes, also other sums */
    tempsum[0] = pkp; tempsum[1] = c_tr; 
    tempsum[2] = c_tt; tempsum[3] = oldrsq;
    g_vecdoublesum( tempsum, 4 );
    pkp = tempsum[0]; c_tr = tempsum[1]; 
    c_tt = tempsum[2]; oldrsq = tempsum[3];
#else
    g_doublesum( &pkp );
#endif
    iteration++;
    total_iters++;
    
    a = (Real) (-rsq/pkp);
    
    /* dest <- dest + a*cg_p */
    /* resid <- resid + a*ttt */
#ifdef FEWSUMS
    actual_rsq=0.0;
#else
    rsq=0.0;
#endif
    FORSOMEPARITYDOMAIN(i,s,parity){
      if( i < loopend-FETCH_UP ){
	prefetch_VVVV( &t_dest[i+FETCH_UP], 
		       &cg_p[i+FETCH_UP], 
		       &resid[i+FETCH_UP], 
		       &ttt[i+FETCH_UP] );
      }
      scalar_mult_add_su3_vector( &t_dest[i], &cg_p[i], a, &t_dest[i] );
      scalar_mult_add_su3_vector( &resid[i], &ttt[i], a, &resid[i]);
#ifdef FEWSUMS
      actual_rsq += (double)magsq_su3vec( &resid[i] );
#else
      rsq += (double)magsq_su3vec( &resid[i] );
#endif
    } END_LOOP
#ifdef FEWSUMS
    /**printf("XXX:  node %d\t%e\t%e\t%e\n",this_node,oldrsq,c_tr,c_tt);**/
    rsq = oldrsq + 2.0*a*c_tr + a*a*c_tt; /*TEST - should equal actual_rsq */
    /**c_tt = actual_rsq;**/ /* TEMP for test */
    /**g_doublesum(&c_tt);**/ /* TEMP true value for rsq */
    /**node0_printf("RSQTEST: %e\t%e\t%e\n",rsq,c_tt,rsq-c_tt);**/
#else
    g_doublesum(&rsq);
#endif	

    if(relrsqmin > 0)
      relrsq = relative_residue(resid, t_dest, parity);
    
    qic->size_r    = (Real)rsq/source_norm;
    qic->size_relr = (Real)relrsq;

#ifdef CG_DEBUG
    if(mynode()==0){printf("iter=%d, rsq/src= %e, relrsq= %e, pkp=%e\n",
			   iteration,(double)qic->size_r,
			   (double)qic->size_relr,
			   (double)pkp);fflush(stdout);}
#endif
    
    b = (Real)rsq/oldrsq;
    /* cg_p  <- resid + b*cg_p */
    FORSOMEPARITY(i,s,parity){
      scalar_mult_add_su3_vector( &resid[i], &cg_p[i], b, &cg_p[i]);
    } END_LOOP
  }

  if(nrestart == max_restarts || iteration == max_cg){
    node0_printf("ks_congrad: CG not converged after %d iterations and %d restarts, \n",
		 iteration, nrestart);
    node0_printf("rsq. = %e wanted %e relrsq = %e wanted %e\n",
		 qic->final_rsq,rsqmin,qic->final_relrsq,relrsqmin);
    fflush(stdout);
  }

  if(special_started==1) {
    cleanup_gathers(tags1,tags2);
    special_started = 0;
  }
  cleanup_dslash_temps();

  free(ttt); free(cg_p); free(resid);
  return iteration;
}

/* API for field arguments */

int ks_congrad_field( su3_vector *src, su3_vector *dest, 
		      quark_invert_control *qic, Real mass,
		      fn_links_t *fn, ks_action_paths *ap)
{
  int iters = 0;
  double dtimec;
  double nflop = 1187;
  int parity = qic->parity;

  if(parity==EVENANDODD)nflop *=2;

  load_fn_links(fn, ap);  /* Do this here so the link build time is not
		      counted in the CG time */
  dtimec = -dclock(); 

  if(parity == EVEN || parity == EVENANDODD){
    qic->parity = EVEN;
    iters += ks_congrad_parity(src, dest, qic, mass, fn, ap);
  }
  if(parity == ODD || parity == EVENANDODD){
    qic->parity = ODD;
    iters += ks_congrad_parity(src, dest, qic, mass, fn, ap);
  }

  qic->parity = parity;

  dtimec += dclock();
#ifdef CGTIME
  if(this_node==0){
    printf("CONGRAD5: time = %e (fn %s) masses = 1 iters = %d mflops = %e\n",
	   dtimec, prec_label[PRECISION-1], iters, 
	   (double)(nflop*volume*iters/(1.0e6*dtimec*numnodes())) );
    fflush(stdout);}
#endif

  return iters;
}

/* New API for site arguments */

int ks_congrad_site( field_offset src, field_offset dest, 
		     quark_invert_control *qic, Real mass,
		     fn_links_t *fn, ks_action_paths *ap)
{
  int i;
  site *s;
  int iters = 0;
  su3_vector *t_src, *t_dest;
  double dtimec;
  double nflop = 1187;
  int parity = qic->parity;

  if(qic->parity==EVENANDODD)nflop *=2;

  load_fn_links(fn, ap);  /* Do this here so the link build time is not
		      counted in the CG time */
  dtimec = -dclock(); 

  /* Map src and dest from site to field of correct precision */
  
  t_src  = (su3_vector *)malloc(sizeof(su3_vector)*sites_on_node);
  t_dest = (su3_vector *)malloc(sizeof(su3_vector)*sites_on_node);
  if(t_src == NULL || t_dest == NULL){
    printf("ks_congrad_site(%d): No room for temporaries\n",this_node);
    terminate(1);
  }

  FORALLSITES(i,s){
    t_src[i]  = *((su3_vector *)F_PT(s,src) );
    t_dest[i] = *((su3_vector *)F_PT(s,dest));
  }

  if(parity == EVEN || parity == EVENANDODD){
    qic->parity = EVEN;
    iters += ks_congrad_parity(t_src, t_dest, qic, mass, fn, ap );
  }
  if(parity == ODD || parity == EVENANDODD){
    qic->parity = ODD;
    iters += ks_congrad_parity(t_src, t_dest, qic, mass, fn, ap );
  }

  /* Map solution to site structure */

  FORALLSITES(i,s){
    *((su3_vector *)F_PT(s,dest)) = t_dest[i];
  }

  qic->parity = parity;

  free(t_src); free(t_dest);

  dtimec += dclock();
#ifdef CGTIME
  if(this_node==0){
    printf("CONGRAD5: time = %e (fn %s) masses = 1 iters = %d mflops = %e\n",
	   dtimec, prec_label[PRECISION-1], iters, 
	   (double)(nflop*volume*iters/(1.0e6*dtimec*numnodes())) );
    fflush(stdout);}
#endif

  return iters;
}

/* Traditional MILC API for site arguments and no relative residual test */

int ks_congrad( field_offset src, field_offset dest, Real mass,
		int niter, int nrestart, Real rsqmin, int prec,
		int parity, Real *final_rsq,
		fn_links_t *fn, ks_action_paths *ap){
  int iters;
  quark_invert_control qic;

  /* Pack structure */
  qic.prec      = prec;  /* Currently ignored */
  qic.parity    = parity;
  qic.max       = niter;
  qic.nrestart  = nrestart;
  qic.resid     = rsqmin;
  qic.relresid  = 0;     /* Suppresses this test */

  /* Solve the system */
  iters = ks_congrad_site( src, dest, &qic, mass, fn, ap );

  /* Unpack the results */
  *final_rsq    = qic.final_rsq;

  return iters;
}
