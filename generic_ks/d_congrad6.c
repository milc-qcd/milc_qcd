/******* d_congrad6.c - conjugate gradient for SU3/fermions ****/
/* MIMD version 7 */
/* Kogut-Susskind fermions */

/* 4/03/00 combined with Schroedinger functional version - CD */
/* 4/03/00 Modified calling sequence to permit general src dest
   and changed name from congrad to ks_congrad CD */
/* 9/30/01 standardized prefetching CD */

/* Formerly d_congrad5.t3e.c */
/* The dslash code in this version waits for gathers from both
   positive and negative directions before computing, thereby combining
   two lattice loops in an attempt to gain prefetching time for
   sub_four_su3_vecs */

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

#include "generic_ks_includes.h"
#define LOOPEND
#include "../include/loopend.h"
#include "../include/prefetch.h"
#define FETCH_UP 1

/* Timing */
/*#define CGTIME*/
/*#define M4TIME*/
double dtimem;
int dtimem_iters;
/* #define DSLASHTIME */
/* #define DSLASHTIMES */

/* Note: precision argument "prec" is ignored here */

int ks_congrad( field_offset src, field_offset dest, Real mass,
		int niter, int nrestart, Real rsqmin, int prec, 
		int parity, Real *final_rsq_ptr ){
  register int i;
  register site *s;
  int iteration;	/* counter for iterations */
  Real a,b;	/* Sugar's a,b */
  double rsq,oldrsq,pkp;	/* resid**2,last resid*2,pkp = cg_p.K.cg_p */
  Real msq_x4;	/* 4*mass*mass */
  double source_norm;	/* squared magnitude of source vector */
  double rsqstop;	/* stopping residual normalized by source norm */
  int l_parity;	/* parity we are currently doing */
  int l_otherparity;	/* the other parity */
  msg_tag * tags1[8], *tags2[8];	/* tags for gathers to parity and opposite */
  int special_started;	/* 1 if dslash_site_special has been called */

/* Timing */

#ifdef CGTIME
double dtimed,dtimec;
#endif
double nflop;
dtimem=0.0;
dtimem_iters=0;

/* debug */
#ifdef CGTIME
 dtimec = -dclock(); 
#endif

nflop = 606;
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

	/* initialization process */
start:
	/**if(this_node==0)printf("CONGRAD: start, parity = %d\n",parity);**/
        /* ttt <-  (-1)*M_adjoint*M*dest
           resid,cg_p <- src + ttt
           rsq = |resid|^2
           source_norm = |src|^2
        */
	if(special_started==1){	/* clean up gathers */
	    for(i=XUP;i<=TUP;i++){
		cleanup_gather( tags1[i] );
		cleanup_gather( tags1[OPP_DIR(i)] );
		cleanup_gather( tags2[i] );
		cleanup_gather( tags2[OPP_DIR(i)] );
	    }
	    special_started=0;
	}
	/**if(this_node==0)if(iteration>1)printf("CONGRAD: start rsq = %.10e\n",rsq);**/
        rsq = source_norm = 0.0;
	dslash_site( dest, F_OFFSET(ttt),l_otherparity);
	dslash_site(F_OFFSET(ttt),F_OFFSET(ttt),l_parity);
	/* ttt  <- ttt - msq_x4*src	(msq = mass squared) */
	FORSOMEPARITYDOMAIN(i,s,l_parity){
	    scalar_mult_add_su3_vector( &(s->ttt), (su3_vector *)F_PT(s,dest),
					-msq_x4, &(s->ttt) );
	    add_su3_vector( (su3_vector *)F_PT(s,src), 
			    &(s->ttt), &(s->resid) );
		/* remember ttt contains -M_adjoint*M*src */
	    s->cg_p = s->resid;
	    source_norm += (double)magsq_su3vec( (su3_vector *)F_PT(s,src) );
            rsq += (double)magsq_su3vec( &(s->resid) );
	} END_LOOP
	g_doublesum( &source_norm );
        g_doublesum( &rsq );
	/**if(this_node==0)if(iteration>1)printf("CONGRAD: start rsq = %.10e\n",rsq);**/
        iteration++ ;  /* iteration counts number of multiplications
                           by M_adjoint*M */
	total_iters++;
	rsqstop = rsqmin * source_norm;
	/**if(this_node==0)printf("congrad: source_norm = %e\n",
	   (double)source_norm);**/
        if( rsq <= rsqstop ){
    	    /* if parity==EVENANDODD, set up to do odd sites and go back */
            if(parity == EVENANDODD) {
		l_parity=ODD; l_otherparity=EVEN;
		parity=EVEN;	/* so we won't loop endlessly */
		iteration = 0;
		/**if(this_node==0)printf("instant goto start\n"); **/
		goto start;
	    }
            *final_rsq_ptr=(Real)rsq;
	    /**if(this_node==0)printf("instant return\n"); fflush(stdout);**/
             return (iteration);
        }
	/**pkp=0.0;
	if(mynode()==0)printf("iter=%d, rsq= %e, pkp=%e\n",
	iteration,(double)rsq,(double)pkp);**/

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
#ifdef DSLASHTIME
dtimed = -dclock();
#endif

	if(special_started==0){
	    /**printf("CONGRAD%: calling dslash_site_special - start\n");**/
	    dslash_site_special(F_OFFSET(cg_p),F_OFFSET(ttt),l_otherparity, tags2,1);
	    dslash_site_special(F_OFFSET(ttt),F_OFFSET(ttt),l_parity,tags1,1);
	    special_started=1;
	}
	else {
	    /**printf("CONGRAD%: calling dslash_site_special - restart\n");**/
	    dslash_site_special(F_OFFSET(cg_p),F_OFFSET(ttt),l_otherparity,tags2,0);
	    dslash_site_special(F_OFFSET(ttt),F_OFFSET(ttt),l_parity,tags1,0);
	}

#ifdef DSLASHTIME
 dtimed += dclock();
if(this_node==0){printf("DSLASH: time = %e iters = %d mflops = %e\n",
dtimed,iteration , (double)(570.0*volume/(1.0e6*dtimed*numnodes())) );
fflush(stdout);} 
#endif

	/* finish computation of M_adjoint*m*p and p*M_adjoint*m*Kp */
	/* ttt  <- ttt - msq_x4*cg_p	(msq = mass squared) */
	/* pkp  <- cg_p.(ttt - msq*cg_p) */
	pkp = 0.0;
	FORSOMEPARITYDOMAIN(i,s,l_parity){
	  if(i < loopend-FETCH_UP){
	    prefetch_VV( &((s+FETCH_UP)->ttt),
			 &((s+FETCH_UP)->cg_p) );
	  }
	    scalar_mult_add_su3_vector( &(s->ttt), &(s->cg_p), -msq_x4,
		&(s->ttt) );
	    pkp += (double)su3_rdot( &(s->cg_p), &(s->ttt) );
	} END_LOOP
	g_doublesum( &pkp );
	iteration++;
	total_iters++;

	a = (Real)(-rsq/pkp);

	/* dest <- dest - a*cg_p */
	/* resid <- resid - a*ttt */
	rsq=0.0;
	FORSOMEPARITYDOMAIN(i,s,l_parity){
	  if(i < loopend-FETCH_UP){
	    prefetch_VVVV( (su3_vector *)F_PT(s+FETCH_UP,dest),
			   &((s+FETCH_UP)->cg_p),
			   &((s+FETCH_UP)->resid), 
			   &((s+FETCH_UP)->ttt) );
	  }
	    scalar_mult_add_su3_vector( (su3_vector *)F_PT(s,dest), 
					&(s->cg_p), a, 
					(su3_vector *)F_PT(s,dest) );
	    scalar_mult_add_su3_vector( &(s->resid), &(s->ttt), a, &(s->resid));
	    rsq += (double)magsq_su3vec( &(s->resid) );
	} END_LOOP
	g_doublesum(&rsq);
	/**if(mynode()==0)printf("iter=%d, rsq= %e, pkp=%e\n",
	   iteration,(double)rsq,(double)pkp);**/

        if( rsq <= rsqstop ){
    	    /* if parity==EVENANDODD, set up to do odd sites and go back */
            if(parity == EVENANDODD) {
		l_parity=ODD; l_otherparity=EVEN;
		parity=EVEN;	/* so we won't loop endlessly */
		iteration = 0;
		/**if(this_node==0)printf("normal goto start\n"); **/
		goto start;
	    }
            *final_rsq_ptr=(Real)rsq;
	    if(special_started==1){	/* clean up gathers */
		for(i=XUP;i<=TUP;i++){
		    cleanup_gather( tags1[i] );
		    cleanup_gather( tags1[OPP_DIR(i)] );
		    cleanup_gather( tags2[i] );
		    cleanup_gather( tags2[OPP_DIR(i)] );
		}
	    }
	    /**if(this_node==0)printf("normal return\n"); fflush(stdout);**/
#ifdef CGTIME
 dtimec += dclock();
if(this_node==0){printf("CONGRAD5: time = %e (naive) masses = 1 iters = %d mflops = %e\n",
dtimec,iteration,(double)(nflop*volume*iteration/(1.0e6*dtimec*numnodes())) );
fflush(stdout);}
#endif
#ifdef M4TIME
if(this_node==0){printf("M4: time = %e iters = %d mflops = %e\n",
dtimem,dtimem_iters,(double)(264*(volume/2)*dtimem_iters/(1.0e6*dtimem*numnodes())) );
fflush(stdout);}
#endif
             return (iteration);
        }

	b = (Real)(rsq/oldrsq);
	/* cg_p  <- resid + b*cg_p */
	scalar_mult_add_latvec( F_OFFSET(resid), F_OFFSET(cg_p),
	    b, F_OFFSET(cg_p), l_parity);

    } while( iteration%niter != 0);

    if( iteration < nrestart*niter ){
	/**if(this_node==0)printf("tryagain goto start\n");**/
	 goto start;
    }

    /* if parity==EVENANDODD, set up to do odd sites and go back */
    if(parity == EVENANDODD) {
	l_parity=ODD; l_otherparity=EVEN;
	parity=EVEN;	/* so we won't loop endlessly */
	iteration = 0;
	goto start;
    }

    *final_rsq_ptr=(Real)rsq;
    if(special_started==1){	/* clean up gathers */
	for(i=XUP;i<=TUP;i++){
	    cleanup_gather( tags1[i] );
	    cleanup_gather( tags1[OPP_DIR(i)] );
	    cleanup_gather( tags2[i] );
	    cleanup_gather( tags2[OPP_DIR(i)] );
	}
    }
    if(this_node==0)printf(
        "ks_congrad: CG not converged after %d iterations, res. = %e wanted %e\n",
        iteration,rsq,rsqstop);
    fflush(stdout);
    return(iteration);
}

/* D_slash routine - sets dest. on each site equal to sum of
   sources parallel transported to site, with minus sign for transport
   from negative directions */
void dslash_site( field_offset src, field_offset dest, int parity ){
register int i;
register site *s;
register int dir,otherparity;
msg_tag *tag[8];
register su3_vector *a,*b1,*b2,*b3,*b4;

    switch(parity){
	case EVEN:	otherparity=ODD; break;
	case ODD:	otherparity=EVEN; break;
	case EVENANDODD:	otherparity=EVENANDODD; break;
    }

    /* Start gathers from positive directions */
    for(dir=XUP; dir<=TUP; dir++){
	tag[dir] = start_gather_site( src, sizeof(su3_vector), dir, parity,
	    gen_pt[dir] );
    }

#ifdef M4TIME
dtimem -= dclock();
#endif
    /* Multiply by adjoint matrix at other sites */
    FORSOMEPARITYDOMAIN(i,s,otherparity){
      if(i < loopend-FETCH_UP){
        prefetch_4MV4V( &((s+FETCH_UP)->link[XUP]),
			(su3_vector *)F_PT((s+FETCH_UP),src),
			&((s+FETCH_UP)->tempvec));
      }
	mult_adj_su3_mat_vec_4dir( s->link,
	    (su3_vector *)F_PT(s,src), s->tempvec );
    } END_LOOP
#ifdef M4TIME
dtimem += dclock();
if(otherparity==EVENANDODD)dtimem_iters +=2; else dtimem_iters++;
#endif

    /* Start gathers from negative directions */
    for( dir=XUP; dir <= TUP; dir++){
	tag[OPP_DIR(dir)] = start_gather_site( F_OFFSET(tempvec[dir]),
	    sizeof(su3_vector), OPP_DIR( dir), parity,
	    gen_pt[OPP_DIR(dir)] );
    }

    /* Wait gathers from positive directions, multiply by matrix and
	accumulate */
    for(dir=XUP; dir<=TUP; dir++){
	wait_gather(tag[dir]);
    }
    /*----------------- (departure from d_congrad.c) -------------*/
    /* Wait gathers from negative directions, accumulate (negative) */
    for(dir=XUP; dir<=TUP; dir++){
	wait_gather(tag[OPP_DIR(dir)]);
    }
    /*------------------------------------------------------------*/
#ifdef SCHROED_FUN
    FORSOMEPARITY(i,s,parity) if(s->t > 0){
      if( i < loopend-FETCH_UP ){
	prefetch_V((su3_vector *)F_PT(s+FETCH_UP,dest));
	prefetch_4MVVVV(&((s+FETCH_UP)->link[XUP]),
			(su3_vector *)gen_pt[XUP][i+FETCH_UP],
			(su3_vector *)gen_pt[YUP][i+FETCH_UP],
			(su3_vector *)gen_pt[ZUP][i+FETCH_UP],
			(su3_vector *)gen_pt[TUP][i+FETCH_UP] );
      }
      /* Don't include last time-link in computation */
      if(s->t == (nt-1)){
	    mult_su3_mat_vec( &(s->link[XUP]),
		(su3_vector *)(gen_pt[XUP][i]), (su3_vector *)F_PT(s,dest));
	    for(dir=YUP; dir<TUP; dir++){
		mult_su3_mat_vec_sum( &(s->link[dir]),
		    (su3_vector *)(gen_pt[dir][i]), (su3_vector *)F_PT(s,dest));
	    }
	}
	else{
#else
    FORSOMEPARITY(i,s,parity){
      if( i < loopend-FETCH_UP ){
	prefetch_V((su3_vector *)F_PT(s+FETCH_UP,dest));
	prefetch_4MVVVV(&((s+FETCH_UP)->link[XUP]),
			(su3_vector *)gen_pt[XUP][i+FETCH_UP],
			(su3_vector *)gen_pt[YUP][i+FETCH_UP],
			(su3_vector *)gen_pt[ZUP][i+FETCH_UP],
			(su3_vector *)gen_pt[TUP][i+FETCH_UP] );
      }
#endif
      mult_su3_mat_vec_sum_4dir( s->link,
          (su3_vector *)gen_pt[XUP][i], (su3_vector *)gen_pt[YUP][i],
          (su3_vector *)gen_pt[ZUP][i], (su3_vector *)gen_pt[TUP][i],
          (su3_vector *)F_PT(s,dest));
    /*------------------------------------------------------------*/

        if( i < loopend-FETCH_UP ){
	  prefetch_VVVV((su3_vector *)gen_pt[XDOWN][i+FETCH_UP],
			(su3_vector *)gen_pt[YDOWN][i+FETCH_UP],
			(su3_vector *)gen_pt[ZDOWN][i+FETCH_UP],
			(su3_vector *)gen_pt[TDOWN][i+FETCH_UP] );
	}
#ifndef INLINE
      sub_four_su3_vecs( (su3_vector *)F_PT(s,dest),
	    (su3_vector *)(gen_pt[XDOWN][i]),
	    (su3_vector *)(gen_pt[YDOWN][i]),
	    (su3_vector *)(gen_pt[ZDOWN][i]),
	    (su3_vector *)(gen_pt[TDOWN][i]) );

#else

      /* Inline version */
      a =  (su3_vector *)F_PT(s,dest);
      b1 = (su3_vector *)(gen_pt[XDOWN][i]);
      b2 = (su3_vector *)(gen_pt[YDOWN][i]);
      b3 = (su3_vector *)(gen_pt[ZDOWN][i]);
      b4 = (su3_vector *)(gen_pt[TDOWN][i]);

      CSUB(a->c[0], b1->c[0], a->c[0]);
      CSUB(a->c[1], b1->c[1], a->c[1]);
      CSUB(a->c[2], b1->c[2], a->c[2]);
      
      CSUB(a->c[0], b2->c[0], a->c[0]);
      CSUB(a->c[1], b2->c[1], a->c[1]);
      CSUB(a->c[2], b2->c[2], a->c[2]);
      
      CSUB(a->c[0], b3->c[0], a->c[0]);
      CSUB(a->c[1], b3->c[1], a->c[1]);
      CSUB(a->c[2], b3->c[2], a->c[2]);
      
      CSUB(a->c[0], b4->c[0], a->c[0]);
      CSUB(a->c[1], b4->c[1], a->c[1]);
      CSUB(a->c[2], b4->c[2], a->c[2]);
#endif

    /*------------------------------------------------------------*/
#ifdef SCHROED_FUN
	}
#endif
    } END_LOOP

    /*------------------------------------------------------------*/
    /* free up the buffers */
    for(dir=XUP; dir<=TUP; dir++){
	cleanup_gather(tag[dir]);
	cleanup_gather(tag[OPP_DIR(dir)]);
    }
}

/* Special dslash_site for use by congrad.  Uses restart_gather_site() when
  possible. Last argument is an array of message tags, to be set
  if this is the first use, otherwise reused. If start=1,use
  start_gather_site, otherwise use restart_gather_site. 
  The calling program must clean up the gathers! */
void dslash_site_special( field_offset src, field_offset dest,
    int parity, msg_tag **tag, int start ){
register int i;
register site *s;
register int dir,otherparity;
register su3_vector *a,*b1,*b2,*b3,*b4;

#ifdef DSLASHTIMES
double dtime0,dtime1,dtime2,dtime3,dtime4,dtime5,dtime6,dclock();
#endif
    switch(parity){
	case EVEN:	otherparity=ODD; break;
	case ODD:	otherparity=EVEN; break;
	case EVENANDODD:	otherparity=EVENANDODD; break;
    }

#ifdef DSLASHTIMES
 dtime0 = -dclock(); 
#endif

    /* Start gathers from positive directions */
    for(dir=XUP; dir<=TUP; dir++){
/**printf("dslash_site_special: up gathers, start=%d\n",start);**/
	if(start==1) tag[dir] = start_gather_site( src, sizeof(su3_vector),
	    dir, parity, gen_pt[dir] );
	else restart_gather_site( src, sizeof(su3_vector),
	    dir, parity, gen_pt[dir] , tag[dir] );
    }

#ifdef DSLASHTIMES
    dtime1 = -dclock(); 
    dtime0 -= dtime1;
#endif


#ifdef M4TIME
dtimem -= dclock();
#endif
    /* Multiply by adjoint matrix at other sites */
    FORSOMEPARITYDOMAIN(i,s,otherparity){
      if( i < loopend-FETCH_UP ){
	prefetch_4MV4V( &((s+FETCH_UP)->link[XUP]),
			(su3_vector *)F_PT((s+FETCH_UP),src),
			&((s+FETCH_UP)->tempvec));
      }
	mult_adj_su3_mat_vec_4dir( s->link,
	    (su3_vector *)F_PT(s,src), s->tempvec );
    } END_LOOP

#ifdef M4TIME
dtimem += dclock();
if(otherparity==EVENANDODD)dtimem_iters +=2; else dtimem_iters++;
#endif

#ifdef DSLASHTIMES
    dtime2 = -dclock();
    dtime1 -= dtime2;
#endif


    /* Start gathers from negative directions */
    for( dir=XUP; dir <= TUP; dir++){
/**printf("dslash_site_special: down gathers, start=%d\n",start);**/
	if (start==1) tag[OPP_DIR(dir)] = start_gather_site( F_OFFSET(tempvec[dir]),
	    sizeof(su3_vector), OPP_DIR( dir), parity, gen_pt[OPP_DIR(dir)] );
	else restart_gather_site( F_OFFSET(tempvec[dir]), sizeof(su3_vector),
	    OPP_DIR( dir), parity, gen_pt[OPP_DIR(dir)] , tag[OPP_DIR(dir)] );
    }


#ifdef DSLASHTIMES
   dtime3 = -dclock(); 
   dtime2 -= dtime3;
#endif

    /* Wait gathers from positive directions, multiply by matrix and
	accumulate */
    for(dir=XUP; dir<=TUP; dir++){
	wait_gather(tag[dir]);
    }
    /*------------------------------------------------------------*/
    /* Wait gathers from negative directions, accumulate (negative) */
    for(dir=XUP; dir<=TUP; dir++){
	wait_gather(tag[OPP_DIR(dir)]);
    }
#ifdef DSLASHTIMES
    dtime4 = -dclock(); 
    dtime3 -= dtime4;
#endif


    /*------------------------------------------------------------*/
#ifdef SCHROED_FUN
    FORSOMEPARITY(i,s,parity) if(s->t > 0){
	if(s->t == (nt-1)){
	    mult_su3_mat_vec( &(s->link[XUP]),
		(su3_vector *)(gen_pt[XUP][i]), (su3_vector *)F_PT(s,dest));
	    for(dir=YUP; dir<TUP; dir++){
		mult_su3_mat_vec_sum( &(s->link[dir]),
		    (su3_vector *)(gen_pt[dir][i]), (su3_vector *)F_PT(s,dest));
	    }
	}
	else{
#else
    FORSOMEPARITY(i,s,parity){
      if( i < loopend-FETCH_UP ){
	prefetch_VVVV( (su3_vector *)gen_pt[XUP][i+FETCH_UP],
		       (su3_vector *)gen_pt[YUP][i+FETCH_UP],
		       (su3_vector *)gen_pt[ZUP][i+FETCH_UP],
		       (su3_vector *)gen_pt[TUP][i+FETCH_UP] );
      }
#endif
      mult_su3_mat_vec_sum_4dir( s->link,
	(su3_vector *)gen_pt[XUP][i], (su3_vector *)gen_pt[YUP][i],
	(su3_vector *)gen_pt[ZUP][i], (su3_vector *)gen_pt[TUP][i],
	(su3_vector *)F_PT(s,dest));

    /*------------------------------------------------------------*/
	if( i < loopend-FETCH_UP ){
	    prefetch_VVVV((su3_vector *)gen_pt[XDOWN][i+FETCH_UP],
			  (su3_vector *)gen_pt[YDOWN][i+FETCH_UP],
			  (su3_vector *)gen_pt[ZDOWN][i+FETCH_UP],
			  (su3_vector *)gen_pt[TDOWN][i+FETCH_UP] );
	}

#ifndef INLINE

      sub_four_su3_vecs( (su3_vector *)F_PT(s,dest),
	    (su3_vector *)(gen_pt[XDOWN][i]),
	    (su3_vector *)(gen_pt[YDOWN][i]),
	    (su3_vector *)(gen_pt[ZDOWN][i]),
	    (su3_vector *)(gen_pt[TDOWN][i]) ); 
#else

      /* Inline version */
      a =  (su3_vector *)F_PT(s,dest);
      b1 = (su3_vector *)(gen_pt[XDOWN][i]);
      b2 = (su3_vector *)(gen_pt[YDOWN][i]);
      b3 = (su3_vector *)(gen_pt[ZDOWN][i]);
      b4 = (su3_vector *)(gen_pt[TDOWN][i]);

      CSUB(a->c[0], b1->c[0], a->c[0]);
      CSUB(a->c[1], b1->c[1], a->c[1]);
      CSUB(a->c[2], b1->c[2], a->c[2]);
      
      CSUB(a->c[0], b2->c[0], a->c[0]);
      CSUB(a->c[1], b2->c[1], a->c[1]);
      CSUB(a->c[2], b2->c[2], a->c[2]);
      
      CSUB(a->c[0], b3->c[0], a->c[0]);
      CSUB(a->c[1], b3->c[1], a->c[1]);
      CSUB(a->c[2], b3->c[2], a->c[2]);
      
      CSUB(a->c[0], b4->c[0], a->c[0]);
      CSUB(a->c[1], b4->c[1], a->c[1]);
      CSUB(a->c[2], b4->c[2], a->c[2]);
#endif
    /*------------------------------------------------------------*/
#ifdef SCHROED_FUN
	}
#endif
    } END_LOOP

#ifdef DSLASHTIMES
    dtime4 +=dclock();
    if(this_node==0){printf("DSLASHTIMES = %e %e %e %e %e\n",
dtime0,dtime1,dtime2,dtime3,dtime4);
fflush(stdout);} 
#endif

} /* d_congrad6.c */
