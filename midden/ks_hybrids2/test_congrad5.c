/******* test_congrad5.c - conjugate gradient for SU3/fermions ****/
/* MIMD version 6 */
/* NOT MAINTAINED.  TEST BEFORE USE! */

/* Kogut-Susskind fermions */

/* This version looks at the initial vector every "niter" passes */
/* The source vector is in "phi", and the initial guess and answer
   in "xxx".  "resid" is the residual vector, and "cg_p" and "ttt" are
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
#include "ks_hyb_includes.h"

int congrad(niter,rsqmin,parity,final_rsq_ptr) 
int niter,parity; Real rsqmin,*final_rsq_ptr;
{
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
void clear_latvec(), copy_latvec(), scalar_mult_add_latvec();
Real magsq_latvec(), realdot_latvec();
void dslash();
register int first,count;
msg_tag * tags1[8], *tags2[8];	/* tags for gathers to parity and opposite */
int special_started;	/* 1 if dslash_special has been called */
void dslash_special();

/**double dtime,dclock();
dtime = -dclock();**/
	
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
	/**if(this_node==0)printf("congrad4: start, parity = %d\n",parity);**/
        /* ttt <-  (-1)*M_adjoint*M*xxx
           resid,cg_p <- phi + ttt
           rsq = |resid|^2
           source_norm = |phi|^2
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
	dslash(F_OFFSET(xxx),F_OFFSET(ttt),l_otherparity);
	dslash(F_OFFSET(ttt),F_OFFSET(ttt),l_parity);
	/* ttt  <- ttt - msq_x4*phi	(msq = mass squared) */
	FORSOMEPARITY(i,s,l_parity){
	    register int j;
	    scalar_mult_add_su3_vector( &(s->ttt), &(s->xxx), -msq_x4,
		&(s->ttt) );
	    add_su3_vector( &(s->phi), &(s->ttt), &(s->resid) );
		/* remember ttt contains -M_adjoint*M*phi */
	    s->cg_p = s->resid;
	    source_norm += (double)magsq_su3vec( &(s->phi) );
            rsq += (double)magsq_su3vec( &(s->resid) );
	}
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
            *final_rsq_ptr=rsq;
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
           xxx <- xxx + a*cg_p
           resid <- resid + a*ttt
           rsq <- |resid|^2
           b <- rsq/oldrsq
           cg_p <- resid + b*cg_p
        */
    do{
        oldrsq = rsq;
        pkp = 0.0;
	/* sum of neighbors */
/**dtime = -dclock();**/
	if(special_started==0){
	    /**printf("CONGRAD%: calling dslash_special - start\n");**/
	    dslash_special(F_OFFSET(cg_p),F_OFFSET(ttt),l_otherparity, tags2,1);
	    dslash_special(F_OFFSET(ttt),F_OFFSET(ttt),l_parity,tags1,1);
	    special_started=1;
	}
	else {
	    /**printf("CONGRAD%: calling dslash_special - restart\n");**/
	    dslash_special(F_OFFSET(cg_p),F_OFFSET(ttt),l_otherparity,tags2,0);
	    dslash_special(F_OFFSET(ttt),F_OFFSET(ttt),l_parity,tags1,0);
	}
/**dtime += dclock();
if(this_node==0)printf("DSLASH: time = %e iters = %d mflops = %e\n",
dtime,iteration , (double)(570.0*volume/(1.0e6*dtime*numnodes())) );**/
	/* finish computation of M_adjoint*m*p and p*M_adjoint*m*Kp */
	/* ttt  <- ttt - msq_x4*cg_p	(msq = mass squared) */
	/* pkp  <- cg_p.(ttt - msq*cg_p) */
	pkp = 0.0;
	FORSOMEPARITY(i,s,l_parity){
	    register int j;
	    scalar_mult_add_su3_vector( &(s->ttt), &(s->cg_p), -msq_x4,
		&(s->ttt) );
	    pkp += (double)su3_rdot( &(s->cg_p), &(s->ttt) );
	}
	g_doublesum( &pkp );
	iteration++;
	total_iters++;

	a = (Real)(-rsq/pkp);

	/* xxx <- xxx - a*cg_p */
	/* resid <- resid - a*ttt */
	rsq=0.0;
	FORSOMEPARITY(i,s,l_parity){
	    /**/
	    scalar_mult_add_su3_vector( &(s->xxx), &(s->cg_p), a, &(s->xxx) );
	    scalar_mult_add_su3_vector( &(s->resid), &(s->ttt), a, &(s->resid));
	    /**/
	    /**
	    scalar_mult_sum_su3_vector( &(s->xxx), &(s->cg_p), a);
	    scalar_mult_sum_su3_vector( &(s->resid), &(s->ttt), a);
	    **/
	    rsq += (double)magsq_su3vec( &(s->resid) );
	}
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
            *final_rsq_ptr=rsq;
	    if(special_started==1){	/* clean up gathers */
		for(i=XUP;i<=TUP;i++){
		    cleanup_gather( tags1[i] );
		    cleanup_gather( tags1[OPP_DIR(i)] );
		    cleanup_gather( tags2[i] );
		    cleanup_gather( tags2[OPP_DIR(i)] );
		}
	    }
	    /**if(this_node==0)printf("normal return\n"); fflush(stdout);**/
/**dtime += dclock();
if(this_node==0)printf("CONGRAD5: time = %e iters = %d mflops = %e\n",
dtime,iteration,(double)(606.0*volume*iteration/(1.0e6*dtime*numnodes())) );**/
             return (iteration);
        }

	b = (Real)(rsq/oldrsq);
	/* cg_p  <- resid + b*cg_p */
	scalar_mult_add_latvec( F_OFFSET(resid), F_OFFSET(cg_p),
	    b, F_OFFSET(cg_p), l_parity);

    } while( iteration%niter != 0);

    if( iteration < 1*niter ){
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

    *final_rsq_ptr=rsq;
    if(special_started==1){	/* clean up gathers */
	for(i=XUP;i<=TUP;i++){
	    cleanup_gather( tags1[i] );
	    cleanup_gather( tags1[OPP_DIR(i)] );
	    cleanup_gather( tags2[i] );
	    cleanup_gather( tags2[OPP_DIR(i)] );
	}
    }
    if(this_node==0)printf(
        "CG not converged after %d iterations, res. = %e wanted %e\n",
        iteration,rsq,rsqstop);
    fflush(stdout);
    return(iteration);
}

/* clear an su3_vector in the lattice */
void clear_latvec(v,parity) field_offset v; int parity;{
register int i,j;
register site *s;
register su3_vector *vv;
    switch(parity){
	case EVEN: FOREVENSITES(i,s){
		vv = (su3_vector *)F_PT(s,v);
		for(j=0;j<3;j++){ vv->c[j].real = vv->c[j].imag = 0.0; }
	    } break;
	case ODD: FORODDSITES(i,s){
		vv = (su3_vector *)F_PT(s,v);
		for(j=0;j<3;j++){ vv->c[j].real = vv->c[j].imag = 0.0; }
	    } break;
	case EVENANDODD: FORALLSITES(i,s){
		vv = (su3_vector *)F_PT(s,v);
		for(j=0;j<3;j++){ vv->c[j].real = vv->c[j].imag = 0.0; }
	    } break;
    } 
}

/* copy an su3_vector in the lattice */
void copy_latvec(src,dest,parity) field_offset src,dest; int parity;{
register int i;
register site *s;
register su3_vector *spt,*dpt;
    switch(parity){
	case EVEN: FOREVENSITES(i,s){
		s = &(lattice[i]);
		spt = (su3_vector *)F_PT(s,src);
		dpt = (su3_vector *)F_PT(s,dest);
		*dpt = *spt;
	    } break;
	case ODD: FORODDSITES(i,s){
		s = &(lattice[i]);
		spt = (su3_vector *)F_PT(s,src);
		dpt = (su3_vector *)F_PT(s,dest);
		*dpt = *spt;
	    } break;
	case EVENANDODD: FORALLSITES(i,s){
		s = &(lattice[i]);
		spt = (su3_vector *)F_PT(s,src);
		dpt = (su3_vector *)F_PT(s,dest);
		*dpt = *spt;
	    } break;
    } 
}

/* scalar multiply and add an SU3 vector in the lattice */
void scalar_mult_add_latvec(src1,src2,scalar,dest,parity)
field_offset src1,src2,dest;
Real scalar;
int parity;
{
register int i;
register site *s;
register su3_vector *spt1,*spt2,*dpt;
    switch(parity){
	case EVEN: FOREVENSITES(i,s){
		spt1 = (su3_vector *)F_PT(s,src1);
		spt2 = (su3_vector *)F_PT(s,src2);
		dpt = (su3_vector *)F_PT(s,dest);
		scalar_mult_add_su3_vector( spt1 , spt2 , scalar , dpt );
	    } break;
	case ODD: FORODDSITES(i,s){
		spt1 = (su3_vector *)F_PT(s,src1);
		spt2 = (su3_vector *)F_PT(s,src2);
		dpt = (su3_vector *)F_PT(s,dest);
		scalar_mult_add_su3_vector( spt1 , spt2 , scalar , dpt );
	    } break;
	case EVENANDODD: FORALLSITES(i,s){
		spt1 = (su3_vector *)F_PT(s,src1);
		spt2 = (su3_vector *)F_PT(s,src2);
		dpt = (su3_vector *)F_PT(s,dest);
		scalar_mult_add_su3_vector( spt1 , spt2 , scalar , dpt );
	    } break;
	} 
}

/* scalar multiply an SU3 vector in the lattice */
void scalar_mult_latvec(src,scalar,dest,parity)
field_offset src,dest;
Real scalar;
int parity;
{
register int i;
register site *s;
register su3_vector *spt,*dpt;
    switch(parity){
	case EVEN: FOREVENSITES(i,s){
		spt = (su3_vector *)F_PT(s,src);
		dpt = (su3_vector *)F_PT(s,dest);
		scalar_mult_su3_vector( spt , scalar , dpt );
	    } break;
	case ODD: FORODDSITES(i,s){
		spt = (su3_vector *)F_PT(s,src);
		dpt = (su3_vector *)F_PT(s,dest);
		scalar_mult_su3_vector( spt , scalar , dpt );
	    } break;
	case EVENANDODD: FORALLSITES(i,s){
		spt = (su3_vector *)F_PT(s,src);
		dpt = (su3_vector *)F_PT(s,dest);
		scalar_mult_su3_vector( spt , scalar , dpt );
	    } break;
    } 
}

/* D_slash routine - sets dest. on each site equal to sum of
   sources parallel transported to site, with minus sign for transport
   from negative directions */
void dslash(src,dest,parity) field_offset src,dest; int parity; {
register int i;
register site *s;
register int dir,otherparity;
msg_tag *tag[8];
register int first,count;

    switch(parity){
	case EVEN:	otherparity=ODD; break;
	case ODD:	otherparity=EVEN; break;
	case EVENANDODD:	otherparity=EVENANDODD; break;
    }

    /* Start gathers from positive directions */
    for(dir=XUP; dir<=TUP; dir++){
	tag[dir] = start_gather( src, sizeof(su3_vector), dir, parity,
	    gen_pt[dir] );
    }

    /* Multiply by adjoint matrix at other sites */
    FORSOMEPARITY(i,s,otherparity){
	/**
	for(dir=XUP; dir<=TUP; dir++){
	    mult_adj_su3_mat_vec( &(s->link[dir]),
		(su3_vector *)F_PT(s,src), &(s->tempvec[dir]) );
	}
	**/
	mult_adj_su3_mat_vec_4dir( s->link,
	    (su3_vector *)F_PT(s,src), s->tempvec );
    }

    /* Start gathers from negative directions */
    for( dir=XUP; dir <= TUP; dir++){
	tag[OPP_DIR(dir)] = start_gather( F_OFFSET(tempvec[dir]),
	    sizeof(su3_vector), OPP_DIR( dir), parity,
	    gen_pt[OPP_DIR(dir)] );
    }

    /* Wait gathers from positive directions, multiply by matrix and
	accumulate */
    for(dir=XUP; dir<=TUP; dir++){
	wait_gather(tag[dir]);
    }
    FORSOMEPARITY(i,s,parity){
	/**
	mult_su3_mat_vec( &(s->link[XUP]),
	    (su3_vector *)(gen_pt[XUP][i]), (su3_vector *)F_PT(s,dest));
	for(dir=YUP; dir<=TUP; dir++){
	    mult_su3_mat_vec_sum( &(s->link[dir]),
	       (su3_vector *)(gen_pt[dir][i]), (su3_vector *)F_PT(s,dest));
	}
	**/
	mult_su3_mat_vec_sum_4dir( s->link,
	    gen_pt[XUP][i], gen_pt[YUP][i], gen_pt[ZUP][i], gen_pt[TUP][i],
	    (su3_vector *)F_PT(s,dest));
    }

    /* Wait gathers from negative directions, accumulate (negative) */
    for(dir=XUP; dir<=TUP; dir++){
	wait_gather(tag[OPP_DIR(dir)]);
    }
    FORSOMEPARITY(i,s,parity){
	sub_four_su3_vecs( (su3_vector *)F_PT(s,dest),
	    (su3_vector *)(gen_pt[XDOWN][i]),
	    (su3_vector *)(gen_pt[YDOWN][i]),
	    (su3_vector *)(gen_pt[ZDOWN][i]),
	    (su3_vector *)(gen_pt[TDOWN][i]) );
    }

    /* free up the buffers */
    for(dir=XUP; dir<=TUP; dir++){
	cleanup_gather(tag[dir]);
	cleanup_gather(tag[OPP_DIR(dir)]);
    }
}

/* Special dslash for use by congrad.  Uses restart_gather() when
  possible. Last argument is an array of message tags, to be set
  if this is the first use, otherwise reused. If start=1,use
  start_gather, otherwise use restart_gather. 
  The calling program must clean up the gathers! */
void dslash_special(src,dest,parity,tag,start)
field_offset src,dest; int parity; msg_tag **tag; int start; {
register int i;
register site *s;
register int dir,otherparity;
register int first,count;

    switch(parity){
	case EVEN:	otherparity=ODD; break;
	case ODD:	otherparity=EVEN; break;
	case EVENANDODD:	otherparity=EVENANDODD; break;
    }

    /* Start gathers from positive directions */
    for(dir=XUP; dir<=TUP; dir++){
/**printf("dslash_special: up gathers, start=%d\n",start);**/
	if(start==1) tag[dir] = start_gather( src, sizeof(su3_vector),
	    dir, parity, gen_pt[dir] );
	else restart_gather( src, sizeof(su3_vector),
	    dir, parity, gen_pt[dir] , tag[dir] );
    }

    /* Multiply by adjoint matrix at other sites */
    FORSOMEPARITY(i,s,otherparity){
	/**
	for(dir=XUP; dir<=TUP; dir++){
	    mult_adj_su3_mat_vec( &(s->link[dir]),
		(su3_vector *)F_PT(s,src), &(s->tempvec[dir]) );
	}
	**/
	mult_adj_su3_mat_vec_4dir( s->link,
	    (su3_vector *)F_PT(s,src), s->tempvec );
    }

    /* Start gathers from negative directions */
    for( dir=XUP; dir <= TUP; dir++){
/**printf("dslash_special: down gathers, start=%d\n",start);**/
	if (start==1) tag[OPP_DIR(dir)] = start_gather( F_OFFSET(tempvec[dir]),
	    sizeof(su3_vector), OPP_DIR( dir), parity, gen_pt[OPP_DIR(dir)] );
	else restart_gather( F_OFFSET(tempvec[dir]), sizeof(su3_vector),
	    OPP_DIR( dir), parity, gen_pt[OPP_DIR(dir)] , tag[OPP_DIR(dir)] );
    }

    /* Wait gathers from positive directions, multiply by matrix and
	accumulate */
    for(dir=XUP; dir<=TUP; dir++){
	wait_gather(tag[dir]);
    }
    FORSOMEPARITY(i,s,parity){
	/**
	mult_su3_mat_vec( &(s->link[XUP]),
	    (su3_vector *)(gen_pt[XUP][i]), (su3_vector *)F_PT(s,dest));
	for(dir=YUP; dir<=TUP; dir++){
	    mult_su3_mat_vec_sum( &(s->link[dir]),
	       (su3_vector *)(gen_pt[dir][i]), (su3_vector *)F_PT(s,dest));
	}
	**/
	mult_su3_mat_vec_sum_4dir( s->link,
	    gen_pt[XUP][i], gen_pt[YUP][i], gen_pt[ZUP][i], gen_pt[TUP][i],
	    (su3_vector *)F_PT(s,dest));
    }

    /* Wait gathers from negative directions, accumulate (negative) */
    for(dir=XUP; dir<=TUP; dir++){
	wait_gather(tag[OPP_DIR(dir)]);
    }
    FORSOMEPARITY(i,s,parity){
	sub_four_su3_vecs( (su3_vector *)F_PT(s,dest),
	    (su3_vector *)(gen_pt[XDOWN][i]),
	    (su3_vector *)(gen_pt[YDOWN][i]),
	    (su3_vector *)(gen_pt[ZDOWN][i]),
	    (su3_vector *)(gen_pt[TDOWN][i]) );
    }

}

