/************************** w_loop1.c *******************************/
/* MIMD version 6 */
/* This version uses gathers to get the neighbors */
/* version of 3/14/94 by UMH */
/* 2/19/98 Version 5 port CD */

/* Computes time-like, planar Wilson loops on gauge configuration
   in axial gauge with time-like gauge fields of last time slice in
   all other time slices as well, instead of the unit matrix! */

#include "string_break_includes.h"

void w_loop1(int tot_smear) {

register int i,dir,r,t;
int nth,nxh;
register site *s;
su3_matrix tmat1,tmat2;
msg_tag *mtag[4];
Real *wils_loop1,ftmp;

    if( nx != ny || nx != ny){
	if(this_node == 0)printf("w_loop1 gives wrong results for nx!=ny!=nz");
	return;
    }

    nth = nt/2;  nxh = nx/2;
    wils_loop1 = (Real *)malloc(nth*nxh*sizeof(Real));
    for(t=0;t<nth;t++) for(r=0;r<nxh;r++){
	wils_loop1[r+nxh*t] = 0.0;
    }

    for(dir=XUP;dir<=ZUP;dir++){

	FORALLSITES(i,s){
	    su3mat_copy( &(s->link[dir]), &(s->s_link));
	    su3mat_copy( &(s->link[TUP]), &(s->t_link_f));
	}

	/* Start gather of forward time-like links */
	mtag[0] = start_gather_site( F_OFFSET(t_link_f), sizeof(su3_matrix),
	    dir, EVENANDODD, gen_pt[0] );

	/* Recursively construct the space-like segments and compute
	   the Wilson loops with that segment */

	for(r=0;r<nxh;r++){

	    if( r>0 ){
		wait_gather( mtag[1]);
		FORALLSITES(i,s){
		    su3mat_copy( (su3_matrix *)(gen_pt[1][i]), &(s->staple));
		}
		FORALLSITES(i,s){
		    mult_su3_nn( &(s->link[dir]), &(s->staple), &(s->s_link));
		}
	    }

	    FORALLSITES(i,s){
		su3mat_copy( &(s->s_link), &(s->s_link_f));
	    }

	    /* Start gather of forward space-like segments */
	    mtag[TUP] = start_gather_site( F_OFFSET(s_link_f), sizeof(su3_matrix),
		TUP, EVENANDODD, gen_pt[TUP] );

	    /* Inbetween gather space-links for next r, if still needed. */
	    if( r==0 ){
		mtag[1] = start_gather_site( F_OFFSET(s_link), sizeof(su3_matrix),
		    dir, EVENANDODD, gen_pt[1] );
	    }
	    else if( r<(nxh-1) ){
		restart_gather_site( F_OFFSET(s_link), sizeof(su3_matrix),
		    dir, EVENANDODD, gen_pt[1], mtag[1] );
	    }
	    else{
		cleanup_gather( mtag[1]);
	    }

	    /* Collect forward time-like links. */
	    wait_gather( mtag[0]);
	    FORALLSITES(i,s){
		su3mat_copy( (su3_matrix *)(gen_pt[0][i]), &(s->staple));
	    }
	    FORALLSITES(i,s){
		su3mat_copy( &(s->staple), &(s->t_link_f));
	    }

	    /* Recursively compute the Wilson loops of different time extent */
	    for(t=0;t<nth;t++){

		/* Collect forward space-like segments */
		wait_gather( mtag[TUP]);
		FORALLSITES(i,s){
		    su3mat_copy( (su3_matrix *)(gen_pt[TUP][i]), &(s->staple));
		}
		FORALLSITES(i,s){
		    su3mat_copy( &(s->staple), &(s->s_link_f));
		}

		/* Start gather for next t, if still needed. */
		if( t<(nth-1) ){
		    restart_gather_site( F_OFFSET(s_link_f), sizeof(su3_matrix),
			TUP, EVENANDODD, gen_pt[TUP], mtag[TUP] );
		}
		else{
		    cleanup_gather( mtag[TUP]);
		}

		/* Finally, compute the Wilson loops. */
		FORALLSITES(i,s){
		    if( ((s->t)+t+1)>=nt ){
			mult_su3_nn( &(s->link[TUP]), &(s->s_link_f), &tmat1);
			mult_su3_na( &tmat1, &(s->t_link_f), &tmat2);
			wils_loop1[r+nxh*t] +=
			    realtrace_su3( &tmat2, &(s->s_link));
		    }
		    else
		    {
			wils_loop1[r+nxh*t] +=
			    realtrace_su3( &(s->s_link_f), &(s->s_link));
		    }
		}

	    } /* end loop over t */

	    /* Start gather of forward time-like links for next r. */
	    if( r<(nxh-1) ){
		restart_gather_site( F_OFFSET(t_link_f), sizeof(su3_matrix),
		    dir, EVENANDODD, gen_pt[0], mtag[0] );
	    }
	    else{
		cleanup_gather( mtag[0]);
	    }

	} /* end loop over r */

    } /* end loop over dir */

    /* Normalize and print the Wilson loops */
    for(t=0;t<nth;t++) for(r=0;r<nxh;r++){
	ftmp = wils_loop1[r+nxh*t];
	g_floatsum( &ftmp);
	ftmp /= (Real)(9*volume);
	if(this_node == 0)printf("WILS_LOOP1_%d  %d  %d  %e\n",
	    tot_smear, r, t, (double)ftmp);
    }

    free( wils_loop1);

} /* w_loop1 */

