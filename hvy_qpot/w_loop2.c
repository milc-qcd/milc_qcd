/************************** w_loop2.c *******************************/
/* MIMD version 7 */
/* This version uses gathers to get the neighbors */
/* version of 3/16/94 by UMH */
/* 2/19/98 Version 5 port CD */

/* Computes time-like, off-axis Wilson loops on gauge configuration
   in axial gauge with time-like gauge fields of last time slice in
   all other time slices as well, instead of the unit matrix! */

#include "hvy_qpot_includes.h"

void w_loop2(int tot_smear) {

  register int i,dir1,dir2,dir3,r,t,r_off;
  int nth,nxh,nrmax;
  register site *s;
  su3_matrix tmat1,tmat2;
  msg_tag *mtag[8],*gmtag;
  Real *wils_loop2,ftmp;
  int disp[4];    /* displacement vector for general gather */
  su3_matrix *s_link, *s_link_f, *t_link_f;

    if( nx != ny || nx != ny){
	if(this_node == 0)printf("w_loop2 gives wrong results for nx!=ny!=nz");
        return;
    }

    nth = nt/2;  nxh = nx/2;  nrmax = 2*nxh + nxh/2;
    wils_loop2 = (Real *)malloc(nth*nrmax*sizeof(Real));
    for(t=0;t<nth;t++) for(r=0;r<nrmax;r++){
	wils_loop2[r+nrmax*t] = 0.0;
    }

    /* Allocate space for space-link product */
    s_link = 
      (su3_matrix *)malloc(sites_on_node*sizeof(su3_matrix));
    if(s_link == NULL){
      fprintf(stderr,"hybrid_loop1: CAN'T MALLOC s_link\n");
      fflush(stderr);
      terminate(1);
    }

    /* Allocate space for space-link product */
    s_link_f = 
      (su3_matrix *)malloc(sites_on_node*sizeof(su3_matrix));
    if(s_link_f == NULL){
      fprintf(stderr,"hybrid_loop1: CAN'T MALLOC s_link_f\n");
      fflush(stderr);
      terminate(1);
    }

    /* Allocate space for time-link */ /* NOT NEEDED */
    t_link_f = 
      (su3_matrix *)malloc(sites_on_node*sizeof(su3_matrix));
    if(t_link_f == NULL){
      fprintf(stderr,"hybrid_loop1: CAN'T MALLOC t_link_f\n");
      fflush(stderr);
      terminate(1);
    }

    for(dir1=XUP;dir1<=YUP;dir1++) for(dir2=dir1+1;dir2<=ZUP;dir2++){

	/* Do off-axis "sqrt(2)" loops in (dir1,dir2)-plane */

	/* First construct the "diagonal" link in the (dir1,dir2) direction */

	mtag[dir1] = start_gather_site( F_OFFSET(link[dir2]), sizeof(su3_matrix),
	    dir1, EVENANDODD, gen_pt[dir1] );
	mtag[dir2] = start_gather_site( F_OFFSET(link[dir1]), sizeof(su3_matrix),
	    dir2, EVENANDODD, gen_pt[dir2] );

	FORALLSITES(i,s){
	    su3mat_copy( &(s->link[TUP]), t_link_f+i);
	}

	/* Concurrently start gather of time-like links across the
	   diagonal. */
	for(i=XUP;i<=TUP;i++)disp[i]=0;
	disp[dir1] = 1;
	disp[dir2] = 1;
	gmtag = start_general_gather_field(
		       (void *)t_link_f, sizeof(su3_matrix),
		       disp, EVENANDODD, gen_pt[4] );

	wait_gather( mtag[dir1]);
	wait_gather( mtag[dir2]);
	ftmp = 0.5;
	FORALLSITES(i,s){
	    mult_su3_nn( &(s->link[dir1]), (su3_matrix *)(gen_pt[dir1][i]),
		&(s->diag));
	    mult_su3_nn( &(s->link[dir2]), (su3_matrix *)(gen_pt[dir2][i]),
		&tmat1);
	    add_su3_matrix( &(s->diag), &tmat1, &(s->diag));
	    scalar_mult_su3_matrix( &(s->diag), ftmp, &(s->diag));
	}
	cleanup_gather( mtag[dir1]);
	cleanup_gather( mtag[dir2]);


	/* Recursively construct the space-like segments and compute
	   the Wilson loops with that segment */

	for(r=0;r<nxh;r++){

	    if( r==0 ){
		FORALLSITES(i,s){
		    su3mat_copy( &(s->diag), s_link+i);
		}
	    }
	    else{
		wait_general_gather( gmtag);
		FORALLSITES(i,s){
		    su3mat_copy( (su3_matrix *)(gen_pt[4][i]), &(s->staple));
		}
		cleanup_general_gather( gmtag);

		/* Inbetween gather time-like links across the diagonal. */
		gmtag = start_general_gather_field(
		       (void *)t_link_f,
		       sizeof(su3_matrix), disp, EVENANDODD, gen_pt[4] );

		FORALLSITES(i,s){
		    mult_su3_nn( &(s->diag), &(s->staple), s_link+i);
		}
	    }

	    FORALLSITES(i,s){
		su3mat_copy( s_link+i, s_link_f+i);
	    }

	    /* Start gather of forward space-like segments */
	    mtag[TUP] = start_gather_field( 
		       (void *)s_link_f, sizeof(su3_matrix),
		       TUP, EVENANDODD, gen_pt[TUP] );

	    /* Collect forward time-like links. */
	    wait_general_gather( gmtag);
	    FORALLSITES(i,s){
		su3mat_copy( (su3_matrix *)(gen_pt[4][i]), &(s->staple));
	    }
	    FORALLSITES(i,s){
		su3mat_copy( &(s->staple), t_link_f+i);
	    }
	    cleanup_general_gather( gmtag);

	    /* Inbetween gather space-links across the diagonal for next r. */
	    if( r<(nxh-1) ){
		gmtag = start_general_gather_field( 
 	         (void *)s_link,
	         sizeof(su3_matrix), disp, EVENANDODD, gen_pt[4] );
	    }


	    /* Recursively compute the Wilson loops of different time extent */
	    for(t=0;t<nth;t++){

		/* Collect forward space-like segments */
		wait_gather( mtag[TUP]);
		FORALLSITES(i,s){
		    su3mat_copy( (su3_matrix *)(gen_pt[TUP][i]), &(s->staple));
		}
		FORALLSITES(i,s){
		    su3mat_copy( &(s->staple), s_link_f+i);
		}

		/* Start gather for next t, if still needed. */
		if( t<(nth-1) ){
		    restart_gather_field( 
		     (void *)s_link_f, sizeof(su3_matrix),
		     TUP, EVENANDODD, gen_pt[TUP], mtag[TUP] );
		}
		else{
		    cleanup_gather( mtag[TUP]);
		}

		/* Finally, compute the Wilson loops. */
		FORALLSITES(i,s){
		    if( ((s->t)+t+1)>=nt ){
			mult_su3_nn( &(s->link[TUP]), s_link_f+i, &tmat1);
			mult_su3_na( &tmat1, t_link_f+i, &tmat2);
			wils_loop2[r+nrmax*t] +=
			    realtrace_su3( &tmat2, s_link+i);
		    }
		    else
		    {
			wils_loop2[r+nrmax*t] +=
			    realtrace_su3( s_link_f+i, s_link+i);
		    }
		}

	    } /* end loop over t */

	} /* end loop over r */

	/* Do off-axis "sqrt(2)" loops in (dir1,-dir2)-plane */

	/* First construct the "diagonal" link in the (dir1,-dir2) direction */

	/* Gather dir2-link from across the diagonal */
	for(i=XUP;i<=TUP;i++)disp[i]=0;
	disp[dir1] = 1;
	disp[dir2] = -1;
	gmtag = start_general_gather_site( F_OFFSET(link[dir2]), sizeof(su3_matrix),
	    disp, EVENANDODD, gen_pt[4] );

	/* Multiply one corner and then gather it */
	FORALLSITES(i,s){
	    mult_su3_an( &(s->link[dir2]), &(s->link[dir1]), &(s->staple));
	}
	mtag[dir2] = start_gather_site( F_OFFSET(staple), sizeof(su3_matrix),
	    OPP_DIR(dir2), EVENANDODD, gen_pt[dir2] );

	FORALLSITES(i,s){
	    su3mat_copy( &(s->link[TUP]), t_link_f+i);
	}

	/* Make second corner */
	wait_general_gather( gmtag);
	FORALLSITES(i,s){
	    mult_su3_na( &(s->link[dir1]), (su3_matrix *)(gen_pt[4][i]),
		&(s->diag));
	}
	cleanup_general_gather( gmtag);

	/* Collect first corner and add to second */
	wait_gather( mtag[dir2]);
	ftmp = 0.5;
	FORALLSITES(i,s){
	    add_su3_matrix( &(s->diag), (su3_matrix *)(gen_pt[dir2][i]),
		&(s->diag));
	    scalar_mult_su3_matrix( &(s->diag), ftmp, &(s->diag));
	}
	cleanup_gather( mtag[dir2]);

	/* Start gather of time-like links across the diagonal. */
	gmtag = start_general_gather_field( 
	       (void *)t_link_f, sizeof(su3_matrix),
	       disp, EVENANDODD, gen_pt[4] );


	/* Recursively construct the space-like segments and compute
	   the Wilson loops with that segment */

	for(r=0;r<nxh;r++){

	    if( r==0 ){
		FORALLSITES(i,s){
		    su3mat_copy( &(s->diag), s_link+i);
		}
	    }
	    else{
		wait_general_gather( gmtag);
		FORALLSITES(i,s){
		    su3mat_copy( (su3_matrix *)(gen_pt[4][i]), &(s->staple));
		}
		cleanup_general_gather( gmtag);

		/* Inbetween gather time-like links across the diagonal. */
		gmtag = start_general_gather_field(
		       (void *)t_link_f,
		       sizeof(su3_matrix), disp, EVENANDODD, gen_pt[4] );

		FORALLSITES(i,s){
		    mult_su3_nn( &(s->diag), &(s->staple), s_link+i);
		}
	    }

	    FORALLSITES(i,s){
		su3mat_copy( s_link+i, s_link_f+i);
	    }

	    /* Start gather of forward space-like segments */
	    mtag[TUP] = start_gather_field( 
		       (void *)s_link_f, sizeof(su3_matrix),
		       TUP, EVENANDODD, gen_pt[TUP] );

	    /* Collect forward time-like links. */
	    wait_general_gather( gmtag);
	    FORALLSITES(i,s){
		su3mat_copy( (su3_matrix *)(gen_pt[4][i]), &(s->staple));
	    }
	    FORALLSITES(i,s){
		su3mat_copy( &(s->staple), t_link_f+i);
	    }
	    cleanup_general_gather( gmtag);

	    /* Inbetween gather space-links across the diagonal for next r. */
	    if( r<(nxh-1) ){
		gmtag = start_general_gather_field( 
		       (void *)s_link,
		       sizeof(su3_matrix), disp, EVENANDODD, gen_pt[4] );
	    }


	    /* Recursively compute the Wilson loops of different time extent */
	    for(t=0;t<nth;t++){

		/* Collect forward space-like segments */
		wait_gather( mtag[TUP]);
		FORALLSITES(i,s){
		    su3mat_copy( (su3_matrix *)(gen_pt[TUP][i]), &(s->staple));
		}
		FORALLSITES(i,s){
		    su3mat_copy( &(s->staple), s_link_f+i);
		}

		/* Start gather for next t, if still needed. */
		if( t<(nth-1) ){
		  restart_gather_field(
			   (void *)s_link_f, sizeof(su3_matrix),
			   TUP, EVENANDODD, gen_pt[TUP], mtag[TUP] );
		}
		else{
		    cleanup_gather( mtag[TUP]);
		}

		/* Finally, compute the Wilson loops. */
		FORALLSITES(i,s){
		    if( ((s->t)+t+1)>=nt ){
			mult_su3_nn( &(s->link[TUP]), s_link_f+i, &tmat1);
			mult_su3_na( &tmat1, t_link_f+i, &tmat2);
			wils_loop2[r+nrmax*t] +=
			    realtrace_su3( &tmat2, s_link+i);
		    }
		    else
		    {
			wils_loop2[r+nrmax*t] +=
			    realtrace_su3( s_link_f+i, s_link+i);
		    }
		}

	    } /* end loop over t */

	} /* end loop over r */

    } /* end loop over dir1 < dir2 */

    /* Off-set for next bunch of Wilson loops */
    r_off = nxh;

    for(dir1=XUP;dir1<=ZUP;dir1++) for(dir2=XUP;dir2<=ZUP;dir2++)
    if( dir1 != dir2){

	/* Do off-axis "sqrt(5)" loops in (dir1,dir2)-plane */

	/* First construct the "diagonal" link in the (2*dir1,dir2) direction */

	mtag[dir1] = start_gather_site( F_OFFSET(link[dir1]), sizeof(su3_matrix),
	    dir1, EVENANDODD, gen_pt[dir1] );

	/* Start gather of dir2-link from "2*dir1" */
	for(i=XUP;i<=TUP;i++)disp[i]=0;
	disp[dir1] = 2;
	gmtag = start_general_gather_site( F_OFFSET(link[dir2]), sizeof(su3_matrix),
	    disp, EVENANDODD, gen_pt[4] );

	FORALLSITES(i,s){
	    su3mat_copy( &(s->link[TUP]), t_link_f+i);
	}

	/* Make double links in dir1 direction */
	wait_gather( mtag[dir1]);
	FORALLSITES(i,s){
	    mult_su3_nn( &(s->link[dir1]), (su3_matrix *)(gen_pt[dir1][i]),
		s_link+i);
	}
	cleanup_gather( mtag[dir1]);

	/* Gather the double links from dir2 direction */
	mtag[dir2] = start_gather_field( 
		    (void *)s_link, sizeof(su3_matrix),
		    dir2, EVENANDODD, gen_pt[dir2] );

	/* Make first corner */
	wait_general_gather( gmtag);
	FORALLSITES(i,s){
	    mult_su3_nn( s_link+i, (su3_matrix *)(gen_pt[4][i]),
		&(s->diag));
	}
	cleanup_general_gather( gmtag);

	/* Make second corner and add to first */
	wait_gather( mtag[dir2]);
	ftmp = 0.5;
	FORALLSITES(i,s){
	    mult_su3_nn( &(s->link[dir2]), (su3_matrix *)(gen_pt[dir2][i]),
		&tmat1);
	    add_su3_matrix( &(s->diag), &tmat1, &(s->diag));
	    scalar_mult_su3_matrix( &(s->diag), ftmp, &(s->diag));
	}
	cleanup_gather( mtag[dir2]);

	/* Start gather of time-like links across the diagonal. */
	for(i=XUP;i<=TUP;i++)disp[i]=0;
	disp[dir1] = 2;
	disp[dir2] = 1;
	gmtag = start_general_gather_field(
		       (void *)t_link_f, sizeof(su3_matrix),
		       disp, EVENANDODD, gen_pt[4] );


	/* Recursively construct the space-like segments and compute
	   the Wilson loops with that segment */

	for(r=0;r<nxh/2;r++){

	    if( r==0 ){
		FORALLSITES(i,s){
		    su3mat_copy( &(s->diag), s_link+i);
		}
	    }
	    else{
		wait_general_gather( gmtag);
		FORALLSITES(i,s){
		    su3mat_copy( (su3_matrix *)(gen_pt[4][i]), &(s->staple));
		}
		cleanup_general_gather( gmtag);

		/* Inbetween gather time-like links across the diagonal. */
		gmtag = start_general_gather_field( 
		       (void *)t_link_f,
		       sizeof(su3_matrix), disp, EVENANDODD, gen_pt[4] );

		FORALLSITES(i,s){
		    mult_su3_nn( &(s->diag), &(s->staple), s_link+i);
		}
	    }

	    FORALLSITES(i,s){
		su3mat_copy( s_link+i, s_link_f+i);
	    }

	    /* Start gather of forward space-like segments */
	    mtag[TUP] = start_gather_field(
		       (void *)s_link_f, sizeof(su3_matrix),
		       TUP, EVENANDODD, gen_pt[TUP] );

	    /* Collect forward time-like links. */
	    wait_general_gather( gmtag);
	    FORALLSITES(i,s){
		su3mat_copy( (su3_matrix *)(gen_pt[4][i]), &(s->staple));
	    }
	    FORALLSITES(i,s){
		su3mat_copy( &(s->staple), t_link_f+i);
	    }
	    cleanup_general_gather( gmtag);

	    /* Inbetween gather space-links across the diagonal for next r. */
	    if( r<(nxh/2-1) ){
		gmtag = start_general_gather_field(
		       (void *)s_link,
		       sizeof(su3_matrix), disp, EVENANDODD, gen_pt[4] );
	    }


	    /* Recursively compute the Wilson loops of different time extent */
	    for(t=0;t<nth;t++){

		/* Collect forward space-like segments */
		wait_gather( mtag[TUP]);
		FORALLSITES(i,s){
		    su3mat_copy( (su3_matrix *)(gen_pt[TUP][i]), &(s->staple));
		}
		FORALLSITES(i,s){
		    su3mat_copy( &(s->staple), s_link_f+i);
		}

		/* Start gather for next t, if still needed. */
		if( t<(nth-1) ){
		    restart_gather_field( 
		      (void *)s_link_f, sizeof(su3_matrix),
		      TUP, EVENANDODD, gen_pt[TUP], mtag[TUP] );
		}
		else{
		    cleanup_gather( mtag[TUP]);
		}

		/* Finally, compute the Wilson loops. */
		FORALLSITES(i,s){
		    if( ((s->t)+t+1)>=nt ){
			mult_su3_nn( &(s->link[TUP]), s_link_f+i, &tmat1);
			mult_su3_na( &tmat1, t_link_f+i, &tmat2);
			wils_loop2[r+r_off+nrmax*t] +=
			    realtrace_su3( &tmat2, s_link+i);
		    }
		    else
		    {
			wils_loop2[r+r_off+nrmax*t] +=
			    realtrace_su3( s_link_f+i, s_link+i);
		    }
		}

	    } /* end loop over t */

	} /* end loop over r */

	/* Do off-axis "sqrt(5)" loops in (dir1,-dir2)-plane */

	/* First construct the "diagonal" link in the (2*dir1,-dir2) dir */

	mtag[dir1] = start_gather_site( F_OFFSET(link[dir1]), sizeof(su3_matrix),
	    dir1, EVENANDODD, gen_pt[dir1] );

	/* Gather dir2-link from across the diagonal */
	for(i=XUP;i<=TUP;i++)disp[i]=0;
	disp[dir1] = 2;
	disp[dir2] = -1;
	gmtag = start_general_gather_site( F_OFFSET(link[dir2]), sizeof(su3_matrix),
	    disp, EVENANDODD, gen_pt[4] );

	FORALLSITES(i,s){
	    su3mat_copy( &(s->link[TUP]), t_link_f+i);
	}

	/* Make double links in dir1 direction */
	wait_gather( mtag[dir1]);
	FORALLSITES(i,s){
	    mult_su3_nn( &(s->link[dir1]), (su3_matrix *)(gen_pt[dir1][i]),
		s_link+i);
	}
	cleanup_gather( mtag[dir1]);

	/* Make one corner and then gather it */
	FORALLSITES(i,s){
	    mult_su3_an( &(s->link[dir2]), s_link+i, &(s->staple));
	}
	mtag[dir2] = start_gather_site( F_OFFSET(staple), sizeof(su3_matrix),
	    OPP_DIR(dir2), EVENANDODD, gen_pt[dir2] );

	/* Make second corner */
	wait_general_gather( gmtag);
	FORALLSITES(i,s){
	    mult_su3_na( s_link+i, (su3_matrix *)(gen_pt[4][i]),
		&(s->diag));
	}
	cleanup_general_gather( gmtag);

	/* Collect first corner and add to second */
	wait_gather( mtag[dir2]);
	ftmp = 0.5;
	FORALLSITES(i,s){
	    add_su3_matrix( &(s->diag), (su3_matrix *)(gen_pt[dir2][i]),
		&(s->diag));
	    scalar_mult_su3_matrix( &(s->diag), ftmp, &(s->diag));
	}
	cleanup_gather( mtag[dir2]);

	/* Start gather of time-like links across the diagonal. */
	gmtag = start_general_gather_field(
	       (void *)t_link_f, sizeof(su3_matrix),
	       disp, EVENANDODD, gen_pt[4] );


	/* Recursively construct the space-like segments and compute
	   the Wilson loops with that segment */

	for(r=0;r<nxh/2;r++){

	    if( r==0 ){
		FORALLSITES(i,s){
		    su3mat_copy( &(s->diag), s_link+i);
		}
	    }
	    else{
		wait_general_gather( gmtag);
		FORALLSITES(i,s){
		    su3mat_copy( (su3_matrix *)(gen_pt[4][i]), &(s->staple));
		}
		cleanup_general_gather( gmtag);

		/* Inbetween gather time-like links across the diagonal. */
		gmtag = start_general_gather_field( 
		       (void *)t_link_f,
		       sizeof(su3_matrix), disp, EVENANDODD, gen_pt[4] );

		FORALLSITES(i,s){
		    mult_su3_nn( &(s->diag), &(s->staple), s_link+i);
		}
	    }

	    FORALLSITES(i,s){
		su3mat_copy( s_link+i, s_link_f+i);
	    }

	    /* Start gather of forward space-like segments */
	    mtag[TUP] = start_gather_field( 
		       (void *)s_link_f, sizeof(su3_matrix),
		       TUP, EVENANDODD, gen_pt[TUP] );

	    /* Collect forward time-like links. */
	    wait_general_gather( gmtag);
	    FORALLSITES(i,s){
		su3mat_copy( (su3_matrix *)(gen_pt[4][i]), &(s->staple));
	    }
	    FORALLSITES(i,s){
		su3mat_copy( &(s->staple), t_link_f+i);
	    }
	    cleanup_general_gather( gmtag);

	    /* Inbetween gather space-links across the diagonal for next r. */
	    if( r<(nxh/2-1) ){
		gmtag = start_general_gather_field(
		       (void *)s_link,
		       sizeof(su3_matrix), disp, EVENANDODD, gen_pt[4] );
	    }


	    /* Recursively compute the Wilson loops of different time extent */
	    for(t=0;t<nth;t++){

		/* Collect forward space-like segments */
		wait_gather( mtag[TUP]);
		FORALLSITES(i,s){
		    su3mat_copy( (su3_matrix *)(gen_pt[TUP][i]), &(s->staple));
		}
		FORALLSITES(i,s){
		    su3mat_copy( &(s->staple), s_link_f+i);
		}

		/* Start gather for next t, if still needed. */
		if( t<(nth-1) ){
		    restart_gather_field(
			     (void *)s_link_f, sizeof(su3_matrix),
			     TUP, EVENANDODD, gen_pt[TUP], mtag[TUP] );
		}
		else{
		    cleanup_gather( mtag[TUP]);
		}

		/* Finally, compute the Wilson loops. */
		FORALLSITES(i,s){
		    if( ((s->t)+t+1)>=nt ){
			mult_su3_nn( &(s->link[TUP]), s_link_f+i, &tmat1);
			mult_su3_na( &tmat1, t_link_f+i, &tmat2);
			wils_loop2[r+r_off+nrmax*t] +=
			    realtrace_su3( &tmat2, s_link+i);
		    }
		    else
		    {
			wils_loop2[r+r_off+nrmax*t] +=
			    realtrace_su3( s_link_f+i, s_link+i);
		    }
		}

	    } /* end loop over t */

	} /* end loop over r */

    } /* end loop over dir1 != dir2 */

    /* Off-set for next bunch of Wilson loops */
    r_off = nxh + nxh/2;

    dir1 = XUP; dir2 = YUP; dir3 = ZUP;

	/* Do off-axis "sqrt(3)" loops in (dir1,dir2,dir3)-space */

	/* First construct the "body diagonal" link in the (dir1,dir2,dir3)
	   direction */

	/* Gather for first "plaquette" */
	mtag[0] = start_gather_site( F_OFFSET(link[dir2]), sizeof(su3_matrix),
	    dir1, EVENANDODD, gen_pt[0] );
	mtag[1] = start_gather_site( F_OFFSET(link[dir1]), sizeof(su3_matrix),
	    dir2, EVENANDODD, gen_pt[1] );

	/* Gather for second "plaquette" */
	mtag[7] = start_gather_site( F_OFFSET(link[dir3]), sizeof(su3_matrix),
	    dir1, EVENANDODD, gen_pt[7] );
	mtag[5] = start_gather_site( F_OFFSET(link[dir1]), sizeof(su3_matrix),
	    dir3, EVENANDODD, gen_pt[5] );

	FORALLSITES(i,s){
	    su3mat_copy( &(s->link[TUP]), t_link_f+i);
	}

	/* Make diagonal link in (dir1,dir2) direction and gather it */
	wait_gather( mtag[0]);
	wait_gather( mtag[1]);
	FORALLSITES(i,s){
	    mult_su3_nn( &(s->link[dir1]), (su3_matrix *)(gen_pt[0][i]),
		s_link+i);
	    mult_su3_nn( &(s->link[dir2]), (su3_matrix *)(gen_pt[1][i]),
		&tmat1);
	    add_su3_matrix( s_link+i, &tmat1, s_link+i);
	}
	cleanup_gather( mtag[0]);
	cleanup_gather( mtag[1]);
	mtag[2] = start_gather_field(s_link, sizeof(su3_matrix),
	    dir3, EVENANDODD, gen_pt[2] );

	/* Gather for third "plaquette" */
	mtag[1] = start_gather_site( F_OFFSET(link[dir3]), sizeof(su3_matrix),
	    dir2, EVENANDODD, gen_pt[1] );
	mtag[3] = start_gather_site( F_OFFSET(link[dir2]), sizeof(su3_matrix),
	    dir3, EVENANDODD, gen_pt[3] );

	/* Make diagonal link in (dir1,dir3) direction and gather it */
	wait_gather( mtag[7]);
	wait_gather( mtag[5]);
	FORALLSITES(i,s){
	    mult_su3_nn( &(s->link[dir1]), (su3_matrix *)(gen_pt[7][i]),
		s_link_f+i);
	    mult_su3_nn( &(s->link[dir3]), (su3_matrix *)(gen_pt[5][i]),
		&tmat1);
	    add_su3_matrix( s_link_f+i, &tmat1, s_link_f+i);
	}
	cleanup_gather( mtag[7]);
	cleanup_gather( mtag[5]);
	mtag[6] = start_gather_field(s_link_f, sizeof(su3_matrix),
	    dir2, EVENANDODD, gen_pt[6] );

	/* Make first body diagonal */
	wait_gather( mtag[2]);
	FORALLSITES(i,s){
	    mult_su3_nn( &(s->link[dir3]), (su3_matrix *)(gen_pt[2][i]),
		&(s->diag));
	}
	cleanup_gather( mtag[2]);

	/* Make diagonal link in (dir2,dir3) direction and gather it */
	wait_gather( mtag[1]);
	wait_gather( mtag[3]);
	FORALLSITES(i,s){
	    mult_su3_nn( &(s->link[dir2]), (su3_matrix *)(gen_pt[1][i]),
		s_link+i);
	    mult_su3_nn( &(s->link[dir3]), (su3_matrix *)(gen_pt[3][i]),
		&tmat1);
	    add_su3_matrix( s_link+i, &tmat1, s_link+i);
	}
	cleanup_gather( mtag[1]);
	cleanup_gather( mtag[3]);
	mtag[0] = start_gather_field( s_link, sizeof(su3_matrix),
	    dir1, EVENANDODD, gen_pt[0] );

	/* Make second body diagonal and add to first */
	wait_gather( mtag[6]);
	FORALLSITES(i,s){
	    mult_su3_nn( &(s->link[dir2]), (su3_matrix *)(gen_pt[6][i]),
		&tmat1);
	    add_su3_matrix( &(s->diag), &tmat1, &(s->diag));
	}
	cleanup_gather( mtag[6]);

	/* Make third body diagonal and add */
	wait_gather( mtag[0]);
	ftmp = 1.0 / 6.0;
	FORALLSITES(i,s){
	    mult_su3_nn( &(s->link[dir1]), (su3_matrix *)(gen_pt[0][i]),
		&tmat1);
	    add_su3_matrix( &(s->diag), &tmat1, &(s->diag));
	    scalar_mult_su3_matrix( &(s->diag), ftmp, &(s->diag));
	}
	cleanup_gather( mtag[0]);

	/* Start gather of time-like links across the body diagonal. */
	disp[dir1] = 1;
	disp[dir2] = 1;
	disp[dir3] = 1;
	disp[TUP] = 0;
	gmtag = start_general_gather_field(
	       (void *)t_link_f, sizeof(su3_matrix),
	       disp, EVENANDODD, gen_pt[4] );


	/* Recursively construct the space-like segments and compute
	   the Wilson loops with that segment */

	for(r=0;r<nxh;r++){

	    if( r==0 ){
		FORALLSITES(i,s){
		    su3mat_copy( &(s->diag), s_link+i);
		}
	    }
	    else{
		wait_general_gather( gmtag);
		FORALLSITES(i,s){
		    su3mat_copy( (su3_matrix *)(gen_pt[4][i]), &(s->staple));
		}
		cleanup_general_gather( gmtag);

		/* Inbetween gather time-like links across the diagonal. */
		gmtag = start_general_gather_field( 
		       (void *)t_link_f,
		       sizeof(su3_matrix), disp, EVENANDODD, gen_pt[4] );

		FORALLSITES(i,s){
		    mult_su3_nn( &(s->diag), &(s->staple), s_link+i);
		}
	    }

	    FORALLSITES(i,s){
		su3mat_copy( s_link+i, s_link_f+i);
	    }

	    /* Start gather of forward space-like segments */
	    mtag[TUP] = start_gather_field(
		       (void *)s_link_f, sizeof(su3_matrix),
		       TUP, EVENANDODD, gen_pt[TUP] );

	    /* Collect forward time-like links. */
	    wait_general_gather( gmtag);
	    FORALLSITES(i,s){
		su3mat_copy( (su3_matrix *)(gen_pt[4][i]), &(s->staple));
	    }
	    FORALLSITES(i,s){
		su3mat_copy( &(s->staple), t_link_f+i);
	    }
	    cleanup_general_gather( gmtag);

	    /* Inbetween gather space-links across the diagonal for next r. */
	    if( r<(nxh-1) ){
		gmtag = start_general_gather_field(
		       (void *)s_link,
		       sizeof(su3_matrix), disp, EVENANDODD, gen_pt[4] );
	    }


	    /* Recursively compute the Wilson loops of different time extent */
	    for(t=0;t<nth;t++){

		/* Collect forward space-like segments */
		wait_gather( mtag[TUP]);
		FORALLSITES(i,s){
		    su3mat_copy( (su3_matrix *)(gen_pt[TUP][i]), &(s->staple));
		}
		FORALLSITES(i,s){
		    su3mat_copy( &(s->staple), s_link_f+i);
		}

		/* Start gather for next t, if still needed. */
		if( t<(nth-1) ){
		    restart_gather_field( 
		     (void *)s_link_f, sizeof(su3_matrix),
		     TUP, EVENANDODD, gen_pt[TUP], mtag[TUP] );
		}
		else{
		    cleanup_gather( mtag[TUP]);
		}

		/* Finally, compute the Wilson loops. */
		FORALLSITES(i,s){
		    if( ((s->t)+t+1)>=nt ){
			mult_su3_nn( &(s->link[TUP]), s_link_f+i, &tmat1);
			mult_su3_na( &tmat1, t_link_f+i, &tmat2);
			wils_loop2[r+r_off+nrmax*t] +=
			    realtrace_su3( &tmat2, s_link+i);
		    }
		    else
		    {
			wils_loop2[r+r_off+nrmax*t] +=
			    realtrace_su3( s_link_f+i, s_link+i);
		    }
		}

	    } /* end loop over t */

	} /* end loop over r */

	/* Do off-axis "sqrt(3)" loops in (dir1,dir2,-dir3)-space */

	/* First construct the "body diagonal" link in the (dir1,dir2,-dir3)
	   direction */

	/* Gather for first "plaquette" */
	mtag[0] = start_gather_site( F_OFFSET(link[dir2]), sizeof(su3_matrix),
	    dir1, EVENANDODD, gen_pt[0] );
	mtag[1] = start_gather_site( F_OFFSET(link[dir1]), sizeof(su3_matrix),
	    dir2, EVENANDODD, gen_pt[1] );

	/* Start one corner for second "plaquette" and gather */
	FORALLSITES(i,s){
	    mult_su3_an( &(s->link[dir3]), &(s->link[dir1]), &(s->staple));
	}
	for(i=XUP;i<=TUP;i++)disp[i]=0;
	disp[dir1] = 1;
	disp[dir3] = -1;
	gmtag = start_general_gather_site( F_OFFSET(link[dir3]), sizeof(su3_matrix),
	    disp, EVENANDODD, gen_pt[4] );
	mtag[5] = start_gather_site( F_OFFSET(staple), sizeof(su3_matrix),
	    OPP_DIR(dir3), EVENANDODD, gen_pt[5] );

	/* Make diagonal link in (dir1,dir2) direction, multiply to get
	   first body diagonal and gather it */
	wait_gather( mtag[0]);
	wait_gather( mtag[1]);
	FORALLSITES(i,s){
	    mult_su3_nn( &(s->link[dir1]), (su3_matrix *)(gen_pt[0][i]),
		&tmat1);
	    mult_su3_nn( &(s->link[dir2]), (su3_matrix *)(gen_pt[1][i]),
		&tmat2);
	    add_su3_matrix( &tmat1, &tmat2, &tmat1);
	    mult_su3_an( &(s->link[dir3]), &tmat1, s_link+i);
	}
	cleanup_gather( mtag[0]);
	cleanup_gather( mtag[1]);
	mtag[2] = start_gather_field(
		 (void *)s_link, sizeof(su3_matrix),
		 OPP_DIR(dir3), EVENANDODD, gen_pt[2] );

	/* Make diagonal link in (dir1,-dir3) direction and gather it */
	wait_general_gather( gmtag);
	wait_gather( mtag[5]);
	FORALLSITES(i,s){
	    mult_su3_na( &(s->link[dir1]), (su3_matrix *)(gen_pt[4][i]),
		s_link_f+i);
	    add_su3_matrix( s_link_f+i, (su3_matrix *)(gen_pt[5][i]),
		s_link_f+i);
	}
	cleanup_general_gather( gmtag);
	cleanup_gather( mtag[5]);
	mtag[6] = start_gather_field(
		 (void *)s_link_f, sizeof(su3_matrix),
		 dir2, EVENANDODD, gen_pt[6] );

	/* Start one corner for third "plaquette" and gather */
	FORALLSITES(i,s){
	    mult_su3_an( &(s->link[dir3]), &(s->link[dir2]), &(s->diag));
	}
	for(i=XUP;i<=TUP;i++)disp[i]=0;
	disp[dir2] = 1;
	disp[dir3] = -1;
	gmtag = start_general_gather_site( F_OFFSET(link[dir3]), sizeof(su3_matrix),
	    disp, EVENANDODD, gen_pt[4] );
	mtag[3] = start_gather_site( F_OFFSET(diag), sizeof(su3_matrix),
	    OPP_DIR(dir3), EVENANDODD, gen_pt[3] );

	FORALLSITES(i,s){
	    su3mat_copy( &(s->link[TUP]), t_link_f+i);
	}

	/* Make diagonal link in (dir2,-dir3) direction and gather it */
	wait_general_gather( gmtag);
	wait_gather( mtag[3]);
	FORALLSITES(i,s){
	    mult_su3_na( &(s->link[dir2]), (su3_matrix *)(gen_pt[4][i]),
		&(s->staple));
	    add_su3_matrix( &(s->staple), (su3_matrix *)(gen_pt[3][i]),
		&(s->staple));
	}
	cleanup_general_gather( gmtag);
	cleanup_gather( mtag[3]);
	mtag[0] = start_gather_site( F_OFFSET(staple), sizeof(su3_matrix),
	    dir1, EVENANDODD, gen_pt[0] );

	/* Make second body diagonal and add to it gathered first */
	wait_gather( mtag[2]);
	wait_gather( mtag[6]);
	FORALLSITES(i,s){
	    mult_su3_nn( &(s->link[dir2]), (su3_matrix *)(gen_pt[6][i]),
		&(s->diag));
	    add_su3_matrix( &(s->diag), (su3_matrix *)(gen_pt[2][i]),
		&(s->diag));
	}
	cleanup_gather( mtag[2]);
	cleanup_gather( mtag[6]);

	/* Make third body diagonal and add */
	wait_gather( mtag[0]);
	ftmp = 1.0 / 6.0;
	FORALLSITES(i,s){
	    mult_su3_nn( &(s->link[dir1]), (su3_matrix *)(gen_pt[0][i]),
		&tmat1);
	    add_su3_matrix( &(s->diag), &tmat1, &(s->diag));
	    scalar_mult_su3_matrix( &(s->diag), ftmp, &(s->diag));
	}
	cleanup_gather( mtag[0]);

	/* Start gather of time-like links across the body diagonal. */
	disp[dir1] = 1;
	disp[dir2] = 1;
	disp[dir3] = -1;
	disp[TUP] = 0;
	gmtag = start_general_gather_field(
		       (void *)t_link_f, sizeof(su3_matrix),
		       disp, EVENANDODD, gen_pt[4] );


	/* Recursively construct the space-like segments and compute
	   the Wilson loops with that segment */

	for(r=0;r<nxh;r++){

	    if( r==0 ){
		FORALLSITES(i,s){
		    su3mat_copy( &(s->diag), s_link+i);
		}
	    }
	    else{
		wait_general_gather( gmtag);
		FORALLSITES(i,s){
		    su3mat_copy( (su3_matrix *)(gen_pt[4][i]), &(s->staple));
		}
		cleanup_general_gather( gmtag);

		/* Inbetween gather time-like links across the diagonal. */
		gmtag = start_general_gather_field( 
		       (void *)t_link_f,
		       sizeof(su3_matrix), disp, EVENANDODD, gen_pt[4] );

		FORALLSITES(i,s){
		    mult_su3_nn( &(s->diag), &(s->staple), s_link+i);
		}
	    }

	    FORALLSITES(i,s){
		su3mat_copy( s_link+i, s_link_f+i);
	    }

	    /* Start gather of forward space-like segments */
	    mtag[TUP] = start_gather_field(
		       (void *)s_link_f, sizeof(su3_matrix),
		       TUP, EVENANDODD, gen_pt[TUP] );

	    /* Collect forward time-like links. */
	    wait_general_gather( gmtag);
	    FORALLSITES(i,s){
		su3mat_copy( (su3_matrix *)(gen_pt[4][i]), &(s->staple));
	    }
	    FORALLSITES(i,s){
		su3mat_copy( &(s->staple), t_link_f+i);
	    }
	    cleanup_general_gather( gmtag);

	    /* Inbetween gather space-links across the diagonal for next r. */
	    if( r<(nxh-1) ){
		gmtag = start_general_gather_field(
		       (void *)s_link,
		       sizeof(su3_matrix), disp, EVENANDODD, gen_pt[4] );
	    }


	    /* Recursively compute the Wilson loops of different time extent */
	    for(t=0;t<nth;t++){

		/* Collect forward space-like segments */
		wait_gather( mtag[TUP]);
		FORALLSITES(i,s){
		    su3mat_copy( (su3_matrix *)(gen_pt[TUP][i]), &(s->staple));
		}
		FORALLSITES(i,s){
		    su3mat_copy( &(s->staple), s_link_f+i);
		}

		/* Start gather for next t, if still needed. */
		if( t<(nth-1) ){
		    restart_gather_field(
			     (void *)s_link_f, sizeof(su3_matrix),
			     TUP, EVENANDODD, gen_pt[TUP], mtag[TUP] );
		}
		else{
		    cleanup_gather( mtag[TUP]);
		}

		/* Finally, compute the Wilson loops. */
		FORALLSITES(i,s){
		    if( ((s->t)+t+1)>=nt ){
			mult_su3_nn( &(s->link[TUP]), s_link_f+i, &tmat1);
			mult_su3_na( &tmat1, t_link_f+i, &tmat2);
			wils_loop2[r+r_off+nrmax*t] +=
			    realtrace_su3( &tmat2, s_link+i);
		    }
		    else
		    {
			wils_loop2[r+r_off+nrmax*t] +=
			    realtrace_su3( s_link_f+i, s_link+i);
		    }
		}

	    } /* end loop over t */

	} /* end loop over r */

	/* Do off-axis "sqrt(3)" loops in (dir1,-dir2,dir3)-space */

	/* First construct the "body diagonal" link in the (dir1,-dir2,dir3)
	   direction */

	/* Gather for first "plaquette" */
	mtag[0] = start_gather_site( F_OFFSET(link[dir3]), sizeof(su3_matrix),
	    dir1, EVENANDODD, gen_pt[0] );
	mtag[2] = start_gather_site( F_OFFSET(link[dir1]), sizeof(su3_matrix),
	    dir3, EVENANDODD, gen_pt[2] );

	/* Start one corner for second "plaquette" and gather */
	FORALLSITES(i,s){
	    mult_su3_an( &(s->link[dir2]), &(s->link[dir1]), &(s->staple));
	}
	for(i=XUP;i<=TUP;i++)disp[i]=0;
	disp[dir1] = 1;
	disp[dir2] = -1;
	gmtag = start_general_gather_site( F_OFFSET(link[dir2]), sizeof(su3_matrix),
	    disp, EVENANDODD, gen_pt[4] );
	mtag[6] = start_gather_site( F_OFFSET(staple), sizeof(su3_matrix),
	    OPP_DIR(dir2), EVENANDODD, gen_pt[6] );

	/* Make diagonal link in (dir1,dir3) direction, multiply to get
	   first body diagonal and gather it */
	wait_gather( mtag[0]);
	wait_gather( mtag[2]);
	FORALLSITES(i,s){
	    mult_su3_nn( &(s->link[dir1]), (su3_matrix *)(gen_pt[0][i]),
		&tmat1);
	    mult_su3_nn( &(s->link[dir3]), (su3_matrix *)(gen_pt[2][i]),
		&tmat2);
	    add_su3_matrix( &tmat1, &tmat2, &tmat1);
	    mult_su3_an( &(s->link[dir2]), &tmat1, s_link+i);
	}
	cleanup_gather( mtag[0]);
	cleanup_gather( mtag[2]);
	mtag[1] = start_gather_field(
		 (void *)s_link, sizeof(su3_matrix),
		 OPP_DIR(dir2), EVENANDODD, gen_pt[1] );

	/* Make diagonal link in (dir1,-dir2) direction and gather it */
	wait_general_gather( gmtag);
	wait_gather( mtag[6]);
	FORALLSITES(i,s){
	    mult_su3_na( &(s->link[dir1]), (su3_matrix *)(gen_pt[4][i]),
		s_link_f+i);
	    add_su3_matrix( s_link_f+i, (su3_matrix *)(gen_pt[6][i]),
		s_link_f+i);
	}
	cleanup_general_gather( gmtag);
	cleanup_gather( mtag[6]);
	mtag[5] = start_gather_field( 
		 (void *)s_link_f, sizeof(su3_matrix),
		 dir3, EVENANDODD, gen_pt[5] );

	/* Start one corner for third "plaquette" and gather */
	FORALLSITES(i,s){
	    mult_su3_an( &(s->link[dir2]), &(s->link[dir3]), &(s->diag));
	}
	for(i=XUP;i<=TUP;i++)disp[i]=0;
	disp[dir2] = -1;
	disp[dir3] = 1;
	gmtag = start_general_gather_site( F_OFFSET(link[dir2]), sizeof(su3_matrix),
	    disp, EVENANDODD, gen_pt[4] );
	mtag[3] = start_gather_site( F_OFFSET(diag), sizeof(su3_matrix),
	    OPP_DIR(dir2), EVENANDODD, gen_pt[3] );

	FORALLSITES(i,s){
	    su3mat_copy( &(s->link[TUP]), t_link_f+i);
	}

	/* Make diagonal link in (-dir2,dir3) direction and gather it */
	wait_general_gather( gmtag);
	wait_gather( mtag[3]);
	FORALLSITES(i,s){
	    mult_su3_na( &(s->link[dir3]), (su3_matrix *)(gen_pt[4][i]),
		&(s->staple));
	    add_su3_matrix( &(s->staple), (su3_matrix *)(gen_pt[3][i]),
		&(s->staple));
	}
	cleanup_general_gather( gmtag);
	cleanup_gather( mtag[3]);
	mtag[0] = start_gather_site( F_OFFSET(staple), sizeof(su3_matrix),
	    dir1, EVENANDODD, gen_pt[0] );

	/* Make second body diagonal and add to it gathered first */
	wait_gather( mtag[1]);
	wait_gather( mtag[5]);
	FORALLSITES(i,s){
	    mult_su3_nn( &(s->link[dir3]), (su3_matrix *)(gen_pt[5][i]),
		&(s->diag));
	    add_su3_matrix( &(s->diag), (su3_matrix *)(gen_pt[1][i]),
		&(s->diag));
	}
	cleanup_gather( mtag[1]);
	cleanup_gather( mtag[5]);

	/* Make third body diagonal and add */
	wait_gather( mtag[0]);
	ftmp = 1.0 / 6.0;
	FORALLSITES(i,s){
	    mult_su3_nn( &(s->link[dir1]), (su3_matrix *)(gen_pt[0][i]),
		&tmat1);
	    add_su3_matrix( &(s->diag), &tmat1, &(s->diag));
	    scalar_mult_su3_matrix( &(s->diag), ftmp, &(s->diag));
	}
	cleanup_gather( mtag[0]);

	/* Start gather of time-like links across the body diagonal. */
	disp[dir1] = 1;
	disp[dir2] = -1;
	disp[dir3] = 1;
	disp[TUP] = 0;
	gmtag = start_general_gather_field(
	       (void *)t_link_f, sizeof(su3_matrix),
	       disp, EVENANDODD, gen_pt[4] );


	/* Recursively construct the space-like segments and compute
	   the Wilson loops with that segment */

	for(r=0;r<nxh;r++){

	    if( r==0 ){
		FORALLSITES(i,s){
		    su3mat_copy( &(s->diag), s_link+i);
		}
	    }
	    else{
		wait_general_gather( gmtag);
		FORALLSITES(i,s){
		    su3mat_copy( (su3_matrix *)(gen_pt[4][i]), &(s->staple));
		}
		cleanup_general_gather( gmtag);

		/* Inbetween gather time-like links across the diagonal. */
		gmtag = start_general_gather_field(
		       (void *)t_link_f,
		       sizeof(su3_matrix), disp, EVENANDODD, gen_pt[4] );

		FORALLSITES(i,s){
		    mult_su3_nn( &(s->diag), &(s->staple), s_link+i);
		}
	    }

	    FORALLSITES(i,s){
		su3mat_copy( s_link+i, s_link_f+i);
	    }

	    /* Start gather of forward space-like segments */
	    mtag[TUP] = start_gather_field(
	        (void *)s_link_f, sizeof(su3_matrix),
		TUP, EVENANDODD, gen_pt[TUP] );

	    /* Collect forward time-like links. */
	    wait_general_gather( gmtag);
	    FORALLSITES(i,s){
		su3mat_copy( (su3_matrix *)(gen_pt[4][i]), &(s->staple));
	    }
	    FORALLSITES(i,s){
		su3mat_copy( &(s->staple), t_link_f+i);
	    }
	    cleanup_general_gather( gmtag);

	    /* Inbetween gather space-links across the diagonal for next r. */
	    if( r<(nxh-1) ){
		gmtag = start_general_gather_field(
		       (void *)s_link,
		       sizeof(su3_matrix), disp, EVENANDODD, gen_pt[4] );
	    }


	    /* Recursively compute the Wilson loops of different time extent */
	    for(t=0;t<nth;t++){

		/* Collect forward space-like segments */
		wait_gather( mtag[TUP]);
		FORALLSITES(i,s){
		    su3mat_copy( (su3_matrix *)(gen_pt[TUP][i]), &(s->staple));
		}
		FORALLSITES(i,s){
		    su3mat_copy( &(s->staple), s_link_f+i);
		}

		/* Start gather for next t, if still needed. */
		if( t<(nth-1) ){
		    restart_gather_field(
		     (void *)s_link_f, sizeof(su3_matrix),
		     TUP, EVENANDODD, gen_pt[TUP], mtag[TUP] );
		}
		else{
		    cleanup_gather( mtag[TUP]);
		}

		/* Finally, compute the Wilson loops. */
		FORALLSITES(i,s){
		    if( ((s->t)+t+1)>=nt ){
			mult_su3_nn( &(s->link[TUP]), s_link_f+i, &tmat1);
			mult_su3_na( &tmat1, t_link_f+i, &tmat2);
			wils_loop2[r+r_off+nrmax*t] +=
			    realtrace_su3( &tmat2, s_link+i);
		    }
		    else
		    {
			wils_loop2[r+r_off+nrmax*t] +=
			    realtrace_su3( s_link_f+i, s_link+i);
		    }
		}

	    } /* end loop over t */

	} /* end loop over r */

	/* Do off-axis "sqrt(3)" loops in (dir1,-dir2,-dir3)-space */

	/* First construct the "body diagonal" link in the (dir1,-dir2,-dir3)
	   direction */

	/* Gather for first "plaquette" */
	mtag[1] = start_gather_site( F_OFFSET(link[dir3]), sizeof(su3_matrix),
	    dir2, EVENANDODD, gen_pt[1] );
	mtag[2] = start_gather_site( F_OFFSET(link[dir2]), sizeof(su3_matrix),
	    dir3, EVENANDODD, gen_pt[2] );

	/* Start one corner for second "plaquette" and gather */
	FORALLSITES(i,s){
	    mult_su3_an( &(s->link[dir2]), &(s->link[dir1]), &(s->staple));
	}
	for(i=XUP;i<=TUP;i++)disp[i]=0;
	disp[dir1] = 1;
	disp[dir2] = -1;
	gmtag = start_general_gather_site( F_OFFSET(link[dir2]), sizeof(su3_matrix),
	    disp, EVENANDODD, gen_pt[4] );
	mtag[6] = start_gather_site( F_OFFSET(staple), sizeof(su3_matrix),
	    OPP_DIR(dir2), EVENANDODD, gen_pt[6] );

	/* Make diagonal link in (dir2,dir3) direction */
	wait_gather( mtag[1]);
	wait_gather( mtag[2]);
	FORALLSITES(i,s){
	    mult_su3_nn( &(s->link[dir2]), (su3_matrix *)(gen_pt[1][i]),
		s_link+i);
	    mult_su3_nn( &(s->link[dir3]), (su3_matrix *)(gen_pt[2][i]),
		&tmat1);
	    add_su3_matrix( s_link+i, &tmat1, s_link+i);
	}
	cleanup_gather( mtag[1]);
	cleanup_gather( mtag[2]);

	/* Make diagonal link in (dir1,-dir2) direction, multiply to get
	   second body diagonal and gather it */
	wait_general_gather( gmtag);
	wait_gather( mtag[6]);
	FORALLSITES(i,s){
	    mult_su3_na( &(s->link[dir1]), (su3_matrix *)(gen_pt[4][i]),
		&tmat1);
	    add_su3_matrix( &tmat1, (su3_matrix *)(gen_pt[6][i]), &tmat1);
	    mult_su3_an( &(s->link[dir3]), &tmat1, s_link_f+i);
	}
	cleanup_general_gather( gmtag);
	cleanup_gather( mtag[6]);
	mtag[2] = start_gather_field(
		 (void *)s_link_f, sizeof(su3_matrix),
		 OPP_DIR(dir3), EVENANDODD, gen_pt[2] );

	/* Start one corner for third "plaquette" and gather */
	FORALLSITES(i,s){
	    mult_su3_an( &(s->link[dir3]), &(s->link[dir1]), &(s->diag));
	}
	for(i=XUP;i<=TUP;i++)disp[i]=0;
	disp[dir1] = 1;
	disp[dir3] = -1;
	gmtag = start_general_gather_site( F_OFFSET(link[dir3]), sizeof(su3_matrix),
	    disp, EVENANDODD, gen_pt[4] );
	mtag[3] = start_gather_site( F_OFFSET(diag), sizeof(su3_matrix),
	    OPP_DIR(dir3), EVENANDODD, gen_pt[3] );

	FORALLSITES(i,s){
	    su3mat_copy( &(s->link[TUP]), t_link_f+i);
	}

	/* Make diagonal link in (dir1,-dir3) direction, multiply to get
	   third body diagonal and gather it */
	wait_general_gather( gmtag);
	wait_gather( mtag[3]);
	FORALLSITES(i,s){
	    mult_su3_na( &(s->link[dir1]), (su3_matrix *)(gen_pt[4][i]),
		&tmat1);
	    add_su3_matrix( &tmat1, (su3_matrix *)(gen_pt[3][i]), &tmat1);
	    mult_su3_an( &(s->link[dir2]), &tmat1, &(s->staple));
	}
	cleanup_general_gather( gmtag);
	cleanup_gather( mtag[3]);
	mtag[1] = start_gather_site( F_OFFSET(staple), sizeof(su3_matrix),
	    OPP_DIR(dir2), EVENANDODD, gen_pt[1] );

	/* Finally gather first "plaquette" */
	disp[dir1] = 1;
	disp[dir2] = -1;
	disp[dir3] = -1;
	disp[TUP] = 0;
	gmtag = start_general_gather_field(
	       (void *)s_link, sizeof(su3_matrix),
	       disp, EVENANDODD, gen_pt[4] );

	/* Collect second and third body diagonal and add them */
	wait_gather( mtag[1]);
	wait_gather( mtag[2]);
	FORALLSITES(i,s){
	    add_su3_matrix( (su3_matrix *)(gen_pt[1][i]),
		(su3_matrix *)(gen_pt[2][i]), &(s->diag));
	}
	cleanup_gather( mtag[1]);
	cleanup_gather( mtag[2]);

	/* Make first body diagonal and add */
	wait_general_gather( gmtag);
	ftmp = 1.0 / 6.0;
	FORALLSITES(i,s){
	    mult_su3_na( &(s->link[dir1]), (su3_matrix *)(gen_pt[4][i]),
		&tmat1);
	    add_su3_matrix( &(s->diag), &tmat1, &(s->diag));
	    scalar_mult_su3_matrix( &(s->diag), ftmp, &(s->diag));
	}
	cleanup_general_gather( gmtag);

	/* Start gather of time-like links across the body diagonal. */
	gmtag = start_general_gather_field(
	       (void *)t_link_f, sizeof(su3_matrix),
	       disp, EVENANDODD, gen_pt[4] );


	/* Recursively construct the space-like segments and compute
	   the Wilson loops with that segment */

	for(r=0;r<nxh;r++){

	    if( r==0 ){
		FORALLSITES(i,s){
		    su3mat_copy( &(s->diag), s_link+i);
		}
	    }
	    else{
		wait_general_gather( gmtag);
		FORALLSITES(i,s){
		    su3mat_copy( (su3_matrix *)(gen_pt[4][i]), &(s->staple));
		}
		cleanup_general_gather( gmtag);

		/* Inbetween gather time-like links across the diagonal. */
		gmtag = start_general_gather_field(
		       (void *)t_link_f,
		       sizeof(su3_matrix), disp, EVENANDODD, gen_pt[4] );

		FORALLSITES(i,s){
		    mult_su3_nn( &(s->diag), &(s->staple), s_link+i);
		}
	    }

	    FORALLSITES(i,s){
		su3mat_copy( s_link+i, s_link_f+i);
	    }

	    /* Start gather of forward space-like segments */
	    mtag[TUP] = start_gather_field(
	        (void *)s_link_f, sizeof(su3_matrix),
		TUP, EVENANDODD, gen_pt[TUP] );

	    /* Collect forward time-like links. */
	    wait_general_gather( gmtag);
	    FORALLSITES(i,s){
		su3mat_copy( (su3_matrix *)(gen_pt[4][i]), &(s->staple));
	    }
	    FORALLSITES(i,s){
		su3mat_copy( &(s->staple), t_link_f+i);
	    }
	    cleanup_general_gather( gmtag);

	    /* Inbetween gather space-links across the diagonal for next r. */
	    if( r<(nxh-1) ){
		gmtag = start_general_gather_field(
		       (void *)s_link,
		       sizeof(su3_matrix), disp, EVENANDODD, gen_pt[4] );
	    }


	    /* Recursively compute the Wilson loops of different time extent */
	    for(t=0;t<nth;t++){

		/* Collect forward space-like segments */
		wait_gather( mtag[TUP]);
		FORALLSITES(i,s){
		    su3mat_copy( (su3_matrix *)(gen_pt[TUP][i]), &(s->staple));
		}
		FORALLSITES(i,s){
		    su3mat_copy( &(s->staple), s_link_f+i);
		}

		/* Start gather for next t, if still needed. */
		if( t<(nth-1) ){
		    restart_gather_field(
			     (void *)s_link_f, sizeof(su3_matrix),
			     TUP, EVENANDODD, gen_pt[TUP], mtag[TUP] );
		}
		else{
		    cleanup_gather( mtag[TUP]);
		}

		/* Finally, compute the Wilson loops. */
		FORALLSITES(i,s){
		    if( ((s->t)+t+1)>=nt ){
			mult_su3_nn( &(s->link[TUP]), s_link_f+i, &tmat1);
			mult_su3_na( &tmat1, t_link_f+i, &tmat2);
			wils_loop2[r+r_off+nrmax*t] +=
			    realtrace_su3( &tmat2, s_link+i);
		    }
		    else
		    {
			wils_loop2[r+r_off+nrmax*t] +=
			    realtrace_su3( s_link_f+i, s_link+i);
		    }
		}

	    } /* end loop over t */

	} /* end loop over r */

    /* end of "sqrt(3)" Wilson loops with dir1 < dir2 < dir3 */

    /* Normalize and print the Wilson loops */
    for(t=0;t<nth;t++){
	for(r=0;r<nxh;r++){
	    ftmp = wils_loop2[r+nrmax*t];
	    g_floatsum( &ftmp);
	    ftmp /= (Real)(18*volume);
	    if(this_node == 0)printf("WILS_LOOP2_%d  %d  %d  %e\n",
		 tot_smear, r, t, (double)ftmp);
	}

	r_off = nxh;
	for(r=0;r<nxh/2;r++){
	    ftmp = wils_loop2[r+r_off+nrmax*t];
	    g_floatsum( &ftmp);
	    ftmp /= (Real)(36*volume);
	    if(this_node == 0)printf("WILS_LOOP2_%d  %d  %d  %e\n",
		 tot_smear, r+r_off, t, (double)ftmp);
	}

	r_off = nxh + nxh/2;
	for(r=0;r<nxh;r++){
	    ftmp = wils_loop2[r+r_off+nrmax*t];
	    g_floatsum( &ftmp);
	    ftmp /= (Real)(12*volume);
	    if(this_node == 0)printf("WILS_LOOP2_%d  %d  %d  %e\n",
		 tot_smear, r+r_off, t, (double)ftmp);
	}
    }

    free( wils_loop2);
    free(s_link);
    free(s_link_f);
    free(t_link_f);

} /* w_loop2 */

