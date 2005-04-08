/************************** w_loop2_red1.c *******************************/
/* MIMD version 6 */
/* This version uses gathers to get the neighbors */
/* version of 3/16/94 by UMH */
/* 05/24/99 Version 5 port CD */

/* Computes time-like, off-axis Wilson loops on gauge configuration
   in axial gauge with time-like gauge fields of last time slice in
   all other time slices as well, instead of the unit matrix! 

   
   05/24/99:
   NOTE: here ONLY sqrt(5) and sqrt(8) correlations calculated
         Also: nrmax = 2 + 2

*/

#include "string_break_includes.h"

void w_loop2(int tot_smear) {

register int i,dir1,dir2,dir3,r,t,r_off;
int nth,nxh,nrmax;
register site *s;
su3_matrix tmat1,tmat2;
msg_tag *mtag[8],*gmtag;
Real *wils_loop2,ftmp;
int disp[4];    /* displacement vector for general gather */

    if( nx != ny || nx != ny){
	if(this_node == 0)printf("w_loop2 gives wrong results for nx!=ny!=nz");
        return;
    }

    nth = nt/2;
    nxh = 2;
    nrmax = 4;
    wils_loop2 = (Real *)malloc(nth*nrmax*sizeof(Real));
    for(t=0;t<nth;t++) for(r=0;r<nrmax;r++){
	wils_loop2[r+nrmax*t] = 0.0;
    }


    /* "sqrt(5)" loops                             */
    /* ----------------------------------------    */

    /* Off-set for next bunch of Wilson loops */
    r_off = 0;

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
	    su3mat_copy( &(s->link[TUP]), &(s->t_link_f));
	}

	/* Make double links in dir1 direction */
	wait_gather( mtag[dir1]);
	FORALLSITES(i,s){
	    mult_su3_nn( &(s->link[dir1]), (su3_matrix *)(gen_pt[dir1][i]),
		&(s->s_link));
	}
	cleanup_gather( mtag[dir1]);

	/* Gather the double links from dir2 direction */
	mtag[dir2] = start_gather_site( F_OFFSET(s_link), sizeof(su3_matrix),
	    dir2, EVENANDODD, gen_pt[dir2] );

	/* Make first corner */
	wait_general_gather( gmtag);
	FORALLSITES(i,s){
	    mult_su3_nn( &(s->s_link), (su3_matrix *)(gen_pt[4][i]),
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
	gmtag = start_general_gather_site( F_OFFSET(t_link_f), sizeof(su3_matrix),
	    disp, EVENANDODD, gen_pt[4] );


	/* Recursively construct the space-like segments and compute
	   the Wilson loops with that segment */

	for(r=0;r<2;r++){

	    if( r==0 ){
		FORALLSITES(i,s){
		    su3mat_copy( &(s->diag), &(s->s_link));
		}
	    }
	    else{
		wait_general_gather( gmtag);
		FORALLSITES(i,s){
		    su3mat_copy( (su3_matrix *)(gen_pt[4][i]), &(s->staple));
		}
		cleanup_general_gather( gmtag);

		/* Inbetween gather time-like links across the diagonal. */
		gmtag = start_general_gather_site( F_OFFSET(t_link_f),
		    sizeof(su3_matrix), disp, EVENANDODD, gen_pt[4] );

		FORALLSITES(i,s){
		    mult_su3_nn( &(s->diag), &(s->staple), &(s->s_link));
		}
	    }

	    FORALLSITES(i,s){
		su3mat_copy( &(s->s_link), &(s->s_link_f));
	    }

	    /* Start gather of forward space-like segments */
	    mtag[TUP] = start_gather_site( F_OFFSET(s_link_f), sizeof(su3_matrix),
		TUP, EVENANDODD, gen_pt[TUP] );

	    /* Collect forward time-like links. */
	    wait_general_gather( gmtag);
	    FORALLSITES(i,s){
		su3mat_copy( (su3_matrix *)(gen_pt[4][i]), &(s->staple));
	    }
	    FORALLSITES(i,s){
		su3mat_copy( &(s->staple), &(s->t_link_f));
	    }
	    cleanup_general_gather( gmtag);

	    /* Inbetween gather space-links across the diagonal for next r. */
	    if( r<(nxh-1) ){
		gmtag = start_general_gather_site( F_OFFSET(s_link),
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
			wils_loop2[r+r_off+nrmax*t] +=
			    realtrace_su3( &tmat2, &(s->s_link));
		    }
		    else
		    {
			wils_loop2[r+r_off+nrmax*t] +=
			    realtrace_su3( &(s->s_link_f), &(s->s_link));
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
	    su3mat_copy( &(s->link[TUP]), &(s->t_link_f));
	}

	/* Make double links in dir1 direction */
	wait_gather( mtag[dir1]);
	FORALLSITES(i,s){
	    mult_su3_nn( &(s->link[dir1]), (su3_matrix *)(gen_pt[dir1][i]),
		&(s->s_link));
	}
	cleanup_gather( mtag[dir1]);

	/* Make one corner and then gather it */
	FORALLSITES(i,s){
	    mult_su3_an( &(s->link[dir2]), &(s->s_link), &(s->staple));
	}
	mtag[dir2] = start_gather_site( F_OFFSET(staple), sizeof(su3_matrix),
	    OPP_DIR(dir2), EVENANDODD, gen_pt[dir2] );

	/* Make second corner */
	wait_general_gather( gmtag);
	FORALLSITES(i,s){
	    mult_su3_na( &(s->s_link), (su3_matrix *)(gen_pt[4][i]),
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
	gmtag = start_general_gather_site( F_OFFSET(t_link_f), sizeof(su3_matrix),
	    disp, EVENANDODD, gen_pt[4] );


	/* Recursively construct the space-like segments and compute
	   the Wilson loops with that segment */

	for(r=0;r<2;r++){

	    if( r==0 ){
		FORALLSITES(i,s){
		    su3mat_copy( &(s->diag), &(s->s_link));
		}
	    }
	    else{
		wait_general_gather( gmtag);
		FORALLSITES(i,s){
		    su3mat_copy( (su3_matrix *)(gen_pt[4][i]), &(s->staple));
		}
		cleanup_general_gather( gmtag);

		/* Inbetween gather time-like links across the diagonal. */
		gmtag = start_general_gather_site( F_OFFSET(t_link_f),
		    sizeof(su3_matrix), disp, EVENANDODD, gen_pt[4] );

		FORALLSITES(i,s){
		    mult_su3_nn( &(s->diag), &(s->staple), &(s->s_link));
		}
	    }

	    FORALLSITES(i,s){
		su3mat_copy( &(s->s_link), &(s->s_link_f));
	    }

	    /* Start gather of forward space-like segments */
	    mtag[TUP] = start_gather_site( F_OFFSET(s_link_f), sizeof(su3_matrix),
		TUP, EVENANDODD, gen_pt[TUP] );

	    /* Collect forward time-like links. */
	    wait_general_gather( gmtag);
	    FORALLSITES(i,s){
		su3mat_copy( (su3_matrix *)(gen_pt[4][i]), &(s->staple));
	    }
	    FORALLSITES(i,s){
		su3mat_copy( &(s->staple), &(s->t_link_f));
	    }
	    cleanup_general_gather( gmtag);

	    /* Inbetween gather space-links across the diagonal for next r. */
	    if( r<(nxh-1) ){
		gmtag = start_general_gather_site( F_OFFSET(s_link),
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
			wils_loop2[r+r_off+nrmax*t] +=
			    realtrace_su3( &tmat2, &(s->s_link));
		    }
		    else
		    {
			wils_loop2[r+r_off+nrmax*t] +=
			    realtrace_su3( &(s->s_link_f), &(s->s_link));
		    }
		}

	    } /* end loop over t */

	} /* end loop over r */

    } /* end loop over dir1 != dir2 */



    /* "sqrt(8)" loops                             */
    /* ----------------------------------------    */
    r_off = 2;

    for(dir1=XUP;dir1<=YUP;dir1++) for(dir2=dir1+1;dir2<=ZUP;dir2++){

	/* Do off-axis "sqrt(8)" loops in (dir1,dir2)-plane */
        /* ------------------------------------------------ */

	/* First construct the "diagonal" link in (2*dir1,2*dir2) direction */
	/* First construct the "double links" */
	mtag[dir1] = start_gather_site( F_OFFSET(link[dir1]), sizeof(su3_matrix),
	    dir1, EVENANDODD, gen_pt[dir1] );
	mtag[dir2] = start_gather_site( F_OFFSET(link[dir2]), sizeof(su3_matrix),
	    dir2, EVENANDODD, gen_pt[dir2] );

	wait_gather( mtag[dir1]);
	FORALLSITES(i,s){
	    mult_su3_nn( &(s->link[dir1]), (su3_matrix *)(gen_pt[dir1][i]),
		&(s->s_link));
	}
	cleanup_gather( mtag[dir1]);

	/* Start gather of dir1-double-link from site "2*dir2" */
	for(i=XUP;i<=TUP;i++)disp[i]=0;
	disp[dir2] = 2;
	gmtag = start_general_gather_site( F_OFFSET(s_link), sizeof(su3_matrix),
	    disp, EVENANDODD, gen_pt[4] );

	wait_gather( mtag[dir2]);
	FORALLSITES(i,s){
	    mult_su3_nn( &(s->link[dir2]), (su3_matrix *)(gen_pt[dir2][i]),
		&(s->s_link_f));
	}
	cleanup_gather( mtag[dir2]);

	/* Make first corner */
	wait_general_gather( gmtag);
        FORALLSITES(i,s){
	    mult_su3_nn( &(s->s_link_f), (su3_matrix *)(gen_pt[4][i]),
		&(s->diag));
        }
	cleanup_general_gather( gmtag);

	/* Start gather of dir2-double-link from site "2*dir1" */
	for(i=XUP;i<=TUP;i++)disp[i]=0;
	disp[dir1] = 2;
	gmtag = start_general_gather_site( F_OFFSET(s_link_f), sizeof(su3_matrix),
	    disp, EVENANDODD, gen_pt[4] );

	FORALLSITES(i,s){
	    su3mat_copy( &(s->link[TUP]), &(s->t_link_f));
	}

	/* Make second corner and add to first */
	ftmp = 0.5;
	wait_general_gather( gmtag);
	FORALLSITES(i,s){
	    mult_su3_nn( &(s->s_link), (su3_matrix *)(gen_pt[4][i]),
		&tmat1);
	    add_su3_matrix( &(s->diag), &tmat1, &(s->diag));
	    scalar_mult_su3_matrix( &(s->diag), ftmp, &(s->diag));
	}
	cleanup_general_gather( gmtag);


	/* Start gather of time-like links across the diagonal. */
	for(i=XUP;i<=TUP;i++)disp[i]=0;
	disp[dir1] = 2;
	disp[dir2] = 2;
	gmtag = start_general_gather_site( F_OFFSET(t_link_f), sizeof(su3_matrix),
	    disp, EVENANDODD, gen_pt[4] );


	/* Recursively construct the space-like segments and compute
	   the Wilson loops with that segment */

	for(r=0;r<2;r++){

	    if( r==0 ){
		FORALLSITES(i,s){
		    su3mat_copy( &(s->diag), &(s->s_link));
		}
	    }
	    else{
		wait_general_gather( gmtag);
		FORALLSITES(i,s){
		    su3mat_copy( (su3_matrix *)(gen_pt[4][i]), &(s->staple));
		}
		cleanup_general_gather( gmtag);

		/* Inbetween gather time-like links across the diagonal. */
		gmtag = start_general_gather_site( F_OFFSET(t_link_f),
		    sizeof(su3_matrix), disp, EVENANDODD, gen_pt[4] );

		FORALLSITES(i,s){
		    mult_su3_nn( &(s->diag), &(s->staple), &(s->s_link));
		}
	    }

	    FORALLSITES(i,s){
		su3mat_copy( &(s->s_link), &(s->s_link_f));
	    }

	    /* Start gather of forward space-like segments */
	    mtag[TUP] = start_gather_site( F_OFFSET(s_link_f), sizeof(su3_matrix),
		TUP, EVENANDODD, gen_pt[TUP] );

	    /* Collect forward time-like links. */
	    wait_general_gather( gmtag);
	    FORALLSITES(i,s){
		su3mat_copy( (su3_matrix *)(gen_pt[4][i]), &(s->staple));
	    }
	    FORALLSITES(i,s){
		su3mat_copy( &(s->staple), &(s->t_link_f));
	    }
	    cleanup_general_gather( gmtag);

	    /* Inbetween gather space-links across the diagonal for next r. */
	    if( r<(nxh-1) ){
		gmtag = start_general_gather_site( F_OFFSET(s_link),
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
			wils_loop2[r+r_off+nrmax*t] +=
			    realtrace_su3( &tmat2, &(s->s_link));
		    }
		    else
		    {
			wils_loop2[r+r_off+nrmax*t] +=
			    realtrace_su3( &(s->s_link_f), &(s->s_link));
		    }
		}

	    } /* end loop over t */

	} /* end loop over r */


	/* Do off-axis "sqrt(8)" loops in (dir1,-dir2)-plane */
        /* ------------------------------------------------ */

	/* First construct the "diagonal" link in the (2*dir1,-2*dir2) dir */
	/* First construct the "double links" */
	mtag[dir2] = start_gather_site( F_OFFSET(link[dir2]), sizeof(su3_matrix),
	    dir2, EVENANDODD, gen_pt[dir2] );
	mtag[dir1] = start_gather_site( F_OFFSET(link[dir1]), sizeof(su3_matrix),
	    dir1, EVENANDODD, gen_pt[dir1] );

	wait_gather( mtag[dir2]);
	FORALLSITES(i,s){
	    mult_su3_nn( &(s->link[dir2]), (su3_matrix *)(gen_pt[dir2][i]),
		&(s->s_link_f));
	}
	cleanup_gather( mtag[dir2]);

	/* Start gather of dir2-double-link from site "(2*dir1,-2*dir2)" */
	for(i=XUP;i<=TUP;i++)disp[i]=0;
	disp[dir1] = 2;
	disp[dir2] = -2;
	gmtag = start_general_gather_site( F_OFFSET(s_link_f), sizeof(su3_matrix),
	    disp, EVENANDODD, gen_pt[4] );

	wait_gather( mtag[dir1]);
	FORALLSITES(i,s){
	    mult_su3_nn( &(s->link[dir1]), (su3_matrix *)(gen_pt[dir1][i]),
		&(s->s_link));
	}
	cleanup_gather( mtag[dir1]);

	/* Make first corner */
        FORALLSITES(i,s){
	    mult_su3_an( &(s->s_link_f), &(s->s_link), &(s->staple));
        }

	/* Make second corner */
	wait_general_gather( gmtag);
        FORALLSITES(i,s){
	    mult_su3_na( &(s->s_link), (su3_matrix *)(gen_pt[4][i]),
		&(s->diag));
        }
	cleanup_general_gather( gmtag);

	/* Start gather first corner from site "-2*dir2" */
	for(i=XUP;i<=TUP;i++)disp[i]=0;
	disp[dir2] = -2;
	gmtag = start_general_gather_site( F_OFFSET(staple), sizeof(su3_matrix),
	    disp, EVENANDODD, gen_pt[4] );

	FORALLSITES(i,s){
	    su3mat_copy( &(s->link[TUP]), &(s->t_link_f));
	}

	/* Collect first corner and add to second */
	ftmp = 0.5;
	wait_general_gather( gmtag);
	FORALLSITES(i,s){
	    add_su3_matrix( &(s->diag), (su3_matrix *)(gen_pt[4][i]),
		 &(s->diag));
	    scalar_mult_su3_matrix( &(s->diag), ftmp, &(s->diag));
	}
	cleanup_general_gather( gmtag);


	/* Start gather of time-like links across the diagonal. */
	for(i=XUP;i<=TUP;i++)disp[i]=0;
	disp[dir1] =  2;
	disp[dir2] = -2;
	gmtag = start_general_gather_site( F_OFFSET(t_link_f), sizeof(su3_matrix),
	    disp, EVENANDODD, gen_pt[4] );

	/* Recursively construct the space-like segments and compute
	   the Wilson loops with that segment */


	for(r=0;r<2;r++){

	    if( r==0 ){
		FORALLSITES(i,s){
		    su3mat_copy( &(s->diag), &(s->s_link));
		}
	    }
	    else{
		wait_general_gather( gmtag);
		FORALLSITES(i,s){
		    su3mat_copy( (su3_matrix *)(gen_pt[4][i]), &(s->staple));
		}
		cleanup_general_gather( gmtag);

		/* Inbetween gather time-like links across the diagonal. */
		gmtag = start_general_gather_site( F_OFFSET(t_link_f),
		    sizeof(su3_matrix), disp, EVENANDODD, gen_pt[4] );

		FORALLSITES(i,s){
		    mult_su3_nn( &(s->diag), &(s->staple), &(s->s_link));
		}
	    }

	    FORALLSITES(i,s){
		su3mat_copy( &(s->s_link), &(s->s_link_f));
	    }

	    /* Start gather of forward space-like segments */
	    mtag[TUP] = start_gather_site( F_OFFSET(s_link_f), sizeof(su3_matrix),
		TUP, EVENANDODD, gen_pt[TUP] );

	    /* Collect forward time-like links. */
	    wait_general_gather( gmtag);
	    FORALLSITES(i,s){
		su3mat_copy( (su3_matrix *)(gen_pt[4][i]), &(s->staple));
	    }
	    FORALLSITES(i,s){
		su3mat_copy( &(s->staple), &(s->t_link_f));
	    }
	    cleanup_general_gather( gmtag);

	    /* Inbetween gather space-links across the diagonal for next r. */
	    if( r<(nxh-1) ){
		gmtag = start_general_gather_site( F_OFFSET(s_link),
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
			wils_loop2[r+r_off+nrmax*t] +=
			    realtrace_su3( &tmat2, &(s->s_link));
		    }
		    else
		    {
			wils_loop2[r+r_off+nrmax*t] +=
			    realtrace_su3( &(s->s_link_f), &(s->s_link));
		    }
		}

	    } /* end loop over t */

	} /* end loop over r */

	} /* end loop over dir1 /dir2 */



    /* Normalize and print the Wilson loops */
    for(t=0;t<nth;t++){
	for(r=0;r<2;r++){
	    ftmp = wils_loop2[r+nrmax*t];
	    g_floatsum( &ftmp);
	    ftmp /= (Real)(36*volume);
	    if(this_node == 0)printf("WILS_LOOP2_%d  %d  %d  %e\n",
		 tot_smear, r, t, (double)ftmp);
	}

	r_off = 2;
	for(r=0;r<2;r++){
	    ftmp = wils_loop2[r+r_off+nrmax*t];
	    g_floatsum( &ftmp);
	    ftmp /= (Real)(18*volume);
	    if(this_node == 0)printf("WILS_LOOP2_%d  %d  %d  %e\n",
		 tot_smear, r+r_off, t, (double)ftmp);
	}

    }

    free( wils_loop2);

} /* w_loop2 */

