/****** update_h.c  -- update the momentum matrices ******************/
/* MIMD version 6 */

/* UMH: Combined with Schroedinger functional version, Jan 2000 */

#include "generic_pg_includes.h"

/* update the momenta with the gauge force */
void update_h(Real eps) {
register int i,dir1,dir2;
register site *st;
msg_tag *tag0,*tag1,*tag2;
int start;
su3_matrix tmat1,tmat2;
register Real eb3;


    eb3 = eps*beta/3.0;
    /* Loop over directions, update mom[dir1] */
    for(dir1=XUP; dir1<=TUP; dir1++){
	/* Loop over other directions, computing force from plaquettes in
	   the dir1,dir2 plane */
	start=1; /* indicates staple sum not initialized */
	for(dir2=XUP;dir2<=TUP;dir2++)if(dir2 != dir1){

	    /* get link[dir2] from direction dir1 */
	    tag0 = start_gather_site( F_OFFSET(link[dir2]), sizeof(su3_matrix),
		dir1, EVENANDODD, gen_pt[0] );

	    /* Start gather for the "upper staple" */
	    tag2 = start_gather_site( F_OFFSET(link[dir1]), sizeof(su3_matrix),
		dir2, EVENANDODD, gen_pt[2] );

	    /* begin the computation "at the dir2DOWN point", we will
		later gather the intermediate result "to the home point" */

	    /* Note: For SCHROED_FUN we don't care if we get this wrong for
		     dir1<TUP and t=0, since then the staple will not be used,
		     as those links and momenta are frozen */

	    wait_gather(tag0);
	    FORALLSITES(i,st){
#ifdef SCHROED_FUN
		if(st->t==(nt-1) && dir1==TUP){
		    mult_su3_an( &(st->link[dir2]), &(st->link[dir1]), &tmat1);
		    mult_su3_nn( &tmat1, &(st->boundary[dir2]),
			&(st->tempmat1));
		}
		else{
		    mult_su3_an( &(st->link[dir2]), &(st->link[dir1]), &tmat1);
		    mult_su3_nn( &tmat1, (su3_matrix *)gen_pt[0][i],
			&(st->tempmat1));
		}
#else
		mult_su3_an( &(st->link[dir2]), &(st->link[dir1]), &tmat1);
		mult_su3_nn( &tmat1, (su3_matrix *)gen_pt[0][i],
		    &(st->tempmat1));
#endif
	    }

	    /* Gather this partial result "up to home site" */
	    tag1 = start_gather_site( F_OFFSET(tempmat1), sizeof(su3_matrix),
		OPP_DIR(dir2), EVENANDODD, gen_pt[1] );

	    /* begin the computation of the "upper" staple.  Note that
		one of the links has already been gathered, since it
		was used in computing the "lower" staple of the site
		above us (in dir2) */
	    wait_gather(tag2);
	    if(start){	/* this is the first contribution to staple */
#ifdef SCHROED_FUN
		FORALLSITES(i,st) if(dir1==TUP || st->t>0){
		    if(st->t==(nt-1) && dir2==TUP){
			mult_su3_nn( &(st->link[dir2]), &(st->boundary[dir1]),
			    &tmat1);
			mult_su3_na( &tmat1, (su3_matrix *)gen_pt[0][i],
			    &(st->staple));
		    }
		    else if(st->t==(nt-1) && dir1==TUP){
			mult_su3_nn( &(st->link[dir2]),
			    (su3_matrix *)gen_pt[2][i], &tmat1);
			mult_su3_na( &tmat1, &(st->boundary[dir2]),
			    &(st->staple));
		    }
		    else{
			mult_su3_nn( &(st->link[dir2]),
			    (su3_matrix *)gen_pt[2][i], &tmat1);
			mult_su3_na( &tmat1, (su3_matrix *)gen_pt[0][i],
			    &(st->staple));
		    }
#else
		FORALLSITES(i,st){
		    mult_su3_nn( &(st->link[dir2]),
			(su3_matrix *)gen_pt[2][i], &tmat1);
		    mult_su3_na( &tmat1, (su3_matrix *)gen_pt[0][i],
			&(st->staple));
#endif
		}
		start=0;
	    }
	    else{
#ifdef SCHROED_FUN
		FORALLSITES(i,st) if(dir1==TUP || st->t>0){
		    if(st->t==(nt-1) && dir2==TUP){
			mult_su3_nn( &(st->link[dir2]), &(st->boundary[dir1]),
			    &tmat1);
			mult_su3_na( &tmat1, (su3_matrix *)gen_pt[0][i],
			    &tmat2);
		    }
		    else if(st->t==(nt-1) && dir1==TUP){
			mult_su3_nn( &(st->link[dir2]),
			    (su3_matrix *)gen_pt[2][i], &tmat1);
			mult_su3_na( &tmat1, &(st->boundary[dir2]), &tmat2);
		    }
		    else{
			mult_su3_nn( &(st->link[dir2]),
			    (su3_matrix *)gen_pt[2][i], &tmat1);
			mult_su3_na( &tmat1, (su3_matrix *)gen_pt[0][i],
			    &tmat2);
		    }
#else
		FORALLSITES(i,st){
		    mult_su3_nn( &(st->link[dir2]),
			(su3_matrix *)gen_pt[2][i], &tmat1);
		    mult_su3_na( &tmat1, (su3_matrix *)gen_pt[0][i],  &tmat2);
#endif
		    add_su3_matrix( &(st->staple),&tmat2,&(st->staple));
		}
	    }

	    wait_gather(tag1);
#ifdef SCHROED_FUN
	    FORALLSITES(i,st) if(dir1==TUP || st->t>0){
#else
	    FORALLSITES(i,st){
#endif
		add_su3_matrix( &(st->staple), (su3_matrix *)gen_pt[1][i],
		    &(st->staple));
	    }
	    cleanup_gather(tag0);
	    cleanup_gather(tag1);
	    cleanup_gather(tag2);
	}
	/* Now multiply the staple sum by the link, then update momentum */
#ifdef SCHROED_FUN
	FORALLSITES(i,st) if(dir1==TUP || st->t>0){
#else
	FORALLSITES(i,st){
#endif
	    mult_su3_na( &(st->link[dir1]), &(st->staple), &tmat1 );
	    uncompress_anti_hermitian( &(st->mom[dir1]), &tmat2 );
	    scalar_mult_sub_su3_matrix( &tmat2, &tmat1,
		eb3, &(st->staple) );
	    make_anti_hermitian( &(st->staple), &(st->mom[dir1]) );
	}
    }
}
