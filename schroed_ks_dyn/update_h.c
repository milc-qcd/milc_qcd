/****** update_h.c  -- update the momentum matrices ******************/
/* MIMD version 7 */
/* THIS CODE NEEDS UPGRADING, NOW */

#include "schroed_ks_includes.h"

void update_h( Real eps ) {
    /* gauge field force */
    gauge_force(eps);
    /* fermionic force */
    /* First compute M*xxx in temporary vector ttt */
	/* The diagonal term in M doesn't matter */
    dslash_site( F_OFFSET(xxx), F_OFFSET(ttt), ODD );
    fermion_force(eps);
} /* update_h */

/* update the momenta with the gauge force */
void gauge_force( Real eps ) {
register int i,dir1,dir2;
register site *st;
msg_tag *tag0,*tag1,*tag2;
#ifdef REWEIGH
 msg_tag *tag3,*tag4;
#endif
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

#ifdef REWEIGH
	    /* get derivative of boundary[dirx] */
	    if(dir1 == TUP)
		tag3 = start_gather_site( F_OFFSET(boundary[dir2]),
		    sizeof(su3_matrix), OPP_DIR(TUP), EVENANDODD, gen_pt[3] );
	    else if(dir2 == TUP)
		tag3 = start_gather_site( F_OFFSET(boundary[dir1]),
		    sizeof(su3_matrix), OPP_DIR(TUP), EVENANDODD, gen_pt[3] );
#endif

	    /* Start gather for the "upper staple" */
	    tag2 = start_gather_site( F_OFFSET(link[dir1]), sizeof(su3_matrix),
		dir2, EVENANDODD, gen_pt[2] );

	    /* begin the computation "at the dir2DOWN point", we will
		later gather the intermediate result "to the home point" */

	    /* Note: we don't care if we get this wrong for dir1<TUP and t=0,
		     since then the staple will not be used, as those links
		     and momenta are frozen */

	    wait_gather(tag0);
#ifdef REWEIGH
	    if(dir1 == TUP) wait_gather(tag3);
#endif
	    FORALLSITES(i,st){
		if(st->t==(nt-1) && dir1==TUP){
		    mult_su3_an( &(st->link[dir2]), &(st->link[dir1]), &tmat1);
		    mult_su3_nn( &tmat1, &(st->boundary[dir2]),
			&(st->tempmat1));
#ifdef REWEIGH
		    mult_su3_nn( &tmat1, (su3_matrix *)gen_pt[3][i],
			&(st->tempmat2));
#endif
		}
		else{
		    mult_su3_an( &(st->link[dir2]), &(st->link[dir1]), &tmat1);
		    mult_su3_nn( &tmat1, (su3_matrix *)gen_pt[0][i],
			&(st->tempmat1));
		}
#ifdef REWEIGH
		if(st->t==0 && dir1==TUP){
		    mult_su3_an( &(st->boundary[dir2]),
			&(st->link[dir1]), &tmat1);
		    mult_su3_nn( &tmat1, (su3_matrix *)gen_pt[0][i],
			&(st->tempmat2));
		}
		if(st->t==0 && dir2==TUP){
		    mult_su3_an( &(st->link[dir2]),
			&(st->boundary[dir1]), &tmat1);
		    mult_su3_nn( &tmat1, (su3_matrix *)gen_pt[0][i],
			&(st->tempmat2));
		}
#endif
	    }

	    /* Gather this partial result "up to home site" */
	    tag1 = start_gather_site( F_OFFSET(tempmat1), sizeof(su3_matrix),
		OPP_DIR(dir2), EVENANDODD, gen_pt[1] );
#ifdef REWEIGH
	    if(dir1==TUP || dir2==TUP)
		tag4 = start_gather_site( F_OFFSET(tempmat2), sizeof(su3_matrix),
		    OPP_DIR(dir2), EVENANDODD, gen_pt[4] );
#endif

	    /* begin the computation of the "upper" staple.  Note that
		one of the links has already been gathered, since it
		was used in computing the "lower" staple of the site
		above us (in dir2) */
	    wait_gather(tag2);
	    if(start){	/* this is the first contribution to staple */
	        FORALLSITES(i,st) if(dir1==TUP || st->t>0){
		    if(st->t==(nt-1) && dir2==TUP){
			mult_su3_nn( &(st->link[dir2]), &(st->boundary[dir1]),
			    &tmat1);
			mult_su3_na( &tmat1, (su3_matrix *)gen_pt[0][i],
			    &(st->staple) );
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
			    &(st->staple) );
		    }
		}
		start=0;
	    }
	    else{
	        FORALLSITES(i,st) if(dir1==TUP || st->t>0){
		    if(st->t==(nt-1) && dir2==TUP){
			mult_su3_nn( &(st->link[dir2]), &(st->boundary[dir1]),
			    &tmat1);
			mult_su3_na( &tmat1, (su3_matrix *)gen_pt[0][i],
			    &tmat2 );
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
			    &tmat2 );
		    }
		    add_su3_matrix( &(st->staple),&tmat2,&(st->staple));
	        }
	    }

	    wait_gather(tag1);
	    FORALLSITES(i,st) if(dir1==TUP || st->t>0){
		add_su3_matrix( &(st->staple), (su3_matrix *)gen_pt[1][i],
		    &(st->staple));
	    }

#ifdef REWEIGH
	    /* Now finish the boundary derivative part */
	    if(dir1 == TUP){
		wait_gather(tag4);
		FORALLSITES(i,st){
		    if(st->t==0){
			mult_su3_nn( &(st->boundary[dir2]),
			    (su3_matrix *)gen_pt[2][i], &tmat1);
			mult_su3_na( &tmat1, (su3_matrix *)gen_pt[0][i],
			    &tmat2);
			add_su3_matrix( &tmat2, (su3_matrix *)gen_pt[4][i],
			    &tmat2);
			scalar_mult_add_su3_matrix( &(st->staple), &tmat2,
			    gamma_rv, &(st->staple));
		    }
		    else if(st->t==(nt-1)){
			mult_su3_nn( &(st->link[dir2]),
			    (su3_matrix *)gen_pt[2][i], &tmat1);
			mult_su3_na( &tmat1, (su3_matrix *)gen_pt[3][i],
			    &tmat2);
			add_su3_matrix( &tmat2, (su3_matrix *)gen_pt[4][i],
			    &tmat2);
			scalar_mult_add_su3_matrix( &(st->staple), &tmat2,
			    gamma_rv, &(st->staple));
		    }
	        }
		cleanup_gather(tag3);
		cleanup_gather(tag4);
	    }
	    else if(dir2 == TUP){
		wait_gather(tag3);
		wait_gather(tag4);
		FORALLSITES(i,st){
		    if(st->t==1){
			scalar_mult_add_su3_matrix( &(st->staple),
			    (su3_matrix *)gen_pt[4][i], gamma_rv, &(st->staple));
		    }
		    else if(st->t==(nt-1)){
			mult_su3_nn( &(st->link[dir2]),
			    (su3_matrix *)gen_pt[3][i], &tmat1);
			mult_su3_na( &tmat1, (su3_matrix *)gen_pt[0][i],
			    &tmat2 );
			scalar_mult_add_su3_matrix( &(st->staple), &tmat2,
			    gamma_rv, &(st->staple));
		    }
	        }
		cleanup_gather(tag3);
		cleanup_gather(tag4);
	    }
#endif

	    cleanup_gather(tag0);
	    cleanup_gather(tag1);
	    cleanup_gather(tag2);
	}

	/* Now multiply the staple sum by the link, then update momentum */
	FORALLSITES(i,st) if(dir1==TUP || st->t>0){
	    mult_su3_na( &(st->link[dir1]), &(st->staple), &tmat1 );
	    uncompress_anti_hermitian( &(st->mom[dir1]), &tmat2 );
	    scalar_mult_add_su3_matrix( &tmat2, &tmat1,
		eb3, &(st->staple) );
	    make_anti_hermitian( &(st->staple), &(st->mom[dir1]) );
	}
    }
}

/* update the  momenta with the fermion force */
/* Assumes that the conjugate gradient has been run, with the answer in
   xxx, and dslash_site(xxx,ttt) has been run. */
void fermion_force( Real eps ) {
register int i,dir;
register site *st;
msg_tag *tag0,*tag1;
su3_vector tvec;
su3_matrix temp1,temp2,temp3;
Real ferm_epsilon;

    ferm_epsilon = (nflavors/2.0)*eps;
    /* For even sites, gather ttt  get first one before entering loop */
    tag0 = start_gather_site( F_OFFSET(ttt), sizeof(su3_vector), XUP, EVEN,
	gen_pt[0] );

    for(dir=XUP;dir<=TUP;dir++){

	/* For odd sites, gather xxx */
	tag1 = start_gather_site( F_OFFSET(xxx), sizeof(su3_vector), dir, ODD,
	    gen_pt[1] );

	/* Note: time-like momenta at t=0 and t=nt-1 do not get
		 fermionic contributions. */

	wait_gather(tag0);
	FOREVENSITES(i,st) if(st->t > 0 && (dir<TUP || st->t<(nt-1))){
	    mult_su3_mat_vec( &(st->link[dir]), (su3_vector *)gen_pt[0][i],
		&tvec);
	    su3_projector( &tvec, &(st->xxx), &temp1);
	    uncompress_anti_hermitian( &(st->mom[dir]), &temp2);
	    scalar_mult_add_su3_matrix( &temp2, &temp1, ferm_epsilon, &temp3);
	    make_anti_hermitian( &temp3, &(st->mom[dir]));
	}
	cleanup_gather(tag0);

	/* For even sites, gather ttt */
	if(dir<TUP){
	    tag0 = start_gather_site( F_OFFSET(ttt), sizeof(su3_vector),
		dir+1, EVEN, gen_pt[0] );
	}

	wait_gather(tag1);
	FORODDSITES(i,st) if(st->t > 0 && (dir<TUP || st->t<(nt-1))){
	    mult_su3_mat_vec( &(st->link[dir]), (su3_vector *)gen_pt[1][i],
		&tvec);
	    su3_projector( &(st->ttt), &tvec, &temp1);
	    uncompress_anti_hermitian( &(st->mom[dir]), &temp2);
	    scalar_mult_add_su3_matrix( &temp2, &temp1, ferm_epsilon, &temp3);
	    make_anti_hermitian( &temp3, &(st->mom[dir]));
	}
	cleanup_gather(tag1);
    }
}

