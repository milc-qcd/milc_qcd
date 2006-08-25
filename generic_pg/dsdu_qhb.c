/****** dsdu_qhb.c  -- compute the staple ******************/
/* MIMD version 6 */
/* UMH: Combined with Schroedinger functional version, Jan 2000 */

#include "generic_pg_includes.h"
#include "../include/loopend.h"

/* The following redefinition of FORSOMEPARITY eliminates some useless
   looping, so saves a tiny bit of time.  It could be extended to all
   architectures, provided we were willing to compile with EVENFIRST
   always. */

void dsdu_qhb(int dir1,int parity)
{
register int i,dir2,otherparity=0;
register site *st;
msg_tag *tag0,*tag1,*tag2,*tag3;
int start;
su3_matrix tmat1,tmat2;

    switch(parity) {
	case EVEN:		otherparity=ODD;	break;
	case ODD:		otherparity=EVEN;	break;
	case EVENANDODD:	otherparity=EVENANDODD;	break;
    }

    /* Loop over other directions, computing force from plaquettes in
       the dir1,dir2 plane */
    start=1; /* indicates staple sum not initialized */

    for(dir2=XUP;dir2<=TUP;dir2++)if(dir2 != dir1)
    {
	/* get link[dir2] from direction dir1 on other parity */
	tag0 = start_gather_site( F_OFFSET(link[dir2]), sizeof(su3_matrix),
	    dir1, otherparity, gen_pt[0] );

	/* get link[dir2] from direction dir1 */
	tag1 = start_gather_site( F_OFFSET(link[dir2]), sizeof(su3_matrix),
	    dir1, parity, gen_pt[1] );

	/* get link[dir1] from direction dir2 */
	tag2 = start_gather_site( F_OFFSET(link[dir1]), sizeof(su3_matrix),
	    dir2, parity, gen_pt[2] );

	/* Lower staple (computed at backward site) */
	/* Note: For SCHROED_FUN we don't care if we get this wrong for
	         dir1<TUP and t=0, since then the staple will not be used,
	         as those links are frozen */
	wait_gather(tag0);
	FORSOMEPARITY(i,st,otherparity){
#ifdef SCHROED_FUN
	    if(st->t==(nt-1) && dir1==TUP){
		mult_su3_an( &(st->link[dir2]), &(st->link[dir1]), &tmat1);
		mult_su3_nn( &tmat1, &(st->boundary[dir2]), &(st->tempmat1));
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
	} END_LOOP
	cleanup_gather(tag0);
	/* get tempmat1 from direction -dir2 */
	tag3 = start_gather_site( F_OFFSET(tempmat1), sizeof(su3_matrix),
	    OPP_DIR(dir2), parity, gen_pt[3] );

	/* Upper staple */
	wait_gather(tag1);
	wait_gather(tag2);
        if(start){  /* this is the first contribution to staple */
#ifdef SCHROED_FUN
	    FORSOMEPARITY(i,st,parity) if(dir1==TUP || st->t>0){
		if(st->t==(nt-1) && dir2==TUP){
		    mult_su3_nn( &(st->link[dir2]), &(st->boundary[dir1]),
			&tmat1 );
		    mult_su3_na( &tmat1, (su3_matrix *)gen_pt[1][i],
			&(st->staple) );
		}
		else if(st->t==(nt-1) && dir1==TUP){
		    mult_su3_nn( &(st->link[dir2]),
			(su3_matrix *)gen_pt[2][i], &tmat1 );
		    mult_su3_na( &tmat1, &(st->boundary[dir2]), &(st->staple) );
		}
		else{
		    mult_su3_nn( &(st->link[dir2]),
			(su3_matrix *)gen_pt[2][i], &tmat1 );
		    mult_su3_na( &tmat1, (su3_matrix *)gen_pt[1][i],
			&(st->staple) );
		}
	    } END_LOOP
#else
	    FORSOMEPARITY(i,st,parity){
		mult_su3_nn( &(st->link[dir2]),
		    (su3_matrix *)gen_pt[2][i], &tmat1 );
		mult_su3_na( &tmat1, (su3_matrix *)gen_pt[1][i],
		    &(st->staple) );
	    } END_LOOP
#endif
	    start=0; 
	}
	else{
#ifdef SCHROED_FUN
	    FORSOMEPARITY(i,st,parity) if(dir1==TUP || st->t>0){
		if(st->t==(nt-1) && dir2==TUP){
		    mult_su3_nn( &(st->link[dir2]), &(st->boundary[dir1]),
			&tmat1 );
		    mult_su3_na( &tmat1, (su3_matrix *)gen_pt[1][i], &tmat2 );
		}
		else if(st->t==(nt-1) && dir1==TUP){
		    mult_su3_nn( &(st->link[dir2]),
			(su3_matrix *)gen_pt[2][i], &tmat1 );
		    mult_su3_na( &tmat1, &(st->boundary[dir2]), &tmat2 );
		}
		else{
		    mult_su3_nn( &(st->link[dir2]),
			(su3_matrix *)gen_pt[2][i], &tmat1 );
		    mult_su3_na( &tmat1, (su3_matrix *)gen_pt[1][i], &tmat2 );
		}
		add_su3_matrix( &(st->staple), &tmat2, &(st->staple));
	    } END_LOOP
#else
	    FORSOMEPARITY(i,st,parity){
		mult_su3_nn( &(st->link[dir2]),
		    (su3_matrix *)gen_pt[2][i], &tmat1 );
		mult_su3_na( &tmat1, (su3_matrix *)gen_pt[1][i], &tmat2 );
		add_su3_matrix( &(st->staple), &tmat2, &(st->staple));
	    } END_LOOP
#endif
	} /* upper staple */
	cleanup_gather(tag1);
	cleanup_gather(tag2);

	/* Add lower staple */
	wait_gather(tag3);
#ifdef SCHROED_FUN
	FORSOMEPARITY(i,st,parity) if(dir1==TUP || st->t>0){
#else
	FORSOMEPARITY(i,st,parity){
#endif
	    add_su3_matrix( &(st->staple), (su3_matrix *)gen_pt[3][i],
		&(st->staple));
	}  END_LOOP	/* lower staple */
	cleanup_gather(tag3);
    }
}

