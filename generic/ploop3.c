/****************** ploop3.c ************************************/
/* MIMD version 6 */
/* evaluate the Polyakov loops.  This version uses general_gathers. */
/* It assumes that nt is even.  Actually, all the code does.  */
/* DT 12/97 use local matrix "tmat" instead of "tempmat2" */
/* Uses only site structure member tempmat1 */

/* Macros ...

   BPCORR
   saves the Polyakov loop value for each site
   in the site member "ploop" for later use 

 */

#include "generic_includes.h"

complex ploop() {
register int i,t;
register site *st;
msg_tag *tag;
complex sum;
complex plp;
su3_matrix tmat;
int d[4];

    sum = cmplx(0.0,0.0);
    d[XUP] = d[YUP] = d[ZUP] = 0;
    /* First multiply the link on every even site by the link above it */
    /* We will compute the Polyakov loop "at" the even sites in the 
	first two time slices. */
    tag=start_gather_site( F_OFFSET(link[TUP]), sizeof(su3_matrix),
	TUP, EVEN, gen_pt[0] );
    wait_gather(tag);
    FOREVENSITES(i,st){
	mult_su3_nn( &(st->link[TUP]), (su3_matrix *)gen_pt[0][i],
	    &(st->tempmat1));
    }
    cleanup_gather(tag);

    for(t=2;t<nt;t+=2){
	d[TUP] = t;	/* distance from which to gather */
	tag=start_general_gather( F_OFFSET(tempmat1), sizeof(su3_matrix),
	    d, EVEN, gen_pt[0] );
	wait_general_gather(tag);
        FOREVENSITES(i,st){
	    if( st->t > 1 )continue;  /* only compute on first two slices */
	    mult_su3_nn( &(st->tempmat1), (su3_matrix *)gen_pt[0][i], &tmat );
	    lattice[i].tempmat1 = tmat;
	    /* We overwrite tempmat1 on the first two time slices,
		leaving the others undisturbed so we can still gather
		them. */
	}
	cleanup_general_gather(tag);
    }
    FOREVENSITES(i,st){
	if( st->t > 1 )continue;
	plp = trace_su3( &(st->tempmat1) );
	/* Save for later correlation measurements */
	CSUM(sum,plp);
#ifdef BPCORR
	/* Save for subsequent correlation measurements */
	/* Note the results are saved on even sites in
	   slices 0 and 1 */
	st->ploop = plp;
#endif
    }
    g_complexsum( &sum );
    plp.real = sum.real /((Real)(nx*ny*nz));
    plp.imag = sum.imag /((Real)(nx*ny*nz));
    return(plp);
}

