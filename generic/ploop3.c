/****************** ploop3.c ************************************/
/* MIMD version 7 */
/* evaluate the Polyakov loops.  This version uses general_gathers. */
/* It assumes that nt is even.  Actually, all the code does.  */
/* DT 12/97 use local matrix "tmat" instead of "tempmat2" */

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
su3_matrix tmat, *tempmat1;
int d[4];

    tempmat1 = (su3_matrix *)malloc(sites_on_node*sizeof(su3_matrix));
    if(tempmat1 == NULL){
      printf("ploop3: Can't malloc temporary\n");
      terminate(1);
  }

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
	    &tempmat1[i]);
    }
    cleanup_gather(tag);

    for(t=2;t<nt;t+=2){
	d[TUP] = t;	/* distance from which to gather */
	tag=start_general_gather_field( tempmat1, sizeof(su3_matrix),
	    d, EVEN, gen_pt[0] );
	wait_general_gather(tag);
        FOREVENSITES(i,st){
	    if( st->t > 1 )continue;  /* only compute on first two slices */
	    mult_su3_nn( &tempmat1[i], (su3_matrix *)gen_pt[0][i], &tmat );
	    tempmat1[i] = tmat;
	    /* We overwrite tempmat1 on the first two time slices,
		leaving the others undisturbed so we can still gather
		them. */
	}
	cleanup_general_gather(tag);
    }
    FOREVENSITES(i,st){
	if( st->t > 1 )continue;
	plp = trace_su3( &tempmat1[i] );
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
    free(tempmat1);
    return(plp);
}

