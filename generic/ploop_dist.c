/****************** ploop_dist.c ************************************/
/* MIMD version 7 */
/* evaluate the Polyakov loops.  This version uses general_gathers. */
/* It assumes that nt is even.  Actually, all the code does.  */
/* Modified C. DeTar July 20, 1993 */
/* Reports the distribution of ploop values as a function of each x y z */
/* Similar to ploop3.c, except doesn't use tmat */

/* Macros ...

   PLOOPDIST
   Turns on reporting of the ploop values

 */

#include "generic_includes.h"

complex ploop() {
register int i,t;
register site *st;
msg_tag *tag;
complex sum;
complex plp;
int d[4];
#ifdef PLOOPDIST
int x,y,z;
#endif

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
	tag=start_general_gather_site( F_OFFSET(tempmat1), sizeof(su3_matrix),
	    d, EVEN, gen_pt[0] );
	wait_general_gather(tag);
        FOREVENSITES(i,st){
	    if( st->t > 1 )continue;  /* only compute on first two slices */
	    mult_su3_nn( &(st->tempmat1), (su3_matrix *)gen_pt[0][i], 
			 &(st->tempmat2));
	    lattice[i].tempmat1 = lattice[i].tempmat2;
	    /* We overwrite tempmat1 on the first two time slices,
		leaving the others undisturbed so we can still gather
		them. */
	}
	cleanup_general_gather(tag);
    }
    FOREVENSITES(i,st){
	if( st->t > 1 )continue;
	plp = trace_su3( &(st->tempmat1) );
#ifdef PLOOPDIST
	/* Save result in tempmat1 for t = 0 or 1 slice */
	st->tempmat1.e[0][0] = plp;
#endif
	CSUM(sum,plp);
    }

#ifdef PLOOPDIST
    /* Report ploop distribution */
    for(x=0;x<nx;x++)for(y=0;y<ny;y++)for(z=0;z<nz;z++)
      {
	t = (x+y+z)%2;
	i=node_index(x,y,z,t);
	if( node_number(x,y,z,t) != mynode() )plp = cmplx(0.,0.);
	else                                plp = lattice[i].tempmat1.e[0][0];
	g_sync();
	g_complexsum(&plp);
	if(mynode()==0)printf("PLOOP %d %d %d %.8e %.8e\n",
			      x,y,z,plp.real,plp.imag);
      }
#endif

    g_complexsum( &sum );
    plp.real = sum.real /((Real)(nx*ny*nz));
    plp.imag = sum.imag /((Real)(nx*ny*nz));
    return(plp);
}
