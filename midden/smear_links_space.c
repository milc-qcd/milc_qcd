/************************** smear_links_space.c *******************************/
/* MIMD version 6 */
/* Doug Toussaint 1/30/96 */
/* 7/17/96 version4 to smear only in space directions */

/* Construct smeared links in dest[dir1].   Smeared link is sum
   of link plus all staples:

		---------		-------> dir1
		|	|
		|	|
		|	|		^
		X--------		| dir2
		|	|		|
		|	|
		|	|
		---------
   dest is probably secretly defined to overlap conjugate gradient
   vectors in the Wilson code.  Also need a temporary matrix hyb_tempmat1.

   For the moment, we add the simple link in with relative weight one,
   arbitrary normalization factor to keep F_mu_nu reasonable size.
   */
#define SPACE_SIMPLE_WEIGHT 1.0	/* relative weight of spatial simple link */
#define SPACE_NORM_FACTOR 0.3	/* arbitrary rescaling of result */
#define TIME_SIMPLE_WEIGHT 1.0	/* relative weight of time simple link */
#define TIME_NORM_FACTOR 0.175	/* arbitrary rescaling of result */
	/* SPACE_NORM_FACTOR = 1 / (SIMPLE_WEIGHT+3) is right for ordered
	   lattice.  Something larger for disordered.  Perhaps
	   1  / ( SIMPLE_WEIGHT + 3*u_0^2 ), where u_0 is tadpole
	   improvement parameter, or perhaps determine NORM_FACTOR
	   by trial at beginning of run. */

#include "cl_hyb_includes.h"

void smear_links( field_offset src, field_offset dest ){
register int i,dir1,dir2;
register site *s;
su3_matrix tmat1,tmat2;
msg_tag *mtag0,*mtag1;
Real simple_weight, norm_factor;

    if(this_node==0){
	printf(
	"Smearing: space_simple_weight = %e , space_norm_factor = %e\n",
	space_simple_weight,space_norm_factor);
	printf(
	"           time_simple_weight = %e ,  time_norm_factor = %e\n",
	time_simple_weight,time_norm_factor);
    }
    for(dir1=XUP;dir1<=TUP;dir1++){
	if(dir1==TUP){
	    simple_weight=time_simple_weight;
	    norm_factor=time_norm_factor;
	} else {
	    simple_weight=space_simple_weight;
	    norm_factor=space_norm_factor;
	}
	FORALLSITES(i,s){
	    scalar_mult_su3_matrix(
	       &(((su3_matrix *)F_PT(s,src))[dir1]), simple_weight,
	       &(((su3_matrix *)F_PT(s,dest))[dir1]) );
	}
	for(dir2=XUP;dir2<TUP;dir2++)if(dir2!=dir1){ /* dir2 only spatial */

	    /* Upper staple, and simple link */
	    mtag0 = start_gather( src+dir2*sizeof(su3_matrix),
		sizeof(su3_matrix), dir1, EVENANDODD, gen_pt[0] );
	    mtag1 = start_gather( src+dir1*sizeof(su3_matrix),
		 sizeof(su3_matrix), dir2, EVENANDODD, gen_pt[1] );
	    wait_gather(mtag0);
	    wait_gather(mtag1);

	    FORALLSITES(i,s){
		mult_su3_na( (su3_matrix *)gen_pt[1][i],
			     (su3_matrix *)gen_pt[0][i], &tmat1 );
		mult_su3_nn( &(((su3_matrix *)F_PT(s,src))[dir2]),
		    &tmat1, &tmat2 );
		add_su3_matrix(
		    &(((su3_matrix *)F_PT(s,dest))[dir1]), &tmat2,
		    &(((su3_matrix *)F_PT(s,dest))[dir1]) );
	    }
	    cleanup_gather(mtag0);
	    cleanup_gather(mtag1);

	    /* lower staple */
	    mtag0 = start_gather( src+dir2*sizeof(su3_matrix),
		 sizeof(su3_matrix), dir1, EVENANDODD, gen_pt[0] );
	    wait_gather(mtag0);
	    FORALLSITES(i,s){
		mult_su3_nn( &(((su3_matrix *)F_PT(s,src))[dir1]),
		    (su3_matrix *)gen_pt[0][i], &tmat1 );
		mult_su3_an( &(((su3_matrix *)F_PT(s,src))[dir2]),
		    &tmat1, &(s->hyb_tempmat1) );
	    }
	    cleanup_gather(mtag0);
	    mtag1 = start_gather( F_OFFSET(hyb_tempmat1), sizeof(su3_matrix),
                OPP_DIR(dir2), EVENANDODD, gen_pt[1] );
	    wait_gather(mtag1);
	    FORALLSITES(i,s){
		add_su3_matrix( &(((su3_matrix *)F_PT(s,dest))[dir1]),
		    (su3_matrix *)gen_pt[1][i],
		    &(((su3_matrix *)F_PT(s,dest))[dir1]) );
	    }
	    cleanup_gather(mtag1);

	} /* dir2 loop */
	FORALLSITES(i,s){
	    scalar_mult_su3_matrix( &(((su3_matrix *)F_PT(s,dest))[dir1]),
		norm_factor, &(((su3_matrix *)F_PT(s,dest))[dir1]) );
	}
    } /*dir1 loop */
} /* smear_links */

