/************************** smear_links.c *******************************/
/* MIMD version 6 */
/* Doug Toussaint 1/30/96 */
/* modification 7/96 to take staple weight as arguments, normalization
	fixed so sum of weights=1 DT.
   4/97 don't automatically copy result back to source.
*/

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
   vectors.  Also need a temporary matrix tempmat1.

*/

#include "ks_hyb_includes.h"

void smear_links( field_offset src, field_offset dest, Real staple_weight ){
register int i,dir1,dir2;
register site *s;
su3_matrix tmat1,tmat2;
msg_tag *mtag0,*mtag1;
Real simple_weight, norm_factor;

    simple_weight =  1.0/staple_weight;
    norm_factor = 1.0/(6.0+simple_weight);

    if(this_node==0)printf(
	"Smearing: simple_weight = %e , norm_factor = %e\n",
	simple_weight,norm_factor);
    for(dir1=XUP;dir1<=TUP;dir1++){
	FORALLSITES(i,s){
	    scalar_mult_su3_matrix(
	       &(((su3_matrix *)F_PT(s,src))[dir1]), simple_weight,
	       &(((su3_matrix *)F_PT(s,dest))[dir1]) );
	}
	for(dir2=XUP;dir2<=TUP;dir2++)if(dir2!=dir1){

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
		    &tmat1, &(s->tempmat1) );
	    }
	    cleanup_gather(mtag0);
	    mtag1 = start_gather( F_OFFSET(tempmat1), sizeof(su3_matrix),
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

