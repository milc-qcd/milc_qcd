/************************** smear_links.c *******************************/
/* MIMD version 6 */
/* Doug Toussaint 1/30/96 */
/* modification 7/96 to take staple weight as arguments, normalization
	fixed so sum of weights=1 DT.
   4/97 don't automatically copy result back to source.
   8/30 add provisions for ape_projection and SPACE_ONLY, DT
*/
#define SPACE_ONLY 1
#define APE_PROJECT
#define NHIT 10

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

   If SPACE_ONLY is defined to be 1, only smear in spatial directions
   If it's defined to be zero, smear isotropically
*/

#include "generic_ks_includes.h"
#include "../include/field_strength.h"
void ape_project_su3mat( su3_matrix *matp );

void smear_links( field_offset src, field_offset dest, Real staple_weight ){
register int i,dir1,dir2;
register site *s;
su3_matrix tmat1,tmat2;
msg_tag *mtag0,*mtag1;
Real simple_weight, norm_factor;

    simple_weight =  1.0/staple_weight;

    node0_printf( "Smearing: staple_weight = %e\n", staple_weight);
#ifdef APE_PROJECT
    node0_printf( "APE projection is on, NHIT = %d\n",NHIT);
#endif
    for(dir1=XUP;dir1<=TUP;dir1++){
	FORALLSITES(i,s){
	    scalar_mult_su3_matrix(
	       &(((su3_matrix *)F_PT(s,src))[dir1]), simple_weight,
	       &(((su3_matrix *)F_PT(s,dest))[dir1]) );
	}
	if( SPACE_ONLY==1 && dir1 != TUP) norm_factor = 1.0/(4.0*u0*u0+simple_weight);
	else norm_factor = 1.0/(6.0*u0*u0+simple_weight);
	for(dir2=XUP;dir2<=(SPACE_ONLY==1?ZUP:TUP);dir2++)if(dir2!=dir1){

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

        /* if desired, APE project to unitary matrix */
#ifdef APE_PROJECT
        FORALLSITES(i,s){
            ape_project_su3mat( &(((su3_matrix *)F_PT(s,dest))[dir1]) );
	}
#endif

    } /*dir1 loop */

} /* smear_links */

/* "APE project" a single su3_matrix to SU3 element, works in place */
void ape_project_su3mat( su3_matrix *matp ){
    su3_matrix tmat2;
    register int ii;

    /* copy original matrix to temp space, initial unitary guess
       from reunit_su3() */
    tmat2 = *matp;
    ii = reunit_su3( matp );
    project_su3(matp, &tmat2, 3*NHIT);

} /* ape_project_su3mat() */

