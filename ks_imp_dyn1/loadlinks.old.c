/******* loadlinks.old.c - load fancy links for improved fermion action ****/
/* MIMD version 6 */
/* NOT MAINTAINED.  TEST BEFORE USE! */
/* Kogut-Susskind fermions  -- improved */
/* three straight links (Naik), fat links (link plus staple) */

/* Jim Hetrick, Kari Rummukainen, Doug Toussaint */

#include "ks_imp_includes.h"	/* definitions files and prototypes */

/* long link calculating routine */
void load_longlinks_old() {
  register int i;
  register site *s;
  register int dir,otherparity;
  msg_tag *tag[4];

  /**node0_printf("Loading long links\n");**/
  /** gather from all forward directions **/
  for (dir=XUP; dir<=TUP; dir++){
    tag[dir] = start_gather(F_OFFSET(link[dir]), sizeof(su3_matrix), dir,
	EVENANDODD, gen_pt[dir] );
  }

  /* wait gather, multiply links and use longlink as a temp array */
  for (dir=XUP; dir<=TUP; dir++) {
    wait_gather(tag[dir]);
    FORALLSITES(i,s) {
      mult_su3_nn( &(s->link[dir]), (su3_matrix *)gen_pt[dir][i],
	&(s->longlink[dir]) );
    }

    cleanup_gather(tag[dir]);
    /* needs to be done -- restart does not work here? */

    /* start the new gather immediately */
    tag[dir] = start_gather(F_OFFSET(longlink[dir]), sizeof(su3_matrix), dir,
	EVENANDODD, gen_pt[dir] );
  }
  
  /* wait for gathers, multiply -- result can't be put directly to longlink,
   * because longlink was used as a scratch */
  for (dir=XUP; dir<=TUP; dir++){
    
    wait_gather(tag[dir]);

    FORALLSITES(i,s) {
      mult_su3_nn( &(s->link[dir]), (su3_matrix *)gen_pt[dir][i],
	&(s->tempmat1) );      
    } 
    /* do not combine these two loops */
    FORALLSITES(i,s) {
      scalar_mult_su3_matrix( &(s->tempmat1), c3, &(s->longlink[dir]) );
    }

    cleanup_gather(tag[dir]);
 
  }

  /** set longlink flag, and return **/
  valid_longlinks = 1;
}
    


/* Fat link calculating routine */
/* Doug Toussaint 1/30/96 

   Construct smeared links in fatlink[dir1].   Smeared link is sum
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
   Uses a temporary matrix tempmat1.

   Note staples are subtracted, since KS phases are in, so each
   plaquette has one minus sign.

*/

void load_fatlinks_old() {
    register int i,dir1,dir2;
    register site *s;
    su3_matrix tmat1,tmat2;
    msg_tag *mtag0,*mtag1;
    Real simple_weight, norm_factor;

    /**node0_printf("Loading fat links, c1=%e, w3 = %e\n",c1,w3);**/
    for(dir1=XUP;dir1<=TUP;dir1++){
	FORALLSITES(i,s){
	    scalar_mult_su3_matrix(
	       &(s->link[dir1]), c1, &(s->fatlink[dir1]) );
	}
	for(dir2=XUP;dir2<=TUP;dir2++)if(dir2!=dir1){

	    /* Upper staple, and simple link */
	    mtag0 = start_gather( F_OFFSET(link[dir2]),
		sizeof(su3_matrix), dir1, EVENANDODD, gen_pt[0] );
	    mtag1 = start_gather( F_OFFSET(link[dir1]),
		 sizeof(su3_matrix), dir2, EVENANDODD, gen_pt[1] );
	    wait_gather(mtag0);
	    wait_gather(mtag1);

	    FORALLSITES(i,s){
		mult_su3_na( (su3_matrix *)gen_pt[1][i],
		    (su3_matrix *)gen_pt[0][i], &tmat1 );
		mult_su3_nn( &(s->link[dir2]), &tmat1, &tmat2 );
		scalar_mult_add_su3_matrix( &(s->fatlink[dir1]), &tmat2, -w3,
		    &(s->fatlink[dir1]) );
	    }
	    cleanup_gather(mtag0);
	    cleanup_gather(mtag1);

	    /* lower staple */
	    mtag0 = start_gather( F_OFFSET(link[dir2]),
		 sizeof(su3_matrix), dir1, EVENANDODD, gen_pt[0] );
	    wait_gather(mtag0);
	    FORALLSITES(i,s){
		mult_su3_nn( &(s->link[dir1]),
		    (su3_matrix *)gen_pt[0][i], &tmat1 );
		mult_su3_an( &(s->link[dir2]), &tmat1, &(s->tempmat1) );
	    }
	    cleanup_gather(mtag0);
	    mtag1 = start_gather( F_OFFSET(tempmat1), sizeof(su3_matrix),
                OPP_DIR(dir2), EVENANDODD, gen_pt[1] );
	    wait_gather(mtag1);
	    FORALLSITES(i,s){
		scalar_mult_add_su3_matrix( &(s->fatlink[dir1]),
		    (su3_matrix *)gen_pt[1][i], -w3, &(s->fatlink[dir1]) );
	    }
	    cleanup_gather(mtag1);

	} /* dir2 loop */
    } /*dir1 loop */

  /** set fatlink flag, and return **/
  valid_fatlinks = 1;
} /* load_fatlinks */

