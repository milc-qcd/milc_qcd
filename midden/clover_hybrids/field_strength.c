/************************** field_strength.c *******************************/
/* MIMD version 6 */
/* started 8/4/95 DT */
/* removed trace 9/3/95 CD */

/* Compute the field strength tensor.  Can add smearing if desired 

This cartoon shows how the plaquettes are calculated. The "O" at the
corner indicates where the path begins and ends.  F_mu_nu is the sum
of these plaquettes minus their adjoints.


  ^		--------<--------       --------<-------
  |dir1		|		|	|		|
  |		|		|	|		|
  |		|		|	|		|
  |		|		^	|		^
  -------->	|		|	|		|
    dir0	|		|	|		|
		|		|	|		|
		|		|	|		|
		-------->-------O       O------->--------


		--------<-------O       O-------<-------
		|		|	|		|
		|		|	|		|
		|		|	|		|
		|		^	|		^
		|		|	|		|
		|		|	|		|
		|		|	|		|
		|		|	|		|
		-------->--------       -------->--------



convention: try to use gen_pt[0] and mtag0 for links in direction dir0
(which are gathered from direction +-dir1), etc.
*/

#include "cl_hyb_includes.h"

/* We might want to make F_mu_nu out of smeared link fields */
#ifdef SMEAR
#define link smearlink
#endif

void make_field_strength() {
register int i,component,dir0,dir1;
register site *s;
int j;
su3_matrix tmat1,tmat2;
complex cc;
msg_tag *mtag0,*mtag1;
    for(component=FS_XY;component<=FS_ZT;component++){
	switch(component){
	    case FS_XY: dir0=XUP; dir1=YUP; break;
	    case FS_XZ: dir0=XUP; dir1=ZUP; break;
	    case FS_YZ: dir0=YUP; dir1=ZUP; break;
	    case FS_XT: dir0=XUP; dir1=TUP; break;
	    case FS_YT: dir0=YUP; dir1=TUP; break;
	    case FS_ZT: dir0=ZUP; dir1=TUP; break;
	}

      /* Plaquette in +dir0 +dir1 direction */
	mtag0 = start_gather( F_OFFSET(link[dir0]), sizeof(su3_matrix),
	    dir1, EVENANDODD, gen_pt[0] );
	mtag1 = start_gather( F_OFFSET(link[dir1]), sizeof(su3_matrix),
	    dir0, EVENANDODD, gen_pt[1] );

	wait_gather(mtag0);
	wait_gather(mtag1);
	FORALLSITES(i,s){
	    mult_su3_nn( &(s->link[dir0]), (su3_matrix *)(gen_pt[1][i]),
		&tmat1 );
	    mult_su3_na( &tmat1, (su3_matrix *)(gen_pt[0][i]), &tmat2 );
	    mult_su3_na( &tmat2, &(s->link[dir1]), &tmat1 );
	    su3_adjoint( &tmat1, &tmat2 );
	    sub_su3_matrix(  &tmat1, &tmat2, &(s->field_strength[component]) );
	}

	/**cleanup_gather(mtag0);   Use same gather in next plaquette**/
	cleanup_gather(mtag1);

      /* Plaquette in -dir0 +dir1 direction */
	/**mtag0 = start_gather( F_OFFSET(link[dir0]), sizeof(su3_matrix),
	    dir1, EVENANDODD, gen_pt[0] );
	wait_gather(mtag0);  Already gathered above**/

	FORALLSITES(i,s){
	    mult_su3_an( &(s->link[dir1]), &(s->link[dir0]), &tmat1 );
	    mult_su3_an( (su3_matrix *)(gen_pt[0][i]), &tmat1, &(s->hyb_tempmat1) );
	}
	mtag1 = start_gather( F_OFFSET(hyb_tempmat1), sizeof(su3_matrix),
	    OPP_DIR(dir0), EVENANDODD, gen_pt[1] );
	wait_gather(mtag1);
	FORALLSITES(i,s){
	    mult_su3_nn( &(s->link[dir1]), (su3_matrix *)(gen_pt[1][i]),
		&tmat1 );
	    su3_adjoint( &tmat1, &tmat2 );
	    add_su3_matrix( &(s->field_strength[component]), &tmat1,
		&(s->field_strength[component]) );
	    sub_su3_matrix( &(s->field_strength[component]), &tmat2,
		&(s->field_strength[component]) );
	}

	cleanup_gather(mtag0);
	cleanup_gather(mtag1);

      /* Plaquette in -dir0 -dir1 direction */
	mtag0 = start_gather( F_OFFSET(link[dir0]), sizeof(su3_matrix),
	    OPP_DIR(dir0), EVENANDODD, gen_pt[0] );
	mtag1 = start_gather( F_OFFSET(link[dir1]), sizeof(su3_matrix),
	    OPP_DIR(dir1), EVENANDODD, gen_pt[1] );
	wait_gather(mtag0);
	wait_gather(mtag1);
	FORALLSITES(i,s){
	    mult_su3_nn( (su3_matrix *)(gen_pt[0][i]), &(s->link[dir1]),
		&(s->hyb_tempmat1) );
	    mult_su3_nn( (su3_matrix *)(gen_pt[1][i]), &(s->link[dir0]),
		&(s->hyb_tempmat2) );
	}
	cleanup_gather(mtag0);
	cleanup_gather(mtag1);
	mtag0 = start_gather( F_OFFSET(hyb_tempmat1), sizeof(su3_matrix),
           OPP_DIR(dir1), EVENANDODD, gen_pt[0] );
        wait_gather(mtag0);
        mtag1 = start_gather( F_OFFSET(hyb_tempmat2), sizeof(su3_matrix),
            OPP_DIR(dir0), EVENANDODD, gen_pt[1] );
        wait_gather(mtag1);
	FORALLSITES(i,s){
	    mult_su3_an( (su3_matrix *)(gen_pt[1][i]), 
		(su3_matrix *)(gen_pt[0][i]), &tmat1 );
	    su3_adjoint( &tmat1, &tmat2 );
	    add_su3_matrix( &(s->field_strength[component]), &tmat1,
		&(s->field_strength[component]) );
	    sub_su3_matrix( &(s->field_strength[component]), &tmat2,
		&(s->field_strength[component]) );
	}
	cleanup_gather(mtag0);
	cleanup_gather(mtag1);

      /* Plaquette in +dir0 -dir1 direction */
	mtag1 = start_gather( F_OFFSET(link[dir1]), sizeof(su3_matrix),
	    dir0, EVENANDODD, gen_pt[1] );
	wait_gather(mtag1);
	FORALLSITES(i,s){
	    mult_su3_an( &(s->link[dir1]), &(s->link[dir0]), &tmat1);
	    mult_su3_nn( &tmat1, (su3_matrix *)(gen_pt[1][i]), 
		&(s->hyb_tempmat1) );
	}
	cleanup_gather(mtag1);
	mtag0 = start_gather( F_OFFSET(hyb_tempmat1), sizeof(su3_matrix),
           OPP_DIR(dir1), EVENANDODD, gen_pt[0] );
        wait_gather(mtag0);
	FORALLSITES(i,s){
	    mult_su3_na( (su3_matrix *)(gen_pt[0][i]), &(s->link[dir0]),
		&tmat1 );
	    su3_adjoint( &tmat1, &tmat2 );
	    add_su3_matrix( &(s->field_strength[component]), &tmat1,
		&(s->field_strength[component]) );
	    sub_su3_matrix( &(s->field_strength[component]), &tmat2,
		&(s->field_strength[component]) );
	}
	cleanup_gather(mtag0);
	/* Make traceless */
	FORALLSITES(i,s){
	  cc = trace_su3(&(s->field_strength[component]));
	  CMULREAL(cc,0.33333333333333333,cc);
	  for(j=0;j<3;j++)
	    CSUB(s->field_strength[component].e[j][j],cc,
		 s->field_strength[component].e[j][j]);
	}
      } /* component */
} /* make_field_strength() */


