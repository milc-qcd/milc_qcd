/************************** field_strength.c *******************************/
/* MIMD version 7 */
/* started 8/4/95 DT */
/* removed trace 9/3/95 CD */
/* temporary storage malloc'ed and source and sink passed as args
   11/26/01 CD */
#include "generic_includes.h"
#include "../include/field_strength.h"

#define LINK_OFFSET(dir) link_src+sizeof(su3_matrix)*(dir)
#define LINK(dir) (((su3_matrix *)F_PT(s,link_src))[dir])
#define FIELD_STRENGTH(component) (((su3_matrix *)F_PT(s,field_dest))[component])

/* Compute the field strength tensor.

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


void make_field_strength(
  field_offset link_src,       /* field offset for su3_matrix[4] type 
				  for the source link matrices */
  field_offset field_dest      /* field offset for su3_matrix[6] type
				  for the resulting field strength */
  )
{
  register int i,component,dir0=-99,dir1=-99;
  register site *s;
  int j;
  su3_matrix tmat1,tmat2;
  su3_matrix *temp1,*temp2;
  complex cc;
  msg_tag *mtag0,*mtag1;
  
  /* Allocate temporary space for two su3_matrix fields */
  temp1 = (su3_matrix *)malloc(sites_on_node*sizeof(su3_matrix));
  if(temp1 == NULL){
    printf("field_strength: No room for temp1\n");
    terminate(1);
  }
  
  temp2 = (su3_matrix *)malloc(sites_on_node*sizeof(su3_matrix));
  if(temp2 == NULL){
    printf("field_strength: No room for temp2\n");
    terminate(1);
  }
  
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
    mtag0 = start_gather_site( LINK_OFFSET(dir0), sizeof(su3_matrix),
			  dir1, EVENANDODD, gen_pt[0] );
    mtag1 = start_gather_site( LINK_OFFSET(dir1), sizeof(su3_matrix),
			  dir0, EVENANDODD, gen_pt[1] );
    
    wait_gather(mtag0);
    wait_gather(mtag1);
    FORALLSITES(i,s){
      mult_su3_nn( &LINK(dir0), 
		   (su3_matrix *)(gen_pt[1][i]), &tmat1 );
      mult_su3_na( &tmat1, (su3_matrix *)(gen_pt[0][i]), &tmat2 );
      mult_su3_na( &tmat2, &LINK(dir1), &tmat1 );
      su3_adjoint( &tmat1, &tmat2 );
      sub_su3_matrix(  &tmat1, &tmat2, &FIELD_STRENGTH(component) );
    }
    
    /**cleanup_gather(mtag0);   Use same gather in next plaquette**/
    cleanup_gather(mtag1);
    
    /* Plaquette in -dir0 +dir1 direction */
    /**mtag0 = start_gather_site( LINK_OFFSET(dir0), 
       sizeof(su3_matrix), dir1, EVENANDODD, gen_pt[0] );
       wait_gather(mtag0);  Already gathered above**/
    
    FORALLSITES(i,s){
      mult_su3_an( &LINK(dir1), 
		   &LINK(dir0), &tmat1 );
      mult_su3_an( (su3_matrix *)(gen_pt[0][i]), &tmat1, &temp1[i] );
    }
    mtag1 = start_gather_field( temp1, sizeof(su3_matrix),
				    OPP_DIR(dir0), EVENANDODD, gen_pt[1] );
    wait_gather(mtag1);
    FORALLSITES(i,s){
      mult_su3_nn( &LINK(dir1), (su3_matrix *)(gen_pt[1][i]), &tmat1 );
      su3_adjoint( &tmat1, &tmat2 );
      add_su3_matrix( &FIELD_STRENGTH(component), &tmat1,
		      &FIELD_STRENGTH(component) );
      sub_su3_matrix( &FIELD_STRENGTH(component), &tmat2,
		      &FIELD_STRENGTH(component) );
    }
    
    cleanup_gather(mtag0);
    cleanup_gather(mtag1);
    
    /* Plaquette in -dir0 -dir1 direction */
    mtag0 = start_gather_site( LINK_OFFSET(dir0), sizeof(su3_matrix),
			  OPP_DIR(dir0), EVENANDODD, gen_pt[0] );
    mtag1 = start_gather_site( LINK_OFFSET(dir1), sizeof(su3_matrix),
			  OPP_DIR(dir1), EVENANDODD, gen_pt[1] );
    wait_gather(mtag0);
    wait_gather(mtag1);
    FORALLSITES(i,s){
      mult_su3_nn( (su3_matrix *)(gen_pt[0][i]), &LINK(dir1), &temp1[i] );
      mult_su3_nn( (su3_matrix *)(gen_pt[1][i]), &LINK(dir0), &temp2[i] );
    }
    cleanup_gather(mtag0);
    cleanup_gather(mtag1);
    mtag0 = start_gather_field( temp1, sizeof(su3_matrix),
				    OPP_DIR(dir1), EVENANDODD, gen_pt[0] );
    wait_gather(mtag0);
    mtag1 = start_gather_field( temp2, sizeof(su3_matrix),
				    OPP_DIR(dir0), EVENANDODD, gen_pt[1] );
    wait_gather(mtag1);
    FORALLSITES(i,s){
      mult_su3_an( (su3_matrix *)(gen_pt[1][i]), 
		   (su3_matrix *)(gen_pt[0][i]), &tmat1 );
      su3_adjoint( &tmat1, &tmat2 );
      add_su3_matrix( &FIELD_STRENGTH(component), &tmat1,
		      &FIELD_STRENGTH(component) );
      sub_su3_matrix( &FIELD_STRENGTH(component), &tmat2,
		      &FIELD_STRENGTH(component) );
    }
    cleanup_gather(mtag0);
    cleanup_gather(mtag1);
    
    /* Plaquette in +dir0 -dir1 direction */
    mtag1 = start_gather_site( LINK_OFFSET(dir1), sizeof(su3_matrix),
			  dir0, EVENANDODD, gen_pt[1] );
    wait_gather(mtag1);
    FORALLSITES(i,s){
      mult_su3_an( &LINK(dir1), &LINK(dir0), &tmat1);
      mult_su3_nn( &tmat1, (su3_matrix *)(gen_pt[1][i]), &temp1[i] );
    }
    cleanup_gather(mtag1);
    mtag0 = start_gather_field( temp1, sizeof(su3_matrix),
				    OPP_DIR(dir1), EVENANDODD, gen_pt[0] );
    wait_gather(mtag0);
    FORALLSITES(i,s){
      mult_su3_na( (su3_matrix *)(gen_pt[0][i]), &LINK(dir0), &tmat1 );
      su3_adjoint( &tmat1, &tmat2 );
      add_su3_matrix( &FIELD_STRENGTH(component), &tmat1,
		      &FIELD_STRENGTH(component) );
      sub_su3_matrix( &FIELD_STRENGTH(component), &tmat2,
		      &FIELD_STRENGTH(component) );
    }
    cleanup_gather(mtag0);
    /* Make traceless */
    FORALLSITES(i,s){
      cc = trace_su3(&FIELD_STRENGTH(component));
      CMULREAL(cc,0.33333333333333333,cc);
      for(j=0;j<3;j++)
	CSUB(FIELD_STRENGTH(component).e[j][j],cc,
	     FIELD_STRENGTH(component).e[j][j]);
    }
  } /* component */
  
  free(temp1);
  free(temp2);
  
} /* make_field_strength() */


