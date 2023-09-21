/************************** field_strength_region.c *******************************/
/* MIMD version 7 */
/* started 8/4/95 DT */
/* removed trace 9/3/95 CD */
/* temporary storage malloc'ed and source and sink passed as args
   11/26/01 CD */
#include "wilson_flow_includes.h"
#include "../include/field_strength.h"
#include <string.h>

#define LINK_OFFSET(dir) link_src+sizeof(su3_matrix)*(dir)
#define LINK(dir) (((su3_matrix *)F_PT(s,link_src))[dir])

#define LINK_HALF_OFFSET(dir) link_half_src+sizeof(su3_matrix)*(dir)
#define LINK_HALF(dir) (((su3_matrix *)F_PT(s,link_half_src))[dir])

#define FIELD_STRENGTH(component) (((su3_matrix *)F_PT(s,field_dest))[component])
#define FIELD_STRENGTH_HALF(component) (((su3_matrix *)F_PT(s,field_half_dest))[component])
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

void clear_fieldstrength ( su3_matrix **fieldstrength, size_t ncomp ) {

  if ( fieldstrength == NULL ) {
    node0_printf("Tried to clear unallocated fieldstrength\n"); 
    terminate(1);
  }
  if ( fieldstrength[0] == NULL ) {
    node0_printf("Tried to clear unallocated fieldstrength\n"); 
    terminate(1);
  }
  memset( fieldstrength[0], '\0', ncomp * sizeof(su3_matrix) * sites_on_node );
}

void destroy_fieldstrength ( su3_matrix **fieldstrength ) {

  if ( fieldstrength == NULL ) {
    node0_printf("Tried to clear unallocated fieldstrength\n"); 
    terminate(1);
  }
  if ( fieldstrength[0] == NULL ) {
    node0_printf("Tried to free unallocated fieldstrength\n"); 
    terminate(1);
  }

  free( fieldstrength[0] );
  free( fieldstrength );
  fieldstrength = NULL;
}

su3_matrix ** new_fieldstrength ( size_t ncomp ) {

  int i;
  su3_matrix **this = (su3_matrix **)malloc( ncomp * sizeof(su3_matrix*) );
  if(this == NULL) {
    printf( "new_fieldstrength: can't malloc this\n" );
    fflush(stdout); terminate(1);
  }
  this[0] = (su3_matrix *)malloc( ncomp * sizeof(su3_matrix) * sites_on_node );
  if(this[0] == NULL) {
    printf( "new_fieldstrength: can't malloc this[0]\n" );
    fflush(stdout); terminate(1);
  }
  memset( this[0], '\0', ncomp * sizeof(su3_matrix) * sites_on_node );
  for ( i = 1; i < ncomp; i++ ) 
    this[i] = this[0] + i * sites_on_node;

  return ( this );
}


static void 
make_spatial_fieldstrength_region ( int region_flag, 
                                         su3_matrix **link, 
                                         su3_matrix **fieldstrength ) {

  register int i, j, icomp;
  register int stride = block_stride;
  #define NGATHER 2
  int dir[NGATHER] = {NODIR,NODIR}, 
      dirb[NGATHER] = {NODIR,NODIR}, 
      diro[NGATHER] = {NODIR,NODIR};
  register int ig;
  register site *s;
  msg_tag *tag[NGATHER];
  #define NTEMP_STORAGE (NGATHER)
  su3_matrix *temp[NTEMP_STORAGE];
  for( ig=0; ig<NTEMP_STORAGE; ig++ ) {
    temp[ig] = (su3_matrix *)malloc( sizeof(su3_matrix)*sites_on_node );
    if( temp[ig] == NULL ) {
      printf( "make_spatial_fieldstrength_region: can't malloc temp[%d]\n", ig );
      fflush(stdout); terminate(1);
    }
  }
  su3_matrix tmat1,tmat2;
  complex cc;

  /* Calculate the spatial/magnetic components */
  for(icomp=FS_XY;icomp<=FS_YZ;icomp++){
    switch(icomp){
      case FS_XY: dir[0]=XUP; dir[1]=YUP; break;
      case FS_XZ: dir[0]=XUP; dir[1]=ZUP; break;
      case FS_YZ: dir[0]=YUP; dir[1]=ZUP; break;
    }
    switch( block_stride ) {
      case (1):
        dirb[0] =      dir[0];  diro[0] = OPP_DIR(dir[0]);
        dirb[1] =      dir[1];  diro[1] = OPP_DIR(dir[1]);
        break;
      case (2):
        dirb[0] = DIR2(dir[0]); diro[0] = DIR2(OPP_DIR(dir[0]));
        dirb[1] = DIR2(dir[1]); diro[1] = DIR2(OPP_DIR(dir[1]));
        break;
      case (4):
        dirb[0] = DIR4(dir[0]); diro[0] = DIR4(OPP_DIR(dir[0]));
        dirb[1] = DIR4(dir[1]); diro[1] = DIR4(OPP_DIR(dir[1]));
        break;
      default:
        node0_printf("Not implemented for this BLOCKING level %d\n",block_stride );
        terminate(1);
        break;
    }
    #define LINK0 link[dir[0]]
    #define LINK1 link[dir[1]]

    /* Consistency for improved readability: gathers in 
       +- dir[0] direction have tag[0] and gen_pt[0], use temp[0] 
       +- dir[1] direction have tag[1] and gen_pt[1], use temp[1] */ 

    /* Plaquette in +dir[0] +dir[1] direction */
    ig = 0;
    // request LINK1 from direction dirb[0]
    tag[ig] = start_gather_field( LINK1 , sizeof(su3_matrix),
                                  dirb[ig], EVENANDODD, gen_pt[ig]);

    ig = 1;
    // request LINK0 from direction dirb[1]
    tag[ig] = start_gather_field( LINK0 , sizeof(su3_matrix),
                                  dirb[ig], EVENANDODD, gen_pt[ig]);

    wait_gather(tag[0]);
    wait_gather(tag[1]);
    FORALLSITES(i,s) 
    IF_BLOCKED(s, block_stride)       
    IF_REGION(s, region_flag) 
    {
      mult_su3_nn( &(LINK0[i]), (su3_matrix *)(gen_pt[0][i]), &tmat1 );
      mult_su3_na( &tmat1, (su3_matrix *)(gen_pt[1][i]), &tmat2 );
      mult_su3_na( &tmat2, &(LINK1[i]), &tmat1 );
      su3_adjoint( &tmat1, &tmat2 );
      sub_su3_matrix(  &tmat1, &tmat2, &(fieldstrength[icomp][i]) );
    }
    /**cleanup_gather(tag[1]);   Use same gather in next plaquette**/
    cleanup_gather(tag[0]);
  
    /* Plaquette in -dir[0] +dir[1] direction */
    /** tag[1] = start_gather_field( LINK0 , 
        sizeof(su3_matrix), dirb[1], EVENANDODD, gen_pt[1]);
       wait_gather(tag[1]);  Already gathered above**/
    
    FORALLSITES(i,s) 
    IF_BLOCKED(s, block_stride)       
    IF_REGION(s, region_flag) 
    {
      mult_su3_an( &(LINK1[i]), &(LINK0[i]), &tmat1 );
      mult_su3_an( (su3_matrix *)(gen_pt[1][i]), &tmat1, &temp[0][i] );
    }

    ig = 0;
    // request request dir0-dir1-dir0 staple from direction diro[0]
    tag[ig] = start_gather_field( temp[ig], sizeof(su3_matrix),
                                  diro[ig], EVENANDODD, gen_pt[ig] );
    wait_gather(tag[ig]);
 
    FORALLSITES(i,s) 
    IF_BLOCKED(s, block_stride)       
    IF_REGION(s, region_flag) 
    {
      mult_su3_nn( &(LINK1[i]), (su3_matrix *)(gen_pt[0][i]), &tmat1 );
      su3_adjoint( &tmat1, &tmat2 );
      add_su3_matrix( &(fieldstrength[icomp][i]), &tmat1,
                      &(fieldstrength[icomp][i]) );
      sub_su3_matrix( &(fieldstrength[icomp][i]), &tmat2,
                      &(fieldstrength[icomp][i]) );
    }
    cleanup_gather(tag[0]);
    cleanup_gather(tag[1]);

    /* Plaquette in -dir[0] -dir[1] direction */
    ig = 0;
    // request LINK0 from direction diro[0]
    tag[ig] = start_gather_field( LINK0 , sizeof(su3_matrix),
                                  diro[ig], EVENANDODD, gen_pt[ig]);
    ig = 1;
    // request LINK1 from direction diro[1]
    tag[ig] = start_gather_field( LINK1 , sizeof(su3_matrix),
                                  diro[ig], EVENANDODD, gen_pt[ig]);

    wait_gather(tag[0]);
    wait_gather(tag[1]);

    FORALLSITES(i,s) 
    IF_BLOCKED(s, block_stride)       
    IF_REGION(s, region_flag) 
    {
      mult_su3_nn( (su3_matrix *)(gen_pt[0][i]), &(LINK1[i]), &(temp[1][i]) );
      mult_su3_nn( (su3_matrix *)(gen_pt[1][i]), &(LINK0[i]), &(temp[0][i]) );
    }

    cleanup_gather(tag[0]);
    cleanup_gather(tag[1]);

    ig = 0;
    // request hook [-dir0,dir1] from direction diro[0]
    tag[ig] = start_gather_field( temp[ig] , sizeof(su3_matrix),
                                  diro[ig], EVENANDODD, gen_pt[ig]);
    ig = 1;
    // request hook [-dir1,dir0] from direction diro[1]
    tag[ig] = start_gather_field( temp[ig] , sizeof(su3_matrix),
                                  diro[ig], EVENANDODD, gen_pt[ig]);

    wait_gather(tag[0]);
    wait_gather(tag[1]);

    FORALLSITES(i,s) 
    IF_BLOCKED(s, block_stride)       
    IF_REGION(s, region_flag) 
    {
      mult_su3_an( (su3_matrix *)(gen_pt[0][i]), 
                   (su3_matrix *)(gen_pt[1][i]), &tmat1 );
      su3_adjoint( &tmat1, &tmat2 );
      add_su3_matrix( &(fieldstrength[icomp][i]), &tmat1,
                      &(fieldstrength[icomp][i]) );
      sub_su3_matrix( &(fieldstrength[icomp][i]), &tmat2,
                      &(fieldstrength[icomp][i]) );
    }

    cleanup_gather(tag[0]);
    cleanup_gather(tag[1]);

    /* Plaquette in +dir0 -dir1 direction */
    ig = 0;
    // request LINK1 from direction dirb[0]
    tag[ig] = start_gather_field( LINK1 , sizeof(su3_matrix),
                                  dirb[ig], EVENANDODD, gen_pt[ig]);

    wait_gather(tag[ig]);

    FORALLSITES(i,s) 
    IF_BLOCKED(s, block_stride)       
    IF_REGION(s, region_flag) 
    {
      mult_su3_an( &(LINK1[i]), &(LINK0[i]), &tmat1);
      mult_su3_nn( &tmat1, (su3_matrix *)(gen_pt[ig][i]), &temp[1][i] );
    }

    cleanup_gather(tag[ig]);

    ig = 1;
    // request LINK1 from direction diro[0]
    tag[ig] = start_gather_field( temp[ig] , sizeof(su3_matrix),
                                  diro[ig], EVENANDODD, gen_pt[ig]);

    wait_gather(tag[ig]);

    FORALLSITES(i,s) 
    IF_BLOCKED(s, block_stride)       
    IF_REGION(s, region_flag) 
    {
      mult_su3_na( (su3_matrix *)(gen_pt[ig][i]), &(LINK0[i]), &tmat1 );
      su3_adjoint( &tmat1, &tmat2 );
      add_su3_matrix( &(fieldstrength[icomp][i]), &tmat1,
                      &(fieldstrength[icomp][i]) );
      sub_su3_matrix( &(fieldstrength[icomp][i]), &tmat2,
                      &(fieldstrength[icomp][i]) );
    }
    cleanup_gather(tag[ig]);


    /* Make traceless */
    FORALLSITES(i,s) 
    IF_BLOCKED(s, block_stride)       
    IF_REGION(s, region_flag) 
    {
      cc = trace_su3( &(fieldstrength[icomp][i]) );
      CMULREAL(cc,0.33333333333333333,cc);
      for(j=0;j<3;j++)
        CSUB( fieldstrength[icomp][i].e[j][j],cc,
              fieldstrength[icomp][i].e[j][j]);
    }

    #undef LINK0
    #undef LINK1
  } /* component */

  // deallocate temporary storage
  for( ig=0; ig < NTEMP_STORAGE; ig++ ) 
    free( temp[ig] );
  #undef NTEMP_STORAGE
  #undef NGATHER
}

static void 
make_temporal_fieldstrength_full ( su3_matrix **link_s, 
                                        su3_matrix **link_t, 
                                        su3_matrix **fieldstrength ) {

  register int i, j, icomp;
  register int stride = block_stride;
  #define NGATHER 2
  int dir[NGATHER] = {NODIR,NODIR}, 
      dirb[NGATHER] = {NODIR, NODIR}, 
      diro[NGATHER] = {NODIR,NODIR};
  register int ig;
  register site *s;
  msg_tag *tag[NGATHER];
  #define NTEMP_STORAGE (NGATHER)
  su3_matrix *temp[NTEMP_STORAGE];
  for( ig=0; ig<NTEMP_STORAGE; ig++ ) {
    temp[ig] = (su3_matrix *)malloc( sizeof(su3_matrix)*sites_on_node );
    if( temp[ig] == NULL ) {
      printf( "make_temporal_fieldstrength_full: can't malloc temp[%d]\n", ig );
      fflush(stdout); terminate(1);
    }
  }
  su3_matrix tmat1,tmat2;
  complex cc;

  /* Calculate the temporal/electric components */  
  for(icomp=FS_XT;icomp<=FS_ZT;icomp++){
    switch(icomp){
      case FS_XT: dir[0]=XUP; dir[1]=TUP; break;
      case FS_YT: dir[0]=YUP; dir[1]=TUP; break;
      case FS_ZT: dir[0]=ZUP; dir[1]=TUP; break;
    }
    dirb[1] =  dir[1]; diro[1] = OPP_DIR(dir[1]);
    switch( block_stride ) {
      case (1):
        dirb[0] =      dir[0];  diro[0] = OPP_DIR(dir[0]);
        break;
      case (2):
        dirb[0] = DIR2(dir[0]); diro[0] = DIR2(OPP_DIR(dir[0]));
        break;
      case (4):
        dirb[0] = DIR4(dir[0]); diro[0] = DIR4(OPP_DIR(dir[0]));
        break;
      default:
        node0_printf("Not implemented for this BLOCKING level %d\n",block_stride );
        terminate(1);
        break;
    }
    #define LINK0 link_s[dir[0]]
    #define LINK1 link_t[0]

    /* Consistency for improved readability: gathers in 
       +- dir[0] direction have tag[0] and gen_pt[0], use temp[0] 
       +- dir[1] direction have tag[1] and gen_pt[1], use temp[1] */ 

    /* Plaquette in +dir[0] +dir[1] direction */
    ig = 0;
    // request LINK1 from direction dirb[0]
    tag[ig] = start_gather_field( LINK1 , sizeof(su3_matrix),
                                  dirb[ig], EVENANDODD, gen_pt[ig]);

    ig = 1;
    // request LINK0 from direction dirb[1]
    tag[ig] = start_gather_field( LINK0 , sizeof(su3_matrix),
                                  dirb[ig], EVENANDODD, gen_pt[ig]);

    wait_gather(tag[0]);
    wait_gather(tag[1]);
    FORALLSITES(i,s) 
    IF_BLOCKED(s, block_stride)       
    {
      mult_su3_nn( &(LINK0[i]), (su3_matrix *)(gen_pt[0][i]), &tmat1 );
      mult_su3_na( &tmat1, (su3_matrix *)(gen_pt[1][i]), &tmat2 );
      mult_su3_na( &tmat2, &(LINK1[i]), &tmat1 );
      su3_adjoint( &tmat1, &tmat2 );
      sub_su3_matrix(  &tmat1, &tmat2, &(fieldstrength[icomp][i]) );
    }
    /**cleanup_gather(tag[1]);   Use same gather in next plaquette**/
    cleanup_gather(tag[0]);
  
    /* Plaquette in -dir[0] +dir[1] direction */
    /** tag[1] = start_gather_field( LINK0, 
        sizeof(su3_matrix), dirb[1], EVENANDODD, gen_pt[1]);
       wait_gather(tag[1]);  Already gathered above**/
    
    FORALLSITES(i,s) 
    IF_BLOCKED(s, block_stride)       
    {
      mult_su3_an( &(LINK1[i]), &(LINK0[i]), &tmat1 );
      mult_su3_an( (su3_matrix *)(gen_pt[1][i]), &tmat1, &temp[0][i] );
    }

    ig = 0;
    // request request dir0-dir1-dir0 staple from direction diro[0]
    tag[ig] = start_gather_field( temp[ig], sizeof(su3_matrix),
                                  diro[ig], EVENANDODD, gen_pt[ig] );
    wait_gather(tag[ig]);
 
    FORALLSITES(i,s) 
    IF_BLOCKED(s, block_stride)       
    {
      mult_su3_nn( &(LINK1[i]), (su3_matrix *)(gen_pt[0][i]), &tmat1 );
      su3_adjoint( &tmat1, &tmat2 );
      add_su3_matrix( &(fieldstrength[icomp][i]), &tmat1,
                      &(fieldstrength[icomp][i]) );
      sub_su3_matrix( &(fieldstrength[icomp][i]), &tmat2,
                      &(fieldstrength[icomp][i]) );
    }
    cleanup_gather(tag[0]);
    cleanup_gather(tag[1]);

    /* Plaquette in -dir[0] -dir[1] direction */
    ig = 0;
    // request LINK0 from direction diro[0]
    tag[ig] = start_gather_field( LINK0 , sizeof(su3_matrix),
                                  diro[ig], EVENANDODD, gen_pt[ig]);
    ig = 1;
    // request LINK1 from direction diro[1]
    tag[ig] = start_gather_field( LINK1 , sizeof(su3_matrix),
                                  diro[ig], EVENANDODD, gen_pt[ig]);

    wait_gather(tag[0]);
    wait_gather(tag[1]);

    FORALLSITES(i,s) 
    IF_BLOCKED(s, block_stride)       
    {
      mult_su3_nn( (su3_matrix *)(gen_pt[0][i]), &(LINK1[i]), &(temp[1][i]) );
      mult_su3_nn( (su3_matrix *)(gen_pt[1][i]), &(LINK0[i]), &(temp[0][i]) );
    }

    cleanup_gather(tag[0]);
    cleanup_gather(tag[1]);

    ig = 0;
    // request hook [-dir0,dir1] from direction diro[0]
    tag[ig] = start_gather_field( temp[ig] , sizeof(su3_matrix),
                                  diro[ig], EVENANDODD, gen_pt[ig]);
    ig = 1;
    // request hook [-dir1,dir0] from direction diro[1]
    tag[ig] = start_gather_field( temp[ig] , sizeof(su3_matrix),
                                  diro[ig], EVENANDODD, gen_pt[ig]);

    wait_gather(tag[0]);
    wait_gather(tag[1]);

    FORALLSITES(i,s) 
    IF_BLOCKED(s, block_stride)       
    {
      mult_su3_an( (su3_matrix *)(gen_pt[0][i]), 
                   (su3_matrix *)(gen_pt[1][i]), &tmat1 );
      su3_adjoint( &tmat1, &tmat2 );
      add_su3_matrix( &(fieldstrength[icomp][i]), &tmat1,
                      &(fieldstrength[icomp][i]) );
      sub_su3_matrix( &(fieldstrength[icomp][i]), &tmat2,
                      &(fieldstrength[icomp][i]) );
    }

    cleanup_gather(tag[0]);
    cleanup_gather(tag[1]);

    /* Plaquette in +dir0 -dir1 direction */
    ig = 0;
    // request LINK1 from direction dirb[0]
    tag[ig] = start_gather_field( LINK1 , sizeof(su3_matrix),
                                  dirb[ig], EVENANDODD, gen_pt[ig]);

    wait_gather(tag[ig]);

    FORALLSITES(i,s) 
    IF_BLOCKED(s, block_stride)       
    {
      mult_su3_an( &(LINK1[i]), &(LINK0[i]), &tmat1);
      mult_su3_nn( &tmat1, (su3_matrix *)(gen_pt[ig][i]), &temp[1][i] );
    }

    cleanup_gather(tag[ig]);

    ig = 1;
    // request LINK1 from direction diro[0]
    tag[ig] = start_gather_field( temp[ig] , sizeof(su3_matrix),
                                  diro[ig], EVENANDODD, gen_pt[ig]);

    wait_gather(tag[ig]);

    FORALLSITES(i,s) 
    IF_BLOCKED(s, block_stride)       
    {
      mult_su3_na( (su3_matrix *)(gen_pt[ig][i]), &(LINK0[i]), &tmat1 );
      su3_adjoint( &tmat1, &tmat2 );
      add_su3_matrix( &(fieldstrength[icomp][i]), &tmat1,
                      &(fieldstrength[icomp][i]) );
      sub_su3_matrix( &(fieldstrength[icomp][i]), &tmat2,
                      &(fieldstrength[icomp][i]) );
    }
    cleanup_gather(tag[ig]);

    /* Make traceless */
    FORALLSITES(i,s) 
    IF_BLOCKED(s, block_stride)       
    {
      cc = trace_su3( &(fieldstrength[icomp][i]) );
      CMULREAL(cc,0.33333333333333333,cc);
      for(j=0;j<3;j++)
        CSUB( fieldstrength[icomp][i].e[j][j],cc,
              fieldstrength[icomp][i].e[j][j]);
    }

    #undef LINK0
    #undef LINK1
  } /* component */

  // deallocate temporary storage
  for( ig=0; ig < NTEMP_STORAGE; ig++ ) 
    free( temp[ig] );
  #undef NTEMP_STORAGE
  #undef NGATHER
}

static void 
make_temporal_fieldstrength_lwr_bulk ( su3_matrix **link_s, 
                                            su3_matrix **link_t, 
                                            su3_matrix **fieldstrength ) {

  register int i, j, icomp;
  register int stride = block_stride;
  #define NGATHER 2
  int dir[NGATHER] = {NODIR,NODIR},
      dirb[NGATHER] = {NODIR, NODIR}, 
      diro[NGATHER] = {NODIR,NODIR};
  register int ig;
  register site *s;
  msg_tag *tag[NGATHER];
  #define NTEMP_STORAGE (NGATHER)
  su3_matrix *temp[NTEMP_STORAGE];
  for( ig=0; ig<NTEMP_STORAGE; ig++ ) {
    temp[ig] = (su3_matrix *)malloc( sizeof(su3_matrix)*sites_on_node );
    if( temp[ig] == NULL ) {
      printf( "make_temporal_fieldstrength_lwr_bulk: can't malloc temp[%d]\n", ig );
      fflush(stdout); terminate(1);
    }
  }
  su3_matrix tmat1,tmat2;
  complex cc;

  /* Calculate the temporal/electric components */  
  for(icomp=FS_XT;icomp<=FS_ZT;icomp++){
    switch(icomp){
      case FS_XT: dir[0]=XUP; dir[1]=TUP; break;
      case FS_YT: dir[0]=YUP; dir[1]=TUP; break;
      case FS_ZT: dir[0]=ZUP; dir[1]=TUP; break;
    }
    dirb[1] =  dir[1]; diro[1] = OPP_DIR(dir[1]);
    switch( block_stride ) {
      case (1):
        dirb[0] =      dir[0];  diro[0] = OPP_DIR(dir[0]);
        break;
      case (2):
        dirb[0] = DIR2(dir[0]); diro[0] = DIR2(OPP_DIR(dir[0]));
        break;
      case (4):
        dirb[0] = DIR4(dir[0]); diro[0] = DIR4(OPP_DIR(dir[0]));
        break;
      default:
        node0_printf("Not implemented for this BLOCKING level %d\n",block_stride );
        terminate(1);
        break;
    }
    #define LINK0 link_s[dir[0]]
    #define LINK1 link_t[0]

    /* Consistency for improved readability: gathers in 
       +- dir[0] direction have tag[0] and gen_pt[0], use temp[0] 
       +- dir[1] direction have tag[1] and gen_pt[1], use temp[1] */ 

    /* Plaquette in +dir[0] +dir[1] direction */
    ig = 0;
    // request LINK1 from direction dirb[0]
    tag[ig] = start_gather_field( LINK1 , sizeof(su3_matrix),
                                  dirb[ig], EVENANDODD, gen_pt[ig]);

    ig = 1;
    // request LINK0 from direction dirb[1]
    tag[ig] = start_gather_field( LINK0 , sizeof(su3_matrix),
                                  dirb[ig], EVENANDODD, gen_pt[ig]);

    wait_gather(tag[0]);
    wait_gather(tag[1]);
    FORALLSITES(i,s) 
    IF_BLOCKED(s, block_stride)       
    IF_LOWER_BULK(s) 
    {
      mult_su3_nn( &(LINK0[i]), (su3_matrix *)(gen_pt[0][i]), &tmat1 );
      mult_su3_na( &tmat1, (su3_matrix *)(gen_pt[1][i]), &tmat2 );
      mult_su3_na( &tmat2, &(LINK1[i]), &tmat1 );
      su3_adjoint( &tmat1, &tmat2 );
      sub_su3_matrix(  &tmat1, &tmat2, &(fieldstrength[icomp][i]) );
    }
    /**cleanup_gather(tag[1]);   Use same gather in next plaquette**/
    cleanup_gather(tag[0]);

    /* Plaquette in -dir[0] +dir[1] direction */
    /** tag[1] = start_gather_field( LINK0, 
        sizeof(su3_matrix), dirb[1], EVENANDODD, gen_pt[1]);
       wait_gather(tag[1]);  Already gathered above**/
    
    FORALLSITES(i,s) 
    IF_BLOCKED(s, block_stride)       
    IF_LOWER_BULK(s) 
    {
      mult_su3_an( &(LINK1[i]), &(LINK0[i]), &tmat1 );
      mult_su3_an( (su3_matrix *)(gen_pt[1][i]), &tmat1, &temp[0][i] );
    }

    ig = 0;
    // request request dir0-dir1-dir0 staple from direction diro[0]
    tag[ig] = start_gather_field( temp[ig], sizeof(su3_matrix),
                                  diro[ig], EVENANDODD, gen_pt[ig] );
    wait_gather(tag[ig]);
 
    FORALLSITES(i,s) 
    IF_BLOCKED(s, block_stride)       
    IF_LOWER_BULK(s) 
    {
      mult_su3_nn( &(LINK1[i]), (su3_matrix *)(gen_pt[0][i]), &tmat1 );
      su3_adjoint( &tmat1, &tmat2 );
      add_su3_matrix( &(fieldstrength[icomp][i]), &tmat1,
                      &(fieldstrength[icomp][i]) );
      sub_su3_matrix( &(fieldstrength[icomp][i]), &tmat2,
                      &(fieldstrength[icomp][i]) );
    }
    cleanup_gather(tag[0]);
    cleanup_gather(tag[1]);
 
    /* Make traceless */
    FORALLSITES(i,s) 
    IF_BLOCKED(s, block_stride)       
    IF_LOWER_BULK(s) 
    {
      cc = trace_su3( &(fieldstrength[icomp][i]) );
      CMULREAL(cc,0.33333333333333333,cc);
      for(j=0;j<3;j++)
        CSUB( fieldstrength[icomp][i].e[j][j],cc,
              fieldstrength[icomp][i].e[j][j]);
    }

    #undef LINK0
    #undef LINK1
  } /* component */

  // deallocate temporary storage
  for( ig=0; ig < NTEMP_STORAGE; ig++ ) 
    free( temp[ig] );
  #undef NTEMP_STORAGE
  #undef NGATHER
}

static void 
make_dropped_temporal_fieldstrength_lwr_bulk ( su3_matrix **link, 
                                                    su3_matrix **fieldstrength ) {

  register int i, j, icomp;
  register int stride = block_stride;
  #define NGATHER 2
  int dir[NGATHER] = {NODIR,NODIR}, 
      dirb[NGATHER] = {NODIR, NODIR}, 
      diro[NGATHER] = {NODIR,NODIR};
  register int ig;
  register site *s;
  msg_tag *tag[NGATHER];
  #define NTEMP_STORAGE (NGATHER)
  su3_matrix *temp[NTEMP_STORAGE];
  for( ig=0; ig<NTEMP_STORAGE; ig++ ) {
    temp[ig] = (su3_matrix *)malloc( sizeof(su3_matrix)*sites_on_node );
    if( temp[ig] == NULL ) {
      printf( "make_dropped_temporal_fieldstrength_lwr_bulk: can't malloc temp[%d]\n", ig );
      fflush(stdout); terminate(1);
    }
  }
  su3_matrix tmat1,tmat2;
  complex cc;

  /* Calculate the temporal/electric components */  
  for(icomp=FS_XT;icomp<=FS_ZT;icomp++){
    switch(icomp){
      case FS_XT: dir[0]=XUP; dir[1]=TUP; break;
      case FS_YT: dir[0]=YUP; dir[1]=TUP; break;
      case FS_ZT: dir[0]=ZUP; dir[1]=TUP; break;
    }
    dirb[1] =  dir[1]; diro[1] = OPP_DIR(dir[1]);
    switch( block_stride ) {
      case (1):
        dirb[0] =      dir[0];  diro[0] = OPP_DIR(dir[0]);
        break;
      case (2):
        dirb[0] = DIR2(dir[0]); diro[0] = DIR2(OPP_DIR(dir[0]));
        break;
      case (4):
        dirb[0] = DIR4(dir[0]); diro[0] = DIR4(OPP_DIR(dir[0]));
        break;
      default:
        node0_printf("Not implemented for this BLOCKING level %d\n",block_stride );
        terminate(1);
        break;
    }
    #define LINK0 link[dir[0]]

    /* Consistency for improved readability: gathers in 
       +- dir[0] direction have tag[0] and gen_pt[0], use temp[0] 
       +- dir[1] direction have tag[1] and gen_pt[1], use temp[1] */ 

    /* Plaquette in +dir[0] +dir[1] direction */
    ig = 1;
    // request LINK0 from direction dirb[1]
    tag[ig] = start_gather_field( LINK0 , sizeof(su3_matrix),
                                  dirb[ig], EVENANDODD, gen_pt[ig]);

    wait_gather(tag[1]);
    FORALLSITES(i,s) 
    IF_BLOCKED(s, block_stride)       
    IF_LOWER_BULK(s) 
    {
      mult_su3_na( &(LINK0[i]), (su3_matrix *)(gen_pt[1][i]), &tmat1 );
      su3_adjoint( &tmat1, &tmat2 );
      sub_su3_matrix(  &tmat1, &tmat2, &(fieldstrength[icomp][i]) );
    }
    /**cleanup_gather(tag[1]);   Use same gather in next plaquette**/

    /* Plaquette in -dir[0] +dir[1] direction */
    /** tag[1] = start_gather_field( LINK0 , 
        sizeof(su3_matrix), dirb[1], EVENANDODD, gen_pt[1]);
       wait_gather(tag[1]);  Already gathered above**/
    
    FORALLSITES(i,s) 
    IF_BLOCKED(s, block_stride)       
    IF_LOWER_BULK(s) 
      mult_su3_an( (su3_matrix *)(gen_pt[1][i]), &(LINK0[i]), &temp[0][i] );

    ig = 0;
    // request request dir0-dir1-dir0 staple from direction diro[0]
    tag[ig] = start_gather_field( temp[ig], sizeof(su3_matrix),
                                  diro[ig], EVENANDODD, gen_pt[ig] );
    wait_gather(tag[ig]);
 
    FORALLSITES(i,s) 
    IF_BLOCKED(s, block_stride)       
    IF_LOWER_BULK(s) 
    {
      su3_adjoint( (su3_matrix *)(gen_pt[0][i]), &tmat2 );
      add_su3_matrix( &(fieldstrength[icomp][i]), (su3_matrix *)(gen_pt[0][i]),
                      &(fieldstrength[icomp][i]) );
      sub_su3_matrix( &(fieldstrength[icomp][i]), &tmat2,
                      &(fieldstrength[icomp][i]) );
    }
    cleanup_gather(tag[0]);
    cleanup_gather(tag[1]);

    /* Make traceless */
    FORALLSITES(i,s) 
    IF_BLOCKED(s, block_stride)       
    IF_LOWER_BULK(s) 
    {
      cc = trace_su3( &(fieldstrength[icomp][i]) );
      CMULREAL(cc,0.33333333333333333,cc);
      for(j=0;j<3;j++)
        CSUB( fieldstrength[icomp][i].e[j][j],cc,
              fieldstrength[icomp][i].e[j][j]);
    }

    #undef LINK0
  } /* component */

  // deallocate temporary storage
  for( ig=0; ig < NTEMP_STORAGE; ig++ ) 
    free( temp[ig] );
  #undef NTEMP_STORAGE
  #undef NGATHER
}

static void 
make_dropped_temporal_fieldstrength_half ( su3_matrix **link, 
                                                su3_matrix **link_half, 
                                                su3_matrix **fieldstrength ) {

  register int i, j, icomp, jcomp;
  register int stride = block_stride;
  #define NGATHER 2
  int dir[NGATHER] = {NODIR,NODIR}, 
      dirb[NGATHER] = {NODIR, NODIR}, 
      diro[NGATHER] = {NODIR,NODIR};
  register int ig;
  register site *s;
  msg_tag *tag[NGATHER];
  #define NTEMP_STORAGE (NGATHER)
  su3_matrix *temp[NTEMP_STORAGE];
  for( ig=0; ig<NTEMP_STORAGE; ig++ ) {
    temp[ig] = (su3_matrix *)malloc( sizeof(su3_matrix)*sites_on_node );
    if( temp[ig] == NULL ) {
      printf( "make_spatial_fieldstrength_region: can't malloc temp[%d]\n", ig );
      fflush(stdout); terminate(1);
    }
  }
  su3_matrix tmat1,tmat2;
  complex cc;

  /* Calculate the temporal/electric components */  
  for(icomp=FS_XT;icomp<=FS_ZT;icomp++){
    switch(icomp){
      case FS_XT: dir[0]=XUP; dir[1]=TUP; break;
      case FS_YT: dir[0]=YUP; dir[1]=TUP; break;
      case FS_ZT: dir[0]=ZUP; dir[1]=TUP; break;
    }
    jcomp = icomp + 3;
    dirb[1] =  dir[1]; diro[1] = OPP_DIR(dir[1]);
    switch( block_stride ) {
      case (1):
        dirb[0] =      dir[0];  diro[0] = OPP_DIR(dir[0]);
        break;
      case (2):
        dirb[0] = DIR2(dir[0]); diro[0] = DIR2(OPP_DIR(dir[0]));
        break;
      case (4):
        dirb[0] = DIR4(dir[0]); diro[0] = DIR4(OPP_DIR(dir[0]));
        break;
      default:
        node0_printf("Not implemented for this BLOCKING level %d\n",block_stride );
        terminate(1);
        break;
    }
    #define LINKF link[dir[0]]
    #define LINKH link_half[dir[0]]

    /* Consistency for improved readability: gathers in 
       +- dir[0] direction have tag[0] and gen_pt[0], use temp[0] 
       +- dir[1] direction have tag[1] and gen_pt[1], use temp[1] */ 

    /* Plaquette in +dir[0] +dir[1] direction, first at integer, then half-integer time */
    ig = 1;
    // request LINKF from direction dirb[1]
    tag[ig] = start_gather_field( LINKF , sizeof(su3_matrix),
                                  dirb[ig], EVENANDODD, gen_pt[ig]);

    FORALLSITES(i,s) 
    IF_BLOCKED(s, block_stride)       
    IF_LOWER_BULK(s) 
    {
      mult_su3_na( &(LINKF[i]), &(LINKH[i]), &tmat1 );
      mult_su3_an( &(LINKH[i]), &(LINKF[i]), &(temp[0][i]) );
      su3_adjoint( &tmat1, &tmat2 );
      sub_su3_matrix(  &tmat1, &tmat2, &(fieldstrength[icomp][i]) );
    }

    ig = 0;
    // request ccw PLAQ from direction diro[0]
    tag[ig] = start_gather_field( temp[0] , sizeof(su3_matrix),
                                  diro[ig], EVENANDODD, gen_pt[ig]);

    wait_gather(tag[1]);

    FORALLSITES(i,s) 
    IF_BLOCKED(s, block_stride)       
    IF_LOWER_BULK(s) 
    {
      mult_su3_na( &(LINKH[i]), (su3_matrix *)(gen_pt[1][i]), &tmat1 );
      su3_adjoint( &tmat1, &tmat2 );
      sub_su3_matrix(  &tmat1, &tmat2, &(fieldstrength[jcomp][i]) );
    }
    /**cleanup_gather(tag[1]);   Use same gather in next plaquette**/

    /* Plaquette in -dir[0] +dir[1] direction */
    /** tag[1] = start_gather_field( LINKF , 
        sizeof(su3_matrix), dirb[1], EVENANDODD, gen_pt[1]);
       wait_gather(tag[1]);  Already gathered above**/
    
    wait_gather(tag[0]);

    FORALLSITES(i,s) 
    IF_BLOCKED(s, block_stride)       
    IF_LOWER_BULK(s) 
    {
      su3_adjoint( (su3_matrix *)(gen_pt[0][i]), &tmat2 );
      add_su3_matrix( &(fieldstrength[icomp][i]), (su3_matrix *)(gen_pt[0][i]),
                      &(fieldstrength[icomp][i]) );
      sub_su3_matrix( &(fieldstrength[icomp][i]), &tmat2,
                      &(fieldstrength[icomp][i]) );
    }
    cleanup_gather(tag[0]);

    FORALLSITES(i,s) 
    IF_BLOCKED(s, block_stride)       
    IF_LOWER_BULK(s) 
      mult_su3_an( (su3_matrix *)(gen_pt[1][i]), &(LINKH[i]), &(temp[0][i]) );

    ig = 0;
    // request request dir0-dir1-dir0 staple from direction diro[0]
    tag[ig] = start_gather_field( temp[ig], sizeof(su3_matrix),
                                  diro[ig], EVENANDODD, gen_pt[ig] );
    wait_gather(tag[ig]);
 
    FORALLSITES(i,s) 
    IF_BLOCKED(s, block_stride)       
    IF_LOWER_BULK(s) 
    {
      su3_adjoint( (su3_matrix *)(gen_pt[0][i]), &tmat2 );
      add_su3_matrix( &(fieldstrength[jcomp][i]), (su3_matrix *)(gen_pt[0][i]),
                      &(fieldstrength[jcomp][i]) );
      sub_su3_matrix( &(fieldstrength[jcomp][i]), &tmat2,
                      &(fieldstrength[jcomp][i]) );
    }
    cleanup_gather(tag[0]);
    cleanup_gather(tag[1]);

    /* Make traceless */
    FORALLSITES(i,s) 
    IF_BLOCKED(s, block_stride)       
    IF_LOWER_BULK(s) 
    {
      cc = trace_su3( &(fieldstrength[icomp][i]) );
      CMULREAL(cc,0.33333333333333333,cc);
      for(j=0;j<3;j++)
        CSUB( fieldstrength[icomp][i].e[j][j],cc,
              fieldstrength[icomp][i].e[j][j]);

      cc = trace_su3( &(fieldstrength[jcomp][i]) );
      CMULREAL(cc,0.33333333333333333,cc);
      for(j=0;j<3;j++)
        CSUB( fieldstrength[jcomp][i].e[j][j],cc,
              fieldstrength[jcomp][i].e[j][j]);
    }

    #undef LINKF
    #undef LINKH
  } /* component */

  // deallocate temporary storage
  for( ig=0; ig < NTEMP_STORAGE; ig++ ) 
    free( temp[ig] );
  #undef NTEMP_STORAGE
  #undef NGATHER
}

static void 
make_dropped_temporal_fieldstrength_bdry ( su3_matrix **link, 
                                        su3_matrix ***link_lf, 
                                        su3_matrix **fieldstrength ) {

  register int i, j, icomp;
  register int stride = block_stride;
  #define NGATHER 1
  int dir[NGATHER] = {NODIR}, 
      dirb[NGATHER] = {NODIR}, 
      diro[NGATHER] = {NODIR};
  register int ig;
  register site *s;
  msg_tag *tag[NGATHER];
  #define NTEMP_STORAGE (NGATHER)
  su3_matrix *temp[NTEMP_STORAGE];
  for( ig=0; ig<NTEMP_STORAGE; ig++ ) {
    temp[ig] = (su3_matrix *)malloc( sizeof(su3_matrix)*sites_on_node );
    if( temp[ig] == NULL ) {
      printf( "make_spatial_fieldstrength_region: can't malloc temp[%d]\n", ig );
      fflush(stdout); terminate(1);
    }
  }
  su3_matrix tmat1,tmat2;
  complex cc;

  /* Calculate the temporal/electric components */  
  for(icomp=FS_XT;icomp<=FS_ZT;icomp++){
    switch(icomp){
      case FS_XT: dir[0]=XUP; dir[1]=TUP; break;
      case FS_YT: dir[0]=YUP; dir[1]=TUP; break;
      case FS_ZT: dir[0]=ZUP; dir[1]=TUP; break;
    }
    dirb[1] =  dir[1]; diro[1] = OPP_DIR(dir[1]);
    switch( block_stride ) {
      case (1):
        dirb[0] =      dir[0];  diro[0] = OPP_DIR(dir[0]);
        break;
      case (2):
        dirb[0] = DIR2(dir[0]); diro[0] = DIR2(OPP_DIR(dir[0]));
        break;
      case (4):
        dirb[0] = DIR4(dir[0]); diro[0] = DIR4(OPP_DIR(dir[0]));
        break;
      default:
        node0_printf("Not implemented for this BLOCKING level %d\n",block_stride );
        terminate(1);
        break;
    }
    #define LINK0 link[dir[0]]
    #define LINK1 link_lf[0][dir[0]]
    #define LINK2 link_lf[1][dir[0]]

    /* Consistency for improved readability: gathers in 
       +- dir[0] direction have tag[0] and gen_pt[0], use temp[0] 
       +- dir[1] direction have tag[1] and gen_pt[1], use temp[1] */ 

    /* Plaquette in +dir[0] +dir[1] direction */
    FORALLSITES(i,s) 
    IF_BLOCKED(s, block_stride)       
    IF_BOUNDARY(s) 
    {
      mult_su3_na( &(LINK1[i]), &(LINK0[i]), &tmat1 );
      su3_adjoint( &tmat1, &tmat2 );
      sub_su3_matrix(  &tmat1, &tmat2, &(fieldstrength[icomp][i]) );
    }

    /* Plaquette in -dir[0] +dir[1] direction */
    
    FORALLSITES(i,s) 
    IF_BLOCKED(s, block_stride)       
    IF_BOUNDARY(s) 
      mult_su3_an(  &(LINK0[i]), &(LINK1[i]), &temp[0][i] );

    ig = 0;
    // request request dir0-dir1-dir0 staple from direction diro[0]
    tag[ig] = start_gather_field( temp[ig], sizeof(su3_matrix),
                                  diro[ig], EVENANDODD, gen_pt[ig] );
    wait_gather(tag[ig]);
 
    FORALLSITES(i,s) 
    IF_BLOCKED(s, block_stride)       
    IF_BOUNDARY(s) 
    {
      su3_adjoint( (su3_matrix *)(gen_pt[0][i]), &tmat2 );
      add_su3_matrix( &(fieldstrength[icomp][i]), (su3_matrix *)(gen_pt[0][i]),
                      &(fieldstrength[icomp][i]) );
      sub_su3_matrix( &(fieldstrength[icomp][i]), &tmat2,
                      &(fieldstrength[icomp][i]) );
    }
    cleanup_gather(tag[0]);

    /* Consistency for improved readability: gathers in 
       +- dir[0] direction have tag[0] and gen_pt[0], use temp[0] 
       +- dir[1] direction have tag[1] and gen_pt[1], use temp[1] */ 

    /* Plaquette in +dir[0] -dir[1] direction */
    FORALLSITES(i,s) 
    IF_BLOCKED(s, block_stride)       
    IF_BOUNDARY(s) 
    {
      mult_su3_na( &(LINK2[i]), &(LINK1[i]), &tmat1 );
      su3_adjoint( &tmat1, &tmat2 );
      add_su3_matrix( &(fieldstrength[icomp][i]), &tmat1,
                      &(fieldstrength[icomp][i]) );
      sub_su3_matrix( &(fieldstrength[icomp][i]), &tmat2,
                      &(fieldstrength[icomp][i]) );
    }

    /* Plaquette in -dir[0] -dir[1] direction */
    
    FORALLSITES(i,s) 
    IF_BLOCKED(s, block_stride)       
    IF_BOUNDARY(s) 
      mult_su3_an(  &(LINK1[i]), &(LINK2[i]), &temp[0][i] );

    ig = 0;
    // request request dir0-dir1-dir0 staple from direction diro[0]
    tag[ig] = start_gather_field( temp[ig], sizeof(su3_matrix),
                                  diro[ig], EVENANDODD, gen_pt[ig] );
    wait_gather(tag[ig]);
 
    FORALLSITES(i,s) 
    IF_BLOCKED(s, block_stride)       
    IF_BOUNDARY(s) 
    {
      su3_adjoint( (su3_matrix *)(gen_pt[0][i]), &tmat2 );
      add_su3_matrix( &(fieldstrength[icomp][i]), (su3_matrix *)(gen_pt[0][i]),
                      &(fieldstrength[icomp][i]) );
      sub_su3_matrix( &(fieldstrength[icomp][i]), &tmat2,
                      &(fieldstrength[icomp][i]) );
    }
    cleanup_gather(tag[0]);

    /* Make traceless */
    FORALLSITES(i,s) 
    IF_BLOCKED(s, block_stride)       
    IF_BOUNDARY(s) 
    {
      cc = trace_su3( &(fieldstrength[icomp][i]) );
      CMULREAL(cc,0.33333333333333333,cc);
      for(j=0;j<3;j++)
        CSUB( fieldstrength[icomp][i].e[j][j],cc,
              fieldstrength[icomp][i].e[j][j]);

    }

    #undef LINK0
    #undef LINK1
    #undef LINK2
  } /* component */

  // deallocate temporary storage
  for( ig=0; ig < NTEMP_STORAGE; ig++ ) 
    free( temp[ig] );
  #undef NTEMP_STORAGE
  #undef NGATHER
}


void 
make_fieldstrength_region ( int region_flag, 
                            field_offset link_src, 
                            su3_matrix **fieldstrength ) {

  register int i, dir;
  register site *s;
  #define NLINK 4
  su3_matrix **link = new_links( region_flag, NLINK );

  make_spatial_fieldstrength_region( region_flag, link, fieldstrength );

  if ( region_flag == FULLVOL ) {
    make_temporal_fieldstrength_full ( link, &(link[TUP]), fieldstrength );
  }
  if ( region_flag == LOWER_BULK ) {
    #ifdef DROP_TIME_LINKS
    make_dropped_temporal_fieldstrength_lwr_bulk ( link, fieldstrength );
    #else
    make_temporal_fieldstrength_lwr_bulk ( link, &(link[TUP]), fieldstrength );
    #endif
  }

  // deallocate temporary spatial links
  destroy_links( link );
  #undef NLINK
}

void 
make_fieldstrength_bdry ( field_offset link_src, 
                          su3_matrix ***link_last_flow, 
                          su3_matrix **fieldstrength ) {

  register int i, dir;
  register site *s;
  #define NLINK 3
  su3_matrix **link = new_links( ACTIVE, NLINK );

  make_spatial_fieldstrength_region( BOUNDARY, link_last_flow[0], fieldstrength );

  make_dropped_temporal_fieldstrength_bdry ( link, link_last_flow, fieldstrength );

  // deallocate temporary spatial links
  destroy_links( link );
  #undef NLINK
}

void 
make_fieldstrength_bulk ( field_offset link_src, 
                          su3_matrix **fieldstrength ) {

  register int i, dir;
  register site *s;
  #ifndef DROP_TIME_LINKS
  #define NLINK 4
  #else
  #define NLINK 3
  #endif
  su3_matrix **link = new_links( ACTIVE, NLINK );

  make_spatial_fieldstrength_region( ACTIVE, link, fieldstrength );

  // #undef DROP_TIME_LINKS
  #ifdef DROP_TIME_LINKS
  make_dropped_temporal_fieldstrength_lwr_bulk ( link, fieldstrength );
  #else
  make_temporal_fieldstrength_lwr_bulk ( link, &(link[TUP]), fieldstrength );
  #endif

  // deallocate temporary spatial links
  destroy_links( link );
  #undef NLINK
}

void 
make_fieldstrength_half ( field_offset link_src, 
                          su3_matrix **fieldstrength ) {

  register int i, dir;
  register site *s;
  #define NLINK 3
  su3_matrix **link = new_links( ACTIVE, NLINK );
  su3_matrix **link_half = new_half_links();

  make_spatial_fieldstrength_region( LOWER_BULK, link_half, fieldstrength );

  make_dropped_temporal_fieldstrength_half ( link, link_half, fieldstrength );

  // deallocate temporary spatial links
  destroy_links( link_half );
  destroy_links( link );
  #undef NLINK
}

void 
make_fieldstrength_full ( field_offset link_src, 
                          su3_matrix **fieldstrength ) {

  register int i, dir;
  register site *s;
  #define NLINK 4
  su3_matrix **link = new_links( FULLVOL, NLINK );

  make_spatial_fieldstrength_region( FULLVOL, link, fieldstrength );

  make_temporal_fieldstrength_full ( link, &(link[TUP]), fieldstrength );

  // deallocate temporary spatial links
  destroy_links( link );
  #undef NLINK
}

///////////////////////////////////

void make_field_strength_bulk   ( 
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

  /* Calculate the spatial/magnetic components at integer time-steps */
  for(component=FS_XY;component<=FS_YZ;component++){
    switch(component){
    case FS_XY: dir0=XUP; dir1=YUP; break;
    case FS_XZ: dir0=XUP; dir1=ZUP; break;
    case FS_YZ: dir0=YUP; dir1=ZUP; break;
    }
    
    /* Plaquette in +dir0 +dir1 direction, counted towards integer time */
    mtag0 = start_gather_site( LINK_OFFSET(dir0), sizeof(su3_matrix),
                          dir1, EVENANDODD, gen_pt[0] );
    mtag1 = start_gather_site( LINK_OFFSET(dir1), sizeof(su3_matrix),
                          dir0, EVENANDODD, gen_pt[1] );
    
    wait_gather(mtag0);
    wait_gather(mtag1);
    FORALLSITES(i,s){ 
    IF_ACTIVE(s) {
      mult_su3_nn( &LINK(dir0), 
                   (su3_matrix *)(gen_pt[1][i]), &tmat1 );
      mult_su3_na( &tmat1, (su3_matrix *)(gen_pt[0][i]), &tmat2 );
      mult_su3_na( &tmat2, &LINK(dir1), &tmat1 );
      su3_adjoint( &tmat1, &tmat2 );
      sub_su3_matrix(  &tmat1, &tmat2, &FIELD_STRENGTH(component) );
    }
    }
    /**cleanup_gather(mtag0);   Use same gather in next plaquette**/
    cleanup_gather(mtag1);
    
    /* Plaquette in -dir0 +dir1 direction, counted towards integer time */
    /**mtag0 = start_gather_site( LINK_OFFSET(dir0), 
       sizeof(su3_matrix), dir1, EVENANDODD, gen_pt[0] );
       wait_gather(mtag0);  Already gathered above**/
    
    FORALLSITES(i,s){
    IF_ACTIVE(s) {
      mult_su3_an( &LINK(dir1), 
                   &LINK(dir0), &tmat1 );
      mult_su3_an( (su3_matrix *)(gen_pt[0][i]), &tmat1, &temp1[i] );
    }
    }
    mtag1 = start_gather_field( temp1, sizeof(su3_matrix),
                                    OPP_DIR(dir0), EVENANDODD, gen_pt[1] );
    wait_gather(mtag1);
    FORALLSITES(i,s){
    IF_ACTIVE(s) {
      mult_su3_nn( &LINK(dir1), (su3_matrix *)(gen_pt[1][i]), &tmat1 );
      su3_adjoint( &tmat1, &tmat2 );
      add_su3_matrix( &FIELD_STRENGTH(component), &tmat1,
                      &FIELD_STRENGTH(component) );
      sub_su3_matrix( &FIELD_STRENGTH(component), &tmat2,
                      &FIELD_STRENGTH(component) );
    }
    }
    cleanup_gather(mtag0);
    cleanup_gather(mtag1);

    /* Plaquette in -dir0 -dir1 direction, counted towards integer time */
    mtag0 = start_gather_site( LINK_OFFSET(dir0), sizeof(su3_matrix),
                          OPP_DIR(dir0), EVENANDODD, gen_pt[0] );
    mtag1 = start_gather_site( LINK_OFFSET(dir1), sizeof(su3_matrix),
                          OPP_DIR(dir1), EVENANDODD, gen_pt[1] );
    wait_gather(mtag0);
    wait_gather(mtag1);
    FORALLSITES(i,s){
    IF_ACTIVE(s) {
      mult_su3_nn( (su3_matrix *)(gen_pt[0][i]), &LINK(dir1), &temp1[i] );
      mult_su3_nn( (su3_matrix *)(gen_pt[1][i]), &LINK(dir0), &temp2[i] );
    }
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
    IF_ACTIVE(s) {
      mult_su3_an( (su3_matrix *)(gen_pt[1][i]), 
                   (su3_matrix *)(gen_pt[0][i]), &tmat1 );

      su3_adjoint( &tmat1, &tmat2 );
      add_su3_matrix( &FIELD_STRENGTH(component), &tmat1,
                      &FIELD_STRENGTH(component) );
      sub_su3_matrix( &FIELD_STRENGTH(component), &tmat2,
                      &FIELD_STRENGTH(component) );
    }
    }
    cleanup_gather(mtag0);
    cleanup_gather(mtag1);

    /* Plaquette in +dir0 -dir1 direction, counted towards integer time */
    mtag1 = start_gather_site( LINK_OFFSET(dir1), sizeof(su3_matrix),
                          dir0, EVENANDODD, gen_pt[1] );
    wait_gather(mtag1);
    FORALLSITES(i,s){
    IF_ACTIVE(s) {
      mult_su3_an( &LINK(dir1), &LINK(dir0), &tmat1);
      mult_su3_nn( &tmat1, (su3_matrix *)(gen_pt[1][i]), &temp1[i] );
    }
    }
    cleanup_gather(mtag1);
    mtag0 = start_gather_field( temp1, sizeof(su3_matrix),
                                    OPP_DIR(dir1), EVENANDODD, gen_pt[0] );
    wait_gather(mtag0);
    FORALLSITES(i,s){
    IF_ACTIVE(s) {
      mult_su3_na( (su3_matrix *)(gen_pt[0][i]), &LINK(dir0), &tmat1 );
      su3_adjoint( &tmat1, &tmat2 );
      add_su3_matrix( &FIELD_STRENGTH(component), &tmat1,
                      &FIELD_STRENGTH(component) );
      sub_su3_matrix( &FIELD_STRENGTH(component), &tmat2,
                      &FIELD_STRENGTH(component) );
    }
    }
    cleanup_gather(mtag0);

    /* Make traceless */
    FORALLSITES(i,s)
    IF_ACTIVE(s) 
    {
      cc = trace_su3(&FIELD_STRENGTH(component));
      CMULREAL(cc,0.33333333333333333,cc);
      for(j=0;j<3;j++)
        CSUB(FIELD_STRENGTH(component).e[j][j],cc,
             FIELD_STRENGTH(component).e[j][j]);

    }

    // if (component == FS_YZ) dumpmat( &(lattice[0].fieldstrength[component]) );

  } /* component */

  /* Calculate the temporal/electric components at integer time-steps */  
  for(component=FS_XT;component<=FS_ZT;component++){
    switch(component){
    case FS_XT: dir0=XUP; dir1=TUP; break;
    case FS_YT: dir0=YUP; dir1=TUP; break;
    case FS_ZT: dir0=ZUP; dir1=TUP; break;
    }
    memset( &FIELD_STRENGTH(component), '\0',sizeof(su3_matrix) );
    
    #ifndef DROP_TIME_LINKS
    /* Plaquette in +dir0 +dir1 direction, counted towards integer time */
    mtag0 = start_gather_site( LINK_OFFSET(dir0), sizeof(su3_matrix),
                          dir1, EVENANDODD, gen_pt[0] );
    mtag1 = start_gather_site( LINK_OFFSET(dir1), sizeof(su3_matrix),
                          dir0, EVENANDODD, gen_pt[1] );
    wait_gather(mtag0);
    wait_gather(mtag1);
    FORALLSITES(i,s){ 
    IF_ACTIVE(s) {
      mult_su3_nn( &LINK(dir0), 
                   (su3_matrix *)(gen_pt[1][i]), &tmat1 );
      mult_su3_na( &tmat1, (su3_matrix *)(gen_pt[0][i]), &tmat2 );
      mult_su3_na( &tmat2, &LINK(dir1), &tmat1 );
      su3_adjoint( &tmat1, &tmat2 );
      sub_su3_matrix(  &tmat1, &tmat2, &FIELD_STRENGTH(component) );
    }
    }
    /**cleanup_gather(mtag0);   Use same gather in next plaquette**/
    cleanup_gather(mtag1);
    
    /* Plaquette in -dir0 +dir1 direction, counted towards integer time */
    /**mtag0 = start_gather_site( LINK_OFFSET(dir0), 
       sizeof(su3_matrix), dir1, EVENANDODD, gen_pt[0] );
       wait_gather(mtag0);  Already gathered above**/
    
    FORALLSITES(i,s){
    IF_ACTIVE(s) {
      mult_su3_an( &LINK(dir1), 
                   &LINK(dir0), &tmat1 );
      mult_su3_an( (su3_matrix *)(gen_pt[0][i]), &tmat1, &temp1[i] );
    }
    }
    mtag1 = start_gather_field( temp1, sizeof(su3_matrix),
                                    OPP_DIR(dir0), EVENANDODD, gen_pt[1] );
    wait_gather(mtag1);
    FORALLSITES(i,s){
    IF_ACTIVE(s) {
      mult_su3_nn( &LINK(dir1), (su3_matrix *)(gen_pt[1][i]), &tmat1 );
      su3_adjoint( &tmat1, &tmat2 );
      add_su3_matrix( &FIELD_STRENGTH(component), &tmat1,
                      &FIELD_STRENGTH(component) );
      sub_su3_matrix( &FIELD_STRENGTH(component), &tmat2,
                      &FIELD_STRENGTH(component) );
    }
    }
    cleanup_gather(mtag0);
    cleanup_gather(mtag1);
    
    // /* Plaquette in -dir0 -dir1 direction, not counted towards integer time */
    // mtag0 = start_gather_site( LINK_OFFSET(dir0), sizeof(su3_matrix),
    //                       OPP_DIR(dir0), EVENANDODD, gen_pt[0] );
    // mtag1 = start_gather_site( LINK_OFFSET(dir1), sizeof(su3_matrix),
    //                       OPP_DIR(dir1), EVENANDODD, gen_pt[1] );
    // wait_gather(mtag0);
    // wait_gather(mtag1);
    // FORALLSITES(i,s){
    // IF_ACTIVE(s) {
    //   mult_su3_nn( (su3_matrix *)(gen_pt[0][i]), &LINK(dir1), &temp1[i] );
    //   mult_su3_nn( (su3_matrix *)(gen_pt[1][i]), &LINK(dir0), &temp2[i] );
    // }
    // }
    // cleanup_gather(mtag0);
    // cleanup_gather(mtag1);
    // mtag0 = start_gather_field( temp1, sizeof(su3_matrix),
    //                                 OPP_DIR(dir1), EVENANDODD, gen_pt[0] );
    // wait_gather(mtag0);
    // mtag1 = start_gather_field( temp2, sizeof(su3_matrix),
    //                                 OPP_DIR(dir0), EVENANDODD, gen_pt[1] );
    // wait_gather(mtag1);
    // FORALLSITES(i,s){
    // IF_ACTIVE(s) {
    //   mult_su3_an( (su3_matrix *)(gen_pt[1][i]), 
    //                (su3_matrix *)(gen_pt[0][i]), &tmat1 );
    //   su3_adjoint( &tmat1, &tmat2 );
    //   add_su3_matrix( &FIELD_STRENGTH(component), &tmat1,
    //                   &FIELD_STRENGTH(component) );
    //   sub_su3_matrix( &FIELD_STRENGTH(component), &tmat2,
    //                   &FIELD_STRENGTH(component) );
    // }
    // }
    // cleanup_gather(mtag0);
    // cleanup_gather(mtag1);
    
    // /* Plaquette in +dir0 -dir1 direction, not counted towards integer time */
    // mtag1 = start_gather_site( LINK_OFFSET(dir1), sizeof(su3_matrix),
    //                       dir0, EVENANDODD, gen_pt[1] );
    // wait_gather(mtag1);
    // FORALLSITES(i,s){
    // IF_ACTIVE(s) {
    //   mult_su3_an( &LINK(dir1), &LINK(dir0), &tmat1);
    //   mult_su3_nn( &tmat1, (su3_matrix *)(gen_pt[1][i]), &temp1[i] );
    // }
    // }
    // cleanup_gather(mtag1);
    // mtag0 = start_gather_field( temp1, sizeof(su3_matrix),
    //                                 OPP_DIR(dir1), EVENANDODD, gen_pt[0] );
    // wait_gather(mtag0);
    // FORALLSITES(i,s){
    // IF_ACTIVE(s) {
    //   mult_su3_na( (su3_matrix *)(gen_pt[0][i]), &LINK(dir0), &tmat1 );
    //   su3_adjoint( &tmat1, &tmat2 );
    //   add_su3_matrix( &FIELD_STRENGTH(component), &tmat1,
    //                   &FIELD_STRENGTH(component) );
    //   sub_su3_matrix( &FIELD_STRENGTH(component), &tmat2,
    //                   &FIELD_STRENGTH(component) );
    // }
    // }
    // cleanup_gather(mtag0);
    #else
    /* Plaquette in +dir0 +dir1 direction, counted towards integer time */
    mtag0 = start_gather_site( LINK_OFFSET(dir0), sizeof(su3_matrix),
                          dir1, EVENANDODD, gen_pt[0] );
    wait_gather(mtag0);
    FORALLSITES(i,s){ 
    IF_ACTIVE(s) {
      mult_su3_na( &LINK(dir0), (su3_matrix *)(gen_pt[0][i]), &tmat1 );
      su3_adjoint( &tmat1, &tmat2 );
      sub_su3_matrix(  &tmat1, &tmat2, &FIELD_STRENGTH(component) );
    }
    }
    /**cleanup_gather(mtag0);   Use same gather in next plaquette**/
    cleanup_gather(mtag1);
    
    /* Plaquette in -dir0 +dir1 direction, counted towards integer time */
    /**mtag0 = start_gather_site( LINK_OFFSET(dir0), 
       sizeof(su3_matrix), dir1, EVENANDODD, gen_pt[0] );
       wait_gather(mtag0);  Already gathered above**/
    
    FORALLSITES(i,s){
    IF_ACTIVE(s) {
      mult_su3_an( (su3_matrix *)(gen_pt[0][i]), &LINK(dir0), &temp1[i] );
    }
    }
    mtag1 = start_gather_field( temp1, sizeof(su3_matrix),
                                    OPP_DIR(dir0), EVENANDODD, gen_pt[1] );
    wait_gather(mtag1);
    FORALLSITES(i,s){
    IF_ACTIVE(s) {
      su3_adjoint( (su3_matrix *)(gen_pt[1][i]), &tmat2 );
      add_su3_matrix( &FIELD_STRENGTH(component), (su3_matrix *)(gen_pt[1][i]),
                      &FIELD_STRENGTH(component) );
      sub_su3_matrix( &FIELD_STRENGTH(component), &tmat2,
                      &FIELD_STRENGTH(component) );
    }
    }
    cleanup_gather(mtag0);
    cleanup_gather(mtag1);
    
    // /* Plaquette in -dir0 -dir1 direction, not counted towards integer time */
    // mtag0 = start_gather_site( LINK_OFFSET(dir0), sizeof(su3_matrix),
    //                       OPP_DIR(dir0), EVENANDODD, gen_pt[0] );
    // mtag1 = start_gather_site( LINK_OFFSET(dir1), sizeof(su3_matrix),
    //                       OPP_DIR(dir1), EVENANDODD, gen_pt[1] );
    // wait_gather(mtag0);
    // wait_gather(mtag1);
    // FORALLSITES(i,s){
    // IF_ACTIVE(s) {
    //   mult_su3_nn( (su3_matrix *)(gen_pt[0][i]), &LINK(dir1), &temp1[i] );
    //   mult_su3_nn( (su3_matrix *)(gen_pt[1][i]), &LINK(dir0), &temp2[i] );
    // }
    // }
    // cleanup_gather(mtag0);
    // cleanup_gather(mtag1);
    // mtag0 = start_gather_field( temp1, sizeof(su3_matrix),
    //                                 OPP_DIR(dir1), EVENANDODD, gen_pt[0] );
    // wait_gather(mtag0);
    // mtag1 = start_gather_field( temp2, sizeof(su3_matrix),
    //                                 OPP_DIR(dir0), EVENANDODD, gen_pt[1] );
    // wait_gather(mtag1);
    // FORALLSITES(i,s){
    // IF_ACTIVE(s) {
    //   mult_su3_an( (su3_matrix *)(gen_pt[1][i]), 
    //                (su3_matrix *)(gen_pt[0][i]), &tmat1 );
    //   su3_adjoint( &tmat1, &tmat2 );
    //   add_su3_matrix( &FIELD_STRENGTH(component), &tmat1,
    //                   &FIELD_STRENGTH(component) );
    //   sub_su3_matrix( &FIELD_STRENGTH(component), &tmat2,
    //                   &FIELD_STRENGTH(component) );
    // }
    // }
    // cleanup_gather(mtag0);
    // cleanup_gather(mtag1);
    
    // /* Plaquette in +dir0 -dir1 direction, not counted towards integer time */
    // mtag1 = start_gather_site( LINK_OFFSET(dir1), sizeof(su3_matrix),
    //                       dir0, EVENANDODD, gen_pt[1] );
    // wait_gather(mtag1);
    // FORALLSITES(i,s){
    // IF_ACTIVE(s) {
    //   mult_su3_an( &LINK(dir1), &LINK(dir0), &tmat1);
    //   mult_su3_nn( &tmat1, (su3_matrix *)(gen_pt[1][i]), &temp1[i] );
    // }
    // }
    // cleanup_gather(mtag1);
    // mtag0 = start_gather_field( temp1, sizeof(su3_matrix),
    //                                 OPP_DIR(dir1), EVENANDODD, gen_pt[0] );
    // wait_gather(mtag0);
    // FORALLSITES(i,s){
    // IF_ACTIVE(s) {
    //   mult_su3_na( (su3_matrix *)(gen_pt[0][i]), &LINK(dir0), &tmat1 );
    //   su3_adjoint( &tmat1, &tmat2 );
    //   add_su3_matrix( &FIELD_STRENGTH(component), &tmat1,
    //                   &FIELD_STRENGTH(component) );
    //   sub_su3_matrix( &FIELD_STRENGTH(component), &tmat2,
    //                   &FIELD_STRENGTH(component) );
    // }
    // }
    // cleanup_gather(mtag0);
    #endif

    /* Make traceless */
    FORALLSITES(i,s){
    IF_LOWER_BULK(s) {
      cc = trace_su3(&FIELD_STRENGTH(component));
      CMULREAL(cc,0.33333333333333333,cc);
      for(j=0;j<3;j++)
        CSUB(FIELD_STRENGTH(component).e[j][j],cc,
             FIELD_STRENGTH(component).e[j][j]);
    }
    }
  } /* component */
  
  free(temp1);
  free(temp2);
} /* make_fieldstrength_bulk */
