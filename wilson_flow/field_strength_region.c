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


static void 
make_spatial_fieldstrength_region ( int region_flag, 
                                    su3_matrix **link, 
                                    su3_matrix **fstrength ) {

  register int i, j, icomp;
  register int stride = block_stride;
  int dir[2] = {NODIR,NODIR}, 
      dirb[2] = {NODIR, NODIR}, 
      diro[2] = {NODIR,NODIR};
  register int ig;
  register site *s;
  #define NGATHER 2
  msg_tag *tag[NGATHER];
  #define NTEMP (NGATHER)
  su3_matrix **temp = new_field( NTEMP );

  su3_matrix tmat1,tmat2;
  complex cc;

  /* Calculate the spatial/magnetic components */
  for(icomp=FS_XY;icomp<=FS_YZ;icomp++){
    switch(icomp){
      case FS_XY: dir[0]=XUP; dir[1]=YUP; break;
      case FS_XZ: dir[0]=XUP; dir[1]=ZUP; break;
      case FS_YZ: dir[0]=YUP; dir[1]=ZUP; break;
    }
    setup_blocked_dirs( dir, dirb ,diro );
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
      sub_su3_matrix(  &tmat1, &tmat2, &(fstrength[icomp][i]) );
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
      add_su3_matrix( &(fstrength[icomp][i]), &tmat1, &(fstrength[icomp][i]) );
      sub_su3_matrix( &(fstrength[icomp][i]), &tmat2, &(fstrength[icomp][i]) );
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
      add_su3_matrix( &(fstrength[icomp][i]), &tmat1, &(fstrength[icomp][i]) );
      sub_su3_matrix( &(fstrength[icomp][i]), &tmat2, &(fstrength[icomp][i]) );
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
      add_su3_matrix( &(fstrength[icomp][i]), &tmat1, &(fstrength[icomp][i]) );
      sub_su3_matrix( &(fstrength[icomp][i]), &tmat2, &(fstrength[icomp][i]) );
    }
    cleanup_gather(tag[ig]);


    /* Make traceless */
    FORALLSITES(i,s) 
    IF_BLOCKED(s, block_stride)       
    IF_REGION(s, region_flag) 
    {
      cc = trace_su3( &(fstrength[icomp][i]) );
      CMULREAL(cc,0.33333333333333333,cc);
      for(j=0;j<3;j++)
        CSUB( fstrength[icomp][i].e[j][j],cc, fstrength[icomp][i].e[j][j]);
    }

    #undef LINK0
    #undef LINK1
  } /* component */

  // deallocate temporary storage
  destroy_field( &temp );
  #undef NTEMP
  #undef NGATHER
}

enum edges { Bs=0, Rs=1, Ts=2, Ls=3 };
enum base { P10=0, P11=1, P01=2 };

// returns all four loops U_munu(n) from the same four edges but starting at different points
static void 
// make_clov( su3_matrix edge[4], su3_matrix loop[4], Real coeff ) 
make_clov( su3_matrix *edge, su3_matrix *loop, Real coeff ) 
{
  su3_matrix tmat1, tmat2;
  
  // first set of hooks
  mult_su3_nn( &(edge[Bs]), &(edge[Rs]), &tmat1 );
  mult_su3_nn( &(edge[Ls]), &(edge[Ts]), &tmat2 );
  /* Plaquette in +dir0 +dir1 direction */
  mult_su3_na( &tmat1, &tmat2, &(loop[0]) );
  scalar_mult_su3_matrix( &(loop[0]), coeff, &(loop[0]) );
  /* Plaquette in -dir0 -dir1 direction */
  mult_su3_an( &tmat2, &tmat1, &(loop[2]) );
  scalar_mult_su3_matrix( &(loop[2]), coeff, &(loop[2]) );

  // second set of hooks
  mult_su3_na( &(edge[Rs]), &(edge[Ts]), &tmat1);
  mult_su3_an( &(edge[Ls]), &(edge[Bs]), &tmat2);
  /* Plaquette in -dir0 +dir1 direction */
  mult_su3_nn( &tmat1, &tmat2, &(loop[1]) );
  scalar_mult_su3_matrix( &(loop[1]), coeff, &(loop[1]) );
  /* Plaquette in +dir0 -dir1 direction */
  mult_su3_nn( &tmat2, &tmat1, &(loop[3]) );
  scalar_mult_su3_matrix( &(loop[3]), coeff, &(loop[3]) );
}

/*
make_improved_spatial_fieldstrength_region: original version by
Parikshit Junnarkar (parikshit@theorie.ikp.physik.tu-darmstadt.de)
Simon Stendebach (sstendebach@theorie.ikp.physik.tu-darmstadt.de)
patched up by 
J.H. Weber (dr.rer.nat.weber@gmail.com)
*/
static void
make_improved_spatial_fieldstrength_region ( int region_flag,
                                            su3_matrix **link,
                                            su3_matrix **fstrength) {

  register int i, j, icomp;
  register int stride = block_stride;
  int dir[2] = {NODIR,NODIR}, 
      dirb[2] = {NODIR, NODIR}, 
      diro[2] = {NODIR,NODIR};
  int disp[4];
  register site *s;
  complex cc;

  Real coeff1x1 = 5.0/3.0, coeff1x2 = -1.0/6.0;
  // Real coeff1x1 = 3.0/3.0, coeff1x2 = -0.0/6.0;
  #define NGATHER 8
  msg_tag *tag[NGATHER];
  #define NTEMP (NGATHER)
  su3_matrix **temp = new_field( NTEMP );
  #define NCLOV 3
  su3_matrix **clov = new_field( NCLOV );
  su3_matrix edge[4], loop[4], tmat;

  /* loop over planes */
  for(icomp=FS_XY;icomp<=FS_YZ;icomp++){
    switch(icomp){
    case FS_XY: dir[0]=XUP; dir[1]=YUP; break;
    case FS_XZ: dir[0]=XUP; dir[1]=ZUP; break;
    case FS_YZ: dir[0]=YUP; dir[1]=ZUP; break;
    }
    setup_blocked_dirs( dir, dirb ,diro );
    memset( disp, '\0', 4 * sizeof(int) );
    memset( loop, '\0', 4 * sizeof(su3_matrix) );
    clear_field( &(fstrength[icomp]), 1 ); 
    disp[dir[0]] = -block_stride;
    disp[dir[1]] = -block_stride;
    #define LINK0 link[dir[0]]
    #define LINK1 link[dir[1]]

    // convention opposite of previous function: tag[0] has link in direction 0

    // U_mu(x+nu)
    tag[0] = start_gather_field( LINK0, sizeof(su3_matrix),
                                 dirb[1], EVENANDODD, gen_pt[0]);
  
    // U_nu(x+mu)
    tag[1] = start_gather_field( LINK1, sizeof(su3_matrix),
                                 dirb[0], EVENANDODD, gen_pt[1]);
  
    // U_nu(x+nu)
    tag[2] = start_gather_field( LINK1, sizeof(su3_matrix),
                                 dirb[1], EVENANDODD, gen_pt[2]);
  
    // U_mu(x+mu)
    tag[3] = start_gather_field( LINK0, sizeof(su3_matrix),
                                 dirb[0], EVENANDODD, gen_pt[3]);

    wait_gather(tag[0]); // U_mu @ P01 (x+nu)
    wait_gather(tag[1]); // U_nu @ P10 (x+mu)
    wait_gather(tag[2]); // U_nu @ P01 (x+nu)
    wait_gather(tag[3]); // U_mu @ P10 (x+mu)

    FORALLSITES(i,s)
    IF_BLOCKED(s, block_stride)       
    IF_REGION(s, region_flag)
    {
      su3mat_copy( (su3_matrix *)(gen_pt[0][i]), &(temp[0][i]) ); // U_mu @ P01 (x+nu)
      su3mat_copy( (su3_matrix *)(gen_pt[1][i]), &(temp[1][i]) ); // U_nu @ P10 (x+mu)
      su3mat_copy( (su3_matrix *)(gen_pt[2][i]), &(temp[2][i]) ); // U_nu @ P01 (x+nu)
      su3mat_copy( (su3_matrix *)(gen_pt[3][i]), &(temp[3][i]) ); // U_mu @ P10 (x+mu)
    }

    cleanup_gather(tag[0]); // U_mu @ P01 (x+nu)
    cleanup_gather(tag[1]); // U_nu @ P10 (x+mu)
    cleanup_gather(tag[2]); // U_nu @ P01 (x+nu)
    cleanup_gather(tag[3]); // U_mu @ P10 (x+mu)

    /* now links more than one step removed from site */

    // U_mu(x+2nu)
    tag[0] = start_gather_field( temp[0], sizeof(su3_matrix),
                                 dirb[1], EVENANDODD, gen_pt[0] );
    
    // U_nu(x+mu+nu)
    tag[1] = start_gather_field( temp[1], sizeof(su3_matrix),
                                 dirb[1], EVENANDODD, gen_pt[1] );
    
    // U_mu(x+mu+nu)
    tag[2] = start_gather_field( temp[0], sizeof(su3_matrix),
                                 dirb[0], EVENANDODD, gen_pt[2] );
    
    // U_nu(x+2mu)
    tag[3] = start_gather_field( temp[1], sizeof(su3_matrix),
                                 dirb[0], EVENANDODD, gen_pt[3] );

    wait_gather(tag[0]); // U_mu @ P02 (x+2nu)
    wait_gather(tag[1]); // U_nu @ P11 (x+mu+nu)
    wait_gather(tag[2]); // U_mu @ P11 (x+mu+nu)
    wait_gather(tag[3]); // U_nu @ P20 (x+2mu)

    clear_field( clov, NCLOV );

    FORALLSITES(i,s) 
    IF_BLOCKED(s, block_stride)       
    IF_REGION(s, region_flag)
    {
      /* links for plaquette */
      su3mat_copy( &(LINK0[i]),   &(edge[Bs]) );   // U_mu(x)
      su3mat_copy( &(temp[1][i]), &(edge[Rs]) );   // U_nu @ P10 (x+mu)
      su3mat_copy( &(temp[0][i]), &(edge[Ts]) );   // U_mu @ P01 (x+nu)
      su3mat_copy( &(LINK1[i]),   &(edge[Ls]) );   // U_nu(x)


      /* the plaquettes U_munu(x) starting at different points (counter-clockwise plaqs) */
      make_clov( edge, loop, coeff1x1 );
      add_su3_matrix( &(fstrength[icomp][i]), &(loop[0]), &(fstrength[icomp][i]) ); // plaq starting at x
      add_su3_matrix( &(loop[1]), &(clov[P10][i]), &(clov[P10][i]) ); // plaq starting @ P10 (x+mu)
      add_su3_matrix( &(loop[2]), &(clov[P11][i]), &(clov[P11][i]) ); // plaq starting @ P11 (x+mu+nu)
      add_su3_matrix( &(loop[3]), &(clov[P01][i]), &(clov[P01][i]) ); // plaq starting @ P01 (x+nu)

      /* sides of 1x2 rectangle */
      // U_mu(x) is still in edge[Bs]
      mult_su3_nn( &(temp[1][i]), (su3_matrix *)(gen_pt[1][i]), &(edge[Rs]) ); // U_nu @ P10 (x+mu) * U_nu @ P11 (x+mu+nu)
      su3mat_copy( (su3_matrix *)(gen_pt[0][i]), &(edge[Ts]) ); // U_mu @ P02 (x+2nu)
      mult_su3_nn( &(LINK1[i]), &(temp[2][i]), &(edge[Ls]) ); // U_nu @ (x) * U_nu @ P01 (x+nu)

      /* construct 1x2 rectangles U_munu(x) but "rotated" */
      make_clov( edge, loop, coeff1x2 );
      add_su3_matrix(  &(fstrength[icomp][i]), &(loop[0]), &(fstrength[icomp][i])); // rect starting at x
      add_su3_matrix(  &(loop[1]), &(clov[P10][i]), &(clov[P10][i]) ) ; // rect starting @ P10 (x+mu), gather it from left (-mu) direction later
      su3mat_copy( &(loop[2]), &(temp[4][i]) ); // rect starting @ P12 (x+mu+2nu), store loop[2] in array to gather it from down (-2nu) and left (-mu) direction later
      su3mat_copy( &(loop[3]), &(temp[5][i]) ); // rect starting @ P02 (x+2nu), store loop[3] in array to gather it from down (-2nu) direction later

      /* sides of 2x1 rectangle */
      mult_su3_nn( &(LINK0[i]), &(temp[3][i]), &(edge[Bs]) ); // U_mu @ (x) * U_mu  @ P10 (x+mu)
      su3mat_copy( (su3_matrix *)(gen_pt[3][i]), &(edge[Rs]) ); // U_nu @ P20 (x+2mu)
      mult_su3_nn( &(temp[0][i]), (su3_matrix *)(gen_pt[2][i]), &(edge[Ts]) ); // U_mu @ P01 (x+nu) * U_mu @ P11 (x+mu+nu)
      su3mat_copy( &(LINK1[i]), &(edge[Ls]) ); // U_nu @ (x)

      /* construct 2x1 rectangles U_munu(x) but "rotated" */
      make_clov( edge, loop, coeff1x2 );
      add_su3_matrix(  &(fstrength[icomp][i]), &(loop[0]), &(fstrength[icomp][i]) ); // rect starting at x
      su3mat_copy( &(loop[1]), &(temp[6][i]) );   // rect starting @ P20 (x+2mu), store loop[2] in array to gather it from left (-2mu) direction later
      su3mat_copy( &(loop[2]), &(temp[7][i]) );   // rect starting @ P21 (x+2mu+nu), store loop[3] in array to gather it from down (-nu) and left (-2mu) direction later
      add_su3_matrix(  &(loop[3]), &(clov[P01][i]), &(clov[P01][i]) ); // rect starting @ P01 (x+nu), gather it from down (-nu) direction later 
    } // end FORALLSITES

    cleanup_gather(tag[0]);
    cleanup_gather(tag[1]);
    cleanup_gather(tag[2]);
    cleanup_gather(tag[3]);

    // temp[4]: 1x2, starts @ P12 (x+mu+2nu) upper right corner
    // temp[5]: 1x2, starts @ P02 (x+2nu) upper left corner
    // temp[6]: 2x1, starts @ P20 (x+2mu) lower right corner
    // temp[7]: 2x1, starts @ P21 (x+2mu+nu) upper right corner

    /* at this point, upper left (clov[P10]) and lower right (clov[P01]) parts miss one rectangle,
       lower left (clov[P11]) misses two, so add them now */

    // gather temp[4] from down direction
    tag[4] = start_gather_field( temp[4], sizeof(su3_matrix), 
                                 diro[1], EVENANDODD, gen_pt[4] ); // contributes to clov @ P11 (x+mu+nu)
    // gather temp[5] from down direction
    tag[5] = start_gather_field( temp[5], sizeof(su3_matrix), 
                                 diro[1], EVENANDODD, gen_pt[5] ); // contributes to clov @ P01 (x+nu)
    // gather temp[6] from left direction
    tag[6] = start_gather_field( temp[6], sizeof(su3_matrix), 
                                 diro[0], EVENANDODD, gen_pt[6] ); // contributes to clov @ P10 (x+mu)
    //gather temp[7] from left direction
    tag[7] = start_gather_field( temp[7], sizeof(su3_matrix), 
                                 diro[0], EVENANDODD, gen_pt[7] ); // contributes to clov @ P11 (x+mu+nu)

    wait_gather(tag[4]);
    wait_gather(tag[5]);
    wait_gather(tag[6]);
    wait_gather(tag[7]);

    FORALLSITES(i,s)
    IF_BLOCKED(s, block_stride)       
    IF_REGION(s, region_flag)
    {
      add_su3_matrix( (su3_matrix *)(gen_pt[4][i]), &(clov[P11][i]), &(clov[P11][i]) );
      add_su3_matrix( (su3_matrix *)(gen_pt[5][i]), &(clov[P01][i]), &(clov[P01][i]) );
      add_su3_matrix( (su3_matrix *)(gen_pt[6][i]), &(clov[P10][i]), &(clov[P10][i]) );
      add_su3_matrix( (su3_matrix *)(gen_pt[7][i]), &(clov[P11][i]), &(clov[P11][i]) );
    }

    cleanup_gather(tag[4]);
    cleanup_gather(tag[5]);
    cleanup_gather(tag[6]);
    cleanup_gather(tag[7]);


    /* general gather lower left clover from diagonally across the square */
    tag[P11] = start_general_gather_field( clov[P11], sizeof(su3_matrix), 
                                           disp, EVENANDODD, gen_pt[P11] );

    /* gather upper left clover from proper position */
    tag[P10] = start_gather_field( clov[P10], sizeof(su3_matrix), 
                                   diro[0], EVENANDODD, gen_pt[P10] );
    
    /* gather lower right clover from proper position */
    tag[P01] = start_gather_field( clov[P01], sizeof(su3_matrix), 
                                   diro[1], EVENANDODD, gen_pt[P01] );


    wait_gather(tag[P10]);
    wait_gather(tag[P01]);
    wait_general_gather(tag[P11]);


    FORALLSITES(i,s)
    IF_BLOCKED(s, block_stride)       
    IF_REGION(s, region_flag) 
    {
      add_su3_matrix( &(fstrength[icomp][i]), (su3_matrix *)(gen_pt[P10][i]), &(fstrength[icomp][i]) );
      add_su3_matrix( &(fstrength[icomp][i]), (su3_matrix *)(gen_pt[P11][i]), &(fstrength[icomp][i]) );
      add_su3_matrix( &(fstrength[icomp][i]), (su3_matrix *)(gen_pt[P01][i]), &(fstrength[icomp][i]) );
    }
    
    cleanup_general_gather(tag[P11]);
    cleanup_gather(tag[P10]);
    cleanup_gather(tag[P01]);

    // cleanup_gather(tag[P11]);

    /* Make traceless */
    FORALLSITES(i,s)
    IF_BLOCKED(s, block_stride)       
    IF_REGION(s, region_flag) 
    {
      su3_adjoint( &(fstrength[icomp][i]), &tmat );
      sub_su3_matrix( &(fstrength[icomp][i]), &tmat, &(fstrength[icomp][i]) );
      
      cc = trace_su3(&(fstrength[icomp][i]));
      CMULREAL(cc,0.33333333333333333,cc);
      for(j=0;j<3;j++)
      {
        CSUB( fstrength[icomp][i].e[j][j],cc, fstrength[icomp][i].e[j][j]);
      } 
    }

    #undef LINK0
    #undef LINK1
  } // end loop over components

  // deallocate temporary storage
  destroy_field( &temp );
  destroy_field( &clov ); 
  #undef NTEMP
  #undef NCLOV
  #undef NGATHER
}

static void 
make_temporal_fieldstrength_full ( su3_matrix **link_s, 
                                   su3_matrix **link_t, 
                                   su3_matrix **fstrength ) {

  register int i, j, icomp;
  register int stride = block_stride;
  int dir[2] = {NODIR,NODIR}, 
      dirb[2] = {NODIR, NODIR}, 
      diro[2] = {NODIR,NODIR};
  register int ig;
  register site *s;
  #define NGATHER 2
  msg_tag *tag[NGATHER];
  #define NTEMP (NGATHER)
  su3_matrix **temp = new_field( NTEMP );

  su3_matrix tmat1,tmat2;
  complex cc;

  /* Calculate the temporal/electric components */  
  for(icomp=FS_XT;icomp<=FS_ZT;icomp++){
    switch(icomp){
      case FS_XT: dir[0]=XUP; dir[1]=TUP; break;
      case FS_YT: dir[0]=YUP; dir[1]=TUP; break;
      case FS_ZT: dir[0]=ZUP; dir[1]=TUP; break;
    }
    setup_blocked_dirs( dir, dirb ,diro );
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
      sub_su3_matrix(  &tmat1, &tmat2, &(fstrength[icomp][i]) );
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
      add_su3_matrix( &(fstrength[icomp][i]), &tmat1,
                      &(fstrength[icomp][i]) );
      sub_su3_matrix( &(fstrength[icomp][i]), &tmat2,
                      &(fstrength[icomp][i]) );
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
      add_su3_matrix( &(fstrength[icomp][i]), &tmat1,
                      &(fstrength[icomp][i]) );
      sub_su3_matrix( &(fstrength[icomp][i]), &tmat2,
                      &(fstrength[icomp][i]) );
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
      add_su3_matrix( &(fstrength[icomp][i]), &tmat1,
                      &(fstrength[icomp][i]) );
      sub_su3_matrix( &(fstrength[icomp][i]), &tmat2,
                      &(fstrength[icomp][i]) );
    }
    cleanup_gather(tag[ig]);

    /* Make traceless */
    FORALLSITES(i,s) 
    IF_BLOCKED(s, block_stride)       
    {
      cc = trace_su3( &(fstrength[icomp][i]) );
      CMULREAL(cc,0.33333333333333333,cc);
      for(j=0;j<3;j++)
        CSUB( fstrength[icomp][i].e[j][j],cc,
              fstrength[icomp][i].e[j][j]);
    }

    #undef LINK0
    #undef LINK1
  } /* component */

  // deallocate temporary storage
  destroy_field( &temp );
  #undef NTEMP
  #undef NGATHER
}

static void
make_improved_temporal_fieldstrength_full ( su3_matrix **link_s, 
                                            su3_matrix **link_t,
                                            su3_matrix **fstrength) {

  register int i, j, icomp;
  register int stride = block_stride;
  int dir[2] = {NODIR,NODIR}, 
      dirb[2] = {NODIR, NODIR}, 
      diro[2] = {NODIR,NODIR};
  int disp[4];
  register int ig;
  register site *s;
  complex cc;

  Real coeff1x1 = 5.0/3.0, coeff1x2 = -1.0/6.0;
  // Real coeff1x1 = 0.0/3.0, coeff1x2 = -1.0/6.0;
  #define NGATHER 8
  msg_tag *tag[NGATHER],*mtagx;
  #define NTEMP (NGATHER)
  su3_matrix **temp = new_field( NTEMP );
  #define NCLOV 4
  su3_matrix **clov = new_field( NCLOV );
  su3_matrix edge[4], loop[4], tmat;

  /* Calculate the temporal/electric components */  
  for(icomp=FS_XT;icomp<=FS_ZT;icomp++){
    switch(icomp){
      case FS_XT: dir[0]=XUP; dir[1]=TUP; break;
      case FS_YT: dir[0]=YUP; dir[1]=TUP; break;
      case FS_ZT: dir[0]=ZUP; dir[1]=TUP; break;
    }
    setup_blocked_dirs( dir, dirb ,diro );
    memset( disp, '\0', 4 * sizeof(int) );
    disp[dir[0]] = -block_stride;
    disp[dir[1]] = -1;
    #define LINK0 link_s[dir[0]]
    #define LINK1 link_t[0]


    // convention opposite of previous function: tag[0] has link in direction 0

    // U_mu(x+nu)
    tag[0] = start_gather_field( LINK0, sizeof(su3_matrix),
                                 dirb[1], EVENANDODD, gen_pt[0]);
  
    // U_nu(x+mu)
    tag[1] = start_gather_field( LINK1, sizeof(su3_matrix),
                                 dirb[0], EVENANDODD, gen_pt[1]);
  
    // U_nu(x+nu)
    tag[2] = start_gather_field( LINK1, sizeof(su3_matrix),
                                 dirb[1], EVENANDODD, gen_pt[2]);
  
    // U_mu(x+mu)
    tag[3] = start_gather_field( LINK0, sizeof(su3_matrix),
                                 dirb[0], EVENANDODD, gen_pt[3]);

    wait_gather(tag[0]); // U_mu @ P01 (x+nu)
    wait_gather(tag[1]); // U_nu @ P10 (x+mu)
    wait_gather(tag[2]); // U_nu @ P01 (x+nu)
    wait_gather(tag[3]); // U_mu @ P10 (x+mu)

    FORALLSITES(i,s)
    IF_BLOCKED(s, block_stride)       
    // IF_REGION(s, region_flag)
    {
      su3mat_copy( (su3_matrix *)(gen_pt[0][i]), &(temp[0][i]) ); // U_mu @ P01 (x+nu)
      su3mat_copy( (su3_matrix *)(gen_pt[1][i]), &(temp[1][i]) ); // U_nu @ P10 (x+mu)
      su3mat_copy( (su3_matrix *)(gen_pt[2][i]), &(temp[2][i]) ); // U_nu @ P01 (x+nu)
      su3mat_copy( (su3_matrix *)(gen_pt[3][i]), &(temp[3][i]) ); // U_mu @ P10 (x+mu)
    }

    cleanup_gather(tag[0]); // U_mu @ P01 (x+nu)
    cleanup_gather(tag[1]); // U_nu @ P10 (x+mu)
    cleanup_gather(tag[2]); // U_nu @ P01 (x+nu)
    cleanup_gather(tag[3]); // U_mu @ P10 (x+mu)

    /* now links more than one step removed from site */

    // U_mu(x+2nu)
    tag[0] = start_gather_field( temp[0], sizeof(su3_matrix),
                                 dirb[1], EVENANDODD, gen_pt[0] );
    
    // U_nu(x+mu+nu)
    tag[1] = start_gather_field( temp[1], sizeof(su3_matrix),
                                 dirb[1], EVENANDODD, gen_pt[1] );
    
    // U_mu(x+mu+nu)
    tag[2] = start_gather_field( temp[0], sizeof(su3_matrix),
                                 dirb[0], EVENANDODD, gen_pt[2] );
    
    // U_nu(x+2mu)
    tag[3] = start_gather_field( temp[1], sizeof(su3_matrix),
                                 dirb[0], EVENANDODD, gen_pt[3] );

    wait_gather(tag[0]); // U_mu @ P02 (x+2nu)
    wait_gather(tag[1]); // U_nu @ P11 (x+mu+nu)
    wait_gather(tag[2]); // U_mu @ P11 (x+mu+nu)
    wait_gather(tag[3]); // U_nu @ P20 (x+2mu)

    clear_field( clov, NCLOV );

    FORALLSITES(i,s) 
    IF_BLOCKED(s, block_stride)       
    // IF_REGION(s, region_flag)
    {
      /* links for plaquette */
      su3mat_copy( &(LINK0[i]),   &(edge[Bs]) );   // U_mu(x)
      su3mat_copy( &(temp[1][i]), &(edge[Rs]) );   // U_nu @ P10 (x+mu)
      su3mat_copy( &(temp[0][i]), &(edge[Ts]) );   // U_mu @ P01 (x+nu)
      su3mat_copy( &(LINK1[i]),   &(edge[Ls]) );   // U_nu(x)


      /* the plaquettes U_munu(x) starting at different points (counter-clockwise plaqs) */
      make_clov( edge, loop, coeff1x1 );
      add_su3_matrix( &(fstrength[icomp][i]), &(loop[0]), &(fstrength[icomp][i]) ); // plaq starting at x
      add_su3_matrix( &(loop[1]), &(clov[P10][i]), &(clov[P10][i]) ); // plaq starting @ P10 (x+mu)
      add_su3_matrix( &(loop[2]), &(clov[P11][i]), &(clov[P11][i]) ); // plaq starting @ P11 (x+mu+nu)
      add_su3_matrix( &(loop[3]), &(clov[P01][i]), &(clov[P01][i]) ); // plaq starting @ P01 (x+nu)

      /* sides of 1x2 rectangle */
      // U_mu(x) is still in edge[Bs]
      mult_su3_nn( &(temp[1][i]), (su3_matrix *)(gen_pt[1][i]), &(edge[Rs]) ); // U_nu @ P10 (x+mu) * U_nu @ P11 (x+mu+nu)
      su3mat_copy( (su3_matrix *)(gen_pt[0][i]), &(edge[Ts]) ); // U_mu @ P02 (x+2nu)
      mult_su3_nn( &(LINK1[i]), &(temp[2][i]), &(edge[Ls]) ); // U_nu @ (x) * U_nu @ P01 (x+nu)

      /* construct 1x2 rectangles U_munu(x) but "rotated" */
      make_clov( edge, loop, coeff1x2 );
      add_su3_matrix(  &(fstrength[icomp][i]), &(loop[0]), &(fstrength[icomp][i])); // rect starting at x
      add_su3_matrix(  &(loop[1]), &(clov[P10][i]), &(clov[P10][i]) ) ; // rect starting @ P10 (x+mu), gather it from left (-mu) direction later
      su3mat_copy( &(loop[2]), &(temp[4][i]) ); // rect starting @ P12 (x+mu+2nu), store loop[2] in array to gather it from down (-2nu) and left (-mu) direction later
      su3mat_copy( &(loop[3]), &(temp[5][i]) ); // rect starting @ P02 (x+2nu), store loop[3] in array to gather it from down (-2nu) direction later

      /* sides of 2x1 rectangle */
      mult_su3_nn( &(LINK0[i]), &(temp[3][i]), &(edge[Bs]) ); // U_mu @ (x) * U_mu  @ P10 (x+mu)
      su3mat_copy( (su3_matrix *)(gen_pt[3][i]), &(edge[Rs]) ); // U_nu @ P20 (x+2mu)
      mult_su3_nn( &(temp[0][i]), (su3_matrix *)(gen_pt[2][i]), &(edge[Ts]) ); // U_mu @ P01 (x+nu) * U_mu @ P11 (x+mu+nu)
      su3mat_copy( &(LINK1[i]), &(edge[Ls]) ); // U_nu @ (x)

      /* construct 2x1 rectangles U_munu(x) but "rotated" */
      make_clov( edge, loop, coeff1x2 );
      add_su3_matrix(  &(fstrength[icomp][i]), &(loop[0]), &(fstrength[icomp][i]) ); // rect starting at x
      su3mat_copy( &(loop[1]), &(temp[6][i]) );   // rect starting @ P20 (x+2mu), store loop[2] in array to gather it from left (-2mu) direction later
      su3mat_copy( &(loop[2]), &(temp[7][i]) );   // rect starting @ P21 (x+2mu+nu), store loop[3] in array to gather it from down (-nu) and left (-2mu) direction later
      add_su3_matrix(  &(loop[3]), &(clov[P01][i]), &(clov[P01][i]) ); // rect starting @ P01 (x+nu), gather it from down (-nu) direction later 
    } // end FORALLSITES

    cleanup_gather(tag[0]);
    cleanup_gather(tag[1]);
    cleanup_gather(tag[2]);
    cleanup_gather(tag[3]);

    // temp[4]: 1x2, starts @ P12 (x+mu+2nu) upper right corner
    // temp[5]: 1x2, starts @ P02 (x+2nu) upper left corner
    // temp[6]: 2x1, starts @ P20 (x+2mu) lower right corner
    // temp[7]: 2x1, starts @ P21 (x+2mu+nu) upper right corner

    /* at this point, upper left (clov[P10]) and lower right (clov[P01]) parts miss one rectangle,
       lower left (clov[P11]) misses two, so add them now */

    // gather temp[4] from down direction
    tag[4] = start_gather_field( temp[4], sizeof(su3_matrix), 
                                 diro[1], EVENANDODD, gen_pt[4] ); // contributes to clov @ P11 (x+mu+nu)
    // gather temp[5] from down direction
    tag[5] = start_gather_field( temp[5], sizeof(su3_matrix), 
                                 diro[1], EVENANDODD, gen_pt[5] ); // contributes to clov @ P01 (x+nu)
    // gather temp[6] from left direction
    tag[6] = start_gather_field( temp[6], sizeof(su3_matrix), 
                                 diro[0], EVENANDODD, gen_pt[6] ); // contributes to clov @ P10 (x+mu)
    //gather temp[7] from left direction
    tag[7] = start_gather_field( temp[7], sizeof(su3_matrix), 
                                 diro[0], EVENANDODD, gen_pt[7] ); // contributes to clov @ P11 (x+mu+nu)

    wait_gather(tag[4]);
    wait_gather(tag[5]);
    wait_gather(tag[6]);
    wait_gather(tag[7]);

    FORALLSITES(i,s)
    IF_BLOCKED(s, block_stride)       
    // IF_REGION(s, region_flag)
    {
      add_su3_matrix( (su3_matrix *)(gen_pt[4][i]), &(clov[P11][i]), &(clov[P11][i]) );
      add_su3_matrix( (su3_matrix *)(gen_pt[5][i]), &(clov[P01][i]), &(clov[P01][i]) );
      add_su3_matrix( (su3_matrix *)(gen_pt[6][i]), &(clov[P10][i]), &(clov[P10][i]) );
      add_su3_matrix( (su3_matrix *)(gen_pt[7][i]), &(clov[P11][i]), &(clov[P11][i]) );
    }

    cleanup_gather(tag[4]);
    cleanup_gather(tag[5]);
    cleanup_gather(tag[6]);
    cleanup_gather(tag[7]);


    /* general gather lower left clover from diagonally across the square */
    tag[P11] = start_general_gather_field( clov[P11], sizeof(su3_matrix), 
                                           disp, EVENANDODD, gen_pt[P11] );

    /* gather upper left clover from proper position */
    tag[P10] = start_gather_field( clov[P10], sizeof(su3_matrix), 
                                   diro[0], EVENANDODD, gen_pt[P10] );
    
    /* gather lower right clover from proper position */
    tag[P01] = start_gather_field( clov[P01], sizeof(su3_matrix), 
                                   diro[1], EVENANDODD, gen_pt[P01] );


    wait_gather(tag[P10]);
    wait_gather(tag[P01]);
    wait_general_gather(tag[P11]);


    FORALLSITES(i,s)
    IF_BLOCKED(s, block_stride)       
    // IF_REGION(s, region_flag) 
    {
      add_su3_matrix( &(fstrength[icomp][i]), (su3_matrix *)(gen_pt[P10][i]), &(fstrength[icomp][i]) );
      add_su3_matrix( &(fstrength[icomp][i]), (su3_matrix *)(gen_pt[P11][i]), &(fstrength[icomp][i]) );
      add_su3_matrix( &(fstrength[icomp][i]), (su3_matrix *)(gen_pt[P01][i]), &(fstrength[icomp][i]) );
    }
    
    cleanup_general_gather(tag[P11]);
    cleanup_gather(tag[P10]);
    cleanup_gather(tag[P01]);

    // cleanup_gather(tag[P11]);

    /* Make traceless */
    FORALLSITES(i,s)
    IF_BLOCKED(s, block_stride)       
    // IF_REGION(s, region_flag) 
    {
      su3_adjoint( &(fstrength[icomp][i]), &tmat );
      sub_su3_matrix( &(fstrength[icomp][i]), &tmat, &(fstrength[icomp][i]) );
      
      cc = trace_su3(&(fstrength[icomp][i]));
      CMULREAL(cc,0.33333333333333333,cc);
      for(j=0;j<3;j++)
      {
        CSUB( fstrength[icomp][i].e[j][j],cc, fstrength[icomp][i].e[j][j]);
      } 
    }

  #undef LINK0
  #undef LINK1
  } // end loop over components

  // deallocate temporary storage
  destroy_field( &temp );
  destroy_field( &clov ); 
  #undef NTEMP
  #undef NCLOV

  #undef NGATHER
}


static void 
make_temporal_fieldstrength_lwr_bulk ( su3_matrix **link_s, 
                                       su3_matrix **link_t, 
                                       su3_matrix **fstrength ) {

  register int i, j, icomp;
  register int stride = block_stride;
  int dir[2] = {NODIR,NODIR}, 
      dirb[2] = {NODIR, NODIR}, 
      diro[2] = {NODIR,NODIR};
  register int ig;
  register site *s;
  #define NGATHER 2
  msg_tag *tag[NGATHER];
  #define NTEMP (NGATHER)
  su3_matrix **temp = new_field( NTEMP );

  su3_matrix tmat1,tmat2;
  complex cc;

  /* Calculate the temporal/electric components */  
  for(icomp=FS_XT;icomp<=FS_ZT;icomp++){
    switch(icomp){
      case FS_XT: dir[0]=XUP; dir[1]=TUP; break;
      case FS_YT: dir[0]=YUP; dir[1]=TUP; break;
      case FS_ZT: dir[0]=ZUP; dir[1]=TUP; break;
    }
    setup_blocked_dirs( dir, dirb ,diro );
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
      sub_su3_matrix(  &tmat1, &tmat2, &(fstrength[icomp][i]) );
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
      add_su3_matrix( &(fstrength[icomp][i]), &tmat1,
                      &(fstrength[icomp][i]) );
      sub_su3_matrix( &(fstrength[icomp][i]), &tmat2,
                      &(fstrength[icomp][i]) );
    }
    cleanup_gather(tag[0]);
    cleanup_gather(tag[1]);
 
    /* Make traceless */
    FORALLSITES(i,s) 
    IF_BLOCKED(s, block_stride)       
    IF_LOWER_BULK(s) 
    {
      cc = trace_su3( &(fstrength[icomp][i]) );
      CMULREAL(cc,0.33333333333333333,cc);
      for(j=0;j<3;j++)
        CSUB( fstrength[icomp][i].e[j][j],cc,
              fstrength[icomp][i].e[j][j]);
    }

    #undef LINK0
    #undef LINK1
  } /* component */

  // deallocate temporary storage
  destroy_field( &temp );
  #undef NTEMP
  #undef NGATHER
}

static void 
make_dropped_temporal_fieldstrength_lwr_bulk ( su3_matrix **link, 
                                               su3_matrix **fstrength ) {

  register int i, j, icomp;
  register int stride = block_stride;
  int dir[2] = {NODIR,NODIR}, 
      dirb[2] = {NODIR, NODIR}, 
      diro[2] = {NODIR,NODIR};
  register int ig;
  register site *s;
  #define NGATHER 2
  msg_tag *tag[NGATHER];
  #define NTEMP (NGATHER)
  su3_matrix **temp = new_field( NTEMP );

  su3_matrix tmat1,tmat2;
  complex cc;

  /* Calculate the temporal/electric components */  
  for(icomp=FS_XT;icomp<=FS_ZT;icomp++){
    switch(icomp){
      case FS_XT: dir[0]=XUP; dir[1]=TUP; break;
      case FS_YT: dir[0]=YUP; dir[1]=TUP; break;
      case FS_ZT: dir[0]=ZUP; dir[1]=TUP; break;
    }
    setup_blocked_dirs( dir, dirb ,diro );
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
      sub_su3_matrix(  &tmat1, &tmat2, &(fstrength[icomp][i]) );
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
      add_su3_matrix( &(fstrength[icomp][i]), (su3_matrix *)(gen_pt[0][i]),
                      &(fstrength[icomp][i]) );
      sub_su3_matrix( &(fstrength[icomp][i]), &tmat2,
                      &(fstrength[icomp][i]) );
    }
    cleanup_gather(tag[0]);
    cleanup_gather(tag[1]);

    /* Make traceless */
    FORALLSITES(i,s) 
    IF_BLOCKED(s, block_stride)       
    IF_LOWER_BULK(s) 
    {
      cc = trace_su3( &(fstrength[icomp][i]) );
      CMULREAL(cc,0.33333333333333333,cc);
      for(j=0;j<3;j++)
        CSUB( fstrength[icomp][i].e[j][j],cc,
              fstrength[icomp][i].e[j][j]);
    }

    #undef LINK0
  } /* component */

  // deallocate temporary storage
  destroy_field( &temp );
  #undef NTEMP
  #undef NGATHER
}

void 
make_fieldstrength_region ( int region_flag, 
                            su3_matrix **fstrength ) {

  register int i, dir;
  register site *s;
  #define NLINK 4
  su3_matrix **link = new_links_from_site( region_flag, NLINK );

  make_spatial_fieldstrength_region( region_flag, link, fstrength );
  make_improved_spatial_fieldstrength_region ( region_flag, link, fstrength + 6 );

  if ( region_flag == FULLVOL ) {
    make_temporal_fieldstrength_full ( link, &(link[TUP]), fstrength );
    make_improved_temporal_fieldstrength_full ( link, &(link[TUP]), fstrength + 6 );
  }
  if ( region_flag == LOWER_BULK ) {
    #ifdef DROP_TIME_LINKS
      make_dropped_temporal_fieldstrength_lwr_bulk ( link, fstrength );
    #else
      make_temporal_fieldstrength_lwr_bulk ( link, &(link[TUP]), fstrength );
    #endif
  }

  // deallocate temporary spatial links
  destroy_field( &link );
  #undef NLINK
}

void 
make_fieldstrength_bulk ( su3_matrix **fstrength ) {

  register int i, dir;
  register site *s;
  #ifndef DROP_TIME_LINKS
  #define NLINK 4
  #else
  #define NLINK 3
  #endif
  su3_matrix **link = new_links_from_site( ACTIVE, NLINK );

  make_spatial_fieldstrength_region( ACTIVE, link, fstrength );
  make_improved_spatial_fieldstrength_region ( ACTIVE, link, fstrength + 6 );

  // #undef DROP_TIME_LINKS
  #ifdef DROP_TIME_LINKS
    make_dropped_temporal_fieldstrength_lwr_bulk ( link, fstrength );
  #else
    make_temporal_fieldstrength_lwr_bulk ( link, &(link[TUP]), fstrength );
  #endif

  // deallocate temporary spatial links
  destroy_field( &link );
  #undef NLINK
}

void 
make_fieldstrength_full ( su3_matrix **fstrength ) {

  register int i, dir;
  register site *s;
  #define NLINK 4
  su3_matrix **link = new_links_from_site( FULLVOL, NLINK );

  make_spatial_fieldstrength_region( FULLVOL, link, fstrength );
  make_improved_spatial_fieldstrength_region ( FULLVOL, link, fstrength + 6 );

  make_temporal_fieldstrength_full ( link, &(link[TUP]), fstrength );
  make_improved_temporal_fieldstrength_full ( link, &(link[TUP]), fstrength + 6 );

  // deallocate temporary spatial links
  destroy_field( &link );
  #undef NLINK
}


#ifdef SPHALERON
static void 
make_dropped_temporal_fieldstrength_half ( su3_matrix **link, 
                                           su3_matrix **link_half, 
                                           su3_matrix **fstrength ) {

  register int i, j, icomp, jcomp;
  register int stride = block_stride;
  int dir[2] = {NODIR,NODIR}, 
      dirb[2] = {NODIR, NODIR}, 
      diro[2] = {NODIR,NODIR};
  register int ig;
  register site *s;
  #define NGATHER 2
  msg_tag *tag[NGATHER];
  #define NTEMP (NGATHER)
  su3_matrix **temp = new_field( NTEMP );
  
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
    setup_blocked_dirs( dir, dirb ,diro );
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
      sub_su3_matrix(  &tmat1, &tmat2, &(fstrength[icomp][i]) );
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
      sub_su3_matrix(  &tmat1, &tmat2, &(fstrength[jcomp][i]) );
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
      add_su3_matrix( &(fstrength[icomp][i]), (su3_matrix *)(gen_pt[0][i]),
                      &(fstrength[icomp][i]) );
      sub_su3_matrix( &(fstrength[icomp][i]), &tmat2,
                      &(fstrength[icomp][i]) );
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
      add_su3_matrix( &(fstrength[jcomp][i]), (su3_matrix *)(gen_pt[0][i]),
                      &(fstrength[jcomp][i]) );
      sub_su3_matrix( &(fstrength[jcomp][i]), &tmat2,
                      &(fstrength[jcomp][i]) );
    }
    cleanup_gather(tag[0]);
    cleanup_gather(tag[1]);

    /* Make traceless */
    FORALLSITES(i,s) 
    IF_BLOCKED(s, block_stride)       
    IF_LOWER_BULK(s) 
    {
      cc = trace_su3( &(fstrength[icomp][i]) );
      CMULREAL(cc,0.33333333333333333,cc);
      for(j=0;j<3;j++)
        CSUB( fstrength[icomp][i].e[j][j],cc, fstrength[icomp][i].e[j][j]);

      cc = trace_su3( &(fstrength[jcomp][i]) );
      CMULREAL(cc,0.33333333333333333,cc);
      for(j=0;j<3;j++)
        CSUB( fstrength[jcomp][i].e[j][j],cc, fstrength[jcomp][i].e[j][j]);
    }

    #undef LINKF
    #undef LINKH
  } /* component */

  // deallocate temporary storage
  destroy_field( &temp );
  #undef NTEMP
  #undef NGATHER
}

static void 
make_dropped_temporal_fieldstrength_bdry ( su3_matrix **link, 
                                           su3_matrix **link_lf, 
                                           su3_matrix **fstrength ) {

  register int i, j, icomp;
  register int stride = block_stride;
  #define NGATHER 1
  int dir[2] = {NODIR,NODIR}, 
      dirb[2] = {NODIR,NODIR}, 
      diro[2] = {NODIR,NODIR};
  register int ig;
  register site *s;
  msg_tag *tag[NGATHER];
  #define NTEMP (NGATHER)
  su3_matrix **temp = new_field( NTEMP );

  su3_matrix tmat1,tmat2;
  complex cc;

  /* Calculate the temporal/electric components */  
  for(icomp=FS_XT;icomp<=FS_ZT;icomp++){
    switch(icomp){
      case FS_XT: dir[0]=XUP; dir[1]=TUP; break;
      case FS_YT: dir[0]=YUP; dir[1]=TUP; break;
      case FS_ZT: dir[0]=ZUP; dir[1]=TUP; break;
    }
    setup_blocked_dirs( dir, dirb ,diro );
    #define LINK0 link[dir[0]]
    #define LINK1 link_lf[dir[0]]

    /* Consistency for improved readability: gathers in 
       +- dir[0] direction have tag[0] and gen_pt[0], use temp[0] 
       +- dir[1] direction have tag[1] and gen_pt[1], use temp[1] */ 

    /* Plaquette in +dir[0] +dir[1] direction */
    FORALLSITES(i,s) 
    IF_BLOCKED(s, block_stride)       
    IF_BOUNDARY(s) 
    {
      mult_su3_na( &(LINK0[i]), &(LINK1[i]), &tmat1 );
      su3_adjoint( &tmat1, &tmat2 );
      sub_su3_matrix(  &tmat1, &tmat2, &(fstrength[icomp][i]) );
    }

    /* Plaquette in -dir[0] +dir[1] direction */
    
    FORALLSITES(i,s) 
    IF_BLOCKED(s, block_stride)       
    IF_BOUNDARY(s) 
      mult_su3_an(  &(LINK1[i]), &(LINK0[i]), &temp[0][i] );

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
      add_su3_matrix( &(fstrength[icomp][i]), (su3_matrix *)(gen_pt[0][i]),
                      &(fstrength[icomp][i]) );
      sub_su3_matrix( &(fstrength[icomp][i]), &tmat2,
                      &(fstrength[icomp][i]) );
    }
    cleanup_gather(tag[0]);

    /* Consistency for improved readability: gathers in 
       +- dir[0] direction have tag[0] and gen_pt[0], use temp[0] 
       +- dir[1] direction have tag[1] and gen_pt[1], use temp[1] */ 

    /* Make traceless */
    FORALLSITES(i,s) 
    IF_BLOCKED(s, block_stride)       
    IF_BOUNDARY(s) 
    {
      cc = trace_su3( &(fstrength[icomp][i]) );
      CMULREAL(cc,0.33333333333333333,cc);
      for(j=0;j<3;j++)
        CSUB( fstrength[icomp][i].e[j][j],cc,
              fstrength[icomp][i].e[j][j]);

    }

    #undef LINK0
    #undef LINK1
  } /* component */

  // deallocate temporary storage
  destroy_field( &temp );
  #undef NTEMP
  #undef NGATHER
}

static void 
make_dropped_temporal_fieldstrength_lwr_bdry ( su3_matrix **link, 
                                           su3_matrix **link_lf, 
                                           su3_matrix **fstrength ) {

  register int i, j, icomp;
  register int stride = block_stride;
  #define NGATHER 1
  int dir[2] = {NODIR,NODIR}, 
      dirb[2] = {NODIR,NODIR}, 
      diro[2] = {NODIR,NODIR};
  register int ig;
  register site *s;
  msg_tag *tag[NGATHER];
  #define NTEMP (NGATHER)
  su3_matrix **temp = new_field( NTEMP );

  su3_matrix tmat1,tmat2;
  complex cc;

  /* Calculate the temporal/electric components */  
  for(icomp=FS_XT;icomp<=FS_ZT;icomp++){
    switch(icomp){
      case FS_XT: dir[0]=XUP; dir[1]=TUP; break;
      case FS_YT: dir[0]=YUP; dir[1]=TUP; break;
      case FS_ZT: dir[0]=ZUP; dir[1]=TUP; break;
    }
    setup_blocked_dirs( dir, dirb ,diro );
    #define LINK0 link[dir[0]]
    #define LINK1 link_lf[dir[0]]

    /* Consistency for improved readability: gathers in 
       +- dir[0] direction have tag[0] and gen_pt[0], use temp[0] 
       +- dir[1] direction have tag[1] and gen_pt[1], use temp[1] */ 

    /* Plaquette in +dir[0] +dir[1] direction */
    FORALLSITES(i,s) 
    IF_BLOCKED(s, block_stride)       
    IF_LWR_BDRY(s) 
    {
      mult_su3_na( &(LINK0[i]), &(LINK1[i]), &tmat1 );
      su3_adjoint( &tmat1, &tmat2 );
      sub_su3_matrix(  &tmat1, &tmat2, &(fstrength[icomp][i]) );
    }

    /* Plaquette in -dir[0] +dir[1] direction */
    
    FORALLSITES(i,s) 
    IF_BLOCKED(s, block_stride)       
    IF_LWR_BDRY(s) 
      mult_su3_an(  &(LINK1[i]), &(LINK0[i]), &temp[0][i] );

    ig = 0;
    // request request dir0-dir1-dir0 staple from direction diro[0]
    tag[ig] = start_gather_field( temp[ig], sizeof(su3_matrix),
                                  diro[ig], EVENANDODD, gen_pt[ig] );
    wait_gather(tag[ig]);
 
    FORALLSITES(i,s) 
    IF_BLOCKED(s, block_stride)       
    IF_LWR_BDRY(s) 
    {
      su3_adjoint( (su3_matrix *)(gen_pt[0][i]), &tmat2 );
      add_su3_matrix( &(fstrength[icomp][i]), (su3_matrix *)(gen_pt[0][i]),
                      &(fstrength[icomp][i]) );
      sub_su3_matrix( &(fstrength[icomp][i]), &tmat2,
                      &(fstrength[icomp][i]) );
    }
    cleanup_gather(tag[0]);

    /* Consistency for improved readability: gathers in 
       +- dir[0] direction have tag[0] and gen_pt[0], use temp[0] 
       +- dir[1] direction have tag[1] and gen_pt[1], use temp[1] */ 

    /* Make traceless */
    FORALLSITES(i,s) 
    IF_BLOCKED(s, block_stride)       
    IF_LWR_BDRY(s) 
    {
      cc = trace_su3( &(fstrength[icomp][i]) );
      CMULREAL(cc,0.33333333333333333,cc);
      for(j=0;j<3;j++)
        CSUB( fstrength[icomp][i].e[j][j],cc,
              fstrength[icomp][i].e[j][j]);

    }

    #undef LINK0
    #undef LINK1
  } /* component */

  // deallocate temporary storage
  destroy_field( &temp );
  #undef NTEMP
  #undef NGATHER
}

static void 
make_dropped_temporal_fieldstrength_upr_bdry ( su3_matrix **link, 
                                           su3_matrix **link_lf, 
                                           su3_matrix **fstrength ) {

  register int i, j, icomp;
  register int stride = block_stride;
  #define NGATHER 1
  int dir[2] = {NODIR,NODIR}, 
      dirb[2] = {NODIR,NODIR}, 
      diro[2] = {NODIR,NODIR};
  register int ig;
  register site *s;
  msg_tag *tag[NGATHER];
  #define NTEMP (NGATHER)
  su3_matrix **temp = new_field( NTEMP );

  su3_matrix tmat1,tmat2;
  complex cc;

  /* Calculate the temporal/electric components */  
  for(icomp=FS_XT;icomp<=FS_ZT;icomp++){
    switch(icomp){
      case FS_XT: dir[0]=XUP; dir[1]=TUP; break;
      case FS_YT: dir[0]=YUP; dir[1]=TUP; break;
      case FS_ZT: dir[0]=ZUP; dir[1]=TUP; break;
    }
    setup_blocked_dirs( dir, dirb ,diro );
    #define LINK0 link[dir[0]]
    #define LINK1 link_lf[dir[0]]

    /* Consistency for improved readability: gathers in 
       +- dir[0] direction have tag[0] and gen_pt[0], use temp[0] 
       +- dir[1] direction have tag[1] and gen_pt[1], use temp[1] */ 

    /* Plaquette in +dir[0] +dir[1] direction */
    FORALLSITES(i,s) 
    IF_BLOCKED(s, block_stride)       
    IF_UPR_BDRY(s) 
    {
      mult_su3_na( &(LINK0[i]), &(LINK1[i]), &tmat1 );
      su3_adjoint( &tmat1, &tmat2 );
      sub_su3_matrix(  &tmat1, &tmat2, &(fstrength[icomp][i]) );
    }

    /* Plaquette in -dir[0] +dir[1] direction */
    
    FORALLSITES(i,s) 
    IF_BLOCKED(s, block_stride)       
    IF_UPR_BDRY(s) 
      mult_su3_an(  &(LINK1[i]), &(LINK0[i]), &temp[0][i] );

    ig = 0;
    // request request dir0-dir1-dir0 staple from direction diro[0]
    tag[ig] = start_gather_field( temp[ig], sizeof(su3_matrix),
                                  diro[ig], EVENANDODD, gen_pt[ig] );
    wait_gather(tag[ig]);
 
    FORALLSITES(i,s) 
    IF_BLOCKED(s, block_stride)       
    IF_UPR_BDRY(s) 
    {
      su3_adjoint( (su3_matrix *)(gen_pt[0][i]), &tmat2 );
      add_su3_matrix( &(fstrength[icomp][i]), (su3_matrix *)(gen_pt[0][i]),
                      &(fstrength[icomp][i]) );
      sub_su3_matrix( &(fstrength[icomp][i]), &tmat2,
                      &(fstrength[icomp][i]) );
    }
    cleanup_gather(tag[0]);

    /* Consistency for improved readability: gathers in 
       +- dir[0] direction have tag[0] and gen_pt[0], use temp[0] 
       +- dir[1] direction have tag[1] and gen_pt[1], use temp[1] */ 

    /* Make traceless */
    FORALLSITES(i,s) 
    IF_BLOCKED(s, block_stride)       
    IF_UPR_BDRY(s) 
    {
      cc = trace_su3( &(fstrength[icomp][i]) );
      CMULREAL(cc,0.33333333333333333,cc);
      for(j=0;j<3;j++)
        CSUB( fstrength[icomp][i].e[j][j],cc,
              fstrength[icomp][i].e[j][j]);

    }

    #undef LINK0
    #undef LINK1
  } /* component */

  // deallocate temporary storage
  destroy_field( &temp );
  #undef NTEMP
  #undef NGATHER
}

void 
make_fieldstrength_bdry ( field_offset link_src, 
                          su3_matrix **link_last_flow, 
                          su3_matrix **fstrength ) {

  register int i, dir;
  register site *s;
  #define NLINK 3
  su3_matrix **link = new_links_from_site( ACTIVE, NLINK );

  make_spatial_fieldstrength_region( BOUNDARY, link, fstrength );
  make_spatial_fieldstrength_region( BOUNDARY, link_last_flow, fstrength + 1 * N_LAST_FLOW );

  make_dropped_temporal_fieldstrength_bdry ( link, link_last_flow, fstrength + 1 * N_LAST_FLOW );

  make_improved_spatial_fieldstrength_region( BOUNDARY, link, fstrength + 3 * N_LAST_FLOW );
  make_improved_spatial_fieldstrength_region( BOUNDARY, link_last_flow, fstrength + 4* N_LAST_FLOW  );

  // deallocate temporary spatial links
  destroy_field( &link );
  #undef NLINK
}

void 
make_fieldstrength_lwr_bdry ( field_offset link_src, 
                          su3_matrix **link_last_flow, 
                          su3_matrix **fstrength ) {

  register int i, dir;
  register site *s;
  #define NLINK 3
  su3_matrix **link = new_links_from_site( ACTIVE, NLINK );

  make_spatial_fieldstrength_region( LOWER_BOUNDARY, link, fstrength );
  make_spatial_fieldstrength_region( LOWER_BOUNDARY, link_last_flow, fstrength + 1 * N_LAST_FLOW );

  make_dropped_temporal_fieldstrength_lwr_bdry ( link, link_last_flow, fstrength + 1 * N_LAST_FLOW );

  make_improved_spatial_fieldstrength_region( LOWER_BOUNDARY, link, fstrength + 3 * N_LAST_FLOW );
  make_improved_spatial_fieldstrength_region( LOWER_BOUNDARY, link_last_flow, fstrength + 4* N_LAST_FLOW  );

  // deallocate temporary spatial links
  destroy_field( &link );
  #undef NLINK
}

void 
make_fieldstrength_upr_bdry ( field_offset link_src, 
                          su3_matrix **link_last_flow, 
                          su3_matrix **fstrength ) {

  register int i, dir;
  register site *s;
  #define NLINK 3
  su3_matrix **link = new_links_from_site( ACTIVE, NLINK );

  make_spatial_fieldstrength_region( UPPER_BOUNDARY, link, fstrength );
  make_spatial_fieldstrength_region( UPPER_BOUNDARY, link_last_flow, fstrength + 1 * N_LAST_FLOW );

  make_dropped_temporal_fieldstrength_upr_bdry ( link, link_last_flow, fstrength + 1 * N_LAST_FLOW );

  make_improved_spatial_fieldstrength_region( UPPER_BOUNDARY, link, fstrength + 3 * N_LAST_FLOW );
  make_improved_spatial_fieldstrength_region( UPPER_BOUNDARY, link_last_flow, fstrength + 4* N_LAST_FLOW  );

  // deallocate temporary spatial links
  destroy_field( &link );
  #undef NLINK
}

void 
make_fieldstrength_half ( field_offset link_src, 
                          su3_matrix **fstrength ) {

  register int i, dir;
  register site *s;
  #define NLINK 3
  su3_matrix **link = new_links_from_site( ACTIVE, NLINK );
  su3_matrix **link_half = new_half_links();

  make_spatial_fieldstrength_region( LOWER_BULK, link_half, fstrength );
  make_improved_spatial_fieldstrength_region( LOWER_BULK, link_half, fstrength + 9 );

  make_dropped_temporal_fieldstrength_half ( link, link_half, fstrength );

  // deallocate temporary spatial links
  destroy_field( &link_half );
  destroy_field( &link );
  #undef NLINK
}
#endif