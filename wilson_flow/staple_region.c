/************************ staple.c ********************************/
/* Computes the wilson and symanzik staples at each link          */

/* Definitions, files, and prototypes */
#include "wilson_flow_includes.h"
#define SU3_PT(field) (su3_matrix*)F_PT(s,field)

#ifdef SPHALERON

static int region = FULLVOL;

// Computes 1x1 staples
void
wilson_staple_region(int dir0, int dir1, field_offset lnk0,
              field_offset lnk1, Real coeff, field_offset stp) {
  register int i;
  register site *s;
  msg_tag *tag0, *tag1, *tag2;
  su3_matrix tmat1, tmat2;

  int dir[2] = {dir0,dir1}, 
      dirb[2] = {NODIR,NODIR}, 
      diro[2] = {NODIR,NODIR};

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
  if ( dir[0] == TUP ) {
    dirb[0] =      dir[0];  diro[0] = OPP_DIR(dir[0]);
  }
  if ( dir[1] == TUP ) {
    dirb[1] =      dir[1];  diro[1] = OPP_DIR(dir[1]);
  }

  /* Consistency for improved readability: gathers in 
     +- dir[0] direction have tag[0] and gen_pt[0], use temp[0] 
     +- dir[1] direction have tag[1] and gen_pt[1], use temp[1] */ 

  // Get blocked_link[dir[1]] from direction dirb[0]
  tag0 = start_gather_site( lnk1, sizeof(su3_matrix), 
                            dirb[0], EVENANDODD, gen_pt[0]);

  // Get blocked_link[dir[0]] from direction dirb[1]
  tag1 = start_gather_site( lnk0, sizeof(su3_matrix), 
                            dirb[1], EVENANDODD, gen_pt[1]);

  // Start working on the lower staple while we wait for the gathers
  // The lower staple is prepared at x-dirb[1] and stored in tempmat[0],
  // then gathered to x
  FORALLSITES(i, s)
  IF_BLOCKED(s, block_stride)       
  // IF_REGION(s, region)
    mult_su3_an(SU3_PT(lnk1), SU3_PT(lnk0), tempmat[0]+i);

  wait_gather(tag0);

  // Finish lower staple
  FORALLSITES(i, s)
  IF_BLOCKED(s, block_stride)       
  // IF_REGION(s, region) 
  {
    mult_su3_nn(tempmat[0]+i, (su3_matrix *)gen_pt[0][i], &tmat1);
    su3mat_copy(&tmat1, tempmat[0]+i);
  }

  // Gather staple from direction diro[1] to "home" site
  tag2 = start_gather_field(tempmat[0], sizeof(su3_matrix),
                            diro[1], EVENANDODD, gen_pt[2]);

  wait_gather(tag1);

  // Calculate upper staple, add it
  FORALLSITES(i, s)
  IF_BLOCKED(s, block_stride)       
  IF_REGION(s, region) 
  {
    mult_su3_nn(SU3_PT(lnk1), (su3_matrix *)gen_pt[1][i], &tmat1);
    mult_su3_na(&tmat1, (su3_matrix *)gen_pt[0][i], &tmat2);
    scalar_mult_add_su3_matrix(SU3_PT(stp), &tmat2, coeff, SU3_PT(stp));
  }

  // Add the lower staple and scale the entire 1x1 contribution
  wait_gather(tag2);
  FORALLSITES(i, s) 
  IF_BLOCKED(s, block_stride)       
  IF_REGION(s, region)
    scalar_mult_add_su3_matrix(SU3_PT(stp), (su3_matrix *)gen_pt[2][i],
                               coeff, SU3_PT(stp));
  cleanup_gather(tag0);
  cleanup_gather(tag1);
  cleanup_gather(tag2);
}

// Computes 1x2 staples
#define TM(tagi) tempmat[tagi-2]
void
symanzik1x2_staple_region(int dir0, int dir1, field_offset lnk0,
                   field_offset lnk1, Real coeff, field_offset stp) {
  register int i;
  register site *s;
  msg_tag *tag0, *tag1, *tag2, *tag3, *tag4, *tag5, *tag6, *tag7, *tag8;
  su3_matrix tmat1, tmat2;

  int dir[2] = {dir0,dir1}, 
      dirb[2] = {NODIR,NODIR}, 
      diro[2] = {NODIR,NODIR};

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
  if ( dir[0] == TUP ) {
    dirb[0] =      dir[0];  diro[0] = OPP_DIR(dir[0]);
  }
  if ( dir[1] == TUP ) {
    dirb[1] =      dir[1];  diro[1] = OPP_DIR(dir[1]);
  }

  /* Consistency for improved readability: gathers in 
     +- dir[0] direction have tag[0] and gen_pt[0], use temp[0] 
     +- dir[1] direction have tag[1] and gen_pt[1], use temp[1] */ 

  // Get blocked_link[dir[1]] from direction dirb[0]
  tag0 = start_gather_site( lnk1, sizeof(su3_matrix), 
                            dirb[0], EVENANDODD, gen_pt[0]);

  // Get blocked_link[dir[0]] from direction dirb[1]
  tag1 = start_gather_site( lnk0, sizeof(su3_matrix), 
                            dirb[1], EVENANDODD, gen_pt[1]);

  wait_gather(tag0);
  wait_gather(tag1);

  //Calculate upper 1x1 staple
  FORALLSITES(i, s)
  IF_BLOCKED(s, block_stride)
  // IF_REGION(s, region) 
  {
    mult_su3_nn(SU3_PT(lnk1), (su3_matrix *)gen_pt[1][i], &tmat1);
    mult_su3_na(&tmat1, (su3_matrix *)gen_pt[0][i], TM(2)+i);
  }

  //Gather upper staple from direction dirb[1]
  tag2 = start_gather_field(TM(2), sizeof(su3_matrix), 
                            dirb[1], EVENANDODD, gen_pt[2]);

  //Calculate lower 1x1 staple
  FORALLSITES(i, s)
  IF_BLOCKED(s, block_stride)
  // IF_REGION(s, region) 
  {
    mult_su3_an(SU3_PT(lnk1), SU3_PT(lnk0), &tmat1);
    mult_su3_nn(&tmat1, (su3_matrix *)gen_pt[0][i], TM(3)+i);
  }

  //Gather lower staple from direction diro[1]
  tag3 = start_gather_field(TM(3), sizeof(su3_matrix), 
                            diro[1], EVENANDODD, gen_pt[3]);

  //Start construction of 2(dirb[0])x1(dirb[1]) staple protruding in dirb[0]
  FORALLSITES(i, s)
  IF_BLOCKED(s, block_stride)
  // IF_REGION(s, region) 
  {
    mult_su3_nn(SU3_PT(lnk0), (su3_matrix *)gen_pt[0][i], &tmat1);
    mult_su3_na(&tmat1, (su3_matrix *)gen_pt[1][i], TM(4)+i);
  }

  //Gather dirb[0] 1x1 protrusion from direction dirb[0]
  tag4 = start_gather_field(TM(4), sizeof(su3_matrix), 
                            dirb[0], EVENANDODD, gen_pt[4]);

  //Start construction of 2(dirb[0])x1(dirb[1]) staple protruding in diro[0]
  FORALLSITES(i, s)
  IF_BLOCKED(s, block_stride)
  // IF_REGION(s, region) 
  {
    mult_su3_nn(SU3_PT(lnk1), (su3_matrix *)gen_pt[1][i], &tmat1);
    mult_su3_an(SU3_PT(lnk0), &tmat1, TM(5)+i);
  }

  //Gather diro[0] 1x1 protrusion from direction diro[0]
  tag5 = start_gather_field(TM(5), sizeof(su3_matrix), 
                            diro[0], EVENANDODD, gen_pt[5]);

  //Continue construction of lower 1(dirb[0])x2(dirb[1]) staple
  wait_gather(tag3);
  FORALLSITES(i, s)
  IF_BLOCKED(s, block_stride)
  // IF_REGION(s, region) 
  {
    mult_su3_an(SU3_PT(lnk1), (su3_matrix *)gen_pt[3][i], &tmat1);
    mult_su3_nn(&tmat1, (su3_matrix *)gen_pt[0][i], TM(6)+i);
  }

  //Gather 1(dirb[0])x2(dirb[1]) lower staple from direction diro[1]
  tag6 = start_gather_field(TM(6), sizeof(su3_matrix), 
                            diro[1], EVENANDODD, gen_pt[6]);

  //Continue construction of 2(dirb[0])x1(dirb[1]) lower staple protruding in dirb[0]
  wait_gather(tag4);
  FORALLSITES(i, s)
  IF_BLOCKED(s, block_stride)
  // IF_REGION(s, region) 
  {
    mult_su3_nn(SU3_PT(lnk0), (su3_matrix *)gen_pt[4][i], &tmat1);
    mult_su3_an(SU3_PT(lnk1), &tmat1, TM(7)+i);
  }

  //Gather 2(dirb[0])x1(dirb[1]) right, lower staple from direction diro[1]
  tag7 = start_gather_field(TM(7), sizeof(su3_matrix), 
                            diro[1], EVENANDODD, gen_pt[7]);

  //Continue construction of 2(dirb[0])x1(dirb[1]) lower staple protruding in diro[0]
  wait_gather(tag5);
  FORALLSITES(i, s)
  IF_BLOCKED(s, block_stride)
  // IF_REGION(s, region) 
  {
    mult_su3_an((su3_matrix *)gen_pt[5][i], SU3_PT(lnk0), &tmat1);
    mult_su3_nn(&tmat1, (su3_matrix *)gen_pt[0][i], TM(8)+i);
  }

  //Gather 2(dirb[0])x1(dirb[1]) left, lower staple from direction diro[1]
  tag8 = start_gather_field(TM(8), sizeof(su3_matrix), 
                            diro[1], EVENANDODD, gen_pt[8]);

  wait_gather(tag2);
  FORALLSITES(i, s)
  IF_BLOCKED(s, block_stride)
  IF_REGION(s, region) 
  {
    //Finish and add upper 1(dirb[0])x2(dirb[1]) staple
    mult_su3_nn(SU3_PT(lnk1), (su3_matrix *)gen_pt[2][i], &tmat1);
    mult_su3_na(&tmat1, (su3_matrix *)gen_pt[0][i], &tmat2);
    scalar_mult_add_su3_matrix(SU3_PT(stp), &tmat2, coeff, SU3_PT(stp));

    //Finish and add upper 2(dirb[0])x1(dirb[1]) staple protruding in dirb[0]
    mult_su3_nn(SU3_PT(lnk1), (su3_matrix *)gen_pt[1][i], &tmat1);
    mult_su3_na(&tmat1, (su3_matrix *)gen_pt[4][i], &tmat2);
    scalar_mult_add_su3_matrix(SU3_PT(stp), &tmat2, coeff, SU3_PT(stp));

    //Finish and add upper 2(dirb[0])x1(dirb[1]) staple protruding in diro[0]
    mult_su3_nn((su3_matrix *)gen_pt[5][i], (su3_matrix *)gen_pt[1][i], &tmat1);
    mult_su3_na(&tmat1, (su3_matrix *)gen_pt[0][i], &tmat2);
    scalar_mult_add_su3_matrix(SU3_PT(stp), &tmat2, coeff, SU3_PT(stp));
  }

  //Finish lower 1(dir1)x2(dir2) staple
  wait_gather(tag6);
  FORALLSITES(i, s)
  IF_BLOCKED(s, block_stride)
  IF_REGION(s, region)
    scalar_mult_add_su3_matrix(SU3_PT(stp), (su3_matrix *)gen_pt[6][i], coeff,
                               SU3_PT(stp));

  //Add lower 2(dir1)x1(dir2) staple protruding in +dir1
  wait_gather(tag7);
  FORALLSITES(i, s)
  IF_BLOCKED(s, block_stride)
  IF_REGION(s, region)
    scalar_mult_add_su3_matrix(SU3_PT(stp), (su3_matrix *)gen_pt[7][i], coeff,
                               SU3_PT(stp));

  //Add lower 2(dir1)x1(dir2) staple protruding in -dir1
  wait_gather(tag8);
  FORALLSITES(i, s)
  IF_BLOCKED(s, block_stride)
  IF_REGION(s, region)
    scalar_mult_add_su3_matrix(SU3_PT(stp), (su3_matrix *)gen_pt[8][i], coeff,
                               SU3_PT(stp));

  //clean up
  cleanup_gather(tag0);
  cleanup_gather(tag1);
  cleanup_gather(tag2);
  cleanup_gather(tag3);
  cleanup_gather(tag4);
  cleanup_gather(tag5);
  cleanup_gather(tag6);
  cleanup_gather(tag7);
  cleanup_gather(tag8);
}

/* Correction to Symanzik flow -- Zeuthen flow,
   Ramos, Sint, 1408.05552 */
#define TMZ(tagi) tempmat[tagi]
void
zeuthen_correction_region( int dir0, field_offset lnk0,
                    Real coeff, field_offset stp ) {

  register int i;
  register site *s;
  Real weight = 1.0 - 2.0 * coeff;
  msg_tag *tag0, *tag1;
  su3_matrix tmat1, tmat2;

  int dir[1] = {dir0}, 
      dirb[1] = {NODIR}, 
      diro[1] = {NODIR};

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
  if ( dir[0] == TUP ) {
    dirb[0] =      dir[0];  diro[0] = OPP_DIR(dir[0]);
  }

  // Get blocked_staple from direction dirb[0]
  tag0 = start_gather_site( stp, sizeof(su3_matrix), 
                            dirb[0], EVENANDODD, gen_pt[0]);

  // Calculate blocked laplacian's downward contribution at diro[0]
  FORALLSITES(i, s)
  IF_BLOCKED(s, block_stride)
  IF_REGION(s, region) 
  {
    mult_su3_an( SU3_PT(lnk0), SU3_PT(stp), &tmat1);
    mult_su3_nn(&tmat1, SU3_PT(lnk0), TMZ(1)+i);
  }

  // Get blocked laplacian's downward contribution from direction diro[1]
  tag1 = start_gather_field( TMZ(1), sizeof(su3_matrix), 
                            diro[0], EVENANDODD, gen_pt[1]);

  wait_gather(tag0);

  // Calculate blocked laplacian's upward contribution
  FORALLSITES(i, s)
  IF_BLOCKED(s, block_stride)
  // IF_REGION(s, region) 
  {
    mult_su3_nn( SU3_PT(lnk0), (su3_matrix *)gen_pt[0][i], &tmat1);
    mult_su3_na(&tmat1, SU3_PT(lnk0), TMZ(0)+i);
  }

  wait_gather(tag1);

  //Combine blocked laplacian's shifted contributions with weighted on-site staple
  FORALLSITES(i, s)
  IF_BLOCKED(s, block_stride)
  IF_REGION(s, region) 
  {
    scalar_mult_su3_matrix( SU3_PT(stp), weight, &tmat1 );
    add_su3_matrix( TMZ(0)+i, (su3_matrix *)gen_pt[1][i], &tmat2);
    scalar_mult_add_su3_matrix( &tmat1, &tmat2, coeff, SU3_PT(stp) );
  }

  cleanup_gather(tag0);
  cleanup_gather(tag1);

}

/* Compiles all contributions to the staple */
void
staple_region( int region_flag )
{
  register int i;
  register site *s;
  int dir0, dir1;
  Real coeff1x1, coeff1x2, coeffz=1.0/12.0;
  field_offset lnk0, lnk1, stp;
#ifdef ANISOTROPY
  Real tmp;
#endif

  /* Pick the coefficients for each loop contributing to the staple */
  if( stapleflag==WILSON ) {
    coeff1x1 = 1.0;
    coeff1x2 = 0.0;
  }
  else { /* stapleflag==SYMANZIK || ZEUTHEN */
    coeff1x1 = 5.0/3.0;
    coeff1x2 = -1.0/12.0;
  }

  /* Pick the direction for the missing link */
  FORALLDIRSUP(dir0, region) {
    // node0_printf("dirs %d\n",dir1);

    /* Clear any previously stored staple */
    FORALLSITES(i, s)
      clear_su3mat(&(s->staple[dir0]));

    /* Pick the second direction for the staple */
    FORALLUPDIRBUT(dir0, dir1) {
      if ( region != FULLVOL && dir1 == TUP ) 
        continue;
      lnk0 = F_OFFSET(link[dir0]);
      lnk1 = F_OFFSET(link[dir1]);
      stp = F_OFFSET(staple[dir0]);

#ifdef ANISOTROPY
      if(dir1==TUP)
        tmp=ani*ani;
      else
        tmp=1.0;

      /* Adds 1x1 (Wilson) contributions */
      wilson_staple_region(dir0, dir1, lnk0, lnk1, tmp*coeff1x1, stp);

      /* Adds 1x2 (Symanzik tree level) contributions -- both
         for Symanzik and Zeuthen flow */
      if( stapleflag!=WILSON )
        symanzik1x2_staple_region(dir0, dir1, lnk0, lnk1, tmp*coeff1x2, stp);
#else

      /* Adds 1x1 (Wilson) contributions */
      wilson_staple_region(dir0, dir1, lnk0, lnk1, coeff1x1, stp);

      /* Adds 1x2 (Symanzik tree level) contributions -- both
         for Symanzik and Zeuthen flow */
      if( stapleflag!=WILSON )
        symanzik1x2_staple_region(dir0, dir1, lnk0, lnk1, coeff1x2, stp);
#endif

    // node0_printf("dirs %d %d\n",dir1,dir2);
    // dumpmat( &(lattice->staple[XUP]) );
    } /* end: dir1 loop */

    if( stapleflag==ZEUTHEN ) 
      zeuthen_correction_region(dir0,lnk0,coeffz,stp);

  } /* end: dir0 loop */
}

#endif