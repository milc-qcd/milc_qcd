/************************ staple.c ********************************/
/* Computes the wilson and symanzik staples at each link          */

/* Definitions, files, and prototypes */
#include "wilson_flow_includes.h"
#define SU3_PT(field) (su3_matrix*)F_PT(s,field)

// Computes 1x1 staples
void
wilson_staple(int dir1, int dir2, field_offset lnk1,
              field_offset lnk2, Real coeff, field_offset stp)
{
  register int i;
  register site *s;
  msg_tag *tag0, *tag1, *tag2;
  su3_matrix tmat1, tmat2;

  // Get blocked_link[dir2] from direction dir1
  tag0 = start_gather_site(lnk2, sizeof(su3_matrix), dir1,
                           EVENANDODD, gen_pt[0]);

  // Get blocked_link[dir1] from direction dir2
  tag1 = start_gather_site(lnk1, sizeof(su3_matrix), dir2,
                           EVENANDODD, gen_pt[1]);

  // Start working on the lower staple while we wait for the gathers
  // The lower staple is prepared at x-dir2 and stored in tempmat[0],
  // then gathered to x
  FORALLSITES(i, s)
    mult_su3_an(SU3_PT(lnk2), SU3_PT(lnk1), tempmat[0]+i);

  wait_gather(tag0);

  // Finish lower staple
  FORALLSITES(i, s) {
    mult_su3_nn(tempmat[0]+i, (su3_matrix *)gen_pt[0][i], &tmat1);
    su3mat_copy(&tmat1, tempmat[0]+i);
  }

  // Gather staple from direction -dir2 to "home" site
  tag2 = start_gather_field(tempmat[0], sizeof(su3_matrix),
                            OPP_DIR(dir2), EVENANDODD, gen_pt[2]);

  wait_gather(tag1);

  // Calculate upper staple, add it
  FORALLSITES(i, s) {
    mult_su3_nn(SU3_PT(lnk2), (su3_matrix *)gen_pt[1][i], &tmat1);
    mult_su3_na(&tmat1, (su3_matrix *)gen_pt[0][i], &tmat2);
    scalar_mult_add_su3_matrix(SU3_PT(stp), &tmat2, coeff, SU3_PT(stp));
  }

  // Add the lower staple and scale the entire 1x1 contribution
  wait_gather(tag2);
  FORALLSITES(i, s)
    scalar_mult_add_su3_matrix(SU3_PT(stp), (su3_matrix *)gen_pt[2][i],
                               coeff, SU3_PT(stp));

  cleanup_gather(tag0);
  cleanup_gather(tag1);
  cleanup_gather(tag2);
}

// Computes 1x2 staples
#define TM(tagi) tempmat[tagi-2]
void
symanzik1x2_staple(int dir1, int dir2, field_offset lnk1,
                   field_offset lnk2, Real coeff, field_offset stp)
{	
  register int i;
  register site *s;
  msg_tag *tag0, *tag1, *tag2, *tag3, *tag4, *tag5, *tag6, *tag7, *tag8;
  su3_matrix tmat1, tmat2;

  // Get blocked_link[dir2] from direction dir1
  tag0 = start_gather_site(lnk2, sizeof(su3_matrix), dir1,
                           EVENANDODD, gen_pt[0]);

  // Get blocked_link[dir1] from direction dir2
  tag1 = start_gather_site(lnk1, sizeof(su3_matrix), dir2,
                           EVENANDODD, gen_pt[1]);

  wait_gather(tag0);
  wait_gather(tag1);
	
  //Calculate upper 1x1 staple
  FORALLSITES(i, s) {
    mult_su3_nn(SU3_PT(lnk2), (su3_matrix *)gen_pt[1][i], &tmat1);
    mult_su3_na(&tmat1, (su3_matrix *)gen_pt[0][i], TM(2)+i);
  }

  //Gather upper staple from direction dir2
  tag2 = start_gather_field(TM(2), sizeof(su3_matrix), dir2,
                            EVENANDODD, gen_pt[2]);

  //Calculate lower 1x1 staple
  FORALLSITES(i, s) {
    mult_su3_an(SU3_PT(lnk2), SU3_PT(lnk1), &tmat1);
    mult_su3_nn(&tmat1, (su3_matrix *)gen_pt[0][i], TM(3)+i);
  }

  //Gather lower staple from direction -dir2
  tag3 = start_gather_field(TM(3), sizeof(su3_matrix), OPP_DIR(dir2),
                            EVENANDODD, gen_pt[3]);
	
  //Start construction of 2(dir1)x1(dir2) staple protruding in +dir1
  FORALLSITES(i, s) {
    mult_su3_nn(SU3_PT(lnk1), (su3_matrix *)gen_pt[0][i], &tmat1);
    mult_su3_na(&tmat1, (su3_matrix *)gen_pt[1][i], TM(4)+i);
  }

  //Gather +dir1 1x1 protrusion from direction dir1
  tag4 = start_gather_field(TM(4), sizeof(su3_matrix), dir1,
                            EVENANDODD, gen_pt[4]);

  //Start construction of 2(dir1)x1(dir2) staple protruding in -dir1
  FORALLSITES(i, s) {
    mult_su3_nn(SU3_PT(lnk2), (su3_matrix *)gen_pt[1][i], &tmat1);
    mult_su3_an(SU3_PT(lnk1), &tmat1, TM(5)+i);
  }

  //Gather -dir1 1x1 protrusion from direction -dir1
  tag5 = start_gather_field(TM(5), sizeof(su3_matrix), OPP_DIR(dir1),
                            EVENANDODD, gen_pt[5]);

  //Continue construction of lower 1(dir1)x2(dir2) staple
  wait_gather(tag3);
  FORALLSITES(i, s) {
    mult_su3_an(SU3_PT(lnk2), (su3_matrix *)gen_pt[3][i], &tmat1);
    mult_su3_nn(&tmat1, (su3_matrix *)gen_pt[0][i], TM(6)+i);
  }

  //Gather 1(dir1)x2(dir2) lower staple from direction -dir2
  tag6 = start_gather_field(TM(6), sizeof(su3_matrix), OPP_DIR(dir2),
                            EVENANDODD, gen_pt[6]);

  //Continue construction of 2(dir1)x1(dir2) lower staple protruding in +dir1
  wait_gather(tag4);
  FORALLSITES(i, s) {
    mult_su3_nn(SU3_PT(lnk1), (su3_matrix *)gen_pt[4][i], &tmat1);
    mult_su3_an(SU3_PT(lnk2), &tmat1, TM(7)+i);
  }

  //Gather 2(dir1)x1(dir2) right, lower staple from direction -dir2
  tag7 = start_gather_field(TM(7), sizeof(su3_matrix), OPP_DIR(dir2),
                            EVENANDODD, gen_pt[7]);

  //Continue construction of 2(dir1)x1(dir2) lower staple protruding in -dir1
  wait_gather(tag5);
  FORALLSITES(i, s) {
    mult_su3_an((su3_matrix *)gen_pt[5][i], SU3_PT(lnk1), &tmat1);
    mult_su3_nn(&tmat1, (su3_matrix *)gen_pt[0][i], TM(8)+i);
  }
	
  //Gather 2(dir1)x1(dir2) left, lower staple from direction -dir2
  tag8 = start_gather_field(TM(8), sizeof(su3_matrix), OPP_DIR(dir2),
                            EVENANDODD, gen_pt[8]);

  wait_gather(tag2);
  FORALLSITES(i, s) {
    //Finish and add upper 1(dir1)x2(dir2) staple
    mult_su3_nn(SU3_PT(lnk2), (su3_matrix *)gen_pt[2][i], &tmat1);
    mult_su3_na(&tmat1, (su3_matrix *)gen_pt[0][i], &tmat2);
    scalar_mult_add_su3_matrix(SU3_PT(stp), &tmat2, coeff, SU3_PT(stp));

    //Finish and add upper 2(dir1)x1(dir2) staple protruding in +dir1
    mult_su3_nn(SU3_PT(lnk2), (su3_matrix *)gen_pt[1][i], &tmat1);
    mult_su3_na(&tmat1, (su3_matrix *)gen_pt[4][i], &tmat2);
    scalar_mult_add_su3_matrix(SU3_PT(stp), &tmat2, coeff, SU3_PT(stp));

    //Finish and add upper 2(dir1)x1(dir2) staple protruding in -dir1
    mult_su3_nn((su3_matrix *)gen_pt[5][i], (su3_matrix *)gen_pt[1][i], &tmat1);
    mult_su3_na(&tmat1, (su3_matrix *)gen_pt[0][i], &tmat2);
    scalar_mult_add_su3_matrix(SU3_PT(stp), &tmat2, coeff, SU3_PT(stp));
  }

  //Finish lower 1(dir1)x2(dir2) staple
  wait_gather(tag6);
  FORALLSITES(i, s)
    scalar_mult_add_su3_matrix(SU3_PT(stp), (su3_matrix *)gen_pt[6][i], coeff, 
                               SU3_PT(stp));

  //Add lower 2(dir1)x1(dir2) staple protruding in +dir1
  wait_gather(tag7);
  FORALLSITES(i, s)
    scalar_mult_add_su3_matrix(SU3_PT(stp), (su3_matrix *)gen_pt[7][i], coeff, 
                               SU3_PT(stp));

  //Add lower 2(dir1)x1(dir2) staple protruding in -dir1
  wait_gather(tag8);
  FORALLSITES(i, s) 
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

/* Compiles all contributions to the staple */
void
staple() 
{
  register int i;
  register site *s;
  int dir1, dir2;
  Real coeff1x1, coeff1x2;
  field_offset lnk1, lnk2, stp;

  /* Pick the coefficients for each loop contributing to the staple */
  if( stapleflag==WILSON ) {
    coeff1x1 = 1.0;
    coeff1x2 = 0.0;
  }
  else { /* stapleflag==SYMANZIK */
    coeff1x1 = 5.0/3.0;
    coeff1x2 = -1.0/12.0;
  }

  /* Pick the direction for the missing link */
  FORALLUPDIR(dir1) {

    /* Clear any previously stored staple */
    FORALLSITES(i, s)
      clear_su3mat(&(s->staple[dir1]));

    /* Pick the second direction for the staple */
    FORALLUPDIRBUT(dir1, dir2) {
      lnk1 = F_OFFSET(link[dir1]);
      lnk2 = F_OFFSET(link[dir2]);
      stp = F_OFFSET(staple[dir1]);

      /* Adds 1x1 (Wilson) contributions */
      wilson_staple(dir1, dir2, lnk1, lnk2, coeff1x1, stp);

      /* Adds 1x2 (Symanzik tree level) contributions */
      if( stapleflag!=WILSON )
        symanzik1x2_staple(dir1, dir2, lnk1, lnk2, coeff1x2, stp);

    } /* end: dir2 loop */
  } /* end: dir1 loop */
}
