/****************** ploop_staple.c ************************************/
/* MIMD version 7 */
/* NOT MAINTAINED.  TEST BEFORE USE! */
/* evaluate fuzzy Polyakov loops (including weighted contributions from
   time-oriented staples).  This version uses general_gathers. */
/* It assumes that nt is even.  Actually, all the code does.  */
/* Mods for fuzzy loops by C. DeTar 7/28/93 */
/* this version saves results in the ploop_fuzz structure */

#include "generic_includes.h"

#define TRUE 1
#define FALSE 0

complex ploop_staple(Real alpha_fuzz) {
  register int i,j,t;
  register site *st;
  msg_tag *tag;
  complex sum;
  complex plp_staple,plp_fuzz;
  int d[4];
  
  /* First we compute the sum of the "left" and "right"
     time-oriented staples.
     This code is cribbed from the gauge_force procedure */
  
  register int dir;
  msg_tag *tag0,*tag1,*tag2;
  int start;
  su3_matrix tmat1,tmat2;
  register Real adoma,oma2nt;
  
  start=TRUE; /* indicates staple sum not initialized */
  adoma = alpha_fuzz/(1.-alpha_fuzz);
  oma2nt = pow((double)(1.-alpha_fuzz),(double)nt);

  for(dir=XUP;dir<=TUP;dir++)if(dir != TUP){

    /*         --[0']-->  --[0]-->
              ^          ^        ^
              |          |        |
       ^    [TUP]        |       [2]
       |      |          |        |
       t      O'-[dir']->O-[dir]->
  
               <-"lower"  "upper"->
     
         dir ->          */


    /* get link[dir] from direction TUP: Pulls [0] down to O */
    tag0 = start_gather_site( F_OFFSET(link[dir]), sizeof(su3_matrix),
			TUP, EVENANDODD, gen_pt[0] );
    
    /* Start gather for the "upper staple" pulls [2] over to O */
    tag2 = start_gather_site( F_OFFSET(link[TUP]), sizeof(su3_matrix),
			dir, EVENANDODD, gen_pt[2] );
    
    /* begin the computation "at the dirDOWN point = O'", we will
       later gather the intermediate result "to the home point" */
    
    /* [1] = [dir']^* [TUP] [0'] on each site */
    wait_gather(tag0);
    FORALLSITES(i,st){
      mult_su3_an( &(st->link[dir]), &(st->link[TUP]), &tmat1 );
      mult_su3_nn( &tmat1, (su3_matrix *)gen_pt[0][i], &(st->tempmat1) );
    }
    
    /* Gather this partial result "up to home site" */
    tag1 = start_gather_site( F_OFFSET(tempmat1), sizeof(su3_matrix),
			OPP_DIR(dir), EVENANDODD, gen_pt[1] );
    
    /* begin the computation of the "upper" staple.  Note that
       one of the links has already been gathered, since it
       was used in computing the "lower" staple of the site
       above us (in dir) */
    
    /* staple += [dir] * [2] * [0]^* */
    wait_gather(tag2);
    if(start){	/* this is the first contribution to staple */
      FORALLSITES(i,st){
	mult_su3_nn( &(st->link[dir]), (su3_matrix *)gen_pt[2][i], &tmat1);
	mult_su3_na( &tmat1, (su3_matrix *)gen_pt[0][i], &(st->staple) );
      }
      start=FALSE;
    }
    else{
      FORALLSITES(i,st){
	mult_su3_nn( &(st->link[dir]), (su3_matrix *)gen_pt[2][i], &tmat1);
	mult_su3_na( &tmat1, (su3_matrix *)gen_pt[0][i], &tmat2 );
	add_su3_matrix( &(st->staple),&tmat2,&(st->staple));
      }
    }
    
    /* Add gathered left (down) staple to right (up) staple */
    /* staple += [1] */
    wait_gather(tag1);
    FORALLSITES(i,st){
      add_su3_matrix( &(st->staple),(su3_matrix *)gen_pt[1][i],&(st->staple));
    }
    cleanup_gather(tag0);
    cleanup_gather(tag1);
    cleanup_gather(tag2);
  }

  /* Now we calculate the link [TUP] + alpha_fuzz/(1 - alpha_fuzz)*[staple]/6 
     to be used in the fuzzy Polyakov loop. (The 1/6 averages over
     the 6 staples)
     We use this form, since there is already a canned routine for
     computing such a linear combination.
     Later on we renormalize Re P by multiplying by (1 - alpha_fuzz)^N_t
     Note: The KS phases are included in the staples.  This means that the 
     staples have an extra minus sign relative to the forward link.
     To compensate, we include a minus sign in the linear combination.

     The result is put in "staple"  */

  FORALLSITES(i,st){
    scalar_mult_add_su3_matrix( &(st->link[TUP]), &(st->staple),
			       -adoma/6., &(st->staple) );
  }
  
  /* Now we proceed with the usual ploop routine, using "staple"
     instead of link[TUP] in the product */
  /* Note that we are assuming that the matrix multplication routines
     do not expect that the matrices are su3, but merely 3x3 complex.
     That is the case for the current coding */

  sum = cmplx(0.0,0.0);
  d[XUP] = d[YUP] = d[ZUP] = 0;
  /* First multiply the link on every even site by the link above it */
  /* We will compute the Polyakov loop "at" the even sites in the 
     first two time slices. */
  tag=start_gather_site( F_OFFSET(staple), sizeof(su3_matrix),
		   TUP, EVEN, gen_pt[0] );
  wait_gather(tag);
  FOREVENSITES(i,st){
    mult_su3_nn( &(st->staple), (su3_matrix *)gen_pt[0][i], &(st->tempmat1));
  }
  cleanup_gather(tag);
  
  for(t=2;t<nt;t+=2){
    d[TUP] = t;	/* distance from which to gather */
    tag=start_general_gather_site( F_OFFSET(tempmat1), sizeof(su3_matrix),
			     d, EVEN, gen_pt[0] );
    wait_general_gather(tag);
    FOREVENSITES(i,st){
      if( st->t > 1 )continue;  /* only compute on first two slices */
      mult_su3_nn( &(st->tempmat1), (su3_matrix *)gen_pt[0][i],
	&(st->tempmat2));
      lattice[i].tempmat1 = lattice[i].tempmat2;
      /* We overwrite tempmat1 on the first two time slices,
	 leaving the others undisturbed so we can still gather
	 them. */
    }
    cleanup_general_gather(tag);
  }
  
  FOREVENSITES(i,st){
    if( st->t > 1 )continue;
    plp_fuzz = trace_su3( &(st->tempmat1) );
    CMULREAL(plp_fuzz,oma2nt,plp_fuzz);
    CSUM(sum,plp_fuzz);
    /* Save for subsequent correlation measurements */
    /* Note the results are saved on even sites in
       slices 0 and 1 */
    st->ploop_fuzz = plp_fuzz;
  }
  g_complexsum( &sum );
  plp_staple.real = sum.real /((Real)(nx*ny*nz));
  plp_staple.imag = sum.imag /((Real)(nx*ny*nz));
  return(plp_staple);
}
