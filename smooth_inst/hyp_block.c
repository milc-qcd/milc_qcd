/******************** hyp_block.c *******************************/
/* Deterministic hypercube blocking ( formerly ape_block_e2.c) */
/* A. Hasenfratz */
/* commented by F. Knechtli - February 2001 */

/* Version 7 port
   Modified to do only one sweep and replace original links
   with result, followed by reunitarization
   Changed entry name from "ape_block_det" to 
   "smooth" C. DeTar - August 2001 */
/*
   the fat links are constructed in three steps.
   1)
   wiggly link(count) = P ( (1-alpha3)*link(dir1) + 
                            (alpha3/2)*staples(dir1,dir2) )
   dir1  = direction of the wiggly link
   dir2  = orthogonal direction of decoration (2 staples)
   count = 3*dir1+dir2; if(dir2>dir1)count=count-1 (0<=count<12)

   the wiggly links are stored in blocked_link[20+count]

   2)
   double link(count) = P ( (1-alpha2)*link(dir1) +
                            (alpha2/4)*wiggly_staples(dir1,NOT dir2) )
   dir1  = direction of the double link
   dir2  = orthogonal directions of NO decoration (4 staples)
   count = 3*dir1+dir2; if(dir2>dir1)count=count-1 (0<=count<12)

   the double links are stored in blocked_link[8+count]

   3)
   fat link(dir)      = P ( (1-alpha)*link(dir) +
                            (alpha/6)*double_staples(all 6) )

   the fat links are stored in blocked_link[0..3]
   a copy of the original thin links is in blocked_link[4..7]
*/

#include "smooth_inst_includes.h"

#define Nc 3

/* deterministic APE blocking */
void smooth()
{
  int parity;
  register int dir,i,dir2;
  register site *st;
  su3_matrix fatq;
  
  void dsdu_ape2(register int dir, int parity);
  void dsdu_ape_ext2(register int dir1, register int dir2, int parity);
  void dsdu_ape_ext1(register int dir1, register int dir2, int parity);
  
  
  FORALLSITES(i,st)for(dir=XUP;dir<=TUP;dir++){
    st->blocked_link[4+dir]= st->link[dir];
  }
  
  /* precompute the wiggly links and store them in blocked_link[20+count] */
  for(dir=XUP;dir<=TUP;dir++)
    for(dir2=XUP;dir2<=TUP;dir2++)if(dir2 != dir){
      dsdu_ape_ext2(dir,dir2,1);
      dsdu_ape_ext2(dir,dir2,2);
    }
  
  /* precompute the double links and store them in blocked_link[8+count] */
  for(dir=XUP;dir<=TUP;dir++)
    for(dir2=XUP;dir2<=TUP;dir2++)if(dir2 != dir){
      dsdu_ape_ext1(dir,dir2,1);
      dsdu_ape_ext1(dir,dir2,2);
    }
  
  for(parity=ODD;parity<=EVEN;parity++){
    for(dir=XUP;dir<=TUP;dir++){
      
      /* compute the decorated staple in tempmat1 */
      dsdu_ape2(dir,parity); 
      
      FORSOMEPARITY(i,st,parity){
	
	st->blocked_link[dir]=st->blocked_link[4+dir];
	
	/* compute the fat link : fatq */
        scalar_mult_su3_matrix( &(st->blocked_link[4+dir]), 
				(1.0-alpha), &fatq );
        scalar_mult_add_su3_matrix( &fatq, 
				    &(st->tempmat1), alpha/6.0, &fatq );
	
	/* Project fatq onto SU(3) - result in blocked_link[dir] */
	project_su3( &(st->blocked_link[dir]), &fatq, hits, 0.);
      } /* site */
    }} /*  direction and parity */
  
  /* Replace original links with result */
  
  FORALLSITES(i,st)for(dir=XUP;dir<=TUP;dir++){
    st->link[dir]= st->blocked_link[dir];
  }
  
  /* reunitarize the gauge field */
  reunitarize();
  
} /* smooth */


/* dsdu_ape2 -- compute the decorated staple, returning the staple in tempmat1 */

void dsdu_ape2(register int dir1, int parity) 
{
  register int i,dir2;
  register site *st;
  msg_tag *tag0,*tag1,*tag2,*tag3,*tag4;
  int start;
  register int count1,count2;
  su3_matrix tmat1,tmat2;
  int disp[4];	/* displacement vector for general gather */
  
  start=1; /* indicates staple sum not initialized */
  /* Loop over other directions, computing force from plaquettes in
     the dir1,dir2 plane */
  for(dir2=XUP;dir2<=TUP;dir2++)if(dir2 != dir1)
    {
      /* displacement vector for bb_link 2 sites away */
      for(i=XUP;i<=TUP;i++)disp[i]=0;
      disp[dir1] = 1;
      disp[dir2] = -1;
      
      /* counters for the double links to be used */
      /* the double links are NOT decorated in the directions dir1 & dir2 */
      count1=3*dir1+dir2;
      if(dir2>dir1)count1=count1-1;
      count2=3*dir2+dir1;
      if(dir1>dir2)count2=count2-1;
      
      /* get blocked_link[8+count2] from direction dir1 */
      tag0 = start_gather_site( F_OFFSET(blocked_link[8+count2]), sizeof(su3_matrix),
			   dir1, parity, gen_pt[0] );
      
      /* get blocked_link[8+count1] from direction dir2 */
      tag1 = start_gather_site( F_OFFSET(blocked_link[8+count1]), sizeof(su3_matrix),
			   dir2, parity, gen_pt[1] );
      
      /* get blocked_link[8+count2] from direction -dir2 */
      tag2 = start_gather_site( F_OFFSET(blocked_link[8+count2]), sizeof(su3_matrix),
			   OPP_DIR(dir2), parity, gen_pt[2] );
      
      /* get blocked_link[8+count1] from direction -dir2 */
      tag3 = start_gather_site( F_OFFSET(blocked_link[8+count1]), sizeof(su3_matrix),
			   OPP_DIR(dir2), parity, gen_pt[3] );
      
      /* get blocked_link[8+count2] from displacement +dir1-dir2 */
      tag4 = start_general_gather_site( F_OFFSET(blocked_link[8+count2]),
				   sizeof(su3_matrix), disp, parity, gen_pt[4] );
      
      /* Upper staple */
      wait_gather(tag0);
      wait_gather(tag1);
      if(start){  /* this is the first contribution to staple */
	FORSOMEPARITY(i,st,parity){
	  mult_su3_nn( &(st->blocked_link[8+count2]), (su3_matrix *)gen_pt[1][i], &tmat1 );
	  mult_su3_na( &tmat1, (su3_matrix *)gen_pt[0][i], &(st->tempmat1) );
	  
	}
	start=0; 
      }
      else{
	FORSOMEPARITY(i,st,parity){
	  mult_su3_nn( &(st->blocked_link[8+count2]), (su3_matrix *)gen_pt[1][i], &tmat1 );
	  mult_su3_na( &tmat1, (su3_matrix *)gen_pt[0][i], &tmat2 );
	  add_su3_matrix( &(st->tempmat1), &tmat2, &(st->tempmat1));
	}
      } /* upper tempmat1 */
      cleanup_gather(tag0);
      cleanup_gather(tag1);
      
      /* Lower tempmat1 */
      wait_gather(tag2);
      wait_gather(tag3);
      wait_general_gather(tag4);
      FORSOMEPARITY(i,st,parity){
	mult_su3_an( (su3_matrix *)gen_pt[2][i], 
		     (su3_matrix *)gen_pt[3][i], &tmat1 );
	mult_su3_nn( &tmat1, (su3_matrix *)gen_pt[4][i], &tmat2 );
	add_su3_matrix( &(st->tempmat1), &tmat2, &(st->tempmat1));
      }  /* lower tempmat1 */
      cleanup_gather(tag2);
      cleanup_gather(tag3);
      cleanup_general_gather(tag4);
    }
}


/* compute the double links in direction dir1 NOT decorated in
   direction dir3 */

void dsdu_ape_ext1(register int dir1, register int dir3, int parity) 
{
  /* dir3 is the direction that's excluded */
  register int i,dir2,dir4 = -1;
  register site *st;
  msg_tag *tag0,*tag1,*tag2,*tag3,*tag4;
  int start;
  register int count;
  register int count1,count2;
  su3_matrix tmat1,tmat2,fatq;
  int disp[4];	/* displacement vector for general gather */
  
  
  start=0; /* indicates staple sum not initialized */
  /* double links are stored in blocked_link[8+count] */
  count=3*dir1+dir3;
  if(dir3>dir1)count=count-1;
  
  /* Loop over other directions, computing force from plaquettes in
     the dir1,dir2 plane */
  for(dir2=XUP;dir2<=TUP;dir2++)if(dir2 != dir1 && dir2 != dir3){
    /* displacement vector for bb_link 2 sites away */
    for(i=XUP;i<=TUP;i++)disp[i]=0;
    disp[dir1] = 1;
    disp[dir2] = -1;
    
    /* dir4 = direction of decoration of the wiggly links */
    for(i=XUP;i<=TUP;i++)if(i != dir1 && i != dir2 && i != dir3)dir4=i;
    count1=3*dir1+dir4;
    if(dir4>dir1)count1=count1-1;
    count2=3*dir2+dir4;
    if(dir4>dir2)count2=count2-1;
    
    /* get blocked_link[20+count2] from direction dir1 */
    tag0 = start_gather_site( F_OFFSET(blocked_link[20+count2]), sizeof(su3_matrix),
			 dir1, parity, gen_pt[0] );
    
    /* get blocked_link[20+count1] from direction dir2 */
    tag1 = start_gather_site( F_OFFSET(blocked_link[20+count1]), sizeof(su3_matrix),
			 dir2, parity, gen_pt[1] );
    
    /* get blocked_link[20+count2] from direction -dir2 */
    tag2 = start_gather_site( F_OFFSET(blocked_link[20+count2]), sizeof(su3_matrix),
			 OPP_DIR(dir2), parity, gen_pt[2] );
    
    /* get blocked_link[20+count1] from direction -dir2 */
    tag3 = start_gather_site( F_OFFSET(blocked_link[20+count1]), sizeof(su3_matrix),
			 OPP_DIR(dir2), parity, gen_pt[3] );
    
    /* get blocked_link[20+count2] from displacement +dir1-dir2 */
    tag4 = start_general_gather_site( F_OFFSET(blocked_link[20+count2]),
				 sizeof(su3_matrix), disp, parity, gen_pt[4] );
    
    /* Upper staple */
    wait_gather(tag0);
    wait_gather(tag1);
    if(start == 0){  /* this is the first contribution to staple */
      FORSOMEPARITY(i,st,parity){
	mult_su3_nn( &(st->blocked_link[20+count2]), (su3_matrix *)gen_pt[1][i], &tmat1 );
	mult_su3_na( &tmat1, (su3_matrix *)gen_pt[0][i], &(st->tempmat1) );
      }
      start++; 
    }
    else{
      FORSOMEPARITY(i,st,parity){
	mult_su3_nn( &(st->blocked_link[20+count2]), (su3_matrix *)gen_pt[1][i], &tmat1 );
	mult_su3_na( &tmat1, (su3_matrix *)gen_pt[0][i], &tmat2 );
	add_su3_matrix( &(st->tempmat1), &tmat2, &(st->tempmat1));
      }
    } /* upper tempmat1 */
    cleanup_gather(tag0);
    cleanup_gather(tag1);
    
    /* Lower tempmat1 */
    wait_gather(tag2);
    wait_gather(tag3);
    wait_general_gather(tag4);
    FORSOMEPARITY(i,st,parity){
      mult_su3_an( (su3_matrix *)gen_pt[2][i], (su3_matrix *)gen_pt[3][i], &tmat1 );
      mult_su3_nn( &tmat1, (su3_matrix *)gen_pt[4][i], &tmat2 );
      add_su3_matrix( &(st->tempmat1), &tmat2, &(st->tempmat1));
    }  /* lower tempmat1 */
    cleanup_gather(tag2);
    cleanup_gather(tag3);
    cleanup_general_gather(tag4);
  }
  FORSOMEPARITY(i,st,parity){
    scalar_mult_su3_matrix( &(st->tempmat1), 
			    (alpha2/4.0), &(st->tempmat1));
    scalar_mult_add_su3_matrix( &(st->tempmat1), 
				&(st->blocked_link[4+dir1]),(1.0-alpha2), &(fatq));
    (st->blocked_link[8+count])= (st->blocked_link[4+dir1]);
    
    /* Project fatq onto SU(3) - result in blocked_link[8+count] */
    project_su3( &(st->blocked_link[8+count]), &fatq, hits, 0.);
    
  }
} /* dsdu_ape_ext1 */


/* compute wiggly links in direction dir1 decorated in direction dir2 */

void dsdu_ape_ext2(register int dir1, register int dir2, int parity) 
{
  register int i;
  register site *st;
  msg_tag *tag0,*tag1,*tag2,*tag3,*tag4;
  int start;
  register int count;
  su3_matrix tmat1,tmat2,fatq;
  int disp[4];	/* displacement vector for general gather */
  
  
  start=1; /* indicates staple sum not initialized */
  /* wiggly links are stored in blocked_link[20+count] */
  count=3*dir1+dir2;
  if(dir2>dir1)count=count-1; 
  /* displacement vector for bb_link 2 sites away */
  for(i=XUP;i<=TUP;i++)disp[i]=0;
  disp[dir1] = 1;
  disp[dir2] = -1;
  
  /* get blocked_link[4+dir2] from direction dir1 */
  tag0 = start_gather_site( F_OFFSET(blocked_link[4+dir2]), sizeof(su3_matrix),
		       dir1, parity, gen_pt[0] );
  
  /* get blocked_link[4+dir1] from direction dir2 */
  tag1 = start_gather_site( F_OFFSET(blocked_link[4+dir1]), sizeof(su3_matrix),
		       dir2, parity, gen_pt[1] );
  
  /* get blocked_link[4+dir2] from direction -dir2 */
  tag2 = start_gather_site( F_OFFSET(blocked_link[4+dir2]), sizeof(su3_matrix),
		       OPP_DIR(dir2), parity, gen_pt[2] );
  
  /* get blocked_link[4+dir1] from direction -dir2 */
  tag3 = start_gather_site( F_OFFSET(blocked_link[4+dir1]), sizeof(su3_matrix),
		       OPP_DIR(dir2), parity, gen_pt[3] );
  
  /* get blocked_link[4+dir2] from displacement +dir1-dir2 */
  tag4 = start_general_gather_site( F_OFFSET(blocked_link[4+dir2]),
			       sizeof(su3_matrix), disp, parity, gen_pt[4] );
  
  /* Upper staple */
  wait_gather(tag0);
  wait_gather(tag1);
  if(start){  /* this is the first contribution to staple */
    FORSOMEPARITY(i,st,parity){
      mult_su3_nn( &(st->blocked_link[4+dir2]), (su3_matrix *)gen_pt[1][i], &tmat1 );
      mult_su3_na( &tmat1, (su3_matrix *)gen_pt[0][i], &(st->tempmat1) );
      
    }
    start=0; 
  }
  else{
    FORSOMEPARITY(i,st,parity){
      mult_su3_nn( &(st->blocked_link[4+dir2]), (su3_matrix *)gen_pt[1][i], &tmat1 );
      mult_su3_na( &tmat1, (su3_matrix *)gen_pt[0][i], &tmat2 );
      add_su3_matrix( &(st->tempmat1), &tmat2, &(st->tempmat1));
    }
  } /* upper tempmat1 */
  cleanup_gather(tag0);
  cleanup_gather(tag1);
  
  /* Lower tempmat1 */
  wait_gather(tag2);
  wait_gather(tag3);
  wait_general_gather(tag4);
  FORSOMEPARITY(i,st,parity){
    mult_su3_an( (su3_matrix *)gen_pt[2][i], (su3_matrix *)gen_pt[3][i], &tmat1 );
    mult_su3_nn( &tmat1, (su3_matrix *)gen_pt[4][i], &tmat2 );
    add_su3_matrix( &(st->tempmat1), &tmat2, &(st->tempmat1));
  }  /* lower tempmat1 */
  cleanup_gather(tag2);
  cleanup_gather(tag3);
  cleanup_general_gather(tag4);
  FORSOMEPARITY(i,st,parity){
    scalar_mult_su3_matrix( &(st->tempmat1), 
			    (alpha3/2.0), &(st->tempmat1));
    scalar_mult_add_su3_matrix( &(st->tempmat1), 
				&(st->blocked_link[4+dir1]),(1.0-alpha3), &(fatq));
    (st->blocked_link[20+count])= (st->blocked_link[4+dir1]);
    
    /* Project fatq onto SU(3) - result in blocked_link[20+count] */
    project_su3( &(st->blocked_link[20+count]), &fatq, hits, 0.);
  }
} /* dsdu_ape_ext2 */
