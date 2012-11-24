/* Routines for field-based HYP-smearing,
   based on HYP-smearing used in smooth_inst application
   A. Bazavov, Oct 2012 */
/* Deterministic hypercube blocking ( formerly ape_block_e2.c) */
/* A. Hasenfratz */
/* commented by F. Knechtli - February 2001 */

/* (from smooth_inst/hyp_block.c):

   the fat links are constructed in three steps.
   1)
   wiggly link(count) = P ( (1-alpha3)*link(dir1) + 
                            (alpha3/2)*staples(dir1,dir2) )
   dir1  = direction of the wiggly link
   dir2  = orthogonal direction of decoration (2 staples)
   count = 3*dir1+dir2; if(dir2>dir1)count=count-1 (0<=count<12)

   the wiggly links are stored in blocked_link[20+count]
   (Wiggly_link array in this code)

   2)
   double link(count) = P ( (1-alpha2)*link(dir1) +
                            (alpha2/4)*wiggly_staples(dir1,NOT dir2) )
   dir1  = direction of the double link
   dir2  = orthogonal directions of NO decoration (4 staples)
   count = 3*dir1+dir2; if(dir2>dir1)count=count-1 (0<=count<12)

   the double links are stored in blocked_link[8+count]
   (Doubly_link array in this code)

   3)
   fat link(dir)      = P ( (1-alpha)*link(dir) +
                            (alpha/6)*double_staples(all 6) )

   the fat links are stored in blocked_link[0..3]
   a copy of the original thin links is in blocked_link[4..7]
   (HYP links are returned in hyp_link array in this code)

   AB: this is modified to handle 3D smearing (spatial but could be any)
       where step 3) is not needed
*/


#include <stdlib.h>
#include <stdio.h>
#include "../include/complex.h"
#include "../include/su3.h"
#include "../include/macros.h"
#include "../include/comdefs.h"
#include "../include/generic.h"
#include <lattice.h>


//AB TEMPORARY HEADER
#include "../include/hyp_coeff.h"


#ifdef QCDOC
#define special_alloc qcdoc_alloc
#define special_free qfree
#else
#define special_alloc malloc
#define special_free free
#endif


// create a field with "n" entities on a site
// (n=4 is a usual link field)
static su3_matrix *create_mn_special(int n){
  char myname[] = "create_mn_special";
  su3_matrix *m;

  m = (su3_matrix *)special_alloc(sites_on_node*n*sizeof(su3_matrix));

  if(m==NULL){
    printf("%s: no room\n",myname);
    terminate(1);
  }

  return m;
}

static void destroy_mn_special(su3_matrix *m){
  special_free(m);
}

/*--------------------------------------------------------------------*/
#if 0
// higher precison su3_matrix dump

static void dumpmat_hp( su3_matrix *m ){
int i,j;
    for(i=0;i<3;i++){
	for(j=0;j<3;j++)printf("(%.15f,%.15f)\t",
	    m->e[i][j].real,m->e[i][j].imag);
	printf("\n");
    }
    printf("\n");
}
#endif

//AB external declaration, NEEDS PUBLIC API(!)
//msg_tag *
//start_general_strided_gather(
//  void *field,          /* source buffer aligned to desired field */
//  int stride,           /* bytes between fields in source buffer */
//  int size,             /* size in bytes of the field (eg sizeof(su3_vector))*/
//  int *displacement,    /* displacement to gather from. four components */
//  int parity,           /* parity of sites to which we gather.
//                           one of EVEN, ODD or EVENANDODD. */
//  char ** dest)         /* one of the vectors of pointers */ ;



/* compute wiggly links in direction dir1 decorated in direction dir2 */
void hyp_block_stage1(register int dir1, register int dir2, int parity,
  su3_matrix *U_link, su3_matrix *Wiggly_link, int dir_exclude,
  hyp_coeffs_t *hc ) {

  register int i;
  register site *st;
  msg_tag *tag0,*tag1,*tag2,*tag3,*tag4;
  int start;
  register int count,nWiggly;
  su3_matrix tmat1,tmat2,fatq;
  su3_matrix *tempmat1;
  int disp[4];	/* displacement vector for general gather */

  /* create temporary storage, one matrix per site */
  tempmat1 = create_mn_special(1);

  start=1; /* indicates staple sum not initialized */

  // array size is fixed for 4D, in 3D not all entries are filled
  nWiggly = 12;

  count=3*dir1+dir2;
  if(dir2>dir1)count=count-1; 

  /* displacement vector for link 2 sites away */
  for(i=XUP;i<=TUP;i++)disp[i]=0;
  disp[dir1] = 1;
  disp[dir2] = -1;
  
  /* get U_link[dir2] from direction dir1 */
  tag0 = declare_strided_gather( U_link + dir2, 4*sizeof(su3_matrix),
           sizeof(su3_matrix), dir1, parity, gen_pt[0] );
  do_gather( tag0 );
  
  /* get U_link[dir1] from direction dir2 */
  tag1 = declare_strided_gather( U_link + dir1, 4*sizeof(su3_matrix),
           sizeof(su3_matrix), dir2, parity, gen_pt[1] );
  do_gather( tag1 );
  
  /* get U_link[dir2] from direction -dir2 */
  tag2 = declare_strided_gather( U_link + dir2, 4*sizeof(su3_matrix),
           sizeof(su3_matrix), OPP_DIR(dir2), parity, gen_pt[2] );
  do_gather( tag2 );
  
  /* get U_link[dir1] from direction -dir2 */
  tag3 = declare_strided_gather( U_link + dir1, 4*sizeof(su3_matrix),
           sizeof(su3_matrix), OPP_DIR(dir2), parity, gen_pt[3] );
  do_gather( tag3 );
  
  /* get U_link[dir2] from displacement +dir1-dir2 */
  tag4 = start_general_strided_gather( (char *)(U_link + dir2), 4*sizeof(su3_matrix),
           sizeof(su3_matrix), disp, parity, gen_pt[4] );

  /* Upper staple */
  wait_gather(tag0);
  wait_gather(tag1);
  if(start){  /* this is the first contribution to staple */
    FORSOMEPARITY(i,st,parity){
      mult_su3_nn( &(U_link[4*i+dir2]), (su3_matrix *)gen_pt[1][i], &tmat1 );
      mult_su3_na( &tmat1, (su3_matrix *)gen_pt[0][i], &(tempmat1[i]) );
      
    }
    start=0; 
  }
  else{
    FORSOMEPARITY(i,st,parity){
      mult_su3_nn( &(U_link[4*i+dir2]), (su3_matrix *)gen_pt[1][i], &tmat1 );
      mult_su3_na( &tmat1, (su3_matrix *)gen_pt[0][i], &tmat2 );
      add_su3_matrix( &(tempmat1[i]), &tmat2, &(tempmat1[i]) );
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
    add_su3_matrix( &(tempmat1[i]), &tmat2, &(tempmat1[i]) );
  }  /* lower tempmat1 */
  cleanup_gather(tag2);
  cleanup_gather(tag3);
  cleanup_general_gather(tag4);
  FORSOMEPARITY(i,st,parity){
    scalar_mult_su3_matrix( &(tempmat1[i]), 
			    (hc->alpha3/2.0), &(tempmat1[i]) );
    scalar_mult_add_su3_matrix( &(tempmat1[i]), 
				&(U_link[4*i+dir1]),(1.0-hc->alpha3), &(fatq));

    switch( hc->proj_method ) {
      case HYP_U3:
        printf( "HYP_U3 projection: NOT READY\n" );
        terminate(1);
        break;
      case HYP_SU3:
        printf( "HYP_SU3 projection: NOT READY\n" );
        terminate(1);
        break;
      case HYP_SU3_TR_MAX:
        /* take original link as guess for projected */
        (Wiggly_link[nWiggly*i+count])=(U_link[4*i+dir1]);
        /* Project fatq onto SU(3) - result in Wiggly_link[count] */
        project_su3( &(Wiggly_link[nWiggly*i+count]), &fatq, hc->hits, 0.);
//AB DEBUGGING
//if(i==0 && dir1==1 && dir2==2 ) {
//  dumpmat_hp( &(Wiggly_link[nWiggly*0+count]) );
//  dumpmat_hp( &fatq );
//}
        break;
      default:
        printf( "Wrong HYP projection method\n" );
        terminate(1);
    }
  }


  /* free temporary storage */
  destroy_mn_special(tempmat1);

} /* hyp_block_stage1 */


/* compute the double links in direction dir1 NOT decorated in
   direction dir3 */
void hyp_block_stage2(register int dir1, register int dir3, int parity,
  su3_matrix *U_link, su3_matrix *Wiggly_link, su3_matrix *Doubly_link,
  hyp_coeffs_t *hc ) {

  /* dir3 is the direction that's excluded by HYP construction,
     dir_exclude is the direction excluded if we want 3D instead of 4D */
  register int i,dir2,dir4 = -1;
  register site *st;
  msg_tag *tag0,*tag1,*tag2,*tag3,*tag4;
  int start;
  register int count;
  register int count1,count2;
  register int nWiggly,nDoubly;
  su3_matrix tmat1,tmat2,fatq;
  su3_matrix *tempmat1;
  int disp[4];	/* displacement vector for general gather */
  
  
  /* create temporary storage, one matrix per site */
  tempmat1 = create_mn_special(1);

  start=0; /* indicates staple sum not initialized */

  // array sizes fixed for 4D, in 3D not all entries are filled
  nWiggly = 12;
  nDoubly = 12;

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

    /* get Wiggly_link[count2] from direction dir1 */
    tag0 = declare_strided_gather( Wiggly_link + count2, nWiggly*sizeof(su3_matrix),
             sizeof(su3_matrix), dir1, parity, gen_pt[0] );
    do_gather( tag0 );

    /* get Wiggly_link[count1] from direction dir2 */
    tag1 = declare_strided_gather( Wiggly_link + count1, nWiggly*sizeof(su3_matrix),
             sizeof(su3_matrix), dir2, parity, gen_pt[1] );
    do_gather( tag1 );

    /* get Wiggly_link[count2] from direction -dir2 */
    tag2 = declare_strided_gather( Wiggly_link + count2, nWiggly*sizeof(su3_matrix),
             sizeof(su3_matrix), OPP_DIR(dir2), parity, gen_pt[2] );
    do_gather( tag2 );

    /* get Wiggly_link[count1] from direction -dir2 */
    tag3 = declare_strided_gather( Wiggly_link + count1, nWiggly*sizeof(su3_matrix),
             sizeof(su3_matrix), OPP_DIR(dir2), parity, gen_pt[3] );
    do_gather( tag3 );
    
    /* get Wiggly_link[count2] from displacement +dir1-dir2 */
    tag4 = start_general_strided_gather( (char *)(Wiggly_link + count2), 
					 nWiggly*sizeof(su3_matrix),
					 sizeof(su3_matrix), disp, parity, gen_pt[4] );

    /* Upper staple */
    wait_gather(tag0);
    wait_gather(tag1);
    if(start == 0){  /* this is the first contribution to staple */
      FORSOMEPARITY(i,st,parity){
	mult_su3_nn( &(Wiggly_link[nWiggly*i+count2]), (su3_matrix *)gen_pt[1][i], &tmat1 );
	mult_su3_na( &tmat1, (su3_matrix *)gen_pt[0][i], &(tempmat1[i]) );
      }
      start++; 
    }
    else{
      FORSOMEPARITY(i,st,parity){
	mult_su3_nn( &(Wiggly_link[nWiggly*i+count2]), (su3_matrix *)gen_pt[1][i], &tmat1 );
	mult_su3_na( &tmat1, (su3_matrix *)gen_pt[0][i], &tmat2 );
	add_su3_matrix( &(tempmat1[i]), &tmat2, &(tempmat1[i]));
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
      add_su3_matrix( &(tempmat1[i]), &tmat2, &(tempmat1[i]));
    }  /* lower tempmat1 */
    cleanup_gather(tag2);
    cleanup_gather(tag3);
    cleanup_general_gather(tag4);
  }
  FORSOMEPARITY(i,st,parity){
    scalar_mult_su3_matrix( &(tempmat1[i]),
			    (hc->alpha2/4.0), &(tempmat1[i]));
    scalar_mult_add_su3_matrix( &(tempmat1[i]),
				&(U_link[4*i+dir1]),(1.0-hc->alpha2), &(fatq));

    switch( hc->proj_method ) {
      case HYP_U3:
        printf( "HYP_U3 projection: NOT READY\n" );
        terminate(1);
        break;
      case HYP_SU3:
        printf( "HYP_SU3 projection: NOT READY\n" );
        terminate(1);
        break;
      case HYP_SU3_TR_MAX:
        /* take original link as guess for projected */
        Doubly_link[nDoubly*i+count]=U_link[4*i+dir1];
        /* Project fatq onto SU(3) - result in blocked_link[8+count] */
        project_su3( &(Doubly_link[nDoubly*i+count]), &fatq, hc->hits, 0.);
//AB DEBUGGING
//if(i==0 && dir1==0 && dir3==3 ) {
//  dumpmat_hp( &(Doubly_link[nDoubly*i+count]) );
//  dumpmat_hp( &fatq );
//}
        break;
      default:
        printf( "Wrong HYP projection method\n" );
        terminate(1);
    }
  }

  /* free temporary storage */
  destroy_mn_special(tempmat1);

} /* hyp_block_stage2 */


/* compute the decorated staple */
void hyp_block_doubly_staple(register int dir1, int parity,
  su3_matrix *Doubly_link, su3_matrix *staple, hyp_coeffs_t *hc ) {
  register int i,dir2;
  register site *st;
  msg_tag *tag0,*tag1,*tag2,*tag3,*tag4;
  int start;
  int nDoubly;
  register int count1,count2;
  su3_matrix tmat1,tmat2;
  int disp[4];	/* displacement vector for general gather */

  nDoubly = 12; // size fixed for 4D, this is not called in 3D

  start=1; /* indicates staple sum not initialized */
  /* Loop over other directions, computing force from plaquettes in
     the dir1,dir2 plane */
  for(dir2=XUP;dir2<=TUP;dir2++)if(dir2 != dir1)
    {
      /* displacement vector for link 2 sites away */
      for(i=XUP;i<=TUP;i++)disp[i]=0;
      disp[dir1] = 1;
      disp[dir2] = -1;
      
      /* counters for the double links to be used */
      /* the double links are NOT decorated in the directions dir1 & dir2 */
      count1=3*dir1+dir2;
      if(dir2>dir1)count1=count1-1;
      count2=3*dir2+dir1;
      if(dir1>dir2)count2=count2-1;

      /* get Doubly_link[count2] from direction dir1 */
      tag0 = declare_strided_gather( Doubly_link + count2, nDoubly*sizeof(su3_matrix),
               sizeof(su3_matrix), dir1, parity, gen_pt[0] );
      do_gather( tag0 );
      
      /* get Doubly_link[count1] from direction dir2 */
      tag1 = declare_strided_gather( Doubly_link + count1, nDoubly*sizeof(su3_matrix),
               sizeof(su3_matrix), dir2, parity, gen_pt[1] );
      do_gather( tag1 );
      
      /* get Doubly_link[count2] from direction -dir2 */
      tag2 = declare_strided_gather( Doubly_link + count2, nDoubly*sizeof(su3_matrix),
               sizeof(su3_matrix), OPP_DIR(dir2), parity, gen_pt[2] );
      do_gather( tag2 );
      
      /* get Doubly_link[count1] from direction -dir2 */
      tag3 = declare_strided_gather( Doubly_link + count1, nDoubly*sizeof(su3_matrix),
               sizeof(su3_matrix), OPP_DIR(dir2), parity, gen_pt[3] );
      do_gather( tag3 );
      
      /* get Doubly_link[count2] from displacement +dir1-dir2 */
      tag4 = start_general_strided_gather( (char *)(Doubly_link + count2), 
					   nDoubly*sizeof(su3_matrix),
					   sizeof(su3_matrix), disp, parity, gen_pt[4] );

     /* Upper staple */
      wait_gather(tag0);
      wait_gather(tag1);
      if(start){  /* this is the first contribution to staple */
        FORSOMEPARITY(i,st,parity){
          mult_su3_nn( &(Doubly_link[nDoubly*i+count2]), (su3_matrix *)gen_pt[1][i], &tmat1 );
          mult_su3_na( &tmat1, (su3_matrix *)gen_pt[0][i], &(staple[i]) );
        }
        start=0; 
      }
      else{
        FORSOMEPARITY(i,st,parity){
          mult_su3_nn( &(Doubly_link[nDoubly*i+count2]), (su3_matrix *)gen_pt[1][i], &tmat1 );
          mult_su3_na( &tmat1, (su3_matrix *)gen_pt[0][i], &tmat2 );
          add_su3_matrix( &(staple[i]), &tmat2, &(staple[i]));
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
        add_su3_matrix( &(staple[i]), &tmat2, &(staple[i]));
      }  /* lower tempmat1 */
      cleanup_gather(tag2);
      cleanup_gather(tag3);
      cleanup_general_gather(tag4);
    }
}


/* HYP smearing
   if dir_exclude==NODIR -- usual 4D HYP smearing
   if dir_exclude==XUP or YUP or ZUP or TUP, 3D smearing
   TODO: VARIOUS PROJECTION METHODS */
void load_hyp_links(su3_matrix *U_link, su3_matrix *hyp_link,
  int dir_exclude, hyp_coeffs_t *hc) {
  int parity;
  register int dir,i,dir2;
  register site *st;
  su3_matrix fatq;
  su3_matrix *Wiggly_link, *Doubly_link, *staple;
  // int nWiggly;
  int nDoubly, count;


  /* dimensions of the arrays are preset for 4D case, where we have
     12 wiggly and 12 doubly links on each site;
     in 3D we have only 6 wiggly and 3 doubly links per site;
     array sizes are fixed but in 3D case not all entries are filled */
  //  nWiggly = 12;
  nDoubly = 12;

  // create temporary storage
  Wiggly_link = create_mn_special(12);
  Doubly_link = create_mn_special(12);
  staple = create_mn_special(1);
  
  /* precompute the wiggly links and store them in Wiggly_link */
  for( dir=XUP; dir<=TUP; dir++ ) {
    if( dir!=dir_exclude ) { // excluded direction for 3D case
      for( dir2=XUP; dir2<=TUP; dir2++ ) {
        if( dir2!=dir && dir2!=dir_exclude ) { // second clause for 3D case
          hyp_block_stage1( dir, dir2, 1, U_link, Wiggly_link, dir_exclude, hc );
          hyp_block_stage1( dir, dir2, 2, U_link, Wiggly_link, dir_exclude, hc );
        }
      }
    }
  }

  /* precompute the double links and store them in Doubly_link */
  for( dir=XUP; dir<=TUP; dir++ ) {
    if( dir!=dir_exclude ) { // excluded direction for 3D case
      for( dir2=XUP; dir2<=TUP; dir2++ ) {
        if( dir2!=dir ) {
          if( NODIR==dir_exclude || 
            ( NODIR!=dir_exclude && dir2==dir_exclude ) ) { // for 3D case,
            // note "==" in second clause: double links are encoded by the second
            // direction where they are NOT decorated, thus in 3D we want those
            // that are decorated in all directions except of dir_exclude,
            // if dir_exclude is not NODIR
            hyp_block_stage2( dir, dir2, 1, U_link, Wiggly_link, Doubly_link, hc );
            hyp_block_stage2( dir, dir2, 2, U_link, Wiggly_link, Doubly_link, hc );
          }
        }
      }
    }
  }


  /* main branching between 4D and 3D case:
     in 3D we are done at this point since Doubly_link contains three links
     smeared in 3D cubes; in 4D we need an extra smearing */
  if( NODIR!=dir_exclude ) { // 3D HYP smearing
    // load Doubly_link into hyp_link for three directions
    // and U_link for the excluded direction (i.e. unsmeared)
    FORALLSITES( i,st ) {
      for( dir=XUP; dir<=TUP; dir++ ) {
        if( dir!=dir_exclude ) { // excluded direction for 3D case
          dir2=dir_exclude;
          // counting scheme needs to be the same as in Doubly_link routine(!)
          count = 3*dir + dir2;
          if( dir2>dir ) count = count-1;
          hyp_link[4*i+dir] = Doubly_link[nDoubly*i+count];
        }
        else { // not smeared in this direction
          hyp_link[4*i+dir] = U_link[4*i+dir];
        }
      }
    }
  }
  else { // 4D HYP smearing
    for( parity=ODD; parity<=EVEN; parity++ ) {
      for( dir=XUP; dir<=TUP; dir++ ) {

        /* compute the decorated staple in "staple" */
        hyp_block_doubly_staple(dir, parity, Doubly_link, staple, hc );
      
        FORSOMEPARITY( i,st,parity ) {


          /* compute the fat link : fatq */
          scalar_mult_su3_matrix( &(U_link[4*i+dir]), 
                                   (1.0-hc->alpha1), &fatq );
          scalar_mult_add_su3_matrix( &fatq, 
                                      &(staple[i]), hc->alpha1/6.0, &fatq );

          switch( hc->proj_method ) {
            case HYP_U3:
              printf( "HYP_U3 projection: NOT READY\n" );
              terminate(1);
              break;
            case HYP_SU3:
              printf( "HYP_SU3 projection: NOT READY\n" );
              terminate(1);
              break;
            case HYP_SU3_TR_MAX:
              /* take original link as guess for projected */
              hyp_link[4*i+dir] = U_link[4*i+dir];
              /* Project fatq onto SU(3) - result in hyp_link[dir] */
              project_su3( &(hyp_link[4*i+dir]), &fatq, hc->hits, 0.);
              break;
            default:
              printf( "Wrong HYP projection method\n" );
              terminate(1);
          }

        } /* site */
      } /* direction */
    } /* parity */
  } // 3D/4D branching

  // free temporary storage
  destroy_mn_special(Wiggly_link);
  destroy_mn_special(Doubly_link);
  destroy_mn_special(staple);

} /* load_hyp_links */


// set HYP coefficients
void set_hyp_coeff( hyp_coeffs_t *hc, Real a1, Real a2, Real a3 ) {
  hc->alpha1 = a1;
  hc->alpha2 = a2;
  hc->alpha3 = a3;
}

// set group projection method: U(3), SU(3), SU(3) by Tr Max
void set_hyp_proj_method( hyp_coeffs_t *hc, int method, int hits ) {
  hc->proj_method = method;
  hc->hits = hits;
}
