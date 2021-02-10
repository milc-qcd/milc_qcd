/************************** w_loop2.c *******************************/
/* MIMD version 7 */
/* This version uses gathers to get the neighbors */
/* version of 3/16/94 by UMH */
/* 2/19/98 Version 5 port CD */

/* Computes time-like, off-axis Wilson loops on gauge configuration
   in axial gauge with time-like gauge fields of last time slice in
   all other time slices as well, instead of the unit matrix! */

//#define ONLYONEPATH  //set this flag to use only one of four possible paths

#include "../generic/generic_includes.h"
#include "hvy_qpot_includes.h"
#include "../include/openmp_defs.h"
#include <assert.h>

static hqp_geom *geom=NULL;
static double *wils_loop2=NULL;

static su3_matrix *s_link=NULL, *s_link_f=NULL, *t_link_f=NULL;

static inline void global_sums( int mysize );
static inline int hqp_buffer_dim (void );
static void hqp_free_buffers( void );
static void hqp_setup_buffers( int mysize );
static inline void output_all( char smtag[], int mi );

/* constructs elementary sqrt(2) transporters, then calls loop_rt */
static void fw_sqrt2_loops ( int mu1, int mu2, double *wils_loop2 );
static void bw_sqrt2_loops ( int mu1, int mu2, double *wils_loop2 );
/* loops over possible sqrt(2) loops */
static int sqrt2_loops( double *wils_loop2 );

/* constructs elementary sqrt(5) transporters, then calls loop_rt */
static void fw_sqrt5_loops ( int mu1, int mu2, double *wils_loop2 );
static void bw_sqrt5_loops ( int mu1, int mu2, double *wils_loop2 );
/* loops over possible sqrt(5) loops, two steps always only forward */
static int sqrt5_loops( double *wils_loop2 );

/* constructs all elementary forward sqrt(3) transporters, then calls loop_rt */
static void fffw_sqrt3_loops ( int mu1, int mu2, int mu3, double *wils_loop2 );
/* constructs all elementary two forward-one backward sqrt(3) transporters, then calls loop_rt */
static void ffbw_sqrt3_loops ( int mu1, int mu2, int mu3, double *wils_loop2 );
static void fbfw_sqrt3_loops ( int mu1, int mu2, int mu3, double *wils_loop2 );
static void bffw_sqrt3_loops ( int mu1, int mu2, int mu3, double *wils_loop2 );
/* not really needed, but good to have for bugfixing the two forward-one backward sqrt(3) transporters */
static void fbbw_sqrt3_loops ( int mu1, int mu2, int mu3, double *wils_loop2 );
static void bfbw_sqrt3_loops ( int mu1, int mu2, int mu3, double *wils_loop2 );
static void bbfw_sqrt3_loops ( int mu1, int mu2, int mu3, double *wils_loop2 );
/* loops over possible sqrt(3) loops, using only all forward or two forward-one backward transporters */
static int sqrt3_loops( double *wils_loop2 );

static int max11[3], max12[6], max111;
static int nr11,nr12,nr111;
static int mi;

void w_loop2(int tot_smear) {

  register int r,t;
  register int mu, mu1, mu2;
  int nr,avgr=0;
  int base_r=0;
  int disp3[4];

  char myname[] = "new_w_loop2";

  char smtag[MAXSMTAG]="";
#ifdef SMEARING
  sprintf(smtag,"_%d",tot_smear);
#endif

  geom = hqp_geometry( myname );
  mi = hqp_buffer_dim();
  geom->maxlen=maxc[TUP];
  hqp_setup_buffers( geom->maxlen*mi );

  /****************************************************************
   * 
   * Proceed to loop over directions and recursively construct the 
   * space-like segments and compute the Wilson loops with these
   *
   ***************************************************************/

  /* Do off-axis "sqrt(2)" loops */
  base_r+=sqrt2_loops ( wils_loop2+base_r );
  /* Do off-axis "sqrt(5)" loops */
  base_r+=sqrt5_loops ( wils_loop2+base_r );
  /* Do off-axis "sqrt(3)" loops */
  base_r+=sqrt3_loops ( wils_loop2+base_r );
  assert( (base_r==mi) );

  /****************************************************************
   * 
   * Proceed to normalizing and printing of the Wilson loops
   * Remark: do not average the directions (impossible if anisotropic),
   * instead print the direction into one column
   * can activate averaging of directions via compiler macro AVERAGE_DIRECTIONS
   *
   ***************************************************************/

  global_sums( geom->maxlen*mi );

  output_all( smtag, mi );

  free( wils_loop2);

} /* w_loop2 */


static void hqp_setup_buffers( int mysize ) { 

  wils_loop2 = hqp_alloc_dble_buffer( mysize );
}

static void hqp_free_buffers( void ) { 

  hqp_free_dble_buffer( wils_loop2 );
}

static inline void global_sums( int mysize ) {

  g_vecdoublesum(wils_loop2, mysize );
}

static inline int hqp_buffer_dim ( void ) {

  /* Determine number of distinct geometries to consider, ignoring anisotropy here */
  max11[XUP] = ( maxc[YUP]<maxc[ZUP] ? maxc[YUP] : maxc[ZUP] );
  max11[YUP] = ( maxc[ZUP]<maxc[XUP] ? maxc[ZUP] : maxc[XUP] );
  max11[ZUP] = ( maxc[XUP]<maxc[YUP] ? maxc[XUP] : maxc[YUP] );
  if (max11[XUP] > max_r/sqrt(2) ) max11[XUP]=ceil(max_r/sqrt(2)); 
  if (max11[YUP] > max_r/sqrt(2) ) max11[YUP]=ceil(max_r/sqrt(2)); 
  if (max11[ZUP] > max_r/sqrt(2) ) max11[ZUP]=ceil(max_r/sqrt(2)); 
  nr11 = max11[XUP] + max11[YUP] + max11[ZUP];

  /* First  sqrt(5) block -- two steps in lower  direction */
  max12[XUP  ] = ( maxc[YUP]/2<maxc[ZUP] ? maxc[YUP]/2 : maxc[ZUP] );
  max12[YUP  ] = ( maxc[XUP]/2<maxc[ZUP] ? maxc[XUP]/2 : maxc[ZUP] );
  max12[ZUP  ] = ( maxc[XUP]/2<maxc[YUP] ? maxc[XUP]/2 : maxc[YUP] );
  /* Second sqrt(5) block -- two steps in higher direction */
  max12[XUP+3] = ( maxc[YUP]<maxc[ZUP]/2 ? maxc[YUP] : maxc[ZUP]/2 );
  max12[YUP+3] = ( maxc[XUP]<maxc[ZUP]/2 ? maxc[XUP] : maxc[ZUP]/2 );
  max12[ZUP+3] = ( maxc[XUP]<maxc[YUP]/2 ? maxc[XUP] : maxc[YUP]/2 );
  if (max12[XUP  ] > max_r/sqrt(5) ) max12[XUP  ]=ceil(max_r/sqrt(5)); 
  if (max12[YUP  ] > max_r/sqrt(5) ) max12[XUP  ]=ceil(max_r/sqrt(5)); 
  if (max12[ZUP  ] > max_r/sqrt(5) ) max12[XUP  ]=ceil(max_r/sqrt(5)); 
  if (max12[XUP+3] > max_r/sqrt(5) ) max12[XUP+3]=ceil(max_r/sqrt(5)); 
  if (max12[YUP+3] > max_r/sqrt(5) ) max12[XUP+3]=ceil(max_r/sqrt(5)); 
  if (max12[ZUP+3] > max_r/sqrt(5) ) max12[XUP+3]=ceil(max_r/sqrt(5)); 
  nr12 = max12[XUP]   + max12[YUP]   + max12[ZUP]
       + max12[XUP+3] + max12[YUP+3] + max12[ZUP+3];

  max111 = ( maxc[YUP]<maxc[ZUP] ? maxc[XUP]<maxc[YUP] ? maxc[XUP] : maxc[YUP] : maxc[ZUP] );
  if (max111 > max_r/sqrt(3) ) max111=ceil(max_r/sqrt(5)); 
  nr111 = max111;

#ifdef DEBUG
  node0_printf("\n maxc  = %d, %d, %d",maxc[XUP],maxc[YUP],maxc[ZUP]);
  node0_printf("\n max11 = %d, %d, %d",max11[XUP],max11[YUP],max11[ZUP]);
  node0_printf("\n max12 = %d, %d, %d, %d, %d, %d",max12[XUP],max12[YUP],max12[ZUP],max12[XUP+3],max12[YUP+3],max12[ZUP+3]);
  node0_printf("\n max111= %d \n",max111);
  node0_printf(" Calculating %d sqrt(2) loops, %d sqrt(5) loops, and %d sqrt(3) loops\n\
                \n####################################\n\n",nr11,nr12,nr111);
#endif

  return( nr11 + nr12 + nr111 );
}

static inline void output_all( char smtag[], int mi ) {

  int avgr, nr, base_r;
  int *r;
  int disp[4]={0,0,0,0};
  double fact;

  for(disp[TUP]=1;disp[TUP]<=geom->maxlen;disp[TUP]++) { 

    avgr = nr = base_r = 0;

#ifndef ONLYONEPATH
  fact=0.25/3.;
#else
  fact=1./3;
#endif
    /* Normalize and print the sqrt(2) Wilson loops */
    for( int mu1=XUP; mu1<=YUP; mu1++) for( int mu2= mu1+1; mu2<=ZUP; mu2++){ 
      int mu=3-mu1-mu2;
      memset(disp, 0,3*sizeof(int));
      r=&(disp[(mu+1)%3]);
      for(*r=1;*r<=max11[mu];(*r)++){
        disp[(mu+2)%3]=*r;
        if ( hqp_disp_rsq_ok( disp, geom ) == 1 ) 
#ifndef AVERAGE_DIRECTIONS
          hqp_output_corr("WILS_LOOP2",smtag, disp, 
          wils_loop2[((*r-1)+base_r+mi*(disp[TUP]-1))]*fact);
#else 
          wils_loop2[(*r-1)+mi*(disp[TUP]-1)]+=((mu+1)%3>XUP?1.:-2.)* 
            wils_loop2[(*r-1)+base_r+mi*(disp[TUP]-1)]/3. ;
#endif
      }
      base_r += max11[mu];
    }
#ifdef AVERAGE_DIRECTIONS
    memset(disp, 0,3*sizeof(int));
    disp[TUP]--;
    r=&(disp[(XUP)%3]);
    for(*r=1;*r<=max11[XUP];(*r)++)
      if ( hqp_disp_rsq_ok( disp, geom ) == 1 ) {
        (*r)--;
        hqp_output_corr("WILS_LOOP2",smtag, disp, wils_loop2[((*r+avgr)+mi*(disp[TUP]))]*fact);
        (*r)++;
      }
    disp[TUP]++;
#endif

  nr = base_r;
  avgr=nr/3;

    /* Normalize and print the sqrt(5) Wilson loops */
    for( int mu1=XUP; mu1<=YUP; mu1++) for( int mu2= mu1+1; mu2<=ZUP; mu2++){ 
      int mu=3-mu1-mu2;
      memset(disp, 0,3*sizeof(int));
      r=&(disp[(mu1)%3]);
      for(*r=1;*r<=max12[3-mu1-mu2];(*r)++) {
        disp[(mu2)%3]=2*(*r);
        if ( hqp_disp_rsq_ok( disp, geom ) == 1 ) 
#ifndef AVERAGE_DIRECTIONS
          hqp_output_corr("WILS_LOOP2",smtag, disp, 
          wils_loop2[((*r-1)+base_r+mi*(disp[TUP]-1))]*fact);
#else 
          wils_loop2[(*r-1)+nr+mi*(disp[TUP]-1)]+=(mu2>YUP?+1:-5)* 
            wils_loop2[(*r-1)+base_r+mi*(disp[TUP]-1)]/6. ;
#endif
      }
      base_r += max12[3-mu1-mu2];
    }

    for( int mu2=XUP; mu2<=YUP; mu2++) for( int mu1= mu2+1; mu1<=ZUP; mu1++){
      int mu=6-mu1-mu2;
      memset(disp, 0,3*sizeof(int));
      r=&(disp[(mu1)%3]);
      for(*r=1;*r<=max12[6-mu1-mu2];(*r)++) {
        disp[xc[(mu2)%3]]=2*(*r);
#ifndef ONLYONEPATH
        if ( hqp_disp_rsq_ok( disp, geom ) == 1 ) 
#  ifndef AVERAGE_DIRECTIONS
          hqp_output_corr("WILS_LOOP2",smtag, disp, 
            wils_loop2[((*r-1)+base_r+mi*(disp[TUP]-1))]*fact);
#  else 
          wils_loop2[(*r-1)+nr+mi*(disp[TUP]-1)]+= 
            wils_loop2[(*r-1)+base_r+mi*(disp[TUP]-1)]/6. ;
#  endif
#endif
      }
      base_r += max12[6-mu1-mu2];
    }
#ifdef AVERAGE_DIRECTIONS
    memset(disp, 0,3*sizeof(int));
    disp[TUP]--;
    r=&(disp[XUP]);
    for(*r=1;*r<=max12[XUP];(*r)++)
      if ( hqp_disp_rsq_ok( disp, geom ) == 1 ) {
        (*r) += avgr-1;
        hqp_output_corr("WILS_LOOP2",smtag, disp, wils_loop2[((*r+nr-avgr)+mi*(disp[TUP]))]*fact);
        (*r) -= avgr-1;
      }
    disp[TUP]++;
#endif

    nr = base_r;
    avgr += (nr-3*avgr)/6;

    /* Normalize and print the sqrt(3) Wilson loops */
#ifndef ONLYONEPATH
  fact=0.0416666666666666667/3.;
#else
  fact=1./3.;
#endif
    r=&(disp[YUP]);
    for(*r=1;*r<=max111;(*r)++) {
      disp[XUP]=(*r);
      disp[ZUP]=(*r);
      if ( hqp_disp_rsq_ok( disp, geom ) == 1 ) {
#ifndef AVERAGE_DIRECTIONS
        hqp_output_corr("WILS_LOOP2",smtag, disp, 
          wils_loop2[((*r-1)+base_r+mi*(disp[TUP]-1))]*fact);
#else 
        disp[TUP]--;
        disp[XUP]=(*r-1+avgr);
        hqp_output_corr("WILS_LOOP2",smtag, disp, 
          wils_loop2[((*r-1)+base_r+mi*(disp[TUP]))]*fact);
        disp[TUP]++;
#endif
      }
    }
    base_r += max111;
    nr = base_r;
    assert( (base_r==mi && nr==mi) );

  }
}

/********************************************************************************
 *
 * Auxiliary functions for constructing and summing up 
 * forward and backward off-axis Wilson loops
 * for different off-axis displacements
 *
 *******************************************************************************/

static void fw_sqrt2_loops ( int mu1, int mu2, double *wils_loop2 ) {
  register int i;
  register site *s;
  int disp[]={0,0,0,0};
  msg_tag *mtag;

  /* First construct the "diagonal" link in the ( mu1, mu2) direction, then make one corner */
  mtag = start_gather_site( F_OFFSET(link[xc[mu2]]), sizeof(su3_matrix),xc[mu1], EVENANDODD, gen_pt[1] );
  wait_gather( mtag);
  FORALLSITES_OMP(i,s, default(shared) ){ mult_su3_nn( &(s->link[xc[mu1]]), (su3_matrix *)(gen_pt[1][i]), &(s->diag)); } END_LOOP_OMP;
  cleanup_gather( mtag);

  disp[xc[mu1]] = 1;
  disp[xc[mu2]] = 1;
  /* Recursively construct the space-like segments and compute the Wilson loops with that segment */
  loop_rt ( geom, max11[3-mu1-mu2], mi, disp, wils_loop2 );
}

static void bw_sqrt2_loops ( int mu1, int mu2, double *wils_loop2 ) {
  register int i;
  register site *s;
  int disp[]={0,0,0,0};
  int nr;
  msg_tag *mtag=NULL,*gmtag=NULL;

  /* Gather mu2 link across the "diagonal" link in the (mu1,-mu2) or (-mu1,mu2) direction */
  if ( xc[mu1] <= TUP ) { 
    disp[xc[mu1]] =  1; 
    disp[xc[OPP_DIR(mu2)]] = -1;
    nr=max11[3-mu1-OPP_DIR(mu2)];
  }else{
    disp[xc[OPP_DIR(mu1)]] = -1;
    disp[xc[mu2]] =  1; 
    nr=max11[3-mu2-OPP_DIR(mu1)];
  }
  if ( xc[mu1] <= TUP ) {
    /* First general gather, later make one corner */
    gmtag = start_general_gather_site( F_OFFSET(link[xc[OPP_DIR(mu2)]]), sizeof(su3_matrix), disp, EVENANDODD, gen_pt[1] );
  }else{
    /* First make one corner, then axis-gather it */
    FORALLSITES_OMP(i,s, default(shared) ){ mult_su3_an(&(s->link[xc[OPP_DIR(mu1)]]),  &(s->link[xc[mu2]]), &(s->staple)); } END_LOOP_OMP;
    mtag = start_gather_site( F_OFFSET(staple), sizeof(su3_matrix),  OPP_DIR(xc[OPP_DIR(mu1)]), EVENANDODD, gen_pt[1] );
  }

  if ( xc[mu1] <= TUP ) {
    /* wait for general gather, then make one corner */
    wait_general_gather( gmtag);
    FORALLSITES_OMP(i,s, default(shared) ){ mult_su3_na( &(s->link[xc[mu1]]), (su3_matrix *)(gen_pt[1][i]),&(s->diag)); } END_LOOP_OMP;
    cleanup_general_gather( gmtag);
  }else{
    /* wait for axis-gather of corner, then copy */
    wait_gather( mtag);
    FORALLSITES_OMP(i,s, default(shared) ){ su3mat_copy( (su3_matrix *)(gen_pt[1][i]), &(s->diag)); } END_LOOP_OMP;
    cleanup_gather( mtag);
  }

  /* Recursively construct the space-like segments and compute the Wilson loops with that segment */
  loop_rt ( geom, nr, mi, disp, wils_loop2 );
}

static void fw_sqrt5_loops ( int mu1, int mu2, double *wils_loop2 ) {
  register int i;
  register site *s;
  int nr, nu1,nu2,twolow=-1;
  int disp[]={0,0,0,0};
  msg_tag *mtag,*gmtag;

  if ( mu1 >= X3UP ) {
    nu1 = mu1 - X3UP;
    nu2 = mu2;
    twolow = ( nu1 < nu2 ? 1 : 0 );
    /* First construct the "diagonal" link in the (2*mu1,mu2) direction */
    mtag = start_gather_site( F_OFFSET(link[xc[nu1]]), sizeof(su3_matrix), xc[nu1], EVENANDODD, gen_pt[1] );
    /* Start gather of mu2-link from "2*mu1" */
    disp[xc[nu1]] = 2;
    gmtag = start_general_gather_site( F_OFFSET(link[xc[nu2]]), sizeof(su3_matrix), disp, EVENANDODD, gen_pt[4] );
  }else{ 
    /* First construct the "diagonal" link in the (mu1,2*mu2) direction */
    nu1=mu1;
    nu2=mu2-X3UP;
    twolow = ( nu1 < nu2 ? 0 : 1 );
    disp[xc[nu1]] = 1;
    disp[xc[nu2]] = 1;
    gmtag = start_general_gather_site( F_OFFSET(link[xc[nu2]]), sizeof(su3_matrix), disp, EVENANDODD, gen_pt[4] );
    mtag = start_gather_site( F_OFFSET(link[xc[nu2]]), sizeof(su3_matrix), xc[nu1], EVENANDODD, gen_pt[1] );
  }
  assert( (twolow > -1) );
  nr = max12[3-nu1-nu2+3*(1-twolow)];

  /* Wait for both gathers, then make corner */
  wait_gather( mtag);
  wait_general_gather( gmtag);
  FORALLSITES_OMP(i,s, default(shared) ){
    mult_su3_nn( &(s->link[xc[nu1]]), (su3_matrix *)(gen_pt[1][i]),&(s->staple));
    mult_su3_nn( &(s->staple), (su3_matrix *)(gen_pt[4][i]),&(s->diag));
  } END_LOOP_OMP;
  cleanup_gather( mtag);
  cleanup_general_gather( gmtag);
  disp[xc[(mu2%X3UP)]] += 1; 

  /* Recursively construct the space-like segments and compute the Wilson loops with that segment */
  loop_rt ( geom, nr, mi, disp, wils_loop2 );
}

static void bw_sqrt5_loops ( int mu1, int mu2, double *wils_loop2 ) {
  register int i;
  register site *s;
  int disp[]={0,0,0,0};
  int nr,nu1,nu2,twolow=-1;
  msg_tag *mtag[2],*gmtag;

  /* Gather mu2 link across the "diagonal" link in the (mu1,-mu2) direction */
  if ( mu1 >= X3UP ) {
    nu1=mu1-X3UP;
    nu2=OPP_DIR(mu2);
    twolow = ( nu1 < nu2 ? 1 : 0 );
    disp[xc[nu1]] = 2;
    disp[xc[nu2]] = -1;
  }else{
    nu1=OPP_DIR(mu1);
    nu2=mu2-X3UP;
    twolow = ( nu1 < nu2 ? 0 : 1 );
    disp[xc[nu1]] = -1;
    disp[xc[nu2]] = 2;
  }
  assert( (twolow > -1) );
  nr = max12[3-nu1-nu2+3*(1-twolow)];

  if ( mu1 >= X3UP ) {
    /* First construct the "diagonal" link in the (2*mu1,-mu2) direction */
    mtag[1] = start_gather_site( F_OFFSET(link[xc[nu1]]), sizeof(su3_matrix), xc[nu1], EVENANDODD, gen_pt[1] );
    /* First gather mu2-link from across the 2,-1 diagonal, later make one corner */
    gmtag = start_general_gather_site( F_OFFSET(link[xc[nu2]]), sizeof(su3_matrix), disp, EVENANDODD, gen_pt[4] );
  }else{
    /* Construct the "diagonal" link in the (-mu1,2*mu2) direction */
    /* First gather the last mu1 link, then make corner, then gahter corner */
    mtag[1] = start_gather_site( F_OFFSET(link[xc[nu2]]), sizeof(su3_matrix), xc[nu2], EVENANDODD, gen_pt[1] );
    wait_gather( mtag[1]);
    FORALLSITES_OMP(i,s, default(shared) ){
      mult_su3_an( &(s->link[xc[nu1]]), &(s->link[xc[nu2]]), &(s->diag) );
      mult_su3_nn( &(s->diag), (su3_matrix *)(gen_pt[1][i]),&(s->staple) );
    } END_LOOP_OMP;
    cleanup_gather( mtag[1]);
    mtag[0] = start_gather_site( F_OFFSET(staple), sizeof(su3_matrix), OPP_DIR(xc[nu1]), EVENANDODD, gen_pt[0] );
  }

  if ( mu1 >= X3UP ) {
    /* Wait for both gathers, then make corner */
    wait_gather( mtag[1]);
    wait_general_gather( gmtag);
    FORALLSITES_OMP(i,s, default(shared) ){
      mult_su3_nn( &(s->link[xc[nu1]]), (su3_matrix *)(gen_pt[1][i]),&(s->staple) );
      mult_su3_na( &(s->staple) , (su3_matrix *)(gen_pt[4][i]),&(s->diag));
    } END_LOOP_OMP;
    cleanup_gather( mtag[1]);
    cleanup_general_gather( gmtag);
  }else{
    /* wait for axis-gather of corner, then copy */
    wait_gather( mtag[0]);
    FORALLSITES_OMP(i,s, default(shared) ){ su3mat_copy( (su3_matrix *)(gen_pt[0][i]),&(s->diag)); } END_LOOP_OMP;
    cleanup_gather( mtag[0]);
  }

  /* Recursively construct the space-like segments and compute the Wilson loops with that segment */
  loop_rt ( geom, nr, mi, disp, wils_loop2 );
}

static void fffw_sqrt3_loops ( int mu1, int mu2, int mu3, double *wils_loop2 ) {
  register int i;
  register site *s;
  int disp[]={0,0,0,0};
  msg_tag *mtag,*gmtag;

#ifdef DEBUG
  node0_printf("Construct the body diagonal link in the (%d,%d,%d) direction\n",xc[mu1],xc[mu2],xc[mu3]);
#endif

  /* Gather mu3-link from across the mu1,mu2 diagonal */
  disp[xc[mu1]] = 1;
  disp[xc[mu2]] = 1;
  gmtag = start_general_gather_site( F_OFFSET(link[xc[mu3]]), sizeof(su3_matrix),disp, EVENANDODD, gen_pt[2] );

  /* Gather mu2-link along mu1 axis for first corner */
  mtag = start_gather_site( F_OFFSET(link[xc[mu2]]), sizeof(su3_matrix), xc[mu1], EVENANDODD, gen_pt[1] );

  /* Wait for gathers, then make lower 2d-corner in (mu1,mu2) direction, then lower body diagonal mu1,mu2,mu3 link */
  wait_gather( mtag);
  wait_general_gather( gmtag);
  FORALLSITES_OMP(i,s, default(shared) ){
    mult_su3_nn( &(s->link[xc[mu1]]), (su3_matrix *)(gen_pt[1][i]),&(s->staple));
    mult_su3_nn( &(s->staple), (su3_matrix *)(gen_pt[2][i]),&(s->diag));
  } END_LOOP_OMP;
  cleanup_gather( mtag);
  cleanup_general_gather( gmtag);
  disp[xc[mu3]] = 1;

  /* Recursively construct the space-like segments and compute the Wilson loops with that segment */
  loop_rt ( geom, max111, mi, disp, wils_loop2 );
}

static void ffbw_sqrt3_loops ( int mu1, int mu2, int mu3, double *wils_loop2 ) {
  register int i;
  register site *s;
  int disp[]={0,0,0,0};
  msg_tag *mtag,*gmtag;

#ifdef DEBUG
  node0_printf("Construct the body diagonal link in the (%d,%d,-%d) direction\n",xc[mu1],xc[mu2],xc[mu3]);
#endif

  /* Gather mu3-link from across the diagonal and one down along -mu3 */
  disp[xc[mu1]] =  1;
  disp[xc[mu2]] =  1;
  disp[xc[mu3]] = -1;
  gmtag = start_general_gather_site( F_OFFSET(link[xc[mu3]]), sizeof(su3_matrix),disp, EVENANDODD, gen_pt[2] );

  /* Gather for first corner */
  mtag = start_gather_site( F_OFFSET(link[xc[mu2]]), sizeof(su3_matrix), xc[mu1], EVENANDODD, gen_pt[1] );

  /* Make lower corner in (mu1,mu2) direction and then lower body diagonal mu1,mu2,-mu3 link */
  wait_gather( mtag);
  wait_general_gather( gmtag);
  FORALLSITES_OMP(i,s, default(shared) ){
      mult_su3_nn( &(s->link[xc[mu1]]), (su3_matrix *)(gen_pt[1][i]),&(s->staple));
      mult_su3_na( &(s->staple), (su3_matrix *)(gen_pt[2][i]), &(s->diag));
  } END_LOOP_OMP;
  cleanup_gather( mtag);
  cleanup_general_gather( gmtag);

  /* Recursively construct the space-like segments and compute the Wilson loops with that segment */
  loop_rt ( geom, max111, mi, disp, wils_loop2 );
}

static void fbfw_sqrt3_loops ( int mu1, int mu2, int mu3, double *wils_loop2 ) {
  register int i;
  register site *s;
  int disp[]={0,0,0,0};
  msg_tag *gmtag;

#ifdef DEBUG
  node0_printf("Construct the body diagonal link in the (%d,-%d,%d) direction\n",xc[mu1],xc[mu2],xc[mu3]);
#endif

  /* Make second corner in (-mu2,mu3) direction */
  FORALLSITES_OMP(i,s, default(shared) ){ mult_su3_an( &(s->link[xc[mu2]]), &(s->link[xc[mu3]]),&(s->staple)); } END_LOOP_OMP;

  /* Gather second corner from across the diagonal (mu1,-mu2) direction, then make lower body diagonal mu1,-mu2,mu3 link */
  disp[xc[mu1]] =  1;
  disp[xc[mu2]] = -1;
  gmtag = start_general_gather_site( F_OFFSET(staple), sizeof(su3_matrix), disp, EVENANDODD, gen_pt[1] );
  wait_general_gather( gmtag);
  FORALLSITES_OMP(i,s, default(shared) ){ mult_su3_nn( &(s->link[xc[mu1]]), (su3_matrix *)(gen_pt[1][i]),&(s->diag)); } END_LOOP_OMP;
  cleanup_general_gather( gmtag);
  disp[xc[mu3]] = 1;

  /* Recursively construct the space-like segments and compute the Wilson loops with that segment */
  loop_rt ( geom, max111, mi, disp, wils_loop2 );
}

static void bffw_sqrt3_loops ( int mu1, int mu2, int mu3, double *wils_loop2 ) {
  register int i;
  register site *s;
  int disp[]={0,0,0,0};
  msg_tag *mtag,*gmtag;

#ifdef DEBUG
  node0_printf("Construct the body diagonal link in the (-%d,%d,%d) direction\n",xc[mu1],xc[mu2],xc[mu3]);
#endif

  disp[xc[mu1]] = -1;
  disp[xc[mu2]] =  1;
  disp[xc[mu2]] =  1;
  /* Gather mu3-link along mu2 axis for second corner */
  mtag = start_gather_site( F_OFFSET(link[xc[mu3]]), sizeof(su3_matrix), xc[mu2], EVENANDODD, gen_pt[1] );

  /* Make first corner */ 
  FORALLSITES_OMP(i,s, default(shared) ){ mult_su3_an( &(s->link[xc[mu1]]), &(s->link[xc[mu2]]), &(s->diag)); } END_LOOP_OMP;

  /* Wait for gather, then make lower body diagonal -mu1,mu2,mu3 link */
  wait_gather( mtag);
  FORALLSITES_OMP(i,s, default(shared) ){ mult_su3_nn( &(s->diag), (su3_matrix *)(gen_pt[1][i]), &(s->staple) ); } END_LOOP_OMP;
  cleanup_gather( mtag);

  /* Gather body-digonal from across the -mu1 direction */
  mtag = start_gather_site( F_OFFSET(staple), sizeof(su3_matrix), xc[OPP_DIR(mu1)], EVENANDODD, gen_pt[0] );
  wait_gather( mtag);
  FORALLSITES_OMP(i,s, default(shared) ){ su3mat_copy( (su3_matrix *)(gen_pt[0][i]), &(s->diag) ); } END_LOOP_OMP;
  cleanup_gather( mtag);

  /* Recursively construct the space-like segments and compute the Wilson loops with that segment */
  loop_rt ( geom, max111, mi, disp, wils_loop2 );
}

static void fbbw_sqrt3_loops ( int mu1, int mu2, int mu3, double *wils_loop2 ) {
  register int i;
  register site *s;
  int disp[]={0,0,0,0};
  msg_tag *gmtag;

#ifdef DEBUG
  node0_printf("Construct the body diagonal link in the (%d,-%d,-%d) direction\n",xc[mu1],xc[mu2],xc[mu3]);
#endif

  /* Gather mu2-link from across the (mu1,-mu2) diagonal */
  disp[xc[mu1]] =  1;
  disp[xc[mu2]] = -1;
  gmtag = start_general_gather_site( F_OFFSET(link[xc[mu2]]), sizeof(su3_matrix),disp, EVENANDODD, gen_pt[1] );

  /* Make first corner in (mu1,-mu2) direction */
  wait_general_gather( gmtag);
  FORALLSITES_OMP(i,s, default(shared) ){ mult_su3_na( &(s->link[xc[mu1]]), (su3_matrix *)(gen_pt[1][i]),&(s->staple)); } END_LOOP_OMP;
  cleanup_general_gather( gmtag);
  disp[xc[mu3]] = -1;

  /* Gather mu3-link from across the (mu1,-mu2) diagonal and one down along -mu3, make body diagonal mu1,-mu2,-mu3 link */
  gmtag = start_general_gather_site( F_OFFSET(link[xc[mu3]]), sizeof(su3_matrix),disp, EVENANDODD, gen_pt[2] );
  wait_general_gather( gmtag);
  FORALLSITES_OMP(i,s, default(shared) ){ mult_su3_na( &(s->staple), (su3_matrix *)(gen_pt[2][i]),&(s->diag)); } END_LOOP_OMP;
  cleanup_general_gather( gmtag);

  /* Recursively construct the space-like segments and compute the Wilson loops with that segment */
  loop_rt ( geom, max111, mi, disp, wils_loop2 );
}

static void bfbw_sqrt3_loops ( int mu1, int mu2, int mu3, double *wils_loop2 ) {
  register int i;
  register site *s;
  int disp[]={0,0,0,0};
  msg_tag *mtag,*gmtag;

#ifdef DEBUG
  node0_printf("Construct the body diagonal link in the (-%d,%d,-%d) direction\n",xc[mu1],xc[mu2],xc[mu3]);
#endif

  /* General Gather mu3-link along (mu2,-mu3 diagonal) for second corner */
  disp[xc[mu2]] = 1;
  disp[xc[mu3]] = -1;
  gmtag = start_general_gather_site( F_OFFSET(link[xc[mu3]]), sizeof(su3_matrix), disp, EVENANDODD, gen_pt[1] );

  /* Make first corner */
  FORALLSITES_OMP(i,s, default(shared) ){ mult_su3_an( &(s->link[xc[mu1]]), &(s->link[xc[mu2]]), &(s->diag)); } END_LOOP_OMP;

  /* Wait for gather, then make lower body diagonal -mu1,mu2,-mu3 link */
  wait_general_gather( gmtag);
  FORALLSITES_OMP(i,s, default(shared) ){ mult_su3_na( &(s->diag), (su3_matrix *)(gen_pt[1][i]), &(s->staple) ); } END_LOOP_OMP;
  cleanup_general_gather( gmtag);
  disp[xc[mu1]] = -1;

  /* Gather body-digonal from across the -mu1 direction */
  mtag = start_gather_site( F_OFFSET(staple), sizeof(su3_matrix), xc[OPP_DIR(mu1)], EVENANDODD, gen_pt[0] );
  wait_gather( mtag);
  FORALLSITES_OMP(i,s, default(shared) ){ su3mat_copy( (su3_matrix *)(gen_pt[2][i]), &(s->diag) ); } END_LOOP_OMP;
  cleanup_gather( mtag);

  /* Recursively construct the space-like segments and compute the Wilson loops with that segment */
  loop_rt ( geom, max111, mi, disp, wils_loop2 );
}

static void bbfw_sqrt3_loops ( int mu1, int mu2, int mu3, double *wils_loop2 ) {
  register int i;
  register site *s;
  int disp[]={0,0,0,0};
  msg_tag *mtag,*gmtag;

#ifdef DEBUG
  node0_printf("Construct the body diagonal link in the (-%d,-%d,%d) direction\n",xc[mu1],xc[mu2],xc[mu3]);
#endif

  /* Gather mu1 link from -mu1 direction */
  mtag = start_gather_site( F_OFFSET(link[xc[mu1]]), sizeof(su3_matrix), OPP_DIR(xc[mu1]), EVENANDODD, gen_pt[1] );

  /* Make first corner in (-mu2,mu3) direction */
  FORALLSITES_OMP(i,s, default(shared) ){ mult_su3_an( &(s->link[xc[mu2]]), &(s->link[xc[mu3]]), &(s->staple) ); } END_LOOP_OMP;

  /* Gather corner from (-mu1,-mu2) direction */
  disp[xc[mu1]] = -1;
  disp[xc[mu2]] = -1;
  gmtag = start_general_gather_site( F_OFFSET(staple), sizeof(su3_matrix), disp, EVENANDODD, gen_pt[0] );

  /* Make body diagonal -mu1,-mu2,mu3 link */
  wait_gather( mtag);
  wait_general_gather( gmtag);
  FORALLSITES_OMP(i,s, default(shared) ){ mult_su3_an( (su3_matrix *)(gen_pt[1][i]), (su3_matrix *)(gen_pt[0][i]), &(s->diag) ); } END_LOOP_OMP;
  cleanup_gather( mtag);
  cleanup_general_gather( gmtag);
  disp[xc[mu3]] = 1;

  /* Recursively construct the space-like segments and compute the Wilson loops with that segment */
  loop_rt ( geom, max111, mi, disp, wils_loop2 );
}

static int sqrt2_loops ( double * wils_loop2 ) {
  register int mu1, mu2;
  int base_r=0;

  /* Use OPP_DIR(mu2) to indicate backward link */
  for( mu1=XUP; mu1<=YUP; mu1++) for( mu2= mu1+1; mu2<=ZUP; mu2++){
    /* Do off-axis "sqrt(2)" loops in ( mu1, mu2)-direction */
    fw_sqrt2_loops ( mu1,mu2, wils_loop2+base_r );
#ifndef ONLYONEPATH
    /* Do off-axis "sqrt(2)" loops in ( mu2, mu1)-direction */
    fw_sqrt2_loops ( mu2,mu1, wils_loop2+base_r );
    /* Do off-axis "sqrt(2)" loops in ( mu1,-mu2)-direction */
    bw_sqrt2_loops ( mu1,OPP_DIR(mu2), wils_loop2+base_r );
    /* Do off-axis "sqrt(2)" loops in (-mu2, mu1)-direction */
    bw_sqrt2_loops ( OPP_DIR(mu2),mu1, wils_loop2+base_r );
#endif
    /* increment base_r counter for loops in different mu1,mu2 planes */
    base_r += max11[3-mu1-mu2];
  } /* end loop over  mu1 < mu2 */

  return base_r;
}

static int sqrt3_loops ( double * wils_loop2 ) {
/* define four constants that indicate which of the three directions go backwards */
  int base_r=0;
  int i,mu,npaths=6;
  int **paths=NULL;

  /* allocate buffer for looping over different body diagonal paths */
  paths=(int**)malloc(npaths*sizeof(int*));
  assert(paths!=NULL);
  for ( i = 0; i < npaths; i++) { 
    paths[i]=(int*)malloc(3*sizeof(int)); 
    assert(paths[i]!=NULL);
    for ( mu=XUP; mu<=ZUP; mu++ ) paths[i][mu] = (i + ( i < 3 ? +mu : -mu ) )%3;
#if DEBUG
    node0_printf("path[%d] %d %d %d\n",i,paths[i][XUP],paths[i][YUP],paths[i][ZUP]);
#endif
  }

  for ( i = 0; i < npaths; i++) { 
    /* Do off-axis "sqrt(3)" loops in ( XUP, YUP, ZUP)-space */
    fffw_sqrt3_loops ( paths[i][XUP], paths[i][YUP], paths[i][ZUP], wils_loop2+base_r );
#ifndef ONLYONEPATH
    /* Do off-axis "sqrt(3)" loops in ( XUP, YUP, -ZUP)-space */
    ffbw_sqrt3_loops ( paths[i][XUP], paths[i][YUP], paths[i][ZUP], wils_loop2+base_r ); 
    /* Same contribution as off-axis "sqrt(3)" loops in ( XUP, -YUP, -ZUP)-space */ // fbbw_sqrt3_loops ( paths[i][XUP], paths[i][YUP], paths[i][ZUP], wils_loop2+base_r ); 
    /* Do off-axis "sqrt(3)" loops in ( XUP, -YUP, ZUP)-space */
    fbfw_sqrt3_loops ( paths[i][XUP], paths[i][YUP], paths[i][ZUP], wils_loop2+base_r ); 
    /* Same contribution as off-axis "sqrt(3)" loops in ( -XUP, YUP, -ZUP)-space */ // bfbw_sqrt3_loops ( paths[i][XUP], paths[i][YUP], paths[i][ZUP], wils_loop2+base_r ); 
    /* Do off-axis "sqrt(3)" loops in ( XUP, -YUP, ZUP)-space */
    bffw_sqrt3_loops ( paths[i][XUP], paths[i][YUP], paths[i][ZUP], wils_loop2+base_r );
    /* Same contribution as off-axis "sqrt(3)" loops in ( -XUP, -YUP, ZUP)-space */ // bbfw_sqrt3_loops ( paths[i][XUP], paths[i][YUP], paths[i][ZUP], wils_loop2+base_r );
#else
    break;
#endif
  }
  for ( i = 0; i < npaths; i++) { free(paths[i]); }
  free(paths);

  base_r=max111;
  return base_r;
}

static int sqrt5_loops ( double * wils_loop2 ) {
  register int mu1, mu2;
  int base_r=0;
  /* Split this into two separate loops, mu1<mu2, or mu2<mu1 */

  /* Dirty trick: use DIR3(mu1) to indicate double link for the subroutines */
  /* Use OPP_DIR(mu2) to indicate backward link */
  for( mu1=XUP; mu1<=YUP; mu1++) for( mu2= mu1+1; mu2<=ZUP; mu2++){
    /* Do off-axis "sqrt(5)" loops in ( 2*mu1, mu2)-direction */
    fw_sqrt5_loops ( DIR3(mu1),mu2, wils_loop2+base_r );
#ifndef ONLYONEPATH
    /* Do off-axis "sqrt(5)" loops in ( mu2, 2*mu1)-direction */
    fw_sqrt5_loops ( mu2,DIR3(mu1), wils_loop2+base_r );
    /* Do off-axis "sqrt(5)" loops in ( 2*mu1,-mu2)-direction */
    bw_sqrt5_loops ( DIR3(mu1),OPP_DIR(mu2) , wils_loop2+base_r );
    /* Do off-axis "sqrt(5)" loops in (-mu2, 2*mu1)-direction */
    bw_sqrt5_loops ( OPP_DIR(mu2),DIR3(mu1) , wils_loop2+base_r );
#endif
    /* increment base_r counter for loops in different mu1,mu2 planes */
    base_r += max12[3-mu1-mu2];
  } /* end loop over  mu1 <  mu2 */

  for( mu2=XUP; mu2<=YUP; mu2++) for( mu1= mu2+1; mu1<=ZUP; mu1++){
#ifndef ONLYONEPATH
    /* Do off-axis "sqrt(5)" loops in ( 2*mu1, mu2)-direction */
    fw_sqrt5_loops ( DIR3(mu1),mu2, wils_loop2+base_r );
    /* Do off-axis "sqrt(5)" loops in ( mu2, 2*mu1)-direction */
    fw_sqrt5_loops ( mu2,DIR3(mu1), wils_loop2+base_r );
    /* Do off-axis "sqrt(5)" loops in ( 2*mu1,-mu2)-direction */
    bw_sqrt5_loops ( DIR3(mu1),OPP_DIR(mu2) , wils_loop2+base_r );
    /* Do off-axis "sqrt(5)" loops in (-mu2, 2*mu1)-direction */
    bw_sqrt5_loops ( OPP_DIR(mu2),DIR3(mu1) , wils_loop2+base_r );
#endif
    /* increment base_r counter for loops in different mu1,mu2 planes */
    base_r += max12[6-mu1-mu2];
  } /* end loop over  mu2 <  mu1 */

  return base_r;
}

