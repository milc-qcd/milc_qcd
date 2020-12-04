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

/* makes Wilson loops and adds them up into the array */
static void make_loops ( int t,
                         Real *wils_loop2, 
                         su3_matrix *s_link, 
                         su3_matrix *s_link_f, 
                         su3_matrix *t_link_f ) ; 

/* loops over r and t, extends Wilson lines as needed, then calls make_loops */
static void loop_rt ( int nr, int disp[], Real *wils_loop2 ) ;

/* constructs elementary sqrt(2) transporters, then calls loop_rt */
static void fw_sqrt2_loops ( int mu1, int mu2, Real *wils_loop2 );
static void bw_sqrt2_loops ( int mu1, int mu2, Real *wils_loop2 );
/* loops over possible sqrt(2) loops */
static int sqrt2_loops( Real *wils_loop2 );

/* constructs elementary sqrt(5) transporters, then calls loop_rt */
static void fw_sqrt5_loops ( int mu1, int mu2, Real *wils_loop2 );
static void bw_sqrt5_loops ( int mu1, int mu2, Real *wils_loop2 );
/* loops over possible sqrt(5) loops, two steps always only forward */
static int sqrt5_loops( Real *wils_loop2 );

/* constructs all elementary forward sqrt(3) transporters, then calls loop_rt */
static void fffw_sqrt3_loops ( int mu1, int mu2, int mu3, Real *wils_loop2 );
/* constructs all elementary two forward-one backward sqrt(3) transporters, then calls loop_rt */
static void ffbw_sqrt3_loops ( int mu1, int mu2, int mu3, Real *wils_loop2 );
static void fbfw_sqrt3_loops ( int mu1, int mu2, int mu3, Real *wils_loop2 );
static void bffw_sqrt3_loops ( int mu1, int mu2, int mu3, Real *wils_loop2 );
/* not really needed, but good to have for bugfixing the two forward-one backward sqrt(3) transporters */
static void fbbw_sqrt3_loops ( int mu1, int mu2, int mu3, Real *wils_loop2 );
static void bfbw_sqrt3_loops ( int mu1, int mu2, int mu3, Real *wils_loop2 );
static void bbfw_sqrt3_loops ( int mu1, int mu2, int mu3, Real *wils_loop2 );
/* loops over possible sqrt(3) loops, using only all forward or two forward-one backward transporters */
static int sqrt3_loops( Real *wils_loop2 );

static int max11[3], max12[6], max111;
static int nrmax;
static int nct;

void w_loop2(int tot_smear) {

  register int r,t;
  register int mu, mu1, mu2;
  double max_r2;
  int nr11,nr12,nr111;
  int nr,avgr=0;
  int base_r=0;
  int disp3[4];
  Real *wils_loop2=NULL;
  char myname[] = "new_w_loop2";

  /*******************************************************
   *
   * Start of header:
   *
   * Text output concerning the geometry and its consistency.
   * In particular, the correlation directions 'xc[mu]' are
   * permuted from the original lattice directions 'mu'
   * to permit static quark correlations with arbitrary 
   * correlation time direction and arbitrary anisotropic 
   * direction.
   *
   ******************************************************/

#ifdef AVERAGE_DIRECTIONS
#ifdef ANISOTROPY
  assert( (ani_dir==xc[TUP]) );
#endif
  assert( (maxc[XUP]==maxc[YUP] && maxc[XUP] == maxc[ZUP] ) );
#endif

  node0_printf("%s with %d times smeared links\n",myname,tot_smear);
#ifdef DEBUG
  node0_printf(" [min,max](CT)[%d] = [%d,%d], \
                \n max(C0)[%x] = %d, max(C1)[%x] = %d, max(C2)[%x] = %d, \
                \n max_r = %.3g\n",cor_dir, min_ct,maxc[TUP], xc[XUP],maxc[XUP] ,xc[YUP],maxc[YUP] ,xc[ZUP],maxc[ZUP] ,max_r);
#endif
  assert ( maxc[XUP] <= nc[XUP]/2 && maxc[YUP] <= nc[YUP]/2 && maxc[ZUP] <= nc[ZUP]/2 );  
  assert ( min_ct <= maxc[TUP] && maxc[TUP] <= nc[TUP] );
  assert ( max_r >= 0 );  
  if ( max_r == 0 ) {   
#ifndef ANISOTROPY
    max_r=sqrt(maxc[XUP]*maxc[XUP] + maxc[YUP]*maxc[YUP] + maxc[ZUP]*maxc[ZUP] ) + 1e-6;
#else
#ifndef ONEDIM_ANISO_TEST
    max_r=sqrt( ( ani_dir != cor_dir ? ani_xiq*ani_xiq : 1 ) * maxc[XUP]*maxc[XUP] + maxc[YUP]*maxc[YUP] + maxc[ZUP]*maxc[ZUP] ) + 1e-6;
#else
    max_r=sqrt( ( ani_dir != cor_dir ? ani_xiq*ani_xiq : iso_xiq ) * maxc[XUP]*maxc[XUP] + (iso_xiq*iso_xiq)*(maxc[YUP]*maxc[YUP] + maxc[ZUP]*maxc[ZUP]) ) + 1e-6;
#endif
#endif
    node0_printf("max_r being automatically  reset to = %g\n",max_r);
  }
  for (mu=XUP; mu <TUP; mu++) {
#ifndef ANISOTROPY
    if ( max_r < maxc[mu] ) { 
      maxc[mu] = ceil(max_r);
      node0_printf("Decrease to max(C%d) = %d\n",mu,maxc[mu]);
#else
    if ( mu==XUP && max_r < ( ani_dir != cor_dir ? ani_xiq : 1 ) * maxc[mu] ) { 
      maxc[mu] = ceil(max_r/ani_xiq);
#ifndef ONEDIM_ANISO_TEST
      if ( max_r < ( mu == XUP && ani_dir != cor_dir ? ani_xiq : iso_xiq ) * maxc[mu] ) { 
        maxc[mu] = ceil(max_r/( mu == XUP && ani_dir != cor_dir ? ani_xiq : iso_xiq )); 
        node0_printf("Decrease to max(C%d) = %d\n",mu,maxc[mu]);
      } // END if ( max_r < ( mu == XUP && ani_dir != cor_dir ? ani_xiq : iso_xiq ) * maxc[mu] )
#else
      node0_printf("Decrease to max(C%d) = %d\n",mu,maxc[mu]);
#endif
#endif
    } // END if ( max_r < maxc[mu] ) or if ( mu==XUP && max_r < ( ani_dir != cor_dir ? ani_xiq : 1 ) * maxc[mu] )
  } // END for (mu=XUP; mu <TUP; mu++) 
  max_r2 = max_r*max_r;

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
  nrmax = nr11 + nr12 + nr111;
  nct = maxc[TUP];

#ifdef DEBUG
  node0_printf("\n maxc  = %d, %d, %d",maxc[XUP],maxc[YUP],maxc[ZUP]);
  node0_printf("\n max11 = %d, %d, %d",max11[XUP],max11[YUP],max11[ZUP]);
  node0_printf("\n max12 = %d, %d, %d, %d, %d, %d",max12[XUP],max12[YUP],max12[ZUP],max12[XUP+3],max12[YUP+3],max12[ZUP+3]);
  node0_printf("\n max111= %d \n",max111);
  node0_printf(" Calculating %d sqrt(2) loops, %d sqrt(5) loops, and %d sqrt(3) loops\n\
                \n####################################\n\n",nr11,nr12,nr111);
#endif

  /****************************************************************
   * 
   * End of header
   * Setup of buffers and prefabrication of 
   * elements to insert in static FORCE calculation (not active)
   * 
   ***************************************************************/

  wils_loop2 = (Real *)malloc(nct*nrmax*sizeof(Real));
  assert(wils_loop2!=NULL);
  memset (wils_loop2,0,nct*nrmax*sizeof(Real));

  /****************************************************************
   * 
   * End of setup 
   * Proceed to loop over directions and 
   * recursively construct the space-like segments and compute
   * the Wilson loops with that segment 
   *
   ***************************************************************/

  /* Do off-axis "sqrt(2)" loops */
  base_r+=sqrt2_loops ( wils_loop2+base_r );
  /* Do off-axis "sqrt(5)" loops */
  base_r+=sqrt5_loops ( wils_loop2+base_r );
  /* Do off-axis "sqrt(3)" loops */
  base_r+=sqrt3_loops ( wils_loop2+base_r );
  assert( (base_r==nrmax) );

  /****************************************************************
   * 
   * End of calculation 
   * Proceed to normalizing and printing of the Wilson loops
   * Remark: do not average the directions (impossible if anisotropic),
   * instead print the direction into one column
   * can activate averaging of directions via compiler macro AVERAGE_DIRECTIONS
   *
   ***************************************************************/

  g_vecfloatsum(wils_loop2,nct*nrmax);
  /* Normalize and print the sqrt(2) Wilson loops */
  for(t=0;t<nct;t++) {
    avgr = nr = base_r = 0;
    for( mu1=XUP; mu1<=YUP; mu1++) for( mu2= mu1+1; mu2<=ZUP; mu2++){ mu=3-mu1-mu2;
      for(r=0;r<max11[mu];r++) {
#ifndef ONLYONEPATH
        wils_loop2[r+base_r+nrmax*t] = ((double)wils_loop2[r+base_r+nrmax*t]/(double)(12*volume));
#else
        wils_loop2[r+base_r+nrmax*t] = ((double)wils_loop2[r+base_r+nrmax*t]/(double)( 3*volume));
#endif
#ifndef AVERAGE_DIRECTIONS
        memset(disp3, 0,4*sizeof(int));
        disp3[xc[(mu+1)%3]]=r+1;
        disp3[xc[(mu+2)%3]]=r+1;
        node0_printf("WILS_LOOP2_%d   %d  %d  %d   %d \t%e\n", tot_smear, 
                      disp3[xc[XUP]],disp3[xc[YUP]],disp3[xc[ZUP]], t+1, (double)wils_loop2[r+base_r+nrmax*t]);
#else
        wils_loop2[r+nr+nrmax*t] += ( mu > XUP ? wils_loop2[r+base_r+nrmax*t] : 0 );
#endif
      }
      base_r += max11[mu];
    }
#ifdef AVERAGE_DIRECTIONS
    for(r=0;r<max11[XUP];r++){
      node0_printf("WILS_LOOP2_%d  %d  %d \t%e\n", tot_smear, r+avgr, t, (double)wils_loop2[r+nr+nrmax*t]/3.);
    }
#endif
    nr = base_r;
    avgr=nr/3;
    /* Normalize and print the sqrt(5) Wilson loops */
    for( mu1=XUP; mu1<=YUP; mu1++) for( mu2= mu1+1; mu2<=ZUP; mu2++){ mu=3-mu1-mu2;
      for(r=0;r<max12[3-mu1-mu2];r++) {
#ifndef ONLYONEPATH
        wils_loop2[r+base_r+nrmax*t] = ((double)wils_loop2[r+base_r+nrmax*t]/(double)(12*volume));
#else
        wils_loop2[r+base_r+nrmax*t] = ((double)wils_loop2[r+base_r+nrmax*t]/(double)( 3*volume));
#endif
#ifndef AVERAGE_DIRECTIONS
        memset(disp3, 0,4*sizeof(int));
        disp3[xc[(mu1)%3]]=r+1;
        disp3[xc[(mu2)%3]]=2*(r+1);
        node0_printf("WILS_LOOP2_%d   %d  %d  %d   %d \t%e\n", tot_smear, 
                      disp3[xc[XUP]],disp3[xc[YUP]],disp3[xc[ZUP]], t+1, (double)wils_loop2[r+base_r+nrmax*t]);
#else
        wils_loop2[r+nr+nrmax*t] += ( mu2 > YUP ? wils_loop2[r+base_r+nrmax*t] : 0 );
#endif
      }
      base_r += max12[3-mu1-mu2];
    }

    for( mu2=XUP; mu2<=YUP; mu2++) for( mu1= mu2+1; mu1<=ZUP; mu1++){ mu=6-mu1-mu2;
      for(r=0;r<max12[6-mu1-mu2];r++) {
#ifndef ONLYONEPATH
        wils_loop2[r+base_r+nrmax*t] = ((double)wils_loop2[r+base_r+nrmax*t]/(double)(12*volume));
#ifndef AVERAGE_DIRECTIONS
        memset(disp3, 0,4*sizeof(int));
        disp3[xc[(mu1)%3]]=r+1;
        disp3[xc[(mu2)%3]]=2*(r+1);
        node0_printf("WILS_LOOP2_%d   %d  %d  %d   %d \t%e\n", tot_smear, 
                      disp3[xc[XUP]],disp3[xc[YUP]],disp3[xc[ZUP]], t+1, (double)wils_loop2[r+base_r+nrmax*t]);
#else
        wils_loop2[r+nr+nrmax*t] += wils_loop2[r+base_r+nrmax*t];
#endif
#endif
      }
      base_r += max12[6-mu1-mu2];
    }
#ifdef AVERAGE_DIRECTIONS
    for(r=0;r<max12[XUP];r++){
      node0_printf("WILS_LOOP2_%d  %d  %d \t%e\n", tot_smear, r+avgr, t, (double)wils_loop2[r+nr+nrmax*t]/6.);
    }
#endif
    nr = base_r;
    avgr += (nr-3*avgr)/6;

    /* Normalize and print the sqrt(3) Wilson loops */
    for(r=0;r<max111;r++) {
#ifndef ONLYONEPATH
        wils_loop2[r+base_r+nrmax*t] = ((double)wils_loop2[r+base_r+nrmax*t]/(double)(72*volume));
#else
        wils_loop2[r+base_r+nrmax*t] = ((double)wils_loop2[r+base_r+nrmax*t]/(double)( 3*volume));
#endif
#ifdef AVERAGE_DIRECTIONS
        node0_printf("WILS_LOOP2_%d  %d  %d \t%e\n", tot_smear,
                       r+avgr, t, (double)wils_loop2[r+base_r+nrmax*t]);
#else
        node0_printf("WILS_LOOP2_%d   %d  %d  %d   %d \t%e\n", tot_smear,
                       r+1,r+1,r+1, t+1, (double)wils_loop2[r+base_r+nrmax*t]);
#endif
    }
    base_r += max111;
    nr = base_r;
    assert( (base_r==nrmax && nr==nrmax) );
  }

  free( wils_loop2);

} /* w_loop2 */


/********************************************************************************
 *
 * Auxiliary functions for constructing and summing up 
 * forward and backward off-axis Wilson loops
 * for different off-axis displacements
 *
 *******************************************************************************/

static void make_loops ( int t, 
                  Real *wils_loop2,
                  su3_matrix *s_link,
                  su3_matrix *s_link_f,
                  su3_matrix *t_link_f ) {
  register int i;
  register site *s;
  su3_matrix tmat1,tmat2;

  FORALLSITES_OMP(i,s, private(tmat1,tmat2) reduction(+:wils_loop2[nrmax*t]) ){
    /* If the loop extends past t = nc[TUP] - 1 the temporal axial gauge link is nontrivial */
    if( (site_coord(s,xc[TUP])+t+1)>=nc[TUP] ){
      mult_su3_nn( &(s->link[xc[TUP]]), s_link_f+i, &tmat1);
      mult_su3_na( &tmat1, t_link_f+i, &tmat2);
      wils_loop2[nrmax*t] += realtrace_su3( &tmat2, s_link+i);
    }else{
      wils_loop2[nrmax*t] += realtrace_su3( s_link_f+i, s_link+i);
    }
  } END_LOOP_OMP;

} /* make_loops */


static void loop_rt ( int nr, int disp[], Real *wils_loop2 ) {
  register int i,t,r;
  register site *s;
  su3_matrix *s_link=NULL, *s_link_f=NULL, *t_link_f=NULL;
  msg_tag *mtag, *gmtag;

  /* Allocate space for space-link product and shifted time-links */
  s_link   = (su3_matrix *)malloc(sites_on_node*sizeof(su3_matrix));
  s_link_f = (su3_matrix *)malloc(sites_on_node*sizeof(su3_matrix));
  t_link_f = (su3_matrix *)malloc(sites_on_node*sizeof(su3_matrix));
  assert( ( s_link  !=NULL ) && ( s_link_f!=NULL ) && ( t_link_f!=NULL ) );
  memset (s_link  ,0,sites_on_node*sizeof(su3_matrix));
  memset (s_link_f,0,sites_on_node*sizeof(su3_matrix));
  memset (t_link_f,0,sites_on_node*sizeof(su3_matrix));

  for(r=0;r<nr;r++){

    if( r==0 ){
      /* Start gather of time-like links across the diagonal. */
      FORALLSITES_OMP(i,s, default(shared) ){ su3mat_copy( &(s->link[xc[TUP]]), t_link_f+i); } END_LOOP_OMP;
      gmtag = start_general_gather_field( (void *)t_link_f, sizeof(su3_matrix), disp, EVENANDODD, gen_pt[4] );
      FORALLSITES_OMP(i,s, default(shared) ){ su3mat_copy( &(s->diag), s_link+i); } END_LOOP_OMP;
    }else{
      wait_general_gather( gmtag);
      FORALLSITES_OMP(i,s, default(shared) ){ su3mat_copy( (su3_matrix *)(gen_pt[4][i]), &(s->staple)); } END_LOOP_OMP;

      /* Inbetween gather time-like links across the diagonal. */
      cleanup_general_gather( gmtag);
      gmtag = start_general_gather_field( (void *)t_link_f, sizeof(su3_matrix), disp, EVENANDODD,gen_pt[4]);

      FORALLSITES_OMP(i,s, default(shared) ){ mult_su3_nn( &(s->diag), &(s->staple), s_link+i); } END_LOOP_OMP;
    }
    FORALLSITES_OMP(i,s, default(shared) ){ su3mat_copy( s_link+i, s_link_f+i); } END_LOOP_OMP;

    /* Start gather of forward space-like segments */
    mtag = start_gather_field( (void *)s_link_f, sizeof(su3_matrix), xc[TUP], EVENANDODD, gen_pt[0] );

    /* Collect forward time-like links. */
    wait_general_gather( gmtag);
    FORALLSITES_OMP(i,s, default(shared) ){ su3mat_copy( (su3_matrix *)(gen_pt[4][i]), &(s->staple)); } END_LOOP_OMP;
    FORALLSITES_OMP(i,s, default(shared) ){ su3mat_copy( &(s->staple), t_link_f+i); } END_LOOP_OMP;

    /* Inbetween gather space-links across the diagonal for next r. */
    cleanup_general_gather( gmtag);
    if( r<(nr-1) ){
      gmtag = start_general_gather_field( (void *)s_link, sizeof(su3_matrix), disp, EVENANDODD, gen_pt[4] );
    }

    /* Recursively compute the Wilson loops of different time extent */
    for(t=0;t<nct;t++){

      /* Collect forward space-like segments */
      wait_gather( mtag);
      FORALLSITES_OMP(i,s, default(shared) ){ su3mat_copy( (su3_matrix *)(gen_pt[0][i]), &(s->staple)); } END_LOOP_OMP;
      FORALLSITES_OMP(i,s, default(shared) ){ su3mat_copy( &(s->staple), s_link_f+i); } END_LOOP_OMP;

      /* Start gather for next t, if still needed. */
      if( t<(nct-1) ){
        restart_gather_field( (void *)s_link_f, sizeof(su3_matrix), xc[TUP], EVENANDODD, gen_pt[0], mtag );
      }else{
        cleanup_gather( mtag);
      }

      /* Finally, compute the Wilson loops. */
      make_loops (t,(wils_loop2+r), s_link, s_link_f, t_link_f );

    } /* end loop over t */
  } /* end loop over r */

  free(s_link  );
  free(s_link_f);
  free(t_link_f);
} /* loop_rt */

static void fw_sqrt2_loops ( int mu1, int mu2, Real *wils_loop2 ) {
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
  loop_rt ( max11[3-mu1-mu2], disp, wils_loop2 );
}

static void bw_sqrt2_loops ( int mu1, int mu2, Real *wils_loop2 ) {
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
  loop_rt ( nr, disp, wils_loop2 );
}

static void fw_sqrt5_loops ( int mu1, int mu2, Real *wils_loop2 ) {
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
  loop_rt ( nr, disp, wils_loop2 );
}

static void bw_sqrt5_loops ( int mu1, int mu2, Real *wils_loop2 ) {
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
  loop_rt ( nr, disp, wils_loop2 );
}

static void fffw_sqrt3_loops ( int mu1, int mu2, int mu3, Real *wils_loop2 ) {
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
  loop_rt ( max111, disp, wils_loop2 );
}

static void ffbw_sqrt3_loops ( int mu1, int mu2, int mu3, Real *wils_loop2 ) {
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
  loop_rt ( max111, disp, wils_loop2 );
}

static void fbfw_sqrt3_loops ( int mu1, int mu2, int mu3, Real *wils_loop2 ) {
  register int i;
  register site *s;
  int disp[]={0,0,0,0};
  msg_tag *gmtag;

#ifdef DEBUG
  node0_printf("Construct the body diagonal link in the (%d,-%d,%d) direction\n",xc[mu1],xc[mu2],xc[mu3]);
#endif

  /* Make second corner in (-mu2,mu3) direction */
  FORALLSITES_OMP(i,s, default(shared) ){ mult_su3_an( &(s->link[mu2]), &(s->link[mu3]),&(s->staple)); } END_LOOP_OMP;

  /* Gather second corner from across the diagonal (mu1,-mu2) direction, then make lower body diagonal mu1,-mu2,mu3 link */
  disp[xc[mu1]] =  1;
  disp[xc[mu2]] = -1;
  gmtag = start_general_gather_site( F_OFFSET(staple), sizeof(su3_matrix), disp, EVENANDODD, gen_pt[1] );
  wait_general_gather( gmtag);
  FORALLSITES_OMP(i,s, default(shared) ){ mult_su3_nn( &(s->link[xc[mu1]]), (su3_matrix *)(gen_pt[1][i]),&(s->diag)); } END_LOOP_OMP;
  cleanup_general_gather( gmtag);
  disp[xc[mu3]] = 1;

  /* Recursively construct the space-like segments and compute the Wilson loops with that segment */
  loop_rt ( max111, disp, wils_loop2 );
}

static void bffw_sqrt3_loops ( int mu1, int mu2, int mu3, Real *wils_loop2 ) {
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
  loop_rt ( max111, disp, wils_loop2 );
}

static void fbbw_sqrt3_loops ( int mu1, int mu2, int mu3, Real *wils_loop2 ) {
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
  loop_rt ( max111, disp, wils_loop2 );
}

static void bfbw_sqrt3_loops ( int mu1, int mu2, int mu3, Real *wils_loop2 ) {
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
  loop_rt ( max111, disp, wils_loop2 );
}

static void bbfw_sqrt3_loops ( int mu1, int mu2, int mu3, Real *wils_loop2 ) {
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
  loop_rt ( max111, disp, wils_loop2 );
}

static int sqrt2_loops ( Real * wils_loop2 ) {
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

static int sqrt3_loops ( Real * wils_loop2 ) {
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

static int sqrt5_loops ( Real * wils_loop2 ) {
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

