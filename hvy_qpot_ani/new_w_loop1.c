/************************** w_loop1.c *******************************/
/* MIMD version 7 */
/* This version uses gathers to get the neighbors */
/* version of 3/14/94 by UMH */
/* 2/19/98 Version 5 port CD */

/* Computes time-like, planar Wilson loops on gauge configuration
   in axial gauge with time-like gauge fields of last time slice in
   all other time slices as well, instead of the unit matrix! */

#include "../generic/generic_includes.h"        /* definitions files and prototypes */
#include "./hvy_qpot_includes.h"
#include "../include/openmp_defs.h"
#include <assert.h>

static void make_loops ( int t,
                         Real *wils_loop1, 
                         su3_matrix *s_link, 
                         su3_matrix *s_link_f, 
                         su3_matrix *t_link_f ) ;

static int onaxis_loops( Real *wils_loop1 );

static int nrmax,nct;

void w_loop1(int tot_smear) {
  register int mu,r,t;
  int base_r=0;
  int disp3[4];
  double max_r2;
  Real *wils_loop1=NULL,ftmp;
  char myname[] = "new_w_loop1";

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
  nrmax = maxc[XUP]+maxc[YUP]+maxc[ZUP];
  nct=maxc[TUP];

//  node0_printf("####################################\n\n");

  /****************************************************************
   * 
   * End of header
   * Setup of buffers and prefabrication of 
   * elements to insert in static FORCE calculation (not active)
   * 
   ***************************************************************/

  wils_loop1 = (Real *)malloc(nct*nrmax*sizeof(Real));
  assert (wils_loop1!=NULL );
  memset (wils_loop1,0,nct*nrmax*sizeof(Real));

  /****************************************************************
  * 
   * End of setup 
   * Proceed to loop over directions and 
   * recursively construct the space-like segments and compute
   * the Wilson loops with that segment 
   *
   ***************************************************************/

  base_r += onaxis_loops ( wils_loop1+base_r ) ;
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

  g_vecfloatsum(wils_loop1,nct*nrmax);
  for(t=0;t<nct;t++){ 
    for( mu=XUP, base_r=0; mu<=ZUP; base_r+=maxc[mu],mu++){
      for(r=0;r<maxc[mu];r++){
        wils_loop1[r+base_r+nrmax*t] = ((double)wils_loop1[r+base_r+nrmax*t]/(double)(3*volume));
#ifndef AVERAGE_DIRECTIONS
        memset(disp3, 0,4*sizeof(int));
        disp3[xc[(mu)%3]]=r+1;
        node0_printf("WILS_LOOP1_%d   %d  %d  %d   %d \t%e\n", tot_smear, disp3[xc[XUP]],disp3[xc[YUP]],disp3[xc[ZUP]], t+1, (double)wils_loop1[r+base_r+nrmax*t]);
#else 
        wils_loop1[r+nrmax*t] += ( mu > XUP ? wils_loop1[r+base_r+nrmax*t] : 0 );
#endif
      }
    }
#ifdef AVERAGE_DIRECTIONS
    for(r=0;r<maxc[XUP];r++){
      node0_printf("WILS_LOOP1_%d  %d   %d \t%e\n", tot_smear, r, t, (double)wils_loop1[r+nrmax*t]/3.);
    }
#endif
  }

  free( wils_loop1);

} /* w_loop1 */

static void make_loops ( int t, 
                  Real *wils_loop1, 
		  su3_matrix *s_link, 
		  su3_matrix *s_link_f, 
		  su3_matrix *t_link_f ) {
  register int i;
  register site *s;
  su3_matrix tmat1,tmat2;

  FORALLSITES_OMP(i,s, private(tmat1,tmat2) reduction(+:wils_loop1[nrmax*t]) ){
    /* If the loop extends past t = nc[TUP] - 1 the temporal axial gauge link is nontrivial */
    if( (site_coord(s,xc[TUP])+t+1)>=nc[TUP] ){
      mult_su3_nn( &(s->link[xc[TUP]]), s_link_f+i, &tmat1);
      mult_su3_na( &tmat1, t_link_f+i, &tmat2);
      wils_loop1[nrmax*t] += realtrace_su3( &tmat2, s_link+i);
    }else{
      wils_loop1[nrmax*t] += realtrace_su3( s_link_f+i, s_link+i);
    }
  } END_LOOP_OMP;
} /* make_loops */

static int onaxis_loops ( Real * wils_loop1 ) {
  register int i,mu,r,t;
  register site *s=NULL;
  int base_r=0;
  msg_tag *mtag[3];
  su3_matrix *s_link=NULL, *s_link_f=NULL, *t_link_f=NULL;

  /* Allocate space for space-link product */
  s_link   = (su3_matrix *)malloc(sites_on_node*sizeof(su3_matrix));
  s_link_f = (su3_matrix *)malloc(sites_on_node*sizeof(su3_matrix));
  t_link_f = (su3_matrix *)malloc(sites_on_node*sizeof(su3_matrix));
  assert( ( s_link!=NULL && s_link_f!=NULL && t_link_f!=NULL ) );
  memset (s_link,0,sites_on_node*sizeof(su3_matrix));
  memset (s_link_f,0,sites_on_node*sizeof(su3_matrix));
  memset (t_link_f,0,sites_on_node*sizeof(su3_matrix));

  for( mu=XUP, base_r=0; mu<=ZUP; base_r+=maxc[mu],mu++){

    FORALLSITES_OMP(i,s, default(shared) ){
      su3mat_copy( &(s->link[xc[mu]]), s_link+i);
      su3mat_copy( &(s->link[xc[TUP]]), t_link_f+i);
    } END_LOOP_OMP;

    /* Start gather of forward time-like links */
    mtag[2] = start_gather_field( (void *)t_link_f, sizeof(su3_matrix), xc[mu], EVENANDODD, gen_pt[2] );

    /* Recursively construct the space-like segments and compute the Wilson loops with that segment */
    for(r=0;r<maxc[mu];r++){

      if( r>0 ){
        /* Collect the space-like segment and extend it by one link */
        wait_gather( mtag[1]);
        FORALLSITES_OMP(i,s, default(shared) ){ su3mat_copy( (su3_matrix *)(gen_pt[1][i]), &(s->staple)); } END_LOOP_OMP;
        FORALLSITES_OMP(i,s, default(shared) ){ mult_su3_nn( &(s->link[xc[mu]]), &(s->staple), s_link+i); } END_LOOP_OMP;
      }

      /* Prepare the forward space-like segment for parallel transport in time */
      FORALLSITES_OMP(i,s, default(shared) ){ su3mat_copy( s_link+i, s_link_f+i); } END_LOOP_OMP;

      /* Start gather of forward space-like segments for next t */
      mtag[0] = start_gather_field( (void *)s_link_f, sizeof(su3_matrix), xc[TUP], EVENANDODD, gen_pt[0] );

      /* Concurrently gather space-like segment for next r, if still needed. */
      if( r==0 ){
        mtag[1] = start_gather_field( (void *)s_link, sizeof(su3_matrix), xc[mu], EVENANDODD, gen_pt[1] );
      }else{ 
        if( r<(maxc[mu]-1) ){
          restart_gather_field( (void *)s_link, sizeof(su3_matrix), xc[mu], EVENANDODD, gen_pt[1], mtag[1] );
        }else{
          cleanup_gather( mtag[1]);
        }
      }

      /* Collect forward time-like links. */
      wait_gather( mtag[2]);
      FORALLSITES_OMP(i,s, default(shared) ){ su3mat_copy( (su3_matrix *)(gen_pt[2][i]), &(s->staple)); } END_LOOP_OMP;
      FORALLSITES_OMP(i,s, default(shared) ){ su3mat_copy( &(s->staple), t_link_f+i); } END_LOOP_OMP;

      /* Recursively compute the Wilson loops of different time extent */
      for(t=0;t<nct;t++){

        /* Collect forward space-like segments */
        wait_gather( mtag[0]);
        FORALLSITES_OMP(i,s, default(shared) ){ su3mat_copy( (su3_matrix *)(gen_pt[0][i]), &(s->staple)); } END_LOOP_OMP;
        FORALLSITES_OMP(i,s, default(shared) ){ su3mat_copy( &(s->staple), s_link_f+i); } END_LOOP_OMP;

        /* Start gather for next t, if still needed. */
        if( t<(nct-1) ){
          restart_gather_field( (void *)s_link_f, sizeof(su3_matrix), xc[TUP], EVENANDODD, gen_pt[0], mtag[0] );
        } else{
          cleanup_gather( mtag[0]);
        }

        /* Finally, compute the Wilson loops. */
        make_loops (t,(wils_loop1+r+base_r), s_link, s_link_f, t_link_f );

      } /* end loop over t */

      /* Start gather of forward time-like links for next r. */
      if( r<(maxc[mu]-1) ){
        restart_gather_field( (void *)t_link_f, sizeof(su3_matrix), xc[mu], EVENANDODD, gen_pt[2], mtag[2] );
      }else{
        cleanup_gather( mtag[2]);
      }
    } /* end loop over r */
  } /* end loop over  mu */

  free(s_link);
  free(s_link_f);
  free(t_link_f);

  return base_r;
} /* onaxis_loops */

