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


static hqp_geom *geom=NULL;
static double *wils_loop1=NULL;

static su3_matrix *s_link=NULL, *s_link_f=NULL, *t_link_f=NULL;

static inline void global_sums( int mysize );
static inline int hqp_buffer_dim (void );
static void hqp_free_buffers( void );
static void hqp_setup_buffers( int mysize );
static inline void output_all( char smtag[], int mi );

static int onaxis_loops( double *wils_loop1, int mysize );

void w_loop1(int tot_smear) {
  int base_r=0;
  char myname[] = "new_w_loop1";

  char smtag[MAXSMTAG]="";
#ifdef SMEARING
  sprintf(smtag,"_%d",tot_smear);
#endif

  geom = hqp_geometry( myname );
  int mi = hqp_buffer_dim();
  hqp_setup_buffers( geom->maxlen*mi );

  /****************************************************************
  * 
   * Proceed to loop over directions and recursively construct the 
   * space-like segments and compute the Wilson loops with these 
   *
   ***************************************************************/

  for( int mu=XUP; mu<=ZUP; mu++) {
    int disp[4]={0,0,0,0}; 
    disp[mu] = 1;
    base_r += loop_rt( geom, maxc[mu], mi, disp, wils_loop1+base_r );
  }
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

  hqp_free_buffers( );
  hqp_free_geometry( geom );

} /* w_loop1 */


static inline void global_sums( int mysize ) {

  g_vecdoublesum(wils_loop1, mysize );
}

static inline int hqp_buffer_dim (void ) {

  return ((maxc[XUP])+(maxc[YUP])+(maxc[ZUP]));
}

static void hqp_setup_buffers( int mysize ) { 

  wils_loop1 = hqp_alloc_dble_buffer( mysize );
}

static void hqp_free_buffers( void ) { 

  hqp_free_dble_buffer( wils_loop1 );
}

static inline void output_all( char smtag[], int mi ) {

  int mu, base_r;
  int disp[4]={0,0,0,0};
  int *r;
  for(disp[TUP]=1;disp[TUP]<=geom->maxlen;disp[TUP]++) { 
    for( mu=XUP, base_r=0; mu<=ZUP; base_r+=maxc[mu],mu++) {
      memset(disp, 0,3*sizeof(int));
      r=&(disp[mu]);
      for((*r)=1;(*r)<=maxc[mu];(*r)++){
        if ( hqp_disp_rsq_ok( disp, geom ) == 1 ) 
#ifndef AVERAGE_DIRECTIONS
          hqp_output_corr("WILS_LOOP1",smtag, disp, wils_loop1[(((*r)-1)+base_r+mi*(disp[TUP]-1))]/3.);
#else 
          wils_loop1[((*r)-1)+mi*(disp[TUP]-1)]+=(mu>XUP?1:-2)* 
            wils_loop1[((*r)-1)+base_r+mi*(disp[TUP]-1)]/3. ;
#endif
      }
    }
#ifdef AVERAGE_DIRECTIONS
    memset(disp, 0,3*sizeof(int));
    disp[TUP]--;
    r=&(disp[XUP]);
    for((*r)=0;*r<maxc[XUP];(*r)++)
      if ( hqp_disp_rsq_ok( disp, geom ) == 1 ) 
        hqp_output_corr("WILS_LOOP1",smtag, disp, wils_loop1[((*r)+mi*(disp[TUP]))]/3.);
    disp[TUP]++;
#endif
  }

}
