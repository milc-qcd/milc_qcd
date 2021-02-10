/****** hqp_helpers.c  -- ******************/

/***********************************************************************
* Auxiliary functions for heavy quark potential
* that are shared between Coulomb gauge correlators 
* and Wilson loops
* MIMD version 7
*
************************************************************************/

#include "../generic/generic_includes.h"	/* definitions files and prototypes */
#include "./hvy_qpot_includes.h"
#include "../include/openmp_defs.h"	
#include <assert.h>


double *hqp_alloc_dble_buffer( int mysize ) {
  double *buf = (double *)malloc(mysize*sizeof(double));
  assert( (buf!=NULL) );
  memset (buf,'\0',mysize*sizeof(double));
  return( buf );
}

su3_matrix *hqp_alloc_su3mat_buffer( int mysize ) {
  su3_matrix *buf = (su3_matrix *)malloc(mysize*sites_on_node*sizeof(su3_matrix));
  assert( (buf!=NULL) );
  memset (buf,'\0',mysize*sites_on_node*sizeof(su3_matrix));
  return( buf );
}

int hqp_disp_rsq_ok ( int disp[], hqp_geom *hqp ) {

  register int mu;
  int dok=1;
#ifdef ANISOTROPY
  Real rsq=0.;
#else
  int rsq=0;
#endif
  for (  mu=XUP; mu<=ZUP; mu++) 
    if ( disp[mu] <= maxc[mu] && disp[mu] >= 0 ) { 
      /* add up to displacement radius squared */
      rsq +=
#ifdef ANISOTROPY
#ifdef ONEDIM_ANISO_TEST
              ( mu==XUP && cor_dir != ani_dir ? ani_xiq*ani_xiq : iso_xiq*iso_xiq ) * 
#else
              ( mu==XUP && cor_dir != ani_dir ? ani_xiq*ani_xiq : 1 ) * 
#endif
#endif
              disp[mu]*disp[mu];
    } else 
      dok = 0; 

  if ( rsq > hqp->max_r2 )
    dok = 0;

  return( dok );
}

void hqp_free_dble_buffer( double *buf ) {

  free(buf);
  buf=NULL;
}

void hqp_free_su3mat_buffer( su3_matrix *buf ) {

  free(buf);
  buf=NULL;
}

void hqp_free_geometry( hqp_geom *hqp ) {
  free( hqp );
  hqp = NULL;
}

hqp_geom *hqp_geometry( char myname[] ) { 

  /*******************************************************
   *
   * Start of geometry setup and report (all algorithms):
   *
   * Text output concerning the geometry and its consistency.
   * In particular, the correlation directions 'xc[mu]' are
   * permuted from the original lattice directions 'mu'
   * to permit static quark correlations with arbitrary 
   * correlation time direction and arbitrary anisotropic 
   * direction.
   *
   ******************************************************/

  int mu;
  hqp_geom *hqp = malloc(sizeof(hqp_geom));;

#ifdef AVERAGE_DIRECTIONS
#ifdef ANISOTROPY
  assert( (ani_dir==xc[TUP]) );
#endif
  assert( (maxc[XUP]==maxc[YUP] && maxc[XUP] == maxc[ZUP] ) );
#endif

  node0_printf("%s with %d times smeared links\n",myname,tot_smear);
#if (defined HQPTIME || defined DEBUG)
  node0_printf(\
" [min,max](CT)[%d] = [%d,%d], \
\n max(C0)[%x] = %d, max(C1)[%x] = %d, max(C2)[%x] = %d, \
\n max_r = %.3g\n",cor_dir, min_ct,maxc[TUP], 
xc[XUP],maxc[XUP] ,xc[YUP],maxc[YUP] ,xc[ZUP],maxc[ZUP] ,max_r);
#endif

  assert ( maxc[XUP] <= nc[XUP]/2 && maxc[YUP] <= nc[YUP]/2 && maxc[ZUP] <= nc[ZUP]/2 ); 
  assert ( min_ct <= maxc[TUP] && maxc[TUP] <= nc[TUP] );
#ifdef PLOOPCOR_MEAS
  hqp->maxlen = nc[TUP];
//  assert ( maxc[TUP] == nc[TUP] );
#else 
  hqp->maxlen = maxc[TUP];
#endif

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
    if ( max_r < maxc[mu] )  
      maxc[mu] = floor(max_r); 
#else
#  ifdef ONEDIM_ANISO_TEST
    if ( max_r < ( mu == ani_dir ? ani_xiq : iso_xiq ) * maxc[mu] )  
      maxc[mu] = floor(max_r*( mu == ani_dir ? 1./ani_xiq : 
      (iso_xiq > 0 ? 1./iso_xiq : 0 ))); 
#  else
    if ( max_r < ( mu == ani_dir ? ani_xiq : 1 ) * maxc[mu] )  
      maxc[mu] = floor(max_r*( mu == ani_dir ? 1./ani_xiq : 1 )); 
#  endif
#endif
  }
  hqp->max_r2=floor(max_r*max_r);

  local_lattice_size(hqp->llat);

#ifdef HQPTIME
  int nll[3];
  nll[XUP]=nc[XUP]/hqp->llat[xc[XUP]]; 
  nll[YUP]=nc[YUP]/hqp->llat[xc[YUP]]; 
  nll[ZUP]=nc[ZUP]/hqp->llat[xc[ZUP]];
  if ( numnodes() > 1 ) {
    node0_printf("####################################\
                \n Using a distributed grid of %dx%dx%d(x%d) nodes\n",
                 nll[XUP],nll[YUP],nll[ZUP],nc[TUP]/hqp->llat[xc[TUP]]);
  } 
  node0_printf("####################################\n\n");
#endif
fflush(stdout);

  return (hqp);
}

void hqp_output_corr(char cname[], char smtag[], int disp[], double corr) {

#ifndef AVERAGE_DIRECTIONS
  node0_printf("%s%s:  %d  %d  %d   %d \t%.6e\n",
    cname,smtag, disp[XUP], disp[YUP], disp[ZUP], disp[TUP], 
    corr /(double)(volume)
  );
#else
  node0_printf("%s%s  %d   %d \t%e\n",
    cname,smtag, disp[XUP], disp[TUP], 
    corr /(double)(volume)
  );
#endif
}

/* this routine returns the local lattice size */
void local_lattice_size (int *llat) {
  int i;
  for (i=XUP;i<=TUP;i++) { llat[i]=locx[i]; }
}

#ifdef WLOOP_MEAS
int loop_rt ( hqp_geom *hqp, 
              int mysize, 
              int mi,
              int disp[], 
              double *wils_loop
            ) {
  register int i,t,r;
  register site *s;
  msg_tag *mtag[3], *gmtag;

  su3_matrix *s_link   = hqp_alloc_su3mat_buffer( 1 );
  su3_matrix *s_link_f = hqp_alloc_su3mat_buffer( 1 );
  su3_matrix *t_link_f = hqp_alloc_su3mat_buffer( 1 );

  /* Indicating that the displacement is not onaxis */
  int mu=NODIR;
  /* Three distinct cases for onaxis loops */
  if ( disp[xc[YUP]] == 0 && disp[xc[ZUP]] == 0 ) 
    mu=XUP;
  if ( disp[xc[ZUP]] == 0 && disp[xc[XUP]] == 0 ) 
    mu=YUP;
  if ( disp[xc[XUP]] == 0 && disp[xc[YUP]] == 0 ) 
    mu=ZUP;

  for(r=0;r<mysize;r++){

    if( r==0 ){

      /* Put time links into buffers */
      FORALLSITES_OMP(i,s, default(shared) ){ su3mat_copy( &(s->link[xc[TUP]]), t_link_f+i); } END_LOOP_OMP;

      /* Put spatial links into buffers, initialize onaxis gather to extend */
      if (mu>NODIR) {
        FORALLSITES_OMP(i,s, default(shared) ){ su3mat_copy( &(s->link[xc[mu]]), s_link+i); } END_LOOP_OMP;
        mtag[1] = start_gather_field( (void *)s_link, sizeof(su3_matrix), xc[mu], EVENANDODD, gen_pt[1] );

      } else 
        FORALLSITES_OMP(i,s, default(shared) ){ su3mat_copy( &(s->diag), s_link+i); } END_LOOP_OMP;


      /* Initiate (offaxis general) gather of time-like links */
      if (mu>NODIR) {
        mtag[2] = start_gather_field( (void *)t_link_f, sizeof(su3_matrix), xc[mu], EVENANDODD, gen_pt[2] );

      } else 
        gmtag = start_general_gather_field( (void *)t_link_f, sizeof(su3_matrix), disp, EVENANDODD, gen_pt[4] );


    }else{

      if (mu>NODIR) {

        /* Collect the space-like line and extend by one step */
        wait_gather( mtag[1]);
        FORALLSITES_OMP(i,s, default(shared) ){ su3mat_copy( (su3_matrix *)(gen_pt[1][i]), &(s->staple)); } END_LOOP_OMP;
        FORALLSITES_OMP(i,s, default(shared) ){ mult_su3_nn( &(s->link[xc[mu]]), &(s->staple), s_link+i); } END_LOOP_OMP;

        /* Reinitiate next gather of spatial lines along time */
        if( r < (mysize-1) ) {
          restart_gather_field( (void *)s_link, sizeof(su3_matrix), xc[mu], EVENANDODD, gen_pt[1], mtag[1] );

        } else
          cleanup_gather( mtag[1]);

      } else {
        /* Collect the space-like line */ 
        wait_general_gather( gmtag);
        FORALLSITES_OMP(i,s, default(shared) ){ su3mat_copy( (su3_matrix *)(gen_pt[4][i]), &(s->staple)); } END_LOOP_OMP;

        /* Reinitiate inbetween gather of time-like links */
        cleanup_general_gather( gmtag);
        gmtag = start_general_gather_field( (void *)t_link_f, sizeof(su3_matrix), disp, EVENANDODD,gen_pt[4]);

        /* Extend space-like line by one step */
        FORALLSITES_OMP(i,s, default(shared) ){ mult_su3_nn( &(s->diag), &(s->staple), s_link+i); } END_LOOP_OMP;

      }
    }

    /* Prepare the forward space-like line for parallel transport in time */
    memcpy(s_link_f,s_link,sites_on_node*sizeof(su3_matrix));
    /* Start time gather of forward space-like line */
    mtag[0] = start_gather_field( (void *)s_link_f, sizeof(su3_matrix), xc[TUP], EVENANDODD, gen_pt[0] );


    /* Collect forward time-like links. */
    if (mu>NODIR) {
      wait_gather( mtag[2]);
      FORALLSITES_OMP(i,s, default(shared) ){ su3mat_copy( (su3_matrix *)(gen_pt[2][i]), &(s->staple)); } END_LOOP_OMP;

    } else {
      wait_general_gather( gmtag);
      FORALLSITES_OMP(i,s, default(shared) ){ su3mat_copy( (su3_matrix *)(gen_pt[4][i]), &(s->staple)); } END_LOOP_OMP;

      /* Inbetween gather space-links across the diagonal for next r */
      cleanup_general_gather( gmtag);
      if( r<(mysize-1) ){
        gmtag = start_general_gather_field( (void *)s_link, sizeof(su3_matrix), disp, EVENANDODD, gen_pt[4] );
      }

    }

    /* Put time-like links into buffer */
    FORALLSITES_OMP(i,s, default(shared) ){ su3mat_copy( &(s->staple), t_link_f+i); } END_LOOP_OMP;

    /* Recursively compute the Wilson loops of different time extent */
    for(t=0;t<hqp->maxlen;t++){

      /* Collect forward space-like line */
      wait_gather( mtag[0]);
      FORALLSITES_OMP(i,s, default(shared) ){ su3mat_copy( (su3_matrix *)(gen_pt[0][i]), &(s->staple)); } END_LOOP_OMP;
      FORALLSITES_OMP(i,s, default(shared) ){ su3mat_copy( &(s->staple), s_link_f+i); } END_LOOP_OMP;

      /* Start gather for next t, if still needed. */
      if( t<(hqp->maxlen-1) ){
        restart_gather_field( (void *)s_link_f, sizeof(su3_matrix), xc[TUP], EVENANDODD, gen_pt[0], mtag[0] );
      }else{
        cleanup_gather( mtag[0]);
      }

      /* Finally, compute the Wilson loops. */
      make_loops (t,mi,(wils_loop+r), s_link,s_link_f,t_link_f);

    } /* end loop over t */
  } /* end loop over r */

  hqp_free_su3mat_buffer( s_link   );
  hqp_free_su3mat_buffer( s_link_f );
  hqp_free_su3mat_buffer( t_link_f );

  return (mysize);
} /* loop_rt */


void make_loops ( int t, int mysize, 
                  double *wils_loop, 
		  su3_matrix *s_link, 
		  su3_matrix *s_link_f, 
		  su3_matrix *t_link_f ) {
  register int i;
  register site *s;
  su3_matrix tmat1,tmat2;

  FORALLSITES_OMP(i,s, private(tmat1,tmat2) reduction(+:wils_loop[mysize*t]) ){
    // If the loop extends past t = nc[TUP] - 1 the temporal axial gauge link is nontrivial 
    if( (site_coord(s,xc[TUP])+t+1)>=nc[TUP] ){
      mult_su3_nn( &(s->link[xc[TUP]]), s_link_f+i, &tmat1);
      mult_su3_na( &tmat1, t_link_f+i, &tmat2);
      wils_loop[mysize*t] += (double)realtrace_su3( &tmat2, s_link+i);
    }else{
      wils_loop[mysize*t] += (double)realtrace_su3( s_link_f+i, s_link+i);
    }
  } END_LOOP_OMP;
} /* make_loops */
#endif

