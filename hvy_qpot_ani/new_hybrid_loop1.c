/********************** hybrid_loop1.c *******************************/
/* MIMD version 7 */
/* Computes on-axis Wilson loops appropriate for excited flux-tube
   potentials (including those relevant to hybrids) with symmetry
   \Pi_u, and Delta_g 

   ALSO computes the conventional Wilson loop with symmetry \Sigma_g^+

   For the diatomic molecular notation applied to this problem, see
   K.J. Juge, J. Kuti, and C.J. Morningstar, Lattice '97 Nucl Phys
   (Proc Supp) B 63, 326 (1998)

   In terms of representations of D_4h these are
   \Sigma_g^+ = A_1g
   \Pi_u = E_u
   \Delta_g = B_1g

 */

/* We generate the hybrid loops by doing parallel displacements of the
   linear space-like Wilson line segment.  The displacement is in the
   plane perpendicular to the space-like segment as follows (Q denotes
   the static quark and o the locus of the Wilson line):

   \Sigma_g^+    Qo    (no displacement - conventional Wilson loop)

                                          o
                                          |
                                          |
   \Pi_u         Q----o - o----Q    or    Q - Q
                                              |
                                              |
                                              o

                                    o
                                    |
                                    |
   \Delta_g      Q----o + o----Q -  Q - Q
                                        |
                                        |
                                        o

*/


/* Assumes corr-temporal axial gauge with the nontrivial link copied to
   all corr-time links as in hvy_qpot_ani/new_ax_gauge.c */

/* based on w_loop1 by UMH */

/* 4/16/02 CD */

#define MAX_PATH_LENGTH 4

typedef struct
{
  char name[12];                  /* name of path */
  int length;                     /* length of path */
  int disp[4];                    /* net displacement of path */
  int dir[MAX_PATH_LENGTH];       /* list of directions in path */
} link_path;

#include "../generic/generic_includes.h"        /* definitions files and prototypes */
#include "./hvy_qpot_includes.h"
#include "../include/openmp_defs.h"
#include <assert.h>


#define WILS_LOOP1(i,t,r)   wils_loop1[ (i)*nrmax*nct + (t)*nrmax + r ]
#define TRANS_LINKS(i,j)    (trans_links + (i)*NTRANS + (j))
#define TRANS_LINKS_F(i,j)  (trans_links_f + (i)*NTRANS_F + (j))
/* flux_links must be laid out with site index varying most rapidly since
   we want to gather j's separately */
#define FLUX_LINKS(i,j)     (flux_links + (i) + (j)*sites_on_node)
#define FLUX_LINKS_F(i,j)   (flux_links_f + (i)*NFLUX_F + (j))

/* Names of messages used to index "mtag" and "gen_pt" */
/* NMSGS must not exceed N_POINTERS */
enum{ M_S_LINK, M_F_LINKS_F, M_T_LINKS_F, M_STAP_POS1, 
      M_STAP_NEG1, M_STAP_POS2, M_STAP_NEG2, NMSGS };

/* Names of flux tube shapes built from transverse links 
   and on-axis links used to index "flux_links" */
enum{ S_LINK, STAP_POS1, STAP_NEG1, STAP_POS2, STAP_NEG2, NFLUX };

/* Names of transverse links 
   used to index "trans_path", "trans_links" and "trans_links_f" */
enum{ XX, YY, ZZ, NTRANS };

/* Names of space-shifted transverse links */
/* used to index "trans_links_f" */
enum{ TRANS_PATH1_F, TRANS_PATH2_F, T_LINK_F, NTRANS_F };

/* Names of flux tube shapes transported forward in time */
/* Used to index "flux_links_f" */
enum{ S_LINK_F, STAP_POS1_F, NFLUX_F };

/* Names of loop observables */
/* Used to index "wils_loop1" */
enum{ W_LOOP1, STAP_SIG_GP1, STAP_PI_U1, STAP_DELTA_G1, NWLOOP1 };

static void make_loops ( int t,
                  Real *wils_loop,
                  su3_matrix *flux_links,
                  su3_matrix *flux_links_f,
                  su3_matrix *trans_links, 
                  su3_matrix *trans_links_f );

static int hybrid_loops( link_path trans_path[], Real *wils_loop1 );

static int nrmax,nct;

void hybrid_loop1(int tot_smear) {

  register int i,j, mu,r,t;
  int mu1=0,mu2=0,trans_path1=0,trans_path2=0;
  int **dispc=NULL;
  register site *s;
  double max_r2;
  int  base_r=0;
  int disp3[4];
  Real *wils_loop1=NULL;
  char myname[] = "new_hybrid_loop1";

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

  msg_tag *mtag[NMSGS];

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
  /* Check the rules */
  if(NMSGS > N_POINTERS){
    if(this_node == 0)fprintf(stderr,"%s: Aborted. gen_pt array too small.",myname);
    terminate(1);
  }

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
  
  /* Paths for transverse links */
  dispc=(int **)malloc(3*sizeof(int*));
  assert(dispc!=NULL);
  for ( mu=XUP; mu<=ZUP; mu++ ) {
    dispc[mu] = (int *)malloc(4*sizeof(int));
    assert(dispc[mu]!=NULL);
    memset (dispc[mu],0,4*sizeof(int));
    dispc[mu][xc[mu]] = 2;
  }
  link_path trans_path[NTRANS];
  sprintf(trans_path[XUP].name ,"%s", "XX" );
  sprintf(trans_path[YUP].name ,"%s", "YY" );
  sprintf(trans_path[ZUP].name ,"%s", "ZZ" );
  for ( mu=XUP; mu<=ZUP; mu++ ) {
    trans_path[mu].length = 2;
    for ( mu1=XUP; mu1<=TUP; mu1++ ) {
      trans_path[mu].disp[xc[mu1]] = dispc[mu][xc[mu1]];
    }
    for ( mu1=0; mu1<MAX_PATH_LENGTH; mu1++) {
      trans_path[mu].dir[mu1] = ( mu1 <2 ? xc[mu] : NODIR );
    }
  }
  for ( mu=XUP; mu<=ZUP; mu++ ) { free(dispc[mu]); }
  free(dispc);

  wils_loop1 = (Real *)malloc(nct*nrmax*sizeof(Real)*NWLOOP1);
  assert( (wils_loop1!=NULL) );
  memset (wils_loop1,0,nct*nrmax*sizeof(Real)*NWLOOP1);

  /****************************************************************
   * 
   * End of setup 
   * Proceed to loop over directions and 
   * recursively construct the space-like segments and compute
   * the Wilson loops with that segment 
   *
   ***************************************************************/

  base_r += hybrid_loops( trans_path, wils_loop1+base_r );
  assert( (base_r==nrmax) );
  
  /****************************************************************
   * 
   * End of calculation 
   * Proceed to normalization and printing of the Wilson loops
   * Remark: do not average the directions (impossible if anisotropic),
   * instead print the direction into one column
   * can activate averaging of directions via compiler macro AVERAGE_DIRECTIONS
   *
   ***************************************************************/

  g_vecfloatsum(wils_loop1,nct*nrmax*NWLOOP1);
  for(t=0;t<nct;t++) {
    for( mu=XUP, base_r=0; mu<=ZUP; base_r+=maxc[mu],mu++) {
      for(r=0;r<maxc[mu];r++) {
        WILS_LOOP1(W_LOOP1,t,r+base_r) = ((double)WILS_LOOP1(W_LOOP1,t,r+base_r)/(double)(3*volume));
#ifndef AVERAGE_DIRECTIONS
        memset(disp3, 0,4*sizeof(int));
        disp3[xc[(mu)%3]]=r+1;
        node0_printf("WILS_LOOP1_%d   %d  %d  %d   %d \t%e\n", tot_smear, disp3[xc[XUP]],disp3[xc[YUP]],disp3[xc[ZUP]], t+1,(double)WILS_LOOP1(W_LOOP1,t,r+base_r));
#else
        WILS_LOOP1(W_LOOP1,t,r) += ( base_r > 0 ? WILS_LOOP1(W_LOOP1,t,r+base_r) : 0 );
#endif
      }
    }
#ifdef AVERAGE_DIRECTIONS
    for(r=0;r<maxc[XUP];r++){
      node0_printf("WILS_LOOP1_%d  %d   %d \t%e\n", tot_smear, r, t,  (double)WILS_LOOP1(W_LOOP1,t,r)/3. );
    }
#endif
  }
  
  for(t=0;t<nct;t++) {
    for( mu=XUP, base_r=0; mu<=ZUP; base_r+=maxc[mu],mu++) {
      for(r=0;r<maxc[mu];r++) {
        WILS_LOOP1(STAP_SIG_GP1,t,r+base_r) = ((double)WILS_LOOP1(STAP_SIG_GP1,t,r+base_r)/(double)(3*volume));
#ifndef AVERAGE_DIRECTIONS
        memset(disp3, 0,4*sizeof(int));
        disp3[xc[(mu)%3]]=r+1;
        node0_printf("STAP_SIG_GP1_%d   %d  %d  %d   %d \t%e\n", tot_smear, disp3[xc[XUP]],disp3[xc[YUP]],disp3[xc[ZUP]], t+1,(double)WILS_LOOP1(STAP_SIG_GP1,t,r+base_r));
#else
        WILS_LOOP1(STAP_SIG_GP1,t,r) += ( base_r > 0 ? WILS_LOOP1(STAP_SIG_GP1,t,r+base_r) : 0 );
#endif
      }
    }
#ifdef AVERAGE_DIRECTIONS
    for(r=0;r<maxc[XUP];r++){
      node0_printf("STAP_SIG_GP1_%d  %d   %d \t%e\n", tot_smear, r, t,  (double)WILS_LOOP1(STAP_SIG_GP1,t,r)/3. );
    }
#endif
  }
    
  for(t=0;t<nct;t++) {
    for( mu=XUP, base_r=0; mu<=ZUP; base_r+=maxc[mu],mu++) {
      for(r=0;r<maxc[mu];r++) {
        WILS_LOOP1(STAP_PI_U1,t,r+base_r) = ((double)WILS_LOOP1(STAP_PI_U1,t,r+base_r)/(double)(3*volume));
#ifndef AVERAGE_DIRECTIONS
        memset(disp3, 0,4*sizeof(int));
        disp3[xc[(mu)%3]]=r+1;
        node0_printf("STAP_PI_U1_%d   %d  %d  %d   %d \t%e\n", tot_smear, disp3[xc[XUP]],disp3[xc[YUP]],disp3[xc[ZUP]], t,(double)WILS_LOOP1(STAP_PI_U1,t,r+base_r));
#else
        WILS_LOOP1(STAP_PI_U1,t,r) += ( base_r > 0 ? WILS_LOOP1(STAP_PI_U1,t+1,r+base_r) : 0 );
#endif
      }
    }
#ifdef AVERAGE_DIRECTIONS
    for(r=0;r<maxc[XUP];r++){
      node0_printf("STAP_PI_U1_%d  %d   %d \t%e\n", tot_smear, r, t,  (double)WILS_LOOP1(STAP_PI_U1,t,r)/3. );
    }
#endif
  }
    
  for(t=0;t<nct;t++) {
    for( mu=XUP, base_r=0; mu<=ZUP; base_r+=maxc[mu],mu++) {
      for(r=0;r<maxc[mu];r++) {
        WILS_LOOP1(STAP_DELTA_G1,t,r+base_r) = ((double)WILS_LOOP1(STAP_DELTA_G1,t,r+base_r)/(double)(3*volume));
#ifndef AVERAGE_DIRECTIONS
        memset(disp3, 0,4*sizeof(int));
        disp3[xc[(mu)%3]]=r+1;
        node0_printf("STAP_DELTA_G1_%d   %d  %d  %d   %d \t%e\n", tot_smear, disp3[xc[XUP]],disp3[xc[YUP]],disp3[xc[ZUP]], t+1,(double)WILS_LOOP1(STAP_DELTA_G1,t,r+base_r));
#else
        WILS_LOOP1(STAP_DELTA_G1,t,r) += ( base_r > 0 ? WILS_LOOP1(STAP_DELTA_G1,t,r+base_r) : 0 );
#endif
      }
    }
#ifdef AVERAGE_DIRECTIONS
    for(r=0;r<maxc[XUP];r++){
      node0_printf("STAP_DELTA_G1_%d  %d   %d \t%e\n", tot_smear, r, t,  (double)WILS_LOOP1(STAP_DELTA_G1,t,r)/3. );
    }
#endif
  }
  
  free( wils_loop1 );
  
} /* hybrid_loop1 */

static void make_loops ( int t,
                  Real *wils_loop1,
                  su3_matrix *flux_links,
                  su3_matrix *flux_links_f,
                  su3_matrix *trans_links, 
                  su3_matrix *trans_links_f ) {
  register int i;
  register site *s;
  su3_matrix tmat1,tmat2;
  su3_matrix *tmatp;
  Real wl=0., sig=0., pi=0., delta=0.;

  /* Compute naive Wilson loop term */
  FORALLSITES_OMP(i,s, private(tmat1,tmat2,tmatp) reduction(+:wl) ){
    /* If the loop extends past t = nc[TUP] - 1 the temporal axial gauge link is nontrivial */
    if( (site_coord(s,xc[TUP])+t+1)>=nc[TUP] ){
      mult_su3_nn( &(s->link[xc[TUP]]), FLUX_LINKS_F(i,S_LINK_F), &tmat1);
      mult_su3_na( &tmat1, TRANS_LINKS_F(i,T_LINK_F), &tmat2);
      tmatp = &tmat2;
    }else
      tmatp = FLUX_LINKS_F(i,S_LINK_F);
    
    wl += realtrace_su3( tmatp, FLUX_LINKS(i,S_LINK) );
  } END_LOOP_OMP;
  WILS_LOOP1(W_LOOP1,t,0) = wl;

  /* Compute big staple terms */
  
  /* Do projection for source and use  mu1 for sink */
  FORALLSITES_OMP(i,s, private(tmat1,tmat2,tmatp) reduction(+:sig,pi,delta) ){
    /* If the loop extends past t = nc[TUP] - 1 the temporal axial gauge link is nontrivial */
    if( (site_coord(s,xc[TUP])+t+1)>=nc[TUP] ){
      mult_su3_nn( &(s->link[xc[TUP]]), FLUX_LINKS_F(i,STAP_POS1_F), &tmat1);
      mult_su3_na( &tmat1, TRANS_LINKS_F(i,T_LINK_F), &tmat2);
      tmatp = &tmat2;  /* Use tmat2 to close loop */
    }else
      /* Use stap_pos1_f to close loop */
      tmatp = FLUX_LINKS_F(i,STAP_POS1_F);
	    
    /* Projection of source onto irreps */

    sig +=
      realtrace_su3( tmatp, FLUX_LINKS(i,STAP_POS1) ) +
      realtrace_su3( tmatp, FLUX_LINKS(i,STAP_NEG1) ) +
      realtrace_su3( tmatp, FLUX_LINKS(i,STAP_POS2) ) +
      realtrace_su3( tmatp, FLUX_LINKS(i,STAP_NEG2) );

    pi +=
      realtrace_su3( tmatp, FLUX_LINKS(i,STAP_POS1) ) -
      realtrace_su3( tmatp, FLUX_LINKS(i,STAP_NEG1) );


    delta +=
      realtrace_su3( tmatp, FLUX_LINKS(i,STAP_POS1) ) +
      realtrace_su3( tmatp, FLUX_LINKS(i,STAP_NEG1) ) -
      realtrace_su3( tmatp, FLUX_LINKS(i,STAP_POS2) ) -
      realtrace_su3( tmatp, FLUX_LINKS(i,STAP_NEG2) );
  } END_LOOP_OMP;
  WILS_LOOP1(STAP_SIG_GP1,t,0) = sig;
  WILS_LOOP1(STAP_PI_U1,t,0) = pi;
  WILS_LOOP1(STAP_DELTA_G1,t,0) = delta;

} /* make_loops */

static int hybrid_loops( link_path trans_path[], Real *wils_loop1 ) {
  register int i,j,mu,r,t;
  register site *s=NULL;
  int mu1=0,mu2=0,trans_path1=0,trans_path2=0;
  int base_r=0;
  int disp[4];
  su3_matrix  *flux_links=NULL,  *flux_links_f=NULL;
  su3_matrix *trans_links=NULL, *trans_links_f=NULL;
  su3_matrix *tempmat1=NULL;
  su3_matrix tmat1;
  msg_tag *mtag[NMSGS];
  char myname[] = "hybrid_loops";

  /* Allocate space for transverse link products */
  trans_links   = (su3_matrix *)malloc(NTRANS  *sites_on_node*sizeof(su3_matrix));
  /* Allocate space for shifted auxiliary link products */
  trans_links_f = (su3_matrix *)malloc(NTRANS_F*sites_on_node*sizeof(su3_matrix));
  assert( ( trans_links!=NULL && trans_links_f!=NULL ) );

  tempmat1 = (su3_matrix *)malloc(sites_on_node*sizeof(su3_matrix));
  assert(tempmat1!=NULL);
  /* Compute and store products of transverse links */
  /* trans_links[i][j] is set to the link product for the
     jth path type that ends at site i */
  for(j = 0; j < NTRANS; j++){
    path_product(trans_path[j].dir, trans_path[j].length, tempmat1);
    FORALLSITES_OMP(i,s, default(shared) ){ su3mat_copy( &(tempmat1[i]), TRANS_LINKS(i,j)); } END_LOOP_OMP;
  }
  free(tempmat1);

  /* Allocate space for timeward shifted flux tube shapes */
  flux_links_f = (su3_matrix *)malloc(NFLUX_F*sites_on_node*sizeof(su3_matrix));
  /* Allocate space for flux-tube shapes */
  flux_links   = (su3_matrix *)malloc(NFLUX  *sites_on_node*sizeof(su3_matrix));
  assert( ( flux_links!=NULL && flux_links_f!=NULL ) );

  for( mu=XUP, base_r=0; mu<=ZUP; base_r+=maxc[mu],mu++){
    
    /* Build giant staples */
    switch(mu){
    case XUP:
      mu1 = YUP; mu2 = ZUP; 
      trans_path1 = YY; trans_path2 = ZZ;
      break;
    case YUP:
      mu1 = ZUP; mu2 = XUP; 
      trans_path1 = ZZ; trans_path2 = XX;
      break;
    case ZUP:
      mu1 = XUP; mu2 = YUP; 
      trans_path1 = XX; trans_path2 = YY;
      break;
    default:
      if(this_node == 0)fprintf(stderr,"%s unknown direction %d\n", myname, mu);
      break;
    }
      
    /* Initialization for loop operators */
      
    /* Start building s_link */
    /* Prepare to shift trans_links */
    FORALLSITES_OMP(i,s, default(shared) ){
      su3mat_copy( &(s->link[xc[mu]]), FLUX_LINKS(i,S_LINK) );
      su3mat_copy(TRANS_LINKS(i,trans_path1), TRANS_LINKS_F(i,TRANS_PATH1_F));
      su3mat_copy(TRANS_LINKS(i,trans_path2), TRANS_LINKS_F(i,TRANS_PATH2_F));
      su3mat_copy( &(s->link[xc[TUP]]), TRANS_LINKS_F(i,T_LINK_F));
    } END_LOOP_OMP;
      
    /* Start spatial (mu) transport of all transverse links */
    /* gen_pt[M_T_LINKS_F][i] will point to trans_links_f[i+mu] for site i */
    mtag[M_T_LINKS_F] = start_gather_field( (void *)(trans_links_f), NTRANS_F*sizeof(su3_matrix), xc[mu], EVENANDODD, gen_pt[M_T_LINKS_F] );
      
    /* Recursively construct the space-like segments and compute the Wilson loops with that segment */
    for(r=0;r<maxc[mu];r++){

      if( r>0 ){
        /* Collect the space-like segment and extend it by one link in the direction mu. */
        /* s_link[i] then contains the product of links BEGINNING at site i */
        wait_gather( mtag[M_S_LINK]);
        FORALLSITES_OMP(i,s, default(shared) ){ su3mat_copy( (su3_matrix *)(gen_pt[M_S_LINK][i]), &(s->staple)); } END_LOOP_OMP;
        FORALLSITES_OMP(i,s, default(shared) ){ mult_su3_nn( &(s->link[xc[mu]]), &(s->staple), FLUX_LINKS(i,S_LINK) ); } END_LOOP_OMP;
      }

      /* Shift the space-like segment for next r, if still needed. */
      /* gen_pt[M_S_LINK][i] will point to s_link[i+mu] */
      if( r==0 ){
        mtag[M_S_LINK] = start_gather_field( (void *)FLUX_LINKS(0,S_LINK), sizeof(su3_matrix), xc[mu], EVENANDODD, gen_pt[M_S_LINK] );
      }else 
        if( r<(maxc[mu]-1) ){
          restart_gather_field( (void *)FLUX_LINKS(0,S_LINK), sizeof(su3_matrix), xc[mu], EVENANDODD, gen_pt[M_S_LINK], mtag[M_S_LINK] );
        }else{
          cleanup_gather( mtag[M_S_LINK]);
        }

      /* Collect the transverse links shifted from mu . */
      /* trans_links_f[i][j] <- trans_links_f[i+mu][j] */
      wait_gather( mtag[M_T_LINKS_F]);
      for(j = 0; j < NTRANS_F; j++){ 
        FORALLSITES_OMP(i,s, default(shared) ){ su3mat_copy( (su3_matrix *)(gen_pt[M_T_LINKS_F][i]) + j, &(s->staple)); } END_LOOP_OMP;
        FORALLSITES_OMP(i,s, default(shared) ){ su3mat_copy( &(s->staple), TRANS_LINKS_F(i,j) ); } END_LOOP_OMP;
      }

      /* Start making flux_links[i][STAP_NEG1] <- backward s_link at -disp */
      FORALLUPDIR(j){ disp[j] = -trans_path[ mu1].disp[j]; }
      mtag[M_STAP_NEG1] = start_general_gather_field( (void *)FLUX_LINKS(0,S_LINK), sizeof(su3_matrix), disp, EVENANDODD, gen_pt[M_STAP_NEG1] );

      /* flux_links[i][STAP_POS1] <- forward staple at disp */
      FORALLSITES_OMP(i,s, private(tmat1) ){
        mult_su3_nn( TRANS_LINKS(i,trans_path1), FLUX_LINKS(i,S_LINK), &tmat1 );
        mult_su3_na( &tmat1, TRANS_LINKS_F(i,TRANS_PATH1_F), FLUX_LINKS(i,STAP_POS1) ); 
      } END_LOOP_OMP;

      /* Can't overlap gen_gathers!! */
      wait_general_gather(mtag[M_STAP_NEG1]);

      /* Start making flux_links[i][STAP_NEG2] <- backward s_link at -disp */
      FORALLUPDIR(j){disp[j] = -trans_path[ mu2].disp[j];}
      mtag[M_STAP_NEG2] = start_general_gather_field( (void *)FLUX_LINKS(0,S_LINK), sizeof(su3_matrix), disp, EVENANDODD, gen_pt[M_STAP_NEG2] );

      /* flux_links[i][STAP_POS2] <- forward staple at disp */
      FORALLSITES_OMP(i,s, private(tmat1) ){
        mult_su3_nn( TRANS_LINKS(i,trans_path2), FLUX_LINKS(i,S_LINK), &tmat1 );
        mult_su3_na( &tmat1, TRANS_LINKS_F(i,TRANS_PATH2_F), FLUX_LINKS(i,STAP_POS2) );
      } END_LOOP_OMP;
	
      wait_general_gather(mtag[M_STAP_NEG2]);

      /* Shift staples to the sites where the t_links join them */
      FORALLUPDIR(j){disp[j] = trans_path[ mu1].disp[j];}
      mtag[M_STAP_POS1] = start_general_gather_field( (void *)FLUX_LINKS(0,STAP_POS1), sizeof(su3_matrix), disp, EVENANDODD, gen_pt[M_STAP_POS1] );

      /* Finish flux_links[i][STAP_NEG1] <- backward staple in dir1 */
      FORALLSITES_OMP(i,s, private(tmat1) ){
        mult_su3_an( TRANS_LINKS(i,trans_path1), (su3_matrix *)gen_pt[M_STAP_NEG1][i], &tmat1 );
        mult_su3_nn( &tmat1, TRANS_LINKS_F(i,TRANS_PATH1_F), FLUX_LINKS(i,STAP_NEG1) );
      } END_LOOP_OMP;

      wait_general_gather(mtag[M_STAP_POS1]);

      /* Collect results of the shift of forward staples */
      FORALLSITES_OMP(i,s, default(shared) ){ su3mat_copy( (su3_matrix *)gen_pt[M_STAP_POS1][i], &(s->staple)); } END_LOOP_OMP;
      FORALLSITES_OMP(i,s, default(shared) ){ su3mat_copy( &(s->staple), FLUX_LINKS(i,STAP_POS1) ); } END_LOOP_OMP;

      /* Shift staples to the sites where the t_links join them */
      FORALLUPDIR(j){disp[j] = trans_path[ mu2].disp[j];}
      mtag[M_STAP_POS2] = start_general_gather_field( (void *)FLUX_LINKS(0,STAP_POS2), sizeof(su3_matrix), disp, EVENANDODD, gen_pt[M_STAP_POS2] );

      /* Finish flux_links[i][STAP_NEG2] <- forward staple in dir2 */
      FORALLSITES_OMP(i,s, private(tmat1) ){
        mult_su3_an( TRANS_LINKS(i,trans_path2), (su3_matrix *)gen_pt[M_STAP_NEG2][i], &tmat1 );
        mult_su3_nn( &tmat1, TRANS_LINKS_F(i,TRANS_PATH2_F), FLUX_LINKS(i,STAP_NEG2) );
      } END_LOOP_OMP;

      wait_general_gather(mtag[M_STAP_POS2]);

      /* Collect results of the shift of forward staples */
      FORALLSITES_OMP(i,s, default(shared) ){ su3mat_copy( (su3_matrix *)gen_pt[M_STAP_POS2][i], &(s->staple)); } END_LOOP_OMP;
      FORALLSITES_OMP(i,s, default(shared) ){ su3mat_copy( &(s->staple), FLUX_LINKS(i,STAP_POS2) ); } END_LOOP_OMP;

      /* Prepare the forward space-like segment for parallel transport in time */
      /* Note: we time shift only one giant staple */
      /* s_link_f[i] <- s_link[i] */
      /* stap_pos1_f[i] <- stap_pos[i] */

      FORALLSITES_OMP(i,s, default(shared) ){
        su3mat_copy( FLUX_LINKS(i,S_LINK), FLUX_LINKS_F(i,S_LINK_F) );
        su3mat_copy( FLUX_LINKS(i,STAP_POS1), FLUX_LINKS_F(i,STAP_POS1_F) );
      } END_LOOP_OMP;
	
      cleanup_general_gather(mtag[M_STAP_NEG1]);
      cleanup_general_gather(mtag[M_STAP_POS1]);
      cleanup_general_gather(mtag[M_STAP_NEG2]);
      cleanup_general_gather(mtag[M_STAP_POS2]);

      /* Start gather of forward space-like segments for next t gen_pt[M_F_LINKS_F][i] will point to flux_links_f[i+TUP] */
      mtag[M_F_LINKS_F] = start_gather_field( (void *)flux_links_f, sizeof(su3_matrix)*NFLUX_F, xc[TUP], EVENANDODD, gen_pt[M_F_LINKS_F] );
	
	
      /* Recursively compute the Wilson loops of different time extent for fixed spatial extent */
      for(t=0;t<nct;t++){
 
        /* Collect forward space-like segments */
        wait_gather( mtag[M_F_LINKS_F]);
        FORALLSITES_OMP(i,s, default(shared) ){
          su3mat_copy( (su3_matrix *)(gen_pt[M_F_LINKS_F][i])+S_LINK_F, &(s->staple));
          su3mat_copy( (su3_matrix *)(gen_pt[M_F_LINKS_F][i])+STAP_POS1_F, &(s->diag));
        } END_LOOP_OMP;
        FORALLSITES_OMP(i,s, default(shared) ){
          su3mat_copy( &(s->staple), FLUX_LINKS_F(i,S_LINK_F) );
          su3mat_copy( &(s->diag), FLUX_LINKS_F(i,STAP_POS1_F) );
        } END_LOOP_OMP;
  
        /* Start gather for next t, if still needed. */
        if( t<(nct-1) ){
          restart_gather_field( (void *)flux_links_f, sizeof(su3_matrix)*NFLUX_F, xc[TUP], EVENANDODD, gen_pt[M_F_LINKS_F], mtag[M_F_LINKS_F] );
        }else{
          cleanup_gather( mtag[M_F_LINKS_F]);
        }

        /* Finally, compute the Wilson loops. */
        make_loops ( t, (wils_loop1+r+base_r), 
                     flux_links, flux_links_f,
                    trans_links, trans_links_f);
  
      } /* end loop over t */

      /* Start gather of forward time-like links for next r. */
      if( r<(maxc[mu]-1) ){
        restart_gather_field( (void *)trans_links_f, NTRANS_F*sizeof(su3_matrix), xc[mu], EVENANDODD, gen_pt[M_T_LINKS_F], mtag[M_T_LINKS_F] );
      }else{
        cleanup_gather( mtag[M_T_LINKS_F]);
      }
    } /* end loop over r */
  } /* end loop over  mu */
  assert( (base_r==nrmax) );

  free( trans_links ); 
  free( trans_links_f ); 
  free( flux_links_f );
  free( flux_links );

  return base_r;
} /* hybrid_loops */

