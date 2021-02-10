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

#define WILS_LOOP1(i,t,r)   wils_loop1[ (i)*mi*geom->maxlen + (t)*mi + r ]
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

static int mi;
static hqp_geom *geom=NULL;
static double *wils_loop1=NULL;

/* Constructeor and destructor for paths */
static link_path *alloc_trans_path ( int ntrans );
static void free_link_path ( link_path *path );

static inline void global_sums( int mysize );
static inline int hqp_buffer_dim (void );
static void hqp_free_buffers( void );
static void hqp_setup_buffers( int mysize );
static inline void output_all( char smtag[], int mi );

/* Build the big staples */
static void initialize_trans_links( int ntrans, link_path *trans_path, su3_matrix *trans_links);
/* Needs a special loop over r and t to build all types of Wilson loops at once */
static int hybrid_loops( link_path trans_path[], double *wils_loop1 );
/* Needs a special routine to contract all types of Wilson loops at once */
static void make_hybrid_loops ( int t,
                  double *wils_loop,
                  su3_matrix *flux_links,
                  su3_matrix *flux_links_f,
                  su3_matrix *trans_links, 
                  su3_matrix *trans_links_f );



void hybrid_loop1(int tot_smear) {

  register int i,j, mu,r,t;
  int mu1=0,mu2=0,trans_path1=0,trans_path2=0;
  register site *s;
  int  base_r=0;
  int disp3[4];
  char myname[] = "new_hybrid_loop1";

  char smtag[MAXSMTAG]="";
#ifdef SMEARING
  sprintf(smtag,"_%d",tot_smear);
#endif

  msg_tag *mtag[NMSGS];

  geom = hqp_geometry( myname );
  mi = hqp_buffer_dim();
  hqp_setup_buffers( geom->maxlen*mi*NWLOOP1 );

  /* Paths for transverse links */
  link_path *trans_path = alloc_trans_path ( NTRANS );

  /****************************************************************
   * 
   * Proceed to loop over directions and recursively construct the 
   * space-like segments and compute the Wilson loops with these 
   *
   ***************************************************************/

  base_r += hybrid_loops( trans_path, wils_loop1+base_r );
  assert( (base_r==mi) );
  
  /****************************************************************
   * 
   * Proceed to normalization and printing of the Wilson loops
   * Remark: do not average the directions (impossible if anisotropic),
   * instead print the direction into one column
   * can activate averaging of directions via compiler macro AVERAGE_DIRECTIONS
   *
   ***************************************************************/

  global_sums( geom->maxlen*mi*NWLOOP1 );

  output_all( smtag, mi );
  
  free_link_path ( trans_path );
  hqp_free_buffers( );
  hqp_free_geometry( geom );
  
} /* hybrid_loop1 */

link_path *alloc_trans_path ( int ntrans ) {

  int **dispc;
  /* Paths for transverse links */
  dispc=(int **)malloc(3*sizeof(int*));
  assert(dispc!=NULL);
  for ( int mu=XUP; mu<=ZUP; mu++ ) {
    dispc[mu] = (int *)malloc(4*sizeof(int));
    assert(dispc[mu]!=NULL);
    memset (dispc[mu],0,4*sizeof(int));
    dispc[mu][xc[mu]] = 2;
  }

  link_path *new_path=malloc(ntrans*sizeof(link_path));

  sprintf(new_path[XUP].name ,"%s", "XX" );
  sprintf(new_path[YUP].name ,"%s", "YY" );
  sprintf(new_path[ZUP].name ,"%s", "ZZ" );

  for ( int mu=XUP; mu<=ZUP; mu++ ) {
    new_path[mu].length = 2;
    for ( int mu1=XUP; mu1<=TUP; mu1++ ) 
      new_path[mu].disp[mu1] = dispc[mu][mu1];

    for ( int mu1=0; mu1<MAX_PATH_LENGTH; mu1++)
      new_path[mu].dir[mu1] = ( mu1 <2 ? xc[mu] : NODIR );
 
  }
  for ( int mu=XUP; mu<=ZUP; mu++ )  
    free(dispc[mu]); 
  free(dispc);

  return (new_path);
}

static void free_link_path ( link_path *path ) {

  free(path);
  path=NULL;
}

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

static int hybrid_loops( link_path trans_path[], double *wils_loop1 ) {
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
  trans_links   = hqp_alloc_su3mat_buffer (NTRANS);
  /* Allocate space for shifted auxiliary link products */
  trans_links_f   = hqp_alloc_su3mat_buffer (NTRANS_F);

  /* Compute and store products of transverse links */
  initialize_trans_links( NTRANS, trans_path, trans_links );

  /* Allocate space for timeward shifted flux tube shapes */
  flux_links_f   = hqp_alloc_su3mat_buffer (NFLUX_F);
  /* Allocate space for flux-tube shapes */
  flux_links   = hqp_alloc_su3mat_buffer (NFLUX);

  for( int mu=XUP; mu<=ZUP; base_r+=maxc[mu],mu++){
    
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
      FORALLUPDIR(j){ disp[xc[j]] = -trans_path[ mu1].disp[xc[j]]; }
      mtag[M_STAP_NEG1] = start_general_gather_field( (void *)FLUX_LINKS(0,S_LINK), sizeof(su3_matrix), disp, EVENANDODD, gen_pt[M_STAP_NEG1] );

      /* flux_links[i][STAP_POS1] <- forward staple at disp */
      FORALLSITES_OMP(i,s, private(tmat1) ){
        mult_su3_nn( TRANS_LINKS(i,trans_path1), FLUX_LINKS(i,S_LINK), &tmat1 );
        mult_su3_na( &tmat1, TRANS_LINKS_F(i,TRANS_PATH1_F), FLUX_LINKS(i,STAP_POS1) ); 
      } END_LOOP_OMP;

      /* Can't overlap gen_gathers!! */
      wait_general_gather(mtag[M_STAP_NEG1]);

      /* Start making flux_links[i][STAP_NEG2] <- backward s_link at -disp */
      FORALLUPDIR(j){disp[xc[j]] = -trans_path[ mu2].disp[xc[j]];}
      mtag[M_STAP_NEG2] = start_general_gather_field( (void *)FLUX_LINKS(0,S_LINK), sizeof(su3_matrix), disp, EVENANDODD, gen_pt[M_STAP_NEG2] );

      /* flux_links[i][STAP_POS2] <- forward staple at disp */
      FORALLSITES_OMP(i,s, private(tmat1) ){
        mult_su3_nn( TRANS_LINKS(i,trans_path2), FLUX_LINKS(i,S_LINK), &tmat1 );
        mult_su3_na( &tmat1, TRANS_LINKS_F(i,TRANS_PATH2_F), FLUX_LINKS(i,STAP_POS2) );
      } END_LOOP_OMP;
	
      wait_general_gather(mtag[M_STAP_NEG2]);

      /* Shift staples to the sites where the t_links join them */
      FORALLUPDIR(j){disp[xc[j]] = trans_path[ mu1].disp[xc[j]];}
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
      FORALLUPDIR(j){disp[xc[j]] = trans_path[ mu2].disp[xc[j]];}
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
      for(t=0;t<geom->maxlen;t++){
 
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
        if( t<(geom->maxlen-1) ){
          restart_gather_field( (void *)flux_links_f, sizeof(su3_matrix)*NFLUX_F, xc[TUP], EVENANDODD, gen_pt[M_F_LINKS_F], mtag[M_F_LINKS_F] );
        }else{
          cleanup_gather( mtag[M_F_LINKS_F]);
        }

        /* Finally, compute the Wilson loops. */
        make_hybrid_loops ( t, (wils_loop1+r+base_r), 
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
  assert( (base_r==mi) );

  hqp_free_su3mat_buffer( trans_links ); 
  hqp_free_su3mat_buffer( trans_links_f ); 
  hqp_free_su3mat_buffer( flux_links_f );
  hqp_free_su3mat_buffer( flux_links );

  return base_r;
} /* hybrid_loops */

static void initialize_trans_links( int ntrans, link_path *trans_path, su3_matrix *trans_links ) {

  su3_matrix *tempmat1 = hqp_alloc_su3mat_buffer( 1 );
  /* trans_links[i][j] is set to the link product for the
     jth path type that ends at site i */
  int i; 
  site *s; 
  for(int j = 0; j < ntrans; j++){
    path_product(trans_path[j].dir, trans_path[j].length, tempmat1);
    FORALLSITES_OMP(i,s, default(shared) ){ su3mat_copy( &(tempmat1[i]), TRANS_LINKS(i,j)); } END_LOOP_OMP;
  }
  hqp_free_su3mat_buffer ( tempmat1 );

}


static void make_hybrid_loops ( int t,
                  double *wils_loop1,
                  su3_matrix *flux_links,
                  su3_matrix *flux_links_f,
                  su3_matrix *trans_links, 
                  su3_matrix *trans_links_f ) {
  register int i;
  register site *s;
  su3_matrix tmat1,tmat2;
  su3_matrix *tmatp;
  double wl=0., sig=0., pi=0., delta=0.;

  /* Compute naive Wilson loop term */
  FORALLSITES_OMP(i,s, private(tmat1,tmat2,tmatp) reduction(+:wl) ){
    /* If the loop extends past t = nc[TUP] - 1 the temporal axial gauge link is nontrivial */
    if( (site_coord(s,xc[TUP])+t+1)>=nc[TUP] ){
      mult_su3_nn( &(s->link[xc[TUP]]), FLUX_LINKS_F(i,S_LINK_F), &tmat1);
      mult_su3_na( &tmat1, TRANS_LINKS_F(i,T_LINK_F), &tmat2);
      tmatp = &tmat2;
    }else
      tmatp = FLUX_LINKS_F(i,S_LINK_F);
    
    wl += (double)realtrace_su3( tmatp, FLUX_LINKS(i,S_LINK) );
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

    sig +=(double)(
      realtrace_su3( tmatp, FLUX_LINKS(i,STAP_POS1) ) +
      realtrace_su3( tmatp, FLUX_LINKS(i,STAP_NEG1) ) +
      realtrace_su3( tmatp, FLUX_LINKS(i,STAP_POS2) ) +
      realtrace_su3( tmatp, FLUX_LINKS(i,STAP_NEG2) ));

    pi +=(double)(
      realtrace_su3( tmatp, FLUX_LINKS(i,STAP_POS1) ) -
      realtrace_su3( tmatp, FLUX_LINKS(i,STAP_NEG1) ));


    delta +=(double)(
      realtrace_su3( tmatp, FLUX_LINKS(i,STAP_POS1) ) +
      realtrace_su3( tmatp, FLUX_LINKS(i,STAP_NEG1) ) -
      realtrace_su3( tmatp, FLUX_LINKS(i,STAP_POS2) ) -
      realtrace_su3( tmatp, FLUX_LINKS(i,STAP_NEG2) ));
  } END_LOOP_OMP;
  WILS_LOOP1(STAP_SIG_GP1,t,0) = sig;
  WILS_LOOP1(STAP_PI_U1,t,0) = pi;
  WILS_LOOP1(STAP_DELTA_G1,t,0) = delta;

} /* make_hybrid_loops */


void output_all( char smtag[], int mi ) {

  int disp[4]={0,0,0,0};
  int *r, base_r=0;

  for( disp[TUP]=1;disp[TUP]<=geom->maxlen;disp[TUP]++) {
    for( int mu=XUP, base_r=0; mu<=ZUP; base_r+=maxc[mu],mu++) {
      memset(disp, 0,3*sizeof(int));
      r=&(disp[mu]);
      for((*r)=1;(*r)<=maxc[mu];(*r)++) {
        if ( hqp_disp_rsq_ok( disp, geom ) == 1 ) 
#ifndef AVERAGE_DIRECTIONS
          hqp_output_corr("WILS_LOOP1",smtag, disp, 
          WILS_LOOP1(W_LOOP1,(disp[TUP]-1),((*r)-1)+base_r)/3. );
#else 
          WILS_LOOP1(W_LOOP1,(disp[TUP]-1),((*r)-1)) +=(mu>XUP?1:-2)* 
          WILS_LOOP1(W_LOOP1,(disp[TUP]-1),((*r)-1)+base_r)/3.;
#endif
      }   
    }   
#ifdef AVERAGE_DIRECTIONS
    memset(disp, 0,3*sizeof(int));
    disp[TUP]--;
    r=&(disp[XUP]);
    for((*r)=0;*r<maxc[XUP];(*r)++)
      if ( hqp_disp_rsq_ok( disp, geom ) == 1 ) 
        hqp_output_corr("WILS_LOOP1",smtag, disp, 
        WILS_LOOP1(W_LOOP1,(disp[TUP]),(*r))/3. );
    disp[TUP]++;
#endif
  }
  
  for( disp[TUP]=1;disp[TUP]<=geom->maxlen;disp[TUP]++) {
    for( int mu=XUP, base_r=0; mu<=ZUP; base_r+=maxc[(mu%3)],mu++) {
      memset(disp, 0,3*sizeof(int));
      r=&(disp[mu]);
      for((*r)=1;(*r)<=maxc[mu];(*r)++) {
        if ( hqp_disp_rsq_ok( disp, geom ) == 1 ) 
#ifndef AVERAGE_DIRECTIONS
          hqp_output_corr("STAP_SIG_GP1",smtag, disp, 
          WILS_LOOP1(STAP_SIG_GP1,(disp[TUP]-1),((*r)-1)+base_r)/3. );
#else 
          WILS_LOOP1(STAP_SIG_GP1,(disp[TUP]-1),((*r)-1)) +=(mu>XUP?1:-2)* 
          WILS_LOOP1(STAP_SIG_GP1,(disp[TUP]-1),((*r)-1)+base_r)/3.;
#endif
      }   
    }   
#ifdef AVERAGE_DIRECTIONS
    memset(disp, 0,3*sizeof(int));
    disp[TUP]--;
    r=&(disp[xc[(XUP%3)]]);
    for((*r)=0;*r<maxc[XUP%3];(*r)++)
      if ( hqp_disp_rsq_ok( disp, geom ) == 1 ) 
        hqp_output_corr("STAP_SIG_GP1",smtag, disp, 
        WILS_LOOP1(STAP_SIG_GP1,(disp[TUP]),(*r)+base_r)/3. );
    disp[TUP]++;
#endif
  }
  
  for( disp[TUP]=1;disp[TUP]<=geom->maxlen;disp[TUP]++) {
    for( int mu=XUP, base_r=0; mu<=ZUP; base_r+=maxc[(mu%3)],mu++) {
      memset(disp, 0,3*sizeof(int));
      r=&(disp[mu]);
      for((*r)=1;(*r)<=maxc[mu];(*r)++) {
        if ( hqp_disp_rsq_ok( disp, geom ) == 1 ) 
#ifndef AVERAGE_DIRECTIONS
          hqp_output_corr("STAP_PI_U1",smtag, disp, 
          WILS_LOOP1(STAP_PI_U1,(disp[TUP]-1),((*r)-1)+base_r)/3. );
#else 
          WILS_LOOP1(STAP_PI_U1,(disp[TUP]-1),((*r)-1)) +=(mu>XUP?1:-2)* 
          WILS_LOOP1(STAP_PI_U1,(disp[TUP]-1),((*r)-1)+base_r)/3.;
#endif
      }   
    }   
#ifdef AVERAGE_DIRECTIONS
    memset(disp, 0,3*sizeof(int));
    disp[TUP]--;
    r=&(disp[XUP]);
    for((*r)=0;*r<maxc[XUP];(*r)++)
      if ( hqp_disp_rsq_ok( disp, geom ) == 1 ) 
        hqp_output_corr("STAP_PI_U1",smtag, disp, 
        WILS_LOOP1(STAP_PI_U1,(disp[TUP]),(*r)+base_r)/3. );
    disp[TUP]++;
#endif
  }
  
  for( disp[TUP]=1;disp[TUP]<=geom->maxlen;disp[TUP]++) {
    for( int mu=XUP, base_r=0; mu<=ZUP; base_r+=maxc[(mu%3)],mu++) {
      memset(disp, 0,3*sizeof(int));
      r=&(disp[mu]);
      for((*r)=0;(*r)<=maxc[mu];(*r)++) {
        if ( hqp_disp_rsq_ok( disp, geom ) == 1 ) 
#ifndef AVERAGE_DIRECTIONS
          hqp_output_corr("STAP_DELTA_G1",smtag, disp, 
          WILS_LOOP1(STAP_DELTA_G1,(disp[TUP]-1),((*r)-1)+base_r)/3. );
#else 
          WILS_LOOP1(STAP_DELTA_G1,(disp[TUP]-1),((*r)-1)) +=(mu>XUP?1:-2)* 
          WILS_LOOP1(STAP_DELTA_G1,(disp[TUP]-1),((*r)-1)+base_r)/3.;
#endif
      }   
    }   
#ifdef AVERAGE_DIRECTIONS
    memset(disp, 0,3*sizeof(int));
    disp[TUP]--;
    r=&(disp[XUP]);
    for((*r)=0;*r<maxc[XUP];(*r)++)
      if ( hqp_disp_rsq_ok( disp, geom ) == 1 ) 
        hqp_output_corr("STAP_DELTA_G1",smtag, disp, 
        WILS_LOOP1(STAP_DELTA_G1,(disp[TUP]),(*r)+base_r)/3. );
    disp[TUP]++;
#endif
  }

}
