/****** fermion_force_hisq_multi.c  -- ******************/
/* MIMD version 7 */
/* Multisource fermion force.  Includes various optimization choices.
 * 
 * 1. General force for any FN-type action. Based on fermion_force_eo_milc.c
 *    and optimized to transport only one set of SU(3) matrices.
 *
 * 2. Same as 1 but with indices on CG solutions reversed
 *    to maybe improve cache hits
 *
 * 3. Generalization of fermion_force_asqtad3.c optimized
 *    for Asqtad.  Transports the full array of source vectors.
 *
 * D.T. 12/05 Version  3. created for Asqtad RHMC.
 * D.T.  6/06 Versions 1. and 2. created.
 * CD   10/06 Collected versions in this file.
 * D.T. 6/07 Two levels of smearing
 * A.B. 10/07 unitarization and its derivatives,
 *            "wrapper" that allows arbitraty levels of smearing
 * A.B. 2/08 combine forces from level 2 smearing and then
 *           calculate force from reunitarization and level 1 smearing
 * A.B. 8/08 Rearrange force to reduce to two Asqtad steps plus
 *           whatever correction to Naik term is present (for c-quark)
 *
 */


#include "generic_ks_includes.h"
#define SU3_MAT_OP_NO_STORAGE
#include "../include/su3_mat_op.h"
#include "../include/umethod.h"
#include "../include/fermion_links.h"
#include <string.h>



//AB output for debugging
#ifdef MILC_GLOBAL_DEBUG
#ifdef HISQ_FF_DEBUG
static void dumpt4( su3_tensor4 dwdv ) {
  int i,j,k,l;
  i=j=k=l=0;

  for(i=0;i<3;i++) {
    for(j=0;j<3;j++) {
      for(k=0;k<3;k++) {
        for(l=0;l<3;l++) {
          printf("T4[%d][%d][%d][%d].real=%lf  T4[%d][%d][%d][%d].imag=%lf\n",
            i,j,k,l,dwdv.t4[i][j][k][l].real,
            i,j,k,l,dwdv.t4[i][j][k][l].imag );
        }
      }
    }
  }
}
#endif /* HISQ_FF_DEBUG */
#endif /* MILC_GLOBAL_DEBUG */

/* All routines in this file require the FN flag */
#ifndef FN
BOMB THE COMPILE
#endif

#define GOES_FORWARDS(dir) (dir<=TUP)
#define GOES_BACKWARDS(dir) (dir>TUP)

#ifdef QCDOC
#define special_alloc qcdoc_alloc
#define special_free qfree
#else
#define special_alloc malloc
#define special_free free
#endif


#ifdef MILC_GLOBAL_DEBUG
// dump matrix on this site/link
#define AB_OUT_ON_SITE 0
#define AB_OUT_ON_LINK 0
#endif /* MILC_GLOBAL_DEBUG */


// Forward declarations

static void 
fn_fermion_force_multi_hisq_mx( info_t *info, Real eps, Real *residues, 
				su3_vector **multi_x, int nterms, 
				fermion_links_t *fl );

// In the following routines smearing and reunitarization
// are split into separate routines. This allows for arbitrary levels of
// smearing and/or reunitarization.

// Combine smearing level 2 from different multi_x,
//               use globals from lattice.h
// HISQ wrapper routine: contribution to the force due to
//    0) smearing level 0 (outer product |X><Y|)
//    1) smearing level 2 + smearing level 2 Naik corrected
//    2) reunitarization
//    3) smearing level 1


/* GPU code always uses wrapper */
#if defined(HISQ_FF_MULTI_WRAPPER) || defined(USE_FF_GPU)

#ifndef USE_FF_GPU
static void 
fn_fermion_force_multi_hisq_wrapper_mx_cpu( info_t *info, Real eps, Real *residues,
					su3_vector **multi_x, int nterms,
					fermion_links_t *fl);
// Contribution to the force from 0 level of smearing, force contains
// only |X><Y| terms
static void 
fn_fermion_force_multi_hisq_smearing0( info_t *info, Real eps, Real *residues, 
				       su3_vector **multi_x, 
				       int nterms,
				       su3_matrix *force_accum[4], 
				       su3_matrix *force_accum_naik[4]);

// This routine calculates contribution to the force due to i-th level
// of smearing, as input it uses array of new links and the force from
// (i-1)-th level of smearing
static void 
fn_fermion_force_multi_hisq_smearing( info_t *info, Real eps, 
				      Real *residues, 
				      su3_vector **multi_x, 
				      int nterms, 
				      su3_matrix *force_accum[4], 
				      su3_matrix *force_accum_old[4], 
				      su3_matrix *force_accum_naik_old[4], 
				      su3_matrix *internal_U_link,
				      int internal_num_q_paths, 
				      Q_path *internal_q_paths_sorted, 
				      int *internal_netbackdir_table);
// Contribution to the force due to reunitarization of links:
//   Y=V*(V^+*V)^-1/2, W=Y/(detY)^1/3
static void 
fn_fermion_force_multi_hisq_reunit( info_t *info, su3_matrix *force_accum[4], 
				    su3_matrix *force_accum_old[4], 
				    fermion_links_t *fl);

#else

static void 
fn_fermion_force_multi_hisq_wrapper_mx_gpu( info_t *info, Real eps, Real *residues,
					su3_vector **multi_x, int nterms,
					fermion_links_t *fl);
#endif /* ifndef USE_FF_GPU */


/* Change wrapper call for GPU */

#ifdef USE_FF_GPU
#define fn_fermion_force_multi_hisq_wrapper_mx fn_fermion_force_multi_hisq_wrapper_mx_gpu
#else 
#define fn_fermion_force_multi_hisq_wrapper_mx fn_fermion_force_multi_hisq_wrapper_mx_cpu
#endif

#endif /* HISQ_FF_MULTI_WRAPPER || USE_FF_GPU */

/**********************************************************************/
/*   Wrapper for fermion force routines with multiple sources         */
/**********************************************************************/
#ifdef FFTIME
static const char *qop_prec[2] = {"F", "D"};
#endif

void eo_fermion_force_multi( Real eps, Real *residues, su3_vector **xxx, 
			     int nterms, int prec, fermion_links_t *fl) {
  double dtime = -dclock();
  info_t info = INFO_ZERO;

  if(prec != PRECISION){
    node0_printf("eo_fermion_force_multi: WARNING, precision requests not supported. Using %d.\n",
		 PRECISION);
  }
  switch(KS_MULTIFF){
  case FNMATREV:
    node0_printf("FNMATREV isn't implemented\n"); exit(0);
    break;
  case FNMAT:
#if defined(HISQ_FF_MULTI_WRAPPER) || defined(USE_FF_GPU)
    fn_fermion_force_multi_hisq_wrapper_mx(&info,  eps, residues, xxx, nterms, fl );
#else
    fn_fermion_force_multi_hisq_mx(&info,  eps, residues, xxx, nterms, fl );
#endif
    break;
  default:
    fn_fermion_force_multi_hisq_mx(&info, eps, residues, xxx, nterms, fl );
  }

  dtime += dclock();
  info.final_sec = dtime;
#ifdef FFTIME
  node0_printf("FFTIME:  time = %e (HISQ %s) terms = %d flops/site = %d mflops = %e\n",
	       info.final_sec,qop_prec[PRECISION-1],nterms,
	       (int)(info.final_flop*numnodes()/volume),
	       info.final_flop/(1e6*info.final_sec) );
#endif

  /* Sum over all nodes and update counter */
  g_intsum(&INFO_HISQ_SVD_COUNTER(&info));
  hisq_svd_counter += INFO_HISQ_SVD_COUNTER(&info);
  
  g_intsum(&INFO_HISQ_FORCE_FILTER_COUNTER(&info));
  hisq_force_filter_counter += INFO_HISQ_FORCE_FILTER_COUNTER(&info);

}

/**********************************************************************/
/*   Accessor for string describing the option                        */
/**********************************************************************/
const char *ks_multiff_opt_chr( void )
{
  switch(KS_MULTIFF){
  case ASVEC:
    return "ASVEC";
    break;
  case FNMATREV:
    return "FNMATREV";
    break;
  case FNMAT:
    return "FNMAT";
    break;
  default:
    return "FNMAT";
  }
  return NULL;
}
/**********************************************************************/
/*   General FN Version for "nterms" sources                          */
/**********************************************************************/

/* update the  momenta with the fermion force */
/* Assumes that the multimass conjugate gradient has been run, with the answer in
   multi_x, and dslash_site(multi_x,multi_x,ODD) has been run. (fills in multi_x on odd sites) */
/* SEE LONG COMMENTS AT END */

// For smearing 1
// Quark paths sorted by net displacement and last directions
static Q_path *q_paths_sorted_1 = NULL;	
// table of net path displacements (backwards from usual convention)
static int *netbackdir_table_1 = NULL; 
// table of gather directions to bring start of path to end, for "FN"
// = "fat-Naik" actions
// static int net_back_dirs_1[8] = 
//   { XDOWN, YDOWN, ZDOWN, TDOWN, XUP, YUP, ZUP, TUP };

// For smearing 2
// Quark paths sorted by net displacement and last directions
static Q_path *q_paths_sorted_2 = NULL;
// table of net path displacements (backwards from usual convention)	
static int *netbackdir_table_2 = NULL; 
// table of gather directions to bring start of path to end, for "FN"
// = "fat-Naik" actions
// static int net_back_dirs_2[16] = 
//   { XDOWN, YDOWN, ZDOWN, TDOWN, XUP, YUP, ZUP, TUP, 
//     X3DOWN, Y3DOWN, Z3DOWN, T3DOWN, X3UP, Y3UP, Z3UP, T3UP };
// Quark paths sorted by net displacement and last directions
static Q_path *q_paths_sorted_3 = NULL;	
// table of net path displacements (backwards from usual convention)
static int *netbackdir_table_3 = NULL; 
// table of gather directions to bring start of path to end, for "FN"
// = "fat-Naik" actions
// static int net_back_dirs_3[16] = 
//   { XDOWN, YDOWN, ZDOWN, TDOWN, XUP, YUP, ZUP, TUP, 
//     X3DOWN, Y3DOWN, Z3DOWN, T3DOWN, X3UP, Y3UP, Z3UP, T3UP };
static Q_path *q_paths_sorted_current = NULL;
static int *netbackdir_table_current = NULL;

// table of gather directions to bring start of path to end, for "FN"
// = "fat-Naik" actions
static int net_back_dirs[16] = 
  { XDOWN, YDOWN, ZDOWN, TDOWN, XUP, YUP, ZUP, TUP, 
    X3DOWN, Y3DOWN, Z3DOWN, T3DOWN, X3UP, Y3UP, Z3UP, T3UP };

static int sort_quark_paths( Q_path *src_table, Q_path *dest_table, 
			     int npaths, int num_back_dirs );
static int find_backwards_gather( Q_path *path );
static int first_force=1;	// 1 if force hasn't been called yet






// Fermion force SKETCH:
// 1) loop on different epsilon corrections to Naik term
//   a) calculate force from smearing level 2 on full multi_x array
//      with epsilon correction = 0
//   b) calculate all other smearing level 2 contributions on 
//      corresponding parts of multi_x array
//   c) add all contributions from smearing level 2
// 2) perform reunitarization
// 3) calculate force from smearing level 1
// This procedure is slightly more than 2 Asqtad smearing steps
static void 
fn_fermion_force_multi_hisq_mx( info_t *info, Real eps, Real *residues, 
				su3_vector **multi_x, int nterms, fermion_links_t *fl)
{
  // Get the terms we need from the fermion links structure */
  ks_action_paths_hisq *ap = get_action_paths_hisq(fl);
  hisq_auxiliary_t *aux = get_hisq_auxiliary(fl);
  int n_naiks = fermion_links_get_n_naiks(fl);
  double *eps_naik = fermion_links_get_eps_naik(fl);
  su3_matrix *U_link = aux->U_link;
  su3_matrix *V_link = aux->V_link;
  su3_matrix *W_unitlink = aux->W_unitlink;
  su3_matrix *Y_unitlink = aux->Y_unitlink;
  // note CG_solution and Dslash * solution are combined in "multi_x" 
  // Use forward part of Dslash to get force 
  // see long comment at end 
  // This version is for two stage smearing, includes unitarization step

  // ALGORITHM SKETCH:
  // For each net displacement
    // compute outer product |X><X|
    // For each path in smearing 2 with that displacement
      // parallel transport outer product to each link in path, accumulate for
      // each link
    // For each path in smearing 1
      // parallel transport accumulated outer product from previous step to
      // each link in path, accumulate new one at each link
    // add to momentum

  // We need sorted path tables for both smearings

  /* note CG_solution and Dslash * solution are combined in "multi_x" */
  /* New version 1/21/99.  Use forward part of Dslash to get force */
  /* see long comment at end */
  /* For each link we need multi_x transported from both ends of path. */
  int term;
  register int i,j,k,lastdir=-99,ipath,ilink;
  register site *s;
  int length,dir,odir;
  su3_matrix tmat,tmat2;
  Real ferm_epsilon_1, ferm_epsilon_2, coeff;
  int num_q_paths_1 = ap->p1.num_q_paths;
  Q_path *q_paths_1 = ap->p1.q_paths;
  int num_q_paths_2 = ap->p2.num_q_paths;
  Q_path *q_paths_2 = ap->p2.q_paths;
  int num_q_paths_3 = ap->p3.num_q_paths;
  Q_path *q_paths_3 = ap->p3.q_paths;
  Q_path *this_path;	// pointer to current path
  Q_path *last_path;	// pointer to previous path
  msg_tag *mtag[2];
  su3_matrix *mat_tmp0,*mat_tmp1;
  su3_matrix *oprod_along_path[MAX_PATH_LENGTH+1]; // length N path has N+1 sites!!
  su3_matrix *mats_along_path[MAX_PATH_LENGTH+1]; // 
  su3_matrix *force_accum_1[4];  // accumulate force through smearing 1
  su3_matrix *force_accum_2[4];  // accumulate force through smearing 2 (done first)
  su3_matrix *force_accum_tmp[4];  // to combine contributions at smearing 2
  int netbackdir, last_netbackdir;	// backwards direction for entire path
  su3_vector **internal_multi_x;
  Real *internal_residues;
  int num_q_paths_current,n_orders_naik_current;
  Real coeff_mult;
  su3_tensor4 dydv, dydagdv, dwdv, dwdagdv;
  su3_matrix force_tmp;
  int m, n, l;
  complex ftmp;
  int inaik, n_naik_shift;
  int aux_stride = 4;
  size_t nflops = 0;
  size_t nflops_all_sites = 0;

  /* Assume that for the fermion force calculation the KS phases are never
     defined for a nonzero origin */
  int r0[4] = {0, 0, 0, 0};
//int tempflops = 0; //TEMP


//  int nflop_old = 1847016 + 6264*(n_naiks-1) 
//    + 720*(2*n_order_naik_total-n_orders_naik[0])
//    + 10488; // HISQ action 10/30/07 version of code;

// CD 11/23/2011: With unit_der_analytic and without SVD the tally is closer to 
// 1877780 + 6264*(n_naiks-1) + 720*(2*n_order_naik_total-n_orders_naik[0]) + 10488
// 
  double dtime=-dclock();
  /* node0_printf("STARTING fn_fermion_force_multi_hisq() nterms = %d\n",nterms); */
  if( nterms==0 )return;

  for(i=0;i<=MAX_PATH_LENGTH;i++){
     oprod_along_path[i] = (su3_matrix *) malloc(sites_on_node*sizeof(su3_matrix) );
  }
  for(i=1;i<=MAX_PATH_LENGTH;i++){ // 0 element is never used (it's unit matrix)
     mats_along_path[i] = (su3_matrix *) malloc(sites_on_node*sizeof(su3_matrix) );
  }
  for(i=XUP;i<=TUP;i++){
     force_accum_1[i] = (su3_matrix *) malloc(sites_on_node*sizeof(su3_matrix) );
     force_accum_2[i] = (su3_matrix *) malloc(sites_on_node*sizeof(su3_matrix) );
     force_accum_tmp[i] = (su3_matrix *) malloc(sites_on_node*sizeof(su3_matrix) );
  }
  mat_tmp0 = (su3_matrix *) special_alloc(sites_on_node*sizeof(su3_matrix) );
  mat_tmp1 = (su3_matrix *) special_alloc(sites_on_node*sizeof(su3_matrix) );
  if( mat_tmp1 == NULL ){printf("Node %d NO ROOM\n",this_node); exit(0); }

	
  ferm_epsilon_1 = 2.0*eps; // we only do forward paths, gives factor of 2
  ferm_epsilon_2 =     1.0; // both forward and backward  CHECK THIS!!
  if( first_force==1 ){
    if( q_paths_sorted_1==NULL ) 
      q_paths_sorted_1 = (Q_path *)malloc( num_q_paths_1*sizeof(Q_path) );
    if(netbackdir_table_1==NULL ) 
      netbackdir_table_1 = (int *)malloc( num_q_paths_1*sizeof(int) );
    if( q_paths_sorted_2==NULL ) 
      q_paths_sorted_2 = (Q_path *)malloc( num_q_paths_2*sizeof(Q_path) );
    if(netbackdir_table_2==NULL ) 
      netbackdir_table_2 = (int *)malloc( num_q_paths_2*sizeof(int) );
    if( q_paths_sorted_3==NULL ) 
      q_paths_sorted_3 = (Q_path *)malloc( num_q_paths_3*sizeof(Q_path) );
    if(netbackdir_table_3==NULL ) 
      netbackdir_table_3 = (int *)malloc( num_q_paths_3*sizeof(int) );
    else{ node0_printf("WARNING: remaking sorted path tables\n"); exit(0); }
    // make sorted tables
    sort_quark_paths( q_paths_1, q_paths_sorted_1, num_q_paths_1, 8 );
    for( ipath=0; ipath<num_q_paths_1; ipath++ )
      netbackdir_table_1[ipath] = 
	find_backwards_gather( &(q_paths_sorted_1[ipath]) );
    sort_quark_paths( q_paths_2, q_paths_sorted_2, num_q_paths_2, 16 );
    for( ipath=0; ipath<num_q_paths_2; ipath++ )
      netbackdir_table_2[ipath] = 
	find_backwards_gather( &(q_paths_sorted_2[ipath]) );
    sort_quark_paths( q_paths_3, q_paths_sorted_3, num_q_paths_3, 16 );
    for( ipath=0; ipath<num_q_paths_3; ipath++ )
      netbackdir_table_3[ipath] = 
	find_backwards_gather( &(q_paths_sorted_3[ipath]) );
    first_force=0;
  }

  // clear force accumulators
  for(dir=XUP;dir<=TUP;dir++)FORALLFIELDSITES(i){
    clear_su3mat( &(force_accum_1[dir][i]) );
    clear_su3mat( &(force_accum_2[dir][i]) );
  }



  // loop on different naik epsilons
  n_naik_shift = 0;
  for( inaik=0; inaik<n_naiks; inaik++ ) {
    // clear temporary force accumulator
    for(dir=XUP;dir<=TUP;dir++)FORALLFIELDSITES(i){
      clear_su3mat( &(force_accum_tmp[dir][i]) );
    }
#ifdef MILC_GLOBAL_DEBUG
    node0_printf("fn_fermion_force_multi_hisq_mx: eps_naik[%d]=%f\n",inaik,eps_naik[inaik]);
#endif /* MILC_GLOBAL_DEBUG */


#ifdef MILC_GLOBAL_DEBUG
    {
    double tabletime=-dclock();
#endif /* MILC_GLOBAL_DEBUG */
#ifdef MILC_GLOBAL_DEBUG
    tabletime+=dclock();
    node0_printf("fftime: time to load fermion links %e\n",tabletime);
    }
#endif /* MILC_GLOBAL_DEBUG */


    internal_multi_x = multi_x + n_naik_shift;
    internal_residues = residues + n_naik_shift;


  // BEGIN SMEARING TWO PART
  /* loop over paths, and loop over links in path */
  last_netbackdir = NODIR;
  last_path = NULL;
  if( 0==inaik ) {
    q_paths_sorted_current = q_paths_sorted_2;
    num_q_paths_current = num_q_paths_2;
    netbackdir_table_current = netbackdir_table_2;
  }
  else {
    q_paths_sorted_current = q_paths_sorted_3;
    num_q_paths_current = num_q_paths_3;
    netbackdir_table_current = netbackdir_table_3;
  }
  for( ipath=0; ipath<num_q_paths_current; ipath++ ){
    this_path = &(q_paths_sorted_current[ipath]);
    /* skip backwards dslash - only in smearing two*/
    if(this_path->forwback== -1)continue;	

    length = this_path->length;
    // find gather to bring multi_x[term] from "this site" to end of path
    netbackdir = netbackdir_table_current[ipath];
    // and bring multi_x to end - no gauge transformation !!
    // resulting outer product matrix has gauge transformation
    // properties of a connection from start to end of path
    // Don't need to repeat this if same net disp. as last path
    if( netbackdir != last_netbackdir){ 
	k=0; // which gather we are using
        mtag[k] = start_gather_field( internal_multi_x[0], sizeof(su3_vector),
           netbackdir, EVENANDODD, gen_pt[k] );
        FORALLFIELDSITES(i){
	  clear_su3mat(  &oprod_along_path[0][i] ); // actually last site in path
        }
        if( 0==inaik ) {
        	// for first set of X links full multi_x is used
          n_orders_naik_current = n_order_naik_total;
        }
        else {
          // for other sets of X links only their part of multi_x is used
          n_orders_naik_current = n_orders_naik[inaik];
        }
        for(term=0;term<n_orders_naik_current;term++){
          if(term<n_orders_naik_current-1)
	    mtag[1-k] = start_gather_field( internal_multi_x[term+1],
              sizeof(su3_vector), netbackdir, EVENANDODD, gen_pt[1-k] );
          wait_gather(mtag[k]);
          FORALLFIELDSITES(i){
	    su3_projector( &internal_multi_x[term][i], 
			   (su3_vector *)gen_pt[k][i], &tmat );
	    scalar_mult_add_su3_matrix( &oprod_along_path[0][i], &tmat, 
	         internal_residues[term], &oprod_along_path[0][i] );
          }
	  nflops += 54 + 36;
          cleanup_gather(mtag[k]);
	  k=1-k; // swap 0 and 1
        } /* end loop over terms in rational function expansion */
//tempflops+=54*n_orders_naik_current;
//tempflops+=36*n_orders_naik_current;
    }

    /* path transport the outer product, or projection matrix, of multi_x[term]
       (EVEN sites)  and Dslash*multi_x[term] (ODD sites) from far end.

       maintain a matrix of the outer product transported backwards
	along the path to all sites on the path.
	If new "net displacement", need to completely recreate it.
	Otherwise, use as much of the previous path as possible 

	Note this array is indexed backwards - the outer product transported
	to site number n along the path is in oprod_along_path[length-n].
	This makes reusing it for later paths easier.

	Sometimes we need this at the start point of the path, and sometimes
	one link into the path, so don't always have to do the last link. */

    // figure out how much of the outer products along the path must be
    // recomputed. j is last one needing recomputation. k is first one.
    j=length-1; // default is recompute all
    if( netbackdir == last_netbackdir )
      while ( j>0 && this_path->dir[j] == 
	      last_path->dir[j+last_path->length-length] ) j--;
    if( GOES_BACKWARDS(this_path->dir[0]) ) k=1; else k=0;

    for(ilink=j;ilink>=k;ilink--){

      link_transport_connection_hisq( oprod_along_path[length-ilink-1], W_unitlink,
      oprod_along_path[length-ilink], mat_tmp0, this_path->dir[ilink]  );

      nflops += 198;
//tempflops+=9*22;
    }


   /* maintain an array of transports "to this point" along the path.
	Don't recompute beginning parts of path if same as last path */
    ilink=0; // first link where new transport is needed
    if( last_path != NULL )
      while( this_path->dir[ilink] == last_path->dir[ilink] ) ilink++ ;
    // Sometimes we don't need the matrix for the last link
    if( GOES_FORWARDS(this_path->dir[length-1]) ) k=length-1; else k=length;

    for( ; ilink<k; ilink++ ){
      if( ilink==0 ){
        dir = this_path->dir[0];
          if( GOES_FORWARDS(dir) ){
	    mtag[1] = declare_strided_gather( W_unitlink + dir,
		      aux_stride*sizeof(su3_matrix), sizeof(su3_matrix),
		      OPP_DIR(dir), EVENANDODD, gen_pt[1] );
	    prepare_gather(mtag[1]);
	    do_gather(mtag[1]);
	    //            mtag[1] = start_gather_field( W_unitlink[dir], sizeof(su3_matrix),
	    //                   OPP_DIR(dir), EVENANDODD, gen_pt[1] );
            wait_gather(mtag[1]);
            FORALLFIELDSITES(i){ 
	      su3_adjoint( (su3_matrix *)gen_pt[1][i], 
			   &(mats_along_path[1][i]) ); 
	    }
            cleanup_gather(mtag[1]);
          }
          else{
            FORALLFIELDSITES(i){ 
	      mats_along_path[1][i] = W_unitlink[4*i+OPP_DIR(dir)]; 
	    }
          }
      }
      else { // ilink != 0
        dir = OPP_DIR(this_path->dir[ilink]);
        link_transport_connection_hisq( mats_along_path[ilink], W_unitlink,
        mats_along_path[ilink+1], mat_tmp0, dir  );
	nflops += 198;
//tempflops+=9*22;
      }
    } // end loop over links

    /* A path has (length+1) points, counting the ends.  At first
	 point, no "down" direction links have their momenta "at this
	 point". At last, no "up" ... */
    if( GOES_FORWARDS(this_path->dir[length-1]) ) k=length-1; else k=length;
    for( ilink=0; ilink<=k; ilink++ ){
      if(ilink<length)dir = this_path->dir[ilink];
      else dir=NODIR;
      coeff = ferm_epsilon_2*this_path->coeff;
      if( (ilink%2)==1 )coeff = -coeff;

      /* add in contribution to the force */
      if( ilink<length && GOES_FORWARDS(dir) ){
	link_gather_connection_hisq( oprod_along_path[length-ilink-1], 
				     mat_tmp1, NULL, dir );
	if ( ilink==0 ){
	   FORALLFIELDSITES(i){ mat_tmp0[i] = mat_tmp1[i]; }
	}
	else {
          FORALLFIELDSITES(i){
	     mult_su3_na( &(mat_tmp1[i]),  
			  &(mats_along_path[ilink][i]), &(mat_tmp0[i]) );
          }
	  nflops += 198;
//tempflops+=9*22;
	}

        FORALLFIELDSITES(i){
	  // sign fixed AB 9/20/07
	  scalar_mult_add_su3_matrix( &(force_accum_tmp[dir][i]), 
				      &(mat_tmp0[i]), coeff, 
				      &(force_accum_tmp[dir][i]) );
	}
	nflops += 36;
//tempflops+=36;
      }
      if( ilink>0 && GOES_BACKWARDS(lastdir) ){
	odir = OPP_DIR(lastdir);
        if ( ilink==1 ){
	    FORALLFIELDSITES(i){ 
	      su3_adjoint( &(oprod_along_path[length-ilink][i]),
		&( mat_tmp0[i] ) ); 
	    }
	}
	else {
	    link_gather_connection_hisq( mats_along_path[ilink-1], mat_tmp1, 
					 NULL, odir );
            FORALLFIELDSITES(i){
              mult_su3_na( &(mat_tmp1[i]), 
			   &(oprod_along_path[length-ilink][i]),  
			   &(mat_tmp0[i]) );
            }
	    nflops += 198;
//tempflops+=9*22;
	}
        FORALLFIELDSITES(i){
	  // sign fixed AB 9/20/07
	  scalar_mult_add_su3_matrix( &(force_accum_tmp[odir][i]), 
				      &(mat_tmp0[i]), +coeff, 
				      &(force_accum_tmp[odir][i]) );
        }
	nflops += 36;
//tempflops+=36;
      }

      lastdir = dir;
    } /* end loop over links in path */
    last_netbackdir = netbackdir;
    last_path = &(q_paths_sorted_current[ipath]);
  } /* end loop over paths */

    if( 0==inaik ) {
      coeff_mult = 1.0;
    }
    else {
      coeff_mult = eps_naik[inaik];
    }
    for(dir=XUP;dir<=TUP;dir++) {
      FORALLFIELDSITES(i) {
//        add_su3_matrix( &( force_accum_2[dir][i] ), &( force_accum_tmp[dir][i] ), &( force_accum_2[dir][i] ) );
        scalar_mult_add_su3_matrix( &( force_accum_2[dir][i] ), 
            &( force_accum_tmp[dir][i] ), coeff_mult,
            &( force_accum_2[dir][i] ) );
      }
      nflops += 36;
//tempflops+=36;
    }
    n_naik_shift += n_orders_naik[inaik];
  } /* end loop over naik epsilons */

  // END SMEARING TWO PART

#ifdef FFTIME
#ifdef MILC_GLOBAL_DEBUG
  {
  double tmptime=dtime+dclock();
  node0_printf("fftime(hisq-smearing1):  time = %e \n",tmptime);
  }
#endif /* MILC_GLOBAL_DEBUG */
#endif


#ifdef MILC_GLOBAL_DEBUG
#ifdef HISQ_FF_DEBUG
  {
  anti_hermitmat ahmat_tmp;
  su3_matrix ttttmat;
printf("Dump W_unitlink[%d][%d]\n",AB_OUT_ON_LINK,AB_OUT_ON_SITE);
dumpmat( &(W_unitlink[AB_OUT_ON_LINK][AB_OUT_ON_SITE]) );
  printf("HISQ FORCE SMEARING 2\n");
  dumpmat( &(force_accum_2[AB_OUT_ON_LINK][AB_OUT_ON_SITE]) );
  printf("HISQ FORCE SMEARING 2 ANTIHERMITIAN\n");
  make_anti_hermitian( &(force_accum_2[AB_OUT_ON_LINK][AB_OUT_ON_SITE]), 
		       &ahmat_tmp );
  uncompress_anti_hermitian( &ahmat_tmp, &(ttttmat) );
  dumpmat( &(ttttmat) );
  }
#endif /* HISQ_FF_DEBUG */
#endif /* MILC_GLOBAL_DEBUG */


  // UNITARIZATION PART

  switch(ap->umethod){
  case UNITARIZE_ROOT:
  case UNITARIZE_RATIONAL:
  case UNITARIZE_ANALYTIC:

    // rephase (out) V, Y
    rephase_field_offset( V_link, OFF, NULL, r0 );
    rephase_field_offset( Y_unitlink, OFF, NULL, r0 );
    
    FORALLFIELDSITES(i){
      for(dir=XUP;dir<=TUP;dir++){
	// derivatives: dY/dV, dY^+/dV
	switch (ap->umethod){
	case UNITARIZE_ROOT:
	  u3_unit_der( &( V_link[4*i+dir] ), &dydv, &dydagdv );
	  break;
	case UNITARIZE_RATIONAL:
	  u3_unit_der_rational( &( V_link[4*i+dir] ), &dydv, &dydagdv );
	  break;
	case UNITARIZE_ANALYTIC:
	  u3_unit_der_analytic( info, &( V_link[4*i+dir] ), &dydv, &dydagdv );
          // +9192 flops
	  break;
	default:
	  node0_printf("Unknown unitarization method\n");
	  terminate(1);
	}
	
  switch(ap->ugroup){
  case UNITARIZE_SU3:
	// piece of derivative from projecting U(3) to SU(3)
	su3_adjoint( &( Y_unitlink[4*i+dir] ), &tmat );
	su3_unit_der_spec( &( V_link[4*i+dir] ), 
			   &( Y_unitlink[4*i+dir] ), &tmat,
			   &dydv, &dydagdv, &dwdv, &dwdagdv );
	  break;
  case UNITARIZE_U3:
	su3t4_copy( &dydv, &dwdv );
	su3t4_copy( &dydagdv, &dwdagdv );
	  break;
	default:
	  node0_printf("Unknown unitarization group. Set UNITARIZE_U3 or UNITARIZE_SU3\n");
	  terminate(1);
	}

	
	// adjoint piece of the force
	su3_adjoint( &( force_accum_2[dir][i] ), &tmat );

	/* LONG COMMENT ABOUT INDEXING CONVENTION
	   f - scalar function, U_{ab} - matrix,
	   derivative df/dU_{ab} is a matrix with indeces transposed,
	   i.e. df/dU_{ab}=D_{ba}
	   In (for example) Wong/Woloshyn notes indeces are NOT transposed,
	   what is called f_{mn} there should be f_{nm}.
	   This code stores f_{mn} (as in notes) because this makes handling
	   smearing part easier (one has several matrix multiplications).
	   However, for the unitarization part matrix derivative is
	   constructed as rank 4-tensor, following usual convention:
	   dV_{kl}/dU_{mn}=T_{knml} (i.e. indeces of the bottom matrix are transposed).
	   Derivative tensor is then contracted with f_{kl} to get f'_{mn}.
	   Usual convention: f'_{nm}=T_{knml}f_{kl},
	   translates into: f'_{mn}=T_{knml}f_{lk} */
	clear_su3mat( &force_tmp );
	for( m=0; m<3; m++) {
	  for( n=0; n<3; n++) {
	    for( k=0; k<3; k++) {
	      for( l=0; l<3; l++) {
		CMUL( dwdv.t4[k][m][n][l], force_accum_2[dir][i].e[l][k], 
		      ftmp );
		CSUM( force_tmp.e[n][m], ftmp );
		/* CAREFUL: adjoint part here! */
		CMUL( dwdagdv.t4[k][m][n][l], tmat.e[l][k], ftmp );
		CSUM( force_tmp.e[n][m], ftmp );
	      }
	    }
	  }
	}
	nflops_all_sites += 81*16;
	// PUT THE RESULT INTO force_accum_2, THIS MAY NEED CHANGES
	su3mat_copy( &force_tmp, &( force_accum_2[dir][i] ) );
      }
    }

    nflops += nflops_all_sites*numnodes()/volume;

//tempflops+=9192+81*16;

    // rephase (in) V, Y
    rephase_field_offset( V_link, ON, NULL, r0 );
    rephase_field_offset( Y_unitlink, ON, NULL, r0 );
    break;

  case UNITARIZE_NONE:
    //  printf("HISQ: Proceeding with UNITARIZE_NONE\n");
    break;
  default:
    printf("HISQ: Only UNITARIZE_ROOT and UNITARIZE_RATIONAL supported so far\n");
    printf("HISQ: Proceeding without reunitarization\n");
  } /* ap->umethod */

#ifdef FFTIME
#ifdef MILC_GLOBAL_DEBUG
  {
  double tmptime=dtime+dclock();
  node0_printf("fftime(hisq-smearing1+reunitarization):  time = %e \n",tmptime);
  }
#endif /* MILC_GLOBAL_DEBUG */
#endif


#ifdef MILC_GLOBAL_DEBUG
#ifdef HISQ_FF_DEBUG
  {
  anti_hermitmat ahmat_tmp;
  su3_matrix ttttmat;
  printf("HISQ FORCE REUNITARIZATION\n");
  dumpmat( &(force_accum_2[AB_OUT_ON_LINK][AB_OUT_ON_SITE]) );
  printf("HISQ FORCE REUNITARIZATION ANTIHERMITIAN\n");
  make_anti_hermitian( &(force_accum_2[AB_OUT_ON_LINK][AB_OUT_ON_SITE]), &ahmat_tmp );
  uncompress_anti_hermitian( &ahmat_tmp, &(ttttmat) );
  dumpmat( &(ttttmat) );
  }
#endif /* HISQ_FF_DEBUG */
#endif /* MILC_GLOBAL_DEBUG */


  //SMEARING ONE PART


  // loop over paths, and loop over links in path 
  // Here each path is one "V link" and each V link should already have a corresponding
  // force_accum_2 computed
  last_netbackdir = NODIR;
  last_path = NULL;
  for( ipath=0; ipath<num_q_paths_1; ipath++ ){
    this_path = &(q_paths_sorted_1[ipath]);
if(this_path->forwback== -1)continue;	// CHECK THIS!!!!

    length = this_path->length;
    // find gather to bring multi_x[term] from "this site" to end of path
    netbackdir = netbackdir_table_1[ipath];
    // and bring force_accum_2 to end - gauge transform one index only!
    // resulting outer product matrix has gauge transformation properties of a connection
    // from start to end of path
// use force_accum_2
// Note if path goes backwards you may need adjoint.  Have to think about
// what its gauge transformation properties are
    if( netbackdir != last_netbackdir){ // don't need to repeat this if same net disp. as last path
	if( GOES_BACKWARDS(netbackdir) ){
	   link_gather_connection_hisq( force_accum_2[OPP_DIR(netbackdir)],
		oprod_along_path[0], NULL, netbackdir );
	}
	else {
	   link_gather_connection_hisq( force_accum_2[netbackdir],
		oprod_along_path[0], NULL, netbackdir );
	}
    }

    /* path transport the outer product force_accum_2(connection) from far end.

       maintain a matrix of the "outer product" transported backwards
	along the path to all sites on the path.
	If new "net displacement", need to completely recreate it.
	Otherwise, use as much of the previous path as possible 

	Note this array is indexed backwards - the outer product transported
	to site number n along the path is in oprod_along_path[length-n].
	This makes reusing it for later paths easier.

	Sometimes we need this at the start point of the path, and sometimes
	one link into the path, so don't always have to do the last link. */

    // figure out how much of the outer products along the path must be
    // recomputed. j is last one needing recomputation. k is first one.
    j=length-1; // default is recompute all
    if( netbackdir == last_netbackdir )
      while ( j>0 && this_path->dir[j] == last_path->dir[j+last_path->length-length] ) j--;
    if( GOES_BACKWARDS(this_path->dir[0]) ) k=1; else k=0;

    for(ilink=j;ilink>=k;ilink--){
      link_transport_connection_hisq( oprod_along_path[length-ilink-1], U_link,
          oprod_along_path[length-ilink], mat_tmp0, this_path->dir[ilink]  );
      nflops += 198;
//tempflops+=9*22;
    }

   /* maintain an array of transports "to this point" along the path.
	Don't recompute beginning parts of path if same as last path */
    ilink=0; // first link where new transport is needed
    if( last_path != NULL )while( this_path->dir[ilink] == last_path->dir[ilink] ) ilink++ ;
    // Sometimes we don't need the matrix for the last link
    if( GOES_FORWARDS(this_path->dir[length-1]) ) k=length-1; else k=length;

    for( ; ilink<k; ilink++ ){
      if( ilink==0 ){
        dir = this_path->dir[0];
          if( GOES_FORWARDS(dir) ){
	    mtag[1] = declare_strided_gather( U_link + dir,
		      aux_stride*sizeof(su3_matrix), sizeof(su3_matrix),
		      OPP_DIR(dir), EVENANDODD, gen_pt[1] );
	    prepare_gather(mtag[1]);
	    do_gather(mtag[1]);
	    //            mtag[1] = start_gather_field( U_link[dir], sizeof(su3_matrix),
	    //                   OPP_DIR(dir), EVENANDODD, gen_pt[1] );
            wait_gather(mtag[1]);
            FORALLFIELDSITES(i){ su3_adjoint( (su3_matrix *)gen_pt[1][i], &(mats_along_path[1][i]) ); }
            cleanup_gather(mtag[1]);
          }
          else{
	    //            FORALLFIELDSITES(i){ mats_along_path[1][i] = U_link[OPP_DIR(dir)][i]; }
            FORALLFIELDSITES(i){ mats_along_path[1][i] = U_link[4*i+OPP_DIR(dir)]; }
          }
      }
      else { // ilink != 0
        dir = OPP_DIR(this_path->dir[ilink]);
        link_transport_connection_hisq( mats_along_path[ilink], U_link,
        mats_along_path[ilink+1], mat_tmp0, dir  );
	nflops += 198;
//tempflops+=9*22;
      }
    } // end loop over links

    /* A path has (length+1) points, counting the ends.  At first
	 point, no "down" direction links have their momenta "at this
	 point". At last, no "up" ... */
    if( GOES_FORWARDS(this_path->dir[length-1]) ) k=length-1; else k=length;
    for( ilink=0; ilink<=k; ilink++ ){
      if(ilink<length)dir = this_path->dir[ilink];
      else dir=NODIR;
      coeff = ferm_epsilon_1*this_path->coeff;
      if( (ilink%2)==1 )coeff = -coeff;

      if(ilink==0 && GOES_FORWARDS(dir) ) FORALLFIELDSITES(i){
	mat_tmp0[i] = oprod_along_path[length][i]; 
      }
      else if( ilink>0){
	FORALLFIELDSITES(i){
        mult_su3_na( &(oprod_along_path[length-ilink][i]),  
		     &(mats_along_path[ilink][i]), &(mat_tmp0[i]) );
	}
	nflops += 198;
      }
//if(ilink>0)tempflops+=9*22;

      /* add in contribution to the force */
      if( ilink<length && GOES_FORWARDS(dir) ){
        FOREVENFIELDSITES(i){
	    scalar_mult_add_su3_matrix( &(force_accum_1[dir][i]), 
					&(mat_tmp0[i]),
					coeff, &(force_accum_1[dir][i]) );
	}
	FORODDFIELDSITES(i){
	    scalar_mult_add_su3_matrix( &(force_accum_1[dir][i]), 
					&(mat_tmp0[i]), -coeff, 
					&(force_accum_1[dir][i]) );
        }
	nflops += 36;
//tempflops+=36;
      }
      if( ilink>0 && GOES_BACKWARDS(lastdir) ){
	odir = OPP_DIR(lastdir);
        FOREVENFIELDSITES(i){
	    scalar_mult_add_su3_matrix( &(force_accum_1[odir][i]), 
					&(mat_tmp0[i]), -coeff, 
					&(force_accum_1[odir][i]) );
	}
        FORODDFIELDSITES(i){
	    scalar_mult_add_su3_matrix( &(force_accum_1[odir][i]), 
					&(mat_tmp0[i]), coeff, 
					&(force_accum_1[odir][i]) );
	}
	nflops += 36;
//tempflops+=36;
      }

      lastdir = dir;
    } /* end loop over links in path */
    last_netbackdir = netbackdir;
    last_path = &(q_paths_sorted_1[ipath]);
  } /* end loop over paths */

  // END SMEARING ONE PART

#ifdef MILC_GLOBAL_DEBUG
#ifdef HISQ_REUNITARIZATION_DEBUG
  {
  anti_hermitmat ahmat_tmp;
  su3_matrix ff_interm;
  double icounter=0;
  double max_force=0;
  Real force_interm;
  FORALLFIELDSITES(i) {
    for( dir=XUP; dir<=TUP; dir++ ) {
      if( fabs(s->RoS[dir])>1.0 ) icounter+=1;

      /* find maximum force */
      make_anti_hermitian( &(force_accum_1[dir][i]), &ahmat_tmp );
      uncompress_anti_hermitian( &ahmat_tmp, &ff_interm );
      force_interm = su3_norm_frob( &ff_interm );
      if( force_interm > max_force ) {
        max_force=force_interm;
      }
    }
  }
  g_doublemax( &max_force );
  g_doublesum( &icounter );
  node0_printf( "REUNITARIZATION: fabs(RoS) > 1.0 occured on %lf link(s)\n", icounter );
  node0_printf( "HISQ FORCE: maximum norm is %lf\n", max_force );
  }
#endif /* HISQ_REUNITARIZATION_DEBUG */
#ifdef HISQ_FF_DEBUG
  {
  anti_hermitmat ahmat_tmp;
  printf("HISQ FORCE\n");
  dumpmat( &(force_accum_1[AB_OUT_ON_LINK][AB_OUT_ON_SITE]) );
  printf("HISQ FORCE ANTIHERMITIAN\n");
  for(dir=XUP; dir<=TUP; dir++) FORALLFIELDSITES(i){
    make_anti_hermitian( &(force_accum_1[dir][i]), &ahmat_tmp );
    uncompress_anti_hermitian( &ahmat_tmp, &(force_accum_1[dir][i]) );
  }
  dumpmat( &(force_accum_1[AB_OUT_ON_LINK][AB_OUT_ON_SITE]) );
  }
#endif /* HISQ_FF_DEBUG */
#endif /* MILC_GLOBAL_DEBUG */


  /* Put antihermitian traceless part into momentum */
  // add force to momentum
  for(dir=XUP; dir<=TUP; dir++)FORALLSITES(i,s){
     uncompress_anti_hermitian( &(s->mom[dir]), &tmat2 );
     add_su3_matrix( &tmat2, &(force_accum_1[dir][i]), &tmat2 );
     make_anti_hermitian( &tmat2, &(s->mom[dir]) );
  }
  nflops += 144;
//tempflops+=4*18;
//tempflops+=4*18;
	
  free( mat_tmp0 );
  free( mat_tmp1 );
  for(i=0;i<=MAX_PATH_LENGTH;i++){
     free( oprod_along_path[i] );
  }
  for(i=1;i<=MAX_PATH_LENGTH;i++){
     free( mats_along_path[i] );
  }
  for(i=XUP;i<=TUP;i++){
     free( force_accum_1[i] );
     free( force_accum_2[i] );
     free( force_accum_tmp[i] );
  }

  info->final_flop = ((double)nflops)*volume/numnodes();
  dtime += dclock();
  info->final_sec = dtime;

//printf("FF flops = %d\n",tempflops);
} //fn_fermion_force_multi_hisq_mx



#if defined(HISQ_FF_MULTI_WRAPPER) || defined(USE_FF_GPU)

#ifndef USE_FF_GPU

/********************************************************************/
/* Non-GPU wrapper support                                          */
/********************************************************************/

// HISQ wrapper routine: contribution to the force due to
//    0) smearing level 0 (outer product |X><Y|)
//    1) smearing level 2 + smearing level 2 Naik corrected
//    2) reunitarization
//    3) smearing level 1

static void 
fn_fermion_force_multi_hisq_wrapper_mx_cpu( info_t *info, Real eps, Real *residues, 
				     su3_vector **multi_x, int nterms,
				     fermion_links_t *fl)
{
  // Get the terms we need from the fermion links structure */
  ks_action_paths_hisq *ap = get_action_paths_hisq(fl);
  hisq_auxiliary_t *aux = get_hisq_auxiliary(fl);
  int n_naiks = fermion_links_get_n_naiks(fl);
  double *eps_naik = fermion_links_get_eps_naik(fl);
  su3_matrix *U_link = aux->U_link;
  su3_matrix *W_unitlink = aux->W_unitlink;
  su3_matrix *force_accum_0[4];  // accumulate force, smearing 0
  su3_matrix *force_accum_0_naik[4];  // accumulate force, smearing 0, Naik
  su3_matrix *force_accum_1[4];  // accumulate force, smearing 1
  su3_matrix *force_accum_1u[4];  // accumulate force, reunitarization
  su3_matrix *force_accum_2[4];  // accumulate force, smearing 2
  su3_matrix *force_final[4];  // accumulate final force, multiply by time step
  su3_matrix tmat2;
  anti_hermitmat *ahmat_tmp;
  register site *s;
  int dir;
  int i;
  int ipath;
  int inaik;
  int n_naik_shift;
  int num_q_paths_current,n_orders_naik_current;
  Real coeff_mult;
  size_t nflops = 0;
  double dtime = -dclock();
  double final_flop = 0;


//  node0_printf("Entering fn_fermion_force_multi_hisq_wrapper_mx()\n");
//  for(i=0;i<n_naiks;i++)
//    node0_printf("wrapper_mx: orders[%d]=%d\n",i,n_orders_naik[i]);
//  node0_printf("wrapper_mx: nterms=%d\n",nterms);

  int num_q_paths_1 = ap->p1.num_q_paths;
  Q_path *q_paths_1 = ap->p1.q_paths;
  int num_q_paths_2 = ap->p2.num_q_paths;
  Q_path *q_paths_2 = ap->p2.q_paths;
  int num_q_paths_3 = ap->p3.num_q_paths;
  Q_path *q_paths_3 = ap->p3.q_paths;


  if( first_force==1 ){
    if( q_paths_sorted_1==NULL ) 
      q_paths_sorted_1 = (Q_path *)malloc( num_q_paths_1*sizeof(Q_path) );
    if(netbackdir_table_1==NULL ) 
      netbackdir_table_1 = (int *)malloc( num_q_paths_1*sizeof(int) );
    if( q_paths_sorted_2==NULL ) 
      q_paths_sorted_2 = (Q_path *)malloc( num_q_paths_2*sizeof(Q_path) );
    if(netbackdir_table_2==NULL ) 
      netbackdir_table_2 = (int *)malloc( num_q_paths_2*sizeof(int) );
    if( q_paths_sorted_3==NULL ) 
      q_paths_sorted_3 = (Q_path *)malloc( num_q_paths_3*sizeof(Q_path) );
    if(netbackdir_table_3==NULL ) 
      netbackdir_table_3 = (int *)malloc( num_q_paths_3*sizeof(int) );
    else{ node0_printf("WARNING: remaking sorted path tables\n"); exit(0); }
    // make sorted tables
    sort_quark_paths( q_paths_1, q_paths_sorted_1, num_q_paths_1, 8 );
    for( ipath=0; ipath<num_q_paths_1; ipath++ )
      netbackdir_table_1[ipath] = 
	find_backwards_gather( &(q_paths_sorted_1[ipath]) );
    sort_quark_paths( q_paths_2, q_paths_sorted_2, num_q_paths_2, 16 );
    for( ipath=0; ipath<num_q_paths_2; ipath++ )
      netbackdir_table_2[ipath] = 
	find_backwards_gather( &(q_paths_sorted_2[ipath]) );
    sort_quark_paths( q_paths_3, q_paths_sorted_3, num_q_paths_3, 16 );
    for( ipath=0; ipath<num_q_paths_3; ipath++ )
      netbackdir_table_3[ipath] = 
	find_backwards_gather( &(q_paths_sorted_3[ipath]) );
    first_force=0;
  }


  for(i=XUP;i<=TUP;i++){
     force_accum_0[i] = 
       (su3_matrix *) malloc(sites_on_node*sizeof(su3_matrix) );
     force_accum_0_naik[i] = 
       (su3_matrix *) malloc(sites_on_node*sizeof(su3_matrix) );
     force_accum_1[i] = 
       (su3_matrix *) malloc(sites_on_node*sizeof(su3_matrix) );
     force_accum_1u[i] = 
       (su3_matrix *) malloc(sites_on_node*sizeof(su3_matrix) );
     force_accum_2[i] = 
       (su3_matrix *) malloc(sites_on_node*sizeof(su3_matrix) );
     force_final[i] = 
       (su3_matrix *) malloc(sites_on_node*sizeof(su3_matrix) );
  }
  ahmat_tmp = 
    (anti_hermitmat *) special_alloc(sites_on_node*sizeof(anti_hermitmat) );

  //  node0_printf("UNITARIZATION_METHOD=%d\n",ap->umethod);

  for(dir=XUP;dir<=TUP;dir++) {
    FORALLFIELDSITES(i) {
      clear_su3mat( &( force_accum_2[dir][i] ) );
    }
  }

  // loop on different naik epsilons
  n_naik_shift = 0;
  for( inaik=0; inaik<n_naiks; inaik++ ) {
#ifdef MILC_GLOBAL_DEBUG
    node0_printf("wrapper_mx: eps_naik[%d]=%f\n",inaik,eps_naik[inaik]);
#endif /* MILC_GLOBAL_DEBUG */

    // smearing level 0
    if( 0==inaik ) {
      n_orders_naik_current = n_order_naik_total;
    }
    else {
      n_orders_naik_current = n_orders_naik[inaik];
    }
    fn_fermion_force_multi_hisq_smearing0(info, eps, residues+n_naik_shift, 
	 multi_x+n_naik_shift, n_orders_naik_current, force_accum_0, 
	   force_accum_0_naik );
    final_flop += info->final_flop;
#ifdef MILC_GLOBAL_DEBUG
#ifdef HISQ_FF_DEBUG
  {
  printf("WRAPPER FORCE EXP level 0\n");
  dumpmat( &(force_accum_0[AB_OUT_ON_LINK][AB_OUT_ON_SITE]) );
  printf("WRAPPER FORCE ANTIHERMITIAN EXP level 0\n");
  for(dir=XUP; dir<=TUP; dir++) FORALLFIELDSITES(i){
    make_anti_hermitian( &(force_accum_0[dir][i]), &(ahmat_tmp[i]) );
    uncompress_anti_hermitian( &(ahmat_tmp[i]), &(force_accum_1u[dir][i]) );
  }
  dumpmat( &(force_accum_1u[AB_OUT_ON_LINK][AB_OUT_ON_SITE]) );
  su3_matrix su3ab, su3ab2;
  clear_su3mat( &su3ab );
  clear_su3mat( &su3ab2 );
  for(dir=XUP; dir<=TUP; dir++) FORALLFIELDSITES(i){
    add_su3_matrix( &su3ab, &(force_accum_0[dir][i]), &su3ab );
    add_su3_matrix( &su3ab, &(force_accum_0_naik[dir][i]), &su3ab );
    add_su3_matrix( &su3ab2, &( W_unitlink[4*i+dir] ), &su3ab2 );
  }
  printf("WRAPPER FORCE EXP level 0 GLOBAL SUM\n");
  dumpmat( &su3ab );
  printf("WRAPPER FORCE EXP W_unitlink GLOBAL SUM\n");
  dumpmat( &su3ab2 );
  }
#endif /* HISQ_FF_DEBUG */
#endif /* MILC_GLOBAL_DEBUG */
    // smearing level 2
    if( 0==inaik ) {
      q_paths_sorted_current = q_paths_sorted_2;
      num_q_paths_current = num_q_paths_2;
      netbackdir_table_current = netbackdir_table_2;
    }
    else {
      q_paths_sorted_current = q_paths_sorted_3;
      num_q_paths_current = num_q_paths_3;
      netbackdir_table_current = netbackdir_table_3;
    }
    fn_fermion_force_multi_hisq_smearing( info, eps, residues+n_naik_shift, 
        multi_x+n_naik_shift, n_orders_naik_current, force_accum_1, 
        force_accum_0, force_accum_0_naik, W_unitlink, 
	num_q_paths_current, q_paths_sorted_current, 
        netbackdir_table_current );
    final_flop += info->final_flop;
#ifdef MILC_GLOBAL_DEBUG
#ifdef HISQ_FF_DEBUG
  {
  printf("WRAPPER FORCE EXP level 2\n");
  dumpmat( &(force_accum_1[AB_OUT_ON_LINK][AB_OUT_ON_SITE]) );
  printf("WRAPPER FORCE ANTIHERMITIAN EXP level 2\n");
  for(dir=XUP; dir<=TUP; dir++) FORALLFIELDSITES(i){
    make_anti_hermitian( &(force_accum_1[dir][i]), &(ahmat_tmp[i]) );
    uncompress_anti_hermitian( &(ahmat_tmp[i]), &(force_accum_1u[dir][i]) );
  }
  dumpmat( &(force_accum_1u[AB_OUT_ON_LINK][AB_OUT_ON_SITE]) );
  }
#endif /* HISQ_FF_DEBUG */
#endif /* MILC_GLOBAL_DEBUG */
    if( 0==inaik ) {
      coeff_mult = 1.0;
    }
    else {
      coeff_mult = eps_naik[inaik];
    }
    for(dir=XUP;dir<=TUP;dir++) {
      FORALLFIELDSITES(i) {
//        add_su3_matrix( &( force_accum_2[dir][i] ), &( force_accum_1[dir][i] ), &( force_accum_2[dir][i] ) );
        scalar_mult_add_su3_matrix( &( force_accum_2[dir][i] ), 
          &( force_accum_1[dir][i] ), coeff_mult, &( force_accum_2[dir][i] ) );
      }
      nflops += 36;
    }
    n_naik_shift += n_orders_naik[inaik];
  }


  switch(ap->umethod){
  case UNITARIZE_NONE:
    // smearing level 1
    fn_fermion_force_multi_hisq_smearing( info, eps, residues, multi_x, nterms, 
        force_accum_1, force_accum_2, NULL,
        U_link, num_q_paths_1, q_paths_sorted_1, netbackdir_table_1 );
    final_flop += info->final_flop;
    break;
  case UNITARIZE_ROOT:
  case UNITARIZE_RATIONAL:
  case UNITARIZE_ANALYTIC:
    // reunitarization
    fn_fermion_force_multi_hisq_reunit(info, force_accum_1u, force_accum_2, fl );
    final_flop += info->final_flop;
#ifdef MILC_GLOBAL_DEBUG
#ifdef HISQ_FF_DEBUG
  {
  printf("WRAPPER FORCE EXP REUNIT\n");
  dumpmat( &(force_accum_1u[AB_OUT_ON_LINK][AB_OUT_ON_SITE]) );
  }
#endif /* HISQ_FF_DEBUG */
#endif /* MILC_GLOBAL_DEBUG */
    // smearing level 1
    fn_fermion_force_multi_hisq_smearing(info, eps, residues, multi_x, 
          nterms, force_accum_1, force_accum_1u, NULL, 
	  U_link, num_q_paths_1, q_paths_sorted_1, netbackdir_table_1 );
    final_flop += info->final_flop;
    break;
  default:
    node0_printf("Unknown unitarization method\n");
    terminate(1);
  }

  // contraction with the link in question should be done here,
  // after contributions from all levels of smearing are taken into account
  for(dir=XUP;dir<=TUP;dir++){
    FORALLSITES(i,s){
      mult_su3_nn( &( s->link[dir] ), &( force_accum_1[dir][i] ), 
		   &( force_final[dir][i] ) );
    }
    nflops += 198;
  }

  // multiply by the time step "eps" twice (only forward paths are considered)
  // take into account even/odd parity (it is NOT done in "smearing" routine)
  for(dir=XUP;dir<=TUP;dir++){
    FOREVENFIELDSITES(i){
      scalar_mult_su3_matrix( &(force_final[dir][i]), 2.0*eps, 
			      &(force_final[dir][i]) );
    }
    FORODDFIELDSITES(i){
      scalar_mult_su3_matrix( &(force_final[dir][i]), -2.0*eps, 
			      &(force_final[dir][i]) );
    }
    nflops += 36;
  }

#ifdef MILC_GLOBAL_DEBUG
#ifdef HISQ_FF_DEBUG
  {
  printf("WRAPPER FORCE\n");
  dumpmat( &(force_final[AB_OUT_ON_LINK][AB_OUT_ON_SITE]) );
  printf("WRAPPER FORCE ANTIHERMITIAN\n");
  for(dir=XUP; dir<=TUP; dir++) FORALLFIELDSITES(i){
    make_anti_hermitian( &(force_final[dir][i]), &(ahmat_tmp[i]) );
    uncompress_anti_hermitian( &(ahmat_tmp[i]), &(force_final[dir][i]) );
  }
  dumpmat( &(force_final[AB_OUT_ON_LINK][AB_OUT_ON_SITE]) );
  }
#endif /* HISQ_FF_DEBUG */
#endif /* MILC_GLOBAL_DEBUG */

  /* Put antihermitian traceless part into momentum */
  // add force to momentum
  for(dir=XUP; dir<=TUP; dir++)FORALLSITES(i,s){
     uncompress_anti_hermitian( &(s->mom[dir]), &tmat2 );
     add_su3_matrix( &tmat2, &(force_final[dir][i]), &tmat2 );
     make_anti_hermitian( &tmat2, &(s->mom[dir]) );
  }

  nflops += 18;


  for(i=XUP;i<=TUP;i++){
     free( force_accum_0[i] );
     free( force_accum_0_naik[i] );
     free( force_accum_1[i] );
     free( force_accum_1u[i] );
     free( force_accum_2[i] );
     free( force_final[i] );
  }
  free(ahmat_tmp);

  dtime += dclock();
  info->final_sec = dtime;
  info->final_flop = final_flop + ((double)nflops)*volume/numnodes();

} //fn_fermion_force_multi_hisq_wrapper_mx


// Contribution to the force from 0 level of smearing, force contains only |X><Y| terms
static void 
fn_fermion_force_multi_hisq_smearing0( info_t *info, Real eps, Real *residues, 
       su3_vector **multi_x, int nterms,
       su3_matrix *force_accum[4], su3_matrix *force_accum_naik[4] )
{
  /* note CG_solution and Dslash * solution are combined in "multi_x" */
  /* New version 1/21/99.  Use forward part of Dslash to get force */
  /* see long comment at end */
  /* For each link we need multi_x transported from both ends of path. */
  int term;
  register int i,k;
  int dir;
  su3_matrix tmat;
  Real coeff;
  size_t nflops = 0;

  msg_tag *mtag[2];
  su3_matrix *mat_tmp0;
  // length N path has N+1 sites!!
  su3_matrix *oprod_along_path[MAX_PATH_LENGTH+1];


  if( nterms==0 )return;

#ifdef MILC_GLOBAL_DEBUG
#ifdef HISQ_FF_DEBUG
  {
  printf("WRAPPER FORCE smearing0 pseudofermion vector at site (0,0,0,0)\n");
  dumpvec( &(multi_x[0][0]) );
  printf("\n");
  su3_vector su3ab;
  clearvec( &su3ab );
  int nnn;
  for( nnn=0; nnn<nterms;nnn++) {
    FORALLFIELDSITES(i){
      add_su3_vector( &su3ab, &(multi_x[nnn][i]), &su3ab );
    }
    nflops += 6;
  }
  printf("WRAPPER FORCE smearing0 pseudofermion vector GLOBAL SUM\n");
  dumpvec( &su3ab );
  printf("\n");
  for( nnn=0; nnn<nterms;nnn++) {
      printf("residues[%d]=%lf\n",nnn,residues[nnn]);
    }
  }
#endif /* HISQ_FF_DEBUG */
#endif /* MILC_GLOBAL_DEBUG */


  for(i=0;i<=MAX_PATH_LENGTH;i++){
     oprod_along_path[i] = 
       (su3_matrix *) malloc(sites_on_node*sizeof(su3_matrix) );
  }
  mat_tmp0 = (su3_matrix *) special_alloc(sites_on_node*sizeof(su3_matrix) );
  if( mat_tmp0 == NULL ){printf("Node %d NO ROOM\n",this_node); exit(0); }



  // clear force accumulators
  for(dir=XUP;dir<=TUP;dir++)
    FORALLFIELDSITES(i)clear_su3mat( &(force_accum[dir][i]) );

  //AB loop on directions, path table is not needed
  for(dir=XUP;dir<=TUP;dir++){
    k=0; // which gather we are using
    //AB netbackdir is just the opposite of dir for 1x1 case
    mtag[k] = start_gather_field( multi_x[0], sizeof(su3_vector),
      OPP_DIR(dir), EVENANDODD, gen_pt[k] );
    FORALLFIELDSITES(i){
      clear_su3mat(  &(oprod_along_path[0][i]) ); // actually last site in path
    }
    for(term=0;term<nterms;term++){
      if(term<nterms-1)mtag[1-k] = start_gather_field( multi_x[term+1],
      sizeof(su3_vector), OPP_DIR(dir), EVENANDODD, gen_pt[1-k] );
      wait_gather(mtag[k]);

      FORALLFIELDSITES(i){
        // build projector as usual
        su3_projector( &multi_x[term][i], (su3_vector *)gen_pt[k][i], &tmat );
        // multiply by alpha_l in rational function expansion
        scalar_mult_add_su3_matrix( &(oprod_along_path[0][i]), &tmat, 
				    residues[term], &(oprod_along_path[0][i]) );
      }
      nflops += 54 + 36;
      cleanup_gather(mtag[k]);
      k=1-k; // swap 0 and 1
    } /* end loop over terms in rational function expansion */


    link_gather_connection_hisq( oprod_along_path[0], oprod_along_path[1], 
				 mat_tmp0, dir );

    coeff = 1; // fermion_eps is outside this routine in "wrapper" routine
    FORALLFIELDSITES(i){
      scalar_mult_add_su3_matrix( &(force_accum[dir][i]), 
				  &(oprod_along_path[1][i]),
				  coeff, &(force_accum[dir][i]) );
    }
    nflops += 36;

  } /* end of loop on directions */


  /* *** Naik part *** */
  // clear force accumulators
  for(dir=XUP;dir<=TUP;dir++)
    FORALLFIELDSITES(i)clear_su3mat( &(force_accum_naik[dir][i]) );

  for(dir=XUP;dir<=TUP;dir++){ //AB loop on directions, path table is not needed
    k=0; // which gather we are using
    mtag[k] = start_gather_field( multi_x[0], sizeof(su3_vector),
      OPP_3_DIR( DIR3(dir) ), EVENANDODD, gen_pt[k] );
    FORALLFIELDSITES(i){
      clear_su3mat(  &(oprod_along_path[0][i]) );
    }
    for(term=0;term<nterms;term++){
      if(term<nterms-1)mtag[1-k] = start_gather_field( multi_x[term+1],
      sizeof(su3_vector), OPP_3_DIR( DIR3(dir) ), EVENANDODD, gen_pt[1-k] );
      wait_gather(mtag[k]);

      FORALLFIELDSITES(i){
        // build projector as usual
        su3_projector( &multi_x[term][i], (su3_vector *)gen_pt[k][i], &tmat );
        // multiply by alpha_l in rational function expansion
        scalar_mult_add_su3_matrix( &(oprod_along_path[0][i]), &tmat, 
				    residues[term], &(oprod_along_path[0][i]) );
      }
      nflops += 54 + 36;
      cleanup_gather(mtag[k]);
      k=1-k; // swap 0 and 1
    } /* end loop over terms in rational function expansion */

    link_gather_connection_hisq( oprod_along_path[0], oprod_along_path[1], 
				 mat_tmp0, DIR3(dir) );

    coeff = 1; // fermion_eps is outside this routine in "wrapper" routine
    FORALLFIELDSITES(i){
      scalar_mult_add_su3_matrix( &(force_accum_naik[dir][i]), 
				  &(oprod_along_path[1][i]),
				  coeff, &(force_accum_naik[dir][i]) );
    }
    nflops += 36;

  } /* end of loop on directions */


  free( mat_tmp0 );
  for(i=0;i<=MAX_PATH_LENGTH;i++){
     free( oprod_along_path[i] );
  }

  info->final_flop = ((double)nflops)*volume/numnodes();
} //fn_fermion_force_multi_hisq_smearing0

// This routine calculates contribution to the force due to i-th level
// of smearing,
// as input it uses array of new links and the force from (i-1)-th
// level of smearing
static void 
fn_fermion_force_multi_hisq_smearing( info_t *info,
				      Real eps, 
				      Real *residues, 
				      su3_vector **multi_x, 
				      int nterms, 
				      su3_matrix *force_accum[4], 
				      su3_matrix *force_accum_old[4], 
				      su3_matrix *force_accum_naik_old[4], 
				      su3_matrix *internal_U_link,
				      int internal_num_q_paths, 
				      Q_path *internal_q_paths_sorted, 
				      int *internal_netbackdir_table)
{
  /* note CG_solution and Dslash * solution are combined in "multi_x" */
  /* New version 1/21/99.  Use forward part of Dslash to get force */
  /* see long comment at end */
  /* For each link we need multi_x transported from both ends of path. */
  register int i,j,k,lastdir=-99,ipath,ilink;
  int length,dir,odir;
  Real coeff;
  Q_path *this_path;	// pointer to current path
  msg_tag *mtag[2];
  su3_matrix *mat_tmp0;
  su3_matrix *mat_tmp1;
  anti_hermitmat *ahmat_tmp;
  su3_matrix *oprod_along_path[MAX_PATH_LENGTH+1]; // length N path has N+1 sites!!
  su3_matrix *mats_along_path[MAX_PATH_LENGTH+1]; // 
  int netbackdir;	// backwards direction for entire path
  int aux_stride = 4;
  size_t nflops = 0;

  if( nterms==0 )return;

  for(i=0;i<=MAX_PATH_LENGTH;i++){
     oprod_along_path[i] = (su3_matrix *) malloc(sites_on_node*sizeof(su3_matrix) );
  }
  for(i=1;i<=MAX_PATH_LENGTH;i++){ // 0 element is never used (it's unit matrix)
     mats_along_path[i] = 
       (su3_matrix *) malloc(sites_on_node*sizeof(su3_matrix) );
  }

  mat_tmp0 = (su3_matrix *) special_alloc(sites_on_node*sizeof(su3_matrix) );
  mat_tmp1 = (su3_matrix *) special_alloc(sites_on_node*sizeof(su3_matrix) );
  if( mat_tmp0 == NULL ){printf("Node %d NO ROOM\n",this_node); exit(0); }
  if( mat_tmp1 == NULL ){printf("Node %d NO ROOM\n",this_node); exit(0); }

  ahmat_tmp = (anti_hermitmat *) 
    special_alloc(sites_on_node*sizeof(anti_hermitmat) );

	

#ifdef MILC_GLOBAL_DEBUG
#ifdef HISQ_FF_DEBUG
  {
  su3_matrix su3ab, su3ab2;
  clear_su3mat( &su3ab );
  clear_su3mat( &su3ab2 );
  for(dir=XUP; dir<=TUP; dir++) {
    FORALLFIELDSITES(i){
      add_su3_matrix( &su3ab, &(force_accum_old[dir][i]), &su3ab );
      if( NULL!=force_accum_naik_old ) {
	add_su3_matrix( &su3ab, &(force_accum_naik_old[dir][i]), &su3ab );
      }
      add_su3_matrix( &su3ab2, &( internal_U_link[4*i+dir] ), &su3ab2 );
    }
  }
  nflops += 36*4;
  if( NULL!=force_accum_naik_old )nflops += 18*4;
  printf("WRAPPER FORCE smearing input force GLOBAL SUM\n");
  dumpmat( &su3ab );
  printf("WRAPPER FORCE smearing input links GLOBAL SUM\n");
  dumpmat( &su3ab2 );
  }
#endif /* HISQ_FF_DEBUG */
#endif /* MILC_GLOBAL_DEBUG */


  // clear force accumulators
  for(dir=XUP;dir<=TUP;dir++)
    FORALLFIELDSITES(i)clear_su3mat( &(force_accum[dir][i]) );

  /* loop over paths, and loop over links in path */
  for( ipath=0; ipath<internal_num_q_paths; ipath++ ){
    this_path = &(internal_q_paths_sorted[ipath]);
    if(this_path->forwback== -1)continue;	/* skip backwards dslash */

    length = this_path->length;
    netbackdir = internal_netbackdir_table[ipath];

    /* move f(i-1) force from current site in positive direction.  This
       corresponds to outer product |X><Y| calculated at the endpoint
       of the path */
    if( netbackdir<8) { // Not a Naik path
      link_gather_connection_hisq( force_accum_old[OPP_DIR(netbackdir)], 
				   oprod_along_path[0], NULL, netbackdir );
    }
    else { // Naik path
      if( NULL==force_accum_naik_old ) {
        node0_printf( "fn_fermion_force_multi_hisq_smearing:  mismatch:\n" );
        node0_printf( "  force_accum_naik_old is NULL, but path table contains Naik paths(!)\n" );
        terminate(1);
      }
      // CONVERSION FROM 3-LINK DIRECTION TO 1-LINK DIRECTION
      link_gather_connection_hisq( force_accum_naik_old[OPP_DIR(netbackdir-8)], 
				   oprod_along_path[0], NULL, netbackdir );
    }


    // figure out how much of the outer products along the path must be
    // recomputed. j is last one needing recomputation. k is first one.
    j=length-1; // default is recompute all
    if( GOES_BACKWARDS(this_path->dir[0]) ) k=1; else k=0;


    for(ilink=j;ilink>=k;ilink--){
        link_transport_connection_hisq( oprod_along_path[length-ilink-1], 
		internal_U_link, oprod_along_path[length-ilink], mat_tmp0, 
		this_path->dir[ilink]  );
	nflops += 198;
    }

   /* maintain an array of transports "to this point" along the path.
	Don't recompute beginning parts of path if same as last path */
    ilink=0; // first link where new transport is needed
    // Sometimes we don't need the matrix for the last link
    if( GOES_FORWARDS(this_path->dir[length-1]) ) k=length-1; else k=length;

    for( ; ilink<k; ilink++ ){
      if( ilink==0 ){
        dir = this_path->dir[0];
          if( GOES_FORWARDS(dir) ){
	    mtag[1] = declare_strided_gather( internal_U_link + dir,
		      aux_stride*sizeof(su3_matrix), sizeof(su3_matrix),
		      OPP_DIR(dir), EVENANDODD, gen_pt[1] );
	    prepare_gather(mtag[1]);
	    do_gather(mtag[1]);
	    //            mtag[1] = start_gather_field( internal_U_link[dir], 
	    //					  sizeof(su3_matrix),
	    //					  OPP_DIR(dir), EVENANDODD, gen_pt[1] );
            wait_gather(mtag[1]);
            FORALLFIELDSITES(i){ su3_adjoint( (su3_matrix *)gen_pt[1][i], 
					   &(mats_along_path[1][i]) ); }
            cleanup_gather(mtag[1]);
          }
          else{
//            FORALLSITES(i,s){ mats_along_path[1][i] = s->link[OPP_DIR(dir)]; }
            FORALLFIELDSITES(i){ 
	      mats_along_path[1][i] = internal_U_link[4*i+OPP_DIR(dir)]; 
	    }
          }
      }
      else { // ilink != 0
        dir = OPP_DIR(this_path->dir[ilink]);
        link_transport_connection_hisq( mats_along_path[ilink], 
		internal_U_link, mats_along_path[ilink+1], mat_tmp0, dir  );
	nflops += 198;
      }
    } // end loop over links



    /* A path has (length+1) points, counting the ends.  At first
	 point, no "down" direction links have their momenta "at this
	 point". At last, no "up" ... */
    if( GOES_FORWARDS(this_path->dir[length-1]) ) k=length-1; else k=length;
    for( ilink=0; ilink<=k; ilink++ ){
      if(ilink<length)dir = this_path->dir[ilink];
      else dir=NODIR;

      coeff = this_path->coeff;
      if( (ilink%2)==1 )coeff = -coeff;

      /* add in contribution to the force */
      if( ilink<length && GOES_FORWARDS(dir) ){

        link_gather_connection_hisq( oprod_along_path[length-ilink-1], 
				     mat_tmp1, NULL, dir );

        if(ilink==0) FORALLFIELDSITES(i){
          mat_tmp0[i] = mat_tmp1[i]; 
        }
        else{
          FORALLFIELDSITES(i){
            mult_su3_na( &( mat_tmp1[i] ),  &( mats_along_path[ilink][i] ), &(mat_tmp0[i]) );
          }
	  nflops += 198;
        }
        FORALLFIELDSITES(i){
	    scalar_mult_add_su3_matrix( &(force_accum[dir][i]), &(mat_tmp0[i]),
                coeff, &(force_accum[dir][i]) );
	}
	nflops += 36;
      }
      if( ilink>0 && GOES_BACKWARDS(lastdir) ){
	odir = OPP_DIR(lastdir);
        if( ilink==1 )FORALLFIELDSITES(i){
          mat_tmp0[i]=oprod_along_path[length-ilink][i];
          su3_adjoint ( &( mat_tmp0[i] ), &( mat_tmp1[i] ) );
        }
        else{
          link_gather_connection_hisq( mats_along_path[ilink-1], mat_tmp1, 
				       NULL, odir );
          FORALLFIELDSITES(i){
            mult_su3_na( &( oprod_along_path[length-ilink][i] ), 
			 &( mat_tmp1[i] ), &( mat_tmp0[i] ) );
            su3_adjoint ( &( mat_tmp0[i] ), &( mat_tmp1[i] ) );
          }
	  nflops += 198;
        }
        FORALLFIELDSITES(i){
          scalar_mult_add_su3_matrix( &(force_accum[odir][i]), &(mat_tmp1[i]),
                +coeff, &(force_accum[odir][i]) );
        }

	nflops += 36;
      }
      lastdir = dir;
    } /* end loop over links in path */
  } /* end loop over paths */

  free( mat_tmp0 );
  free( mat_tmp1 );
  free( ahmat_tmp );
  for(i=0;i<=MAX_PATH_LENGTH;i++){
     free( oprod_along_path[i] );
  }
  for(i=1;i<=MAX_PATH_LENGTH;i++){
     free( mats_along_path[i] );
  }

  info->final_flop = ((double)nflops)*volume/numnodes();
}// fn_fermion_force_multi_hisq_smearing


// Contribution to the force due to reunitarization of links:
//   Y=V*(V^+*V)^-1/2, W=Y/(detY)^1/3
static void 
fn_fermion_force_multi_hisq_reunit( info_t *info, su3_matrix *force_accum[4], 
				    su3_matrix *force_accum_old[4], 
				    fermion_links_t *fl)
{
  ks_action_paths_hisq *ap = get_action_paths_hisq(fl);
  hisq_auxiliary_t *aux = get_hisq_auxiliary(fl);
  su3_matrix *internal_V_link = aux->V_link;
  su3_matrix *internal_Y_link = aux->Y_unitlink;
  int i, k, l, m, n, dir;
  complex ftmp;
  su3_matrix tmat;
  su3_tensor4 dwdv, dwdagdv;
  su3_tensor4 dwdvs, dwdagdvs;
  /* Assume that for the fermion force calculation the KS phases are never
     defined for a nonzero origin */
  int r0[4] = {0, 0, 0, 0};
  size_t nflops_all_sites = 0;

  info->final_flop = 0.0;
  
  /* phases out */
  rephase_field_offset( internal_V_link, OFF, NULL , r0);
  rephase_field_offset( internal_Y_link, OFF, NULL , r0);
  
  FORALLFIELDSITES(i){
    for(dir=XUP;dir<=TUP;dir++){
      // calculate derivatives
      if(ap->umethod == UNITARIZE_ROOT )
	u3_unit_der( &( internal_V_link[4*i+dir] ), &dwdv, &dwdagdv );
      else if(ap->umethod == UNITARIZE_RATIONAL )
	u3_unit_der_rational( &( internal_V_link[4*i+dir] ), &dwdv, &dwdagdv );
      else if(ap->umethod == UNITARIZE_ANALYTIC ) {
	u3_unit_der_analytic( info, &( internal_V_link[4*i+dir] ), &dwdv, &dwdagdv );
#ifdef MILC_GLOBAL_DEBUG
#ifdef HISQ_FF_DEBUG
  {
    if( i==0 ) {
      printf("WRAPPER FORCE REUNIT V link at site 0\n");
      dumpmat( &(internal_V_link[4*i+dir]) );
      printf("WRAPPER FORCE REUNIT derivative at site 0\n");
      dumpt4( dwdv );
      dumpt4( dwdagdv );
    }
  }
#endif /* HISQ_FF_DEBUG */
#endif /* MILC_GLOBAL_DEBUG */
      }
      else if(ap->umethod == UNITARIZE_NONE ) {
	node0_printf("fn_fermion_force_multi_hisq_reunit should not be called for UNITARIZE_NONE\n");
	terminate(1);
      }
      else{
	node0_printf("Unknown unitarization method\n");
	terminate(1);
      }

  switch(ap->ugroup){
  case UNITARIZE_SU3:
      /* piece of derivative from projecting U(3) to SU(3) */
      su3_adjoint( &( internal_Y_link[4*i+dir] ), &tmat );
      su3_unit_der_spec( &( internal_V_link[4*i+dir] ), 
			 &( internal_Y_link[4*i+dir] ), &tmat,
                         &dwdv, &dwdagdv, &dwdvs, &dwdagdvs );
	  break;
  case UNITARIZE_U3:
      su3t4_copy( &dwdv, &dwdvs );
      su3t4_copy( &dwdagdv, &dwdagdvs );
	  break;
	default:
	  node0_printf("Unknown unitarization group. Set UNITARIZE_U3 or UNITARIZE_SU3\n");
	  terminate(1);
	}

      // adjoint piece of force from the previous level
      su3_adjoint( &( force_accum_old[dir][i] ), &tmat );

      // see LONG COMMENT in fermion_force_fn_multi_hisq
      clear_su3mat( &( force_accum[dir][i] ) );
      for( m=0; m<3; m++) {
        for( n=0; n<3; n++) {
          for( k=0; k<3; k++) {
            for( l=0; l<3; l++) {
              CMUL( dwdvs.t4[k][m][n][l], 
		    force_accum_old[dir][i].e[l][k], ftmp );
              CSUM( force_accum[dir][i].e[n][m], ftmp );
              /* CAREFUL with the adjoint part here! */
              CMUL( dwdagdvs.t4[k][m][n][l], tmat.e[l][k], ftmp );
              CSUM( force_accum[dir][i].e[n][m], ftmp );
            }
          }
        }
      }
      nflops_all_sites += 81*16;
      

#ifdef MILC_GLOBAL_DEBUG
#ifdef HISQ_FF_DEBUG
  {
    if( i==0 ) {
      printf("WRAPPER FORCE REUNIT old force\n");
      dumpmat( &( force_accum_old[dir][i] ) );
      printf("WRAPPER FORCE REUNIT new force\n");
      dumpmat( &( force_accum[dir][i] ) );
    }
  }
#endif /* HISQ_FF_DEBUG */
#endif /* MILC_GLOBAL_DEBUG */
    }
  }

    /* Sum over all nodes and update counters */

  // printf("nflops_all_sites = %lu\n", nflops_all_sites);

  info->final_flop += nflops_all_sites;

  /* phases in */
  rephase_field_offset( internal_V_link, ON, NULL , r0);
  rephase_field_offset( internal_Y_link, ON, NULL , r0);
} // fn_fermion_force_multi_hisq_reunit

#else /* ifndef USE_FF_GPU */

/********************************************************************/
/* GPU support                                                      */
/********************************************************************/

#include <quda_milc_interface.h>
#include "../include/generic_quda.h"

static void outer_product_append( Real one_hop_coeff, Real three_hop_coeff, 
				  Real *residues, 
				  su3_vector **multi_x, int nterms,
				  su3_matrix *one_hop_oprod[4], 
				  su3_matrix *three_hop_oprod[4] )
{
  /* note CG_solution and Dslash * solution are combined in "multi_x" */
  /* New version 1/21/99.  Use forward part of Dslash to get force */
  /* see long comment at end */
  /* For each link we need multi_x transported from both ends of path. */
  int term;
  int i,k;
  int dir;
  su3_matrix tmat;
  
  msg_tag *mtag[2];
  su3_matrix *mat_tmp0;
  // length N path has N+1 sites!!
  su3_matrix *oprod_along_path[MAX_PATH_LENGTH+1];
  
  
  if( nterms==0 )return;
  
#ifdef MILC_GLOBAL_DEBUG
#ifdef HISQ_FF_DEBUG
  {
    printf("WRAPPER FORCE smearing0 pseudofermion vector at site (0,0,0,0)\n");
    dumpvec( &(multi_x[0][0]) );
    printf("\n");
    su3_vector su3ab;
    clearvec( &su3ab );
    int nnn;
    for( nnn=0; nnn<nterms;nnn++) {
      FORALLFIELDSITES(i){
	add_su3_vector( &su3ab, &(multi_x[nnn][i]), &su3ab );
      }
    }
    printf("WRAPPER FORCE smearing0 pseudofermion vector GLOBAL SUM\n");
    dumpvec( &su3ab );
    printf("\n");
    for( nnn=0; nnn<nterms;nnn++) {
      printf("residues[%d]=%lf\n",nnn,residues[nnn]);
    }
  }
#endif /* HISQ_FF_DEBUG */
#endif /* MILC_GLOBAL_DEBUG */
  
  
  for(i=0;i<=MAX_PATH_LENGTH;i++){
    oprod_along_path[i] = 
      (su3_matrix *) malloc(sites_on_node*sizeof(su3_matrix) );
  }
  mat_tmp0 = (su3_matrix *) special_alloc(sites_on_node*sizeof(su3_matrix) );
  if( mat_tmp0 == NULL ){printf("Node %d NO ROOM\n",this_node); exit(0); }
  
  
  
  
  // clear force accumulators
  //  for(dir=XUP;dir<=TUP;dir++)
  //    FORALLFIELDSITES(i)clear_su3mat( &(one_hop_coeff[dir][i]) );
  
  //AB loop on directions, path table is not needed
  for(dir=XUP;dir<=TUP;dir++){
    k=0; // which gather we are using
    //AB netbackdir is just the opposite of dir for 1x1 case
    mtag[k] = start_gather_field( multi_x[0], sizeof(su3_vector),
				  OPP_DIR(dir), EVENANDODD, gen_pt[k] );
    FORALLFIELDSITES(i){
      clear_su3mat(  &(oprod_along_path[0][i]) ); // actually last site in path
    }
    for(term=0;term<nterms;term++){
      if(term<nterms-1)mtag[1-k] = start_gather_field( multi_x[term+1],
			       sizeof(su3_vector), OPP_DIR(dir), EVENANDODD, gen_pt[1-k] );
      wait_gather(mtag[k]);
      
      FORALLFIELDSITES(i){
        // build projector as usual
        su3_projector( &multi_x[term][i], (su3_vector *)gen_pt[k][i], &tmat );
        // multiply by alpha_l in rational function expansion
        scalar_mult_add_su3_matrix( &(oprod_along_path[0][i]), &tmat, 
				    residues[term], &(oprod_along_path[0][i]) );
      }
      cleanup_gather(mtag[k]);
      k=1-k; // swap 0 and 1
    } /* end loop over terms in rational function expansion */
    
    
    link_gather_connection_hisq( oprod_along_path[0], oprod_along_path[1], 
				 mat_tmp0, dir );
    
    FORALLFIELDSITES(i){
      scalar_mult_add_su3_matrix( &(one_hop_oprod[dir][i]), 
				  &(oprod_along_path[1][i]),
				  one_hop_coeff, &(one_hop_oprod[dir][i]) );
    }
    
  } /* end of loop on directions */
  
  
  /* *** Naik part *** */
  // clear force accumulators
  //  for(dir=XUP;dir<=TUP;dir++)
  //    FORALLFIELDSITES(i)clear_su3mat( &(three_hop_coeff[dir][i]) );
  
  for(dir=XUP;dir<=TUP;dir++){ //AB loop on directions, path table is not needed
    k=0; // which gather we are using
    mtag[k] = start_gather_field( multi_x[0], sizeof(su3_vector),
				  OPP_3_DIR( DIR3(dir) ), EVENANDODD, gen_pt[k] );
    FORALLFIELDSITES(i){
      clear_su3mat(  &(oprod_along_path[0][i]) );
    }
    for(term=0;term<nterms;term++){
      if(term<nterms-1)mtag[1-k] = start_gather_field( multi_x[term+1],
	       sizeof(su3_vector), OPP_3_DIR( DIR3(dir) ), EVENANDODD, gen_pt[1-k] );
      wait_gather(mtag[k]);
      
      FORALLFIELDSITES(i){
        // build projector as usual
        su3_projector( &multi_x[term][i], (su3_vector *)gen_pt[k][i], &tmat );
        // multiply by alpha_l in rational function expansion
        scalar_mult_add_su3_matrix( &(oprod_along_path[0][i]), &tmat, 
				    residues[term], &(oprod_along_path[0][i]) );
      }
      cleanup_gather(mtag[k]);
      k=1-k; // swap 0 and 1
    } /* end loop over terms in rational function expansion */
    
    link_gather_connection_hisq( oprod_along_path[0], oprod_along_path[1], 
				 mat_tmp0, DIR3(dir) );
    
    FORALLFIELDSITES(i){
      scalar_mult_add_su3_matrix( &(three_hop_oprod[dir][i]), 
				  &(oprod_along_path[1][i]),
				  three_hop_coeff, &(three_hop_oprod[dir][i]) );
    }
  } /* end of loop on directions */
  
  
  free( mat_tmp0 );
  for(i=0;i<=MAX_PATH_LENGTH;i++){
    free( oprod_along_path[i] );
  }
  return;
} // outer_product_append


static void outer_product_create( Real one_hop_coeff, Real three_hop_coeff, 
	  Real *residues, 
	  su3_vector **multi_x, int nterms,
	  su3_matrix *one_hop_oprod[4], su3_matrix *three_hop_oprod[4] )
{
  int dir, i;
  for(dir=XUP;dir<=TUP;dir++){
    FORALLFIELDSITES(i){ 
      clear_su3mat( &(one_hop_oprod[dir][i]) );
      clear_su3mat( &(three_hop_oprod[dir][i]));
    }
  }
  
  outer_product_append(one_hop_coeff, three_hop_coeff,
		       residues, 
		       multi_x, 
		       nterms,
		       one_hop_oprod,
		       three_hop_oprod);
  return;
} // outer_product_create




static void write_path_coeffs_to_array(ks_component_paths paths, double array[6])
{
  array[0] = paths.act_path_coeff.one_link; 
  array[1] = paths.act_path_coeff.naik;
  array[2] = paths.act_path_coeff.three_staple;
  array[3] = paths.act_path_coeff.five_staple;
  array[4] = paths.act_path_coeff.seven_staple;
  array[5] = paths.act_path_coeff.lepage;
  return;
}


static void 
fn_fermion_force_multi_hisq_wrapper_mx_gpu(info_t* info, Real eps, Real *residues, 
					   su3_vector **multi_x, int nterms,
					   fermion_links_t *fl)
{
  // Get the terms we need from the fermion links structure */
  ks_action_paths_hisq *ap = get_action_paths_hisq(fl);
  hisq_auxiliary_t *aux = get_hisq_auxiliary(fl);
  int n_naiks = fermion_links_get_n_naiks(fl);
  double *eps_naik = fermion_links_get_eps_naik(fl);
  su3_matrix *U_link = aux->U_link;
  su3_matrix *V_link = aux->V_link;
  su3_matrix *W_unitlink = aux->W_unitlink;
  site *s;
  int i, j, dir;
  int inaik;
  int n_naik_shift;
  int n_orders_naik_current;
  
  su3_matrix* one_link_oprod[4];
  su3_matrix* three_link_oprod[4];
  su3_matrix* staple_oprod[4];

  initialize_quda();
  
  // First construct the outer products
  for(i=XUP; i<=TUP; ++i){
    one_link_oprod[i] = (su3_matrix*)malloc(sites_on_node*sizeof(su3_matrix));
    three_link_oprod[i] = (su3_matrix*)malloc(sites_on_node*sizeof(su3_matrix));
    staple_oprod[i] = (su3_matrix*)malloc(sites_on_node*sizeof(su3_matrix));
  }
  
  n_naik_shift = 0;
  n_orders_naik_current = n_order_naik_total;
  Real one_hop_coeff = 2.0*eps;
  Real three_hop_coeff = ap->p2.act_path_coeff.naik*2.0*eps;   
  
  outer_product_create(one_hop_coeff, three_hop_coeff, residues, 
		       multi_x, n_orders_naik_current, staple_oprod, three_link_oprod);
  
  for(dir=XUP; dir<=TUP; ++dir){
    FORALLFIELDSITES(i){
      scalar_mult_su3_matrix(&(staple_oprod[dir][i]), ap->p2.act_path_coeff.one_link, 
			     &(one_link_oprod[dir][i]));
    }
  }
  
  n_naik_shift = n_orders_naik[0];
  for(inaik=1; inaik<n_naiks; ++inaik){
    n_orders_naik_current = n_orders_naik[inaik];
    
    one_hop_coeff = ap->p3.act_path_coeff.one_link*eps_naik[inaik]*2.0*eps;
    three_hop_coeff = ap->p3.act_path_coeff.naik*eps_naik[inaik]*2.0*eps;
    
    outer_product_append(one_hop_coeff, 
			 three_hop_coeff, 
			 residues+n_naik_shift,
			 multi_x+n_naik_shift, n_orders_naik_current, 
			 one_link_oprod,
			 three_link_oprod);
    
    n_naik_shift += n_orders_naik[inaik];  
  } // end loop over inaik
  // done constructing the outer products
  
  Real* momentum = (Real*)special_alloc(sites_on_node*4*sizeof(anti_hermitmat));
  
  
  double level2_coeff[6];
  write_path_coeffs_to_array(ap->p2, level2_coeff);
  double fat7_coeff[6];
  write_path_coeffs_to_array(ap->p1, fat7_coeff);
  
  // The following lines need to be deleted
  const int dim[4] = {nx,ny,nz,nt};
  QudaLayout_t layout;
  layout.latsize = dim;
  qudaInit(layout);
  // end delete
  
  // The following lines of code are optional
  // If the parameters are not set explicitly here, 
  // they are set to default values the first time 
  // qudaHisqForce is called.
  QudaHisqParams_t params;
  params.reunit_allow_svd = 1;
  params.reunit_svd_only = 0;
  params.reunit_svd_abs_error = 1e-8;
  params.reunit_svd_rel_error = 1e-7;
  params.force_filter = 5e-5;

  qudaHisqParamsInit(params);
  // end optional code
		
  qudaHisqForce(PRECISION, level2_coeff, fat7_coeff, 
		(const void* const*)staple_oprod, 
		(const void* const*)one_link_oprod, 
		(const void* const*)three_link_oprod, 
		W_unitlink, V_link, U_link, momentum);
  
  // append result
  FORALLSITES(i,s){
    for(dir=0; dir<4; ++dir){
      for(j=0; j<10; ++j){
	*((Real*)(&(s->mom[dir])) + j) += *(momentum + (4*i+ dir)*10 + j);
      }
    }
  }
  
  for(i=XUP; i<=TUP; ++i){
    free(one_link_oprod[i]);
    free(three_link_oprod[i]);
    free(staple_oprod[i]);
  }
  free(momentum);
  
  return;
} //fn_fermion_force_multi_hisq_wrapper_mx_gpu

#endif /* ifndef USE_FF_GPU */

#endif  /* HISQ_FF_MULTI_WRAPPER || USE_FF_GPU */

int 
find_backwards_gather( Q_path *path )
{
    int disp[4], i;
    /* compute total displacement of path */
    for(i=XUP;i<=TUP;i++)disp[i]=0;
    for( i=0; i<path->length; i++){
	if( GOES_FORWARDS(path->dir[i]) )
	    disp[        path->dir[i]  ]++;
	else
	    disp[OPP_DIR(path->dir[i]) ]--;
    }

   // There must be an elegant way??
   if( disp[XUP]==+1 && disp[YUP]== 0 && disp[ZUP]== 0 && disp[TUP]== 0 )
     return(XDOWN);
   if( disp[XUP]==-1 && disp[YUP]== 0 && disp[ZUP]== 0 && disp[TUP]== 0 )
     return(XUP);
   if( disp[XUP]== 0 && disp[YUP]==+1 && disp[ZUP]== 0 && disp[TUP]== 0 )
     return(YDOWN);
   if( disp[XUP]== 0 && disp[YUP]==-1 && disp[ZUP]== 0 && disp[TUP]== 0 )
     return(YUP);
   if( disp[XUP]== 0 && disp[YUP]== 0 && disp[ZUP]==+1 && disp[TUP]== 0 )
     return(ZDOWN);
   if( disp[XUP]== 0 && disp[YUP]== 0 && disp[ZUP]==-1 && disp[TUP]== 0 )
     return(ZUP);
   if( disp[XUP]== 0 && disp[YUP]== 0 && disp[ZUP]== 0 && disp[TUP]==+1 )
     return(TDOWN);
   if( disp[XUP]== 0 && disp[YUP]== 0 && disp[ZUP]== 0 && disp[TUP]==-1 )
     return(TUP);

   if( disp[XUP]==+3 && disp[YUP]== 0 && disp[ZUP]== 0 && disp[TUP]== 0 )
     return(X3DOWN);
   if( disp[XUP]==-3 && disp[YUP]== 0 && disp[ZUP]== 0 && disp[TUP]== 0 )
     return(X3UP);
   if( disp[XUP]== 0 && disp[YUP]==+3 && disp[ZUP]== 0 && disp[TUP]== 0 )
     return(Y3DOWN);
   if( disp[XUP]== 0 && disp[YUP]==-3 && disp[ZUP]== 0 && disp[TUP]== 0 )
     return(Y3UP);
   if( disp[XUP]== 0 && disp[YUP]== 0 && disp[ZUP]==+3 && disp[TUP]== 0 )
     return(Z3DOWN);
   if( disp[XUP]== 0 && disp[YUP]== 0 && disp[ZUP]==-3 && disp[TUP]== 0 )
     return(Z3UP);
   if( disp[XUP]== 0 && disp[YUP]== 0 && disp[ZUP]== 0 && disp[TUP]==+3 )
     return(T3DOWN);
   if( disp[XUP]== 0 && disp[YUP]== 0 && disp[ZUP]== 0 && disp[TUP]==-3 )
     return(T3UP);
   node0_printf("OOOPS: NODIR\n"); exit(0);
   return( NODIR );
} //find_backwards_gather


// Make a new path table.  Sorted principally by total displacement of path.
// Below that, sort by direction of first link
// Below that, sort by direction of second link - note special case of one link paths
int 
sort_quark_paths( Q_path *src_table, Q_path *dest_table, int npaths, 
		  int num_back_dirs )
{
    int netdir,dir0,dir1,dir1tmp,num_new,i,j;

    num_new=0; // number of paths in sorted table
    for( i=0; i<num_back_dirs; i++ ){ // loop over net_back_dirs
        netdir = net_back_dirs[i]; // table of possible displacements for Fat-Naik
	for( dir0=0; dir0<=7; dir0++){ // XUP ... TDOWN
	  for( dir1=-1; dir1<=7; dir1++){ // NODIR, XUP ... TDOWN
	    if( dir1==-1 )dir1tmp=NODIR; else dir1tmp=dir1;
	    for( j=0; j<npaths; j++ ){ // pick out paths with right net displacement
	        if( find_backwards_gather( &(src_table[j]) ) == netdir && 
			src_table[j].dir[0]==dir0 &&
			src_table[j].dir[1]==dir1tmp ){
		    dest_table[num_new] = src_table[j];
		    num_new++;
	        }
	    } // loop over paths
	  } //dir1
	} //dir0
    }
    if( num_new!=npaths){ node0_printf("OOPS: path table error, num_new=%d, npaths=%d\n",num_new,npaths); exit(0); }
    return 0;
}



/**********************************************************************/
/*   Version for HISQ.  Parallel transport nterms source vectors    */
/**********************************************************************/

void u_shift_veclist_fermion(veclist *src, veclist *dest, 
	     int dir, msg_tag **mtag, veclist *tmpvec, int listlength );
void add_3f_force_to_mom_list(veclist *back,
	 veclist *forw, int dir, Real coeff[2], int listlength ) ;
void side_link_3f_force_list(int mu, int nu, Real *coeff, 
	veclist *Path, veclist *Path_nu, 
	veclist *Path_mu, veclist *Path_numu, int listlength ) ;

static su3_matrix *backwardlink[4];
static anti_hermitmat *tempmom[4];

void u_shift_veclist_fermion(veclist *src, veclist *dest, 
	int dir, msg_tag **mtag, veclist *tmpvec, int listlength ) {
#if 0
  site *s ;
  int i ;
#endif

#ifdef FFSTIME
  double time0, time1;
  time1 = -dclock();
#endif
  /* Copy source to temporary vector */
  memcpy( (char *)tmpvec, (char *)src, 
	 sites_on_node*sizeof(veclist) );

  if(*mtag == NULL)
    *mtag = start_gather_field(tmpvec, sizeof(veclist), 
			       dir, EVENANDODD, gen_pt[dir]);
  else
    restart_gather_field(tmpvec, sizeof(veclist), 
			 dir, EVENANDODD, gen_pt[dir], *mtag);
  wait_gather(*mtag);

#ifdef FFSTIME
  time0 = -dclock();
  time1 -= time0;
#endif

  if(GOES_FORWARDS(dir)) /* forward shift */
    {
#if 0
      FORALLSITES(i,s)
	mult_su3_mat_hwvec( &(s->link[dir]), 
			    (half_wilson_vector *)gen_pt[dir][i], 
			    dest + i );
#else
      mult_su3_sitelink_latveclist(dir, (veclist **)gen_pt[dir], dest, listlength );
#endif
    }
  else /* backward shift */
    {
#if 0
      FORALLFIELDSITES(i)
	mult_adj_su3_mat_hwvec( backwardlink[OPP_DIR(dir)] + i, 
				(half_wilson_vector *)gen_pt[dir][i], 
				dest + i );
#else
     // mult_adj_su3_fieldlink_lathwvec( backwardlink[OPP_DIR(dir)],
				       //(half_wilson_vector **)gen_pt[dir],
				       //dest);
      mult_adj_su3_fieldlink_latveclist( backwardlink[OPP_DIR(dir)],
	(veclist **)gen_pt[dir], dest, listlength );
#endif
    }

#ifdef FFSTIME
  time0 += dclock();
  node0_printf("FFSHIFT time0 = %e\nFFSHIFT time1 = %e\n",time0,time1);
#endif
}


void scalar_mult_add_lathwvec_proj(anti_hermitmat *mom, half_wilson_vector *back, 
				   half_wilson_vector *forw, Real coeff[2]);


void add_3f_force_to_mom_list(veclist *back,
	veclist *forw, int dir, Real *coeff, int listlength ) {
#if 0
  register site *s ;
  register int i;
  su3_matrix tmat, *tmat2;
  Real *tmp_coeff = (Real *)malloc(listlength*sizeof(Real) );
#endif
  int j ;  
  Real *my_coeff = (Real *)malloc(listlength*sizeof(Real) );
  int mydir;
#ifdef FFSTIME
  double time;
  time = -dclock();
#endif

  if(GOES_BACKWARDS(dir)) {
    mydir = OPP_DIR(dir); for( j=0; j<listlength; j++ ){my_coeff[j] = -coeff[j];}
  }
  else {
    mydir = dir; for( j=0; j<listlength; j++ ){my_coeff[j] = coeff[j];}
  }
#if 0
  FORALLFIELDSITES(i){
    if(s->parity==ODD)
      {	tmp_coeff[0] = -my_coeff[0] ; tmp_coeff[1] = -my_coeff[1] ; }
    else
      { tmp_coeff[0] = my_coeff[0] ; tmp_coeff[1] = my_coeff[1] ; }
    tmat2 = tempmom[mydir] + i;
    su3_projector(&(back[i].h[0]), &(forw[i].h[0]), &tmat);
    scalar_mult_add_su3_matrix(tmat2, &tmat,  tmp_coeff[0], tmat2 );
    su3_projector(&(back[i].h[1]), &(forw[i].h[1]), &tmat);
    scalar_mult_add_su3_matrix(tmat2, &tmat,  tmp_coeff[1], tmat2 );
  }
#else
  scalar_mult_add_latveclist_proj( tempmom[mydir], back,
	forw, my_coeff, listlength );
#endif

#ifdef FFSTIME
  time += dclock();
  node0_printf("FFSHIFT time2 = %e\n", time);
#endif
  free(my_coeff); 
#if 0
  free(tmp_coeff);
#endif
}

 

void side_link_3f_force_list(int mu, int nu, Real *coeff, 
	veclist *Path, veclist *Path_nu, veclist *Path_mu, veclist *Path_numu,
	int listlength ) {

  Real *m_coeff;
  int j;

  m_coeff = (Real *)malloc( listlength*sizeof(Real) );
  for( j=0; j<listlength; j++){
     m_coeff[j] = -coeff[j] ;
  }

  if(GOES_FORWARDS(mu))
    {
      /*                    nu           * 
       * Add the force :  +----+         *
       *               mu |    |         *
       *                  x    (x)       *
       *                  o    o         */
      if(GOES_FORWARDS(nu))
	add_3f_force_to_mom_list( Path_numu, Path, mu, coeff, listlength ) ;
      else
	add_3f_force_to_mom_list( Path, Path_numu,OPP_DIR(mu),m_coeff, listlength);/* ? extra - */
    }
  else /*GOES_BACKWARDS(mu)*/
    {
      /* Add the force :  o    o         *
       *               mu |    |         *
       *                  x    (x)       *
       *                  +----+         *
       *                    nu           */ 
      if(GOES_FORWARDS(nu))
	add_3f_force_to_mom_list( Path_nu, Path_mu, mu, m_coeff, listlength) ; /* ? extra - */
      else
	add_3f_force_to_mom_list( Path_mu, Path_nu, OPP_DIR(mu), coeff, listlength) ;
    }
    free(m_coeff);
}


/* LONG COMMENTS
   Here we have combined "xxx", (offset "x_off")  which is
(M_adjoint M)^{-1} phi, with Dslash times this vector, which goes in the
odd sites of xxx.  Recall that phi is defined only on even sites.  In
computing the fermion force, we are looking at

< X |  d/dt ( Dslash_eo Dslash_oe ) | X >
=
< X | d/dt Dslash_eo | T > + < T | d/dt Dslash_oe | X >
where T = Dslash X.

The subsequent manipulations to get the coefficent of H, the momentum
matrix, in the simulation time derivative above look the same for
the two terms, except for a minus sign at the end, if we simply stick
T, which lives on odd sites, into the odd sites of X

 Each path in the action contributes terms when any link of the path
is the link for which we are computing the force.  We get a minus sign
for odd numbered links in the path, since they connect sites of the
opposite parity from what it would be for an even numbered link.
Minus signs from "going around" plaquette - ie KS phases, are supposed
to be already encoded in the path coefficients.
Minus signs from paths that go backwards are supposed to be already
encoded in the path coefficients.

Here, for example, are comments reproduced from the force routine for
the one-link plus Naik plus single-staple-fat-link action:

 The three link force has three contributions, where the link that
was differentiated is the first, second, or third link in the 3-link
path, respectively.  Diagramatically, where "O" represents the momentum,
the solid line the link corresponding to the momentum, and the dashed
lines the other links:
 

	O______________ x ............ x ...............
+
	x..............O______________x.................
+
	x..............x..............O________________
Think of this as
	< xxx | O | UUUxxx >		(  xxx, UUUX_p3 )
+
	< xxx U | O | UUxxx >		( X_m1U , UUX_p2 )
+
	< xxx U U | O | Uxxx >		( X_m2UU , UX_p1 )
where "U" indicates parallel transport, "X_p3" is xxx displaced
by +3, etc.
Note the second contribution has a relative minus sign
because it effectively contributes to the <odd|even>, or M_adjoint,
part of the force when we work on an even site. i.e., for M on
an even site, this three link path begins on an odd site.

The staple force has six contributions from each plane containing the
link direction:
Call these diagrams A-F:


	x...........x		O____________x
		    .			     .
		    .			     .
		    .			     .
		    .			     .
		    .			     .
	O___________x		x............x
	   (A)			    (B)



	x	    x		O____________x
	.	    .		.	     .
	.	    .		.	     .
	.	    .		.	     .
	.	    .		.	     .
	.	    .		.	     .
	O___________x		x	     x
	   (C)			    (D)



	x...........x		O____________x
	.			.
	.			.
	.			.
	.			.
	.			.
	O___________x		x............x
	   (E)			    (F)

As with the Naik term, diagrams C and D have a relative minus
sign because they connect sites of the other parity.

Also note an overall minus sign in the staple terms relative to the
one link term because, with the KS phase factors included, the fat
link is  "U - w3 * UUU", or the straight link MINUS w3 times the staples.

Finally, diagrams B and E get one more minus sign because the link
we are differentiating is in the opposite direction from the staple
as a whole.  You can think of this as this "U" being a correction to
a "U_adjoint", but the derivative of U is iHU and the derivative
of U_adjoint is -iHU_adjoint.

*/
/* LONG COMMENT on sign conventions
In most of the program, the KS phases and antiperiodic boundary
conditions are absorbed into the link matrices.  This greatly simplfies
multiplying by the fermion matrix.  However, it requires care in
specifying the path coefficients.  Remember that each time you
encircle a plaquette, you pick up a net minus sign from the KS phases.
Thus, when you have more than one path to the same point, you generally
have a relative minus sign for each plaquette in a surface bounded by
this path and the basic path for that displacement.

Examples:
  Fat Link:
    Positive:	X-------X

    Negative     --------
	 	|	|
		|	|
		X	X

  Naik connection, smeared
    Positive:	X-------x-------x-------X

    Negative:	---------
		|	|
		|	|
		X	x-------x-------X

    Positive:	--------x--------
		|		|
		|		|
		X		x-------X

    Negative:	--------x-------x-------x
		|			|
		|			|
		X			X
*/



/* Comment on acceptable actions.
   We construct the backwards part of dslash by reversing all the
   paths in the forwards part.  So, for example, in the p4 action
   the forwards part includes +X+Y+Y

		X
		|
		|
		X
		|
		|
	X---->--X

  so we put -X-Y-Y in the backwards part.  But this isn't the adjoint
  of U_x(0)U_y(+x)U_y(+x+y).  Since much of the code assumes that the
  backwards hop is the adjoint of the forwards (for example, in
  preventing going to 8 flavors), the code only works for actions
  where this is true.  Roughly, this means that the fat link must
  be symmetric about reflection around its midpoint.  Equivalently,
  the paths in the backwards part of Dslash are translations of the
  paths in the forwards part.  In the case of the "P4" or knight's move
  action, this means that we have to have both paths
   +X+Y+Y and +Y+Y+X to the same point, with the same coefficients.
  Alternatively, we could just use the symmetric path +Y+X+Y.
*/
