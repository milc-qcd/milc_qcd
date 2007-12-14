/****** fermion_force_fn_multi.c  -- ******************/
/* MIMD version 7 */
/* Multisource fermion force.  Includes various optimization choices.
 * 
 * 1. General force for any FN-type action. Based on fermion_force_general.c
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
 * A.B. 10/07 unitarization and its derivatives
 *
 */


#include "generic_ks_includes.h"	/* definitions files and prototypes */
#include "../include/umethod.h"
#include <string.h>
void printpath( int *path, int length );


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


// dump matrix on this site/link
#define AB_OUT_ON_SITE 3
#define AB_OUT_ON_LINK TUP

//AB block 071015
// fermion_force_fn_multi() does calculate conjugate (as in previous version)
#define AB_CHANGED_ASQTAD
// THESE ROUTINES DO NOT HANDLE NAIK LINKS!!!
// IT CAN BE FIXED BY INTRODUCING EXTRA univ0 TYPE ROUTINE

// Forward declarations
// Universal wrapper routine
static void 
fn_fermion_force_multi_wrapper( Real eps, Real *residues, 
				su3_vector **multi_x, int nterms,
				ferm_links_t *fn, ks_action_paths *ap );
// Contribution to the force from 0 level of smearing, force contains only |X><Y| terms
static void 
fn_fermion_force_multi_univ0( Real eps, Real *residues, 
			      su3_vector **multi_x, int nterms,
			      su3_matrix *force_accum[4] );
// This routine calculates contribution to the force due to i-th level of smearing,
// as input it uses array of new links and the forth from (i-1)-th level of smearing
static void 
fn_fermion_force_multi_universal( Real eps, Real *residues, su3_vector **multi_x, 
				  int nterms, su3_matrix *force_accum[4], 
				  su3_matrix *force_accum_old[4], 
				  su3_matrix **internal_U_link,
				  int internal_num_q_paths, 
				  Q_path *internal_q_paths_sorted, 
				  int *internal_netbackdir_table );
// Contribution to force due to reunitarization of links: W=V*(V^+*V)^-1/2
static void 
fn_fermion_force_multi_reunit( su3_matrix *force_accum[4], 
			       su3_matrix *force_accum_old[4], 
			       ferm_links_t *fn, ks_action_paths *ap );
//End of AB block 071015

// More forward declarations
static void 
fn_fermion_force_multi_hisq( Real eps, Real *residues, su3_vector **multi_x, 
			     int nterms, ferm_links_t *fn, 
			     ks_action_paths *ap );


/**********************************************************************/
/*   Wrapper for fermion force routines with multiple sources         */
/**********************************************************************/
void eo_fermion_force_multi( Real eps, Real *residues, su3_vector **xxx, 
			     int nterms, int prec, ferm_links_t *fn,
			     ks_action_paths *ap ) {
  switch(KS_MULTIFF){
  case FNMATREV:
    node0_printf("FNMATREV isn't implemented\n"); exit(0);
    break;
  case FNMAT:
    fn_fermion_force_multi_hisq( eps, residues, xxx, nterms, fn, ap ); 
    //fn_fermion_force_multi_wrapper( eps, residues, xxx, nterms ); //TEMP, DOES NOT HANDLE NAIK!!!
    break;
  default:
    fn_fermion_force_multi_hisq( eps, residues, xxx, nterms, fn, ap ); //TEMP
  }
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
Q_path *q_paths_sorted_1 = NULL;	// Quark paths sorted by net displacement and last directions
int *netbackdir_table_1 = NULL; // table of net path displacements (backwards from usual convention)
// table of gather directions to bring start of path to end, for "FN" = "fat-Naik" actions
int net_back_dirs_1[8] = { XDOWN, YDOWN, ZDOWN, TDOWN, XUP, YUP, ZUP, TUP };
// For smearing 2
Q_path *q_paths_sorted_2 = NULL;	// Quark paths sorted by net displacement and last directions
int *netbackdir_table_2 = NULL; // table of net path displacements (backwards from usual convention)
// table of gather directions to bring start of path to end, for "FN" = "fat-Naik" actions
int net_back_dirs_2[16] = { XDOWN, YDOWN, ZDOWN, TDOWN, XUP, YUP, ZUP, TUP, 
	X3DOWN, Y3DOWN, Z3DOWN, T3DOWN, X3UP, Y3UP, Z3UP, T3UP };

// For old code
Q_path *q_paths_sorted = NULL;	// Quark paths sorted by net displacement and last directions
int *netbackdir_table = NULL; // table of net path displacements (backwards from usual convention)
// table of gather directions to bring start of path to end, for "FN" = "fat-Naik" actions
int net_back_dirs[16] = { XDOWN, YDOWN, ZDOWN, TDOWN, XUP, YUP, ZUP, TUP, 
	X3DOWN, Y3DOWN, Z3DOWN, T3DOWN, X3UP, Y3UP, Z3UP, T3UP };

int sort_quark_paths( Q_path *src_table, Q_path *dest_table, int npaths, int num_back_dirs );
int find_backwards_gather( Q_path *path );
int first_force=1;	// 1 if force hasn't been called yet




//TEMP, AB070910
//Attempt: a universal force routine to calculate the contribution
//to the force for any level of smearing iteratively
// univ0 -- calculate |X><Y| without F, B -- no smearing (or smearing level 0) ??!
// see fn_fermion_force_multi_universal below!
static void 
fn_fermion_force_multi_univ0( Real eps, Real *residues, 
			      su3_vector **multi_x, int nterms,
			      su3_matrix *force_accum[4] )
{
  /* note CG_solution and Dslash * solution are combined in "multi_x" */
  /* New version 1/21/99.  Use forward part of Dslash to get force */
  /* see long comment at end */
  /* For each link we need multi_x transported from both ends of path. */
  int term;
  register int i,k;
  register site *s;
  int dir;
  su3_matrix tmat;
  Real ferm_epsilon, coeff;

  msg_tag *mtag[2];
  su3_matrix *mat_tmp0;
  su3_matrix *oprod_along_path[MAX_PATH_LENGTH+1]; // length N path has N+1 sites!!


  /* node0_printf("STARTING fn_fermion_force_multi() nterms = %d\n",nterms); */
  if( nterms==0 )return;


  for(i=0;i<=MAX_PATH_LENGTH;i++){
     oprod_along_path[i] = (su3_matrix *) malloc(sites_on_node*sizeof(su3_matrix) );
  }
  mat_tmp0 = (su3_matrix *) special_alloc(sites_on_node*sizeof(su3_matrix) );
  if( mat_tmp0 == NULL ){printf("Node %d NO ROOM\n",this_node); exit(0); }


  ferm_epsilon = 2.0*eps; // we only do forward paths, gives factor of 2


  // clear force accumulators
  for(dir=XUP;dir<=TUP;dir++)FORALLSITES(i,s)clear_su3mat( &(force_accum[dir][i]) );

  for(dir=XUP;dir<=TUP;dir++){ //AB loop on directions, path table is not needed
    k=0; // which gather we are using
//AB        mtag[k] = start_gather_field( multi_x[0], sizeof(su3_vector),
//AB           netbackdir, EVENANDODD, gen_pt[k] );
    mtag[k] = start_gather_field( multi_x[0], sizeof(su3_vector),
      OPP_DIR(dir), EVENANDODD, gen_pt[k] ); //AB netbackdir is just the opposite of dir for 1x1 case
    FORALLSITES(i,s){
      clear_su3mat(  &(oprod_along_path[0][i]) ); // actually last site in path
    }
    for(term=0;term<nterms;term++){
      if(term<nterms-1)mtag[1-k] = start_gather_field( multi_x[term+1],
      sizeof(su3_vector), OPP_DIR(dir), EVENANDODD, gen_pt[1-k] );
      wait_gather(mtag[k]);

      FORALLSITES(i,s){
        //AB build projector as usual
        su3_projector( &multi_x[term][i], (su3_vector *)gen_pt[k][i], &tmat );
        //AB multiply by alpha_l in rational function expansion
        scalar_mult_add_su3_matrix( &(oprod_along_path[0][i]), &tmat, residues[term], &(oprod_along_path[0][i]) );
      }
      cleanup_gather(mtag[k]);
      k=1-k; // swap 0 and 1
    } /* end loop over terms in rational function expansion */


    link_gather_connection_hisq( oprod_along_path[0], oprod_along_path[1], mat_tmp0, dir );

    coeff = 1; //AB fermion_eps is outside this routine in "wrapper" routine
/*
    FOREVENSITES(i,s){
      scalar_mult_add_su3_matrix( &(force_accum[dir][i]), &(oprod_along_path[1][i]),
                coeff, &(force_accum[dir][i]) );
    }
    FORODDSITES(i,s){
      scalar_mult_add_su3_matrix( &(force_accum[dir][i]), &(oprod_along_path[1][i]),
                -coeff, &(force_accum[dir][i]) );
    }
*/
    FORALLSITES(i,s){
      scalar_mult_add_su3_matrix( &(force_accum[dir][i]), &(oprod_along_path[1][i]),
                coeff, &(force_accum[dir][i]) );
    }

  } /* end of loop on directions */

	
  free( mat_tmp0 );
  for(i=0;i<=MAX_PATH_LENGTH;i++){
     free( oprod_along_path[i] );
  }
}





//AB: UNIVERSAL_FORCE
static void 
fn_fermion_force_multi_universal( Real eps, Real *residues, su3_vector **multi_x, 
				  int nterms, su3_matrix *force_accum[4], 
				  su3_matrix *force_accum_old[4], 
				  su3_matrix **internal_U_link,
				  int internal_num_q_paths, 
				  Q_path *internal_q_paths_sorted, 
				  int *internal_netbackdir_table )
{
  /* note CG_solution and Dslash * solution are combined in "multi_x" */
  /* New version 1/21/99.  Use forward part of Dslash to get force */
  /* see long comment at end */
  /* For each link we need multi_x transported from both ends of path. */
  register int i,j,k,lastdir=-99,ipath,ilink;
  register site *s;
  int length,dir,odir;
  Real ferm_epsilon, coeff;
  Q_path *this_path;	// pointer to current path
  Q_path *last_path;	// pointer to previous path
  msg_tag *mtag[2];
  su3_matrix *mat_tmp0;
  su3_matrix *mat_tmp1;
  anti_hermitmat *ahmat_tmp;
  su3_matrix *oprod_along_path[MAX_PATH_LENGTH+1]; // length N path has N+1 sites!!
  su3_matrix *mats_along_path[MAX_PATH_LENGTH+1]; // 
//  su3_matrix *force_accum[4];  // accumulate force
//  su3_matrix *force_final[4];  // accumulate force
  int netbackdir, last_netbackdir;	// backwards direction for entire path
//int tempflops = 0; //TEMP

  /* node0_printf("STARTING fn_fermion_force_multi() nterms = %d\n",nterms); */
  if( nterms==0 )return;

  for(i=0;i<=MAX_PATH_LENGTH;i++){
     oprod_along_path[i] = (su3_matrix *) malloc(sites_on_node*sizeof(su3_matrix) );
  }
  for(i=1;i<=MAX_PATH_LENGTH;i++){ // 0 element is never used (it's unit matrix)
     mats_along_path[i] = (su3_matrix *) malloc(sites_on_node*sizeof(su3_matrix) );
  }

  mat_tmp0 = (su3_matrix *) special_alloc(sites_on_node*sizeof(su3_matrix) );
  mat_tmp1 = (su3_matrix *) special_alloc(sites_on_node*sizeof(su3_matrix) );
  if( mat_tmp0 == NULL ){printf("Node %d NO ROOM\n",this_node); exit(0); }
  if( mat_tmp1 == NULL ){printf("Node %d NO ROOM\n",this_node); exit(0); }

  ahmat_tmp = (anti_hermitmat *) special_alloc(sites_on_node*sizeof(anti_hermitmat) );

	
  ferm_epsilon = 2.0*eps; // we only do forward paths, gives factor of 2

  // clear force accumulators
  for(dir=XUP;dir<=TUP;dir++)FORALLSITES(i,s)clear_su3mat( &(force_accum[dir][i]) );

//#ifndef AB_OLD_UNIV
//for(idir=XUP;idir<=TUP;idir++){ /* loop over link for f(0)_\nu */
//#endif /* AB_OLD_UNIV */


  /* loop over paths, and loop over links in path */
  last_netbackdir = NODIR;
  last_path = NULL;
  for( ipath=0; ipath<internal_num_q_paths; ipath++ ){
    this_path = &(internal_q_paths_sorted[ipath]);
    if(this_path->forwback== -1)continue;	/* skip backwards dslash */

    length = this_path->length;
    // find gather to bring multi_x[term] from "this site" to end of path
    //netbackdir = find_backwards_gather( &(q_paths_sorted[ipath]) );
    netbackdir = internal_netbackdir_table[ipath];
//AB    printf("netbackdir[ipath=%d]=%d\n",ipath,netbackdir);

    // and bring multi_x to end - no gauge transformation !!
    // resulting outer product matrix has gauge transformation properties of a connection
    // from start to end of path


    /* move f(i-1) force from current site in positive direction,
       this corresponds to outer product |X><Y| calculated at the endpoint of the path */
    link_gather_connection_hisq( force_accum_old[OPP_DIR(netbackdir)], 
				 oprod_along_path[0], NULL, netbackdir );


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
//ABopt    if( netbackdir == last_netbackdir )
//ABopt      while ( j>0 && this_path->dir[j] == last_path->dir[j+last_path->length-length] ) j--;
    if( GOES_BACKWARDS(this_path->dir[0]) ) k=1; else k=0;


    for(ilink=j;ilink>=k;ilink--){
        link_transport_connection_hisq( oprod_along_path[length-ilink-1], internal_U_link,
        oprod_along_path[length-ilink], mat_tmp0, this_path->dir[ilink]  );
    }

   /* maintain an array of transports "to this point" along the path.
	Don't recompute beginning parts of path if same as last path */
    ilink=0; // first link where new transport is needed
//ABopt    if( last_path != NULL )while( this_path->dir[ilink] == last_path->dir[ilink] ) ilink++ ;
    // Sometimes we don't need the matrix for the last link
    if( GOES_FORWARDS(this_path->dir[length-1]) ) k=length-1; else k=length;

    for( ; ilink<k; ilink++ ){
      if( ilink==0 ){
        dir = this_path->dir[0];
          if( GOES_FORWARDS(dir) ){
            mtag[1] = start_gather_site( F_OFFSET(link[dir]), sizeof(su3_matrix),
                   OPP_DIR(dir), EVENANDODD, gen_pt[1] );
            wait_gather(mtag[1]);
            FORALLSITES(i,s){ su3_adjoint( (su3_matrix *)gen_pt[1][i], &(mats_along_path[1][i]) ); }
            cleanup_gather(mtag[1]);
          }
          else{
            FORALLSITES(i,s){ mats_along_path[1][i] = s->link[OPP_DIR(dir)]; }
          }
      }
      else { // ilink != 0
        dir = OPP_DIR(this_path->dir[ilink]);
        //TEST link_transport_connection( mats_along_path[ilink],
        //TEST mats_along_path[ilink+1], mat_tmp0, dir  );
        link_transport_connection_hisq( mats_along_path[ilink], internal_U_link,
        mats_along_path[ilink+1], mat_tmp0, dir  );
      }
    } // end loop over links



    /* A path has (length+1) points, counting the ends.  At first
	 point, no "down" direction links have their momenta "at this
	 point". At last, no "up" ... */
    if( GOES_FORWARDS(this_path->dir[length-1]) ) k=length-1; else k=length;
    for( ilink=0; ilink<=k; ilink++ ){
      if(ilink<length)dir = this_path->dir[ilink];
      else dir=NODIR;

//      coeff = ferm_epsilon*this_path->coeff;
      coeff = this_path->coeff;
      if( (ilink%2)==1 )coeff = -coeff;


      /* add in contribution to the force */

      if( ilink<length && GOES_FORWARDS(dir) ){

        link_gather_connection_hisq( oprod_along_path[length-ilink-1], 
				     mat_tmp1, NULL, dir );

        if(ilink==0) FORALLSITES(i,s){
          mat_tmp0[i] = mat_tmp1[i]; 
        }
        else{
          FORALLSITES(i,s){
//AB            mult_su3_na( &(mat_tmp1[i]),  &(mats_along_path[ilink][i]), &(mat_tmp0[i]) );
            mult_su3_na( &( mat_tmp1[i] ),  &( mats_along_path[ilink][i] ), &(mat_tmp0[i]) );
          }
        }

/*
        FOREVENSITES(i,s){
	    scalar_mult_add_su3_matrix( &(force_accum[dir][i]), &(mat_tmp0[i]),
                coeff, &(force_accum[dir][i]) );
	}
	FORODDSITES(i,s){
	    scalar_mult_add_su3_matrix( &(force_accum[dir][i]), &(mat_tmp0[i]),
                -coeff, &(force_accum[dir][i]) );
        }
*/
        FORALLSITES(i,s){
	    scalar_mult_add_su3_matrix( &(force_accum[dir][i]), &(mat_tmp0[i]),
                coeff, &(force_accum[dir][i]) );
	}

      }
      if( ilink>0 && GOES_BACKWARDS(lastdir) ){
	odir = OPP_DIR(lastdir);
        if( ilink==1 )FORALLSITES(i,s){
          mat_tmp0[i]=oprod_along_path[length-ilink][i];
          su3_adjoint ( &( mat_tmp0[i] ), &( mat_tmp1[i] ) );
        }
        else{
          link_gather_connection_hisq( mats_along_path[ilink-1], mat_tmp1, 
				       NULL, odir );
          FORALLSITES(i,s){
            mult_su3_na( &( oprod_along_path[length-ilink][i] ), &( mat_tmp1[i] ), &( mat_tmp0[i] ) );
            su3_adjoint ( &( mat_tmp0[i] ), &( mat_tmp1[i] ) );
          }
        }

/*
        FOREVENSITES(i,s){
          scalar_mult_add_su3_matrix( &(force_accum[odir][i]), &(mat_tmp1[i]),
                +coeff, &(force_accum[odir][i]) );
        }
        FORODDSITES(i,s){
          scalar_mult_add_su3_matrix( &(force_accum[odir][i]), &(mat_tmp1[i]),
                -coeff, &(force_accum[odir][i]) );
        }
*/
        FORALLSITES(i,s){
          scalar_mult_add_su3_matrix( &(force_accum[odir][i]), &(mat_tmp1[i]),
                +coeff, &(force_accum[odir][i]) );
        }

      }
      lastdir = dir;
    } /* end loop over links in path */
  } /* end loop over paths */

//#ifndef AB_OLD_UNIV
//} /* end loop over link for f(0)_\nu */
//#endif /* AB_OLD_UNIV */


	
  free( mat_tmp0 );
  free( mat_tmp1 );
  free( ahmat_tmp );
  for(i=0;i<=MAX_PATH_LENGTH;i++){
     free( oprod_along_path[i] );
  }
  for(i=1;i<=MAX_PATH_LENGTH;i++){
     free( mats_along_path[i] );
  }
}


//AB TEMP
// Contribution to force due to reunitarization of links: W=V*(V^+*V)^-1/2
static void 
fn_fermion_force_multi_reunit( su3_matrix *force_accum[4], 
			       su3_matrix *force_accum_old[4], 
			       ferm_links_t *fn, ks_action_paths *ap )
{
  int i, k, l, m, n, dir;
  register site *s;
  complex force, ftmp;
  su3_matrix tmat;
  su3_tensor4 dwdv, dwdagdv, dwdvre;
  su3_tensor4 dwdvs, dwdagdvs;
  su3_matrix **internal_V_link = fn->hl.V_link;
  su3_matrix **internal_Y_link = fn->hl.Y_unitlink;
  
  /* phases out */
  custom_rephase( internal_V_link, OFF, &fn->hl.phases_in_V );
  
  FORALLSITES(i,s){
    for(dir=XUP;dir<=TUP;dir++){
      // calculate derivatives
      if(ap->umethod == UNITARIZE_ROOT )
	su3_unit_der( &( internal_V_link[dir][i] ), &dwdv, &dwdagdv );
      else if(ap->umethod == UNITARIZE_RATIONAL )
	su3_unit_der_rational( &( internal_V_link[dir][i] ), &dwdv, &dwdagdv );
      else{
	node0_printf("Unknown unitarization method\n");
	terminate(1);
      }
      
      /* piece of derivative from projecting U(3) to SU(3) */
      //      su3_unitarize( &( internal_V_link[dir][i] ), &tmat ); //AB here it is done for the second time!!!
      su3_adjoint( &( internal_Y_link[dir][i] ), &tmat );
      su3_unit_der_spec( &( internal_V_link[dir][i] ), &( internal_Y_link[dir][i] ), &tmat,
                         &dwdv, &dwdagdv, &dwdvs, &dwdagdvs );
      
//      su3_spec_unit_der( &( internal_V_link[dir][i] ), &dwdvs, &dwdagdvs ); // SAME??

// this will be used for UNITARIZE_HISQ, later
//      su3_unit_der_rational( &( internal_V_link[dir][i] ), &dwdv, &dwdagdv );


//      su3_unit_der_reim( &( internal_V_link[dir][i] ), &dwdvre, &dwdvim, &dwdagdvre, &dwdagdvim );
//      su3_unit_der_reim_join( &dwdvre, &dwdvim, &dwdagdvre, &dwdagdvim, &dwdv, &dwdagdv );

      if( (dir==AB_OUT_ON_LINK) && (i==AB_OUT_ON_SITE) && (1==0) ) {
        printf("TENSOR4 dump\n");
        dumptensor4( &dwdv );
        printf("TENSOR4 adjoint dump\n");
        dumptensor4( &dwdagdv );
        /* for a unitary matrix dwdagdv should be a direct product of U^+ U^+ */
        su3_adjoint( &( internal_V_link[dir][i] ), &tmat );
        for( k=0; k<3; k++) {
          for( l=0; l<3; l++) {
            for( m=0; m<3; m++) {
              for( n=0; n<3; n++) {
                CMUL( tmat.e[k][m], tmat.e[n][l], dwdvre.t4[k][m][n][l] );
                dwdvre.t4[k][m][n][l].real /= -2;
                dwdvre.t4[k][m][n][l].imag /= -2;
              }
            }
          }
        }
        printf("TENSOR4 adjoint check dump\n");
        dumptensor4( &dwdvre );
      }


      // adjoint piece of force from the previous level
      su3_adjoint( &( force_accum_old[dir][i] ), &tmat );

      // contraction of derivative (rank 4 tensor) with previous force (matrix),
      // the rule: RESULT_{mn} = TENSOR_{kmnl} * MATRIX_{kl}
      // TODO: EXPLAIN indexing convention (transposed to Wong, Woloshyn)
      // TODO: MATCH it to operations in su3_mat_op.c
      // TODO: this may be separated into a routine (if used not only here!)
      clear_su3mat( &( force_accum[dir][i] ) );
      for( m=0; m<3; m++) {
        for( n=0; n<3; n++) {
          force.real=0.0;
          force.imag=0.0;
          for( k=0; k<3; k++) {
            for( l=0; l<3; l++) {
              CMUL( dwdvs.t4[k][m][n][l], force_accum_old[dir][i].e[l][k], ftmp );
//              CSUM( force, ftmp );
              CSUM( force_accum[dir][i].e[n][m], ftmp );
              /* CAREFUL with the adjoint part here! */
              CMUL( dwdagdvs.t4[k][m][n][l], tmat.e[l][k], ftmp );
//              CMUL( dwdagdv.t4[k][m][n][l], force_accum_old[dir][i].e[k][l], ftmp );
//              CSUM( force, ftmp );
              CSUM( force_accum[dir][i].e[n][m], ftmp );
            }
          }
//          force_accum[dir][i].e[m][n] = force;
//          force_accum[dir][i].e[m][n] = force_accum_old[dir][i].e[m][n];
        }
      }
if( (dir==AB_OUT_ON_LINK) && (i==AB_OUT_ON_SITE) ) {
  printf("REUNIT FORCE\n");
  dumpmat( &(force_accum[AB_OUT_ON_LINK][AB_OUT_ON_SITE]) );
}
    }
  }

  /* phases in */
  custom_rephase( internal_V_link, ON, &fn->hl.phases_in_V );
}





//AB WRAPPER_FORCE
static void 
fn_fermion_force_multi_wrapper( Real eps, Real *residues, 
				su3_vector **multi_x, int nterms,
				ferm_links_t *fn, ks_action_paths *ap)
{
  su3_matrix *force_accum_0[4];  // accumulate force, smearing 0
  su3_matrix *force_accum_1[4];  // accumulate force, smearing 1
  su3_matrix *force_accum_1u[4];  // accumulate force, reunitarization
  su3_matrix *force_accum_2[4];  // accumulate force, smearing 2
  su3_matrix *force_final[4];  // accumulate final force, multiply by time step
  su3_matrix **U_link = fn->hl.U_link;
  su3_matrix **W_unitlink = fn->hl.W_unitlink;
  anti_hermitmat *ahmat_tmp;
  register site *s;
  int dir;
  int i;
  int ipath;

  int num_q_paths_1 = ap->p1.num_q_paths;
  Q_path *q_paths_1 = ap->p1.q_paths;
  int num_q_paths_2 = ap->p2.num_q_paths;
  Q_path *q_paths_2 = ap->p2.q_paths;


  if( first_force==1 ){
    if( q_paths_sorted_1==NULL ) q_paths_sorted_1 = (Q_path *)malloc( num_q_paths_1*sizeof(Q_path) );
    if(netbackdir_table_1==NULL ) netbackdir_table_1 = (int *)malloc( num_q_paths_1*sizeof(int) );
    if( q_paths_sorted_2==NULL ) q_paths_sorted_2 = (Q_path *)malloc( num_q_paths_2*sizeof(Q_path) );
    if(netbackdir_table_2==NULL ) netbackdir_table_2 = (int *)malloc( num_q_paths_2*sizeof(int) );
    else{ node0_printf("WARNING: remaking sorted path tables\n"); exit(0); }

// Note change for HISQ - using only one "first_force", so have to allocate
// all tables here
    sort_quark_paths( q_paths_1, q_paths_sorted_1, num_q_paths_1, 8 );
    for( ipath=0; ipath<num_q_paths_1; ipath++ )
	netbackdir_table_1[ipath] = find_backwards_gather( &(q_paths_sorted_1[ipath]) );
    sort_quark_paths( q_paths_2, q_paths_sorted_2, num_q_paths_2, 16 );
    for( ipath=0; ipath<num_q_paths_2; ipath++ )
	netbackdir_table_2[ipath] = find_backwards_gather( &(q_paths_sorted_2[ipath]) );
    first_force=0;
  }


  for(i=XUP;i<=TUP;i++){
     force_accum_0[i] = (su3_matrix *) malloc(sites_on_node*sizeof(su3_matrix) );
     force_accum_1[i] = (su3_matrix *) malloc(sites_on_node*sizeof(su3_matrix) );
     force_accum_1u[i] = (su3_matrix *) malloc(sites_on_node*sizeof(su3_matrix) );
     force_accum_2[i] = (su3_matrix *) malloc(sites_on_node*sizeof(su3_matrix) );
     force_final[i] = (su3_matrix *) malloc(sites_on_node*sizeof(su3_matrix) );
  }
  ahmat_tmp = (anti_hermitmat *) special_alloc(sites_on_node*sizeof(anti_hermitmat) );

  printf("UNITARIZATION_METHOD=%d\n",ap->umethod);
  //AB smearing level 0
  fn_fermion_force_multi_univ0( eps, residues, multi_x, nterms, force_accum_0 );
  //AB smearing level 2
  fn_fermion_force_multi_universal( eps, residues, multi_x, nterms, force_accum_2, force_accum_0,
     W_unitlink, num_q_paths_2, q_paths_sorted_2, netbackdir_table_2 );

  switch(ap->umethod){
  case UNITARIZE_NONE:
    //AB smearing level 1
    fn_fermion_force_multi_universal( eps, residues, multi_x, nterms, force_accum_1, force_accum_2,
				      U_link, num_q_paths_1, q_paths_sorted_1, netbackdir_table_1 );
    break;
  case UNITARIZE_ROOT:
    printf("UNITARIZE_ROOT call in fermion_force_fn_multi\n");
    //AB reunitarization
    fn_fermion_force_multi_reunit( force_accum_1u, force_accum_2, fn, ap );
    //AB smearing level 1
    fn_fermion_force_multi_universal( eps, residues, multi_x, nterms, force_accum_1, force_accum_1u,
				      U_link, num_q_paths_1, q_paths_sorted_1, netbackdir_table_1 );
    break;
  case UNITARIZE_RATIONAL:
    // ON THIS LEVEL THIS IS IDENTICAL TO UNITARIZE_ROOT
    printf("UNITARIZE_RATIONAL call in fermion_force_fn_multi\n");
    //AB reunitarization
    fn_fermion_force_multi_reunit( force_accum_1u, force_accum_2, fn, ap );
    //AB smearing level 1
    fn_fermion_force_multi_universal( eps, residues, multi_x, nterms, force_accum_1, force_accum_1u,
				      U_link, num_q_paths_1, q_paths_sorted_1, netbackdir_table_1 );
    break;
  default:
    node0_printf("Unknown unitarization method\n");
    terminate(1);
  }


#ifdef AB_DEBUG
  fn_fermion_force_multi_univ0( eps, residues, multi_x, nterms, force_accum_0 );
  fn_fermion_force_multi_reunit( force_accum_1, force_accum_0, fn, ap );
//  fn_fermion_force_multi_det( force_accum_1, force_accum_0, fn );
#endif /* AB_DEBUG */


//  custom_rephase( U_link, OFF, &fn->hl.phases_in_U );

  //AB contraction with the link in question should be done here,
  //AB after contributions from all levels of smearing are taken into account
  for(dir=XUP;dir<=TUP;dir++){
    FORALLSITES(i,s){
      mult_su3_nn( &( s->link[dir] ), &( force_accum_1[dir][i] ), &( force_final[dir][i] ) );

//#define ABcheck
#ifdef ABcheck
      complex ft1, ft2;
      su3_matrix tt1, tt2;
      ft1.real = 0.0;
      ft1.imag = 0.0;
      su3mat_copy( &( s->link[dir] ), &tt1 );
      for( k=0; k<3; k++) {
        for( l=0; l<3; l++) {
          CMUL( tt1.e[k][l], force_accum_0[dir][i], ft2 );
          CSUM( ft1, ft2 );
        }
      }
#endif /* ABcheck */



#ifdef ABcheck_unit
// THIS CHECK IS INTENDED TO PERFORM THE SAME OPERATIONS
// AS CALCULATING THE FORCE FROM UNITARIZATION OF A MATRIX
// WHICH IS ALREADY UNITARY
      su3_adjoint( &( force_accum_0[dir][i] ), &( tmp2 ) );
      mult_su3_an( &( U_link[dir][i] ), &( tmp2 ), &( tmp_su3mat ) );
      mult_su3_na( &( tmp_su3mat ), &( U_link[dir][i] ), &( tmp2 ) );
      sub_su3_matrix( &( force_accum_0[dir][i] ), &( tmp2 ), &( tmp_su3mat ) );
      scalar_mult_su3_matrix( &( tmp_su3mat ), 0.5, &( force_accum_2[dir][i] ) );
if( (dir==AB_OUT_ON_LINK) && (i==AB_OUT_ON_SITE) ) {
  printf("WRAPPER Intermediate/Reunit FORCE\n");
  dumpmat( &(force_accum_2[AB_OUT_ON_LINK][AB_OUT_ON_SITE]) );
}
      mult_su3_nn( &( s->link[dir] ), &( force_accum_2[dir][i] ), &( force_final[dir][i] ) );
#endif /* ABcheck_unit */

    }
  }

//  custom_rephase( U_link, ON, &fn->hl.phases_in_U );

  //AB multiply by the time step "eps" twice (only forward paths are considered)
  //AB take into account even/odd parity (it is NOT done in "universal" routine
  for(dir=XUP;dir<=TUP;dir++){
/*
    FORALLSITES(i,s){
      scalar_mult_su3_matrix( &(force_final[dir][i]), 2.0*eps, &(force_final[dir][i]) );
    }
*/
    FOREVENSITES(i,s){
      scalar_mult_su3_matrix( &(force_final[dir][i]), 2.0*eps, &(force_final[dir][i]) );
    }
    FORODDSITES(i,s){
      scalar_mult_su3_matrix( &(force_final[dir][i]), -2.0*eps, &(force_final[dir][i]) );
    }
  }

  printf("WRAPPER FORCE\n");
  dumpmat( &(force_final[AB_OUT_ON_LINK][AB_OUT_ON_SITE]) );
  printf("WRAPPER FORCE ANTIHERMITIAN\n");
  for(dir=XUP; dir<=TUP; dir++) FORALLSITES(i,s){
    make_anti_hermitian( &(force_final[dir][i]), &(ahmat_tmp[i]) );
    uncompress_anti_hermitian( &(ahmat_tmp[i]), &(force_final[dir][i]) );
  }
  dumpmat( &(force_final[AB_OUT_ON_LINK][AB_OUT_ON_SITE]) );


  //AB: here TA projection and addition to momenta should be done
  //AB: once force is calculated correctly
//  printf("*** Evolve with WRAPPER force ***\n");
//  for(dir=XUP; dir<=TUP; dir++)FORALLSITES(i,s){
//     uncompress_anti_hermitian( &(s->mom[dir]), &tmp_su3mat );
//     add_su3_matrix( &tmp_su3mat, &(force_final[dir][i]), &tmp_su3mat );
//     make_anti_hermitian( &tmp_su3mat, &(s->mom[dir]) );
//  }


  for(i=XUP;i<=TUP;i++){
     free( force_accum_0[i] );
     free( force_accum_1[i] );
     free( force_accum_1u[i] );
     free( force_accum_2[i] );
     free( force_final[i] );
  }
  free(ahmat_tmp);
}



void custom_rephase( su3_matrix **internal_links, int flag, int* status_now ) {
  register int i,j,k,dir;
  register site *s;

  // s->phase assumed to be set at this stage (!)
  // check to see if settings are consistent
  if( !( flag==ON && *status_now==OFF )  && !( flag==OFF && *status_now==ON ) ){
    node0_printf("DUMMY: you messed up the phase\n");
    //node0_printf("flag = %d, status=%d, ON=%d, address = %x\n",flag,*status_now,ON,status_now);
    terminate(0);
  }
  FORALLSITES(i,s){
    for(dir=XUP;dir<=TUP;dir++){
      for(j=0;j<3;j++)for(k=0;k<3;k++){
        (internal_links[dir][i].e[j][k]).real *= s->phase[dir];
        (internal_links[dir][i].e[j][k]).imag *= s->phase[dir];
      }
    }
  }
  *status_now=flag;   // keep track of the current phases
}


//TEMP
static void 
fn_fermion_force_multi_hisq( Real eps, Real *residues, su3_vector **multi_x, 
			     int nterms, ferm_links_t *fn, ks_action_paths *ap)
{
  // note CG_solution and Dslash * solution are combined in "multi_x" 
  // Use forward part of Dslash to get force 
  // see long comment at end 
  // This version is temporary.  For two stage smearing, includes unitarization step

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
  Q_path *this_path;	// pointer to current path
  Q_path *last_path;	// pointer to previous path
  msg_tag *mtag[2];
  su3_matrix *mat_tmp0,*mat_tmp1;
  su3_matrix *oprod_along_path[MAX_PATH_LENGTH+1]; // length N path has N+1 sites!!
  su3_matrix *mats_along_path[MAX_PATH_LENGTH+1]; // 
  su3_matrix *force_accum_1[4];  // accumulate force through smearing 1
  su3_matrix *force_accum_2[4];  // accumulate force through smearing 2 (done first)
  int netbackdir, last_netbackdir;	// backwards direction for entire path
  su3_matrix **U_link = fn->hl.U_link;
  su3_matrix **V_link = fn->hl.V_link;
  su3_matrix **W_unitlink = fn->hl.W_unitlink;
  su3_matrix **Y_unitlink = fn->hl.Y_unitlink;

  su3_tensor4 dydv, dydagdv, dwdv, dwdagdv;
  su3_matrix force_tmp;
  int m, n, l;
  complex force, ftmp;

#ifdef FFTIME
  int nflop = 0 + 0*nterms; // HISQ action 10/30/07 version of code;
  double dtime;
#endif
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
  }
  mat_tmp0 = (su3_matrix *) special_alloc(sites_on_node*sizeof(su3_matrix) );
  mat_tmp1 = (su3_matrix *) special_alloc(sites_on_node*sizeof(su3_matrix) );
  if( mat_tmp1 == NULL ){printf("Node %d NO ROOM\n",this_node); exit(0); }

#ifdef FFTIME
	dtime=-dclock();
#endif
	
  ferm_epsilon_1 = 2.0*eps; // we only do forward paths, gives factor of 2
  ferm_epsilon_2 =     1.0; // both forward and backward  CHECK THIS!!
  if( first_force==1 ){
    if( q_paths_sorted_1==NULL ) q_paths_sorted_1 = (Q_path *)malloc( num_q_paths_1*sizeof(Q_path) );
    if(netbackdir_table_1==NULL ) netbackdir_table_1 = (int *)malloc( num_q_paths_1*sizeof(int) );
    if( q_paths_sorted_2==NULL ) q_paths_sorted_2 = (Q_path *)malloc( num_q_paths_2*sizeof(Q_path) );
    if(netbackdir_table_2==NULL ) netbackdir_table_2 = (int *)malloc( num_q_paths_2*sizeof(int) );
    else{ node0_printf("WARNING: remaking sorted path tables\n"); exit(0); }

// Note change for HISQ - using only one "first_force", so have to allocate
// all tables here
    sort_quark_paths( q_paths_1, q_paths_sorted_1, num_q_paths_1, 8 );
    for( ipath=0; ipath<num_q_paths_1; ipath++ )
	netbackdir_table_1[ipath] = find_backwards_gather( &(q_paths_sorted_1[ipath]) );
    sort_quark_paths( q_paths_2, q_paths_sorted_2, num_q_paths_2, 16 );
    for( ipath=0; ipath<num_q_paths_2; ipath++ )
	netbackdir_table_2[ipath] = find_backwards_gather( &(q_paths_sorted_2[ipath]) );
    first_force=0;
  }

  // clear force accumulators
  for(dir=XUP;dir<=TUP;dir++)FORALLSITES(i,s){
    clear_su3mat( &(force_accum_1[dir][i]) );
    clear_su3mat( &(force_accum_2[dir][i]) );
  }

  //BEGIN SMEARING TWO PART
  /* loop over paths, and loop over links in path */
  last_netbackdir = NODIR;
  last_path = NULL;
  for( ipath=0; ipath<num_q_paths_2; ipath++ ){
    this_path = &(q_paths_sorted_2[ipath]);
    if(this_path->forwback== -1)continue;	/* skip backwards dslash - only in smearing two*/

    length = this_path->length;
    // find gather to bring multi_x[term] from "this site" to end of path
    netbackdir = netbackdir_table_2[ipath];
    // and bring multi_x to end - no gauge transformation !!
    // resulting outer product matrix has gauge transformation properties of a connection
    // from start to end of path
    if( netbackdir != last_netbackdir){ // don't need to repeat this if same net disp. as last path
	k=0; // which gather we are using
        mtag[k] = start_gather_field( multi_x[0], sizeof(su3_vector),
           netbackdir, EVENANDODD, gen_pt[k] );
        FORALLSITES(i,s){
	  clear_su3mat(  &oprod_along_path[0][i] ); // actually last site in path
        }
        for(term=0;term<nterms;term++){
          if(term<nterms-1)mtag[1-k] = start_gather_field( multi_x[term+1],
		sizeof(su3_vector), netbackdir, EVENANDODD, gen_pt[1-k] );
          wait_gather(mtag[k]);
          FORALLSITES(i,s){
	    su3_projector( &multi_x[term][i], (su3_vector *)gen_pt[k][i], &tmat );
	    scalar_mult_add_su3_matrix( &oprod_along_path[0][i], &tmat, residues[term], &oprod_along_path[0][i] );
          }
          cleanup_gather(mtag[k]);
	  k=1-k; // swap 0 and 1
        } /* end loop over terms in rational function expansion */

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
      while ( j>0 && this_path->dir[j] == last_path->dir[j+last_path->length-length] ) j--;
    if( GOES_BACKWARDS(this_path->dir[0]) ) k=1; else k=0;

    for(ilink=j;ilink>=k;ilink--){

      link_transport_connection_hisq( oprod_along_path[length-ilink-1], W_unitlink,
      oprod_along_path[length-ilink], mat_tmp0, this_path->dir[ilink]  );
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
            mtag[1] = start_gather_field( W_unitlink[dir], sizeof(su3_matrix),
                   OPP_DIR(dir), EVENANDODD, gen_pt[1] );
            wait_gather(mtag[1]);
            FORALLSITES(i,s){ su3_adjoint( (su3_matrix *)gen_pt[1][i], &(mats_along_path[1][i]) ); }
            cleanup_gather(mtag[1]);
          }
          else{
            FORALLSITES(i,s){ mats_along_path[1][i] = W_unitlink[OPP_DIR(dir)][i]; }
          }
      }
      else { // ilink != 0
        dir = OPP_DIR(this_path->dir[ilink]);
        link_transport_connection_hisq( mats_along_path[ilink], W_unitlink,
        mats_along_path[ilink+1], mat_tmp0, dir  );
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
	   FORALLSITES(i,s){ mat_tmp0[i] = mat_tmp1[i]; }
	}
	else {
          FORALLSITES(i,s){
	     mult_su3_na( &(mat_tmp1[i]),  &(mats_along_path[ilink][i]), &(mat_tmp0[i]) );
          }
	}
        FOREVENSITES(i,s){
	    scalar_mult_add_su3_matrix( &(force_accum_2[dir][i]), &(mat_tmp0[i]),
                coeff, &(force_accum_2[dir][i]) );
	}
	FORODDSITES(i,s){
	    scalar_mult_add_su3_matrix( &(force_accum_2[dir][i]), &(mat_tmp0[i]),
                //-coeff, &(force_accum_2[dir][i]) );
                coeff, &(force_accum_2[dir][i]) ); // Alexei 9/20/07
        }
      }
      if( ilink>0 && GOES_BACKWARDS(lastdir) ){
	odir = OPP_DIR(lastdir);
        if ( ilink==1 ){
	    //FORALLSITES(i,s){ mat_tmp0[i] = oprod_along_path[length-ilink][i]; }
	    FORALLSITES(i,s){ su3_adjoint( &(oprod_along_path[length-ilink][i]),
		&( mat_tmp0[i] ) ); }
	}
	else {
	    link_gather_connection_hisq( mats_along_path[ilink-1], mat_tmp1, 
					 NULL, odir );
            FORALLSITES(i,s){
              //mult_su3_na( &(oprod_along_path[length-ilink][i]),  &(mat_tmp1[i]), &(mat_tmp0[i]) );
              mult_su3_na( &(mat_tmp1[i]), &(oprod_along_path[length-ilink][i]),  &(mat_tmp0[i]) );
            }
	}
        FOREVENSITES(i,s){
	    scalar_mult_add_su3_matrix( &(force_accum_2[odir][i]), &(mat_tmp0[i]),
                //-coeff, &(force_accum_2[odir][i]) );
                +coeff, &(force_accum_2[odir][i]) );
	}
        FORODDSITES(i,s){
	    scalar_mult_add_su3_matrix( &(force_accum_2[odir][i]), &(mat_tmp0[i]),
                //coeff, &(force_accum_2[odir][i]) );
                coeff, &(force_accum_2[odir][i]) ); // Alexei 9/20/07
                //-coeff, &(force_accum_2[odir][i]) );
	}
      }

      lastdir = dir;
    } /* end loop over links in path */
    last_netbackdir = netbackdir;
    last_path = &(q_paths_sorted_2[ipath]);
  } /* end loop over paths */

  // END SMEARING TWO PART

  //UNITARIZATION PART

  switch(ap->umethod){
  case UNITARIZE_ROOT:
  case UNITARIZE_RATIONAL:
  case UNITARIZE_ANALYTIC:

    // WARNING: custom_rephase IS A QUICK FIX, IT DOES NOT SET ANY GLOBAL VARIABLES!
    // rephase (out) V, Y
    custom_rephase( V_link, OFF, &fn->hl.phases_in_V );
    custom_rephase( Y_unitlink, OFF, &fn->hl.phases_in_Y );
    
    FORALLSITES(i,s){
      for(dir=XUP;dir<=TUP;dir++){
	// derivatives: dY/dV, dY^+/dV
	switch (ap->umethod){
	case UNITARIZE_ROOT:
	  su3_unit_der( &( V_link[dir][i] ), &dydv, &dydagdv );
	  break;
	case UNITARIZE_RATIONAL:
	  su3_unit_der_rational( &( V_link[dir][i] ), &dydv, &dydagdv );
	  break;
	case UNITARIZE_ANALYTIC:
	  su3_unit_der_analytic( &( V_link[dir][i] ), &dydv, &dydagdv );
	  break;
	default:
	  node0_printf("Unknown unitarization method\n");
	  terminate(1);
	}
	
	// piece of derivative from projecting U(3) to SU(3)
	su3_adjoint( &( Y_unitlink[dir][i] ), &tmat );
	su3_unit_der_spec( &( V_link[dir][i] ), 
			   &( Y_unitlink[dir][i] ), &tmat,
			   &dydv, &dydagdv, &dwdv, &dwdagdv );
	
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
	    force.real=0.0;
	    force.imag=0.0;
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
	// PUT THE RESULT INTO force_accum_2, THIS MAY NEED CHANGES
	su3mat_copy( &force_tmp, &( force_accum_2[dir][i] ) );
      }
    }
    
    // WARNING: custom_rephase IS A QUICK FIX, IT DOES NOT SET ANY GLOBAL VARIABLES!
    // rephase (in) V, Y
    custom_rephase( V_link, ON, &fn->hl.phases_in_V );
    custom_rephase( Y_unitlink, ON, &fn->hl.phases_in_Y );
    break;

  case UNITARIZE_NONE:
    //  printf("HISQ: Proceeding with UNITARIZE_NONE\n");
    break;
  default:
    printf("HISQ: Only UNITARIZE_ROOT and UNITARIZE_RATIONAL supported so far\n");
    printf("HISQ: Proceeding without reunitarization\n");
  } /* ap->umethod */



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
            mtag[1] = start_gather_field( U_link[dir], sizeof(su3_matrix),
                   OPP_DIR(dir), EVENANDODD, gen_pt[1] );
            wait_gather(mtag[1]);
            FORALLSITES(i,s){ su3_adjoint( (su3_matrix *)gen_pt[1][i], &(mats_along_path[1][i]) ); }
            cleanup_gather(mtag[1]);
          }
          else{
            FORALLSITES(i,s){ mats_along_path[1][i] = U_link[OPP_DIR(dir)][i]; }
          }
      }
      else { // ilink != 0
        dir = OPP_DIR(this_path->dir[ilink]);
        link_transport_connection_hisq( mats_along_path[ilink], U_link,
        mats_along_path[ilink+1], mat_tmp0, dir  );
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

      if(ilink==0 && GOES_FORWARDS(dir) ) FORALLSITES(i,s){
	mat_tmp0[i] = oprod_along_path[length][i]; 
      }
      else if( ilink>0) FORALLSITES(i,s){
        mult_su3_na( &(oprod_along_path[length-ilink][i]),  &(mats_along_path[ilink][i]), &(mat_tmp0[i]) );
      }

      /* add in contribution to the force */
      if( ilink<length && GOES_FORWARDS(dir) ){
        FOREVENSITES(i,s){
	    scalar_mult_add_su3_matrix( &(force_accum_1[dir][i]), &(mat_tmp0[i]),
                coeff, &(force_accum_1[dir][i]) );
	}
	FORODDSITES(i,s){
	    scalar_mult_add_su3_matrix( &(force_accum_1[dir][i]), &(mat_tmp0[i]),
                -coeff, &(force_accum_1[dir][i]) );
        }
      }
      if( ilink>0 && GOES_BACKWARDS(lastdir) ){
	odir = OPP_DIR(lastdir);
        FOREVENSITES(i,s){
	    scalar_mult_add_su3_matrix( &(force_accum_1[odir][i]), &(mat_tmp0[i]),
                -coeff, &(force_accum_1[odir][i]) );
	}
        FORODDSITES(i,s){
	    scalar_mult_add_su3_matrix( &(force_accum_1[odir][i]), &(mat_tmp0[i]),
                coeff, &(force_accum_1[odir][i]) );
	}
      }

      lastdir = dir;
    } /* end loop over links in path */
    last_netbackdir = netbackdir;
    last_path = &(q_paths_sorted_1[ipath]);
  } /* end loop over paths */

  // END SMEARING ONE PART

  /* Put antihermitian traceless part into momentum */
  // add force to momentum
  for(dir=XUP; dir<=TUP; dir++)FORALLSITES(i,s){
     uncompress_anti_hermitian( &(s->mom[dir]), &tmat2 );
     add_su3_matrix( &tmat2, &(force_accum_1[dir][i]), &tmat2 );
     make_anti_hermitian( &tmat2, &(s->mom[dir]) );
  }
	
//printf("FF flops = %d\n",tempflops);
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
  }


#ifdef FFTIME
  dtime += dclock();
  node0_printf("FFTIME:  time = %e (HISQ) terms = %d mflops = %e\n",dtime,nterms,
	     (Real)nflop*volume/(1e6*dtime*numnodes()) );
#endif
} //fn_fermion_force_multi_hisq

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
   if( disp[XUP]==+1 && disp[YUP]== 0 && disp[ZUP]== 0 && disp[TUP]== 0 )return(XDOWN);
   if( disp[XUP]==-1 && disp[YUP]== 0 && disp[ZUP]== 0 && disp[TUP]== 0 )return(XUP);
   if( disp[XUP]== 0 && disp[YUP]==+1 && disp[ZUP]== 0 && disp[TUP]== 0 )return(YDOWN);
   if( disp[XUP]== 0 && disp[YUP]==-1 && disp[ZUP]== 0 && disp[TUP]== 0 )return(YUP);
   if( disp[XUP]== 0 && disp[YUP]== 0 && disp[ZUP]==+1 && disp[TUP]== 0 )return(ZDOWN);
   if( disp[XUP]== 0 && disp[YUP]== 0 && disp[ZUP]==-1 && disp[TUP]== 0 )return(ZUP);
   if( disp[XUP]== 0 && disp[YUP]== 0 && disp[ZUP]== 0 && disp[TUP]==+1 )return(TDOWN);
   if( disp[XUP]== 0 && disp[YUP]== 0 && disp[ZUP]== 0 && disp[TUP]==-1 )return(TUP);

   if( disp[XUP]==+3 && disp[YUP]== 0 && disp[ZUP]== 0 && disp[TUP]== 0 )return(X3DOWN);
   if( disp[XUP]==-3 && disp[YUP]== 0 && disp[ZUP]== 0 && disp[TUP]== 0 )return(X3UP);
   if( disp[XUP]== 0 && disp[YUP]==+3 && disp[ZUP]== 0 && disp[TUP]== 0 )return(Y3DOWN);
   if( disp[XUP]== 0 && disp[YUP]==-3 && disp[ZUP]== 0 && disp[TUP]== 0 )return(Y3UP);
   if( disp[XUP]== 0 && disp[YUP]== 0 && disp[ZUP]==+3 && disp[TUP]== 0 )return(Z3DOWN);
   if( disp[XUP]== 0 && disp[YUP]== 0 && disp[ZUP]==-3 && disp[TUP]== 0 )return(Z3UP);
   if( disp[XUP]== 0 && disp[YUP]== 0 && disp[ZUP]== 0 && disp[TUP]==+3 )return(T3DOWN);
   if( disp[XUP]== 0 && disp[YUP]== 0 && disp[ZUP]== 0 && disp[TUP]==-3 )return(T3UP);
node0_printf("OOOPS: NODIR\n"); exit(0);
   return( NODIR );
} //find_backwards_gather

// Make a new path table.  Sorted principally by total displacement of path.
// Below that, sort by direction of first link
// Below that, sort by direction of second link - note special case of one link paths
int 
sort_quark_paths( Q_path *src_table, Q_path *dest_table, int npaths, int num_back_dirs )
{
    int netdir,dir0,dir1,dir1tmp,thislength,num_new,i,j;

    num_new=0; // number of paths in sorted table
    for( i=0; i<num_back_dirs; i++ ){ // loop over net_back_dirs
        netdir = net_back_dirs[i]; // table of possible displacements for Fat-Naik
	for( dir0=0; dir0<=7; dir0++){ // XUP ... TDOWN
	  for( dir1=-1; dir1<=7; dir1++){ // NODIR, XUP ... TDOWN
	    if( dir1==-1 )dir1tmp=NODIR; else dir1tmp=dir1;
	    for( j=0; j<npaths; j++ ){ // pick out paths with right net displacement
	    thislength = src_table[j].length;
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
/*   Version for asqtad.  Parallel transport nterms source vectors    */
/**********************************************************************/

/*
 * 10/01/02, flopcount for ASQ_OPTIMIZED
 * Fermion force: 21698*nterms
 */

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
      FORALLSITES(i,s)
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
  FORALLSITES(i,s){
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
