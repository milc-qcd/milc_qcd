/****** fermion_force_fn_rhmc.c -- ******************/
/* OBSOLETE - CD */
/* MIMD version 7 */
/* Fermion force for a general quark action.  Not optimized.
* D.T. 1/28/98, starting from gauge_stuff.c
* K.O. 3/99 Added optimized fattening for Asq actions
* D.T. 4/99 Combine force calculations for both mass quarks
* K.O. 4/99 Optimized force for Asq action
* S.G. 7/01, modified to use t_longlink and t_fatlink
* C.D. 10/02, consolidated quark_stuff.c and quark_stuff_tmp.c
* T.B. 11/01, Added d(M)/d(u0) dependencies for equation of state calc's with
*             Asqtad action - ATTN: site structure needs 'dfatlink_du0' in
*             addition to 'fatlink': #define DM_DU0
*
* C.D. 3/05 Separated from quark_stuff.c
* D.T. 6/06 Adapt to rational function code

* assume all paths are one-link or three-lin (Naik_
* imples all paths connect even to odd sites, etc.
* Tabulate "backwards" paths (e.g. "XDOWN" is backward path to "XUP")
* as separate parity transforms of the fundamental paths.  They will
* generally need a negative sign in Dslash.  See bottom for a long
* comment on sign conventions.
*/

/**#define FFTIME**/
/**#define FFSTIME**/

#include "ks_imp_includes.h"	/* definitions files and prototypes */

#define GOES_FORWARDS(dir) (dir<=TUP)
#define GOES_BACKWARDS(dir) (dir>TUP)

/**********************************************************************/
/*   Version for a rational function with "nterms" "masses"           */
/**********************************************************************/

/* update the  momenta with the fermion force */
/* Assumes that the multimass conjugate gradient has been run, with the answer in
   multi_x, and dslash_site(multi_x,multi_x,ODD) has been run. (fills in multi_x on odd sites) */
/* SEE LONG COMMENTS AT END */

Q_path *q_paths_sorted = NULL;	// Quark paths sorted by net displacement and last directions
// table of gather directions to bring start of path to end, for "FN" = "fat-Naik" actions
int net_back_dirs[16] = { XDOWN, YDOWN, ZDOWN, TDOWN, XUP, YUP, ZUP, TUP, 
	X3DOWN, Y3DOWN, Z3DOWN, T3DOWN, X3UP, Y3UP, Z3UP, T3UP };
int sort_quark_paths( Q_path *src_table, Q_path *dest_table, int npaths );
int find_backwards_gather( Q_path *path );
int first_force=1;	// 1 if force hasn't been called yet

void fn_fermion_force_rhmc( Real eps, Real *residues, su3_vector **multi_x, 
			    int nterms, ferm_links_t *fn )
{
  /* note CG_solution and Dslash * solution are combined in "multi_x" */
  /* New version 1/21/99.  Use forward part of Dslash to get force */
  /* see long comment at end */
  /* For each link we need multi_x transported from both ends of path. */
  int term;
  register int i,dir,lastdir=-99,ipath,ilink;
  register site *s;
  int length;
  su3_matrix tmat,tmat2;
  Real ferm_epsilon, coeff;
  int num_q_paths = get_num_q_paths();
  Q_path *q_paths = get_q_paths();
  msg_tag *mtag1;
  su3_matrix *mat_outerprod,*mat_tmp0,*mat_tmp1,*tmp_matpt;
  int netbackdir, last_netbackdir;

#ifdef FFTIME
  int nflop = 0;
  double dtime;
#endif
node0_printf("STARTING fn_fermion_force_rhmc() nterms = %d\n",nterms);
  if( nterms==0 )return;

  mat_outerprod = (su3_matrix *) malloc(sites_on_node*sizeof(su3_matrix) );
  mat_tmp0 = (su3_matrix *) malloc(sites_on_node*sizeof(su3_matrix) );
  mat_tmp1 = (su3_matrix *) malloc(sites_on_node*sizeof(su3_matrix) );
  if( mat_tmp1 == NULL ){printf("Node %d NO ROOM\n",this_node); exit(0); }

#ifdef FFTIME
	dtime=-dclock();
#endif
	
  ferm_epsilon = 2.0*eps;
  if( first_force==1 ){
    if( q_paths_sorted==NULL ) q_paths_sorted = (Q_path *)malloc( num_q_paths*sizeof(Q_path) );
    else{ node0_printf("WARNING: remaking sorted path table\n"); exit(0); }
    sort_quark_paths( q_paths, q_paths_sorted, num_q_paths );
    first_force=0;
  }

	  /* loop over paths, and loop over links in path */
	  last_netbackdir = NODIR;
	  for( ipath=0; ipath<num_q_paths; ipath++ ){
	    if(q_paths_sorted[ipath].forwback== -1)continue;	/* skip backwards dslash */

	    length = q_paths_sorted[ipath].length;
	    // find gather to bring multi_x[term] from "this site" to end of path
	    netbackdir = find_backwards_gather( &(q_paths_sorted[ipath]) );
	    // and bring multi_x to end - no gauge transformation !!
	    // resulting outer product matrix has gauge transformation properties of a connection
	    // from start to end of path
	    if( netbackdir != last_netbackdir){ // don't need to repeat this if same net disp. as last path
	        FORALLSITES(i,s){
		  clear_su3mat(  &mat_outerprod[i] );
	        }
	        for(term=0;term<nterms;term++){
	          mtag1 = start_gather_field( multi_x[term], sizeof(su3_vector),
	             netbackdir, EVENANDODD, gen_pt[1] );
	          wait_gather(mtag1);
	          FORALLSITES(i,s){
		    su3_projector( &multi_x[term][i], (su3_vector *)gen_pt[1][i], &tmat );
		    scalar_mult_add_su3_matrix( &mat_outerprod[i], &tmat, residues[term], &mat_outerprod[i] );
	          }
	          cleanup_gather(mtag1);
	        } /* end loop over terms in rational function expansion */
	    }
	
	    /* path transport the outer product, or projection matrix, of multi_x[term]
	       (EVEN sites)  and Dslash*multi_x[term] (ODD sites) from far end.  Sometimes
		we need them at the start point of the path, and sometimes
		one link into the path - an optimization for later */
	    path_transport_connection( mat_outerprod, mat_tmp0, EVENANDODD, q_paths_sorted[ipath].dir, length );

	    /* A path has (length+1) points, counting the ends.  At first
		 point, no "down" direction links have their momenta "at this
		 point". At last, no "up" ... */
	    for( ilink=0; ilink<=length; ilink++ ){
	      if(ilink<length)dir = q_paths_sorted[ipath].dir[ilink];
	      else dir=NODIR;
	      coeff = ferm_epsilon*q_paths_sorted[ipath].coeff;
	      if( (ilink%2)==1 )coeff = -coeff;
	
	      /* add in contribution to the force */
	      /* Put antihermitian traceless part into momentum */
	      if( ilink<length && GOES_FORWARDS(dir) ){
	        FORALLSITES(i,s){
	          uncompress_anti_hermitian( &(s->mom[dir]), &tmat2 );
		  if( s->parity==EVEN ){
		    scalar_mult_add_su3_matrix(&tmat2, &(mat_tmp0[i]),  coeff, &tmat2 );
		  }
		  else{
		    scalar_mult_add_su3_matrix(&tmat2, &(mat_tmp0[i]), -coeff, &tmat2 );
		  }
	          make_anti_hermitian( &tmat2, &(s->mom[dir]) );
	        }
	      }
	      if( ilink>0 && GOES_BACKWARDS(lastdir) ){
	        FORALLSITES(i,s){
	          uncompress_anti_hermitian( &(s->mom[OPP_DIR(lastdir)]), &tmat2 );
		  if( s->parity==EVEN ){
		    scalar_mult_add_su3_matrix(&tmat2, &(mat_tmp0[i]), -coeff, &tmat2 );
		  }
		  else{
		    scalar_mult_add_su3_matrix(&tmat2, &(mat_tmp0[i]),  coeff, &tmat2 );
		  }
	          make_anti_hermitian( &tmat2, &(s->mom[OPP_DIR(lastdir)]) );
	        }
	      }

	      /* path transport multi_x[term] and Dslash*multi_x[term] to next point */
	      /* sometimes we don't need them at last point */
	      if( ilink<length-1 || 
		(ilink==length-1 && GOES_BACKWARDS(dir)) ){
	
		if( GOES_FORWARDS(dir) ){
		  FORALLSITES(i,s){
	            mult_su3_an( &(s->link[dir]), &(mat_tmp0[i]), &(tmat) );
	            mult_su3_nn( &(tmat), &(s->link[dir]), &(mat_tmp1[i]) );
		  }
		  mtag1 = start_gather_field( mat_tmp1, sizeof(su3_matrix),
	             OPP_DIR(dir), EVENANDODD, gen_pt[1] );
	          wait_gather(mtag1);
	          FORALLSITES(i,s){
		     mat_tmp0[i] = *(su3_matrix *)gen_pt[1][i];
		  }
	          cleanup_gather(mtag1);
		}
		else{   /* GOES_BACKWARDS(dir) */
	          mtag1 = start_gather_field( mat_tmp0, sizeof(su3_matrix),
	               OPP_DIR(dir), EVENANDODD, gen_pt[1] );
	          wait_gather(mtag1);
	          FORALLSITES(i,s){
	            mult_su3_nn( &(s->link[OPP_DIR(dir)]), (su3_matrix *)(gen_pt[1][i]), &(tmat) );
	            mult_su3_na( &(tmat), &(s->link[OPP_DIR(dir)]), &(mat_tmp1[i]) );
		  }
		  tmp_matpt = mat_tmp0; mat_tmp0 = mat_tmp1; mat_tmp1 = tmp_matpt;
	          cleanup_gather(mtag1);
		}
	      }
	      lastdir = dir;
	
	    } /* end loop over links in path */
	    last_netbackdir = netbackdir;
	  } /* end loop over paths */
	
  free( mat_outerprod ); free( mat_tmp0 ); free( mat_tmp1 );
#ifdef FFTIME
  dtime += dclock();
  node0_printf("FFTIME:  time = %e (%d terms) mflops = %e\n",dtime,nterms,
	     (Real)nflop*volume*nterms/(1e6*dtime*numnodes()) );
#endif
}


//version with "X" vectors in "site major" order.  Since "nterms" is variable, use
// single index array
void fn_fermion_force_rhmc_reverse( Real eps, Real *residues, 
				    su3_vector **multi_x, int nterms,
				    ferm_links_t *fn )
{
  /* note CG_solution and Dslash * solution are combined in "multi_x" */
  /* New version 1/21/99.  Use forward part of Dslash to get force */
  /* see long comment at end */
  /* For each link we need multi_x transported from both ends of path. */
  int term;
  register int i,j,dir,lastdir=-99,ipath,ilink;
  register site *s;
  int length;
  su3_matrix tmat,tmat2;
  Real ferm_epsilon, coeff;
  int num_q_paths = get_num_q_paths();
  Q_path *q_paths = get_q_paths();
  msg_tag *mtag1;
  su3_matrix *mat_outerprod,*mat_tmp0,*mat_tmp1,*tmp_matpt;
  su3_matrix *mat_opm1; // "Outer Product at Minus One from end of path"
  int netbackdir, last_netbackdir;	// backwards direction for entire path
  int pathenddir, last_pathenddir;  // final direction in path
  su3_vector *multi_x_rev;

#ifdef FFTIME
  int nflop = 0;
  double dtime;
#endif
node0_printf("STARTING fn_fermion_force_rhmc_reverse() nterms = %d\n",nterms);
  if( nterms==0 )return;

  multi_x_rev = (su3_vector *)malloc( nterms*sites_on_node*sizeof(su3_vector) );
  if( multi_x_rev==NULL ){printf("NO ROOM\n"); exit(0);}
  FORALLSITES(i,s){
        for( j=0; j<nterms; j++){
            multi_x_rev[i*nterms+j] = multi_x[j][i];
        }
  }

  mat_outerprod = (su3_matrix *) malloc(sites_on_node*sizeof(su3_matrix) );
  mat_opm1 = (su3_matrix *) malloc(sites_on_node*sizeof(su3_matrix) );
  mat_tmp0 = (su3_matrix *) malloc(sites_on_node*sizeof(su3_matrix) );
  mat_tmp1 = (su3_matrix *) malloc(sites_on_node*sizeof(su3_matrix) );
  if( mat_tmp1 == NULL ){printf("Node %d NO ROOM\n",this_node); exit(0); }

#ifdef FFTIME
	dtime=-dclock();
#endif
	
  ferm_epsilon = 2.0*eps;
  if( first_force==1 ){
    if( q_paths_sorted==NULL ) q_paths_sorted = (Q_path *)malloc( num_q_paths*sizeof(Q_path) );
    else{ node0_printf("WARNING: remaking sorted path table\n"); exit(0); }
    sort_quark_paths( q_paths, q_paths_sorted, num_q_paths );
    first_force=0;
  }

	  /* loop over paths, and loop over links in path */
	  last_netbackdir = NODIR;
	  last_pathenddir = NODIR;
	  for( ipath=0; ipath<num_q_paths; ipath++ ){
	    if(q_paths_sorted[ipath].forwback== -1)continue;	/* skip backwards dslash */

	    length = q_paths_sorted[ipath].length;
	    pathenddir = q_paths_sorted[ipath].dir[length-1];
	    // find gather to bring multi_x[term] from "this site" to end of path
	    netbackdir = find_backwards_gather( &(q_paths_sorted[ipath]) );
	    // and bring multi_x to end - no gauge transformation !!
	    // resulting outer product matrix has gauge transformation properties of a connection
	    // from start to end of path
	    if( netbackdir != last_netbackdir){ // don't need to repeat this if same net disp. as last path
	        FORALLSITES(i,s){
		  clear_su3mat(  &mat_outerprod[i] );
	        }
	        mtag1 = start_gather_field( multi_x_rev, nterms*sizeof(su3_vector),
	           netbackdir, EVENANDODD, gen_pt[1] );
	        wait_gather(mtag1);
	        FORALLSITES(i,s){
	          for(term=0;term<nterms;term++){
		    su3_projector( &multi_x_rev[nterms*i+term], ((su3_vector *)gen_pt[1][i])+term, &tmat );
		    scalar_mult_add_su3_matrix( &mat_outerprod[i], &tmat, residues[term], &mat_outerprod[i] );
	          } /* end loop over terms in rational function expansion */
	        }
	        cleanup_gather(mtag1);
	    }
	
	    /* path transport the outer product, or projection matrix, of multi_x[term]
	       (EVEN sites)  and Dslash*multi_x[term] (ODD sites) from far end.  Sometimes
		we need them at the start point of the path, and sometimes
		one link into the path - an optimization for later */
//TRY
//node0_printf("PATH, length = %d netbackdir = %d first = %d  last = %d\n",
//length,netbackdir,q_paths_sorted[ipath].dir[0],q_paths_sorted[ipath].dir[length-1]);
	    if( pathenddir != last_pathenddir || netbackdir != last_netbackdir){
	        path_transport_connection( mat_outerprod, mat_opm1, EVENANDODD,
		    &(q_paths_sorted[ipath].dir[length-1]), 1 );
//node0_printf("LAST\n");
	    }
	    path_transport_connection( mat_opm1, mat_tmp0, EVENANDODD, q_paths_sorted[ipath].dir, length-1 );
	    //path_transport_connection( mat_outerprod, mat_tmp0, EVENANDODD, q_paths_sorted[ipath].dir, length );
//ENDTRY

	    /* A path has (length+1) points, counting the ends.  At first
		 point, no "down" direction links have their momenta "at this
		 point". At last, no "up" ... */
	    for( ilink=0; ilink<=length; ilink++ ){
	      if(ilink<length)dir = q_paths_sorted[ipath].dir[ilink];
	      else dir=NODIR;
	      coeff = ferm_epsilon*q_paths_sorted[ipath].coeff;
	      if( (ilink%2)==1 )coeff = -coeff;
	
	      /* add in contribution to the force */
	      /* Put antihermitian traceless part into momentum */
	      if( ilink<length && GOES_FORWARDS(dir) ){
	        FORALLSITES(i,s){
	          uncompress_anti_hermitian( &(s->mom[dir]), &tmat2 );
		  if( s->parity==EVEN ){
		    scalar_mult_add_su3_matrix(&tmat2, &(mat_tmp0[i]),  coeff, &tmat2 );
		  }
		  else{
		    scalar_mult_add_su3_matrix(&tmat2, &(mat_tmp0[i]), -coeff, &tmat2 );
		  }
	          make_anti_hermitian( &tmat2, &(s->mom[dir]) );
	        }
	      }
	      if( ilink>0 && GOES_BACKWARDS(lastdir) ){
	        FORALLSITES(i,s){
	          uncompress_anti_hermitian( &(s->mom[OPP_DIR(lastdir)]), &tmat2 );
		  if( s->parity==EVEN ){
		    scalar_mult_add_su3_matrix(&tmat2, &(mat_tmp0[i]), -coeff, &tmat2 );
		  }
		  else{
		    scalar_mult_add_su3_matrix(&tmat2, &(mat_tmp0[i]),  coeff, &tmat2 );
		  }
	          make_anti_hermitian( &tmat2, &(s->mom[OPP_DIR(lastdir)]) );
	        }
	      }

	      /* path transport multi_x[term] and Dslash*multi_x[term] to next point */
	      /* sometimes we don't need them at last point */
	      if( ilink<length-1 || 
		(ilink==length-1 && GOES_BACKWARDS(dir)) ){
	
		if( GOES_FORWARDS(dir) ){
		  FORALLSITES(i,s){
	            mult_su3_an( &(s->link[dir]), &(mat_tmp0[i]), &(tmat) );
	            mult_su3_nn( &(tmat), &(s->link[dir]), &(mat_tmp1[i]) );
		  }
		  mtag1 = start_gather_field( mat_tmp1, sizeof(su3_matrix),
	             OPP_DIR(dir), EVENANDODD, gen_pt[1] );
	          wait_gather(mtag1);
	          FORALLSITES(i,s){
		     mat_tmp0[i] = *(su3_matrix *)gen_pt[1][i];
		  }
	          cleanup_gather(mtag1);
		}
		else{   /* GOES_BACKWARDS(dir) */
	          mtag1 = start_gather_field( mat_tmp0, sizeof(su3_matrix),
	               OPP_DIR(dir), EVENANDODD, gen_pt[1] );
	          wait_gather(mtag1);
	          FORALLSITES(i,s){
	            mult_su3_nn( &(s->link[OPP_DIR(dir)]), (su3_matrix *)(gen_pt[1][i]), &(tmat) );
	            mult_su3_na( &(tmat), &(s->link[OPP_DIR(dir)]), &(mat_tmp1[i]) );
		  }
		  tmp_matpt = mat_tmp0; mat_tmp0 = mat_tmp1; mat_tmp1 = tmp_matpt;
	          cleanup_gather(mtag1);
		}
	      }
	      lastdir = dir;
	
	    } /* end loop over links in path */
	    last_netbackdir = netbackdir;
	    last_pathenddir = pathenddir;
	  } /* end loop over paths */
	
  free( mat_outerprod ); free( mat_opm1 ); free( mat_tmp0 ); free( mat_tmp1 );
  free( multi_x_rev );
#ifdef FFTIME
  dtime += dclock();
  node0_printf("FFTIME:  time = %e (%d terms) mflops = %e\n",dtime,nterms,
	     (Real)nflop*volume*nterms/(1e6*dtime*numnodes()) );
#endif
}


// OLDER VERSION BEFORE SOME OPTIMIZATIONS.  THIS IS PROBABLY CLEARER AS TO
// WHAT IS GOING ON

void fn_fermion_force_rhmc_june05( Real eps, Real *residues, 
				   su3_vector **multi_x, int nterms,
				   ferm_links_t *fn )
{
  /* note CG_solution and Dslash * solution are combined in "multi_x" */
  /* New version 1/21/99.  Use forward part of Dslash to get force */
  /* see long comment at end */
  /* For each link we need multi_x transported from both ends of path. */
  int term;
  register int i,dir,lastdir=-99,ipath,ilink;
  register site *s;
  int length;
  su3_matrix tmat,tmat2;
  Real ferm_epsilon, coeff;
  int num_q_paths = get_num_q_paths();
  Q_path *q_paths = get_q_paths();
  msg_tag *mtag1;
  su3_matrix *mat_outerprod,*mat_tmp0,*mat_tmp1,*tmp_matpt;
  int netbackdir, last_netbackdir;

#ifdef FFTIME
  int nflop = 0;
  double dtime;
#endif
node0_printf("STARTING fn_fermion_force_rhmc() nterms = %d\n",nterms);
  if( nterms==0 )return;

  mat_outerprod = (su3_matrix *) malloc(sites_on_node*sizeof(su3_matrix) );
  mat_tmp0 = (su3_matrix *) malloc(sites_on_node*sizeof(su3_matrix) );
  mat_tmp1 = (su3_matrix *) malloc(sites_on_node*sizeof(su3_matrix) );
  if( mat_tmp1 == NULL ){printf("Node %d NO ROOM\n",this_node); exit(0); }

#ifdef FFTIME
	dtime=-dclock();
#endif
	
  ferm_epsilon = 2.0*eps;
  if( first_force==1 ){
    if( q_paths_sorted==NULL ) q_paths_sorted = (Q_path *)malloc( num_q_paths*sizeof(Q_path) );
    else{ node0_printf("WARNING: remaking sorted path table\n"); exit(0); }
    sort_quark_paths( q_paths, q_paths_sorted, num_q_paths );
    first_force=0;
  }

	  /* loop over paths, and loop over links in path */
	  last_netbackdir = NODIR;
	  for( ipath=0; ipath<num_q_paths; ipath++ ){
	    if(q_paths_sorted[ipath].forwback== -1)continue;	/* skip backwards dslash */

	    length = q_paths_sorted[ipath].length;
	    // find gather to bring multi_x[term] from "this site" to end of path
	    netbackdir = find_backwards_gather( &(q_paths_sorted[ipath]) );
	    // and bring multi_x to end - no gauge transformation !!
	    // resulting outer product matrix has gauge transformation properties of a connection
	    // from start to end of path
	    if( netbackdir != last_netbackdir){ // don't need to repeat this if same net disp. as last path
	        FORALLSITES(i,s){
		  clear_su3mat(  &mat_outerprod[i] );
	        }
	        for(term=0;term<nterms;term++){
	          mtag1 = start_gather_field( multi_x[term], sizeof(su3_vector),
	             netbackdir, EVENANDODD, gen_pt[1] );
	          wait_gather(mtag1);
	          FORALLSITES(i,s){
		    su3_projector( &multi_x[term][i], (su3_vector *)gen_pt[1][i], &tmat );
		    scalar_mult_add_su3_matrix( &mat_outerprod[i], &tmat, residues[term], &mat_outerprod[i] );
	          }
	          cleanup_gather(mtag1);
	        } /* end loop over terms in rational function expansion */
	    }
	
	    /* path transport the outer product, or projection matrix, of multi_x[term]
	       (EVEN sites)  and Dslash*multi_x[term] (ODD sites) from far end.  Sometimes
		we need them at the start point of the path, and sometimes
		one link into the path - an optimization for later */
	    path_transport_connection( mat_outerprod, mat_tmp0, EVENANDODD, q_paths_sorted[ipath].dir, length );

	    /* A path has (length+1) points, counting the ends.  At first
		 point, no "down" direction links have their momenta "at this
		 point". At last, no "up" ... */
	    for( ilink=0; ilink<=length; ilink++ ){
	      if(ilink<length)dir = q_paths_sorted[ipath].dir[ilink];
	      else dir=NODIR;
	      coeff = ferm_epsilon*q_paths_sorted[ipath].coeff;
	      if( (ilink%2)==1 )coeff = -coeff;
	
	      /* add in contribution to the force */
	      /* Put antihermitian traceless part into momentum */
	      if( ilink<length && GOES_FORWARDS(dir) ){
	        FORALLSITES(i,s){
	          uncompress_anti_hermitian( &(s->mom[dir]), &tmat2 );
		  if( s->parity==EVEN ){
		    scalar_mult_add_su3_matrix(&tmat2, &(mat_tmp0[i]),  coeff, &tmat2 );
		  }
		  else{
		    scalar_mult_add_su3_matrix(&tmat2, &(mat_tmp0[i]), -coeff, &tmat2 );
		  }
	          make_anti_hermitian( &tmat2, &(s->mom[dir]) );
	        }
	      }
	      if( ilink>0 && GOES_BACKWARDS(lastdir) ){
	        FORALLSITES(i,s){
	          uncompress_anti_hermitian( &(s->mom[OPP_DIR(lastdir)]), &tmat2 );
		  if( s->parity==EVEN ){
		    scalar_mult_add_su3_matrix(&tmat2, &(mat_tmp0[i]), -coeff, &tmat2 );
		  }
		  else{
		    scalar_mult_add_su3_matrix(&tmat2, &(mat_tmp0[i]),  coeff, &tmat2 );
		  }
	          make_anti_hermitian( &tmat2, &(s->mom[OPP_DIR(lastdir)]) );
	        }
	      }

	      /* path transport multi_x[term] and Dslash*multi_x[term] to next point */
	      /* sometimes we don't need them at last point */
	      if( ilink<length-1 || 
		(ilink==length-1 && GOES_BACKWARDS(dir)) ){
	
		if( GOES_FORWARDS(dir) ){
		  FORALLSITES(i,s){
	            mult_su3_an( &(s->link[dir]), &(mat_tmp0[i]), &(tmat) );
	            mult_su3_nn( &(tmat), &(s->link[dir]), &(mat_tmp1[i]) );
		  }
		  mtag1 = start_gather_field( mat_tmp1, sizeof(su3_matrix),
	             OPP_DIR(dir), EVENANDODD, gen_pt[1] );
	          wait_gather(mtag1);
	          FORALLSITES(i,s){
		     mat_tmp0[i] = *(su3_matrix *)gen_pt[1][i];
		  }
	          cleanup_gather(mtag1);
		}
		else{   /* GOES_BACKWARDS(dir) */
	          mtag1 = start_gather_field( mat_tmp0, sizeof(su3_matrix),
	               OPP_DIR(dir), EVENANDODD, gen_pt[1] );
	          wait_gather(mtag1);
	          FORALLSITES(i,s){
	            mult_su3_nn( &(s->link[OPP_DIR(dir)]), (su3_matrix *)(gen_pt[1][i]), &(tmat) );
	            mult_su3_na( &(tmat), &(s->link[OPP_DIR(dir)]), &(mat_tmp1[i]) );
		  }
		  tmp_matpt = mat_tmp0; mat_tmp0 = mat_tmp1; mat_tmp1 = tmp_matpt;
	          cleanup_gather(mtag1);
		}
	      }
	      lastdir = dir;
	
	    } /* end loop over links in path */
	    last_netbackdir = netbackdir;
	  } /* end loop over paths */
	
  free( mat_outerprod ); free( mat_tmp0 ); free( mat_tmp1 );
#ifdef FFTIME
  dtime += dclock();
  node0_printf("FFTIME:  time = %e (%d terms) mflops = %e\n",dtime,nterms,
	     (Real)nflop*volume*nterms/(1e6*dtime*numnodes()) );
#endif
}


int find_backwards_gather( Q_path *path ){
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
int sort_quark_paths( Q_path *src_table, Q_path *dest_table, int npaths ){
    int netdir,dir0,dir1,dir1tmp,thislength,num_new,i,j;

    num_new=0; // number of paths in sorted table
    for( i=0; i<16; i++ ){ // loop over net_back_dirs
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
    if( num_new!=npaths){ node0_printf("OOPS: path table error\n"); exit(0); }
    return 0;
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
*/  /* fermion_force_fn_rhmc.c */
