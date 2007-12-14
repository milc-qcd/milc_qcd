/**************** fermion_force_fn_qdp.c ******************************/
/* MIMD version 7 */
/* General force for any FN-type action.  Compile with
 * fermion_force_general.c to cover all cases of one, two, and N terms.

 * Optimized to transport only one set of SU(3) matrices.
 * D.T. 12/05 Version  3. created for improved fermion RHMC.
 * C.D. 10/06 Converted to QDP (original code in fermion_force_fn_multi.c)
 *      Added multi vector version.	
 */

/* External entry points

   fermion_force_fn_multi

 */

#include "generic_ks_includes.h"	/* definitions files and prototypes */
#include "../include/generic_qdp.h"
#include "../include/generic_ks_qdp.h"
#include <lattice_qdp.h>
#include <string.h>

static int net_back_dirs[16] = { XDOWN, YDOWN, ZDOWN, TDOWN, 
				 XUP, YUP, ZUP, TUP, 
				 X3DOWN, Y3DOWN, Z3DOWN, T3DOWN, 
				 X3UP, Y3UP, Z3UP, T3UP };

#ifdef QCDOC
#define special_alloc qcdoc_alloc
#define special_free qfree
#else
#define special_alloc malloc
#define special_free free
#endif

#define GOES_FORWARDS(dir) (dir<=TUP)
#define GOES_BACKWARDS(dir) (dir>TUP)


/**********************************************************************/
/*   Utilities                                                        */
/**********************************************************************/

/* Map netdir to the QDP shift dir */

/* Assumes the ordering
   XUP, YUP, ZUP, TUP, TDOWN, ZDOWN, YDOWN, XDOWN,
   X3UP, Y3UP, Z3UP, T3UP, T3DOWN, Z3DOWN, Y3DOWN, X3DOWN
*/

static QDP_Shift fnshift(int netdir){
  if      (netdir <  4) return shiftdirs[netdir];
  else if (netdir <  8) return shiftdirs[7-netdir];
  else if (netdir < 12) return shiftdirs[netdir-4];
  else                  return shiftdirs[19-netdir];
}

static QDP_ShiftDir fndir(int netdir){
  if      (netdir <  4) return QDP_forward;
  else if (netdir <  8) return QDP_backward;
  else if (netdir < 12) return QDP_forward;
  else                  return QDP_backward;
}

/* special case to transport a "connection" by one link, does both parities */
static void 
link_transport_connection_qdp( QDP_ColorMatrix *dest, QDP_ColorMatrix *src,
			       QDP_ColorMatrix *gf[4], QDP_ColorMatrix *work,
                               QDP_ColorMatrix *st[8], int dir ){
  if( GOES_FORWARDS(dir) ) {
    QDP_M_eq_M(work, src, QDP_all);
    QDP_M_eq_sM(st[dir], work, QDP_neighbor[dir], QDP_forward, QDP_all);
    QDP_M_eq_M_times_M(dest, gf[dir], st[dir], QDP_all);
    QDP_discard_M(st[dir]);
  }
  else { /* GOES_BACKWARDS(dir) */
    QDP_M_eq_Ma_times_M(work, gf[OPP_DIR(dir)], src, QDP_all);
    QDP_M_eq_sM(st[dir], work, QDP_neighbor[OPP_DIR(dir)], 
		QDP_backward,QDP_all);
    QDP_M_eq_M(dest, st[dir], QDP_all);
    QDP_discard_M(st[dir]);
  }
} /* link_transport_connection_qdp */

static int 
find_backwards_gather( Q_path *path ){
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


// Make a new path table.  Sorted principally by total displacement of
// path.  Below that, sort by direction of first link 
// Below that, sort by direction of second link - note special case of
// one link paths
static int 
sort_quark_paths( Q_path *src_table, Q_path *dest_table, int npaths ){
  int netdir,dir0,dir1,dir1tmp,thislength,num_new,i,j;
  
  num_new=0; // number of paths in sorted table
  for( i=0; i<16; i++ ){ // loop over net_back_dirs
    netdir = net_back_dirs[i]; // table of possible displacements for Fat-Naik
    for( dir0=0; dir0<=7; dir0++){ // XUP ... TDOWN
      for( dir1=-1; dir1<=7; dir1++){ // NODIR, XUP ... TDOWN
	if( dir1==-1 )dir1tmp=NODIR; else dir1tmp=dir1;
	for( j=0; j<npaths; j++ ){
	  // pick out paths with right net displacement
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


/**********************************************************************/
/*   General QDP FN Version for "nterms" sources                      */
/**********************************************************************/
static void 
fn_fermion_force_multi_qdp( QDP_ColorMatrix *force[], QDP_ColorMatrix *gf[], 
			    Real eps, QLA_Real *res, 
			    QDP_ColorVector *x[], int nterms,
			    ferm_links_t *fn, ks_action_paths *ap ){

  /* note CG_solution and Dslash * solution are combined in "x" */
  /* New version 1/21/99.  Use forward part of Dslash to get force */
  /* see long comment at end */
  /* For each link we need x transported from both ends of path. */
  int term;
  int i,j,k,lastdir=-99,ipath,ilink;
  int length,dir,odir;
  QDP_ColorMatrix *tmat;
  Real ferm_epsilon, coeff;
  int num_q_paths = ap->num_q_paths;
  Q_path *q_paths = ap->q_paths;
  Q_path *this_path;	// pointer to current path
  Q_path *last_path;	// pointer to previous path
  QDP_ColorMatrix *mat_tmp0, *stmp[8];
  QDP_ColorMatrix *oprod_along_path[MAX_PATH_LENGTH+1];
  QDP_ColorMatrix *mats_along_path[MAX_PATH_LENGTH+1];
  QDP_ColorMatrix *force_accum[4];  // accumulate force
  QDP_ColorVector *vec_tmp[2];
  int netbackdir, last_netbackdir;	// backwards direction for entire path
  //int tempflops = 0; //TEMP
  // Quark paths sorted by net displacement and last directions
  static Q_path *q_paths_sorted = NULL; 
  // table of net path displacements (backwards from usual convention)
  static int *netbackdir_table = NULL;
  static int first_force = 1;  // 1 if force hasn't been called yet
  char myname[] = "fn_fermion_force_multi_qdp";
#ifdef FFTIME
  int nflop = 966456 + 1440*nterms; // Asqtad action 10/5/06 version of code;
#endif
  double dtime;

  /* node0_printf("STARTING fn_fermion_force_multi_qdp() nterms = %d\n",nterms);*/
  if( nterms==0 )return;

  dtime=-dclock();
	
  /* Allocate fields */
  for(i=0;i<=MAX_PATH_LENGTH;i++){
     oprod_along_path[i] = QDP_create_M();
  }
  for(i=1;i<=MAX_PATH_LENGTH;i++){ 
    // 0 element is never used (it's unit matrix)
     mats_along_path[i] = QDP_create_M();
  }
  for(i=XUP;i<=TUP;i++){
     force_accum[i] = QDP_create_M();
  }
  mat_tmp0   = QDP_create_M();
  for(i=0; i<8; i++) stmp[i] = QDP_create_M();
  tmat       = QDP_create_M();
  vec_tmp[0] = QDP_create_V();
  vec_tmp[1] = QDP_create_V();
  if( vec_tmp[1] == NULL ){
    printf("%s(%d) NO ROOM\n",myname,this_node); 
    terminate(1);
  }

  ferm_epsilon = 2.0*eps; // we only do forward paths, gives factor of 2

  /* Sort the paths */
  if( first_force==1 ){
    if( q_paths_sorted==NULL ) 
      q_paths_sorted = (Q_path *)malloc( num_q_paths*sizeof(Q_path) );
    if(netbackdir_table==NULL ) 
      netbackdir_table = (int *)malloc( num_q_paths*sizeof(int) );
    else{ node0_printf("WARNING: remaking sorted path table\n"); exit(0); }
    sort_quark_paths( q_paths, q_paths_sorted, num_q_paths );
    for( ipath=0; ipath<num_q_paths; ipath++ )
      netbackdir_table[ipath] = 
	find_backwards_gather( &(q_paths_sorted[ipath]) );
    first_force=0;
  }
  
  // clear force accumulators
  for(dir=XUP;dir<=TUP;dir++)
    QDP_M_eq_zero(force_accum[dir], QDP_all);

  /* loop over paths, and loop over links in path */
  last_netbackdir = NODIR;
  last_path = NULL;
  for( ipath=0; ipath<num_q_paths; ipath++ ){
    this_path = &(q_paths_sorted[ipath]);
    if(this_path->forwback== -1)continue;	/* skip backwards dslash */
    
    length = this_path->length;
    // find gather to bring x[term] from "this site" to end of path
    //netbackdir = find_backwards_gather( &(q_paths_sorted[ipath]) );
    netbackdir = netbackdir_table[ipath];
    // and bring x to end - no gauge transformation !!
    // The resulting outer product matrix has gauge transformation
    // properties of a connection from start to end of path
    if( netbackdir != last_netbackdir){ 
      // don't need to repeat this if same net disp. as last path
      k=0; // which vec_tmp we are using (0 or 1)
      QDP_V_eq_sV(vec_tmp[k], x[0], 
	  fnshift(netbackdir), fndir(netbackdir), QDP_all);
      // actually last site in path
      QDP_M_eq_zero(oprod_along_path[0], QDP_all);
      
      for(term=0;term<nterms;term++){
	if(term<nterms-1){
	  QDP_V_eq_sV(vec_tmp[1-k], x[term+1], 
		      fnshift(netbackdir), fndir(netbackdir), QDP_all);
	}
	QDP_M_eq_V_times_Va(tmat, x[term], vec_tmp[k], QDP_all);
	QDP_discard_V(vec_tmp[k]);
	QDP_M_peq_r_times_M(oprod_along_path[0], &res[term], tmat, QDP_all);
	k=1-k; // swap 0 and 1
      } /* end loop over terms in rational function expansion */
    }
    //tempflops+=54*nterms;
    //tempflops+=36*nterms;

    /* path transport the outer product, or projection matrix, of x[term]
       (EVEN sites)  and Dslash*x[term] (ODD sites) from far end.
       
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
      while ( j>0 && 
	      this_path->dir[j] == 
	      last_path->dir[j+last_path->length-length] ) j--;
    if( GOES_BACKWARDS(this_path->dir[0]) ) k=1; else k=0;
    
    for(ilink=j;ilink>=k;ilink--){
      link_transport_connection_qdp( oprod_along_path[length-ilink], 
				     oprod_along_path[length-ilink-1], gf,
                                     mat_tmp0, stmp, this_path->dir[ilink] );
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
	  QDP_M_eq_sM(tmat, gf[dir], QDP_neighbor[dir],
		      QDP_backward, QDP_all);
	  QDP_M_eq_Ma(mats_along_path[1], tmat, QDP_all);
	  QDP_discard_M(tmat);
	}
	else{
	  QDP_M_eq_M(mats_along_path[1], gf[OPP_DIR(dir)], QDP_all); 
	}
      }
      else { // ilink != 0
        dir = OPP_DIR(this_path->dir[ilink]);
        link_transport_connection_qdp( mats_along_path[ilink+1], 
				       mats_along_path[ilink], gf,
				       mat_tmp0, stmp, dir );
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
      coeff = ferm_epsilon*this_path->coeff;
      if( (ilink%2)==1 )coeff = -coeff;
      
      if(ilink==0 && GOES_FORWARDS(dir) ) {
	QDP_M_eq_M(mat_tmp0, oprod_along_path[length],QDP_all); 
      } 
      else if( ilink>0) {
	QDP_M_eq_M_times_Ma(mat_tmp0, oprod_along_path[length-ilink],  
			    mats_along_path[ilink], QDP_all);
      }
      //if(ilink>0)tempflops+=9*22;

      /* add in contribution to the force */
      /* Put antihermitian traceless part into momentum */
      if( ilink<length && GOES_FORWARDS(dir) ){
	QDP_M_peq_r_times_M(force_accum[dir], &coeff, mat_tmp0, QDP_even);
	QDP_M_meq_r_times_M(force_accum[dir], &coeff, mat_tmp0, QDP_odd);
	//tempflops+=36;
      }
      if( ilink>0 && GOES_BACKWARDS(lastdir) ){
	odir = OPP_DIR(lastdir);
	QDP_M_meq_r_times_M(force_accum[odir], &coeff, mat_tmp0, QDP_even);
	QDP_M_peq_r_times_M(force_accum[odir], &coeff, mat_tmp0, QDP_odd);
      }
      //tempflops+=36;
      
      lastdir = dir;
    } /* end loop over links in path */
    last_netbackdir = netbackdir;
    last_path = &(q_paths_sorted[ipath]);
  } /* end loop over paths */

  // add force to momentum
  for(dir=XUP; dir<=TUP; dir++){
    QDP_M_eq_antiherm_M(mat_tmp0, force_accum[dir], QDP_all);
    QDP_M_peq_M(force[dir], mat_tmp0, QDP_all);
  }
  //tempflops+=4*18;
  //tempflops+=4*18;
  
  QDP_destroy_V( vec_tmp[0] );
  QDP_destroy_V( vec_tmp[1] );
  QDP_destroy_M( mat_tmp0 );
  for(i=0; i<8; i++) QDP_destroy_M(stmp[i]);
  QDP_destroy_M( tmat );
  for(i=0;i<=MAX_PATH_LENGTH;i++){
     QDP_destroy_M( oprod_along_path[i] );
  }
  for(i=1;i<=MAX_PATH_LENGTH;i++){
     QDP_destroy_M( mats_along_path[i] );
  }
  for(i=XUP;i<=TUP;i++){
     QDP_destroy_M( force_accum[i] );
  }
  dtime += dclock();
#ifdef FFTIME
  node0_printf("FFTIME:  time = %e (FNMAT) terms = %d mflops = %e\n",dtime,nterms,
	       (Real)nflop*volume/(1e6*dtime*numnodes()) );
#endif
  //printf("FF flops = %d\n",tempflops);
}


/**********************************************************************/
/*   General FN Version for "nterms" sources                          */
/**********************************************************************/
void fermion_force_fn_multi( Real eps, Real residues[], 
			     su3_vector **multi_x, int nterms,
			     int prec, ferm_links_t *fn, ks_action_paths *ap )
{
  /* prec is ignored for now */
  int i,dir;
  site *s;
  su3_matrix tmat;
  QLA_ColorMatrix *temp;
  QDP_ColorMatrix *gf[4];
  QDP_ColorMatrix *force[4];
  QDP_ColorVector **x;
  QLA_Real *res;
  double remaptime = -dclock();

  /* Allocate space for gauge field and force */
  FORALLUPDIR(dir){
    gf[dir] = QDP_create_M();
    force[dir]  = QDP_create_M();
  }
  
  /* Allocate space for multisource vector */
  x = (QDP_ColorVector **)malloc(nterms*sizeof(QDP_ColorVector *));
  for(i = 0; i < nterms; i++){
    x[i] = QDP_create_V();
    set_V_from_field(x[i], multi_x[i],EVENANDODD);
  }

  /* Map gauge links to QDP */
  set4_M_from_site(gf, F_OFFSET(link), EVENANDODD);

  /* The force requires a special conversion from the antihermit type */
  FORALLUPDIR(dir){
    temp = QDP_expose_M (force[dir]);
    FORALLSITES(i,s) {
      uncompress_anti_hermitian( &(s->mom[dir]), &tmat );
      memcpy((void *)&temp[i], (void *)&tmat, sizeof(QLA_ColorMatrix));
    }
    QDP_reset_M (force[dir]);
  }

  /* Map the residues */
  res = (QLA_Real *)malloc(nterms*sizeof(QLA_Real));
  for(i = 0; i < nterms; i++)
    res[i] = residues[i];

  /* Evaluate the fermion force */
  remaptime += dclock();
  fn_fermion_force_multi_qdp(force, gf, eps, res, x, nterms, fn, ap);
  remaptime -= dclock();
  
  /* Map the force back to MILC */
  /* It requires a special conversion to the antihermitian type */
  FORALLUPDIR(dir){
    temp = QDP_expose_M (force[dir]);
    FORALLSITES(i,s) {
      memcpy((void *)&tmat, (void *)&temp[i], sizeof(su3_matrix));
      make_anti_hermitian( &tmat, &(s->mom[dir])); 
    }
    QDP_reset_M (force[dir]);
  }
  
  /* Clean up */
  FORALLUPDIR(dir){
    QDP_destroy_M(gf[dir]);
    QDP_destroy_M(force[dir]);
  }
  for(i = 0; i < nterms; i++)
    QDP_destroy_V(x[i]);
  free(res);
  remaptime += dclock();
#ifdef FFTIME
#ifdef REMAP
  node0_printf("FFREMAP:  time = %e\n",remaptime);
#endif
#endif
}

