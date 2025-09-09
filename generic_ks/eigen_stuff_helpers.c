
/****** eigen_stuff_linalg.c  ******************/
/* Linear algebra for eigenvalue and eigevector computation routines.
   These routines work with color vectors

* K.O. 8/99 Started. 
* UH did some stuff to this I think
* EBG 2/15/02 changed phi -> phi1, mass -> mass1 so could be used 
*    w/ ks_imp_dyn2
* EBG 6/2004 fixed some memory leaks
* MIMD version 7
*/

/* This routine is the routine that applies the Matrix, whose eigenvalues    *
 * we want to compute, to a vector. For this specific application it is the  *
 * -D_slash^2 of the KS fermions. We only compute on the "parity" sites.     *
 * Where parity can be EVEN, ODD, or ENENANDODD                              */

#include "generic_ks_includes.h"
#include "../include/dslash_ks_redefine.h"
#include "../include/flavor_ops.h"
#include "../include/openmp_defs.h"

/************************************************************************/
/* Temporary su3_vectors used for the squaring */ 
static su3_vector *temp = NULL ;
/* flag indicating wether to start dslash. */
static int dslash_start = 1 ; /* 1 means start dslash */
/* Message tags to be used for the Matrix Vector multiplication */
static msg_tag *tags1[16],*tags2[16];
/* Other parity */
static int otherparity;

/************************************************************************/
static int 
opposite_parity(int parity){
  
  otherparity = EVENANDODD;
  
  switch(parity){
  case EVEN:
    otherparity = ODD ;
    break ;
  case ODD:
    otherparity = EVEN ;
    break ;
  case EVENANDODD:
    otherparity = EVENANDODD ;
    break ;
  default:
    node0_printf("ERROR: wrong parity in eigen_stuff::opposite_parity\n") ;
    terminate(1) ;
  }
  
  return otherparity;
} 

/************************************************************************/
static void 
DdagD( su3_vector *res, su3_vector *src, int ds, msg_tag *t1[], msg_tag *t2[], 
       ks_eigen_param *eigen_param, imp_ferm_links_t *fn )
{
  int i;
  int parity = eigen_param->parity;
  otherparity = opposite_parity(parity);
  
  dslash_fn_field_special(src , temp, otherparity, t1, ds, fn) ;
  dslash_fn_field_special(temp, res , parity     , t2, ds, fn) ;
  
  FORSOMEFIELDPARITY_OMP(i,parity, ){ 
    scalar_mult_su3_vector( &(res[i]), -1.0, &(res[i])) ;
  } END_LOOP_OMP;
}

/************************************************************************/
/* Chebyshev polynomial                                                 */

static double 
poly( double am, double aM, int p, double x) {
  double delta = 0.5*(aM-am);
  double theta = 0.5*(aM+am);
  
  double d1    = -1.0/theta;
  double sig0  = d1*delta;
  double d2    = 1.0;
  double d3    = 0.0;
  
  double x0    = 1.0;
  double x1    = d1*x + d2;
  double x2    = x1;
  
  double sig1  = sig0;
  double sigma = 0.0;
  
  for (int i=2; i<=p; i++) {
    sigma = 1.0/(2.0/sig0 - sig1);
    
    d1 = 2.0*sigma/delta;
    d2 = -d1*theta;
    d3 = -sigma*sig1;
    
    x2 = d3*x0 + d2*x1 + d1*x1*x;
    
    x0 = x1;
    x1 = x2;
    
    sig1 = sigma;
  }
  
  return x2;
}

/************************************************************************/
/* Matrix vector operators                                              */

/* Use Chebyshev preconditioned operator p(DdagD) vec */

#if defined(HAVE_PRIMME) || defined(HAVE_ARPACK)
static void 
PDdagD( su3_vector *res, su3_vector *src, ks_eigen_param *eigen_param, imp_ferm_links_t *fn )

{
  /* Chebyshev operator */

  int parity = eigen_param->parity;
  double am = eigen_param->poly.minE;   /* Lower focus of window of exclusion */
  double aM = eigen_param->poly.maxE;   /* Upper focus of window of exclusion */
  int p = eigen_param->poly.norder;     /* Degree. Must be >= 1 */
  
  double delta = 0.5*(aM-am);
  double theta = 0.5*(aM+am);
  
  double d1    = -1.0/theta;
  double sig0  = d1*delta;
  double d2    = 1.0;
  double d3    = 0.0;
  msg_tag *tags3[16],*tags4[16];
  /* double x0    = 1; */
  su3_vector *x0 = create_v_field();
  copy_v_field(x0, src);

  /* double x1    = d1*x + d2; */
  su3_vector *x1 = create_v_field();
  su3_vector *y1 = create_v_field();
  DdagD(y1, src, 1, tags3, tags4, eigen_param, fn);
  saxpby_v_field( x1, d1, y1, d2, x0);
  cleanup_gathers(tags3,tags4);

  /* double x2    = x1; */
  su3_vector *x2 = create_v_field();
  copy_v_field(x2, x1);
  
  double sig1  = sig0;
  double sigma = 0.0;

  int dslash_start2 = 1;
  for (int i=2; i<=p; i++) {
    sigma = 1.0/(2.0/sig0 - sig1);
    
    d1 = 2.0*sigma/delta;
    d2 = -d1*theta;
    d3 = -sigma*sig1;
    
    /* x2 = theta*(d3*x0 + d2*x1 + d1*x*x1); */
    DdagD(y1, x1, dslash_start2, tags3, tags4, eigen_param, fn);
    dslash_start2 = 0; /* So we do restart_gathers when we loop again */
    saxpbypcz_v_field( x2, d3, x0, d2, x1, d1, y1);

    /* x0 = x1; */
    copy_v_field(x0, x1);

    /* x1 = x2; */
    copy_v_field(x1, x2);
    
    sig1 = sigma;
  }

  if(!dslash_start2)
    cleanup_gathers(tags3,tags4);

  /* return x2; */
  copy_v_field(res, x2);

  destroy_v_field(x2);
  destroy_v_field(y1);
  destroy_v_field(x1);
  destroy_v_field(x0);
}
#endif
/*****************************************************************************/
/* The Matrix_Vec_mult and cleanup_Matrix() */

void Matrix_Vec_mult(su3_vector *src, su3_vector *res, 
		     ks_eigen_param *eigen_param,
                     imp_ferm_links_t *fn ){

  /* store last source so that we know when to reinitialize the message tags */
  static su3_vector *last_src=NULL ;

  if(temp == NULL){
    temp = create_v_field();
  }

  /*reinitialize the tags if we have a new source */
  if(last_src != src){
    if(!dslash_start) cleanup_gathers(tags1,tags2);
    dslash_start = 1 ;
    last_src = src ;
  }

#ifdef MATVEC_PRECOND
  PDdagD(res, src, eigen_param, fn);
#else
  DdagD(res, src, dslash_start, tags1, tags2, eigen_param, fn);
  dslash_start = 0 ;  /* Signals that tags1 and tags2 are initialized */
#endif

}

/*****************************************************************************/
#if defined(HAVE_PRIMME) || defined(HAVE_ARPACK)

/* The Matrix_Vec_mult and cleanup_Matrix() */

void Precond_Matrix_Vec_mult(su3_vector *src, su3_vector *res, ks_eigen_param *eigen_param,
			     imp_ferm_links_t *fn ){

  /* store last source so that we know when to reinitialize the message tags */
  static su3_vector *last_src=NULL ;

  if(dslash_start){
    temp = (su3_vector *)malloc(sites_on_node*sizeof(su3_vector));
  }

  /*reinitialize the tags if we have a new source */
  if(last_src != src){
    if(!dslash_start) cleanup_gathers(tags1,tags2);
    dslash_start = 1 ;
    last_src = src ;
  }

  PDdagD(res, src, eigen_param, fn);

  dslash_start = 0 ;
}
#endif

/************************************************************************/
/* Deallocates the tags and the temporaries the Matrix_Vec_mult needs */
void cleanup_Matrix(){
  if(!dslash_start) {
    cleanup_gathers(tags1,tags2); 
  }
  cleanup_dslash_temps() ;
  destroy_v_field(temp) ;
  temp = NULL;
  dslash_start = 1 ;
#ifdef DEBUG
  node0_printf("cleanup_Matrix(): done!\n") ; fflush(stdout) ;
#endif
}

/*****************************************************************************/
/* measures the chiraliry of a normalized fermion state */
void measure_chirality(su3_vector *src, double *chirality, int parity){
  register int i;
  register site *s;
  register double cc ;
  su3_vector *tempvec, *ttt;
  int r0[4] = {0,0,0,0};
  complex tmp ;
  char gam5Xgam1[] = "G1-G5";  /* Note, spin_taste_op inteprets this as G5-G1 */

  ttt = create_v_field();
  tempvec = create_v_field();
  copy_v_field(tempvec, src);
  spin_taste_op(spin_taste_index(gam5Xgam1), r0, ttt, tempvec);

  cc = 0.0 ; 
  FORSOMEPARITY_OMP(i,s,parity,private(tmp) reduction(+:cc) ){
    tmp = su3_dot( tempvec+i, ttt+i );
    cc +=  tmp.real ; /* chirality is real since Gamma_5 is hermitian */
  } END_LOOP_OMP;
  *chirality = cc ;
  g_doublesum(chirality);

  destroy_v_field(tempvec);
  destroy_v_field(ttt);
}


/*****************************************************************************/
/* prints the density and chiral density of a normalized fermion state */
void print_densities(su3_vector *src, char *tag, int y,int z,int t, 
		     int parity){

  register int i;
  register site *s;
  su3_vector *tempvec, *ttt;
  int r0[4] = {0,0,0,0};
  complex tmp1,tmp2 ;

  ttt = create_v_field();
  tempvec = create_v_field();
  copy_v_field(tempvec, src);
  spin_taste_op(spin_taste_index("G5-G1"), r0, ttt, tempvec);

  FORSOMEPARITY(i,s,parity){ 
    if((s->y==y)&&(s->z==z)&&(s->t==t)){
      tmp1 = su3_dot( tempvec+i, ttt+i ) ;
      tmp2 = su3_dot( tempvec+i, tempvec+i ) ;
      node0_printf("%s: %i %e %e %e\n",tag,
		   s->x,tmp2.real,tmp1.real,tmp1.imag);
    }
  } END_LOOP;
}

/*****************************************************************************/
/* Restore the ODD (EVEN) part of KS eigenvectors from the EVEN (ODD) part */
/* Here we don't normalize the result */

void restore_eigVec(int Nvecs, Real *eigVal, su3_vector **eigVec, int parity,
		    imp_ferm_links_t *fn)
{
  char myname[] = "restore_eigVec";
  
  su3_vector *ttt = create_v_field();
  
  for(int i=0; i < Nvecs; i++){
    dslash_fn_field(eigVec[i], ttt, parity, fn);
    int j;
    FORSOMEFIELDPARITY(j, parity){
      complex c = cmplx(0.0,1.0/sqrt(eigVal[i]));
      c_scalar_mult_su3vec(ttt + j, &c, eigVec[i] + j) ;
    } END_LOOP;
  }
  
  destroy_v_field(ttt);
}

/*****************************************************************************/
/* Reset the eigenvalues from the eigenvectors using the Rayleigh quotient   */
/* (Needed when the eigenvalues are from a preconditioned operator)          */

void reset_eigenvalues(su3_vector *eigVec[], Real *eigVal,
		       int Nvecs, int parity, imp_ferm_links_t *fn)
{
  int i,j;
  su3_vector *ttt = create_v_field();

  otherparity = opposite_parity(parity);

  node0_printf("Resetting eigenvalues\n");


  for(i = 0; i < Nvecs; i++){
    double norm = 0., expect = 0.;
    dslash_fn_field(eigVec[i], ttt, otherparity, fn);
    dslash_fn_field(ttt, ttt, parity, fn);
    FORSOMEFIELDPARITY_OMP(j,parity,reduction(+:expect,norm)){
      expect += su3_rdot(eigVec[i]+j, ttt+j); 
      norm += magsq_su3vec(eigVec[i]+j);
    } END_LOOP_OMP;

    /* The expectation value of -DdaggerD on the eigenvector. */
    g_doublesum(&expect);
    g_doublesum(&norm);
    node0_printf("eigVal[%d] = %e --> %e\n", i,eigVal[i], -expect/norm);
    eigVal[i] = -expect/norm;
  }
  destroy_v_field(ttt);
}
  
/*****************************************************************************/
/* Returns the dot product of two fermion vectors */
static void dot_product(su3_vector *vec1, su3_vector *vec2, 
			double_complex *dot, int parity) {
  register double re,im ;
  register  int i;
  
  re=im=0.0;
  FORSOMEFIELDPARITY_OMP(i,parity,reduction(+:re,im)){
    complex cc = su3_dot( &(vec1[i]), &(vec2[i]) );
    re += cc.real ;
    im += cc.imag ;
  } END_LOOP_OMP;
  dot->real = re ; 
  dot->imag = im ;
  g_dcomplexsum(dot);
}

/*****************************************************************************/
/* Returns the dot product of two fermion vectors */
static double dot_product_real(su3_vector *vec1, su3_vector *vec2, int parity) {
  double rdot = 0.;
  register  int i;
  
  FORSOMEFIELDPARITY_OMP(i,parity,reduction(+:rdot)){
    rdot += su3_rdot( &(vec1[i]), &(vec2[i]) );
  } END_LOOP_OMP;
  g_doublesum(&rdot);
  return rdot;
}

/*****************************************************************************/
/* Returns vec2 = vec2 - cc*vec1   cc is a double complex   */
static void complex_vec_mult_sub(double_complex *cc, su3_vector *vec1, 
				 su3_vector *vec2, int parity){
  register  int i;
  complex sc ;
  
  sc.real= (Real)(cc->real) ; 
  sc.imag= (Real)(cc->imag) ;

  FORSOMEFIELDPARITY_OMP(i,parity,){
    c_scalar_mult_sub_su3vec(&(vec2[i]), (&sc), &(vec1[i])) ;
  } END_LOOP_OMP;
}

/*****************************************************************************/
/* Use first-order perturbation theory to shift eigenvalues and eigenvectorsa
   from an old to a new Dslash^2 operator */

void perturb_eigpair(su3_vector *eigVec_new[], Real *eigVal_new,
		     su3_vector *eigVec_old[], Real *eigVal_old,
		     int Nvecs, int parity, imp_ferm_links_t *fn_new,
		     imp_ferm_links_t *fn_old){


  int otherparity = opposite_parity(parity);
  su3_vector *ttt = create_v_field();
  for(int n = 0; n < Nvecs; n++){
    /* ttt <- D_new^2 eigVec[n] */
    dslash_fn_field(eigVec_old[n], ttt, otherparity, fn_new);
    dslash_fn_field(ttt, ttt, parity, fn_new);
    copy_v_field(eigVec_new[n], eigVec_old[n]);
    for(int m = 0; m < Nvecs; m++){
      /* Note that eigVal is -1 times the eigenvalue of D^2 */
      /* eigVal_new = - < eigVec D_new^2 eigVec > */
      if(m == n){
	eigVal_new[m] = -dot_product_real(eigVec_old[m], ttt, parity);
      } else {
	/* deigVec_n = sum(m!=n) eigVec_m <eigVec_m D_new^2 eigVec_n> / (-eigVal_n + eigVal_m) */
	double_complex cc;
	/* cc = <eigVec_m D_new^2 eigVec_n> */
	dot_product(eigVec_old[m], ttt, &cc, parity);
	/* cc = <eigVec_m D_new^2 eigVec_n> / (eigVal_n - eigVal_m) */
	CDIVREAL(cc, eigVal_old[n] - eigVal_old[m], cc);
	complex_vec_mult_sub(&cc, eigVec_old[n], eigVec_new[m], parity);
	//	node0_printf("eig(%d) %d |cc| %f\n",n,m,sqrt(cc.real*cc.real+cc.imag*cc.imag));
      }
    }
  }
  destroy_v_field(ttt);
}
/*****************************************************************************/
/* Check the residues and norms of the eigenvectors                          */
/* Hiroshi Ono */

void check_eigres(double *resid, su3_vector *eigVec[], Real *eigVal,
		  int Nvecs, int parity, imp_ferm_links_t *fn){

  su3_vector *ttt;
  int i,j;
  otherparity = opposite_parity(parity);

  /* compute residual and norm of eigenvectors */
  double *norm = (double *)malloc(Nvecs*sizeof(double));
  ttt = create_v_field();
  for(i = 0; i < Nvecs; i++){
    double rsd = (double)0.0;
    double nrm  = (double)0.0;
    dslash_fn_field(eigVec[i], ttt, otherparity, fn);
    dslash_fn_field(ttt, ttt, parity, fn);
    FORSOMEFIELDPARITY_OMP(j,parity, reduction(+:rsd,nrm)){
      scalar_mult_sum_su3_vector(ttt+j, eigVec[i]+j, eigVal[i]);
      rsd += magsq_su3vec(ttt+j);
      nrm += magsq_su3vec(eigVec[i]+j);
    } END_LOOP_OMP;
    resid[i] = rsd;
    norm[i] = nrm;
  }
  destroy_v_field(ttt);
  g_vecdoublesum(resid, Nvecs);
  g_vecdoublesum(norm, Nvecs);

  /* Fix the normalization if needed */
  double tol = 1e-9;
  int count_renorm = 0;
  int count_zero = 0;
  double maxdev = 0;
  for(i = 0; i < Nvecs; i++){
    if(norm[i] == 0.){
      count_zero++;
    } else {
      double dev = fabs(sqrt(norm[i])-1.);
      maxdev = fmax(maxdev, dev);
      if(dev > tol){
	count_renorm++;
	double fact = 1./sqrt(norm[i]);
	FORSOMEFIELDPARITY_OMP(j,parity,){
	  scalar_mult_su3_vector( eigVec[i]+j, fact, eigVec[i]+j );
	} END_LOOP_OMP;
      }
    }
  }

  if(count_zero > 0){
    node0_printf("check_eigres: %d eigVecs have zero norm\n", count_zero);
    terminate(1);
  }

  node0_printf("Checking eigensolutions\n");
  for(i = 0; i < Nvecs; i++){
    resid[i] = sqrt(resid[i]/norm[i]);
    norm[i] = sqrt(norm[i]);
    node0_printf("eigVal[%d] = %e ( resid = %e , |eigVec[%d]|-1 = %e )\n",
		 i, eigVal[i], resid[i], i, norm[i]-1);
  }
  if(count_renorm > 0){
    node0_printf("check_eigres: NOTE: %d eigVec norms deviated from 1 by more than %g with max deviation %g.\n",
		 count_renorm, tol, maxdev);
    node0_printf("They were renormalized to 1\n");
  }
  
  node0_printf("End of eigensolutions\n");
  free(norm);
}

/**********************************************************************************/
/* Returns the squared norm of a list of n su3_vector fields over the
   lattice parity subset */

static void dvecmagsq(double magsq[], su3_vector *vec[], int parity, int n) {
  int i,j;
  
  for(j = 0; j < n; j++){
    double mgsq = 0.;
    FORSOMEFIELDPARITY_OMP(i,parity,reduction(+:mgsq)){
      mgsq += magsq_su3vec( vec[j]+i );
    } END_LOOP_OMP;
    magsq[j] = mgsq;
  }

  g_vecdoublesum(magsq, n);
}

/*****************************************************************************/
/* Returns vec2 = a*vec1   a is a double vec2 can be vec1*/
static void double_vec_mult(double *a, su3_vector *vec1, 
			    su3_vector *vec2, int parity){
  
  register site *s;
  register  int i;
  
  FORSOMEPARITY_OMP(i,s,parity,){ 
    scalar_mult_su3_vector( &(vec1[i]),(Real)*a, &(vec2[i])) ;
  } END_LOOP_OMP;
}

/*****************************************************************************/
/* Construct odd-site (even-site)  eigenvectors from even (odd)              */
/* (Simply multiply even-site vectors by dslash and normalize                */

void construct_eigen_other_parity(su3_vector *eigVec[], Real eigVal[], 
				  ks_eigen_param *eigen_param, imp_ferm_links_t *fn){
  
  char myname[] = "construct_eigen_other_parity";
  //node0_printf("Entered %s\n", myname);
  int i,j;
  double *magsq;
  int Nvecs = eigen_param->Nvecs;
  int otherparity = 0;
  /* If parity = ODD we assume we have EVEN eigenvectors and vice-versa */
  int parity = eigen_param->parity;
  site *s;

  switch(parity){
  case(EVEN): otherparity=ODD; break;
  case(ODD):  otherparity=EVEN; break;
  }

  for(j = 0; j < Nvecs; j++){
    FORSOMEPARITY_OMP(i,s,otherparity,){
      clearvec(eigVec[j]+i);
    } END_LOOP_OMP;
    dslash_field(eigVec[j], eigVec[j], otherparity, fn);
  }

  /* If we calculate the 2-norms all at once we do only one large
     reduction sum */
  magsq = (double *)malloc(sizeof(double)*Nvecs);
  if(magsq == NULL){
    node0_printf("magsq: no room\n");
    terminate(1);
  }

  dvecmagsq(magsq, eigVec, otherparity, Nvecs);

  /* Normalize the odd-site vectors */
  for(j = 0; j < Nvecs; j++){
    double fact;
    if(magsq[j] != 0.){
      fact = 1./sqrt(magsq[j]);
    } else {
      node0_printf("%s: ERROR. zero norm for eigenvector %d.", myname, j);
      fact = 0.;
    }
    double_vec_mult(&fact, eigVec[j], eigVec[j], otherparity);
  }
  
  free(magsq);
}


