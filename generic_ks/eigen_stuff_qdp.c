/*************************** eigen_stuff_qdp.c  ***************************/
/* Eigenvalue and Eigevector computation routines.  QDP version.
 * K.O. 8/99 Started.
 * MIMD version 7
 *
 * These routines are for the computation of the Eigenvalues and Eigevectors
 * of the Kogut-Susskind dslash^2.
 */

//#define DEBUG


#include <float.h> /* To control the tolerance of the JACOBI iteration */
#include <lattice_qdp.h>
#include "generic_ks_includes.h"
#include "../include/generic_qdp.h"
#include "../include/generic_ks_qdp.h"
#include "../include/jacobi.h"

#if ( QDP_Precision == 'F' )

#define dslash_qdp_fn_special2 dslash_qdp_F_fn_special2
#define FATLINKS          fatlinks_F
#define LONGLINKS         longlinks_F
#define SETUP_DSLASH setup_dslash_F

#else

#define dslash_qdp_fn_special2 dslash_qdp_D_fn_special2
#define FATLINKS          fatlinks_D
#define LONGLINKS         longlinks_D
#define SETUP_DSLASH setup_dslash_D

#endif

extern QDP_ColorMatrix **FATLINKS, **LONGLINKS;

#define JACOBI_TOL FLT_EPSILON
#define MINITER 5

/*   If you use strict convergence then you stop when g = |M*v - l*v|     *
 *  becomes less than the eigenvalue tolerance. Otherwise you stop        *
 *  when the estimated error of the eigenvalue gets smaller than          *
 *  the eigenvalue tolerance. The later option takes abvantage of         *
 *  the quadratic convergence of the algorithm. (See Kalkreuter's paper   * 
 *  for details). I prefer the strict convergenve which results in some   *
 *  overkill.                                                             */
#define STRICT_CONVERGENCE

#include <qdp.h>
#include <string.h>
#include <qdp_string.h>

void GramSchmidt_qdp(QDP_ColorVector **vector, int Num, QDP_Subset subset);
void normalize_qdp(QDP_ColorVector *vec, QDP_Subset parity);
void project_out_qdp(QDP_ColorVector *vec, QDP_ColorVector **vector,
		 int Num, QDP_Subset parity);
void constructArray_qdp(QDP_ColorVector **eigVec, Matrix *A, QLA_Real *err,
		    QDP_Subset parity);
void RotateBasis_qdp(QDP_ColorVector **eigVec, Matrix *V, QDP_Subset parity);
void mult_spin_pseudoscalar_qdp(field_offset src, field_offset dest);

/* Temporary QDP_ColorVectors used for the squaring */ 
static QDP_ColorVector *temp0=NULL, *temp1[16], *temp2[16];

/* This is the routine that applies the Matrix, whose eigenvalues            *
 * we want to compute, to a vector. For this specific application it is the  *
 * -D_slash^2 of the KS fermions. We only compute on the "subset" sites.     *
 * Where subset can be EVEN, ODD, or ENENANDODD                              */
void
Matrix_Vec_mult_qdp(QDP_ColorVector *src, QDP_ColorVector *res,
		    QDP_Subset subset)
{
  QDP_Subset othersubset;

#ifdef DEBUG
  if(QDP_this_node==0) printf("begin Matrix_Vec_mult_qdp\n");
#endif

  if(subset==QDP_even) othersubset = QDP_odd;
  else if(subset==QDP_odd) othersubset = QDP_even;
  else othersubset = QDP_all;

  load_ferm_links(&fn_links);
  dslash_qdp_fn_special2(src, temp0, othersubset, temp1);
  dslash_qdp_fn_special2(temp0, res, subset, temp2);
  QDP_V_eqm_V(res, res, subset);

#ifdef DEBUG
  if(QDP_this_node==0) printf("end Matrix_Vec_mult_qdp\n");
#endif
}

/* allocates the tags and the temporaries the Matrix_Vec_mult needs */
void
prepare_Matrix()
{
  int i;

  SETUP_DSLASH();
  temp0 = QDP_create_V();
  for(i=0; i<16; i++) {
    temp1[i] = QDP_create_V();
    temp2[i] = QDP_create_V();
  }
}

/* deallocates the tags and the temporaries the Matrix_Vec_mult needs */
void
cleanup_Matrix()
{
  int i;

  QDP_destroy_V(temp0);
  for(i=0; i<16; i++) {
    QDP_destroy_V(temp1[i]);
    QDP_destroy_V(temp2[i]);
  }
}

/* Projects out the *vectors from the vec. Num is the Number of vectors *
 * and subset is the subset on which we work on.                        *
 * The vectors are assumed to be orthonormal.                           */
void
project_out_qdp(QDP_ColorVector *vec, QDP_ColorVector *vector[], int Num,
		QDP_Subset subset)
{
  QLA_Complex cc1;
  //QLA_Complex cc[Num];
  //QDP_ColorVector *vv[Num];
  int i;

  //for(i=0; i<Num; i++) vv[i] = vec;

  //if(Num>0) {
    //QDP_c_veq_V_dot_V(cc, vector, vv, subset, Num);
    //QDP_c_eq_V_dot_V(&cc1, vector[0], vec, subset);
    //node0_printf("%f %f  %f %f\n", QLA_real(cc1), QLA_imag(cc1), QLA_real(cc[0]), QLA_imag(cc[0]));
    //QDP_V_vmeq_c_times_V(vv, cc, vector, subset, Num);
    //}
  for(i=Num-1; i>=0; i--) {
  //for(i=0; i<Num; i++) {
    QDP_c_eq_V_dot_V(&cc1, vector[i], vec, subset);
    QDP_V_meq_c_times_V(vec, &cc1, vector[i], subset);
    //node0_printf("%f %f  %f %f\n", QLA_real(cc1), QLA_imag(cc1), QLA_real(cc[i]), QLA_imag(cc[i]));
  }
}

/* normalizes the vecror vec. Work only on subset. */
void
normalize_qdp(QDP_ColorVector *vec, QDP_Subset subset)
{
  QLA_Real norm;
  QDP_r_eq_norm2_V(&norm, vec, subset);
  norm = 1.0/sqrt(norm);
  QDP_V_eq_r_times_V(vec, &norm, vec, subset);
}

int
Rayleigh_min_qdp(QDP_ColorVector *vec, QDP_ColorVector **eigVec,
		 Real Tolerance,  Real RelTol, int Nvecs, int MaxIter,
		 int Restart, QDP_Subset subset)
{
  QLA_Complex cc;
  QLA_Real beta, cos_theta, sin_theta;
  QLA_Real quot, P_norm, theta, real_vecMp, pMp;
  QLA_Real g_norm, old_g_norm, start_g_norm;
  QDP_ColorVector *Mvec, *grad, *P, *MP;
  int iter;

#ifdef DEBUG
  if(QDP_this_node==0) printf("begin Rayleigh_min_qdp\n");
#endif

  Mvec = QDP_create_V();
  grad = QDP_create_V();
  //oldgrad = QDP_create_V();
  P = QDP_create_V();
  MP = QDP_create_V();

  project_out_qdp(vec, eigVec, Nvecs, subset);
  normalize_qdp(vec, subset);
  Matrix_Vec_mult_qdp(vec, Mvec, subset);
  project_out_qdp(Mvec, eigVec, Nvecs, subset);

  /* Compute the quotient quot=vev*M*vec */
  QDP_r_eq_re_V_dot_V(&quot, vec, Mvec, subset);
  /* quot is real since M is hermitian. quot = vec*M*vec */
#ifdef DEBUG
  if(QDP_this_node==0) printf("Rayleigh_min: Start -- quot=%g\n", quot);
#endif
  /* Compute the grad=M*vec - quot*vec */
  QDP_V_eq_V(grad, Mvec, QDP_all);
  QDP_V_meq_r_times_V(grad, &quot, vec, subset);
  /* set P (the search direction) equal to grad */
  QDP_V_eq_V(P, grad, QDP_all);
  /* compute the norms of P and grad */
  QDP_r_eq_norm2_V(&P_norm, P, subset);
  P_norm = sqrt(P_norm);
  QDP_r_eq_norm2_V(&g_norm, grad, subset);
  g_norm = sqrt(g_norm);
  start_g_norm = g_norm;
  //QDP_V_eq_V(oldgrad, grad, subset);
#ifdef DEBUG
  if(QDP_this_node==0) printf("Rayleigh_min: Start -- g_norm=%g\n", g_norm);
#endif  

  iter = 0;
  while( (g_norm>Tolerance*quot) &&
	 ( ((iter<MaxIter)&&(g_norm/start_g_norm>RelTol)) || (iter<MINITER) )
	 ) {
    iter++;
    Matrix_Vec_mult_qdp(P, MP, subset);
    QDP_r_eq_re_V_dot_V(&real_vecMp, vec, MP, subset);
    QDP_r_eq_re_V_dot_V(&pMp, P, MP, subset);
    theta = 0.5*atan(2.0*real_vecMp/(quot*P_norm - pMp/P_norm));
    sin_theta = sin(theta);
    cos_theta = cos(theta);
    if(sin_theta*cos_theta*real_vecMp>0) {
      theta = theta - 0.5*M_PI;  /* chose the minimum not the maximum */
      sin_theta = sin(theta);  /* the sin,cos calls can be avoided */
      cos_theta = cos(theta);
    }
    sin_theta = sin_theta/P_norm;
    /* vec = cos(theta)*vec +sin(theta)*P/p_norm */
    //dax_p_by_qdp(cos_theta, vec, sin_theta, P, subset);
    QDP_V_eq_r_times_V(vec, &cos_theta, vec, subset);
    QDP_V_peq_r_times_V(vec, &sin_theta, P, subset);
    /* Mvec = cos(theta)*Mvec +sin(theta)*MP/p_norm */
    //dax_p_by_qdp(cos_theta, Mvec, sin_theta, MP, subset);
    QDP_V_eq_r_times_V(Mvec, &cos_theta, Mvec, subset);
    QDP_V_peq_r_times_V(Mvec, &sin_theta, MP, subset);
    /* renormalize vec ... */
    if( iter%Restart == 0 ) {
#ifdef DEBUG
      {
	QLA_Real vec_norm;
	if(QDP_this_node==0) printf("Renormalizing...");
	QDP_r_eq_norm2_V(&vec_norm, vec, subset);
	if(QDP_this_node==0) printf("  norm: %g\n", sqrt(vec_norm));
      }
#endif
      /* Project vec on the orthogonal complement of eigVec */
      project_out_qdp(vec, eigVec, Nvecs, subset);
      normalize_qdp(vec, subset);
      Matrix_Vec_mult_qdp(vec, Mvec, subset);
      /* Recompute the quotient */
      QDP_r_eq_re_V_dot_V(&quot, vec, Mvec, subset);
      /* Recompute the grad */
      QDP_V_eq_V(grad, Mvec, QDP_all);
      QDP_V_meq_r_times_V(grad, &quot, vec, subset);
      //QDP_r_eq_norm2_V(&g_norm, grad, subset);
      //printf("g_norm = %g\n", g_norm);
      /* Project P on the orthogonal complement of eigVec */
      //QDP_r_eq_norm2_V(&P_norm, P, subset);
      //printf("P_norm = %g\n", P_norm);
      project_out_qdp(P, eigVec, Nvecs, subset);
      //QDP_r_eq_norm2_V(&P_norm, P, subset);
      //printf("P_norm = %g\n", P_norm);
      /* make P orthogonal to vec */
      QDP_c_eq_V_dot_V(&cc, vec, P, subset);
      //printf("cc = %g\n", QLA_real(cc));
      QDP_V_meq_c_times_V(P, &cc, vec, subset);
      //QDP_r_eq_norm2_V(&P_norm, P, subset);
      //printf("P_norm = %g\n", P_norm);
      /* make P orthogonal to grad */
      QDP_c_eq_V_dot_V(&cc, grad, P, subset);
      //printf("cc = %g\n", QLA_real(cc));
      QDP_V_meq_c_times_V(P, &cc, grad, subset);
      QDP_r_eq_norm2_V(&P_norm, P, subset);
      P_norm = sqrt(P_norm);
    }
    QDP_r_eq_re_V_dot_V(&quot, vec, Mvec, subset);
#ifdef DEBUG
    node0_printf("Rayleigh_min: %i, quot=%8g g=%8g b=%6g P:%6g\n",
		 iter, quot, g_norm, beta, P_norm);
#endif
    old_g_norm = g_norm;

    QDP_V_eq_V(grad, Mvec, QDP_all);
    QDP_V_meq_r_times_V(grad, &quot, vec, subset);

    //QDP_V_meq_V(oldgrad, grad, subset);
    //QDP_r_eq_re_V_dot_V(&g_norm, oldgrad, grad, subset);
    //QDP_V_eq_V(oldgrad, grad, subset);

    QDP_r_eq_norm2_V(&g_norm, grad, subset);
    g_norm = sqrt(g_norm);

    beta = cos_theta*g_norm*g_norm/(old_g_norm*old_g_norm);
    if( beta>2.0 ) beta = 2.0;  /* Cut off beta */

    QDP_c_eq_V_dot_V(&cc, vec, P, subset);
    QLA_real(cc) *= beta;
    QLA_imag(cc) *= beta;
    QDP_V_eq_r_times_V_plus_V(P, &beta, P, grad, subset);
    QDP_V_meq_c_times_V(P, &cc, vec, subset);
    QDP_r_eq_norm2_V(&P_norm, P, subset);
    P_norm = sqrt(P_norm);
  }
  project_out_qdp(vec, eigVec, Nvecs, subset);
  normalize_qdp(vec, subset);
  QDP_destroy_V(MP);
  QDP_destroy_V(P);
  //QDP_destroy_V(oldgrad);
  QDP_destroy_V(grad);
  QDP_destroy_V(Mvec);

  iter++;
#ifdef DEBUG
  if(QDP_this_node==0) printf("end Rayleigh_min_qdp\n");
#endif
  return iter;
}

/* Returns the projected matrix A and the error of each eigenvector */
void
constructArray_qdp(QDP_ColorVector **eigVec, Matrix *A, QLA_Real *err,
		   QDP_Subset subset)
{
  QLA_Complex cc;
  QLA_Real rr;
  QDP_ColorVector *res, *grad;
  int i, j, Nvecs;

  Nvecs = A->N;
  res = QDP_create_V();
  grad = QDP_create_V();

  for(i=0; i<Nvecs; i++) {
    Matrix_Vec_mult_qdp(eigVec[i], res, subset);
    QDP_r_eq_re_V_dot_V(&rr, res, eigVec[i], subset);
    A->M[i][i].real = rr;
    A->M[i][i].imag = 0.0;
    QDP_V_eq_V(grad, res, subset);
    QDP_V_meq_r_times_V(grad, &rr, eigVec[i], subset);
    QDP_r_eq_norm2_V(&rr, grad, subset);
    err[i] = sqrt(rr);
    for(j=i+1; j<Nvecs; j++) {
      QDP_c_eq_V_dot_V(&cc, res, eigVec[j], subset);
      A->M[i][j].real = QLA_real(cc);
      A->M[i][j].imag = QLA_imag(cc);
      A->M[j][i].real = QLA_real(cc);
      A->M[j][i].imag = -QLA_imag(cc);
    }
  }

  QDP_destroy_V(grad);
  QDP_destroy_V(res);
}

void
RotateBasis_qdp(QDP_ColorVector **eigVec, Matrix *V, QDP_Subset subset)
{
  QLA_Complex z;
  QDP_ColorVector **Tmp;
  int i, j, N;

  N = V->N;
  /* Allocate the temporary vectors needed */
  Tmp = malloc(N*sizeof(QDP_ColorVector *));
  for(i=0; i<N; i++) Tmp[i] = QDP_create_V();

  for(i=0; i<N; i++) {
    QDP_V_eq_zero(Tmp[i], subset);
    for(j=0; j<N; j++) {
      QLA_real(z) = V->M[j][i].real;
      QLA_imag(z) = V->M[j][i].imag;
      QDP_V_peq_c_times_V(Tmp[i], &z, eigVec[j], subset);
    }
  }

  /* Copy rotated basis to the eigVec and free temporaries */
  for(i=0; i<N; i++) {
    QDP_V_eq_V(eigVec[i], Tmp[i], subset);
    normalize_qdp(eigVec[i], subset);
    QDP_destroy_V(Tmp[i]);
  }
  free(Tmp);
}

int
Kalkreuter_qdp(QDP_ColorVector **eigVec, double *eigVal, Real Tolerance, 
	       Real RelTol, int Nvecs, int MaxIter, int Restart, int Kiters,
	       QDP_Subset subset)
{
  QLA_Real max_error = 1.0e+10;
  QLA_Real min_grad;
  QLA_Real *grad, *err;
  Matrix Array, V;
  QDP_ColorVector *vec;
  int total_iters=0;
  int i, j;
  int iter = 0;

#ifdef DEBUG
  if(QDP_this_node==0) printf("begin Kalkreuter_qdp\n");
#endif

  prepare_Matrix();

  Array = AllocateMatrix(Nvecs);  /* Allocate the array */
  V = AllocateMatrix(Nvecs);      /* Allocate the Eigenvector matrix */

  vec = QDP_create_V();
  grad = malloc(Nvecs*sizeof(QLA_Real));
  err = malloc(Nvecs*sizeof(QLA_Real));

  /* Initiallize all the eigenvectors to a random vector */
  for(j=0; j<Nvecs; j++) {
    grad[j] = 1.0e+10;
    QDP_V_eq_gaussian_S(eigVec[j], rand_state, QDP_all);
    eigVal[j] = 1.0e+16;
    //project_out_qdp(eigVec[j], eigVec, j, subset);
    //normalize_qdp(eigVec[j], subset);
  }

#if 0
  constructArray_qdp(eigVec, &Array, grad, subset);
  Jacobi(&Array, &V, JACOBI_TOL);
  sort_eigenvectors(&Array, &V);
  RotateBasis_qdp(eigVec, &V, subset);
#endif

  while( (max_error>Tolerance) && (iter<Kiters) ) {
    iter++;

    min_grad = grad[0]/eigVal[0];
    for(i=1; i<Nvecs; i++) {
      if(grad[i]<min_grad*eigVal[i]) min_grad = grad[i]/eigVal[i];
    }

    RelTol = 0.3;
    for(j=0; j<Nvecs; j++) {
      if(grad[j]>Tolerance*eigVal[j]) {
	QLA_Real rt;
	rt = RelTol*min_grad*eigVal[j]/grad[j];
	//rt = 1e-5/grad[j];
	if(rt>RelTol) rt = RelTol;
	//rt = RelTol;
	QDP_V_eq_V(vec, eigVec[j], QDP_all);
	total_iters += Rayleigh_min_qdp(vec, eigVec, Tolerance, rt,
					j, MaxIter, Restart, subset);
	QDP_V_eq_V(eigVec[j], vec, QDP_all);
      }
    }
    constructArray_qdp(eigVec, &Array, grad, subset);

    for(i=0; i<Nvecs; i++)
      node0_printf("quot(%i) = %g +/- %8e |grad|=%g\n",
		   i, Array.M[i][i].real, err[i], grad[i]);

#ifdef DEBUG
    node0_printf("Eigenvalues before diagonalization\n");
    for(i=0;i<Nvecs;i++)
      node0_printf("quot(%i) = %g |grad|=%g\n",i,Array.M[i][i].real,grad[i]);
#endif

    Jacobi(&Array, &V, JACOBI_TOL);
    sort_eigenvectors(&Array, &V);
    RotateBasis_qdp(eigVec, &V, subset);
    constructArray_qdp(eigVec, &Array, grad, subset);

    /* find the maximum error */
    max_error = 0.0;
    for(i=0; i<Nvecs; i++) {
      err[i] = eigVal[i];
      eigVal[i] = Array.M[i][i].real;
      err[i] = fabs(err[i] - eigVal[i])/(1.0 - RelTol*RelTol);
      if(eigVal[i]>1e-10) {
#ifndef STRICT_CONVERGENCE
	if(err[i]/eigVal[i]>max_error) max_error = err[i]/eigVal[i];
#else
	if(grad[i]/eigVal[i]>max_error) max_error = grad[i]/eigVal[i];
#endif
      }
    }

    node0_printf("\nEigenvalues after diagonalization at iteration %i\n",iter);
    for(i=0; i<Nvecs; i++)
      node0_printf("quot(%i) = %g +/- %8e |grad|=%g\n",
		   i, eigVal[i], err[i], grad[i]);
  }

  node0_printf("BEGIN RESULTS\n");
  for(i=0;i<Nvecs;i++){
    node0_printf("Eigenvalue(%i) = %g +/- %8e\n",
		 i,eigVal[i],err[i]);
  }

#if 0
  node0_printf("BEGIN EIGENVALUES\n");
  for(i=0; i<Nvecs; i++) {
    double ev, er;
    ev = sqrt(eigVal[i]);
    er = err[i]/(2*ev);
    node0_printf("%.8g\t%g\n", ev, er);
  }
  node0_printf("END EIGENVALUES\n");

  {
    QDP_Writer *qw;
    QDP_String *md;
    char evstring[100], *fn="eigenvecs.out";
    md = QDP_string_create();

    sprintf(evstring, "%i", Nvecs);
    QDP_string_set(md, evstring);
    qw = QDP_open_write(md, fn, QDP_SINGLEFILE);

    for(i=0; i<Nvecs; i++) {
      double ev, er;
      ev = sqrt(eigVal[i]);
      er = err[i]/(2*ev);
      sprintf(evstring, "%.8g\t%g", ev, er);
      QDP_string_set(md, evstring);
      QDP_write_V(qw, md, eigVec[i]);
    }
    QDP_close_write(qw);
    QDP_string_destroy(md);
  }

#endif

  /** Deallocate the arrays **/
  deAllocate(&V) ;
  deAllocate(&Array) ;
  free(err);
  free(grad);
  QDP_destroy_V(vec);
  cleanup_Matrix();
#ifdef DEBUG
  if(QDP_this_node==0) printf("end Kalkreuter_qdp\n");
#endif
  return total_iters;
}

/*
 *  Wrapper for the MILC API
 */

int
Kalkreuter(su3_vector **eigVec, double *eigVal, Real Tolerance,
	   Real RelTol, int Nvecs, int MaxIter,
	   int Restart, int Kiters, int parity,
	   ferm_links_t *fn)
{
  //QLA_Real *ev;
  QDP_ColorVector **vec;
  QDP_Subset subset;
  int i, its;
  su3_matrix *t_fatlink;
  su3_matrix *t_longlink;

#ifdef DEBUG
  if(QDP_this_node==0) printf("begin Kalkreuter\n");
#endif

  if(parity==EVEN) subset = QDP_even;
  else if(parity==ODD) subset = QDP_odd;
  else subset = QDP_all;

  t_longlink = fn->fl.lng;
  t_fatlink = fn->fl.fat;

  set4_M_from_field(FATLINKS, t_fatlink, EVENANDODD);
  set4_M_from_field(LONGLINKS, t_longlink, EVENANDODD);

  //ev = malloc(Nvecs*sizeof(QLA_Real));
  vec = malloc(Nvecs*sizeof(QDP_ColorVector *));
  for(i=0; i<Nvecs; i++) {
    vec[i] = QDP_create_V();
    //ev[i] = eigVal[i];
    set_V_from_field(vec[i], eigVec[i],EVENANDODD);
  }

  its = Kalkreuter_qdp(vec, eigVal, Tolerance, RelTol, Nvecs, MaxIter, Restart,
		       Kiters, subset);

  for(i=0; i<Nvecs; i++) {
    //eigVal[i] = ev[i];
    set_field_from_V(eigVec[i], vec[i],EVENANDODD);
    QDP_destroy_V(vec[i]);
  }
  free(vec);
  //free(ev);
#ifdef DEBUG
  if(QDP_this_node==0) printf("end Kalkreuter\n");
#endif
  return its;
}

/* measures the chirality of a normalized fermion state */
/* DOESN'T APPEAR TO BE USED - CD */

void
measure_chirality_qdp(QDP_ColorVector *src, double *chirality,
		      QDP_Subset subset)
{
  complex tmp;
  double cc;
  int i, parity;
  site *s;

  if(subset==QDP_even) parity = EVEN;
  else if(subset==QDP_odd) parity = ODD;
  else parity = EVENANDODD;

  set_site_from_V(F_OFFSET(tempvec[3]), src,EVENANDODD);
  mult_spin_pseudoscalar(F_OFFSET(tempvec[3]),F_OFFSET(ttt));

  cc = 0.0;
  FORSOMEPARITY(i,s,parity) {
    tmp = su3_dot( &(s->tempvec[3]), &(s->ttt) );
    cc += tmp.real;  /* chirality is real since Gamma_5 is hermitian */
  }
  g_doublesum(&cc);
  *chirality = cc;
}

/* prints the density and chiral density of a normalized fermion state */
void
print_densities_qdp(QDP_ColorVector *src, char *tag, int y, int z, int t,
		    QDP_Subset subset)
{
  complex tmp1, tmp2;
  int i, parity;
  site *s;

  if(subset==QDP_even) parity = EVEN;
  else if(subset==QDP_odd) parity = ODD;
  else parity = EVENANDODD;

  set_site_from_V(F_OFFSET(tempvec[3]), src,EVENANDODD);

  mult_spin_pseudoscalar(F_OFFSET(tempvec[3]),F_OFFSET(ttt));

  FORSOMEPARITY(i,s,parity) {
    if( (s->y==y)&&(s->z==z)&&(s->t==t) ) {
      tmp1 = su3_dot( &(s->tempvec[3]), &(s->ttt) ) ;
      tmp2 = su3_dot( &(s->tempvec[3]), &(s->tempvec[3]) ) ;
      node0_printf("%s: %i %e %e %e\n", tag, s->x, tmp2.real,
		   tmp1.real, tmp1.imag);
    }
  }
}

/* measures the chiraliry of a normalized fermion state */
void measure_chirality(su3_vector *src, double *chirality, int parity)
{
  register int i;
  register site *s;
  register double cc ;
  complex tmp ;
 
  FORSOMEPARITY(i,s,parity){
    su3vec_copy(&src[i],&(s->tempvec[3])) ;
  }
 
  mult_spin_pseudoscalar(F_OFFSET(tempvec[3]),F_OFFSET(ttt)) ;
 
  cc = 0.0 ;
  FORSOMEPARITY(i,s,parity){
    tmp = su3_dot( &(s->tempvec[3]), &(s->ttt) ) ;
    cc +=  tmp.real ; /* chirality is real since Gamma_5 is hermitian */
  }
  *chirality = cc ;
  g_doublesum(chirality);

} /* chirality.c */
