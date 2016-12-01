/******* inc_eigcg.c - Incremental eigCG for SU3/fermions ****/
/* MIMD version 7 */
/* Kogut-Susskind fermions -- this version for "fat plus Naik" quark
   actions.  
*/

/* Hiroshi Ohno */
/* About eigCG algorithm, see A. Stathopoulos and K. Orginos [arXiv:0707.0131] */
#include "generic_ks_includes.h"
#include "../include/blas_lapack.h"
#include "../include/prefetch.h"
#include "../include/openmp_defs.h"
#include <string.h>

#ifdef CGTIME
static const char *prec_label[2] = {"F", "D"};
#endif

#if (PRECISION==1)
#error Requires double precision!
#endif

//#define CG_DEBUG
//#define EIGCG_DEBUG

/* The Fermilab relative residue */
static Real my_relative_residue(su3_vector *p, su3_vector *q, int parity){

  register int i;
  double residue, num, den;
  
  residue = (double)0.0;
  FORSOMEFIELDPARITY_OMP(i, parity, private(num,den) reduction(+:residue)){
    num = (double)magsq_su3vec(p+i);
    den = (double)magsq_su3vec(q+i);
    residue += (den==0) ? 1.0 : (num/den);
  } END_LOOP_OMP;
  g_doublesum(&residue);

  if(parity == EVENANDODD)
    return sqrt(residue/volume);
  else
    return sqrt(2*residue/volume);
}

/* init-CG */
/* dest <- dest + sum_i eigVec[i]*H^{-1}*eigVec[i].resid, 
   resid = src - (-Dslash^2 + 4*mass^2)*dest
*/
static void initCG(su3_vector *src, su3_vector *dest, int Nvecs_curr, int Nvecs_max,
		   su3_vector **eigVec, double_complex *H, Real mass, int parity,
		   imp_ferm_links_t *fn){

  /* Constants */
  int ione = 1;
  int otherparity = (parity == EVEN) ? ODD : EVEN;
  Real msq_x4 = 4.0*mass*mass;
  double dzero = (double)0.0;
  double_complex zzero = dcmplx(dzero, dzero);

  register int i;
  int j, info;
  double_complex cc;
  double_complex *c, *H2;
  su3_vector *resid;

  c = (double_complex *)malloc(Nvecs_curr*sizeof(double_complex));
  resid = (su3_vector *)malloc(sites_on_node*sizeof(su3_vector));

  /* resid <- src - (-Dslash^2 + 4*mass^2)*dest */
  dslash_fn_field(dest, resid, otherparity, fn);
  dslash_fn_field(resid, resid, parity, fn);
  FORSOMEFIELDPARITY_OMP(i, parity, default(shared)){
    scalar_mult_sum_su3_vector(resid+i, dest+i, -msq_x4);
    add_su3_vector(resid+i, src+i, resid+i);
  } END_LOOP_OMP;

  /* c[i] = eigVec[i].resid */
  for(j = 0; j < Nvecs_curr; j++){
//    c[j] = zzero;
//    FORSOMEFIELDPARITY_OMP(i, parity, private(cc) reduction(+:c[j])){
//      cc = su3_dot(eigVec[j]+i, resid+i);
//      CSUM(c[j], cc);
//    } END_LOOP_OMP;

    double cctotr=0., cctoti=0.;
    FORSOMEFIELDPARITY_OMP(i, parity, private(cc) reduction(+:cctotr,cctoti)){
      cc = su3_dot(eigVec[j]+i, resid+i);
      cctotr += cc.real;
      cctoti += cc.imag;
    } END_LOOP_OMP;
    c[j].real = cctotr;
    c[j].imag = cctoti;

  }
  g_vecdcomplexsum(c, Nvecs_curr);

  free(resid);

  H2 = (double_complex *)malloc(Nvecs_curr*Nvecs_curr*sizeof(double_complex));

  /* H2 = H + 4*mass^2*I */
  for(j = 0; j < Nvecs_curr; j++){
    zcopy_(&Nvecs_curr, H+Nvecs_max*j, &ione, H2+Nvecs_curr*j, &ione);
    CSUM(H2[(Nvecs_curr+1)*j], dcmplx(msq_x4, dzero));
  }

  /* Compute H^{-1}*c = H^{-1}*eigVec.resid with Cholesky decomposition */
  zpotrf_("U", &Nvecs_curr, H2, &Nvecs_curr, &info);
  zpotrs_("U", &Nvecs_curr, &ione, H2, &Nvecs_curr, c, &Nvecs_curr, &info);

  free(H2);

  /* dest <- dest + sum_j c[i]*eigVec[i] = dest + sum_i eigVec[i]*H^{-1}*eigVec[i].resid */ 
  for(j = 0; j < Nvecs_curr; j++){
    FORSOMEFIELDPARITY_OMP(i, parity, default(shared)){
      c_scalar_mult_add_su3vec(dest+i, c+j, eigVec[j]+i);
    } END_LOOP_OMP;
  }

  free(c);
}

/* Orthogonalize eigVec by the (modified) Gram-Schmidt method */
/* Assume the first Nvecs_curr eigenvectors have been already orthonormalized.
   If norm of an eigenvector is less than ORTHO_EPS, remove it.
   Rturn the number of new eigenvectors to be added.
*/
static int orthogonalize(int Nvecs, int Nvecs_curr, su3_vector **eigVec, int parity){

  register int i;
  int j, k, Nvecs_add, n;
  double norm;
  double_complex cc;
  double_complex *c;

  j = Nvecs_curr;
  Nvecs_add = Nvecs;
  n = Nvecs_curr + Nvecs_add;

  c = (double_complex *)malloc(n*sizeof(double_complex));

  while(j < n){
    /* Modified Gram-Schmidt
       Orthogonality is better but more communications are needed */
    for(k = 0; k < j; k++){
//      c[k] = dcmplx((double)0.0,(double)0.0);
//      FORSOMEFIELDPARITY_OMP(i, parity, private(cc) reduction(+:c[k])){
//	cc = su3_dot(eigVec[k]+i, eigVec[j]+i);
//	CSUM(c[k], cc);
//      } END_LOOP_OMP;

      double cctotr=0., cctoti=0.;
      FORSOMEFIELDPARITY_OMP(i, parity, private(cc) reduction(+:cctotr,cctoti)){
	cc = su3_dot(eigVec[k]+i, eigVec[j]+i);
	cctotr += cc.real;
	cctoti += cc.imag;
      } END_LOOP_OMP;
      c[k].real = cctotr;
      c[k].imag = cctoti;

      g_dcomplexsum(c+k);
      FORSOMEFIELDPARITY_OMP(i, parity, default(shared)){
	c_scalar_mult_sub_su3vec(eigVec[j]+i, c+k, eigVec[k]+i);
      } END_LOOP_OMP;
    }
    /* Gram-Schmidt
       Less communications but
       poor orthogonality might happen if the number of vectors is too large. */
    /*
    for(k = 0; k < j; k++){
      c[k] = dcmplx((double)0.0,(double)0.0);
      FORSOMEFIELDPARITY_OMP(i, parity, private(cc) reduction(+:c[k])){
	cc = su3_dot(eigVec[k]+i, eigVec[j]+i);
	CSUM(c[k], cc);
      } END_LOOP_OMP;
    }
    g_vecdcomplexsum(c, j);
    for(k = 0; k < j; k++){
      FORSOMEFIELDPARITY_OMP(i, parity, default(shared)){
	c_scalar_mult_sub_su3vec(eigVec[j]+i, c+k, eigVec[k]+i);
      } END_LOOP_OMP;
    }
    */
    norm = (double)0.0;
    FORSOMEFIELDPARITY_OMP(i, parity, reduction(+:norm)){
      norm += magsq_su3vec(eigVec[j]+i);
    } END_LOOP_OMP;
    g_doublesum(&norm);
    norm = sqrt(norm);
    if( norm < ORTHO_EPS ){
      Nvecs_add--;
      n--;
      for(k = j; k < n; k++){
	FORSOMEFIELDPARITY_OMP(i, parity, default(shared)){
	  eigVec[k][i] = eigVec[k+1][i];
	} END_LOOP_OMP;
      }
    }
    else{
      norm = 1.0/norm;
      FORSOMEFIELDPARITY_OMP(i, parity, default(shared)){
	scalar_mult_su3_vector(eigVec[j]+i, norm, eigVec[j]+i);
      } END_LOOP_OMP;
      j++;
    }
  }

  free(c);

  return Nvecs_add;
}

/* Extend the matrix H = -U^+ Dslash^2 U */
/* H_new = | H         W^+U|
           |U^+W  W^+V|,
   where U and V are current and new Ritz vectors,respectively, 
   and W = -Dslash^2*V.
*/
static void extend_H(int Nvecs, int Nvecs_curr, int Nvecs_max, su3_vector **eigVec,
		     double_complex *H, int parity, imp_ferm_links_t *fn){

  /* Constants */
  int otherparity = (parity == EVEN) ? ODD : EVEN;

  register int i;
  int j, k, kk;
  double_complex cc;
  su3_vector *ttt;

  ttt = (su3_vector *)malloc(sites_on_node*sizeof(su3_vector));

  for(j = Nvecs_curr; j < Nvecs_curr+Nvecs; j++){
    /* ttt <- Dslash^2*eigVec[i] */
    dslash_fn_field(eigVec[j], ttt, otherparity, fn);
    dslash_fn_field(ttt, ttt, parity, fn);
    for(k = 0; k < Nvecs_curr+Nvecs; k++){
      /* H_{k,j} = -eigVec[k].ttt  */
      kk = k + Nvecs_max*j;
//      H[kk] = dcmplx((double)0.0, (double)0.0);
//      FORSOMEFIELDPARITY_OMP(i, parity, private(cc) reduction(+:H[kk])){
//	cc = su3_dot(eigVec[k]+i, ttt+i);
//	CSUM(H[kk], cc);
//      } END_LOOP_OMP;

      double cctotr=0., cctoti=0.;
      FORSOMEFIELDPARITY_OMP(i, parity, private(cc) reduction(+:cctotr,cctoti)){
	cc = su3_dot(eigVec[k]+i, ttt+i);
	cctotr += cc.real;
	cctoti += cc.imag;
      } END_LOOP_OMP;
      H[kk].real = cctotr;
      H[kk].imag = cctoti;


      CNEGATE(H[kk], H[kk]);
    }
    g_vecdcomplexsum(H+Nvecs_max*j, Nvecs_curr+Nvecs);
  }

  free(ttt);
}

/* Rayleigh-Ritz procedure */
/* H = Q^+AQ, Solve HZ = ZM -> eigVal[i] = M_{i,i}, eigVec = QZ */
static void RayleighRitz(int m, int n, double *eigVal, su3_vector **eigVec, double_complex *H,
			 int ldH, int parity){

  /* Constant */
  int lwork = 2*n+1;

  register int i;
  int j, k, info;
  double *rwork;
  double_complex *work;
  su3_vector *ttt;

  rwork = (double *)malloc(3*n*sizeof(double));
  work = (double_complex *)malloc(lwork*sizeof(double_complex));

  /* Solve HZ = ZM -> eigVal[j] = M_{j,j} */
  zheev_("V", "U", &n, H, &ldH, eigVal, work, &lwork, rwork, &info);

  free(rwork); free(work);

  ttt = (su3_vector *)malloc(m*sizeof(su3_vector));

  /* eigVec = QZ */
  FORSOMEFIELDPARITY(i, parity){
    for(j = 0; j < m; j++){
      clearvec( ttt+j );
      for(k = 0; k < n; k++)
	c_scalar_mult_add_su3vec(ttt+j, H+k+ldH*j, eigVec[k]+i);
    }
    for(j = 0; j < m; j++) eigVec[j][i] = ttt[j];
  } END_LOOP;

  free(ttt);
}

/* Calculate eigenpairs */
void calc_eigenpairs(double *eigVal, su3_vector **eigVec, eigcg_params *eigcgp, int parity){

  /* Constants */
  double dzero = (double)0.0;
  double_complex zzero = dcmplx(dzero, dzero);

  int i, j, Nvecs_curr, Nvecs_max;
  double_complex *H;

  Nvecs_curr = eigcgp->Nvecs_curr;
  Nvecs_max = eigcgp->Nvecs_max;
  H = eigcgp->H;

  RayleighRitz(Nvecs_curr, Nvecs_curr, eigVal, eigVec, H, Nvecs_max, parity);

  for(i = 0; i < Nvecs_curr; i++){
    for(j = 0; j < i; j++) H[j + Nvecs_max*i] = zzero;
    H[(Nvecs_max+1)*i] = dcmplx(eigVal[i], dzero);
  }
}

#if 0   /* Obsolete: See eigen_stuff_helpers for supported code */
void calc_eigresid(int Nvecs, double *resid, double *norm, double *eigVal,
		   su3_vector **eigVec, int parity, imp_ferm_links_t *fn){

  /* Constants */
  int otherparity = (parity == EVEN) ? ODD : EVEN;
  double dzero = (double)0.0;

  register int i;
  int j;
  double *dvec;
  su3_vector *ttt;

  ttt = (su3_vector *)malloc(sites_on_node*sizeof(su3_vector));
  dvec = (double *)malloc(2*Nvecs*sizeof(double));

  /* dvec[0] = ||Dslash^2*eigVec + eigVal*eigVec||^2,
     dvec[1] = ||eigVec||^2
  */
  for(j = 0; j < Nvecs; j++){
    dslash_fn_field(eigVec[j], ttt, otherparity, fn);
    dslash_fn_field(ttt, ttt, parity, fn);

    dvec[2*j] = dzero;
    dvec[2*j + 1] = dzero;
    FORSOMEFIELDPARITY_OMP(i, parity, reduction(+:dvec[2*j], dvec[2*j + 1])){
      scalar_mult_sum_su3_vector(ttt+i, eigVec[j]+i, eigVal[j]);
      dvec[2*j] += magsq_su3vec(ttt+i);
      dvec[2*j + 1] += magsq_su3vec(eigVec[j]+i);
    } END_LOOP_OMP;
  }
  g_vecdoublesum(dvec, 2*Nvecs);
  for(j = 0; j < Nvecs; j++){
    resid[j] = sqrt(dvec[2*j]/dvec[2*j + 1]);
    norm[j] = sqrt(dvec[2*j  + 1]);
  }

  free(dvec); free(ttt);
}
#endif

/* eigCG */
/* Solve a linear equation as well as calculate Nvecs eigenpairs of -Dslash^2 simultaniously. */
/* Lanczos part restarts after first m steps and then every (m - 2*Nvecs) steps */
/* This version looks at the initial vector every "niter" passes. */
/* The source vector is in "src", and the initial guess and answer
   in "dest".  "resid" is the residual vector, and "cg_p" and "ttt" are
   working vectors for the conjugate gradient.
   niter = maximum number of iterations before restarting.
   max_restarts = max number of restarts
   rsqmin = desired rsq, quit when we reach rsq <= rsqmin*source_norm.

   reinitialize after niters iterations and try once more.
*/
int ks_eigCG_parity(su3_vector *src, su3_vector *dest, double *eigVal, su3_vector **eigVec,
		    int m, int Nvecs, quark_invert_control *qic, Real mass, imp_ferm_links_t *fn){

  /* constants */
  int ione = 1;
  double dzero = (double)0.0;
  double done = (double)1.0;
  double_complex zzero = dcmplx(dzero, dzero);
  double_complex zone = dcmplx(done, dzero);
  Real msq_x4 = 4.0*mass*mass;
  /*** EigCG ***/
  int Nvecs_x2 = 2*Nvecs;
  int m1 = m - 1;
  int mm = m*m;
  int mNvecs = m*Nvecs;
  int lwork = 2*m;
  /*************/

  register int i;
  int j, iteration;	           /* counter for iterations */
  double a = done, b = dzero; /* Sugar's a,b */
  double rsq, relrsq;              /* resid**2, rel resid*2 */
  double pkp;	           /* pkp = cg_p.K.cg_p */
  double source_norm;              /* squared magnitude of source vector */
  msg_tag * tags1[16], *tags2[16]; /* tags for gathers to parity and opposite */
  int special_started;             /* 1 if dslash_fn_field_special has been called */
  int nrestart;                    /* Restart counter */
  su3_vector *ttt, *cg_p, *resid;
  double dtimec;
  char myname[] = "ks_eigCG_parity";

  /* Unpack structure */
  int niter        = qic->max;                    /* maximum number of iters per restart */
  int max_restarts = qic->nrestart;               /* maximum restarts */
  Real rsqmin      = qic->resid*qic->resid;       /* desired residual - 
						     normalized as sqrt(r*r)/sqrt(src_e*src_e) */
  Real relrsqmin   = qic->relresid*qic->relresid; /* desired relative residual (FNAL)*/
  int parity       = qic->parity;                 /* EVEN, ODD */

  int otherparity = (parity == EVEN) ? ODD : EVEN;
  int max_cg = max_restarts*niter;                /* Maximum number of iterations */

  /*** EigCG ***/
  int jj;
  int k = -1, info;
  int *iwork = NULL, *ifail = NULL;
  double *rwork = NULL;
  double_complex *T = NULL, *T2 = NULL, *Y = NULL, *tau = NULL, *work = NULL;
  su3_vector *ttt2 = NULL, *tmp = NULL;
  double dtimec2 = 0.0, dtimec3;
  /*************/

  if(fn == NULL){
    printf("%s(%d): Called with NULL fn\n", myname, this_node);
    terminate(1);
  }

  dtimec = -dclock(); 

  rsq = dzero;
  relrsq = done;
  special_started = 0;
  nrestart = 0;
  iteration = 0;
  qic->size_r        = 0;
  qic->size_relr     = 1.;
  qic->final_iters   = 0;
  qic->final_restart = 0;
  qic->converged     = 1;
  qic->final_rsq     = 0.;
  qic->final_relrsq  = 0.;

  /* Source norm */
  source_norm = dzero;
  FORSOMEFIELDPARITY_OMP(i, parity, reduction(+:source_norm)){
    source_norm += (double)magsq_su3vec(src+i);
  } END_LOOP_OMP;
  g_doublesum(&source_norm);
#ifdef CG_DEBUG
  node0_printf("congrad: source_norm = %e\n", (double)source_norm); fflush(stdout);
#endif

  /* Provision for trivial solution */
  if(source_norm == dzero){
    /* Zero the solution and return zero iterations */
    FORSOMEFIELDPARITY_OMP(i, parity, default(shared)){
      memset(dest+i, 0, sizeof(su3_vector));
    } END_LOOP_OMP;

    dtimec += dclock();
#ifdef CGTIME
    node0_printf("CONGRAD5_EIGCG: time = %e (fn_eigcg %s) masses = 1 iters = %d\n",
		 dtimec, prec_label[PRECISION-1], qic->final_iters);
#endif

    return 0;
  }

  /* Allocate temporary variables */
  /* PAD may be used to avoid cache trashing */
#define PAD 0
  ttt = (su3_vector *)malloc((sites_on_node+PAD)*sizeof(su3_vector));
  cg_p = (su3_vector *)malloc((sites_on_node+PAD)*sizeof(su3_vector));
  resid = (su3_vector *)malloc((sites_on_node+PAD)*sizeof(su3_vector));

  if(ttt == NULL || cg_p == NULL || resid == NULL){
    printf("%s(%d): No room for temporaries\n",myname,this_node);
  }

  /*** EigCG ***/
  if( Nvecs > 0 ){
    iwork = (int *)malloc(5*m*sizeof(int));
    ifail = (int *)malloc(m*sizeof(int));
    rwork = (double *)malloc(7*m*sizeof(double));
    T = (double_complex *)malloc(mm*sizeof(double_complex));
    T2 = (double_complex *)malloc(mm*sizeof(double_complex));
    Y = (double_complex *)malloc(m*Nvecs_x2*sizeof(double_complex));
    tau = (double_complex *)malloc(Nvecs_x2*sizeof(double_complex));
    work = (double_complex *)malloc(lwork*sizeof(double_complex));
    ttt2 = (su3_vector *)malloc((sites_on_node+PAD)*sizeof(su3_vector));
    tmp = (su3_vector *)malloc(Nvecs_x2*sizeof(su3_vector));

    if(rwork == NULL || T == NULL || T2 == NULL || Y == NULL || tau == NULL ||
	work == NULL || ttt2 == NULL || tmp == NULL){
      printf("%s(%d): No room for temporaries\n", myname, this_node);
    }
  }
  /*************/

  /* Start CG iterations */

#ifdef CG_DEBUG
  node0_printf("rsqmin = %g relmin = %g\n", rsqmin, relrsqmin); fflush(stdout);
#endif
  while(1) {
    /* Check for completion */
    /* Start at niter = 0 and restart at niter intervals.
       But if we have met tolerances with the cumulative residual,
       check the true residuals and quit if they have also met tolerances */
    if((iteration % niter == 0) || 
       ((rsqmin    <= 0 || rsqmin    > qic->size_r) &&
	(relrsqmin <= 0 || relrsqmin > qic->size_relr))){
	
      /* (re)initialization process */
	
      /* Compute true residual and relative residual */
      /* ttt <-  (-1)*M_adjoint*M*dest
	 resid,cg_p <- src + ttt
	 rsq = |resid|^2
	 source_norm = |src|^2
      */
      if(special_started==1) {	/* clean up gathers */
	cleanup_gathers(tags1,tags2);
	special_started=0;
      }
      rsq = dzero;
      dslash_fn_field(dest, ttt, otherparity, fn);
      dslash_fn_field(ttt, ttt, parity, fn);
      /* ttt  <- ttt - msq_x4*src	(msq = mass squared) */
      FORSOMEFIELDPARITY_OMP(i, parity, reduction(+:rsq)){
	scalar_mult_sum_su3_vector(ttt+i, dest+i, -msq_x4);
	add_su3_vector(src+i, ttt+i, resid+i);
	/* remember ttt contains -M_adjoint*M*src */
	cg_p[i] = resid[i];
	rsq += (double)magsq_su3vec(resid+i);
      } END_LOOP_OMP;
      g_doublesum(&rsq);

      if(relrsqmin > 0)
	relrsq = my_relative_residue(resid, dest, parity);

      qic->final_rsq    = (Real)rsq/source_norm;
      qic->final_relrsq = (Real)relrsq;

      iteration++ ; /* iteration counts number of multiplications
		       by M_adjoint*M */
      total_iters++;

#ifdef CG_DEBUG
      node0_printf("CONGRAD: (re)start %d rsq = %.10e relrsq %.10e\n",
		   nrestart, qic->final_rsq, qic->final_relrsq); fflush(stdout);
#endif
      /* Quit when true residual and true relative residual are within
	 tolerance or when we exhaust iterations or restarts */
	
      if(iteration >= max_cg || 
	 nrestart  >= max_restarts ||
	 ((rsqmin    <= 0 || rsqmin    > qic->final_rsq) &&
	  (relrsqmin <= 0 || relrsqmin > qic->final_relrsq))) break;

      /*** EigCG ***/
      if(Nvecs > 0){
	a = done;
	b = dzero;
	k = -1;
        for(j = 0; j < mm; j++) T[j] = zzero;
      }
      /*************/
	
      nrestart++;
    }

    /* ttt <- (-1)*M_adjoint*M*cg_p */
    if(special_started==0){
      dslash_fn_field_special(cg_p, ttt, otherparity, tags2, 1, fn);
      dslash_fn_field_special(ttt, ttt, parity, tags1, 1, fn);
      special_started=1;
    }
    else {
      dslash_fn_field_special(cg_p, ttt, otherparity, tags2, 0, fn);
      dslash_fn_field_special(ttt, ttt, parity, tags1, 0, fn);
    }    
    /* finish computation of M_adjoint*M*cg_p and cg_p*M_adjoint*M*cg_p */
    /* ttt  <- ttt - msq_x4*cg_p	(msq = mass squared) */
    /* pkp  <- cg_p.(ttt - msq*cg_p) = (-1)*cg_p.M_adjoint*M.cg_p */
    pkp = dzero;
    FORSOMEFIELDPARITY_OMP(i, parity, reduction(+:pkp)){
      scalar_mult_sum_su3_vector(ttt+i, cg_p+i, -msq_x4);
      pkp += (double)su3_rdot(cg_p+i, ttt+i);
    } END_LOOP_OMP;
    g_doublesum(&pkp);

    /*** EigCG ***/
    /* Compute Ritz pairs, then restart */
    if(Nvecs > 0){
      dtimec3 = - dclock();

#ifdef CG_DEBUG
      node0_printf("ks_eigCG_parity computing Ritz pairs with k %d and m1 %d\n",k,m1); fflush(stdout);
#endif    

      if(k == m1){
	/* Solve T_m Y1 = Y1 M1 for the lowest Nvecs eigenpairs */
	zcopy_(&mm, T, &ione, T2, &ione);
	zheevx_("V", "I", "U", &m, T2, &m, &dzero, &dzero, &ione, &Nvecs, &dzero, &k,
		eigVal, Y, &m, work, &lwork, rwork, iwork, ifail, &info);

	/* Solve T_{m-1} Y2 = Y2 M2 for the lowest Nvecs eigenpairs */
	/* Now Y = [Y1, Y2] */
	zcopy_(&mm, T, &ione, T2, &ione);
	zheevx_("V", "I", "U", &m1, T2, &m, &dzero, &dzero, &ione, &Nvecs, &dzero, &k,
		eigVal, Y+mNvecs, &m, work, &lwork, rwork, iwork, ifail, &info);
	for(j = Nvecs; j <= Nvecs_x2; j++) Y[m*j-1] = zzero;

	/* Orthogonalize Y -> QR, Y <- Q */
	zgeqrf_(&m, &Nvecs_x2, Y, &m, tau, work, &lwork, &info);
	zungqr_(&m, &Nvecs_x2, &Nvecs_x2, Y, &m, tau, work, &lwork, &info);

	/* T_m <- Q^+ T_m Q */
	zhemm_("L", "U", &m, &Nvecs_x2, &zone, T, &m, Y, &m, &zzero, T2, &m);
	zgemm_("C", "N", &Nvecs_x2, &Nvecs_x2, &m, &zone, Y, &m, T2, &m, &zzero, T,
	       &Nvecs_x2);

	/* Solve T Z = Z M, eigVal[j] <- M_{j,j}, T <- Z */
	zheev_("V", "U", &Nvecs_x2, T, &Nvecs_x2, eigVal, work, &lwork, rwork, &info);

	/* eigVec[j] <- sum_jj eigVec[jj]*(QZ)[jj][j] */
	zgemm_("N", "N", &m, &Nvecs_x2, &Nvecs_x2, &zone, Y, &m, T, &Nvecs_x2, &zzero,
	       T2, &m);

	FORSOMEFIELDPARITY(i, parity){
	  for(j = 0; j < Nvecs_x2; j++){
	    clearvec( tmp+j );
	    for(jj = 0; jj < m; jj++)
	      c_scalar_mult_add_su3vec(tmp+j, T2+jj+m*j, eigVec[jj]+i);
	  }
	  for(j = 0; j < Nvecs_x2; j++) eigVec[j][i] = tmp[j];
	} END_LOOP;
	
	/* T <- M */
	for(j = 0; j < Nvecs_x2; j++){
	  for(jj = 0; jj < j; jj++) T[jj + m*j] = zzero;
	  T[(m+1)*j] = dcmplx(eigVal[j], dzero);
	}

        k = Nvecs_x2 - 1;

	/* ttt2 < ttt2 - ttt, where ttt2 = b*ttt_old */
	FORSOMEFIELDPARITY_OMP(i, parity, default(shared)){
	  sub_su3_vector(ttt2+i, ttt+i, ttt2+i);
        } END_LOOP_OMP;

	/* T_{j,k+1} = eigVec[j].ttt2/sqrt(rsq) */
	for(j = 0; j < Nvecs_x2; j++){

//	  tau[j] = zzero;
//	  FORSOMEFIELDPARITY_OMP(i, parity, private(work[0]) reduction(+:tau[j])){
//	    work[0] = su3_dot(eigVec[j]+i, ttt2+i);
//	    CSUM(tau[j], work[0]);
//	  } END_LOOP_OMP;

	  double cctotr=0., cctoti=0.;
	  double_complex cc;
	  FORSOMEFIELDPARITY_OMP(i, parity, private(cc) reduction(+:cctotr,cctoti)){
	    cc = su3_dot(eigVec[j]+i, ttt2+i);
	    cctotr += cc.real;
	    cctoti += cc.imag;
	  } END_LOOP_OMP;
	  tau[j].real = cctotr;
	  tau[j].imag = cctoti;

	}
	g_vecdcomplexsum(tau, Nvecs_x2);
	for(j = 0; j < Nvecs_x2; j++)
	  CDIVREAL(tau[j], sqrt(rsq), T[j + m*(k+1)]);

#ifdef EIGCG_DEBUG
	if(special_started==1) {	/* clean up gathers */
	  cleanup_gathers(tags1,tags2);
	  special_started=0;
	}
	for(j = 0; j < Nvecs; j++){
	  dslash_fn_field(eigVec[j], ttt2, otherparity, fn);
	  dslash_fn_field(ttt2, ttt2, parity, fn);

	  
	  double rw0 = 0., rw1 = 0.;
	  //rwork[2*j] = dzero;
	  //rwork[2*j + 1] = dzero;
	  //FORSOMEFIELDPARITY_OMP(i, parity, reduction(+:rwork[2*j], rwork[2*j + 1])){
	  FORSOMEFIELDPARITY_OMP(i, parity, reduction(+:rw0, rw1)){
	    scalar_mult_sum_su3_vector(ttt2+i, eigVec[j]+i, eigVal[j] - msq_x4);
	    rw0 += magsq_su3vec(ttt2+i);
	    rw1 += magsq_su3vec(eigVec[j]+i);
	  } END_LOOP_OMP;
	  rwork[2*j] = rw0;
	  rwork[2*j + 1] = rw1;
        }
	g_vecdoublesum(rwork, 2*Nvecs);
	for(j = 0; j < Nvecs; j++){
	  rwork[2*j] = sqrt(rwork[2*j]/rwork[2*j + 1]);
	  rwork[2*j + 1] = sqrt(rwork[2*j + 1]);
	  node0_printf("eigVal[%d] = %e, resid = %e, ||eigVec[%d]|| = %e\n",
		       j, eigVal[j] - msq_x4, rwork[2*j], j, rwork[2*j + 1]);
	}
#endif
      }
      else if( k >= 0 ){
	/* T_{k,k+1} = -sqrt(b)/a */
	T[k + m*(k+1)] = dcmplx(-sqrt(b)/a, dzero);
      }

      k++;

      /* eigVec[k] = resid/sqrt(rsq) */
      FORSOMEFIELDPARITY_OMP(i, parity, default(shared)){
	scalar_mult_su3_vector(resid+i, done/sqrt(rsq), eigVec[k]+i);
      } END_LOOP_OMP;

      /* T_{k,k} = 1/a + b_old/a_old */
      T[(m+1)*k] = dcmplx(b/a, dzero);

      dtimec2 += dclock() + dtimec3;
    }
    /*************/

    /* a <- -rsq/pkp */
    a = (Real)(-rsq/pkp);

    /*
      dest <- dest + a*cg_p
      resid <- resid + a*ttt 
      rsq <- |resid|^2
    */
    b = rsq;
    rsq = dzero;
    FORSOMEFIELDPARITY_OMP(i, parity, reduction(+:rsq)){
      scalar_mult_sum_su3_vector(dest+i, cg_p+i, a);
      scalar_mult_sum_su3_vector(resid+i, ttt+i, a);
      rsq += (double)magsq_su3vec(resid+i);
    } END_LOOP_OMP;
    g_doublesum(&rsq);

    if(relrsqmin > 0)
      relrsq = my_relative_residue(resid, dest, parity);

    iteration++;
    total_iters++;
    
    qic->size_r        = (Real)rsq/source_norm;
    qic->size_relr     = relrsq;
    qic->final_iters   = iteration;
    qic->final_restart = nrestart;
    qic->converged     = 1;

#ifdef CG_DEBUG
    node0_printf("iter=%d, rsq/src= %e, relrsq= %e, pkp=%e\n",
		 iteration, (double)qic->size_r, (double)qic->size_relr, (double)pkp);
    fflush(stdout);
#endif

    /* b <- rsq/oldrsq */
    b = (Real)rsq/b;

    /* cg_p  <- resid + b*cg_p */
    FORSOMEFIELDPARITY_OMP(i, parity, default(shared)){
      scalar_mult_add_su3_vector(resid+i, cg_p+i, b, cg_p+i);
    } END_LOOP_OMP;

    /*** EigCG ***/
    if( Nvecs > 0 ){
      dtimec3 = -dclock();

      /* T_{k,k} = 1/a + b_old/a_old */
      T[(m+1)*k].real += done/a;

      /* ttt2 = b*ttt */
      if(k == m1){
	FORSOMEFIELDPARITY_OMP(i, parity, default(shared)){
	  scalar_mult_su3_vector(ttt+i, b, ttt2+i);
	} END_LOOP_OMP;
      }

      dtimec2 += dclock() + dtimec3;
    }
    /*************/
  }

  if(special_started==1) {
    cleanup_gathers(tags1, tags2);
    special_started = 0;
  }
  cleanup_dslash_temps();

  /*** EigCG ***/
  /* compute final eigenpairs */
  if( Nvecs > 0 ){
    dtimec3 = -dclock();

    k++;

    RayleighRitz(Nvecs, k, eigVal, eigVec, T, m, parity);
    for(j = 0; j < Nvecs; j++) eigVal[j] -= msq_x4;

    dtimec2 += dclock() + dtimec3;

#ifdef EIGCG_DEBUG
    check_eigres( rwork, eigVec, eigVal, Nvecs, parity, fn );
//    calc_eigresid(Nvecs, rwork, rwork+Nvecs, eigVal, eigVec, parity, fn);
//    for(j = 0; j < Nvecs; j++)
//      node0_printf("eigVal[%d] = %e, resid = %e, ||eigVec[%d]|| = %e\n",
//		   j, eigVal[j], rwork[j], j, rwork[j + Nvecs]);
#endif
  }
  /*************/

  if(nrestart == max_restarts || iteration == max_cg){
    qic->converged = 0;
    node0_printf("%s: CG not converged after %d iterations and %d restarts, \n",
		 myname, iteration, nrestart);
    node0_printf("rsq. = %e wanted %e relrsq = %e wanted %e\n",
		 qic->final_rsq, rsqmin,qic->final_relrsq, relrsqmin);
  }

  free(ttt); free(cg_p); free(resid);
  /*** EigCG ***/
  if( Nvecs > 0 ){
    free(rwork); free(T); free(T2); free(Y); free(tau); free(work); free(ttt2); free(tmp);
  }
  /*************/

  dtimec += dclock();
#ifdef CGTIME
  node0_printf("CONGRAD5_EIGCG: time = %e time_eig = %e (fn_eigcg %s) masses = 1 iters = %d\n",
	       dtimec, dtimec2, prec_label[PRECISION-1], qic->final_iters);
#endif

  return iteration;
}

/* Incremental eigCG */
/* !!! This computes eigenpairs of -Dslash^2 !!! */
/* Since this routine only comupte H = -U^+ Dslash^2 U, calc_eigenpairs() must be called
   if eigenpairs are really needed. */
int ks_inc_eigCG_parity( su3_vector *src, su3_vector *dest, double *eigVal,
			 su3_vector **eigVec, eigcg_params *eigcgp, quark_invert_control *qic,
			 Real mass, imp_ferm_links_t *fn ){

  int i, parity, m, Nvecs, Nvecs_curr, Nvecs_max, Nvecs_add, iteration;
  double dtimec, dtimec2=0.0, dtimec3, dtimec4=0.0, dtimec5=0.0;
  double_complex *H;
#ifdef EIGCG_DEBUG
  double *resid_debug, *norm;
#endif

  dtimec = -dclock();

  parity = qic->parity;
  m = eigcgp->m;
  Nvecs = eigcgp->Nvecs;
  Nvecs_curr = eigcgp->Nvecs_curr;
  Nvecs_max = eigcgp->Nvecs_max;

  if(Nvecs_curr == 0){
    H = (double_complex *)malloc(Nvecs_max*Nvecs_max*sizeof(double_complex));
    for(i = 0; i < Nvecs_max*Nvecs_max; i++) H[i] = dcmplx((double)0.0, (double)0.0);
    if(eigcgp->H != NULL) free(eigcgp->H);
    eigcgp->H = H;
  }
  else{
    /* Deflation */
    H = eigcgp->H;
    dtimec2 = -dclock();
    initCG(src, dest, Nvecs_curr, Nvecs_max, eigVec, H, mass, parity, fn);
    dtimec2 += dclock();
  }

#ifdef CG_DEBUG
    node0_printf("Calling ks_eigCG_parity with Nvecs_curr = %d\n", Nvecs_curr); fflush(stdout);
#endif    

  /* Solve a linear equation */
  dtimec3 = -dclock();
  //iteration = ks_eigCG_parity(src, dest, eigVal+Nvecs_curr, eigVec+Nvecs_curr, m, Nvecs,
			       //		      qic, mass, fn);
  iteration = ks_eigCG_parity(src, dest, eigVal+Nvecs_curr, eigVec+Nvecs_curr, m, Nvecs,
			      qic, mass, fn);
  dtimec3 += dclock();

  if(Nvecs > 0){

#ifdef CG_DEBUG
    node0_printf("Orthogonalization step with Nvecs_curr = %d\n", Nvecs_curr); fflush(stdout);
#endif    

    /* Orthogonalize vectors */
    dtimec4 = -dclock();
    Nvecs_add = orthogonalize(Nvecs, Nvecs_curr, eigVec, parity);
    dtimec4 += dclock();

    /* Construct H = -U^+ Dslash^2 U */
    dtimec5 = -dclock();
    extend_H(Nvecs_add, Nvecs_curr, Nvecs_max, eigVec, H, parity, fn);
    dtimec5 += dclock();

    Nvecs_curr += Nvecs_add;
    eigcgp->Nvecs_curr = Nvecs_curr;
    eigcgp->Nvecs = (Nvecs_max - Nvecs_curr < Nvecs) ? (Nvecs_max - Nvecs_curr) : Nvecs;
  }

  dtimec += dclock();

#ifdef EIGCG_DEBUG
  //  norm = (double *)malloc(Nvecs_curr*sizeof(double));
  calc_eigenpairs(eigVal, eigVec, eigcgp, parity);

  resid_debug = (double *)malloc(Nvecs_curr*sizeof(double));
  check_eigres( resid_debug, eigVec, eigVal, Nvecs_curr, parity, fn );
//  calc_eigresid(Nvecs_curr, resid_debug, norm, eigVal, eigVec, parity, fn);
//  for(i = 0; i < Nvecs_curr; i++)
//    node0_printf("eigVal[%d] = %e, resid = %e, ||eigVec[%d]|| = %e\n",
//		 i, eigVal[i], resid[i], i, norm[i]);
  free(resid_debug); // free(norm);
#endif

#ifdef CGTIME
  node0_printf("INC_EIGCG: time for init-CG           %e\n", dtimec2);
  node0_printf("INC_EIGCG: time for eigCG             %e\n", dtimec3);
  node0_printf("INC_EIGCG: time for orthogonalization %e\n", dtimec4);
  node0_printf("INC_EIGCG: time for extending H       %e\n", dtimec5);
  node0_printf("INC_EIGCG: total time                 %e\n", dtimec);
#endif

  return iteration;
}
