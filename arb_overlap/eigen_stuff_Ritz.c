/****** eigen_stuff_Ritz.c  ******************/
/* Eigenvalue and Eigevector computation routines.
* K.O. 8/99 Started. 
* MIMD version 7
*
*  These routines are for the computation of the Eigenvalues and Eigevectors
* of ``any'' Hermetian matrix acting on data type wilson_vector 
*/
/**#define DEBUG**/
#define MINITER 5
#include "arb_ov_includes.h"	/* definitions files and prototypes */

#ifdef SSE
#define SSE_SUBS
#include "../sse/include/inline_sse.h"
#endif

#if 0
/* This routine is the routine that applies the Matrix, whose eigenvalues    *
 * we want to compute, to a vector. */
static void Matrix_Vec_mult(wilson_vector *src, wilson_vector *res)
{  
  register site *s;
  register  int i;
  /* chirality check */
  complex ctmp,cn;
  Real chirality,cd;
  wilson_vector wtmp;

  if(kind_of_h0 == HOVERLAP){
    /* check src for chirality */
    cn=cmplx(0.0,0.0);
    cd=0.0;
    FORALLSITES(i,s){
      mult_by_gamma(&(src[i]),&wtmp,GAMMAFIVE);
      cd += magsq_wvec(&(src[i]));
      ctmp =  wvec_dot(&(src[i]),&wtmp);
      CSUM(cn,ctmp);
    }
    g_complexsum(&cn);
    g_floatsum(&cd);
    chirality=cn.real/cd;
 /*if(fabs(chirality) < 0.95) node0_printf("chirality %e\n",chirality);*/
    if(chirality >0.0) chirality_flag=1;
    if(chirality <0.0) chirality_flag= -1;

    hdelta0_field(src,res,kind_of_h0);
    /*
      node0_printf(" Matrix_Vec_mult: chirality_flag %d\n",chirality_flag);
    */
  }
  else 
  {
      hdelta0_field(src,res,kind_of_h0);
  }

}
#endif

int Rayleigh_min(wilson_vector *vec, wilson_vector **eigVec, Real Tolerance, 
		 Real RelTol, int Nvecs, int MaxIter, int Restart)
{
  int iter ;
  double beta, cos_theta, sin_theta ;
  double quot, old_quot, P_norm, theta, real_vecMp, pMp ;
#ifdef DEBUG
  double vec_norm ;
#endif  
  double g_norm, old_g_norm, start_g_norm ;
  double_complex cc ;
  wilson_vector *Mvec, *grad, *P, *MP ;
  
  Mvec     = (wilson_vector *)malloc(sites_on_node*sizeof(wilson_vector)) ;
  grad     = (wilson_vector *)malloc(sites_on_node*sizeof(wilson_vector)) ;
  P        = (wilson_vector *)malloc(sites_on_node*sizeof(wilson_vector)) ;
  MP       = (wilson_vector *)malloc(sites_on_node*sizeof(wilson_vector)) ;


  old_quot=1.0e+16;
  /*
    printf("\n  In Rayleigh, MaxIter=%d\n",MaxIter);
    */
  project_out(vec, eigVec, Nvecs);
  normalize(vec) ; 
  Matrix_Vec_mult(vec,Mvec) ;
  project_out(Mvec, eigVec, Nvecs);
  
  /* Compute the quotient quot=vev*M*vec */
  dot_product(vec, Mvec, &cc) ;
  /* quot is real since M is hermitian. quot = vec*M*vec */
  quot = cc.real ;
#ifdef DEBUG
  node0_printf("Rayleigh_min: Start -- quot=%g,%g\n", quot,cc.imag) ;
#endif  
  /* Compute the grad=M*vec - quot*vec */
  copy_Vector(Mvec,grad) ;
  double_vec_mult_sub(&quot, vec, grad) ;
  /* set P (the search direction) equal to grad */
  copy_Vector(grad,P) ;
  /* compute the norms of P and grad */
  norm2(P   , &P_norm) ;
  g_norm=P_norm;
  start_g_norm = g_norm ;
#ifdef DEBUG
  node0_printf("Rayleigh_min: Start -- g_norm=%g\n", g_norm) ;
#endif

  iter = 0 ;
  while( (g_norm>Tolerance)&&
         ( ((iter<MaxIter)&&(g_norm/start_g_norm>RelTol)) || (iter<MINITER) )
         )
    {
      iter++ ;
      Matrix_Vec_mult(P,MP) ;
      dot_product(vec, MP, &cc) ; 
      real_vecMp = cc.real ;
      dot_product(P, MP, &cc) ; 
      pMp = cc.real ; /*pMp is real */
      theta = 0.5*atan(2.0*real_vecMp/(quot*P_norm - pMp/P_norm)) ;
      sin_theta = sin(theta) ;
      cos_theta = cos(theta) ;
      if(sin_theta*cos_theta*real_vecMp>0)
	{ 
	  theta = theta - 0.5*PI ; /* chose the minimum not the maximum */
	  sin_theta = sin(theta) ; /* the sin,cos calls can be avoided */
	  cos_theta = cos(theta) ;
	}
      sin_theta = sin_theta/P_norm ;
      /* vec = cos(theta)*vec +sin(theta)*P/p_norm */
      dax_p_by(&cos_theta, vec,&sin_theta, P) ; 
      /* Mvec = cos(theta)*Mvec +sin(theta)*MP/p_norm */
      dax_p_by(&cos_theta,Mvec,&sin_theta,MP) ; 
      /* renormalize vec ... */
      if(iter%Restart == 0 ) {
#ifdef DEBUG
	node0_printf("Renormalizing...");
	norm2(vec,&vec_norm) ;
	node0_printf("  norm: %g\n",1.0/vec_norm);
#endif
	/* Project vec on the orthogonal complement of eigVec */ 
	project_out(vec, eigVec, Nvecs);
	normalize(vec) ; 
	Matrix_Vec_mult(vec,Mvec);
	/* Recompute the quotient */
	dot_product(vec, Mvec, &cc) ;
	/* quot is real since M is hermitian. quot = vec*M*vec */
	quot = cc.real ;
	/* Recompute the grad */
	copy_Vector(Mvec,grad) ;
	double_vec_mult_sub(&quot, vec, grad) ;
	norm2(grad,&g_norm) ; /* recompute the g_norm */
	/* Project P on the orthogonal complement of eigVec */
	project_out(P, eigVec, Nvecs);
	/* make P orthogonal to vec */
	dot_product(vec, P, &cc) ;
	complex_vec_mult_sub(&cc, vec, P);
	/* make P orthogonal to grad */
	dot_product(grad, P, &cc) ;
	complex_vec_mult_sub(&cc, grad, P);
	norm2(P,&P_norm) ; /* recompute the P_norm */
      }
      dot_product(vec, Mvec, &cc) ;
      /* quot is real since M is hermitian. quot = vec*M*vec */
      quot = cc.real ;
#ifdef DEBUG
      node0_printf("Rayleigh_min: %i, quot=%8g g=%8g b=%6g P:%6g\n",
		   iter,cc.real,g_norm,beta,P_norm) ;    
#endif
      old_g_norm = g_norm ;
       
      copy_Vector(Mvec,grad) ;
      double_vec_mult_sub(&quot, vec, grad) ;

      norm2(grad,&g_norm) ;
      beta = cos_theta*g_norm*g_norm/(old_g_norm*old_g_norm) ;
      /* Cut off beta */
      if( beta>2.0 ) 
	beta = 2.0 ;
      dot_product(vec, P, &cc) ;
      cc.real *= beta ; cc.imag *= beta ;
      /*      P_norm = 1.0/P_norm ;
	      double_vec_mult(&P_norm,  P,  P) ; */
      vec_plus_double_vec_mult(grad, &beta, P) ;/* P = grad + beta*P */
      complex_vec_mult_sub(&cc, vec, P) ; /* P = P - cc*vec */
      norm2(P   , &P_norm) ;
/*
  node0_printf("DEBUG: iter %d  g_norm %e Tolerance %e checkratio %e RelTol %e\n",
		     iter,g_norm,Tolerance,g_norm/start_g_norm,RelTol);
*/
  if(fabs(old_quot -quot)< Tolerance/100.0) g_norm= Tolerance/10.0; /*break*/
      old_quot = quot;
    }
  project_out(vec, eigVec, Nvecs);
  normalize(vec) ;
  /*  cleanup_Matrix() ; */
  free(MP) ;
  free(P) ;
  free(grad) ;
  free(Mvec) ;
  iter++ ;
/*
  node0_printf("in Rayleigh %d\n",iter);
*/
  return iter ;
}

int Kalkreuter_Ritz(wilson_vector **eigVec, double *eigVal, Real Tolerance, 
	       Real RelTol, int Nvecs, int MaxIter, 
	       int Restart, int Kiters, int parity)
{
  int total_itns=0 ;
  int j ;
  Matrix Array,V ;
  register  int i ;
  wilson_vector *vec ;

  wilson_vector **MeigVec ;

  double max_error = 1.0e+10 ;
  double *grad, *err ;
  int iter = 0 ;
  int *converged;
  Real ToleranceG;
    /* code to adaptively vary the step function */
  Real resid_inner_hold;

  resid_inner_hold=resid_inner;
  ToleranceG= 10.0*Tolerance;
node0_printf("gradient tolerance %e eigVal change tolerance %e\n",
	     ToleranceG, Tolerance);

  /** Allocate the array **/
  Array = AllocateMatrix(Nvecs) ;  
  /** Allocate the Eigenvector matrix **/
  V = AllocateMatrix(Nvecs) ;  

  vec = (wilson_vector *)malloc(sites_on_node*sizeof(wilson_vector));
  grad = (double *) malloc(Nvecs*sizeof(double)) ;
  err = (double *) malloc(Nvecs*sizeof(double)) ;

  for(j=0;j<Nvecs;j++)err[j]=1.e6;

  converged = (int *)malloc(Nvecs*sizeof(int));
  for(j=0;j<Nvecs;j++)converged[j]=0;

		  
  /* Initiallize all the eigenvectors to a random vector 
  for(j=0;j<Nvecs;j++)
    {
      grad[j] = 1.0e+10 ;
      grsource(EVENANDODD);  
      FORALLSITES(i,s){
	copy_wvec(&(s->g_rand),&(eigVec[j][i]));
      }
      eigVal[j] = 1.0e+16 ;
    }
*/
  for(j=0;j<Nvecs;j++){grad[j] = 1.0e+10 ;/* eigVal[j] = 1.0e+16 ;*/}


  MeigVec = (wilson_vector **)malloc(Nvecs*sizeof(wilson_vector*));
  for(i=0;i<Nvecs;i++)
    MeigVec[i]=
      (wilson_vector*)malloc(sites_on_node*sizeof(wilson_vector));

  while((max_error>Tolerance)&&(iter<Kiters))
    {
      iter++ ;
      for(j=0;j<Nvecs;j++)
	{
	  if(grad[j]>(ToleranceG)) 
	    {
#ifdef MINNI
	      if(kind_of_h0 == HOVERLAP){
        /* adjust resid_inner */
		resid_inner= 0.001*grad[j];
	/* reset the j=0 resid_inner_run */
		if(resid_inner > 1.e-3) resid_inner=1.e-3;
		if(resid_inner < resid_inner_hold) resid_inner=resid_inner_hold;
	/* freeze the CG accuracy for the lowest mode or after some (4) sweeps */
		do_minn=1;
		if(j==0 || iter > 4){
		  do_minn=0;
		  resid_inner=resid_inner_hold;
		  resid_inner_run=resid_inner_hold;
		}
		node0_printf("Residuals:%d  running  %e desired eps  %e\n",
			j,resid_inner_run,resid_inner);
		}
#endif
	      converged[j]=0;
	      copy_Vector(eigVec[j],vec) ;
	      total_itns += Rayleigh_min(vec, eigVec,
					  ToleranceG, RelTol, 
					  j, MaxIter , Restart) ;
	      copy_Vector(vec,eigVec[j]) ;
	    }
	  else{
	    node0_printf("vector %d converged\n",j);
	    converged[j]=1;
	  }
	}

/*
      node0_printf("running Rit %d\n",total_itns);
*/


      /* if you didn't act on eigVec[i] last time, converged[i]=1,
	 and  MeigVec hasn't changed, so don't compute it */
      constructArray(eigVec,MeigVec, &Array, converged) ;

#ifdef DEBUG
      node0_printf("Eigenvalues before diagonalization\n");
      for(i=0;i<Nvecs;i++)
	node0_printf("quot(%i) = %g |grad|=%g\n",i,Array.M[i][i].real,grad[i]);
#endif
      Jacobi(&Array, &V, JACOBI_TOL) ;
      sort_eigenvectors(&Array,&V) ;
      RotateBasis(eigVec,&V) ;
      RotateBasis(MeigVec,&V) ;
     /*  recompute gradient--vec=MeigVec - eigVal*eigVec */
      /* recall we are recycling MeigVec's */

      for(i=0;i<Nvecs;i++) {
	Array.M[i][i].imag=0.0;
	copy_Vector(MeigVec[i],vec) ;
	double_vec_mult_sub(&Array.M[i][i].real,eigVec[i],vec) ;
	norm2(vec, &grad[i]) ;
      }


      /* find the maximum error */
      max_error = 0.0 ;
      for(i=0;i<Nvecs;i++)
        {
          err[i] = eigVal[i] ; /* the previous eigVal */
          eigVal[i] = Array.M[i][i].real ; /* the one from this iteration */
          err[i] = fabs(err[i] - eigVal[i])/(1.0 - RelTol*RelTol) ;
          if(err[i]>max_error) 
            max_error = err[i] ;
        }




      /* old code--see kostas_29_feb.tex
	{
	  if(grad[i]>max_error) 
	    max_error = grad[i] ;
	  err[i] = eigVal[i] ;
	  eigVal[i] = Array.M[i][i].real ;
	  err[i] = fabs(err[i] - eigVal[i])/(1.0 - RelTol*RelTol) ;
	}
	*/



      node0_printf("\nEigenvalues after diagonalization at iteration %i\n",
		   iter);
      for(i=0;i<Nvecs;i++)
	node0_printf("quot(%i) = %g +/- %8e |grad|=%g\n",
		     i,eigVal[i],err[i],grad[i]);

#ifdef MINNI
      if(kind_of_h0 == HOVERLAP)node0_printf("RSD %e\n",resid_inner_run);
#endif
      fflush(stdout);

      /* dump eigenvals and eigenvecs */
      if(kind_of_h0 == HOVERLAP){

	if( out_hov_flag != FORGET ){
	  write_eigen(eigVec,Nvecs,out_hov_flag,out_hov);
	  write_eigenval(eigVal,Nvecs,out_hov);
	}

      }

    }
  /*
  for(i=0;i<Nvecs;i++)
    node0_printf("Eigenvalue(%i) = %g +/- %8e \n",i,eigVal[i],err[i]);
    */
  /** Deallocate the arrays **/

  resid_inner=resid_inner_hold;
#ifdef MINNI
  resid_inner_run=resid_inner;
#endif

  deAllocate(&V) ;
  deAllocate(&Array) ;
  free(err);
  free(grad);
  free(vec);

  free(converged);
  for(i=0;i<Nvecs;i++) free(MeigVec[i]);
  free(MeigVec);

  return total_itns ;
}

