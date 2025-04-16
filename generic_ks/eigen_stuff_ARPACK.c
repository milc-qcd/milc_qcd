/****** eigen_stuff_ARPACK_poly.c   ******************/
/* Eigenvalue and Eigevector computation routines.
 *
 *  Apply polynomial acceleration.
 *
* This version uses ARPACK (July 2017)
* MIMD version 7
*
*  These routines are for the computation of the Eigenvalues and Eigevectors
* of the Kogut-Susskind dslash^2. 
*/


/* If you do not define the USE_DSLASH_SPECIAL then calls to the standard   *
 * dslash are made. If you do define USE_DSLASH_SPECIAL then DSLASH_SPECIAL *
 * is used                                                                  */
#ifdef FN
#define USE_DSLASH_SPECIAL
#endif


enum accel_type { NO_POLY, POLY_X , MSQ_POLY , CHEBY , POLY_INV , NO_POLY_SQRT , CHEBY_ORG  }  ;

/* Include files */

#include "generic_ks_includes.h"
#include "../include/arpack.h"
#include "../include/openmp_defs.h"
#include <string.h>

#ifdef MPI_COMMS
#include <mpi.h>
#endif

/***  coefficients of polynomial interpolataion **/
static Real *poly ;
static Real debug_val = 0.5 ;  /***  debug value **/

void Matrix_Vec_mult_sqrt(su3_vector *src, su3_vector *res, 
			  ks_eigen_param *eigen_param,
 			  imp_ferm_links_t *fn ) ;


/************************************************************************/

/**** some debug routines   ***/

static void constant_vec(su3_vector *v , Real c){
  
  v->c[0].real =  c ;
  v->c[1].real =  c ;
  v->c[2].real =  c ;
  
  v->c[0].imag = 0.0;
  v->c[1].imag = 0.0;
  v->c[2].imag = 0.0;
  
}



static void constant_mult(su3_vector *in, Real c ,su3_vector *out) {
  
  out->c[0].real =  c * in->c[0].real   ;
  out->c[1].real =  c * in->c[1].real   ;
  out->c[2].real =  c * in->c[2].real   ;
  
  out->c[0].imag =  c * in->c[0].imag   ;
  out->c[1].imag =  c * in->c[1].imag   ;
  out->c[2].imag =  c * in->c[2].imag   ;
  
}


static void dump_vec(su3_vector *v){
  
  node0_printf("vec_[0]_re = %g\n" , v->c[0].real ) ;
  node0_printf("vec_[1]_re = %g\n" , v->c[1].real ) ;
  node0_printf("vec_[2]_re = %g\n" , v->c[2].real ) ;
  
}




static void constant_vec_even(su3_vector *v , Real c)
{
  int i ;
  
  FOREVENFIELDSITES_OMP(i,)
  {
    constant_vec(v + i, c) ;
  } END_LOOP_OMP;
  
  
}






static void constant_mult_even(su3_vector *in, Real c ,su3_vector *out)
{
  int i ;
  
  FOREVENFIELDSITES_OMP(i,)
  {
    constant_mult(in + i, c, out + i) ;
  } END_LOOP_OMP;
  
  
}




static void dump_vec_even(su3_vector *v)
{
  int i = 0 ;
  
  node0_printf("Site %d\n" , i );
  dump_vec(v + i) ;
  
  
}

/* Bubble sort eigenvalues and create hash table for eigenvector
   reordering */
static int *sort_eigenvalues(double eigVal[],int Nvecs)
{
  int *hash = malloc(sizeof(int)*Nvecs);
  double fu;
  int i,j,k,m;
  for(i=0; i<Nvecs; i++)
    hash[i] = i;

  for(i=0; i<Nvecs; i++)
    for (j=1; j<Nvecs-i; j++)
      {
	k = j-1 ; 
	if(eigVal[j] < eigVal[k]) /* swicth them */
	  {
	    /* switch the eigenvalues */
	    fu = eigVal[k] ;
	    eigVal[k] = eigVal[j] ;
	    eigVal[j] = fu ;
	    /* switch the eigenvector hash index */
	    m = hash[k] ;
	    hash[k] = hash[j] ;
	    hash[j] = m ;
	  }
      }
  return hash;
}



/**************************************************/



/**
   This routines assumes that the parity is EVEN.
   
**/


static void Matrix_Vec_mult_poly(su3_vector *src, su3_vector *res, 
				 ks_eigen_param *eigen_param,
				 imp_ferm_links_t *fn)
{
  ks_eigen_poly *poly_param = &eigen_param->poly;
  int parity = eigen_param->parity;
  su3_vector *old = create_v_field();
  su3_vector *tmp = create_v_field();
  su3_vector  *sv = create_v_field();
  
  su3_vector *older = create_v_field() ;
  int j , i  ; 
  
  Real a = poly_param->minE  ;
  Real b = poly_param->maxE  ;
  int norder = poly_param->norder ;
  
  clearvec(old) ;  /** old = 0 **/
  clearvec(older) ;  /** older = 0 **/
  
  Real norm  = 4.0 / (b -a)   ;
  Real norm2 = 2.0 / (b -a)   ;
  
  
  //    constant_vec_even(src,1.0 ) ;  /** debug  ***/
  
  /** apply Clenshaw's recursion ***/
  
  for (j=norder-1;j>=1;j--) 
    { 
      /**
	 y=(2.0*x-a-b)/(b-a); 
	 y2=2.0*y  ;
	 sv=d;
	 d=y2*d-dd+c[j];
	 dd=sv;
      **/
      
      
      
      Matrix_Vec_mult(old,tmp,eigen_param,fn);  /**  tmp = M *  old **/
      /**  constant_mult_even(old , debug_val, tmp) ;   ***/
      
      FOREVENFIELDSITES_OMP(i,)
      {
	su3vec_copy(old +i , sv + i ) ;                        /* sv = old  **/
	scalar_mult_su3_vector(tmp+i, norm, old+i) ;           /* old = norm * M * old   ***/
	scalar_mult_sub_su3_vector(old+i, sv+i, 2.0 , old+i) ; /* old = old - 2*old **/
	
	
	sub_su3_vector(old+i, older+i, old+i )  ;              /* old = old - older  **/           
	scalar_mult_add_su3_vector(old+i,src + i, poly[j],  old+i )  ;  /** old = old + poly[j] * src **/
	
	
	su3vec_copy(sv + i , older +i ) ;   /** older = sv ***/
      } END_LOOP_OMP;
      
    }
  
  
  /**  return y*d-dd+0.5*c[0];     **/
  
  Matrix_Vec_mult(old,res,eigen_param,fn);   /**   res = M * old **/
  /**    constant_mult_even(old , debug_val, res) ;  **/
  
  FOREVENFIELDSITES_OMP(i,)
  {
    scalar_mult_su3_vector(res+i, norm2, res+i) ;   /* res = norm2 * M * old   ***/
    sub_su3_vector(res+i, old+i, res+i )  ;              /* res = res - old   **/           
    
    sub_su3_vector(res + i,older + i, res + i)  ; /** res = res - older **/
    scalar_mult_add_su3_vector(res + i, src + i, 0.5 * poly[0], res + i) ; /**  res = res + 0.5*poly[0] * src **/
  } END_LOOP_OMP;
  
  
  /**    dump_vec_even(res) ;   exit(0) ;  **/
  
  /** clean up memory  **/    
  destroy_v_field(old);
  destroy_v_field(tmp);
  destroy_v_field(sv);
  destroy_v_field(older);
  
}




/**
   This routines assumes that the parity is EVEN.
   Apply a Chebyshev polynomial
**/


static void Matrix_Vec_mult_chebyshev_poly(su3_vector *src, su3_vector *res, 
					   ks_eigen_param *eigen_param,
					   imp_ferm_links_t *fn)
{
  ks_eigen_poly *poly_param = &eigen_param->poly;
  int parity = eigen_param->parity;
  Real d1 , d2 , d3 ;
  Real e,c;
  Real sigma,sigma1,sigma_old;
  
  su3_vector *v1 = create_v_field();
  su3_vector *v2 = create_v_field();
  su3_vector  *v3 = create_v_field();
  
  su3_vector *vtmp = create_v_field() ;
  int j , i  ; 
  
  Real a = poly_param->minE  ;
  Real b = poly_param->maxE ;
  int norder = poly_param->norder ;
  
  //    clearvec(v1) ;  /** v1 = 0 **/
  // clearvec(vtmp) ;  /** vtmp = 0 **/
  
  e = (b-a)/2.0;
  c = (b+a)/2.0;
  sigma1 = -e/c;
  
  //   printf("DEBUG a=%f, b=%f ==> e=%f , c=%f \n" , a , b , e , c) ;
  
  //T_0(Q)=1 
  if(norder== 0){
    
    FOREVENFIELDSITES_OMP(i,)
    {
      su3vec_copy(src +i , res + i ) ;                        /* res = src  **/
    } END_LOOP_OMP;
    
    return;
  }
  
  //T_1(Q) = [2/(b-a)*Q - (b+a)/(b-a)]*sigma_1 
  d1 =  sigma1/e;
  d2 =  1.0;
  
  //   d1 =  fabs(sigma1/e) ;  /** debug -- this is s a hack  ***/
  
  
  //      node0_printf("DEBUG (d1,d2) =  (%f,%f) \n", d1 , d2) ;
  
  Matrix_Vec_mult(src,v1,eigen_param,fn);  /**  v1 = M * src **/
  
  // out[s] = d1*v1 + d2*in; 
  FOREVENFIELDSITES_OMP(i,)
  {
    scalar_mult_su3_vector(v1+i, d1, res+i) ;           /* res = d1 * v1   ***/
    scalar_mult_add_su3_vector(res+i, src+i, d2 , res+i) ; /* res = res + d2*src **/
    
  } END_LOOP_OMP;
  
  if(norder== 1){
    return ;
  }
  
  /**   **/
  //degree >=2  
  
  FOREVENFIELDSITES_OMP(i,)
  {
    su3vec_copy(src +i , v1 + i ) ;                        /* v1 = src  **/
    su3vec_copy(res +i , v2 + i ) ;                        /* v2 = res  **/
  } END_LOOP_OMP;
  sigma_old = sigma1;
  
  for (j=2 ; j <=  norder ; ++j)
    { 
      sigma = 1.0/(2.0/sigma1-sigma_old);
      d1 = 2.0*sigma/e;
      d2 = -d1*c;
      d3 = -sigma*sigma_old;
      
      Matrix_Vec_mult(v2,v3,eigen_param,fn);    /** v3 = M * v2  **/
      
      /**  constant_mult_even(v1 , debug_val, v2) ;   ***/
      
      FOREVENFIELDSITES_OMP(i,)
      {
	//  res = d1*v3 + d2*v2 + d3*v1; 
	scalar_mult_su3_vector(v3+i, d1, res+i) ;           /* res = d1  * v3   ***/
	
	scalar_mult_add_su3_vector(res+i, v2+i, d2 , res+i) ; /* res = res + d2*v2 **/
	scalar_mult_add_su3_vector(res+i, v1+i, d3 , res+i) ; /* res = res + d3*v1 **/
	
	su3vec_copy(v2  + i , v1 +i ) ;   /** v1 = v2 ***/
	su3vec_copy(res + i , v2 +i ) ;   /** res = v2 ***/
      } END_LOOP_OMP;
      
      sigma_old  = sigma;
      
    }
  
  
  /**    dump_vec_even(res) ;   exit(0) ;  **/
  
  /** clean up memory  **/    
  destroy_v_field(v1);
  destroy_v_field(v2);
  destroy_v_field(v3);
  destroy_v_field(vtmp);
  
}



/**
   This routines assumes that the parity is EVEN.
   Apply a Chebyshev polynomial
   This is a pure Chebyshev polynomial
**/


static void Matrix_Vec_mult_chebyshev_poly_org(su3_vector *src, su3_vector *res, 
					       ks_eigen_param *eigen_param,
					       imp_ferm_links_t *fn)
{
  ks_eigen_poly *poly_param = &eigen_param->poly;
  int parity = eigen_param->parity;
  Real d1 , d2 , d3 ;
  Real e,c;
  
  su3_vector *v1 = create_v_field();
  su3_vector *v2 = create_v_field();
  su3_vector  *v3 = create_v_field();
  
  su3_vector *vtmp = create_v_field() ;
  int j , i  ; 
  
  Real a = poly_param->minE  ;
  Real b = poly_param->maxE ;
  int norder = poly_param->norder ;
  
  //    clearvec(v1) ;  /** v1 = 0 **/
  // clearvec(vtmp) ;  /** vtmp = 0 **/
  
  e = (b-a)/2.0;
  c = (b+a)/2.0;
  
  
  //T_0(Q)=1 
  if(norder== 0){
    
    FOREVENFIELDSITES_OMP(i,)
    {
      su3vec_copy(src +i , res + i ) ;                        /* res = src  **/
    } END_LOOP_OMP;
    
    return;
  }
  
  //T_1(Q) = [2/(b-a)*Q - (b+a)/(b-a)]*sigma_1 
  d1 =  1.0/e;
  d2 =  1.0;
  
  
  Matrix_Vec_mult(src,v1,eigen_param,fn);  /**  v1 = M * src **/
  
  // out[s] = d1*v1 + d2*in; 
  FOREVENFIELDSITES_OMP(i,)
  {
    scalar_mult_su3_vector(v1+i, d1, res+i) ;           /* res = d1 * v1   ***/
    scalar_mult_add_su3_vector(res+i, src+i, d2 , res+i) ; /* res = res + d2*src **/
    
  } END_LOOP_OMP;
  
  if(norder== 1){
    return ;
  }
  
  /**   **/
  //degree >=2  
  
  FOREVENFIELDSITES_OMP(i,)
  {
    su3vec_copy(src +i , v1 + i ) ;                        /* v1 = src  **/
    su3vec_copy(res +i , v2 + i ) ;                        /* v2 = res  **/
  } END_LOOP_OMP;
  
  for (j=2 ; j <=  norder ; ++j)
    { 
      d1 = 2.0/e;
      d2 = -d1*c;
      d3 = -1.0;
      
      Matrix_Vec_mult(v2,v3,eigen_param,fn);    /** v3 = M * v2  **/
      
      FOREVENFIELDSITES_OMP(i,)
      {
	//  res = d1*v3 + d2*v2 + d3*v1; 
	scalar_mult_su3_vector(v3+i, d1, res+i) ;           /* res = d1  * v3   ***/
	
	scalar_mult_add_su3_vector(res+i, v2+i, d2 , res+i) ; /* res = res + d2*v2 **/
	scalar_mult_add_su3_vector(res+i, v1+i, d3 , res+i) ; /* res = res + d3*v1 **/
	
	su3vec_copy(v2  + i , v1 +i ) ;   /** v1 = v2 ***/
	su3vec_copy(res + i , v2 +i ) ;   /** res = v2 ***/
      } END_LOOP_OMP;
    }
  
  
  /** clean up memory  **/    
  destroy_v_field(v1);
  destroy_v_field(v2);
  destroy_v_field(v3);
  destroy_v_field(vtmp);
  
}








static void Matrix_Vec_mult_sq(su3_vector *src, su3_vector *res, ks_eigen_param *eigen_param,
			       imp_ferm_links_t *fn)
{
  su3_vector *tmp = create_v_field();
  
  Matrix_Vec_mult(src,tmp,eigen_param,fn);  /**  tmp = M * src **/
  Matrix_Vec_mult(tmp,res,eigen_param,fn);   /**   res = M * tmp **/
  
  destroy_v_field(tmp);
}


double chebev(double a, double b, Real  c[], int m, double x)
{
  double d=0.0,dd=0.0,sv,y,y2;
  int j;
  if ((x-a)*(x-b) > 0.0) 
    {
      node0_printf("x not in range (%f,%f) \n", a , b) ;
      exit(1)  ;
    }
  
  /**  y2=2.0*(y=(2.0*x-a-b)/(b-a));   **/
  
  y=(2.0*x-a-b)/(b-a); 
  y2=2.0*y  ;
  
  for (j=m-1;j>=1;j--) 
    { 
      sv=d;
      d=y2*d-dd+c[j];
      dd=sv;
    }
  return y*d-dd+0.5*c[0]; 
  
}


static void cheb_ft (double a, double b , Real  c[] , int n , double (*func)(double, double),
		     double poly_param_1)
{
  const double myPI =  3.1415926535897932385  ;
  double bma ;
  double bpa ;
  double *f ;
  int j , k ;
  double fac ;
  
  f=  (double*)  calloc(n, sizeof(double));
  bma=0.5*(b-a);
  bpa=0.5*(b+a);
  
  for (k=0;k<n;k++) 
    { 
      double y=cos(myPI*(k+0.5)/n); 
      f[k]=(*func)(y*bma+bpa, poly_param_1);
    }
  fac=2.0/n;
  
  for (j=0;j<n;j++) 
    {
      double sum=0.0; 
      for (k=0;k<n;k++) 
	sum +=  f[k]*cos(PI*j*(k+0.5)/n);
      c[j]=fac*sum;
    }
  free(f) ;
  
}

/**
   A simple cut off.
   
**/
double myfunc(double x, double poly_param_1)
{
  double ans ;
  
  if( x < poly_param_1 )
    {
      /**    ans = 1.0 ;   **/
      ans =  100 + 10*x ;
    }
  else
    ans = 0.0 ;
  
  return ans ;
}

/**
   p(x)  = 1 / x
**/

double myfuncA(double x, double poly_param_1)
{
  double ans ;
  
  if( x < poly_param_1 )
    {
      /**    ans = 1.0 ;   **/
      ans =  1/x ;
    }
  else
    ans = 0.0 ;
  
  return ans ;
}


void dump_poloynom(int n, double a, double b, double poly_param_1)
{
  int m = 100 ;
  double del = (b - a) / m ;
  double xx =  a , yy;
  int i ;  
  double del_store = del ; 
  
  del /= 10 ; 
  
  for(i= 0 ; xx <= b ; ++i)
    {
      yy =  chebev(a,b,poly, n, xx) ;
      node0_printf("POLY %f,%f\n", xx , yy) ;
      xx += del ;
      if( xx > 2 )
	{
          del = del_store ; 
	}
      
    }
  
  xx = debug_val ;
  yy =  chebev(a,b,poly, n, xx) ;
  node0_printf("DEBUG_POLY %g %g\n", xx, yy) ;
  
  for(i= 0 ; i < m ; ++i)
    {
      yy =  myfunc(xx, poly_param_1) ;
      node0_printf("FUNC %f,%f\n", xx , yy) ;
      xx += del ;
    }
  
  
}


void dump_poloynomA(int n, double a, double b, double poly_param_1)
{
  int m = 100 ;
  double del = (b - a) / m ;
  double xx =  a , yy;
  int i ;  
  
  for(i= 0 ; i < m ; ++i)
    {
      yy =  chebev(a,b,poly, n, xx) ;
      node0_printf("POLY %f,%f\n", xx , yy) ;
      xx += del ;
    }
  
  xx = debug_val ;
  yy =  chebev(a,b,poly, n, xx) ;
  node0_printf("DEBUG_POLY %g %g\n", xx, yy) ;
  
  for(i= 0 ; i < m ; ++i)
    {
      yy =  myfuncA(xx, poly_param_1) ;
      node0_printf("FUNC %f,%f\n", xx , yy) ;
      xx += del ;
    }
  
  
}




/**
   Reconstruct the eigenvalues from the polynomial
**/


static void my_check_eigen(su3_vector **eigVec, double *eigVal, int Nvecs,
			   int naik_term_epsilon_index, ks_eigen_param *eigen_param,
			   imp_ferm_links_t *fn)
{
  ks_eigen_poly *poly_param = &eigen_param->poly;
  int save_parity = eigen_param->parity;
  int i , ieig ;
  su3_vector *M_eigvec = NULL;
  Real eigVal_est_re ;
  Real eigVal_est_im ;
  Real bot ;
  Real eigen_mod ;
  int k ; 
  
  //      eigen_param->parity = EVENANDODD ;
  eigen_param->parity = EVEN ;

  M_eigvec = create_v_field();
  
  node0_printf("Checking linear eigenvalues\n") ; 
  for(ieig= 0 ; ieig < Nvecs ; ++ieig )
    {
      Matrix_Vec_mult(eigVec[ieig],M_eigvec,eigen_param,fn);  
      
      eigVal_est_re = 0.0 ;
      eigVal_est_im = 0.0 ;
      
      bot = 0.0 ;
      
      
      FOREVENFIELDSITES_OMP(i,reduction(+:bot,eigVal_est_re,eigVal_est_im))
      {
	complex z = su3_dot( eigVec[ieig] + i, M_eigvec +i )  ;
	eigVal_est_re += z.real ; 
	eigVal_est_im += z.imag ; 
	
	bot  += su3_rdot(eigVec[ieig] + i, eigVec[ieig] + i)  ;
      } END_LOOP_OMP;

      g_doublesum(&eigVal_est_re );
      g_doublesum(&eigVal_est_im);
      g_doublesum(&bot);
      
      /**	eigen_mod =  2.0 * pow(eigVal[ieig] , 0.5 )	;  **/
      eigen_mod = 0.0 	;
#if 0
      if( which_poly == POLY_X  || which_poly == POLY_INV  )
	{
	  Real Poly_re = myfunc(eigVal_est_re, poly->poly_param_1)  ;
	  Real a = poly_param->minE  ;
	  Real b = poly_param->maxE ;
	  
	  Real Cheb_re = chebev(a,b,poly, norder, eigVal_est_re) ;  ;
	  
	  node0_printf("EIGVALEST[%d] (%g,%g) norm = %g, Poly(Re_eig) = %g Cheby = %g\n" , 
		       ieig, eigVal_est_re ,eigVal_est_im ,  bot, Poly_re, Cheb_re) ;
	  
	}
      else
	{
	  node0_printf("EIGVALEST[%d] (%g,%g) norm = %g \n" , 
		       ieig, eigVal_est_re ,eigVal_est_im ,  bot) ;
	  
	}
#endif
      
      /**   update the eigenvalue  **/
      eigVal[ieig] =  eigVal_est_re  ;
      
    }  /** end loop over eigenvalues  **/
  
  eigen_param->parity = save_parity ;
  destroy_v_field(M_eigvec); M_eigvec = NULL;
  
}

extern struct
{
  int  logfil, ndigit, mgetv0 ;
  int       msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd;
  int         mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd ;
  int         mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd ;
  
  
  
} debug_ ;





/*****************************************************************************/
int ks_eigensolve_ARPACK(su3_vector **eigVec, double *eigVal, 
			 ks_eigen_param *eigen_param, int init){
  
  int maxnev = eigen_param->Nvecs;    /* number of eigenvalues to compute*/
  Real minE = eigen_param->poly.minE ;
  Real maxE = eigen_param->poly.maxE ;
  int parity = eigen_param->parity;
  ks_eigen_poly *poly_param = &eigen_param->poly; /* Polynomial parameters */
  int which_poly    = poly_param->which_poly;
  Real poly_param_1 = poly_param->poly_param_1;
  Real poly_param_2 = poly_param->poly_param_2;
  int norder        = poly_param->norder;

  int maxn;
  double * evals , *rnorms;    /*work space*/
  double_complex* evecs;    
  
  int i,j,k;
  int ret;
  site* s;
  double *xx;		/* for copying */
  Real *yy;		/* for copying */
  
  /*****  parameters of polynomial  **/
  double a = minE ; 
  double b = maxE ; 
  
  su3_vector tmp1[sites_on_node], tmp2[sites_on_node];
  
  /***  arpack   ***/
  int ido=0;           //control of the action taken by reverse communications
  //set initially to zero
  
  
  char  which[] = "SR" ; 
  //     char  which[] = "LR" ; 
  //"SR" (smallest real part), "LR" (largest real part)
  //"SI" (smallest imaginary part), "LI" (largest imaginary part)
  //"SM" (smallest absolute value), "LM" (largest absolute value)
  
  char bmat[]="I";     /* Specifies that the right hand side matrix
			  should be the identity matrix; this makes
			  the problem a standard eigenvalue problem.
		       */
  
  int nev = eigen_param->Nvecs ;   //number of eigenvectors to be computed (IN) 
  
  double tol = eigen_param->tol  ;   //tolerance for approximating the eigenvalues (IN) 
  int MaxIter = eigen_param->MaxIter ;
  
  int iparam[11] ; 
  int nconv  ;
  
  iparam[0]=1;  //use exact shifts
  //     iparam[0]=0;  //use NO shifts
  iparam[2]= MaxIter ;
  iparam[3]=1;
  iparam[6]=1;
  
  int nsize ; //problem size (IN)
  int info ;
  
  int nkv = eigen_param->nArnoldi  ;    //size of search subspace nkv>nev , nkv<nsize (IN)
  
  int ipntr[14];
  
  dcomplex *resid ; 
  dcomplex *workd ; 
  int lworkl ;
  dcomplex *workl ;
  double *rwork ; 
  dcomplex *workev  ;
  
  
  int rvec=1;       //always call the subroutine that computes orthonormal bais for the eigenvectors
  
  //  int rvec=0;       // eigenvalues only
  char howmany[]="A";     //compute eigenvectors
  
  dcomplex sigma ;
  int *select  ;
  
  dcomplex *eigenVal ; 
  
  int nupdates ;
  int nops ; 
  
  double rrr ; 
  
  //#define PARPACK 1
  
#if defined(MPI_COMMS) && defined(PARPACK)
  MPI_Comm comm; //communicator used when we call PARPACK
  int comm_err ;
#endif
  
  int ierr ;
  
#ifdef PARPACK 
  char tag[] = "P-ARPACK"  ;
#else
  char tag[] = "ARPACK"  ;
#endif
  
  
  /**  debug  ***/
#if 0
  debug_.ndigit = -3 ;
  debug_.logfil = 6 ;
  
  debug_.mcaup2 = 3 ;
  debug_.mcaupd = 3 ;
#endif
  
  /***/
  
#ifdef EIGTIME
  double dtimec;
#endif
  
  imp_ferm_links_t *fn = get_fm_links(fn_links, 0);
  
  if(parity == EVENANDODD){
    maxn=sites_on_node*3;			/*local size of matrix*/
  }
  else
    maxn=sites_on_node*3/2;			/*local size of matrix*/
  
  
  if(   which_poly == NO_POLY_SQRT  ) 
    {
      maxn=sites_on_node*3;			/*local size of matrix*/
      parity = EVENANDODD ;
    }
  
  
  node0_printf("Calculate eigenvalues with %s\n" ,tag ) ;
  node0_printf("%s Number of eigenvales requested = %d\n", tag, maxnev ) ;
  node0_printf("%s size of search subspace (nkv>nev)  = %d\n", tag, nkv ) ;
  node0_printf("%s maximum number of iterations  = %d\n", tag, MaxIter ) ;
  
  
  /** arpack  **/  nsize  = maxn ;
  node0_printf("%s Order of matrix  = %d\n", tag, nsize ) ;
  
  info = 0; //means use a random starting vector with Arnoldi
  
  
  
  /***  arpack memory allocation ****/
  resid = malloc(nsize*sizeof(dcomplex));
  if (resid==NULL) exit(1);
  
  workd = malloc(3*nsize*sizeof(dcomplex));
  if (workd==NULL) exit(1);
  
  lworkl=(3*(nkv)*(nkv)+5*(nkv))*2; //size of some work arrays in arpack
  //  workl = malloc(lworkl*sizeof(dcomplex));
  workl = malloc(10*lworkl*sizeof(dcomplex));
  if (workl==NULL) exit(1);
  
  rwork = malloc(nkv*sizeof(double)) ; 
  if (rwork ==NULL) exit(1);
  
  // workev = malloc(3*nkv*sizeof(dcomplex));
  workev = malloc(30*nkv*sizeof(dcomplex));
  if (workev == NULL) exit(1);
  
  select   = malloc(nkv*sizeof(int));
  if (select == NULL) exit(1);
  
  //eigenVal = malloc(nkv*sizeof(dcomplex));
  eigenVal = malloc(nkv*sizeof(dcomplex));
  if (eigenVal == NULL) exit(1);
  
  evecs = malloc(nsize*nkv*sizeof(dcomplex));
  if (evecs==NULL) exit(1);
  
  
  
  node0_printf("ARPACK eigensolver with polynomial acceleration\n");
  if(   which_poly == NO_POLY ) 
    {
      node0_printf("%s NO POLYNOMIAL ACCELERATION\n" , tag );
    }
  else if(   which_poly == MSQ_POLY ) 
    {
      node0_printf("M**2 POLYNOMIAL ACCELERATION\n");
    }
  else if( which_poly == POLY_X ) 
    {
      node0_printf("CHEBYSHEV POLYNOMIAL INTERPOLATION ACCELERATION used\n");
      node0_printf("Order of polynomial %d\n",norder);
      node0_printf("Cut off with x >  %g, f(x) = 0  \n", poly_param_1);
      node0_printf("Cut off with x <  %g, f(x) = x  \n", poly_param_1);
      node0_printf("Eigenvalue range %g to %g \n",minE , maxE);
      node0_printf("Polynomial parameter [0]  %g \n", poly_param_1);
      node0_printf("Polynomial parameter [1]  %g \n", poly_param_2);
    }
  else if( which_poly == POLY_INV ) 
    {
      node0_printf("CHEBYSHEV POLYNOMIAL INTERPOLATION ACCELERATION used\n");
      node0_printf("Order of polynomial %d\n",norder);
      node0_printf("Cut off with x >  %g, f(x) = 0  \n", poly_param_1);
      node0_printf("Cut off with x <  %g, f(x) = 1/x  \n", poly_param_1);
      node0_printf("Eigenvalue range %g to %g \n",minE , maxE);
      node0_printf("Polynomial parameter [0]  %g \n", poly_param_1);
      node0_printf("Polynomial parameter [1]  %g \n", poly_param_2);
    }
  else if( which_poly == CHEBY ) 
    {
      node0_printf("CHEBYSHEV POLYNOMIAL used for ACCELERATION used\n");
      node0_printf("Order of CHEBYSHEV polynomial %d\n",norder);
      node0_printf("Eigenvalue range %g to %g \n",minE , maxE);
    }
  else if( which_poly == NO_POLY_SQRT  ) 
    {
      node0_printf("(D_operator, sqrt) NO POLYNOMIAL ACCELERATION\n");
    }
  else if( which_poly == CHEBY_ORG  ) 
    {
      node0_printf("CHEBYSHEV ACCELERATION  normal normalization\n");
      node0_printf("Order of CHEBYSHEV polynomial %d\n",norder);
      node0_printf("Eigenvalue range %g to %g \n",minE , maxE);
    }
  else
    {
      node0_printf("ERROR polynomial flag = %d out of range\n",  which_poly  );
      terminate(1);
    }
  
#if 1
  if( poly_param->eigmax == 0 )
    {
      char s_tmp[] = "SR" ;
      strncpy(which,s_tmp,2) ;
      node0_printf("\nMinimum eigenvalues %s\n", which);
    }
  else
    {
      char s_tmp[] = "LR" ;
      strncpy(which,s_tmp,2) ;
      
      node0_printf("\nMaximum eigenvalues %s\n", which);
    }
  
  if( which_poly == NO_POLY_SQRT  )
    {
      char s_tmp[] = "SM" ;
      strncpy(which,s_tmp,2) ;
      
    }
  
#endif
  
  
  /** setup the polynomial    **/
  
  poly = (Real *) calloc(norder, sizeof(double)) ;
  
  if( which_poly == POLY_X )
    {
      cheb_ft(a,b ,poly , norder , myfunc, poly_param_1 ) ;
      dump_poloynom(norder, a,  b, poly_param_1) ;
    }
  else if( which_poly == POLY_INV )
    {
      cheb_ft(a,b ,poly , norder , myfuncA, poly_param_1 ) ;
      dump_poloynomA(norder, a,  b, poly_param_1) ;
    }
  
  
#if defined(MPI_COMMS) && defined(PARPACK)
  /****   parallel stuff  ***********/
  //duplicate the MPI_COMM_WORLD to create a communicator to be used with arpack
  comm_err = MPI_Comm_dup(MPI_COMM_WORLD,&comm); 
  if(comm_err != MPI_SUCCESS) { //error when trying to duplicate the communicator
    node0_printf("MPI_Comm_dup return with an error\n" ) ;
    exit(1);
  }
#endif
  
  /***  arpack-start  ***/
  do
    {
      
#ifdef PARPACK
      
      _AFT(pznaupd)(&comm, &ido,bmat,&nsize,which,&nev,&tol,resid,&nkv,
		    evecs,&nsize,iparam,ipntr,workd, 
		    workl,&lworkl,rwork,&info);
      
      
#else
      _AFT(znaupd)(&ido,bmat,&nsize,which,&nev,&tol,resid,&nkv,
                   evecs,&nsize,iparam,ipntr,workd, 
                   workl,&lworkl,rwork,&info);
      
      
#endif
      
      
      if (ido == 99 || info == 1)
	break;
      
      if ((ido==-1)||(ido==1))
	{
	  
#if 0
	  // unit matrix out = in 
	  xx= &((workd  +ipntr[0] -1 )->real)  ;
	  yy = &((workd  +ipntr[1] -1)->real)  ;
	  FORSOMEPARITY(i,s,parity)
	    {
	      for(k=0;k<6;k++) 
		*(yy++) = *(xx++);
            } END_LOOP;
#endif
	  
#if 1
	  //          printf("debug matrix multiplication called\n");
	  su3_vector *wd = (su3_vector *)&((workd  +ipntr[0] -1 )->real);
	  FORSOMEFIELDPARITY_OMP(i,parity,)
	    {
	      su3vec_copy(wd+i, tmp1+i);
	    } END_LOOP_OMP;
	  
	  /**	  Matrix_Vec_mult(tmp1,tmp2,parity,fn[0]);   **/
	  if( which_poly == NO_POLY  )
	    Matrix_Vec_mult(tmp1,tmp2,eigen_param,fn);  
	  else if ( which_poly == POLY_X  || which_poly == POLY_INV   )
	    Matrix_Vec_mult_poly(tmp1,tmp2,eigen_param,fn)  ;  
	  else if ( which_poly == MSQ_POLY   )
	    Matrix_Vec_mult_sq(tmp1,tmp2,eigen_param,fn)  ;  
	  else if ( which_poly == CHEBY    )
	    Matrix_Vec_mult_chebyshev_poly(tmp1,tmp2,eigen_param,fn)  ;  
	  // 	  else if(   which_poly == NO_POLY_SQRT  ) 
	  // 	    {
	  // 	      Matrix_Vec_mult_sqrt(tmp1,tmp2,parity,fn[0]) ;
	  // 	    }
	  else if(   which_poly == CHEBY_ORG  ) 
	    {
	      Matrix_Vec_mult_chebyshev_poly_org(tmp1,tmp2,eigen_param,fn)  ;  
	    }
	  
	  wd  = (su3_vector *)&((workd  +ipntr[1] -1)->real)  ;
	  FORSOMEFIELDPARITY_OMP(i,parity,)
	    {
	      su3vec_copy(tmp2+i, wd+i);
	    } END_LOOP_OMP;
#endif
	  
	}
      
    } while (ido != 99);
  
  
  if ( info  < 0 ) 
    {
      node0_printf("Error with znaupd, info =%d Check documentation\n" , info) ;
      exit(1) ;
    }
  else 
    {
      nconv  = iparam[4] ;
      
      nupdates = iparam[2];
      nops     = iparam[8]; 
      
      
#ifdef PARPACK
      //       node0_printf("pzneupd called\n") ; 
      _AFT(pzneupd) (&comm, &rvec,howmany,select,eigenVal,evecs,&nsize,&sigma, 
                     workev,bmat,&nsize,which,&nev,&tol,resid,&nkv, 
                     evecs,&nsize,iparam,ipntr,workd,workl,&lworkl, 
                     rwork,&ierr);
      
#else
      node0_printf("zneupd called\n") ; 
      _AFT(zneupd) (&rvec,howmany,select,eigenVal,evecs,&nsize,&sigma, 
		    workev,bmat,&nsize,which,&nev,&tol,resid,&nkv, 
		    evecs,&nsize,iparam,ipntr,workd,workl,&lworkl, 
		    rwork,&ierr);
#endif
      

      /* Sort the eigenvectors according to ascending eigenvalue*/

      /* copy Evalues in global arrays 
       * (convert from double to single precision) */
      for(j=0;j<eigen_param->Nvecs;j++)
	eigVal[j]= eigenVal[j].real ;
      
      int *hash = sort_eigenvalues(eigVal, eigen_param->Nvecs);

      /* copy Evectors in global arrays 
       * (convert from double to single precision if need be) */
      for(int m=0; m<eigen_param->Nvecs; m++) {
	j = hash[m];
	/* Points to the bottom of the jth eigenvector */
	su3_vector *ev = (su3_vector *)(&evecs[0].real + 2*j*maxn);
	FORSOMEFIELDPARITY_OMP(i,parity,){
	  su3vec_copy(ev+i, eigVec[m]+i);
	} END_LOOP_OMP;
      }

      /*report the computed evals and their residuals  **/
      node0_printf("%s number of operations %d\n" , tag, nops);
      node0_printf("ARPACK number of iterations %d (maximum %d) \n" ,
		   iparam[2], MaxIter);
      node0_printf("%s number of converged eigenvalues  %d\n" , tag, nconv);
      node0_printf("Showing eigenvalues, errors and imaginary parts:\n" );
      node0_printf("BEGIN RESULTS\n");
      
      for(k=0; k< nev ; k++)
	{
	  
	  rrr = fabs((workl+ipntr[10]-1+hash[k])->real) ;
	  
	  if( which_poly == NO_POLY_SQRT  )
	    {
	      double eig_sq = eigenVal[hash[k]].real*eigenVal[hash[k]].real ;
	      node0_printf("Eigenvalue(%d)**2 = %g +/- %g Im = %g\n", k, eig_sq, rrr, eigenVal[hash[k]].imag) ; 
	    }
	  else
	    {
	      node0_printf("Eigenvalue(%d) = %g +/- %g Im = %g\n", k, eigenVal[hash[k]].real, rrr, eigenVal[hash[k]].imag) ; 
	    }  
	  
	  //           printf("Eigenvalue(%d) = %g , %g \n", k, eigenVal[krev].real,eigenVal[krev].imag) ; 
	}

      free(hash);
    }  
  

  
#ifdef EIGTIME
  dtimec += dclock();
  node0_printf("KAULKRITER: time = %e iters = %d iters/vec = %e\n",
	       dtimec,total_iters, (double)(total_iters)/eigen_param->Nvecs);
#endif
  

  {
    int naik_term_epsilon_index = 0 ;
    if( which_poly != NO_POLY_SQRT  )
      {
	my_check_eigen(eigVec,eigVal,eigen_param->Nvecs,
		       naik_term_epsilon_index, eigen_param, fn) ;
      }
    
  }

  destroy_fn_links(fn);
  
  /***  free up memory arpack  *************/
  free(resid); 
  free(workd) ;
  free(workl) ;
  free(rwork) ;
  free(workev) ;
  free(select) ;
  free(eigenVal) ;
  
  free(evecs)  ; 
  
  
  cleanup_Matrix() ;
  return iparam[2] ;  /*****  update ************/
}
