/*
 *  This file contains a number of routines to calculate
 *  the smearing wave functions.
 *
 *
 */

#include<stdio.h>
#include<math.h>
#include<string.h>

/*** functions required on this file ********/
Real ran1(int *idum) ;
Real periodic_radius(int x, int y, int z, int nx);
Real periodic_radius_x_squeeze(int x, int y, int z, int nx, Real squeeze) ;

int where(int x, int y, int z ,int nx) ;

/*
 *  Create the local smearing operator.
 *
 *   src(x) = 0   if   x != 0
 *   src(x) = 1   if   x =  0
 *
 */


void create_local_oper(Real *src, int dim, char filename[])
{
  int i ;

  for(i=0 ; i < dim ; ++i)
    *(src + i) = 0.0 ;

  *src = 1.0 ;

  sprintf(filename,"local_smear");
  printf("I have created a LOCAL source \n");

}




/*
 *  Create the wall smearing operator.
 *
 *   src(x) = 1   all x 
 *
 */


void create_wall_oper(Real *src, int dim, char filename[])
{
  int i ;

  for(i=0 ; i < dim ; ++i)
    *(src + i) = 1.0 ;

  printf("I have created a WALL source \n");
  sprintf(filename,"wall_smear");


}




/*
 *  Create a random smearing operator.
 *  Useful for debuggin purposes.
 *
 *   src(x) = [0 .. 1 ]   all x 
 *
 */


void create_random_oper(Real *src, int *seed, int dim, char filename[])
{
  static int rem = 0 ;
  int i ;

  for(i=0 ; i < dim ; ++i)
    *(src + i) = ran1(seed) ;

  printf("I have created a RANDOM source \n");
  sprintf(filename,"random_smear_%d",rem);

  ++rem ;

}



/*
 *  Create an exponential smearing function
 *
 *   src(x) = exp( -alpha * x ) 
 *
 */


void create_expon_oper(Real *src, int nx, char filename[])
{
  int x,y,z ;
  int pt ;
  Real rad ;
  Real src_local ;
  Real decay ;

  /***  input the required source parameters *****/
  printf("Exponential source \n");
  printf("src(r) = exp(-alpha*r) \nInput the alpha parameter\n");
  scanf("%lfHELP",&decay);

  /*** create the source ********/
  
  for(x = 0 ; x < nx ; ++x)
    for(y = 0 ; y < nx ; ++y)
      for(z = 0 ; z < nx ; ++z)
      {

	rad =  periodic_radius(x, y, z, nx) ;

	src_local = (Real) exp( (double)  -decay*rad );

	pt = where(x,y,z,nx) ;
	*(src + pt) = src_local ;
      }


  printf("I have created a source:: src(x) = exp(-%g * x) \n",decay);
  sprintf(filename,"exp_smear_%g",decay);


}





/*
 *  Create an exponential smearing function
 *  with the x direction being squeezed
 *
 *   src(x) = exp( -alpha * r ) 
 *
 *  r = sqrt( (sq*x)**2 + y**2 + z**2 )
 *
 */


void create_squeeze_expon_oper(Real *src, int nx, char filename[])
{
  int x,y,z ;
  int pt ;
  Real rad ;
  Real src_local ;
  Real decay ;
  Real x_squeeze ; 

  /***  input the required source parameters *****/
  printf("Exponential source SQUEEZED in the x-direction\n");

  printf("src(r) = exp(-alpha*r) \nInput the alpha parameter\n");
  scanf("%lfHELP",&decay);

  printf("Input the x squeezing factor \n") ; 
  scanf("%lfHELP",&x_squeeze) ; 

  /*** create the source ********/
  
  for(x = 0 ; x < nx ; ++x)
    for(y = 0 ; y < nx ; ++y)
      for(z = 0 ; z < nx ; ++z)
      {

	rad =  periodic_radius_x_squeeze(x, y, z, nx, x_squeeze) ;

	src_local = (Real) exp( (double)  -decay*rad );

	pt = where(x,y,z,nx) ;
	*(src + pt) = src_local ;
      }


  printf("I have created a source:: src(r) = exp(-%g * r) with x-> %f x\n",decay,x_squeeze);
  sprintf(filename,"exp_smear_%g_squeeze_%g",decay,x_squeeze);


}




/*
 *  Create a gaussian smearing function
 *
 *   src(x) = exp( - alpha * x**2 ) 
 *
 */


void create_gaussian_oper(Real *src, int nx, char filename[])
{
  int x,y,z ;
  int pt ;
  Real rad ;
  Real src_local ;
  Real decay ;

  /***  input the required source parameters *****/
  printf("Gaussian source \n");
  printf("src(r) = exp(-alpha*r*r) \nInput the alpha parameter\n");
  scanf("%lfHELP",&decay);

  /*** create the source ********/
  
  for(x = 0 ; x < nx ; ++x)
    for(y = 0 ; y < nx ; ++y)
      for(z = 0 ; z < nx ; ++z)
      {

	rad =  periodic_radius(x, y, z, nx) ;

	src_local = (Real) exp( (double)  -decay*rad*rad );

	pt = where(x,y,z,nx) ;
	*(src + pt) = src_local ;
      }


  printf("I have created a GAUSSIAN source:: src(x) = exp(-%g * x*x) \n",decay);
  sprintf(filename,"gaussian_smear_%g",decay);


}







/*
 *  Create a gaussian smearing function
 *
 *   src(x) = exp( - x**2/(sigma**2)  ) 
 *
 */


void create_milc_gaussian_oper(Real *src, int nx, char filename[])
{
  int x,y,z ;
  int pt ;
  Real rad ;
  Real src_local ;
  Real decay ;

  /***  input the required source parameters *****/
  printf("MILC Gaussian source \n");
  printf("src(r) = exp(-r*r/(sigma*sigma) ) \nInput the sigma parameter\n");
  scanf("%lfHELP",&decay);

  /*** create the source ********/
  
  for(x = 0 ; x < nx ; ++x)
    for(y = 0 ; y < nx ; ++y)
      for(z = 0 ; z < nx ; ++z)
      {

	rad =  periodic_radius(x, y, z, nx) ;

	src_local = (Real) exp( (double)  -(rad*rad)/(decay*decay)  );

	pt = where(x,y,z,nx) ;
	*(src + pt) = src_local ;
      }


  printf("I have created a GAUSSIAN source:: src(x) = exp(- x*x/(%g)**2) \n",decay);
  sprintf(filename,"gaussian_milc_smear_%g",decay);


}


/*
 *  Create an exponential smearing function
 *  Hydrogenic smearing function
 *
 *   src(x) = exp( -alpha * x /n )  * polynomial (2 *alpha *x/n)
 * 
 *   The polynomial is an asscoaited laguerre ploynomial
 *
 *   Source normalization:   src(0) = 1
 */


void create_hydrogen(Real *src, int quantum_n, Real (*func)(Real),
int nx, char filename[])
{
  int x,y,z ;
  int pt ;
  Real rad ;
  Real src_local ;
  Real nrm = (*func)(0.0);  /** normalization factor at the origin ***/
  Real decay ;
  /***----------------------------------------********/

  if( quantum_n <= 1 )
  {
    printf("quantum n factor %d is out of range",quantum_n );
    exit(2);
  }

  /***  input the required source parameters *****/
  printf("Hydrogen wave function source with quantum number n = %d \n",quantum_n);
  printf("src(x) = exp( -alpha * x /%d )  * polynomial_[ %d ] (2 *alpha *x/%d)\n",quantum_n,quantum_n,quantum_n);

  printf("Input the alpha parameter  \n");
  scanf("%lfHELP",&decay);

  /*** create the source ********/

  sprintf(filename,"hydro_smear_%g_S%d",decay,quantum_n );
  decay /= (Real) quantum_n  ;


  for(x = 0 ; x < nx ; ++x)
    for(y = 0 ; y < nx ; ++y)
      for(z = 0 ; z < nx ; ++z)
      {
	rad =  periodic_radius(x, y, z, nx) ;

	src_local = (Real) exp( (double)  -decay*rad );
	rad *= 2.0*decay ;
	  
	src_local *= (*func)(rad);
	src_local /= nrm ;

	pt = where(x,y,z,nx) ;
	*(src + pt) = src_local ;

      }


  printf("I have created a Hydrogen smearing function source: decay = %g n =  %d S\n",decay,quantum_n);



}


/*  2S model smearing function
 *
 *  Create an exponential times x smearing function
 *
 *   src(x) = exp( -alpha * x ) * (1 - A*x)
 *
 */


void create_2S_polyexpon_oper(Real *src, int nx, char filename[])
{
  int x,y,z ;
  int pt ;
  Real rad ;
  Real src_local ;
  Real decay ;
  Real A ;

  /***  input the required source parameters *****/
  printf("2S Exponential times polynomial source \n");
  printf("src(r) = (1 - A*x) * exp(-alpha*r) \n");

  printf("Input the alpha parameter \n");
  scanf("%lfHELP",&decay);

  printf("Input the A parameter \n");
  scanf("%lfHELP",&A);

  /*** create the source ********/
  
  for(x = 0 ; x < nx ; ++x)
    for(y = 0 ; y < nx ; ++y)
      for(z = 0 ; z < nx ; ++z)
      {

	rad =  periodic_radius(x, y, z, nx) ;

	src_local = (Real) (1.0 - A*rad)*exp( (double)  -decay*rad );

	pt = where(x,y,z,nx) ;
	*(src + pt) = src_local ;
      }


  printf("I have created a source:: src(x) = ( 1  - %g *x ) * exp(-%g * x) \n",A,decay);
  sprintf(filename,"polyexp_smear_%g_x%g",decay,A);


}




/*  2S model smearing function
 *
 *  Create an exponential times x smearing function
 *
 *   src(x) = exp( -alpha * r ) * (1 - A*r)
 *
 *  r = sqrt( (sq*x)**2 + y**2 + z**2 )
 *
 */


void create_2S_polyexpon_oper_x_squeeze(Real *src, int nx, char filename[])
{
  int x,y,z ;
  int pt ;
  Real rad ;
  Real src_local ;
  Real decay ;
  Real A ;
  Real x_squeeze ; 

  /***  input the required source parameters *****/
  printf("2S Exponential times polynomial source SQUEEZED in the x-direction\n");
  printf("src(r) = (1 - A*x) * exp(-alpha*r) \n");

  printf("Input the alpha parameter \n");
  scanf("%lfHELP",&decay);

  printf("Input the A parameter \n");
  scanf("%lfHELP",&A);

  printf("Input the x squeezing factor \n") ; 
  scanf("%lfHELP",&x_squeeze) ; 

  /*** create the source ********/
  
  for(x = 0 ; x < nx ; ++x)
    for(y = 0 ; y < nx ; ++y)
      for(z = 0 ; z < nx ; ++z)
      {
	rad =  periodic_radius_x_squeeze(x, y, z, nx, x_squeeze) ;

	src_local = (Real) (1.0 - A*rad)*exp( (double)  -decay*rad );

	pt = where(x,y,z,nx) ;
	*(src + pt) = src_local ;
      }


  printf("I have created a source:: src(r) = ( 1  - %g *r ) * exp(-%g * r) with x-> %f x\n",
	 A,decay,x_squeeze);
  sprintf(filename,"polyexp_smear_%g_x%g_squeeze_%g",decay,A,x_squeeze);


}



/*  3S model smearing function
 *
 *  Create an exponential times a polynomial smearing function
 *
 *   src(x) = exp( -alpha * x ) * (1 - B*x - D*x**2)
 *
 */


void create_3S_polyexpon_oper(Real *src, int nx, char filename[])
{
  int x,y,z ;
  int pt ;
  Real rad ;
  Real src_local ;
  Real decay ;
  Real B,D ;

  /***  input the required source parameters *****/
  printf("3S Exponential times polynomial source \n");
  printf("src(r) = (1 - B*x -D*x**2 ) * exp(-alpha*r) \n");

  printf("Input the alpha parameter \n");
  scanf("%lfHELP",&decay);

  printf("Input the B parameter \n");
  scanf("%lfHELP",&B);

  printf("Input the D parameter \n");
  scanf("%lfHELP",&D);

  /*** create the source ********/
  
  for(x = 0 ; x < nx ; ++x)
    for(y = 0 ; y < nx ; ++y)
      for(z = 0 ; z < nx ; ++z)
      {

	rad =  periodic_radius(x, y, z, nx) ;

	src_local = (Real) (1.0 - B*rad -D*rad*rad )*exp( (double)  -decay*rad );

	pt = where(x,y,z,nx) ;
	*(src + pt) = src_local ;
      }


  printf("I have created a source:: src(x) = (1 - %g*x -%g*x**2) * exp(-%g * x) \n",B,D,decay);
  sprintf(filename,"polyexp_smear_%g_x%g_xsq%g",decay,B,D);


}




/*  3S model smearing function. SQUEEZED in the X diraction.
 *
 *  Create an exponential times a polynomial smearing function
 *
 *   src(r) = exp( -alpha * r ) * (1 - B*r - D*r**2)
 *
 *  r = sqrt( (sq*x)**2 + y**2 + z**2 )
 */


void create_3S_polyexpon_oper_x_squeeze(Real *src, int nx, char filename[])
{
  int x,y,z ;
  int pt ;
  Real rad ;
  Real src_local ;
  Real decay ;
  Real B,D ;
  Real x_squeeze ; 

  /***  input the required source parameters *****/
  printf("3S Exponential times polynomial source \n");
  printf("src(r) = (1 - B*x -D*x**2 ) * exp(-alpha*r) \n");

  printf("Input the alpha parameter \n");
  scanf("%lfHELP",&decay);

  printf("Input the B parameter \n");
  scanf("%lfHELP",&B);

  printf("Input the D parameter \n");
  scanf("%lfHELP",&D);

  printf("Input the x squeezing factor \n") ; 
  scanf("%lfHELP",&x_squeeze) ; 

  /*** create the source ********/
  
  for(x = 0 ; x < nx ; ++x)
    for(y = 0 ; y < nx ; ++y)
      for(z = 0 ; z < nx ; ++z)
      {

	rad =  periodic_radius_x_squeeze(x, y, z, nx, x_squeeze) ;

	src_local = (Real) (1.0 - B*rad -D*rad*rad )*exp( (double)  -decay*rad );

	pt = where(x,y,z,nx) ;
	*(src + pt) = src_local ;
      }


  printf("I have created a source:: src(r) = (1 - %g*r -%g*r**2) * exp(-%g * r) with x-> %f x\n",
	 B,D,decay,x_squeeze);
  sprintf(filename,"polyexp_smear_%g_x%g_xsq%g_squeeze_%g",decay,B,D,x_squeeze);


}


