/*
  Code for making rational functions for MILC ks_imp_rhmc code
  Version 1   5/10/06 Doug Toussaint
  Based on "test" code by Mike Clark, using Mike Clarks alg_remez.C

  Usage: poly4 < [params_in_file] > [params_rhmc_file]

  where the first line in params_in_file contains the number of lines,
  one per pseudofermion field and each subsequent line contains 13
  values as follows:

  nflavors1 mass1 nflavors2 mass2 nflavors3 mass3 nflavors4 mass4 \
     moldyn_order action_order eval_min eval_max precision

  Example: 2 .0036 -2 .018 0 1 0 1 4 5 1e-15 90 65
  
  Calculate three rational approximations, for molecular dynamics
  evolution, gaussian random heatbath update, and fermion action computation.

  For the moment, specify four flavor numbers and masses (flavor number
    can be negative )

  Molecular dynamics:
	PROD  X^(-nf_i/4) where X = M^dagger M
  Heatbath ("GR")
	PROD  X^(+nf_i/8) where X = M^dagger M
  Fermion Acttion:
	PROD  X^(-nf_i/8) where X = M^dagger M

*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>

#include"alg_remez.h"
int work( int y1,int z1, double m1,  int y2,int z2, double m2,  int y3,int z3, double m3,
    int y4,int z4, double m4, int order, double lambda_low,  double lambda_high, int precision,
    char *tag1, char *tag2 );

int main (int argc, char* argv[]) {

  int i=0;
  int order1; // The degree of the numerator and denominator polynomials for molecular dynamics
  int order2; // The degree of the numerator and denominator polynomials for grsource and
		// fermion action
  int y1; // The numerator of the exponent for flavor 1
  double m1; // The mass for flavor 1
  int y2; // The numerator of the exponent for flavor 2
  double m2; // The mass for flavor 2
  int y3; // The numerator of the exponent for flavor 3
  double m3; // The mass for flavor 3
  int y4; // The numerator of the exponent for flavor 4
  double m4; // The mass for flavor 4
  int precision; // The precision that gmp uses
  double lambda_low, lambda_high; // The bounds of the approximation
  int iphi,n_pseudo;
  double naik_term_epsilon;

  // Read the number of pseudofermion fields
  scanf("%d",&n_pseudo);
  printf("n_pseudo %d\n\n",n_pseudo);

  for(iphi=0; iphi<n_pseudo; iphi++){
    // Set the exponents and masses
    fprintf(stderr,"For pseudofermion %d\n",iphi);
    scanf("%le",&naik_term_epsilon);
    printf("naik_term_epsilon %g\n",naik_term_epsilon);
    scanf("%d",&y1);
    scanf("%le",&m1);
    scanf("%d",&y2);
    scanf("%le",&m2);
    scanf("%d",&y3);
    scanf("%le",&m3);
    scanf("%d",&y4);
    scanf("%le",&m4);
    
    // Set the required degrees of approximation
    scanf("%d",&order1);
    scanf("%d",&order2);
    
    // Set the approximation bounds
    scanf("%le",&lambda_low);
    scanf("%le",&lambda_high);
    
    // Set the precision of the arithmetic
    scanf("%d",&precision);

    // For the MD term we need only the inverse
    work( y1,4,m1,  y2,4,m2,  y3,4,m3,  y4,4,m4, order1, 
	  lambda_low, lambda_high, precision, "OMIT", "MD" );
    // The random source term takes the function and action term,
    // the inverse
    work( y1,8,m1,  y2,8,m2,  y3,8,m3,  y4,8,m4, order2, 
	  lambda_low, lambda_high, precision, "GR", "FA" );
  }
}

void print_check_ratfunc(int y1,int y2,int y3,int y4,
			 int z1,int z2,int z3,int z4,
			 double m1,double m2,double m3,double m4,
			 double lambda_low, int order,
			 double norm,  double *res, double *pole, char *tag){

  if( strcmp(tag,"OMIT") == 0 )return;

  printf("\n\n# Rational function for %s\n",tag);
  printf("y_%s %d %d %d %d\n", tag, y1, y2, y3, y4);
  printf("z_%s %d %d %d %d\n", tag, z1, z2, z3, z4);
  printf("m_%s %f %f %f %f\n", tag, m1, m2, m3, m4);
  printf("order_%s %d\n",tag,order);
  printf("\n");
  printf("res_%s %18.16e\n", tag, norm);
  for (int i = 0; i < order; i++)
    printf("res_%s %18.16e\n", tag, res[i]);
  printf("\n");
  printf("pole_%s 99.9\n", tag); //DUMMY!!
  for (int i = 0; i < order; i++) {
    printf("pole_%s %18.16e\n", tag, pole[i] ); 
  }
  printf("\n");

  // Check - compute the function at the low endpoint of the interval
  double x,sum;
  int ii;
  for ( x = lambda_low, sum=norm, ii = 0; ii < order; ii++) {
    sum += res[ii]/(x+pole[ii]);
  }
  
  double f1 = pow(x+4*m1*m1,((double)y1)/z1);
  double f2 = pow(x+4*m2*m2,((double)y2)/z2);
  double f3 = pow(x+4*m3*m3,((double)y3)/z3);
  double f4 = pow(x+4*m4*m4,((double)y4)/z4);

  printf("# CHECK: f(%e) = %e = %e?\n",x,sum,f1*f2*f3*f4);
}

int work(int y1,int z1, double m1,  
	 int y2,int z2, double m2,  
	 int y3,int z3, double m3,
	 int y4,int z4, double m4, 
	 int order, 
	 double lambda_low,  double lambda_high, 
	 int precision,
	 char *tag1, char *tag2 ){
  // The error from the approximation (the relative error is minimised
  // - if another error minimisation is requried, then line 398 in
  // alg_remez.C is where to change it)
  double error;
  
  double *res = new double[order];
  double *pole = new double[order];
  
  // The partial fraction expansion takes the form 
  // r(x) = norm + sum_{k=1}^{n} res[k] / (x + pole[k])
  double norm;

  double bulk = exp(0.5*(log(lambda_low)+log(lambda_high)));

  // Instantiate the Remez class,
  AlgRemez remez(lambda_low,lambda_high,precision);

  // Generate the required approximation
  fprintf(stderr,
	  "Generating a (%d,%d) rational function using %d digit precision.\n",
	  order,order,precision);
  error = remez.generateApprox(order,order,y1,z1,m1,y2,z2,m2,
			       y3,z3,m3,y4,z4,m4);

  // Find the partial fraction expansion of the approximation 
  // to the function x^{y/z} (this only works currently for 
  // the special case that n = d)
  remez.getPFE(res,pole,&norm);
  
  print_check_ratfunc(y1,y2,y3,y4,z1,z2,z3,z4,m1,m2,m3,m4,
		      lambda_low,order,norm,res,pole,tag1);

  // Find pfe of the inverse function
  remez.getIPFE(res,pole,&norm);

  print_check_ratfunc(-y1,-y2,-y3,-y4,z1,z2,z3,z4,m1,m2,m3,m4,
		      lambda_low,order,norm,res,pole,tag2);

  FILE *error_file = fopen("error.dat", "w");
  for (double x=lambda_low; x<lambda_high; x*=1.01) {
    double f = remez.evaluateFunc(x);
    double r = remez.evaluateApprox(x);
    fprintf(error_file,"%e %e\n", x,  (r - f)/f);
  }
  fclose(error_file);

  delete res;
  delete pole;

}
