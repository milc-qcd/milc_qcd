/*
  Code for making rational functions for MILC ks_imp_rhmc code
  Version 1   5/10/06 Doug Toussaint
  Based on "test" code by Mike Clark, using Mike Clarks alg_remez.C
TEST - generate sqrt of usual functions, use two spinors

  Usage: poly4 nflavors1 mass1 nflavors2 mass2 nflavors3 mass3 nflavors4 mass4 \
     moldyn_order action_order eval_min eval_max precision
  Example: poly4 2 .0036 -2 .018 0 1 0 1 4 5 1e-15 90 65
  
  Calculate three rational approximations, for molecular dynamics
  evolution, gaussian random heatbath update, and fermion action computation.

  For the moment, specify four flavor numbers and masses (flavor number
    can be negative )

  Molecular dynamics:
	PROD  X^(-nf_i/8) where X = M^dagger M
  Heatbath ("GR")
	PROD  X^(+nf_i/16) where X = M^dagger M
  Fermion Acttion:
	PROD  X^(-nf_i/16) where X = M^dagger M

*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

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

  // Set the exponents and masses
  sscanf(argv[++i],"%d",&y1);
  sscanf(argv[++i],"%le",&m1);
  sscanf(argv[++i],"%d",&y2);
  sscanf(argv[++i],"%le",&m2);
  sscanf(argv[++i],"%d",&y3);
  sscanf(argv[++i],"%le",&m3);
  sscanf(argv[++i],"%d",&y4);
  sscanf(argv[++i],"%le",&m4);

  // Set the required degrees of approximation
  sscanf(argv[++i],"%d",&order1);
  sscanf(argv[++i],"%d",&order2);

  // Set the approximation bounds
  sscanf(argv[++i],"%le",&lambda_low);
  sscanf(argv[++i],"%le",&lambda_high);

  // Set the precision of the arithmetic
  sscanf(argv[++i],"%d",&precision);

  printf("#define MOLDYN_ORDER_X %d // order of approximation for updating factor X\n",order1);
  printf("#define GRSOURCE_ORDER_X %d // order for heat bath pseudofermion factor X\n",order2);
  printf("#define ACTION_ORDER_X %d // order for factorting fermion action factor X\n",order2);

// 4->8 and 8>16 IN THE NEXT TWO LINES (TEST)
  work( y1,8,m1,  y2,8,m2,  y3,8,m3,  y4,8,m4, order1, lambda_low, lambda_high, precision,
    "OMIT", "MD" );
  work( y1,16,m1,  y2,16,m2,  y3,16,m3,  y4,16,m4, order2, lambda_low, lambda_high, precision,
    "GR", "FA" );
}

int work( int y1,int z1, double m1,  int y2,int z2, double m2,  int y3,int z3, double m3,
    int y4,int z4, double m4, int order, double lambda_low,  double lambda_high, int precision,
    char *tag1, char *tag2 ){
  // The error from the approximation (the relative error is minimised
  // - if another error minimisation is requried, then line 398 in
  // alg_remez.C is where to change it)
  double error;

  // The partial fraction expansion takes the form 
  // r(x) = norm + sum_{k=1}^{n} res[k] / (x + pole[k])
  double norm;
  double *res = new double[order];
  double *pole = new double[order];

  double bulk = exp(0.5*(log(lambda_low)+log(lambda_high)));

  // Instantiate the Remez class
  AlgRemez remez(lambda_low,lambda_high,precision);

  // Generate the required approximation
  error = remez.generateApprox(order,order,y1,z1,m1,y2,z2,m2,y3,z3,m3,y4,z4,m4);

  // Find the partial fraction expansion of the approximation 
  // to the function x^{y/z} (this only works currently for 
  // the special case that n = d)
  remez.getPFE(res,pole,&norm);
  
  if( strcmp(tag1,"OMIT") != 0 ){
    printf("\n//Approximation to f(x) = (x+4*%f^2)^(%d/%d) (x+4*%f^2)^(%d/%d)\n// (x+4*%f^2)^(%d/%d) (x+4*%f^2)^(%d/%d)",
  	m1, y1, z1, m2, y2, z2, m3, y3, z3, m4, y4, z4);
    printf("//where x = M^dagger M\n",m1);
    printf("//residues:\nReal A_%s_X[%d] = {\n",tag1,order+1);
    printf("%18.16e,\n", norm);
    for (int i = 0; i < order; i++) {
      printf("%18.16e", res[i]);
      if(i<order-1)printf(",\n" ); 
      else printf("\n};\n");
    }
    printf("//(-)roots\nReal B_%s_X[%d] = {\n",tag1,order+1);
    printf("99.9,  //DUMMY\n"); //DUMMY!!
    for (int i = 0; i < order; i++) {
      printf("%18.16e", pole[i] ); 
      if(i<order-1)printf(",\n" ); 
      else printf("\n};\n");
    }
    // Check - compute the function at the endpoints of the interval
    double x,sum;
    int ii;
    for ( x = lambda_low, sum=norm, ii = 0; ii < order; ii++) {
       sum += res[ii]/(x+pole[ii]);
    }
    printf("//CHECK: f(%e) = %e\n",x,sum);
  }


  // Find pfe of inverse function
  remez.getIPFE(res,pole,&norm);

  printf("\n//Approximation to f(x) = (x+4*%f^2)^(-%d/%d) (x+4*%f^2)^(-%d/%d) (x+4*%f^2)^(-%d/%d) (x+4*%f^2)^(-%d/%d)\n\n", m1, y1, z1, m2, y2, z2, m3, y3, z3, m4, y4, z4);
  printf("//where x = M^dagger M\n",m1);
  printf("//residues:\nReal A_%s_X[%d] = {\n",tag2,order+1);
  printf("%18.16e,\n", norm);
  for (int i = 0; i < order; i++) {
    printf("%18.16e", res[i]);
    if(i<order-1)printf(",\n" ); 
    else printf("\n};\n");
  }
  printf("//(-)roots\nReal B_%s_X[%d] = {\n",tag2,order+1);
  printf("99.9, //DUMMY\n"); //DUMMY!!
  for (int i = 0; i < order; i++) {
    printf("%18.16e", pole[i] ); 
    if(i<order-1)printf(",\n" ); 
    else printf("\n};\n");
  }

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
