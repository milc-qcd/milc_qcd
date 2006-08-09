/*
  Mike Clark - 25th May 2005

  Quick test code for calculating optimal rational approximations for
  the functions x^(y/z) with appropriate bounds.

*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include"alg_remez.h"

int main (int argc, char* argv[]) {

  int i=0;
  int n; // The degree of the numerator polynomial
  int d; // The degree of the denominator polynomial
  int y1; // The numerator of the exponent for flavor 1
  int z1; // The denominator of the exponent for flavor 1
  double m1; // The mass for flavor 1
  int y2; // The numerator of the exponent for flavor 2
  int z2; // The denominator of the exponent for flavor 2
  double m2; // The mass for flavor 2
  int y3; // The numerator of the exponent for flavor 3
  int z3; // The denominator of the exponent for flavor 3
  double m3; // The mass for flavor 3
  int y4; // The numerator of the exponent for flavor 4
  int z4; // The denominator of the exponent for flavor 4
  double m4; // The mass for flavor 4
  int precision; // The precision that gmp uses
  double lambda_low, lambda_high; // The bounds of the approximation

  // Set the exponent
  sscanf(argv[++i],"%d",&y1);
  sscanf(argv[++i],"%d",&z1);
  sscanf(argv[++i],"%le",&m1);
  sscanf(argv[++i],"%d",&y2);
  sscanf(argv[++i],"%d",&z2);
  sscanf(argv[++i],"%le",&m2);
  sscanf(argv[++i],"%d",&y3);
  sscanf(argv[++i],"%d",&z3);
  sscanf(argv[++i],"%le",&m3);
  sscanf(argv[++i],"%d",&y4);
  sscanf(argv[++i],"%d",&z4);
  sscanf(argv[++i],"%le",&m4);

  // Set the required degree of approximation
  sscanf(argv[++i],"%d",&n);
  sscanf(argv[++i],"%d",&d);

  // Set the approximation bounds
  sscanf(argv[++i],"%le",&lambda_low);
  sscanf(argv[++i],"%le",&lambda_high);

  // Set the precision of the arithmetic
  sscanf(argv[++i],"%d",&precision);

  // The error from the approximation (the relative error is minimised
  // - if another error minimisation is requried, then line 398 in
  // alg_remez.C is where to change it)
  double error;

  // The partial fraction expansion takes the form 
  // r(x) = norm + sum_{k=1}^{n} res[k] / (x + pole[k])
  double norm;
  double *res = new double[n];
  double *pole = new double[d];

  double bulk = exp(0.5*(log(lambda_low)+log(lambda_high)));

  // Instantiate the Remez class
  AlgRemez remez(lambda_low,lambda_high,precision);

  // Generate the required approximation
  error = remez.generateApprox(n,d,y1,z1,m1,y2,z2,m2,y3,z3,m3,y4,z4,m4);

  FILE *output = fopen("approx.dat", "w");

  fprintf(output, "Approximation to f(x) = (x+4*%f^2)^(%d/%d) (x+4*%f^2)^(%d/%d) (x+4*%f^2)^(%d/%d) (x+4*%f^2)^(%d/%d)\n\n", m1, y1, z1, m2, y2, z2, m3, y3, z3, m4, y4, z4);

  // Find the partial fraction expansion of the approximation 
  // to the function x^{y/z} (this only works currently for 
  // the special case that n = d)
  remez.getPFE(res,pole,&norm);
  
  fprintf(output, "alpha[0] = %18.16e\n", norm);
  for (int i = 0; i < n; i++) {
    fprintf(output, "alpha[%d] = %18.16e, beta[%d] = %18.16e\n", 
	    i+1, res[i], i+1, pole[i]);
  }
  printf("\nApproximation to f(x) = (x+4*%f^2)^(%d/%d) (x+4*%f^2)^(%d/%d)\n (x+4*%f^2)^(%d/%d) (x+4*%f^2)^(%d/%d)",
	m1, y1, z1, m2, y2, z2, m3, y3, z3, m4, y4, z4);
  printf("where x = M^dagger M\n",m1);
  printf("RESIDUES\n");
  printf("%18.16e,\n", norm);
  for (int i = 0; i < n; i++) {
    printf("%18.16e,\n", res[i]);
  }
  printf("ROOTS\n");
  printf("99.9,  //DUMMY\n"); //DUMMY!!
  for (int i = 0; i < n; i++) {
    printf("%18.16e,\n", pole[i] ); 
  }
  // Check - compute the function at the endpoints of the interval
  double x,sum;
  int ii;
  for ( x = lambda_low, sum=norm, ii = 0; ii < n; ii++) {
     sum += res[ii]/(x+pole[ii]);
  }
  printf("CHECK: f(%e) = %e\n",x,sum);


  // Find pfe of inverse function
  remez.getIPFE(res,pole,&norm);
  fprintf(output, "\nApproximation to f(x) = (x+4*%f^2)^(-%d/%d) (x+4*%f^2)^(-%d/%d) (x+4*%f^2)^(-%d/%d) (x+4*%f^2)^(-%d/%d)\n\n", m1, y1, z1, m2, y2, z2, m1, y1, z1, m2, y2, z2, m3, y3, z3, m4, y4, z4);
  fprintf(output, "alpha[0] = %18.16e\n", norm);
  for (int i = 0; i < n; i++) {
    fprintf(output, "alpha[%d] = %18.16e, beta[%d] = %18.16e,\n", 
	    i+1, res[i], i+1, pole[i]);
  }
  printf("\nApproximation to f(x) = (x+4*%f^2)^(-%d/%d) (x+4*%f^2)^(-%d/%d) (x+4*%f^2)^(-%d/%d) (x+4*%f^2)^(-%d/%d)\n\n", m1, y1, z1, m2, y2, z2, m3, y3, z3, m4, y4, z4);
  printf("where x = M^dagger M\n",m1);
  printf("RESIDUES\n");
  printf("%18.16e,\n", norm);
  for (int i = 0; i < n; i++) {
    printf("%18.16e,\n", res[i]);
  }
  printf("ROOTS\n");
  printf("99.9, //DUMMY\n"); //DUMMY!!
  for (int i = 0; i < n; i++) {
    printf("%18.16e,\n", pole[i] ); 
  }

  fclose(output);

  FILE *error_file = fopen("error.dat", "w");
  for (double x=lambda_low; x<lambda_high; x*=1.01) {
    double f = remez.evaluateFunc(x);
    double r = remez.evaluateApprox(x);
    fprintf(error_file,"%e %e\n", x,  (r - f)/f);
  }
  fclose(error_file);

  delete res;
  delete pole;

  exit(0);

}
