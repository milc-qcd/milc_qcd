/*

  Mike Clark - 25th May 2005

  alg_remez.h

  AlgRemez is an implementation of the Remez algorithm, which in this
  case is used for generating the optimal nth root rational
  approximation.

  Note this class requires the gnu multiprecision (GNU MP) library.

*/

#ifndef INCLUDED_ALG_REMEZ_H
#define INCLUDED_ALG_REMEZ_H

#include "bigfloat.h"

#define JMAX 10000 //Maximum number of iterations of Newton's approximation
#define SUM_MAX 10 // Maximum number of terms in exponential

class AlgRemez
{
 private:
  char *cname;

  // The approximation parameters
  bigfloat *param, *roots, *poles;
  bigfloat norm;

  // The numerator and denominator degree (n=d)
  int n, d;
  
  // The bounds of the approximation
  bigfloat apstrt, apwidt, apend;

  // The number of flavors.  This can be 1 to 4 only.
  int num_flav;

  // The numerator and denominator of the power we are approximating for flavor 1.
  unsigned long power_num_flav_1;
  unsigned long power_den_flav_1;
  // The numerator and denominator of the power we are approximating for flavor 2.
  unsigned long power_num_flav_2;
  unsigned long power_den_flav_2;
  // The numerator and denominator of the power we are approximating for flavor 3.
  unsigned long power_num_flav_3;
  unsigned long power_den_flav_3;
  // The numerator and denominator of the power we are approximating for flavor 4.
  unsigned long power_num_flav_4;
  unsigned long power_den_flav_4;

  // The mass for flavor 1.
  bigfloat mass_flav_1;
  // The mass for flavor 2.
  bigfloat mass_flav_2;
  // The mass for flavor 3.
  bigfloat mass_flav_3;
  // The mass for flavor 4.
  bigfloat mass_flav_4;

  // Flag to determine whether the arrays have been allocated
  int alloc;

  // Flag to determine whether the roots have been found
  int foundRoots;

  // Variables used to calculate the approximation
  int nd1, iter;
  bigfloat *xx, *mm, *step;
  bigfloat delta, spread, tolerance;

  // The exponential summation coefficients
  bigfloat *a;
  int *a_power;
  int a_length;

  // The number of equations we must solve at each iteration (n+d+1)
  int neq;

  // The precision of the GNU MP library
  long prec;

  // Initial values of maximal and minmal errors
  void initialGuess();

  // Solve the equations
  void equations();

  // Search for error maxima and minima
  void search(bigfloat *step); 

  // Initialise step sizes
  void stpini(bigfloat *step);

  // Calculate the roots of the approximation
  int root();

  // Evaluate the polynomial
  bigfloat polyEval(bigfloat x, bigfloat *poly, long size);
  //complex_bf polyEval(complex_bf x, complex_bf *poly, long size);

  // Evaluate the differential of the polynomial
  bigfloat polyDiff(bigfloat x, bigfloat *poly, long size);
  //complex_bf polyDiff(complex_bf x, complex_bf *poly, long size);

  // Newton's method to calculate roots
  bigfloat rtnewt(bigfloat *poly, long i, bigfloat x1, bigfloat x2, bigfloat xacc);
  //complex_bf rtnewt(complex_bf *poly, long i, bigfloat x1, bigfloat x2, bigfloat xacc);

  // Evaluate the partial fraction expansion of the rational function
  // with res roots and poles poles.  Result is overwritten on input
  // arrays.
  void pfe(bigfloat *res, bigfloat* poles, bigfloat norm);

  // Calculate function required for the approximation
  bigfloat func(bigfloat x);

  // Compute size and sign of the approximation error at x
  bigfloat getErr(bigfloat x, int *sign);

  // Solve the system AX=B
  int simq(bigfloat *A, bigfloat *B, bigfloat *X, int n);

  // Free memory and reallocate as necessary
  void allocate(int num_degree, int den_degree);

  // Evaluate the rational form P(x)/Q(x) using coefficients from the
  // solution vector param
  bigfloat approx(bigfloat x);

 public:
  
  // Constructor
  AlgRemez(double lower, double upper, long prec);

  // Destructor
  virtual ~AlgRemez();

  // Reset the bounds of the approximation
  void setBounds(double lower, double upper);

  // Generate the rational approximation x^(pnum/pden)
  double generateApprox(int num_degree, int den_degree, 
			unsigned long power_num1, unsigned long power_den1, double m1,
			unsigned long power_num2, unsigned long power_den2, double m2,
			unsigned long power_num3, unsigned long power_den3, double m3,
			unsigned long power_num4, unsigned long power_den4, double m4,
			int a_len, double* a_param, int* a_pow);
  double generateApprox(int num_degree, int den_degree, 
			unsigned long power_num1, unsigned long power_den1, double m1,
			unsigned long power_num2, unsigned long power_den2, double m2,
			unsigned long power_num3, unsigned long power_den3, double m3,
			unsigned long power_num4, unsigned long power_den4, double m4);
  double generateApprox(int degree,
                        unsigned long power_num1, unsigned long power_den1, double m1,
                        unsigned long power_num2, unsigned long power_den2, double m2,
                        unsigned long power_num3, unsigned long power_den3, double m3,
                        unsigned long power_num4, unsigned long power_den4, double m4);

  // Return the partial fraction expansion of the approximation x^(pnum/pden)
  int getPFE(double *res, double *pole, double *norm);

  // Return the partial fraction expansion of the approximation x^(-pnum/pden)
  int getIPFE(double *res, double *pole, double *norm);

  // Evaluate the rational form P(x)/Q(x) using coefficients from the
  // solution vector param
  double evaluateApprox(double x);

  // Evaluate the rational form Q(x)/P(x) using coefficients from the
  // solution vector param
  double evaluateInverseApprox(double x);

  // Calculate function required for the approximation
  double evaluateFunc(double x);

  // Calculate inverse function required for the approximation
  double evaluateInverseFunc(double x);

};

#endif  // Include guard



