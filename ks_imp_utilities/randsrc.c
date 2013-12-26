/**************** randsrc.c ************************************/
/* MIMC version 7 */
/* Make random source for testing inverter and fermion force term */
/* Adapted from normdist.c */

/*   Generate a list of random numbers according to
     a Gaussian normal distribution  
 
      C. DeTar (4/17/91)   */

/*   This program should be compiled with the command 

     cc normdist.c -o normdist -lm
                                      */

/* stdio is needed for fprintf, fprint, scanf, stderr */
#include <stdio.h>

/* stdlib is needed for srand and rand */
#include <stdlib.h>

/* math is needed for sqrt, log, cos */
#include <math.h>

#define PI 3.14159265358979323
#define TINY 1e-307

/* According to the writeup, rand gives integers up to this value */
/* We want to normalize to the real interval [0,1] via rand()/MAXRAND */
#define MAXRAND 32768

int main()
{
  float mean, stdev, x;
  float v,p,r;
  int n, i;
  unsigned int seed;

  mean = 0;
  stdev = 10;
  fprintf(stderr,"Number of random values?\t");
  scanf("%d",&n);
  fprintf(stderr,"Random number seed\t");
  scanf("%d",&seed);

/* Initialize random number seed */
  srand(seed);

/* Generate list */
  for(i = 0; i < n; i++)
    {
      v = (float)rand()/MAXRAND;
      p = (float)rand()*2.*PI/MAXRAND;
      r = sqrt(2.*(-log(v+TINY)));
      x = r*cos(p);
      x = x*stdev + mean;
      fwrite(&x,sizeof(x),1,stdout);
    }
  return 0;
}
