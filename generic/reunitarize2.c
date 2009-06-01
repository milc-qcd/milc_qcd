/*********************** reunitarize2.c ***************************/
/* MIMD version 7 */

/* reunitarize the link matrices */

/* Modifications ...

   09/03/96 Incorporates unitarity checking C.D.
   01/20/00 combined with Schroedinger functional version - UMH
   03/28/00 Modified reunit_su3 to allow calling program to handle errors.
*/

#include "generic_includes.h"

#define TOLERANCE (0.0001)
#define MAXERRCOUNT 100
/**#define UNIDEBUG**/

Real max_deviation;
double av_deviation;

/* canopy qcdlib code - stolen, of course */
#define fixsu3(matrix) \
 { \
    bj0r = (*matrix).e[0][0].real; \
    bj0i = (*matrix).e[0][0].imag; \
    bj1r = (*matrix).e[0][1].real; \
    bj1i = (*matrix).e[0][1].imag; \
    bj2r = (*matrix).e[0][2].real; \
    bj2i = (*matrix).e[0][2].imag; \
    ar = (*matrix).e[1][2].real; \
    ai = (*matrix).e[1][2].imag; \
    tr = bj1r*ar - bj1i*ai; \
    ti = bj1r*ai + bj1i*ar; \
    ar = (*matrix).e[1][1].real; \
    ai = (*matrix).e[1][1].imag; \
    tr = tr - bj2r*ar + bj2i*ai; \
    ti = ti - bj2r*ai - bj2i*ar; \
    (*matrix).e[2][0].real = tr; \
    (*matrix).e[2][0].imag = -ti; \
    ar = (*matrix).e[1][0].real; \
    ai = (*matrix).e[1][0].imag; \
    tr = bj2r*ar - bj2i*ai; \
    ti = bj2r*ai + bj2i*ar; \
    ar = (*matrix).e[1][2].real; \
    ai = (*matrix).e[1][2].imag; \
    tr = tr - bj0r*ar + bj0i*ai; \
    ti = ti - bj0r*ai - bj0i*ar; \
    (*matrix).e[2][1].real = tr; \
    (*matrix).e[2][1].imag = -ti; \
    ar = (*matrix).e[1][1].real; \
    ai = (*matrix).e[1][1].imag; \
    tr = bj0r*ar - bj0i*ai; \
    ti = bj0r*ai + bj0i*ar; \
    ar = (*matrix).e[1][0].real; \
    ai = (*matrix).e[1][0].imag; \
    tr = tr - bj1r*ar + bj1i*ai; \
    ti = ti - bj1r*ai - bj1i*ar; \
    (*matrix).e[2][2].real = tr; \
    (*matrix).e[2][2].imag = -ti; \
 } /* define fixsu3 */

static int check_deviation(Real deviation)
{

  if(max_deviation<deviation) max_deviation=deviation;
  av_deviation += deviation*deviation;

  if (deviation<TOLERANCE){
    return 0;
  }
  else return 1;
} /* check_deviation */

void reunit_report_problem_matrix(su3_matrix *mat, int i,int dir)
{
  int ii,jj;
  union {
    Real fval;
    int ival;
  } ifval;

  printf("Unitarity problem on node %d, site %d, dir %d tolerance=%e\n",
	 mynode(),i,dir,TOLERANCE);
  printf("SU3 matrix:\n");
  for(ii=0;ii<=2;ii++){
    for(jj=0;jj<=2;jj++){
      printf("%f ",(*mat).e[ii][jj].real); 
      printf("%f ",(*mat).e[ii][jj].imag); 
    }
    printf("\n");
  }
  printf("repeat in hex:\n");
  for(ii=0;ii<=2;ii++){
    for(jj=0;jj<=2;jj++){
      ifval.fval = (*mat).e[ii][jj].real; 
      printf("%08x ", ifval.ival); 
      ifval.fval = (*mat).e[ii][jj].imag; 
      printf("%08x ", ifval.ival); 
    }
    printf("\n");
  }
  printf("  \n \n");
  fflush(stdout); 
} /* reunit_report_problem_matrix */

int reunit_su3(su3_matrix *c)
{
     register Real bj0r, bj0i, bj1r, bj1i, bj2r, bj2i;
     register Real c0r, c0i, c1r, c1i, c2r, c2i;
     register Real ar, ai, tr, ti;
     Real deviation;
     int errors;

     errors = 0;
     /* first normalize row 0 */
     ar = (*c).e[0][0].real * (*c).e[0][0].real +    /* sum of squares of row */
          (*c).e[0][0].imag * (*c).e[0][0].imag +
          (*c).e[0][1].real * (*c).e[0][1].real +
          (*c).e[0][1].imag * (*c).e[0][1].imag +
          (*c).e[0][2].real * (*c).e[0][2].real +
          (*c).e[0][2].imag * (*c).e[0][2].imag;
  
     deviation = fabs(ar - 1.);
     errors += check_deviation(deviation);
     
     ar = 1.0 / sqrt( (double)ar);	       /* used to normalize row */
     (*c).e[0][0].real *= ar;
     (*c).e[0][0].imag *= ar;
     (*c).e[0][1].real *= ar;
     (*c).e[0][1].imag *= ar;
     (*c).e[0][2].real *= ar;
     (*c).e[0][2].imag *= ar;

     /* now make row 1 orthogonal to row 0 */
     ar = (*c).e[0][0].real * (*c).e[1][0].real +     /* real part of 0 dot 1 */
          (*c).e[0][0].imag * (*c).e[1][0].imag +
          (*c).e[0][1].real * (*c).e[1][1].real +
          (*c).e[0][1].imag * (*c).e[1][1].imag +
          (*c).e[0][2].real * (*c).e[1][2].real +
          (*c).e[0][2].imag * (*c).e[1][2].imag;
     ai = (*c).e[0][0].real * (*c).e[1][0].imag -     /* imag part of 0 dot 1 */
          (*c).e[0][0].imag * (*c).e[1][0].real +
          (*c).e[0][1].real * (*c).e[1][1].imag -
          (*c).e[0][1].imag * (*c).e[1][1].real +
          (*c).e[0][2].real * (*c).e[1][2].imag -
          (*c).e[0][2].imag * (*c).e[1][2].real;

     deviation = ar*ar + ai*ai;
     errors += check_deviation(deviation);

     /* row 1 -= a * row 0 */
     (*c).e[1][0].real -= ar*(*c).e[0][0].real - ai*(*c).e[0][0].imag;
     (*c).e[1][0].imag -= ar*(*c).e[0][0].imag + ai*(*c).e[0][0].real;
     (*c).e[1][1].real -= ar*(*c).e[0][1].real - ai*(*c).e[0][1].imag;
     (*c).e[1][1].imag -= ar*(*c).e[0][1].imag + ai*(*c).e[0][1].real;
     (*c).e[1][2].real -= ar*(*c).e[0][2].real - ai*(*c).e[0][2].imag;
     (*c).e[1][2].imag -= ar*(*c).e[0][2].imag + ai*(*c).e[0][2].real;

     /* now normalize row 1 */
     ar = (*c).e[1][0].real * (*c).e[1][0].real +    /* sum of squares of row */
          (*c).e[1][0].imag * (*c).e[1][0].imag +
          (*c).e[1][1].real * (*c).e[1][1].real +
          (*c).e[1][1].imag * (*c).e[1][1].imag +
          (*c).e[1][2].real * (*c).e[1][2].real +
          (*c).e[1][2].imag * (*c).e[1][2].imag;
  
     deviation = fabs(ar - 1.);
     errors += check_deviation(deviation);

     ar = 1.0 / sqrt( (double)ar);	       /* used to normalize row */
     (*c).e[1][0].real *= ar;
     (*c).e[1][0].imag *= ar;
     (*c).e[1][1].real *= ar;
     (*c).e[1][1].imag *= ar;
     (*c).e[1][2].real *= ar;
     (*c).e[1][2].imag *= ar;

     /* Save for checking */
     c0r = (*c).e[2][0].real;
     c0i = (*c).e[2][0].imag;
     c1r = (*c).e[2][1].real;
     c1i = (*c).e[2][1].imag;
     c2r = (*c).e[2][2].real;
     c2i = (*c).e[2][2].imag;

     fixsu3(c); /* reconstruct row 2 */

     /* Now check deviation */
     ar = (c0r - (*c).e[2][0].real) * (c0r - (*c).e[2][0].real) +
          (c0i - (*c).e[2][0].imag) * (c0i - (*c).e[2][0].imag) +
          (c1r - (*c).e[2][1].real) * (c1r - (*c).e[2][1].real) +
          (c1i - (*c).e[2][1].imag) * (c1i - (*c).e[2][1].imag) +
          (c2r - (*c).e[2][2].real) * (c2r - (*c).e[2][2].real) +
          (c2i - (*c).e[2][2].imag) * (c2i - (*c).e[2][2].imag);
     deviation = ar;
     errors += check_deviation(deviation);

     return errors;

} /* reunit_su3 */

void reunitarize() {
  register su3_matrix *mat;
  register int i,dir;
  register site *s;
  int errcount = 0;
  int errors;

  max_deviation = 0.;
  av_deviation = 0.;
  
  FORALLSITES(i,s){
#ifdef SCHROED_FUN
  for(dir=XUP; dir<=TUP; dir++ ) if(dir==TUP || s->t>0 ){
#else
    for(dir=XUP; dir<=TUP; dir++ ){
#endif
      mat = (su3_matrix *)&(s->link[dir]);
      errors = reunit_su3( mat );
      errcount += errors;
      if(errors)reunit_report_problem_matrix(mat,i,dir);
      if(errcount > MAXERRCOUNT)
	{
	  printf("Unitarity error count exceeded.\n");
	  terminate(1);
	}
    }
  }

#ifdef UNIDEBUG
  printf("Deviation from unitarity on node %d: max %.3e, avrg %.3e\n",
	 mynode(), max_deviation, av_deviation);
#endif
  if(max_deviation> TOLERANCE) 
    {
      printf("reunitarize: Node %d unitarity problem, maximum deviation=%e\n",
	     mynode(),max_deviation);
      errcount++;
      if(errcount > MAXERRCOUNT)
	{
	  printf("Unitarity error count exceeded.\n");
	  terminate(1);
	}
    }

}  /* reunitarize2 */

