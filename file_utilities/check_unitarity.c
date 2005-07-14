/*********************** check_unitarity.c ***************************/
/* MIMD version 7 */
/* Claude Bernard, original version */
/* Modified 7/96 DT to quit after finding violation */
/* 9/4/96 DT added Urs' row orthogonality checks, if STRONG defined */
#define STRONG	/* check row orthogonality as well as norms */

/* Check unitarity of link matrices, quit if not unitary */
#include <stdio.h>
#include <math.h>
#include "../include/complex.h"
#include "../include/su3.h"
#include <lattice.h>
#include "../include/macros.h"
#include "../include/comdefs.h"
#define TOLERANCE (0.0001)
/*#define UNIDEBUG */
Real check_su3();

void check_unitarity() {
register int i,dir;
int ii,jj;
register site *s;
register su3_matrix *mat;
Real deviation,max_deviation;
double av_deviation;
union {
                Real fval;
                int ival;
                } ifval;

    max_deviation=av_deviation=0;
    FORALLSITES(i,s){
        for(dir=XUP; dir<=TUP; dir++ ){
            mat = (su3_matrix *)&(s->link[dir]);
            deviation=check_su3( mat );
            if (deviation>TOLERANCE){
                printf("Unitarity problem on node %d, site %d, dir %d, deviation=%f\n",
                mynode(),i,dir,deviation);
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
           	 if(max_deviation<deviation) max_deviation=deviation;
	         av_deviation += deviation*deviation;
fflush(stdout); terminate(1);
	    }
          }
     }
        av_deviation = sqrt(av_deviation/(4*i));
#ifdef UNIDEBUG
        printf("Deviation from unitarity on node %d: max %g, avrg %g\n",
               mynode(), max_deviation, av_deviation);
#endif
        if(max_deviation> TOLERANCE) 
          printf("Unitarity problem on node %d, maximum deviation=%f\n",
                mynode(),max_deviation);
}  /*check_unitarity() */

Real check_su3(c) su3_matrix *c; {
     register Real ar, ai, ari, max;
     register int i;

     /* first normalize row */
  for(i=0 , max=0.; i<3; ++i) {
     ar = (*c).e[i][0].real * (*c).e[i][0].real +    /* sum of squares of row */
          (*c).e[i][0].imag * (*c).e[i][0].imag +
          (*c).e[i][1].real * (*c).e[i][1].real +
          (*c).e[i][1].imag * (*c).e[i][1].imag +
          (*c).e[i][2].real * (*c).e[i][2].real +
          (*c).e[i][2].imag * (*c).e[i][2].imag;
     ar =  fabs( sqrt((double)ar) - 1.);
     if(max<ar) max=ar;
  }

#ifdef STRONG

     /* Test orthogonality of row 0 and row 1 */
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

  ari = sqrt((double)(ar*ar + ai*ai));
  if(max<ari) max=ari;

     /* Test orthogonality of row 0 and row 2 */
  ar = (*c).e[0][0].real * (*c).e[2][0].real +     /* real part of 0 dot 1 */
       (*c).e[0][0].imag * (*c).e[2][0].imag +
       (*c).e[0][1].real * (*c).e[2][1].real +
       (*c).e[0][1].imag * (*c).e[2][1].imag +
       (*c).e[0][2].real * (*c).e[2][2].real +
       (*c).e[0][2].imag * (*c).e[2][2].imag;
  ai = (*c).e[0][0].real * (*c).e[2][0].imag -     /* imag part of 0 dot 1 */
       (*c).e[0][0].imag * (*c).e[2][0].real +
       (*c).e[0][1].real * (*c).e[2][1].imag -
       (*c).e[0][1].imag * (*c).e[2][1].real +
       (*c).e[0][2].real * (*c).e[2][2].imag -
       (*c).e[0][2].imag * (*c).e[2][2].real;

  ari = sqrt((double)(ar*ar + ai*ai));
  if(max<ari) max=ari;

     /* Test orthogonality of row 1 and row 2 */
  ar = (*c).e[1][0].real * (*c).e[2][0].real +     /* real part of 0 dot 1 */
       (*c).e[1][0].imag * (*c).e[2][0].imag +
       (*c).e[1][1].real * (*c).e[2][1].real +
       (*c).e[1][1].imag * (*c).e[2][1].imag +
       (*c).e[1][2].real * (*c).e[2][2].real +
       (*c).e[1][2].imag * (*c).e[2][2].imag;
  ai = (*c).e[1][0].real * (*c).e[2][0].imag -     /* imag part of 0 dot 1 */
       (*c).e[1][0].imag * (*c).e[2][0].real +
       (*c).e[1][1].real * (*c).e[2][1].imag -
       (*c).e[1][1].imag * (*c).e[2][1].real +
       (*c).e[1][2].real * (*c).e[2][2].imag -
       (*c).e[1][2].imag * (*c).e[2][2].real;

  ari = sqrt((double)(ar*ar + ai*ai));
  if(max<ari) max=ari;

#endif /*STRONG*/

     return(max);

} /* check_su3 */
