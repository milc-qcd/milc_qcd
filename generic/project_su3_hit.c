/***************** project_su3_hit.c *******************************/
/* Project a single 3x3 complex matrix Q onto SU(3) W by maximizing
   Re Tr (W Q+).  This code uses an iterative method based on hits
   in the various diagonal SU(2) subgroups.  A convergence test
   is also available. 
*/

/* MIMD version 7 */

/* C. DeTar 11/30/97 based on monte.c by T. DeGrand March 1991.
   See cubic.c from Tom's inverse blocking code for a more
   direct method. 
   CD 11/25/01 Included features from Urs Heller's smearing.c */

#include "generic_includes.h"

#define Nc 3

void project_su3(
   su3_matrix *w,         /* input initial guess. output resulting
                             SU(3) matrix */
   su3_matrix *q,         /* 3 x 3 complex matrix to be projected */
   int Nhit,              /* number of SU(2) hits. 0 for no projection */
   Real tol              /* tolerance for SU(3) projection.
			     If nonzero, treat Nhit as a maximum
			     number of hits.  If zero, treat Nhit
			     as a prescribed number of hits. */ 
)
{
   int index1, ina, inb,ii;
   Real v0,v1,v2,v3, vsq;
   Real z;
   su3_matrix action;
   su2_matrix h;
   Real conver, old_tr = 0, new_tr;

   if(tol > 0)
     old_tr = realtrace_su3(w,q)/3.0;
   conver = 1.0;

   /* Do SU(2) hits */
   for(index1=0;index1<Nhit && conver > tol; index1++)
   {
      /*  pick out an SU(2) subgroup */
      ina =  index1    % Nc;
      inb = (index1+1) % Nc;
      if(ina > inb){ ii=ina; ina=inb; inb=ii; }

      mult_su3_na( w, q, &action );

      /* decompose the action into SU(2) subgroups using
         Pauli matrix expansion */
      /* The SU(2) hit matrix is represented as v0 + i *
         Sum j (sigma j * vj)*/
      v0 =  action.e[ina][ina].real + action.e[inb][inb].real;
      v3 =  action.e[ina][ina].imag - action.e[inb][inb].imag;
      v1 =  action.e[ina][inb].imag + action.e[inb][ina].imag;
      v2 =  action.e[ina][inb].real - action.e[inb][ina].real;

      /* Normalize v */
      vsq = v0*v0 + v1*v1 + v2*v2 + v3*v3;
      z = sqrt((double)vsq );
      if(z == 0.){z = 1.;v0 = 1.;}
      else {v0 = v0/z; v1 = v1/z; v2 = v2/z; v3 = v3/z;}

      /* Elements of SU(2) matrix */

      h.e[0][0] = cmplx( v0,-v3);
      h.e[0][1] = cmplx(-v2,-v1);
      h.e[1][0] = cmplx( v2,-v1);
      h.e[1][1] = cmplx( v0, v3);

      /* update the link */
      left_su2_hit_n(&h,ina,inb,w);

      /* convergence measure every third hit */
      if(tol>0 && (index1 % 3) == 2){
	new_tr = realtrace_su3(w,q)/3.;
	conver = (new_tr-old_tr)/old_tr; /* trace always increases */
	old_tr = new_tr;
      }
      
   } /* hits */
   
   if( Nhit > 0 && tol > 0 && conver > tol )
     printf("project_su3: node %d No convergence: conver = %e\n",
	    this_node, (double)conver);
}
