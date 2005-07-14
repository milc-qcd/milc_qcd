/********** boundary_flip.c *************/
/* MIMD version 7 */

/* Flip the time direction boundary conditions.  (Just multiply time
   links on last time slice by -1.)
   argument "sign" is PLUS or MINUS, says what you want the boundary
   conditions to become.  If the boundary conditions are already set
   to "sign", you get a warning. */

#include "generic_wilson_includes.h"

int current_boundary = PLUS;
void boundary_flip( int sign ) {
  register int i,j,k;
  register site *s;
  if(this_node==0 && sign==current_boundary)
    printf("WARNING: you lost track of the boundary conditions!\n");
  if(sign==current_boundary)return;
  
  FORALLSITES(i,s){
    if(s->t != nt-1)continue;	/* hit only last time slice */
    for(j=0;j<3;j++)for(k=0;k<3;k++){
      s->link[TUP].e[j][k].real *= -1.0;
      s->link[TUP].e[j][k].imag *= -1.0;
    }
  }
  
  current_boundary = sign;
}


