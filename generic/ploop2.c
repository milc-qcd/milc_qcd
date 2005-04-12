/****************** ploop2.c ************************************/
/* MIMD version 6 */
/* evaluate the Polyakov loops.  This version uses gathers. */
/* This version uses site structure members tempmat1, tempmat2, staple */

#include "generic_includes.h"

complex ploop() {
  register int i,t;
  register site *st;
  msg_tag *tag;
  complex sum;
  complex plp;
  su3_matrix *tempmat;

  tempmat = (su3_matrix *)malloc(sites_on_node*sizeof(su3_matrix));
  if(tempmat == NULL){
    printf("ploop: Can't allocate temporary\n");
    terminate(1);
  }

  sum = cmplx(0.0,0.0);
  FORALLSITES(i,st){tempmat[i] = lattice[i].link[TUP];}
  for(t=1;t<nt;t++){
    tag=start_gather_field( tempmat, sizeof(su3_matrix),
			    TUP, EVENANDODD, gen_pt[0] );
    wait_gather(tag);
    FORALLSITES(i,st){
      if( st->t != 0 )continue;
      if(t==1){
	mult_su3_nn( &(st->link[TUP]), (su3_matrix *)gen_pt[0][i],
		     &(st->staple));
      }
      else {
	mult_su3_nn( &(st->staple), (su3_matrix *)gen_pt[0][i],
		     &(tempmat[i]));
	lattice[i].staple = tempmat[i];
      }
    }
    FORALLSITES(i,st){
      st->tempmat1 = *(su3_matrix *)(gen_pt[0][i]);
    }
    FORALLSITES(i,st){
     tempmat[i] = st->tempmat1;
    }
    cleanup_gather(tag);
  }
  FORALLSITES(i,st){
    if( st->t != 0 )continue;
    plp = trace_su3( &(st->staple) );
    CSUM(sum,plp);
  }
  g_complexsum( &sum );
  plp.real = sum.real /((Real)(nx*ny*nz));
  plp.imag = sum.imag /((Real)(nx*ny*nz));
  return(plp);
}
