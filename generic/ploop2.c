/****************** ploop2.c ************************************/
/* MIMD version 7 */
/* evaluate the Polyakov loops.  This version uses gathers. */

#include "generic_includes.h"

complex ploop() {
  register int i,t;
  register site *st;
  msg_tag *tag;
  complex sum;
  complex plp;
  su3_matrix *tempmat1, *tempmat2, *staple;

  staple = (su3_matrix *)malloc(sites_on_node*sizeof(su3_matrix));
  if(staple == NULL){
    printf("ploop: Can't allocate temporary\n");
    terminate(1);
  }

  tempmat1 = (su3_matrix *)malloc(sites_on_node*sizeof(su3_matrix));
  if(tempmat1 == NULL){
    printf("ploop: Can't allocate temporary\n");
    terminate(1);
  }

  tempmat2 = (su3_matrix *)malloc(sites_on_node*sizeof(su3_matrix));
  if(tempmat2 == NULL){
    printf("ploop: Can't allocate temporary\n");
    terminate(1);
  }

  sum = cmplx(0.0,0.0);
  FORALLSITES(i,st){tempmat2[i] = lattice[i].link[TUP];}
  for(t=1;t<nt;t++){
    tag=start_gather_field( tempmat2, sizeof(su3_matrix),
			    TUP, EVENANDODD, gen_pt[0] );
    wait_gather(tag);
    FORALLSITES(i,st){
      if( st->t != 0 )continue;
      if(t==1){
	mult_su3_nn( &(st->link[TUP]), (su3_matrix *)gen_pt[0][i],
		     &staple[i]);
      }
      else {
	mult_su3_nn( &staple[i], (su3_matrix *)gen_pt[0][i],
		     &(tempmat2[i]));
	staple[i] = tempmat2[i];
      }
    }
    FORALLSITES(i,st){
      tempmat1[i] = *(su3_matrix *)(gen_pt[0][i]);
    }
    FORALLSITES(i,st){
     tempmat2[i] = tempmat1[i];
    }
    cleanup_gather(tag);
  }
  FORALLSITES(i,st){
    if( st->t != 0 )continue;
    plp = trace_su3( &staple[i] );
    CSUM(sum,plp);
  }
  g_complexsum( &sum );
  plp.real = sum.real /((Real)(nx*ny*nz));
  plp.imag = sum.imag /((Real)(nx*ny*nz));
  free(tempmat1);
  free(tempmat2);
  free(staple);
  return(plp);
}
