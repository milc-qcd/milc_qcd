/******* t3e_vhelp2_SF.c -vector routine for T3E */
/* MIMD version 6 */
/* NOT MAINTAINED.  TEST BEFORE USE! */

/* Kogut-Susskind fermions with Schroedinger Functional bc. */

#include "schroed_ks_includes.h"
#include "../include/prefetch.h"

void V_mult_su3_mat_vec_sum_4dir( su3_matrix * lpt,
    su3_vector ** xpt, su3_vector ** ypt, su3_vector ** zpt, su3_vector ** tpt,
    su3_vector * destpt, int nsites, int stride ){
    register int i;
    register site *s;
    su3_matrix *lpt2;
    nsites--;	/*last site has no prefetch*/
    for(i=0, s=&(lattice[i]); i<nsites; i++, s++){
      if(i != nsites-1){
	prefetch_V( *(xpt+1) );
	prefetch_V( *(ypt+1) );
	prefetch_V( *(zpt+1) );
	prefetch_V( *(tpt+1) );
	prefetch_V( (su3_vector *)( (char *)destpt + stride ) );
      }
	if(s->t > 0){
	    if(s->t == (nt-1)){
		mult_su3_mat_vec( lpt, *xpt, destpt);
		lpt2 = (su3_matrix *)( (char *)lpt + sizeof(su3_matrix) );
		mult_su3_mat_vec_sum( lpt2, *ypt, destpt);
		lpt2 = (su3_matrix *)( (char *)lpt2 + sizeof(su3_matrix) );
		mult_su3_mat_vec_sum( lpt2, *zpt, destpt);
	    }
	    else
		mult_su3_mat_vec_sum_4dir( lpt, *xpt, *ypt, *zpt, *tpt, destpt);
	}
	lpt = (su3_matrix *)( (char *)lpt + stride );
	xpt++; ypt++; zpt++; tpt++;
	destpt = (su3_vector *)( (char *)destpt + stride );
    }
    if(s->t > 0){
	if(s->t == (nt-1)){
	    mult_su3_mat_vec( lpt, *xpt, destpt);
	    lpt2 = (su3_matrix *)( (char *)lpt + sizeof(su3_matrix) );
	    mult_su3_mat_vec_sum( lpt2, *ypt, destpt);
	    lpt2 = (su3_matrix *)( (char *)lpt2 + sizeof(su3_matrix) );
	    mult_su3_mat_vec_sum( lpt2, *zpt, destpt);
	}
	else
	    mult_su3_mat_vec_sum_4dir( lpt, *xpt, *ypt, *zpt, *tpt, destpt);
    }
}
