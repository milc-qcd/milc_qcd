/******* t3e_vhelp1.c -vector routine for T3E */
/* MIMD version 6 */
/* NOT MAINTAINED.  TEST BEFORE USE! */

/* Kogut-Susskind fermions */

#include "schroed_ks_includes.h"
#include "../include/prefetch.h"

void V_mult_adj_su3_mat_vec_4dir( su3_matrix * lpt,
    su3_vector * srcpt, su3_vector * destpt, int n, int stride ){
    register int i;
    n--;	/* last iteration has no prefetch */
    for(i=0;i<n;i++){
      prefetch_V( (su3_vector *)( (char *)srcpt + stride ) );
      mult_adj_su3_mat_vec_4dir( lpt, srcpt, destpt );
      lpt = (su3_matrix *)( (char *)lpt + stride );
      srcpt = (su3_vector *)( (char *)srcpt + stride );
      destpt = (su3_vector *)( (char *)destpt + stride );
    }
    mult_adj_su3_mat_vec_4dir( lpt, srcpt, destpt );
}
