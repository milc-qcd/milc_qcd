/************************** d_plaq6.c *******************************/
/* MIMD version 7 */
/* This version mallocs the temporary su3_matrix */
/* This version uses gathers to get the neighbors */
/* UMH: Combined with Schroedinger functional version, Jan 2000 */
/* JHW: All six plaquettes separated to be useful for anisotropic settings, Feb 2020 */

/* Measure the six average plaquettes Pmunu */

#include "generic_includes.h"
#include "../include/openmp_defs.h"
#include <assert.h>

void d_plaquette6(double plaq[]) {
register int i,mu1,mu2,imunu;
register site *s;
register su3_matrix *m1,*m4;
double sum[]={0.,0.,0.,0.,0.,0.};
msg_tag *mtag0,*mtag1;
/* su3mat is scratch space of size su3_matrix */
su3_matrix mtmp, *su3mat=NULL;

  su3mat = (su3_matrix *)malloc(sizeof(su3_matrix)*sites_on_node);
  assert( (su3mat!=NULL) );

  for(imunu=0,mu1=YUP;mu1<=TUP;mu1++){
    for(mu2=XUP;mu2<mu1;mu2++,imunu++){

      mtag0 = start_gather_site( F_OFFSET(link[mu2]), sizeof(su3_matrix),mu1, EVENANDODD, gen_pt[0] );
      mtag1 = start_gather_site( F_OFFSET(link[mu1]), sizeof(su3_matrix),mu2, EVENANDODD, gen_pt[1] );

      FORALLSITES_OMP(i,s, private(m1,m4) ){
        m1 = &(s->link[mu1]);
        m4 = &(s->link[mu2]);
        mult_su3_an(m4,m1,&(su3mat[i]));
      } END_LOOP_OMP

      wait_gather(mtag0);
      wait_gather(mtag1);
      FORALLSITES_OMP(i,s,private(mtmp) reduction(+:sum[imunu]) ){
#ifdef SCHROED_FUN
        if( mu1 == TUP ){
          if( s->t == (nt-1) ){
            mult_su3_nn( &(su3mat[i]),&(s->boundary[mu2]), &(mtmp));
          }else{
            mult_su3_nn( &(su3mat[i]),(su3_matrix *)(gen_pt[0][i]), &(mtmp));
          }
        }else if(s->t > 0){
          mult_su3_nn( &(su3mat[i]), (su3_matrix *)(gen_pt[0][i]),&(mtmp));
        }
        if( mu1 == TUP || s->t > 0) 
#else
        mult_su3_nn( &(su3mat[i]),(su3_matrix *)(gen_pt[0][i]),&(mtmp) );
#endif
        sum[imunu] += (double)realtrace_su3((su3_matrix *)(gen_pt[1][i]),&(mtmp) );
      } END_LOOP_OMP
      cleanup_gather(mtag0);
      cleanup_gather(mtag1);
#ifdef SCHROED_FUN
      sum[imunu] /= (mu1 < TUP ? nt-1. : nt*1. );
#else 
      sum[imunu] /= (nt*1.);
#endif
    }
  }
  g_vecdoublesum(sum,6);

  for ( imunu=0; imunu<6; imunu++ ) plaq[imunu] = sum[imunu]/((double)(nx*ny*nz));

  free(su3mat);
} /* d_plaquette6 */

