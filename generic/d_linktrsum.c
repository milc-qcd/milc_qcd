/************************** d_linktrsum.c *******************************/
/* MIMD version 7 */
/* Computes the mean global sum of the trace of the gauge links --
   used to aid checking lattice file integrity */

#include "generic_includes.h"

void d_linktrsum(double_complex *linktrsum) {
  int i,dir;
  site *s;
  su3_matrix *a;

  linktrsum->real = 0.;
  linktrsum->imag = 0.;

  FORALLSITES(i,s){
    FORALLUPDIR(dir){
      a = &s->link[dir];
      CSUM(*linktrsum,a->e[0][0]);
      CSUM(*linktrsum,a->e[1][1]);
      CSUM(*linktrsum,a->e[2][2]);
    }
  }

  g_dcomplexsum(linktrsum);
  CDIVREAL(*linktrsum,(4*volume),*linktrsum);

} /* d_linktrsum */

