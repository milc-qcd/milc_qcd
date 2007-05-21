/************************ nersc_cksum.c *******************************/
/* MIMD version 7 */

/* Two utilities used in the NERSC archive formate */

/* Compute the low order 32 bits of the unsigned integer sum of the
   float precision real and complex parts of the elements of the gauge
   matrices.
*/

/* Computes the mean global sum of the trace of the gauge links --
   used to aid checking lattice file integrity */

#include "generic_includes.h"

/**
static int 
my_big_endian() {
  union  {
    long l;
    char c[sizeof (long)];
  } u;
  u.l = 1;
  return (u.c[sizeof (long) - 1] == 1);
}
**/

u_int32type 
nersc_cksum( void ) {
  int i,mu,a,b;
  site *s;
  u_int32type chksum = 0;
  union {
    float       flt;
    u_int32type p32;
  } tmp;
  u_int32type p32;

  FORALLSITES(i,s) {
    for(mu=0; mu<4; ++mu) {
      for(a=0; a<2; a++) for(b=0; b<3; b++) {
	tmp.flt = s->link[mu].e[a][b].real;
	p32 = tmp.p32;
	chksum += p32;
	tmp.flt = s->link[mu].e[a][b].imag;
	p32 = tmp.p32;
	chksum += p32;
      }
    }
  }

  g_uint32sum(&chksum);
  return chksum;

} /* nersc_cksum.c */


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

