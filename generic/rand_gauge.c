/*********************** rand_gauge.c ***************************/
/* MIMD version 7 */

/* original code by UMH */
/* 2/19/98 Version 5 port CD */

/* Makes a random gauge transformation on the gauge fields
   and then reunitarizes them */
/* Warning KS fermion applications: Must be done with KS phases OUT! */
/* Requires site workspace G */

#include "generic_includes.h"

void randomize(field_offset G, Real radius);
void gauge_trans(field_offset G);

void rand_gauge(field_offset G)
	/* G holds the gauge transformation matrices */
{
    randomize(G, 1.0);
    gauge_trans(G);
    reunitarize();
}

void randomize(field_offset G, Real radius)
{
  register int a, b, i;
  site *s;

  FORALLSITES(i,s) {
    for(a=0; a<3; a++) for(b=0; b<3; b++)
      (*(su3_matrix *)F_PT(s,G)).e[a][b] 
             = cmplx(radius*((Real)drand48()-0.5),
                     radius*((Real)drand48()-0.5));
    reunit_su3((su3_matrix *)F_PT(s,G));
  }
}

void gauge_trans(field_offset G)
{
  register int i,mu;
  site *s;
  su3_matrix tmp;
  msg_tag *tag[4];

  FORALLUPDIR(mu) 
    tag[mu] = start_gather_site(G,sizeof(su3_matrix),mu,EVENANDODD,
		       gen_pt[mu]);

  FORALLUPDIR(mu) {
    wait_gather(tag[mu]);
    FORALLSITES(i,s) {

       mult_su3_an((su3_matrix *)F_PT(s,G), &(s->link[mu]), &tmp);
       mult_su3_nn(&tmp, (su3_matrix *)gen_pt[mu][i],
		       &(s->link[mu]));

    }
    cleanup_gather(tag[mu]);
  }
}
