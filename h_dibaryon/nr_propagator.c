/********** nr_propagator.c *************/
/* MIMD version 6 */
/* 4/5/97 C. DeTar */

/* Project Dirac propagator onto nonrelativistic spin state 
   to get Pauli propagator */

#include "h_dibaryon_includes.h"

/* Create non-relativistic quark propagator from Wilson propagator */
/* These are defined as projections onto +1 eigenstates of gamma(TUP)
   as specified in wp_shrink:

 gamma(TUP)			eigenvectors	eigenvalue
 	    0  0  1  0		( 1, 0, 1, 0)	+1
            0  0  0  1		( 0, 1, 0, 1)	+1
            1  0  0  0		( 1, 0,-1, 0)	-1
            0  1  0  0		( 0, 1, 0,-1)	-1

   For forward propagation (forw_back == 1),
       Pauli spin 0 is Dirac (1,0,-1,0)
   and Pauli spin 1 is Dirac (0,1,0,-1)

   For backward propagation (forw_back == 2),
       Pauli spin 0 is Dirac (1,0,1,0)
   and Pauli spin 1 is Dirac (0,1,0,1)


	    */

void nr_propagator(field_offset qk, field_offset nr, int forw_back)
{
  register int i;
  register site *s;
  int ci,cj;

  if(forw_back == 1)
    {
      
      FORALLSITES(i,s)
	{
	  for(ci = 0; ci < 3; ci++) for(cj = 0; cj < 3; cj++)
	    {
	      CSUB(
		   ((wilson_propagator *)F_PT(s,qk))->c[ci].d[0].d[0].c[cj],
		   ((wilson_propagator *)F_PT(s,qk))->c[ci].d[2].d[0].c[cj],
		   ((pauli_propagator  *)F_PT(s,nr))->c[ci].p[0].p[0].c[cj]);
	      CSUB(
		   ((pauli_propagator  *)F_PT(s,nr))->c[ci].p[0].p[0].c[cj],
		   ((wilson_propagator *)F_PT(s,qk))->c[ci].d[0].d[2].c[cj],
		   ((pauli_propagator  *)F_PT(s,nr))->c[ci].p[0].p[0].c[cj]);
	      CSUM(
		   ((pauli_propagator  *)F_PT(s,nr))->c[ci].p[0].p[0].c[cj],
		   ((wilson_propagator *)F_PT(s,qk))->c[ci].d[2].d[2].c[cj]);
	      
	      CSUB(
		   ((wilson_propagator *)F_PT(s,qk))->c[ci].d[1].d[0].c[cj],
		   ((wilson_propagator *)F_PT(s,qk))->c[ci].d[3].d[0].c[cj],
		   ((pauli_propagator  *)F_PT(s,nr))->c[ci].p[1].p[0].c[cj]);
	      CSUB(
		   ((pauli_propagator  *)F_PT(s,nr))->c[ci].p[1].p[0].c[cj],
		   ((wilson_propagator *)F_PT(s,qk))->c[ci].d[1].d[2].c[cj],
		   ((pauli_propagator  *)F_PT(s,nr))->c[ci].p[1].p[0].c[cj]);
	      CSUM(
		   ((pauli_propagator  *)F_PT(s,nr))->c[ci].p[1].p[0].c[cj],
		   ((wilson_propagator *)F_PT(s,qk))->c[ci].d[3].d[2].c[cj]);
	      
	      CSUB(
		   ((wilson_propagator *)F_PT(s,qk))->c[ci].d[0].d[1].c[cj],
		   ((wilson_propagator *)F_PT(s,qk))->c[ci].d[2].d[1].c[cj],
		   ((pauli_propagator  *)F_PT(s,nr))->c[ci].p[0].p[1].c[cj]);
	      CSUB(
		   ((pauli_propagator  *)F_PT(s,nr))->c[ci].p[0].p[1].c[cj],
		   ((wilson_propagator *)F_PT(s,qk))->c[ci].d[0].d[3].c[cj],
		   ((pauli_propagator  *)F_PT(s,nr))->c[ci].p[0].p[1].c[cj]);
	      CSUM(
		   ((pauli_propagator  *)F_PT(s,nr))->c[ci].p[0].p[1].c[cj],
		   ((wilson_propagator *)F_PT(s,qk))->c[ci].d[2].d[3].c[cj]);
	      
	      CSUB(
		   ((wilson_propagator *)F_PT(s,qk))->c[ci].d[1].d[1].c[cj],
		   ((wilson_propagator *)F_PT(s,qk))->c[ci].d[3].d[1].c[cj],
		   ((pauli_propagator  *)F_PT(s,nr))->c[ci].p[1].p[1].c[cj]);
	      CSUB(
		   ((pauli_propagator  *)F_PT(s,nr))->c[ci].p[1].p[1].c[cj],
		   ((wilson_propagator *)F_PT(s,qk))->c[ci].d[1].d[3].c[cj],
		   ((pauli_propagator  *)F_PT(s,nr))->c[ci].p[1].p[1].c[cj]);
	      CSUM(
		   ((pauli_propagator  *)F_PT(s,nr))->c[ci].p[1].p[1].c[cj],
		   ((wilson_propagator *)F_PT(s,qk))->c[ci].d[3].d[3].c[cj]);
	    }
	}
    }
  else
    {
      
      FORALLSITES(i,s)
	{
	  for(ci = 0; ci < 3; ci++) for(cj = 0; cj < 3; cj++)
	    {
	      CADD(
		   ((wilson_propagator *)F_PT(s,qk))->c[ci].d[0].d[0].c[cj],
		   ((wilson_propagator *)F_PT(s,qk))->c[ci].d[2].d[0].c[cj],
		   ((pauli_propagator  *)F_PT(s,nr))->c[ci].p[0].p[0].c[cj]);
	      CADD(
		   ((pauli_propagator  *)F_PT(s,nr))->c[ci].p[0].p[0].c[cj],
		   ((wilson_propagator *)F_PT(s,qk))->c[ci].d[0].d[2].c[cj],
		   ((pauli_propagator  *)F_PT(s,nr))->c[ci].p[0].p[0].c[cj]);
	      CSUM(
		   ((pauli_propagator  *)F_PT(s,nr))->c[ci].p[0].p[0].c[cj],
		   ((wilson_propagator *)F_PT(s,qk))->c[ci].d[2].d[2].c[cj]);
	      
	      CADD(
		   ((wilson_propagator *)F_PT(s,qk))->c[ci].d[1].d[0].c[cj],
		   ((wilson_propagator *)F_PT(s,qk))->c[ci].d[3].d[0].c[cj],
		   ((pauli_propagator  *)F_PT(s,nr))->c[ci].p[1].p[0].c[cj]);
	      CADD(
		   ((pauli_propagator  *)F_PT(s,nr))->c[ci].p[1].p[0].c[cj],
		   ((wilson_propagator *)F_PT(s,qk))->c[ci].d[1].d[2].c[cj],
		   ((pauli_propagator  *)F_PT(s,nr))->c[ci].p[1].p[0].c[cj]);
	      CSUM(
		   ((pauli_propagator  *)F_PT(s,nr))->c[ci].p[1].p[0].c[cj],
		   ((wilson_propagator *)F_PT(s,qk))->c[ci].d[3].d[2].c[cj]);
	      
	      CADD(
		   ((wilson_propagator *)F_PT(s,qk))->c[ci].d[0].d[1].c[cj],
		   ((wilson_propagator *)F_PT(s,qk))->c[ci].d[2].d[1].c[cj],
		   ((pauli_propagator  *)F_PT(s,nr))->c[ci].p[0].p[1].c[cj]);
	      CADD(
		   ((pauli_propagator  *)F_PT(s,nr))->c[ci].p[0].p[1].c[cj],
		   ((wilson_propagator *)F_PT(s,qk))->c[ci].d[0].d[3].c[cj],
		   ((pauli_propagator  *)F_PT(s,nr))->c[ci].p[0].p[1].c[cj]);
	      CSUM(
		   ((pauli_propagator  *)F_PT(s,nr))->c[ci].p[0].p[1].c[cj],
		   ((wilson_propagator *)F_PT(s,qk))->c[ci].d[2].d[3].c[cj]);
	      
	      CADD(
		   ((wilson_propagator *)F_PT(s,qk))->c[ci].d[1].d[1].c[cj],
		   ((wilson_propagator *)F_PT(s,qk))->c[ci].d[3].d[1].c[cj],
		   ((pauli_propagator  *)F_PT(s,nr))->c[ci].p[1].p[1].c[cj]);
	      CADD(
		   ((pauli_propagator  *)F_PT(s,nr))->c[ci].p[1].p[1].c[cj],
		   ((wilson_propagator *)F_PT(s,qk))->c[ci].d[1].d[3].c[cj],
		   ((pauli_propagator  *)F_PT(s,nr))->c[ci].p[1].p[1].c[cj]);
	      CSUM(
		   ((pauli_propagator  *)F_PT(s,nr))->c[ci].p[1].p[1].c[cj],
		   ((wilson_propagator *)F_PT(s,qk))->c[ci].d[3].d[3].c[cj]);
	    }
	}
    }
} /* nr_propagator */

