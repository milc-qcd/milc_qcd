/********** diquarkprop.c *************/
/* MIMD version 6 */
/* 4/5/97 C. DeTar */
/* 04/06/00 Corrected CSUB bug */

/* Construct diquark propagator from nonrelativistic propagator */

#include "h_dibaryon_includes.h"

/* Create diquark hashing tables */
void diquarkhash(int qk2diqk[6][6], int diqk2qk1[15], int diqk2qk2[15])
{

  /* Diquark hashing tables */
  /* Convert from 2 color-spin indices to diquark index */

  qk2diqk[0][1] =  0;
  qk2diqk[0][2] =  1;
  qk2diqk[0][3] =  2;
  qk2diqk[0][4] =  3;
  qk2diqk[0][5] =  4;
  qk2diqk[1][2] =  5;
  qk2diqk[1][3] =  6;
  qk2diqk[1][4] =  7;
  qk2diqk[1][5] =  8;
  qk2diqk[2][3] =  9;
  qk2diqk[2][4] = 10;
  qk2diqk[2][5] = 11;
  qk2diqk[3][4] = 12;
  qk2diqk[3][5] = 13;
  qk2diqk[4][5] = 14;
 
  /* Convert from diquark index to color-spin index */

  diqk2qk1[ 0] = 0;
  diqk2qk1[ 1] = 0;
  diqk2qk1[ 2] = 0;
  diqk2qk1[ 3] = 0;
  diqk2qk1[ 4] = 0;
  diqk2qk1[ 5] = 1;
  diqk2qk1[ 6] = 1;
  diqk2qk1[ 7] = 1;
  diqk2qk1[ 8] = 1;
  diqk2qk1[ 9] = 2;
  diqk2qk1[10] = 2;
  diqk2qk1[11] = 2;
  diqk2qk1[12] = 3;
  diqk2qk1[13] = 3;
  diqk2qk1[14] = 4;

  diqk2qk2[ 0] = 1;
  diqk2qk2[ 1] = 2;
  diqk2qk2[ 2] = 3;
  diqk2qk2[ 3] = 4;
  diqk2qk2[ 4] = 5;
  diqk2qk2[ 5] = 2;
  diqk2qk2[ 6] = 3;
  diqk2qk2[ 7] = 4;
  diqk2qk2[ 8] = 5;
  diqk2qk2[ 9] = 3;
  diqk2qk2[10] = 4;
  diqk2qk2[11] = 5;
  diqk2qk2[12] = 4;
  diqk2qk2[13] = 5;
  diqk2qk2[14] = 5;
}

/* Create diquark propagator from Pauli quark propagator */
void diquarkprop(field_offset qk, field_offset diqk) 
{
  int dqi,dqj;
  int cs1,cs2,cs3,cs4,c1,s1,c2,s2,c3,s3,c4,s4;
  
  register int i;
  register site *s;

  complex diquark_temp;
  int qk2diqk[6][6], diqk2qk1[15], diqk2qk2[15];

  /* Construct hash tables for diquark indices */

  diquarkhash(qk2diqk, diqk2qk1, diqk2qk2);

  FORALLSITES(i,s)
    {
      for(dqi = 0; dqi < 15; dqi++) for(dqj = 0; dqj < 15; dqj++)
	{
	  /* Decode diquark index to get color and spin for 4 quarks */

	  cs1 = diqk2qk1[dqi];
	  cs2 = diqk2qk2[dqi];
	  cs3 = diqk2qk1[dqj];
	  cs4 = diqk2qk2[dqj];

	  /* Six-value color-spin index is color + 3*spin */

	  c1 = cs1 % 3;	  s1 = cs1/3;
	  c2 = cs2 % 3;	  s2 = cs2/3;
	  c3 = cs3 % 3;	  s3 = cs3/3;
	  c4 = cs4 % 3;	  s4 = cs4/3;

	  /* Antisymmetric product of two nonrelativistic quark propagators */

	  CMUL( 
	       ((pauli_propagator *)F_PT(s,qk))->c[c1].p[s1].p[s3].c[c3],
	       ((pauli_propagator *)F_PT(s,qk))->c[c2].p[s2].p[s4].c[c4],
	       ((dipauli_propagator *)F_PT(s,diqk))->q[dqi].q[dqj]);
	  
	  CMUL( 
	       ((pauli_propagator *)F_PT(s,qk))->c[c1].p[s1].p[s4].c[c4],
	       ((pauli_propagator *)F_PT(s,qk))->c[c2].p[s2].p[s3].c[c3],
	       diquark_temp);
	  
	  CSUB(((dipauli_propagator *)F_PT(s,diqk))->q[dqi].q[dqj],
	       diquark_temp,
	       ((dipauli_propagator *)F_PT(s,diqk))->q[dqi].q[dqj]);
	}
    }
} /* diquarkprop */
