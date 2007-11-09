/********************* canopy2weyl_rot.c *********************************/
/* MIMD version 7 */

/* Convert propagator matrix between gamma matrix conventions 
   for MILC (DeGrand-Rossi), FNAL, and Weyl  */

#include "generic_wilson_includes.h"

//V G V+

/* 
     V  = 
           0  1  0 -i
          -1  0  i  0
           0  1  0  i
          -1  0 -i  0

*/

/* Convert a single propagator */
void canopy2weyl(wilson_propagator *src, wilson_propagator *dest)
{
  int s0,c0,c1;
  wilson_propagator qtmp;

  for(c0=0;c0<3;c0++)for(c1=0;c1<3;c1++)
    {
      for(s0=0;s0<4;s0++){
	qtmp.c[c0].d[0].d[s0].c[c1].real =
	  src->c[c0].d[1].d[s0].c[c1].real+
	  src->c[c0].d[3].d[s0].c[c1].imag;
	
	qtmp.c[c0].d[0].d[s0].c[c1].imag =
	  src->c[c0].d[1].d[s0].c[c1].imag-
	  src->c[c0].d[3].d[s0].c[c1].real;
	
	qtmp.c[c0].d[1].d[s0].c[c1].real =
	  -src->c[c0].d[0].d[s0].c[c1].real-
	  src->c[c0].d[2].d[s0].c[c1].imag;
	
	qtmp.c[c0].d[1].d[s0].c[c1].imag =
	  -src->c[c0].d[0].d[s0].c[c1].imag+
	  src->c[c0].d[2].d[s0].c[c1].real;
	
	qtmp.c[c0].d[2].d[s0].c[c1].real =
	  src->c[c0].d[1].d[s0].c[c1].real-
	  src->c[c0].d[3].d[s0].c[c1].imag;
	
	qtmp.c[c0].d[2].d[s0].c[c1].imag =
	  src->c[c0].d[1].d[s0].c[c1].imag+
	  src->c[c0].d[3].d[s0].c[c1].real;
	
	
	qtmp.c[c0].d[3].d[s0].c[c1].real =
	  -src->c[c0].d[0].d[s0].c[c1].real+
	  src->c[c0].d[2].d[s0].c[c1].imag;
	
	qtmp.c[c0].d[3].d[s0].c[c1].imag =
	  -src->c[c0].d[0].d[s0].c[c1].imag-
	  src->c[c0].d[2].d[s0].c[c1].real;
      }
      
      for(s0=0;s0<4;s0++){
	dest->c[c0].d[s0].d[0].c[c1].real =
	  0.5*(qtmp.c[c0].d[s0].d[1].c[c1].real-
	       qtmp.c[c0].d[s0].d[3].c[c1].imag);
	
	dest->c[c0].d[s0].d[0].c[c1].imag =
	  0.5*(qtmp.c[c0].d[s0].d[1].c[c1].imag+
	       qtmp.c[c0].d[s0].d[3].c[c1].real);
	
	dest->c[c0].d[s0].d[1].c[c1].real =
	  0.5*(-qtmp.c[c0].d[s0].d[0].c[c1].real+
	       qtmp.c[c0].d[s0].d[2].c[c1].imag);
	
	dest->c[c0].d[s0].d[1].c[c1].imag =
	  0.5*(-qtmp.c[c0].d[s0].d[0].c[c1].imag-
	       qtmp.c[c0].d[s0].d[2].c[c1].real);
	
	dest->c[c0].d[s0].d[2].c[c1].real =
	  0.5*(qtmp.c[c0].d[s0].d[1].c[c1].real+
	       qtmp.c[c0].d[s0].d[3].c[c1].imag);
	
	
	dest->c[c0].d[s0].d[2].c[c1].imag =
	  0.5*(qtmp.c[c0].d[s0].d[1].c[c1].imag-
	       qtmp.c[c0].d[s0].d[3].c[c1].real);
	
	
	dest->c[c0].d[s0].d[3].c[c1].real =
	  0.5*(-qtmp.c[c0].d[s0].d[0].c[c1].real-
	       qtmp.c[c0].d[s0].d[2].c[c1].imag);
	
	dest->c[c0].d[s0].d[3].c[c1].imag =
	  0.5*(-qtmp.c[c0].d[s0].d[0].c[c1].imag+
	       qtmp.c[c0].d[s0].d[2].c[c1].real);
      }
    }
}

void canopy2weyl_site(field_offset src, field_offset dest)
{
  int i; site *s;
  wilson_propagator *a,*b;

  FORALLSITES(i,s)
    {
      a = (wilson_propagator *)F_PT(s,src);
      b = (wilson_propagator *)F_PT(s,dest);
      canopy2weyl(a, b);
    }

}
void canopy2weyl_field(wilson_propagator *src, wilson_propagator *dest)
{
  int i; site *s;

  FORALLSITES(i,s)
    {
      canopy2weyl(&src[i], &dest[i]);
    }

}

//V+ G V

void weyl2canopy(wilson_propagator *src, wilson_propagator *dest)
{
  int s0,c0,c1;
  wilson_propagator qtmp;

  for(c0=0;c0<3;c0++)for(c1=0;c1<3;c1++)
    {
      for(s0=0;s0<4;s0++){
	qtmp.c[c0].d[0].d[s0].c[c1].real =
	  -1.*(src->c[c0].d[1].d[s0].c[c1].real+
	       src->c[c0].d[3].d[s0].c[c1].real);
	
	qtmp.c[c0].d[0].d[s0].c[c1].imag =
	  -1.*(src->c[c0].d[1].d[s0].c[c1].imag+
	       src->c[c0].d[3].d[s0].c[c1].imag);
	
	qtmp.c[c0].d[1].d[s0].c[c1].real =
	  src->c[c0].d[0].d[s0].c[c1].real+
	  src->c[c0].d[2].d[s0].c[c1].real;
	
	qtmp.c[c0].d[1].d[s0].c[c1].imag =
	  src->c[c0].d[0].d[s0].c[c1].imag+
	  src->c[c0].d[2].d[s0].c[c1].imag;
	
	qtmp.c[c0].d[2].d[s0].c[c1].real =
	  src->c[c0].d[1].d[s0].c[c1].imag-
	  src->c[c0].d[3].d[s0].c[c1].imag;
	
	qtmp.c[c0].d[2].d[s0].c[c1].imag =
	  -src->c[c0].d[1].d[s0].c[c1].real+
	  src->c[c0].d[3].d[s0].c[c1].real;
	
	
	qtmp.c[c0].d[3].d[s0].c[c1].real =
	  -src->c[c0].d[0].d[s0].c[c1].imag+
	  src->c[c0].d[2].d[s0].c[c1].imag;
	
	qtmp.c[c0].d[3].d[s0].c[c1].imag =
	  src->c[c0].d[0].d[s0].c[c1].real-
	  src->c[c0].d[2].d[s0].c[c1].real;
      }
      
      for(s0=0;s0<4;s0++){
	dest->c[c0].d[s0].d[0].c[c1].real =
	  0.5*(-qtmp.c[c0].d[s0].d[1].c[c1].real-
	       qtmp.c[c0].d[s0].d[3].c[c1].real);
	
	dest->c[c0].d[s0].d[0].c[c1].imag =
	  0.5*(-qtmp.c[c0].d[s0].d[1].c[c1].imag-
	       qtmp.c[c0].d[s0].d[3].c[c1].imag);
	
	dest->c[c0].d[s0].d[1].c[c1].real =
	  0.5*(qtmp.c[c0].d[s0].d[0].c[c1].real+
	       qtmp.c[c0].d[s0].d[2].c[c1].real);
	
	dest->c[c0].d[s0].d[1].c[c1].imag =
	  0.5*(qtmp.c[c0].d[s0].d[0].c[c1].imag+
	       qtmp.c[c0].d[s0].d[2].c[c1].imag);
	
	dest->c[c0].d[s0].d[2].c[c1].real =
	  0.5*(-qtmp.c[c0].d[s0].d[1].c[c1].imag+
	       qtmp.c[c0].d[s0].d[3].c[c1].imag);
	
	
	dest->c[c0].d[s0].d[2].c[c1].imag =
	  0.5*(qtmp.c[c0].d[s0].d[1].c[c1].real-
	       qtmp.c[c0].d[s0].d[3].c[c1].real);
	
	
	dest->c[c0].d[s0].d[3].c[c1].real =
	  0.5*(qtmp.c[c0].d[s0].d[0].c[c1].imag-
	       qtmp.c[c0].d[s0].d[2].c[c1].imag);
	
	dest->c[c0].d[s0].d[3].c[c1].imag =
	  0.5*(-qtmp.c[c0].d[s0].d[0].c[c1].real+
	       qtmp.c[c0].d[s0].d[2].c[c1].real);
      }
    }
}

void bj_to_w_rot(wilson_propagator *src, wilson_propagator *dest)
{
  int s0,c0,c1;
  wilson_propagator qtmp;
  
  for(c0=0;c0<3;c0++)for(c1=0;c1<3;c1++)
    {
      for(s0=0;s0<4;s0++){
	qtmp.c[c0].d[0].d[s0].c[c1].real =
	  (src->c[c0].d[1].d[s0].c[c1].imag-
	   src->c[c0].d[3].d[s0].c[c1].imag);
	
	qtmp.c[c0].d[0].d[s0].c[c1].imag =
	  (-src->c[c0].d[1].d[s0].c[c1].real+
	   src->c[c0].d[3].d[s0].c[c1].real);
	
	qtmp.c[c0].d[1].d[s0].c[c1].real =
	  -src->c[c0].d[0].d[s0].c[c1].imag+
	  src->c[c0].d[2].d[s0].c[c1].imag;
	
	qtmp.c[c0].d[1].d[s0].c[c1].imag =
	  src->c[c0].d[0].d[s0].c[c1].real-
	  src->c[c0].d[2].d[s0].c[c1].real;
	
	qtmp.c[c0].d[2].d[s0].c[c1].real =
	  src->c[c0].d[1].d[s0].c[c1].imag+
	  src->c[c0].d[3].d[s0].c[c1].imag;
	
	qtmp.c[c0].d[2].d[s0].c[c1].imag =
	  -src->c[c0].d[1].d[s0].c[c1].real-
	  src->c[c0].d[3].d[s0].c[c1].real;
	
	
	qtmp.c[c0].d[3].d[s0].c[c1].real =
	  -src->c[c0].d[0].d[s0].c[c1].imag-
	  src->c[c0].d[2].d[s0].c[c1].imag;
	
	qtmp.c[c0].d[3].d[s0].c[c1].imag =
	  src->c[c0].d[0].d[s0].c[c1].real+
	  src->c[c0].d[2].d[s0].c[c1].real;
      }
      
      for(s0=0;s0<4;s0++){
	dest->c[c0].d[s0].d[0].c[c1].real =
	  0.5*(-qtmp.c[c0].d[s0].d[1].c[c1].imag+
	       qtmp.c[c0].d[s0].d[3].c[c1].imag);
	
	dest->c[c0].d[s0].d[0].c[c1].imag =
	  0.5*( qtmp.c[c0].d[s0].d[1].c[c1].real-
		qtmp.c[c0].d[s0].d[3].c[c1].real);
	
	dest->c[c0].d[s0].d[1].c[c1].real =
	  0.5*(qtmp.c[c0].d[s0].d[0].c[c1].imag-
	       qtmp.c[c0].d[s0].d[2].c[c1].imag);
	
	dest->c[c0].d[s0].d[1].c[c1].imag =
	  0.5*(-qtmp.c[c0].d[s0].d[0].c[c1].real+
	       qtmp.c[c0].d[s0].d[2].c[c1].real);
	
	dest->c[c0].d[s0].d[2].c[c1].real =
	  0.5*(-qtmp.c[c0].d[s0].d[1].c[c1].imag-
	       qtmp.c[c0].d[s0].d[3].c[c1].imag);
	
	dest->c[c0].d[s0].d[2].c[c1].imag =
	  0.5*(qtmp.c[c0].d[s0].d[1].c[c1].real+
	       qtmp.c[c0].d[s0].d[3].c[c1].real);
	
	dest->c[c0].d[s0].d[3].c[c1].real =
	  0.5*(qtmp.c[c0].d[s0].d[0].c[c1].imag+
	       qtmp.c[c0].d[s0].d[2].c[c1].imag);
	
	dest->c[c0].d[s0].d[3].c[c1].imag =
	  0.5*(-qtmp.c[c0].d[s0].d[0].c[c1].real-
	       qtmp.c[c0].d[s0].d[2].c[c1].real);
      }
    }
}


void weyl2canopy_site(field_offset src, field_offset dest)
{
  int i; site *s;
  wilson_propagator *a,*b;

  FORALLSITES(i,s)
    {
      a = (wilson_propagator *)F_PT(s,src);
      b = (wilson_propagator *)F_PT(s,dest);
      weyl2canopy(a, b);
    }

}
void weyl2canopy_field(wilson_propagator *src, wilson_propagator *dest)
{
  int i; site *s;

  FORALLSITES(i,s)
    {
      weyl2canopy(&src[i], &dest[i]);
    }

}

/* Convert a propagator from FNAL to MILC conventions */
/* Conversion is done in-place in the site structure */

void convert_wprop_fnal_to_milc_site(field_offset wprop)
{
  int i, color;
  site *s;
  spin_wilson_vector tmp;
  wilson_propagator *w;

  /* Rotate FNAL to MILC format  V g0 G g0 V^+ */
  
  FORALLSITES(i,s)
    {	 
      for(color=0;color<3;color++){
	w = (wilson_propagator *)F_PT(s,wprop);
	mult_swv_by_gamma_l( &(w->c[color]), &tmp, TUP );
	mult_swv_by_gamma_r( &tmp, &(w->c[color]), TUP );
      }
    }

  canopy2weyl_site( wprop, wprop);
}

/* Convert a propagator from FNAL to MILC conventions */
/* Conversion is done in-place in the specified field */

void convert_wprop_fnal_to_milc_field(wilson_propagator *wprop)
{
  int i, color;
  site *s;
  spin_wilson_vector tmp;
  wilson_propagator *w;

  /* Rotate FNAL to MILC format  V g0 G g0 V^+ */
  
  FORALLSITES(i,s)
    {	 
      for(color=0;color<3;color++){
	w = &wprop[i];
	mult_swv_by_gamma_l( &(w->c[color]), &tmp, TUP );
	mult_swv_by_gamma_r( &tmp, &(w->c[color]), TUP );
      }
    }

  canopy2weyl_field( wprop, wprop);
}

/* Convert a propagator from FNAL to MILC conventions */
/* Conversion is done in-place in the site structure */

void convert_wprop_milc_to_fnal_site(field_offset wprop)
{
  int i, color;
  site *s;
  spin_wilson_vector tmp;
  wilson_propagator *w;


  /* Rotate MILC to FNAL format  g0 V^+ G V g0 */
  
  weyl2canopy_site( wprop, wprop);

  FORALLSITES(i,s)
    {	 
      for(color=0;color<3;color++){
	w = (wilson_propagator *)F_PT(s,wprop);
	mult_swv_by_gamma_l( &(w->c[color]), &tmp, TUP );
	mult_swv_by_gamma_r( &tmp, &(w->c[color]), TUP );
      }
    }
}

/* Convert a propagator from FNAL to MILC conventions */
/* Conversion is done in-place in the specified field */

void convert_wprop_milc_to_fnal_field(wilson_propagator *wprop)
{
  int i, color;
  site *s;
  spin_wilson_vector tmp;
  wilson_propagator *w;

  /* Rotate MILC to FNAL format  g0 V^+ G V g0 */
  
  weyl2canopy_field( wprop, wprop);

  FORALLSITES(i,s)
    {	 
      for(color=0;color<3;color++){
	w = &wprop[i];
	mult_swv_by_gamma_l( &(w->c[color]), &tmp, TUP );
	mult_swv_by_gamma_r( &tmp, &(w->c[color]), TUP );
      }
    }
}

/* canopy2weyl_rot.c */

/*********************************************************************/
/* Notes

The MILC Dslash operator is defined as:

Dslash * src(x) =   sum_dir=0^3 { 
    ( 1 + gamma[dir] ) * U(x,dir) * src(x+dir)
  + ( 1 - gamma[dir] ) * U_adj(x-dir,dir) * src(x-dir) 
}

where the gamma[dir] are the MILC (Weyl-Degrand-Rossi) gamma matrices.

The Dslash operator in the massive fermion paper,
Phys.Rev.D55:3933-3957,1997, uses the opposite sign convention for
gamma in Dslash:

Dslash * src(x) =  sum_dir=0^3 { 
    ( 1 - gamma[dir] ) * U(x,dir) * src(x+dir)
  + ( 1 + gamma[dir] ) * U_adj(x-dir,dir) * src(x-dir)
}

So the transformation of the propagators has to take into account both
the change of basis and the sign flip.

First, let's see how we convert from the canopy gamma matrices to the
MILC ones.  We do this:

      gamma_m = V gamma_c Vdag

where

      [  0  1  0 -i ]
  V = [ -1  0  i  0 ] * 1/sqrt(2)
      [  0  1  0  i ]
      [ -1  0 -i  0 ]

(As I understand it, the DeGrand-Rossi minus sign for gamma_MILC[2] is
 also used in canopy.  In any case this transformation reproduces 
 that sign.)

When we import an FNAL clover propagator, we do a transformation to
put it into the MILC basis.  To use matrix multiplication notation, we
need to take the transpose of the propagator matrix to put the sink
index first.  Then the propagator matrix S satisfies the Dirac
equation as in

    D transp(S) = I

in either basis. So let's start from the canopy Dirac equation with the
FNAL propagator Sc:

    Dc transp(Sc) = I

and flip the sign of the gamma matrix to prepare the MILC convention

   gamma^5_c Dc_flip gamma^5_c^dag transp(Sc) = I

where gamma^5_c is the canopy gamma^5.  But

   Dc_flip = Vdag Dm V

where Dm is the MILC Dirac matrix.  A little algebra then gives

   Dm V gamma^5_c^dag transp(Sc) gamma^5_c Vdag = I

from which we can read off the MILC propagator matrix:

   transp(Sm) = V gamma^5_c^dag transp(Sc) gamma^5_c Vdag

In the MILC code the comments in generic_wilson/canopy2weyl.c say this
transformation is

   Sm = V gamma^0_m Sc gamma^0_m Vdag

which is, in fact, identical, since

   V^* gamma^5_c = V gamma^0_m

*********************************************************************/
