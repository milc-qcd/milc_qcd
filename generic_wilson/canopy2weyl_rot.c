#include "generic_wilson_includes.h"

//V G V+
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

  /* Rotate FNAL to Milc format  V g0 G g0 V^+ */
  
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

  /* Rotate FNAL to Milc format  V g0 G g0 V^+ */
  
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


  /* Rotate FNAL to MILC format  V g0 G g0 V^+ */
  
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

  /* Rotate FNAL to MILC format  V g0 G g0 V^+ */
  
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

