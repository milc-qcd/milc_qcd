/**************************** wilson_source_routines.c ****************************/
/* MIMD version 7 */
/*
   A collection of routines to create the source for a light
   quark inversion
 */

/* MIMD version 7 */

/* Initialize a source for the inverter */
#include "generic_form_includes.h"

/*
 *  Put the source at the origin of (0,0,0,0)
 *
 *  Function arguments 
 *    src   :: wilson_vector in the site structure wich get loadeed 
 *    color :: source colour of the propagator
 *    spin  :: source spin of the propagator
 */



/*
 *  Set a wilson vector equal to a complex number.
 *  
 *
 *  Function arguments
 *    On input
 *       z  :: pointer to the complex number
 *    On output 
 *      dest
 *
 */
void load_wvec(wilson_vector *dest, complex *z, int spin, int colour)
{

  clear_wvec( dest )    ; 

  dest->d[spin].c[colour].real = z->real ; 
  dest->d[spin].c[colour].imag = z->imag ; 


}







/*
 *  Load a smaering function into the source for a quark propagator
 *  inversion.
 *
 *  dest_[colour][spin] = smearing_func ; for t = 0 
 *  dest = 0  everywhere else
 *
 *  The smearing function has a copy of the spatial
 *  smearing function on every timeslice. This is convenient for 
 *  smearing at the sink. In this function the wislon source
 *  is set equal to the smearing function on the zeroth time slice.
 *
 *  Function arguments 
 *    dest   :: wilson_vector in the site structure which gets loadeed 
 *    src    :: complex smearing function  in the site structure 
 *    color  :: source colour of the propagator
 *    spin   :: source spin of the propagator
 */

void load_wilson_source(field_offset src, field_offset dest,int color,int spin)
{
  register int i;
  register site *s; 

	
  /* Zero the src over the lattice *****/
  FORALLSITES(i,s) 
  {
    if( s->t == 0 )
      load_wvec((wilson_vector *)F_PT(s,dest), (complex *)F_PT(s,src), spin,color); 
    else
      clear_wvec((wilson_vector *)F_PT(s,dest)); 

  }



}  /*** end of the  source loading routine *****/
