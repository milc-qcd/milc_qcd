/********************** bparam_smear_routines.c ****************************/
/* MIMD version 7 */
/*
 *  A collection of routines to smear the light quark 
 *  propagator in the way required for the B parameter
 *  calculation.
 *
 */

#include "w_static_includes.h"

/** Function calls required by this code ****/
void conjg_wvec(wilson_vector *src,wilson_vector *dest );

/*
 * Copy the complex conjugate of a wilson vector of the quark propagator 
 * into another wilson vector.
 *
 */

void dagger_quark_prop(field_offset conj_quark, field_offset quark)
{
  register site *s;
  int i ;


    FORALLSITES(i,s)
    {
      conjg_wvec( ((wilson_vector *)F_PT(s,quark)) , ((wilson_vector *)F_PT(s,conj_quark))  );  
    }


}


/* 
 *  Complex conjugate a Wilson vector
 * 
 *   dest = complex_conjugate (src)
 */

void conjg_wvec(wilson_vector *src,wilson_vector *dest )
{
  int ic,ispin ;

  for(ic =0 ; ic < 3 ;++ic)
    for(ispin = 0 ; ispin < 4 ; ++ispin)
    {
      CONJG( src->d[ispin].c[ic] , dest->d[ispin].c[ic] );

    }


}





