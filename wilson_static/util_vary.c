/**************************** util_vary.c ****************************/
/* MIMD version 7 */
/*
 * Some useful general purpose routines for the
 * static - B parameter calculation.
 */

#include "w_static_includes.h"

/*
 *  Driver routine for setting up the smearing
 *  routines.
 *
 */


void setup_vary(complex *meson, int nodata)
{
  int i ;

  /*** create the smearing functions *****/
  create_smearing_funcs() ;

  /** intialise the Fourier transform routines ****/
  setup_timeslice_fft();


  for(i=0 ; i < nodata ;++i)
  {
    (*(meson + i )).real = 0.0 ;
    (*(meson + i )).imag = 0.0 ;
  }

  /** Zero the store of the quark propagator required for the smearing matrix ***/
  zero_strip_quark();

  /* Calculate the gauge part of the static propagator ***/
  static_prop() ; 



}





/*
 *  Zero the store for the quark popagatots with its
 *  spin indices are stripped off
 *
 *
 */


void zero_strip_quark()
{
  int ic,jc ;
  int i ;
  register site *s ;


    FORALLSITES(i,s)
    {   
       for(ic = 0 ; ic < 3 ;++ic )
          for(jc = 0 ; jc < 3 ;++jc )
          {
            s->strip_quark.e[ic][jc].real = 0.0  ;
            s->strip_quark.e[ic][jc].imag = 0.0  ;
          }

    }


}


