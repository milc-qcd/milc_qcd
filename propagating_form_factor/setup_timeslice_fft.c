/**************************** setup_timeslice_fft.c ****************************/
/* MIMD version 6 */
/*
 *
 *
 */

#include "prop_form_includes.h"

#ifdef  DEBUGDEF
#include "debug_form.h"
#endif


/*
 *  Set up the FFT routines, for a three dimensional
 *  FFT transfrom over the the spatial lattice. For
 *  each time slice.
 *
 * This subroutine is just a driver for MILC's
 * FFT routines.
 */

void setup_timeslice_fft()
{
  int key[4] = {1,1,1,0};
  int l_restrict[4] = {0,0,0,0};

  setup_restrict_fourier(key,l_restrict);

}



