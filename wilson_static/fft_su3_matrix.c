/*  $Header: /lqcdproj/detar/cvsroot/milc_qcd/wilson_static/fft_su3_matrix.c,v 1.1 2005/02/23 00:06:10 detar Exp $
 *  A collection of routines to FFT su3 matrix
 *  objecs, such as the Wilson line, or the quark propagator
 *  which has had its spin indices stripped off.
 *
 *  This is just a driver routine, it loops
 *  over the colour indices of the matrix and FFT's each
 *  complex element.
 *
 *  This routine assumes that the set up routine
 *  has been run.
 *
 */

#include "w_static_includes.h"
#define restrict rstrict /* C-90 T3D cludge */

/* functions required in this file ****/
void setup_restrict_fourier(int *key, int *restrict ) ;

void restrict_fourier(field_offset src,field_offset space,field_offset space2,
int size,int isign) ;



/*
 *  FFT the smeared Wilson line from position
 *  to momentum space.
 *
 *
 *  smear_w_line(p,t) --->  \sum_{x}  e^{ -i p x } smear_w_line(x,t)
 *
 *  To reduce the amount of work space required the routine loops
 *  over the colour components of the su3 matrices.
 *
 */

void fft_smeared_wline()
{
  int sgn = 1 ;

  restrict_fourier(F_OFFSET(smear_w_line) , F_OFFSET(su3fftwk1), F_OFFSET(su3fftwk2), sizeof(su3_matrix), sgn);

}


/*
 *  FFT the light quark propagator stripped of its spin
 *  indices.
 *
 *  The quark propagator for negative momentum is calculated.
 *
 *  strip_quark(-p,t) --->  \sum_{x}  e^{ i p x } strip_quark(x,t)
 *
 *
 *
 *  To reduce the amount of work space required the routine loops
 *  over the colour components of the su3 matrices.
 */

void fft_negmom_quark()
{
  int sgn = -1 ;

  restrict_fourier(F_OFFSET(strip_quark) , F_OFFSET(su3fftwk1), F_OFFSET(su3fftwk2), sizeof(su3_matrix), sgn);


}


/*
 *  Set up the FFT routines, for a three dimensional
 *  FFT transfrom over the the spatial lattice. For
 *  each time slice.
 *
 * This subroutine is just a dricer for Carleton's
 * FFT routines.
 */

void setup_timeslice_fft()
{
  int key[4] = {1,1,1,0};
  int restrict[4] = {0,0,0,0};

  setup_restrict_fourier(key,restrict);

}



