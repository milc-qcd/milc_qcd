/*
 *  FFT a smearing function, the results is written to standard
 *  output. This is really a debug test.
 *
 *  THis function assumes a cubic lattice
 *
 *   array  ::  thre dimensional array pf real data
 *   dim :: linear dimension of the latticex
 */


#include<stdio.h>
#include<math.h>


void fft_smearing_func(Real *array, int dim)
{
  const int sgn = 1 ;
  const Real pi =  3.14159265359   ;
  double fft_array_re , fft_array_im ;
  int where ;
  int nx,ny,nz ;
  int px,py,pz ;
  double phase ;



  printf("src(p) = sum_{x} exp( %d i p x ) src(x) \n",sgn);
  printf("Linear dimension of the box = %d \n",dim);

  for(px = 0 ; px < dim ;  ++px )
    for(py = 0 ; py < dim ;  ++py )
      for(pz = 0 ; pz < dim ;  ++pz )
      {

	fft_array_re  = 0.0 ;
	fft_array_im  = 0.0 ;

	for(nx = 0 ; nx < dim ;  ++nx )
	  for(ny = 0 ; ny < dim ;  ++ny )
	    for(nz = 0 ; nz < dim ;  ++nz )
	    {

	      phase = 2.0*pi*(nx*px + ny*py + nz*pz )/(dim);
	      where = nx + dim*(ny + dim*nz  );

	      fft_array_re += cos(phase)* (*(array + where)); 
	      fft_array_im += sgn*sin(phase)* (*(array + where)); 


	    }
	where = px + dim*(py + dim*pz  );

	printf("f( %d %d %d) org = %g fft=  %g %g  \n",px,py,pz,
	       *(array + where) , fft_array_re,fft_array_im);


      } /*** end of the loop over momentum ****/



}
