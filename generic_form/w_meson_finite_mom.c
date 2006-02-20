/************* w_meson_finite_mom.c **************************/
/* MIMD version 7  ***/

/* Modifications
   Original version UMH April 96  w_meson.c
   cmn Nov 1996: added momentum phase factors, added ANSI protoptypes.
   C. DeTar 4/29/97 Extended source and sink meson types.
                    Switched to more general version of gamma matrix multiplication 
   C. DeTar 4/30/97 Takes a table of source and sink gamma matrices
   C. DeTar 5/24/97 Operator indexing now based on corrlist table
*/

#include "generic_form_includes.h"

/*****************************************

This is a general function that can contract two quark
propagators together to form a meson correlator.

\sum_{x} exp( i p .x } 
    trace( \Gamma_in  \gamma_5  S_1^{\dagger} \gamma_5 \Gamma_out  S_2 )

   where S_1 and S_2 are quark propagators

  Function arguments

     On input 
        src1 :: site structure pouinter to spin_wilson_vector 
        src2 :: site structure pouinter to spin_wilson_vector 
	base_pt:: partially computed offset for storage of result in prop array
	q_stride, op_stride:: stride for q index and operator index in prop array
	q_momstore :: table of momentum values to use
        no_q_values :: number of momentum values in table to use
        gamma_corr :: table of vertex and source gamma matrices to use
        no_gamma_corr :: number of entries in table to use

     On output
         prop :: complex vector to the data correlators


*******************************************/


void meson_cont_mom(complex prop[],
		    field_offset src1,field_offset src2,
		    int base_pt, int q_stride, int op_stride,
		    gamma_corr gamma_table[], int no_gamma_corr)
{
  register int i;
  register site *s; 
  
  double theta ; 
  double factx = 2.0*PI/(1.0*nx) ; 
  double facty = 2.0*PI/(1.0*ny) ; 
  double factz = 2.0*PI/(1.0*nz) ; 
  Real px,py,pz;
  complex phase_fact ; 
  
  int my_t;
  int cf, sf, si;
  int i_gamma_corr,q_pt,prop_pt;
  
  complex g1,g2;
  
  spin_wilson_vector localmat,localmat2;  /* temporary storage */
  spin_wilson_vector quark;               /* temporary storage for quark */
  spin_wilson_vector antiquark;           /* temporary storage for antiquark */
  
  
  FORALLSITES(i,s)
    {
      
      my_t = s->t;
      
      /* copy src2 into quark */
      for(si=0;si<4;si++)
	for(sf=0;sf<4;sf++)
	  for(cf=0;cf<3;cf++)
	    {
	      quark.d[si].d[sf].c[cf] = 
		((spin_wilson_vector *)F_PT(s,src2))->d[si].d[sf].c[cf];
	    }
      
      /* next, construct antiquark from src1 */
      /*first, dirac multiplication by the source gamma matrices (on left) */
      
      /*  antiquark = c.c. of quark propagator */
      for(si=0;si<4;si++)
	for(sf=0;sf<4;sf++)
	  for(cf=0;cf<3;cf++)
	    {
	      
	      CONJG(((spin_wilson_vector *)F_PT(s,src1))->d[si].d[sf].c[cf],
		    antiquark.d[sf].d[si].c[cf]); 
	      
	    }
      
      
      /* left multiply antiquark by source gamma matrices,
	 beginning with gamma_5 for quark -> antiquark */
      mult_sw_by_gamma_l( &antiquark, &localmat, G5);    
      
      /* right dirac multiplication by gamma-5 (finishing up antiquark) */
      mult_sw_by_gamma_r( &localmat, &antiquark, G5);     
      
      /* Run through the table of source-sink gamma matrices */
      for(i_gamma_corr=0; i_gamma_corr<no_gamma_corr; i_gamma_corr++)
	{
	  /* left multiply by the particular source dirac matrices */
	  /* result in localmat2 */
	  
	  mult_sw_by_gamma_l( &antiquark, &localmat, 
			     gamma_table[i_gamma_corr].gin);
	  
	  mult_sw_by_gamma_r( &localmat, &localmat2,
			     gamma_table[i_gamma_corr].gout);
	  
	  /* Run through all sink momenta */
	  for(q_pt=0; q_pt<no_q_values; q_pt++)
	    {
	      px = q_momstore[q_pt][0];
	      py = q_momstore[q_pt][1];
	      pz = q_momstore[q_pt][2];
	      
	      theta = factx*(s->x)*px + facty*(s->y)*py + factz*(s->z)*pz; 
	      phase_fact = cmplx((Real) cos(theta)  , (Real) sin(theta)) ; 
	      
	      prop_pt = my_t + base_pt + q_pt * q_stride + i_gamma_corr * op_stride;
	      
	      /* trace over propagators */
	      for(si=0;si<4;si++)
		for(sf=0;sf<4;sf++)
		  for(cf=0;cf<3;cf++)
		    {
		      g1 = localmat2.d[si].d[sf].c[cf];
		      CMUL( quark.d[sf].d[si].c[cf] , phase_fact,  g2);    
		      
		      prop[prop_pt ].real += 
			(g1.real*g2.real - g1.imag*g2.imag); 
		      prop[prop_pt ].imag += 
			(g1.real*g2.imag + g1.imag*g2.real); 
		    }
	    }
	}
    }  /**** end of the loop over lattice sites ******/
  
  


} /* end of meson_cont_mom function  */


