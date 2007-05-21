/************* twopt_sequential.c **************************/
/* MIMD version 4  **/
/*   
    Tie the sequential source together with the gamma_5
    matrix. This is one of the two point functions required
    in the extraction of the matrix element from the simulation.

 */

/* MIMD version 6 */

/* Initialize a source for the inverter */

#include "hqet_light_includes.h"

#ifdef DEBUGDEF
#include DEBUGDEF
#endif


/*** Function prototypes *******/

void calc_anti_quark(spin_wilson_vector  *anti_quark, spin_wilson_vector *quark) ;

complex trace_hqet_light(spin_wilson_vector *light, su3_matrix *hqet, int color)  ;

void zero_spin_wilson_vector( spin_wilson_vector *in ) ; 
void inc_c_add(spin_wilson_vector *out, spin_wilson_vector *inc, complex z   ) ;

void left_vel_proj(spin_wilson_vector *out, spin_wilson_vector *in,  int v_pt, int sgn) ;

void transpose_spin_wilson_vector( spin_wilson_vector *in ) ;

/*** End of function prototypes *******/



/*
 *
 *  Function arguments 
 *    On entry
 *      v_pt  :: pointer to which velocity to use.
 *      spect_pt :: pointer to which spectator quark
 *      color    :: the source color
 *      spin     :: the source spin 
 *  On return
 *      corr :: the sequantial correlators
 *      hqet_prop :: the HQET propagator. A su3_matrix lattice site pointer.
 *                   This is really just workspace in thsi routine.
 *
 *  Global variables in the site structure
 *
 *      quark_sequential   :: the quark propagator source, as a spinwilson_vector pointer
 *               to a site structure.
 *
 *
 *  This code calculates the smeared-smeared correlators, if the light quark is 
 *  already smeared at the sink.
 */

void twopt_sequential_corr(complex *corr, field_offset hqet_prop,
			   int v_pt, int spect_pt, int color)
{
  register int i;
  int j;
  register site *s; 
  wilson_vector gammafive_quark ; 
  int pt ; 
  int t ; 
  double  t_total ;
  int sgn ; 
  int tsrc , tend ; 
  int spin ; 

  Real scale =1.0/((Real) nx*ny*nz ) ; 
  double  t_smear , t_hqet ; 

  spin_wilson_vector  anti_quark ;

  spin_wilson_vector  anti_quark_proj ;
  complex z ,  z_tmp ; 

  int tcurrent = 2*nt ; /*** dummy variable ****/
  /********** ----------------------------------------*********/

  t_total = -1.0*dclock() ; 

  /**** smear the light quark propagator st the sink *****/

  t_smear = dclock() ; 


  for(spin = 0 ; spin < 4 ; ++spin)
  {
    /**
      quark(\p) = \sum_{\p} exp( -i\p \x )  quark(\x) corresponds to isign = -1 
      **/
    sgn = 1  ;
  
    restrict_fourier_site(F_OFFSET(quark_sequential.d[spin]),
			  sizeof(wilson_vector), sgn);
		   
		   
    FORALLSITES(i,s) 
    {
      z_tmp.real = s->seq_smear_func_fft[v_pt].real * scale ;
      z_tmp.imag = s->seq_smear_func_fft[v_pt].imag * scale ; 

      c_scale_wilson_vector2( (wilson_vector *) &(s->quark_sequential.d[spin] ), &z_tmp ); 
    }


    /**
      FFT the convolution back to real space.

      smeared_quark(x) = \sum_{\p} exp( ip x )  quark(p) f(-p)
      
      corresponds to isign = -1 

      **/


    sgn = -1 ;
    restrict_fourier_site( F_OFFSET(quark_sequential.d[spin]), 
			   sizeof(wilson_vector), sgn); 

  } /*** end of the loop over spin *****/


  IF_VERBOSE_ON(1)
    printf("twopt_sequential_corr::Time to smear the light quark = %g sec\n",dclock() - t_smear );



  /**** .......... generate the HQET propagator ..........********/
  t_hqet = -dclock() ; 

  /**** smear the HQET source *******/
  smear_hqet_prop(hqet_prop, v_pt , v_pt) ; 

  tsrc = 0 ; 
  tend = nt/2  ; 
  generate_hqet_prop(hqet_prop, tsrc, tend, tcurrent, v_pt, v_pt) ;

  /*** generate the backward moving hqet propagator ***/
  tsrc = nt ; 
  tend = nt/2 ;
  generate_hqet_prop_back(hqet_prop, tsrc,tend ,tcurrent, v_pt, v_pt)  ; 

  t_hqet += dclock() ; 

  IF_VERBOSE_ON(1)
    printf("twopt_sequential_corr::Time to generate the HQET propagator = %g\n", t_hqet) ; 

  /**** trace the light quark propagator and the HQET propagator together *****/


  FORALLSITES(i,s) 
  {

    /*
       As this is a pseudoscalar --> pseudoscalar transition, in
       principle there are no gamma_5's involved. Howevee the gamma_5's
       are required to get the standard action with 
       ( 1 \pm \gamma_{\mu} ) / 2 factors.

    */
    


    calc_anti_quark(&anti_quark, (spin_wilson_vector *) &(s->quark_sequential) ) ; 
    transpose_spin_wilson_vector( &anti_quark  ) ; 

    /*
       The transpose is required because the MILC convention for 
       the spin_wilsn_vector data structure is [source_spin][sink_spin],
       which is opposite to the way the propagatos are written down.

       This just be better though out in the future

    */



    if( s->t < nt/2 )
      sgn = 1 ;
    else
      sgn = -1 ;

    left_vel_proj( &anti_quark_proj   , &anti_quark , v_pt , sgn) ; 


    /***** tie the HQET and light quark propagators together *******/
    z = trace_hqet_light( &anti_quark_proj , (su3_matrix *) F_PT(s,hqet_prop)  , color)  ; 


    t = s->t ; 
    pt = TWOPT_SEQ_WHERE(t, spect_pt, v_pt) ; 

    (corr + pt )->real += z.real ; 
    (corr + pt )->imag += z.imag ; 
     
  } /**** end of the loop over the lattice ****/

  t_total  += dclock() ; 

  IF_VERBOSE_ON(1)
    printf("Time in the twopt_sequential_source routine %g sec\n", t_total); 

}




/*
 *   Trace a su3_matrix operator with a 
 *   spin_wilson_vector_object
 *
 *    trace( hqet * light ) 
 *    
 *  Subroutine arguments
 *     color :: the source color
 *     hqet  :: HQET 
 *     light ::
 *
 */ 


complex trace_hqet_light(spin_wilson_vector *light, su3_matrix *hqet, int color) 
{
  int ispin ;  
  int i , j ; 
  complex ans  , z ; 
  

  ans.real = 0.0 ; 
  ans.imag = 0.0 ; 



  for( ispin = 0 ; ispin < 4 ; ++ispin )
    for( i = 0 ; i < 3 ; ++i)  
    {
      CMUL( hqet->e[i][color ], light->d[ispin].d[ispin].c[i] , z  );   


      ans.real += z.real ; 
      ans.imag += z.imag ; 
      
    }


  return ans ;

}




/*
 *                             [dagger]
 *  anti_quark = gamma_5 * quark * gamma_5
 *
 *
 *
 */


void calc_anti_quark(spin_wilson_vector  *anti_quark, spin_wilson_vector *quark)
{
  int sf,si,cf ; 
  spin_wilson_vector localmat_A ;
  spin_wilson_vector localmat_B ;
  


  for(si=0;si<4;si++)
    for(sf=0;sf<4;sf++)
      for(cf=0;cf<3;cf++)
      {
	CONJG( quark->d[si].d[sf].c[cf] , localmat_A.d[sf].d[si].c[cf] ) ;

	/*** note that I have transposed te spin indices ***/

      }

  /* left multiply antiquark by source gamma matrices,
     beginning with gamma_5 for quark -> antiquark */
  mult_sw_by_gamma_l( &localmat_A, &localmat_B, G5);    

  /* right dirac multiplication by gamma-5 (finishing up antiquark) */
  mult_sw_by_gamma_r( &localmat_B, anti_quark, G5);     



}


/*
 * Apply the velocity projection operator to 
 * spin_wilson matrix.
 *
 *  P_[+](v) = 1/2( 1 + v_0 * \gamma_0 - i *  v_k * \gamma_k )

    P_[-](v) = 1/2( 1 - v_0 * \gamma_0 + i *  v_k * \gamma_k )

 Apply the operator to the left.


 */



void left_vel_proj(spin_wilson_vector *out, spin_wilson_vector *in,  int v_pt, int sgn)
{
  complex z ; 
  spin_wilson_vector gamma_in ; 

  zero_spin_wilson_vector(out) ; 



  /**** unit matrix contribution ****/
  z = cmplx( 0.5 , 0.0 ) ;
  inc_c_add(out, in , z) ;  

  /******   gamma_time *********/
  mult_sw_by_gamma_l( in , &gamma_in ,  TUP ); 

  z = cmplx(sgn*0.5*velocity[ v_pt ][ TUP] , 0.0 ) ; 
  inc_c_add(out, &gamma_in , z) ; 


  /******   gamma_X *********/
  mult_sw_by_gamma_l(in ,   &gamma_in , XUP );     
  z = cmplx( 0.0 , -sgn*0.5*velocity[ v_pt ][ XUP] ) ;
  inc_c_add(out, &gamma_in , z) ; 

		    
  /******   gamma_Y *********/
  mult_sw_by_gamma_l(in ,   &gamma_in , YUP );     
  z = cmplx( 0.0 , -sgn*0.5*velocity[ v_pt ][ YUP ] ) ;
  inc_c_add(out, &gamma_in , z) ; 

		    
  /******   gamma_Z *********/
  mult_sw_by_gamma_l( in ,  &gamma_in , ZUP );     
  z = cmplx( 0.0 , -sgn*0.5*velocity[ v_pt ][ ZUP ] ) ;
  inc_c_add(out,  &gamma_in , z) ; 

  /****DEBUG*******DEBUG*******DEBUG*******DEBUG*******DEBUG*****

  IF_MASTER printf("Here is IN\n") ; 
  dump_spin_wilson_vector( in )  ;

  IF_MASTER printf("Here is OUT\n") ; 
  dump_spin_wilson_vector( out )  ;

  IF_MASTER printf("Here is GAMMA_INn") ; 
  dump_spin_wilson_vector( &gamma_in )  ;

  exit(1); 
  ****DEBUG********DEBUG********DEBUG********DEBUG********DEBUG********DEBUG****/

}



/*  For spin_wilson vectors
 *  
 *  out += inc * z 
 *           with z a complex number
 */


void inc_c_add(spin_wilson_vector *out, spin_wilson_vector *inc, complex z   )
{
  int i,j,k ; 
  complex z_tmp ;
  
  for(i = 0 ; i < 4 ; ++i)
    for(j = 0 ; j < 4 ; ++j)
      for(k = 0 ; k < 3 ; ++k)
      {
	z_tmp = cmul(&z, &(inc->d[i].d[j].c[k])  ) ;

	out->d[i].d[j].c[k].real +=  z_tmp.real ; 
	out->d[i].d[j].c[k].imag +=   z_tmp.imag ; 

      }

}

/*
 *  in = 0 
 */

void zero_spin_wilson_vector( spin_wilson_vector *in )
{
  int i , j , k ;

  for(i = 0 ; i < 4 ; ++i)
    for(j = 0 ; j < 4 ; ++j)
      for(k = 0 ; k < 3 ; ++k)
      {
	in->d[i].d[j].c[k].real = 0.0 ; 
	in->d[i].d[j].c[k].imag = 0.0 ; 
      }


}







/*
 *  in =====> in^[ transpose ]
 */

void transpose_spin_wilson_vector( spin_wilson_vector *in )
{
  int i , j , k ;
  spin_wilson_vector tmp ; 

  for(i = 0 ; i < 4 ; ++i)
    for(j = 0 ; j < 4 ; ++j)
      for(k = 0 ; k < 3 ; ++k)
      {
	tmp.d[j].d[i].c[k] =  in->d[i].d[j].c[k] ;
      }


  *in = tmp ; 

}




