/******************* sequential_hqet_source.c *****************************/
/* MIMD version 6 */
/*   Create a sequential source for the hqet propagator, from a light 
     quark propagator.

     THe quark propagator at the position t_f is the source for the HQET 
     propagator inversion. The light quark propagator is smeared 
     at the sink.

  The  source  for the HQET inversion is

   src(x,t) = \sum_{r} \gamma_5  f(r) S(x+r, t = tf)
            = 0   otherwise


  IMPORTANT, this code now assumes that the light quark
  propagator is already smeared at the sink.

 */

/* MIMD version 6 */

/* Initialize a source for the inverter */
#include "hqet_light_includes.h"

/*** Function prototypes *******/


void zero_zu3_matrix(su3_matrix *in) ;

void get_wilson_vector_from_hqet_prop(field_offset out, field_offset in, int spin) ;
void wilson_vector_hqet_src(field_offset out, field_offset in, int spin, int tsrc) ;

/*** End of function prototypes *******/



/*
 *
 *  Function arguments 
 *    On entry
 *      color :: source colour of the propagator
 *      tB    :: the time slice the B meson is at.
 *      ans   :: the quark propagator source, as a wilson_vector pointer
 *               to a site structure.
 *
 *    On return
 *      ans   :: wilson_vector in the site structure which will get loadeed 
 *               This code now assumes that the light quark is alraedy smeared at the sink
 *      hqet_prop :: su3_matrix in the site structure which will get over written
 *
 */

void smeared_sequential_source(field_offset ans,
			       field_offset hqet_prop, 
			       int tB, int v_pt)
{
  complex z_tmp ;
  register int i;
  int j;
  register site *s; 
  complex phase_fact ; 
  double theta ; 
  double fact = 2.0*PI/(1.0*nx) ; 

  int sink_spin ; 
  wilson_vector gammafive_quark ; 
  Real scale =1.0/((Real) nx*ny*nz ) ; 

  double t_total , t_hqet , t_smear ; 
  /********** ----------------------------------------*********/

  t_total = dclock() ; 

  /**** use the input quark propagator as the source ***/

  FORALLSITES(i,s) 
  {
    if( s->t == tB)
    {   
      mult_by_gamma((wilson_vector *)F_PT(s,ans), 
		    &gammafive_quark , GAMMAFIVE );

      *((wilson_vector *)F_PT(s,ans)) = gammafive_quark  ; 
    } 
    else
    {
      clear_wvec((wilson_vector *)F_PT(s,ans)); 
    } 
  }




  t_hqet = dclock() ; 

  /*** generate the HQET propagator over all space ****/
  for(sink_spin = 0 ; sink_spin < 4 ; ++sink_spin )
  {
    /*** create the source for the HQET inversion ********/
    wilson_vector_hqet_src(hqet_prop , ans, sink_spin, tB) ; 

    /**** generate the HQET propagator *****/
    evolve_hqet_forwards_and_backwards(hqet_prop,tB, v_pt);  


    /***** put the HQET propagator in the correct part of the propagator ***/
    get_wilson_vector_from_hqet_prop(ans , hqet_prop, sink_spin);



  } /*** end of the loop over sink_spin   ****/


 IF_VERBOSE_ON(1)
    printf("smeared_sequential_source::Time for HQET inversion = %g sec\n",dclock() - t_hqet );



  /*** multiply the (1 + i v_slash)/2  operator into the propagator ***/
  apply_hqet_proj(ans, v_pt) ;


 IF_VERBOSE_ON(1)
    printf("smeared_sequential_source::Total time for HQET sequential inversion = %g sec\n",dclock() - t_total );



}  /*** end of the HQET sequential source function *****/






/*
 *  Copy a particular spin a wilson_vector into the
 *  diagonal elements of a su3_matrix.
 *
 *   This is only on one time slice 
 *
 *   Subroutine arguments
 *     On entry
 *        in   :: site structure pointer to wilson_matrix
 *        tB :: the time slice of the B meson
 *        spin :: spin
 *
 *     On return
 *        out :: site structure pointer to su3_matrix
 *
 */

void wilson_vector_hqet_src(field_offset out, field_offset in, int spin, int tB)
{
  register int i;
  register site *s; 
  int colour ; 
  double t_start ;


  t_start = dclock() ; 


  FORALLSITES(i,s) 
  {

    /*** zero the hqet source   **/
    zero_zu3_matrix( (su3_matrix *)F_PT(s,out) );

    if( s->t == tB)
    {
      for(colour = 0 ; colour < 3 ; ++colour)
	((su3_matrix *)F_PT(s,out))->e[colour][colour] 
	  = ((wilson_vector *)F_PT(s,in))->d[spin].c[colour] ;

    }



  } /** end of the loop over lattice sites ****/



 IF_VERBOSE_ON(1)
    printf("wilson_vector_hqet_src::Total time in function = %g sec\n",dclock() - t_start) ;



}







/*
 *  Extract a Wilson vector, from a hqet propagator stored in a 
 *  su3_matrix.
 * 
 *
 *
 *   Subroutine arguments
 *     On entry
 *        in   :: site structure pointer to su3_matrix
 *        spin :: spin
 *
 *     On return
 *        out :: site structure pointer to wilson_matrix
 */

void get_wilson_vector_from_hqet_prop(field_offset out, field_offset in, int spin)
{
  register int i;
  register site *s; 
  int colour ; 
  double t_start ;

  t_start = dclock() ; 

  FORALLSITES(i,s) 
  {
    
    for(colour = 0 ; colour < 3 ; ++colour)
    {

      CADD( ((su3_matrix *)F_PT(s,in))->e[colour][ 0 ]  , 
	   ((su3_matrix *)F_PT(s,in))->e[colour][ 1 ]  , 
	   ((wilson_vector *)F_PT(s,out))->d[spin].c[colour] );

      CSUM( ((wilson_vector *)F_PT(s,out))->d[spin].c[colour] ,
	   ((su3_matrix *)F_PT(s,in))->e[colour][ 2 ]  ) ; 

    }

  } /** end of the loop over lattice sites ****/



 IF_VERBOSE_ON(1)
    printf("get_wilson_vector_from_hqet_prop::Total time in function = %g sec\n",dclock() - t_start) ;



}









/*
 *  Set an su3_matrix equal to zero
 *
 *     in = 0 
 */

void zero_zu3_matrix(su3_matrix *in)
{
  int i , j ; 

  for(i=0 ; i < 3 ; ++i)
    for(j=0 ; j < 3 ; ++j)
    {
      in->e[i][j].real = 0.0 ; 
      in->e[i][j].imag = 0.0 ; 
    }


}
