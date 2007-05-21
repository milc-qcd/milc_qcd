/**************************** calc_smear_meson.c ****************************/
/* MIMD version 7 */
/*
 *  This file contains the function that calculates the 
 *  smeared meson operators used as buikding blocks in the
 *  calculation of the static B_B parameter.
 *
 */


#include "w_static_includes.h"

/*
 *  Subroutine arguments
 *
 *   meson       :: array storing the smaered meson correlators
 *   quark       :: wilson vector containg a spin-coloir source element
 *                  of the quark propagator.
 *   quark_smear :: wilson vector storing the smeared quark propagator
 *   colour      :: colour componet of the source of the quark propagator
 *   spin        :: spin componet of the source of the quark propagator
 *
 *  The spin quark propagator indices are not manipulated in this code,
 *  this is done in the Wick contraction part of the calculation
 *  (a separate program).
 *
 */

void calc_smeared_meson(complex *meson, 
field_offset quark,field_offset quark_smear,
int colour, int spin)
{
  int ismear ;
  wilson_vector prod ;
  int ic,ispin ;
  int i ;
  int t;
  register site *s ;
  wilson_vector *smear_qrk_origin ;
  Real ts,te ;
  /****..................................................**/
  ts = dclock();

  if( (smear_qrk_origin = (wilson_vector *) calloc( (size_t) nt, sizeof(wilson_vector) )  ) == NULL )
  {
    printf("ERROR: could not reserve space in calc_smeared_meson\n");
    terminate(1);
  }

  /**** complex conjugate the quark propagator ****/
  dagger_quark_prop(quark_smear , quark);  

  for(ismear=0 ; ismear < nosmear ; ++ismear)
  {

    /*** smear the light quark *************/
    for(t=0 ;t < nt ;++t)
      clear_wvec(&smear_qrk_origin[t]);

    FORALLSITES(i,s)
    {
      c_scalar_mult_add_wvec(&smear_qrk_origin[s->t], (wilson_vector *)F_PT(s,quark_smear) , 
			      &s->smear_func[ismear] ,&smear_qrk_origin[s->t]) ;
    } 


    /*** sum up the smeared quark propagator at the origin over all the nodes *******/

    g_veccomplexsum(&(smear_qrk_origin[0].d[0].c[0]), nt*12 );

     /** calculate the smeared meson operator ****/
     FORALLSITES(i,s)
     {

       if( s->x == 0 && s->y == 0 && s->z == 0 )
       {
	 mult_wilson_vec_matdag(&prod,&smear_qrk_origin[s->t],&(s->w_line));
	 
	 t = s->t ;
       
	 for(ispin = 0 ; ispin < 4 ;++ispin)
	   for(ic=0 ; ic < 3 ;++ic)
	   {
	     CADD(*(meson + MESON_WHERE ),prod.d[ispin].c[ic] , *(meson + MESON_WHERE ));
	   }
	 
       } /** end the loop over the spatial origin ***/
       
     }  /** end the loop over sites *****/

    
  } /** end of the loop over the smearing functions ***/



  free(smear_qrk_origin);


  /*** dump out some timing information *****/
  te = dclock() - ts ;
 IF_MASTER
   printf("Time to calculate the smeared mesons (for B parameter) = %g sec\n",te);

  fflush(stdout);



}

