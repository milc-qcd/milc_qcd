/**************************** smear_hqet_prop.c ****************************/
/* MIMD version 6 */
/*
 *   $Header: /lqcdproj/detar/cvsroot/milc_qcd/hqet_heavy_to_light/smear_hqet_prop.c,v 1.1 2005/02/23 00:05:07 detar Exp $
 */

#include "hqet_light_includes.h"

#ifdef DEBUGDEF
#include DEBUGDEF
#endif


/*
 *  Initialize hqet propagator -- set up
 *
 *  prop(t= 0 , x)  =  f(x)
 *
 *  where f is the smearing function.
 *
 *
 *  In my notes the the source should be multipled by
 *  the inverse of the time component of the velocity.
 *  However, as this is just a constant factor it is more
 *  convent to pull it of the fron of the correlators.
 *
 */

void smear_hqet_prop(field_offset hqet_prop,int which_smear, int vpt)
{
  int i,j,k ;
  register site *s;
  su3_matrix tmp ;


    FORALLSITES(i,s)
    {
      if ( s->t == 0)
      {
	/*** First time slice ***************/
	for(k=0; k < 3 ; ++k)
	{
	  for(j=k+1; j < 3 ; ++j)
	  {
	    tmp.e[k][j].real = 0.0 ;
	    tmp.e[k][j].imag = 0.0 ;
	    
	    tmp.e[j][k].real = 0.0 ;
	    tmp.e[j][k].imag = 0.0 ;
	  }
	  tmp.e[k][k] = s->seq_smear_func[which_smear] ;

	}
	su3mat_copy(&tmp, (su3_matrix *) F_PT(s,hqet_prop) );

      } /**** end of the first time slice ****************/
      else
      {
	for(k=0; k < 3 ; ++k)
	  for(j=0; j < 3 ; ++j)
	  {
	    tmp.e[k][j].real = 0.0 ;
	    tmp.e[k][j].imag = 0.0 ;
	  }

	su3mat_copy(&tmp, (su3_matrix *) F_PT(s,hqet_prop) );

      }  /**** end of middle of the time  ****************/

    } /** end of the loop over the sites ****/

}

