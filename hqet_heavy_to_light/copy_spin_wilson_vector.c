/**************************** copy_spin_wilson_vector.c ****************************/
/* MIMD version 6 */
/*
 *  A couple of routines to copy spin wilson vectors
 *
 *
 *
 */


#include "hqet_light_includes.h"

#ifdef DEBUGDEF
#include DEBUGDEF
#endif









/*
 *  Copy spin wilson vector 
 *    out = in 
 *
 */

void copy_spin_wilson_vector(spin_wilson_vector *in , spin_wilson_vector *out) 
{
  int i , j , k ; 

  for(i = 0 ; i < 4 ; ++i)  
    for(j = 0 ; j < 4 ; ++j)  
      for(k = 0 ; k < 3 ; ++k)  
      {
	out->d[i].d[j].c[k] = in->d[i].d[j].c[k] ; 
      }



}



/*
 *  Copy a spin_wilson_vector in the site structure into 
 *  another spin_wilson_vector
 *
 *
 */


void copy_lattice_spin_wilson_vector(field_offset out, field_offset in)
{
  register int i;
  register site *s; 


  FORALLSITES(i,s) 
  {
    copy_spin_wilson_vector((spin_wilson_vector *) F_PT(s,in) , (spin_wilson_vector *) F_PT(s,out)  )  ;
  }


}
