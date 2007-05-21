/**************************** copy_site_wilson_vector.c ****************************/
/* MIMD version 7 */
#include "generic_wilson_includes.h"

#ifdef DEBUGDEF
  #include DEBUGDEF
#endif



/***

     Copy a Wilson vector data structure

     dest = src 

***/


void copy_site_wilson_vector(field_offset src, field_offset dest)
{
  int i ; 
  register site *s; 
  /****** --------------------------------------------------*****/


  FORALLSITES(i,s)
  {
    copy_wvec( (wilson_vector *)F_PT(s,src)  , (wilson_vector *)F_PT(s,dest)  ) ; 

  }



}

