/**************************** copy_site_spin_wilson_vector.c ****************************/
/* MIMD version 6 */
#include "prop_form_includes.h"

#ifdef DEBUGDEF
#include "debug_form.h"
#endif



/***

     Copy a Wilson vector data structure

     dest = src 

***/


void copy_site_spin_wilson_vector(field_offset src, field_offset dest)
{
  int i ; 
  register site *s; 
  /****** --------------------------------------------------*****/


  FORALLSITES(i,s)
  {
    *((spin_wilson_vector *)F_PT(s,dest)) = *((spin_wilson_vector *)F_PT(s,src))  ; 
  }



}

