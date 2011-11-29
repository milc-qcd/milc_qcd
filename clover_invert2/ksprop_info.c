/*********************** ksprop_info.c *************************/
/* MIMD version 7 */

/* For clover_invert2 */

/* NEEDS WORK: PLEASE REWRITE THIS TO CONFORM TO clov_info.c */

/* Application-dependent routine for writing ksprop info file */
/* This file is an ASCII companion to the KS propagator file
   and contains information about the action used to generate it.
   This information is consistently written in the pattern

       keyword  value

   or

       keyword[n] value1 value2 ... valuen

   where n is an integer.

   To maintain a semblance of consistency, the possible keywords are
   listed in io_ksprop.h. 

   */

#include "cl_inv_includes.h"

#define MAX_XML 2049

/*---------------------------------------------------------------------------*/
/* This routine writes the ASCII info file.  It is called from one of
   the lattice output routines in ../generic_ks/io_prop_ks.c.*/

void write_appl_ksprop_info(FILE *fp)
{

  /* Note that the file has already been opened and
     the required magic number, time stamp, and lattice
     dimensions have already been written */

  /* The rest are optional */

  write_ksprop_info_item(fp,"quark.mass","%f",(char *)&mass,0,0);


}

/* Temporary until we can create XML */
/* TO DO: Proviso for a number of quark action paths different from 6 */
char *create_ks_XML(void)
{
  char *xml;
  char *ac_str = get_action_parameter_string(fn_links);
  char bc[] = "antiperiodic";

  xml = (char *)malloc(MAX_XML);
  
  snprintf(xml,MAX_XML,"\nDerived MILC KS field\ngauge.filename %s\npropagator.boundary_conditions space: periodic time: %s\nasqtad.u0 %7.5f\n%s\ninv_arg.rsqprop %e\n",
	   param.startfile,
	   bc,
	   u0,
	   ac_str,
	   rsqprop);

  xml[MAX_XML-1] = '\0';
  return xml;
}


void free_ks_XML(char *xml){
  if(xml != NULL)free(xml);
}
