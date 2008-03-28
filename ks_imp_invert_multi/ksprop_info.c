/*********************** ksprop_info.c *************************/
/* MIMD version 7 */

/* For ks_imp_invert_multi */

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

#include "ks_imp_includes.h"

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

  write_ksprop_info_item(fp,"quark.mass","%f",(char *)&propmass,0,0);


}

/* Temporary until we can create XML */
/* TO DO: Proviso for a number of quark action paths different from 6 */
char *create_ks_XML()
{
  char *xml;
  Real *act_path_coeff = ks_act_paths.act_path_coeff;
#ifdef PERIODICBC
  char bc[] = "periodic";
#else
  char bc[] = "antiperiodic";
#endif

  xml = (char *)malloc(MAX_XML);
  
  snprintf(xml,MAX_XML,"\nDerived MILC KS field\ngauge.filename %s\npropagator.boundary_conditions space: periodic time: %s\nasqtad.u0 %7.5f\naction.path_coeff[0] %e\naction.path_coeff[1] %e\naction.path_coeff[2] %e\naction.path_coeff[3] %e\naction.path_coeff[4] %e\naction.path_coeff[5] %e\ninv_arg.rsqprop %e\n",
	   startfile,
	   bc,
	   u0,
	   act_path_coeff[0],
	   act_path_coeff[1],
	   act_path_coeff[2],
	   act_path_coeff[3],
	   act_path_coeff[4],
	   act_path_coeff[5],
	   rsqprop);

  xml[MAX_XML-1] = '\0';
  return xml;
}


void free_ks_XML(char *xml){
  if(xml != NULL)free(xml);
}
