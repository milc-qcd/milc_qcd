/*********************** ks_info.c *************************/
/* MIMD version 7 */

/* For ks_imp_utilities */

/* Application-dependent routine for generating metadata for KS
   propagator and source files */

#include "ks_imp_includes.h"
#include <quark_action.h>
#include "../include/io_ksprop.h"
#define MAX_XML 2049

/* Temporary until we can create XML */
/* TO DO: Proviso for a number of quark action paths different from 6 */
char *create_ks_XML()
{
  char *xml;
  char *ac_str = get_action_parameter_string(fn_links);
  
  xml = (char *)malloc(MAX_XML);
  
  snprintf(xml,MAX_XML,"\nDerived MILC KS field\ngauge.filename %s\nlayout.boundary_conditions space: periodic time: antiperiodic\nasqtad.u0 %7.5f\n%s\ninv_arg.mass %g\ninv_arg.rsqprop %e\n",
	   startfile,
	   u0,
	   ac_str,
	   mass,
	   rsqprop);

  xml[MAX_XML-1] = '\0';
  return xml;
}


void free_ks_XML(char *xml){
  if(xml != NULL)free(xml);
}

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
