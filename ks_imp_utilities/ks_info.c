/*********************** ks_info.c *************************/
/* MIMD version 7 */

/* For ks_imp_dyn */

/* Application-dependent routine for generating metadata for KS
   propagator and source files */

#include "ks_imp_includes.h"
#include <quark_action.h>
#define MAX_XML 2049

/* Temporary until we can create XML */
/* TO DO: Proviso for a number of quark action paths different from 6 */
char *create_ks_XML()
{
  char *xml;
  Real *act_path_coeff = get_quark_path_coeff();
  
  xml = (char *)malloc(MAX_XML);
  
  snprintf(xml,MAX_XML,"\nDerived MILC KS field\ngauge.filename %s\nlayout.boundary_conditions space: periodic time: antiperiodic\nasqtad.u0 %7.5f\nasqtad.path_coeff[0] %e\nasqtad.path_coeff[1] %e\nasqtad.path_coeff[2] %e\nasqtad.path_coeff[3] %e\nasqtad.path_coeff[4] %e\nasqtad.path_coeff[5] %e\nasqtad_arg.sign MILC-convention\ninv_arg.mass %g\ninv_arg.rsqprop %e\n",
	  startfile,
	  u0,
	  act_path_coeff[0],
	  act_path_coeff[1],
	  act_path_coeff[2],
	  act_path_coeff[3],
	  act_path_coeff[4],
	  act_path_coeff[5],
	  mass,
	  rsqprop);

  xml[MAX_XML-1] = '\0';
  return xml;
}


void free_ks_XML(char *xml){
  if(xml != NULL)free(xml);
}
