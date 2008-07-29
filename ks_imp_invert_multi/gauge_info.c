/*********************** gauge_info.c *************************/
/* MIMD version 7 */

/* For ks_imp_dyn */

/* Application-dependent routine for writing gauge info file */
/* This file is an ASCII companion to the gauge configuration file
   and contains information about the action used to generate it.
   This information is consistently written in the pattern

       keyword  value

   or

       keyword[n] value1 value2 ... valuen

   where n is an integer.

   To maintain a semblance of consistency, the possible keywords are
   listed in io_lat.h.  Add more as the need arises, but be sure
   to notify the rest of the collaboration.

   */

#include "ks_imp_includes.h"
#define IMP_QUARK_ACTION_INFO_ONLY
#include <quark_action.h>
#include <string.h>
/*---------------------------------------------------------------------------*/
/* This routine writes the ASCII info file.  It is called from one of
   the lattice output routines in io_lat4.c.*/

extern char gauge_action_description[128]; /* in gauge_stuff.c */
extern int gauge_action_nloops,gauge_action_nreps;
void write_appl_gauge_info(FILE *fp, gauge_file *gf)
{
  Real gauge_fix_tol = GAUGE_FIX_TOL;

  /* Write generic information */
  write_generic_gauge_info(fp, gf);

  /* The rest are optional */

  write_gauge_info_item(fp,"action.description","%s",
			"\"Gauge plus fermion\"",0,0);

  write_gauge_info_item(fp,"gauge.description","%s",
			gauge_action_description,0,0);
  write_gauge_info_item(fp,"gauge.nloops","%d",(char *)&gauge_action_nloops,0,0);
  write_gauge_info_item(fp,"gauge.nreps","%d",(char *)&gauge_action_nreps,0,0);
  write_gauge_info_item(fp,"gauge.beta11","%f",(char *)&beta,0,0);
  write_gauge_info_item(fp,"gauge.tadpole.u0","%f",(char *)&u0,0,0);

  write_gauge_info_item(fp,"quark.description","%s",QUARK_ACTION_DESCRIPTION,0,0);
  write_gauge_info_item(fp,"quark.flavors1","%d",(char *)&nflavors1,0,0);
  write_gauge_info_item(fp,"quark.flavors2","%d",(char *)&nflavors2,0,0);
  write_gauge_info_item(fp,"quark.mass1","%f",(char *)&mass1,0,0);
  write_gauge_info_item(fp,"quark.mass2","%f",(char *)&mass2,0,0);


  if(fixflag==COULOMB_GAUGE_FIX)
    {
      write_gauge_info_item(fp,"gauge.fix.description","%s",
                             "\"Coulomb\"",0,0);
      write_gauge_info_item(fp,"gauge.fix.tolerance","%g",
                             (char *)&gauge_fix_tol,0,0);
    }

}

char *create_QCDML(){
  char dummy[] = "<?xml version=\"1.0\" encoding=\"UTF-8\"?><title>Dummy QCDML</title>";
  char *qcdml = (char *)malloc(sizeof(dummy)+1);
  strcpy(qcdml,dummy);
  return qcdml;
}

void free_QCDML(char *qcdml){
  if(qcdml != NULL)free(qcdml);
}
