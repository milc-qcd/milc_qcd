/*********************** gauge_info.c *************************/
/* MIMD version 6 */

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
#include <quark_action.h>
/*---------------------------------------------------------------------------*/
/* This routine writes the ASCII info file.  It is called from one of
   the lattice output routines in io_lat4.c.*/

extern char gauge_action_description[128]; /* in gauge_stuff.c */
extern int gauge_action_nloops,gauge_action_nreps;
void write_appl_gauge_info(FILE *fp)
{

  /* Note that the file has already been opened and
     the required magic number, time stamp, and lattice
     dimensions have already been written */

  /* The rest are optional */

  write_gauge_info_item(fp,"action.description","%s",
			"\"Gauge plus fermion (improved)\"",0,0);

  write_gauge_info_item(fp,"gauge.description","%s",
			gauge_action_description,0,0);
  write_gauge_info_item(fp,"gauge.nloops","%d",(char *)&gauge_action_nloops,0,0);
  write_gauge_info_item(fp,"gauge.nreps","%d",(char *)&gauge_action_nreps,0,0);
  write_gauge_info_item(fp,"gauge.beta11","%f",(char *)&beta,0,0);
  write_gauge_info_item(fp,"gauge.tadpole.u0","%f",(char *)&u0,0,0);

  write_gauge_info_item(fp,"quark.description","%s",quark_action_description,0,0);
  write_gauge_info_item(fp,"quark.flavors","%d",(char *)&nflavors,0,0);
  write_gauge_info_item(fp,"quark.mass","%f",(char *)&mass,0,0);

}
