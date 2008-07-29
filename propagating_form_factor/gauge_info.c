/*********************** gauge_info.c *************************/
/* MIMD version 6 */

/* For propagating_form_factor */

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

#include "prop_form_includes.h"

/*---------------------------------------------------------------------------*/
/* This routine writes the ASCII info file.  It is called from one of
   the lattice output routines in io_lat4.c.*/

void write_appl_gauge_info(FILE *fp, gauge_file *gf)
{
  Real gauge_fix_tol = GAUGE_FIX_TOL;

  /* Write generic information */
  write_generic_gauge_info(fp, gf);

  /* The rest are optional */
  if(fixflag==COULOMB_GAUGE_FIX)
    {
      write_gauge_info_item(fp,"gauge.fix.description","%s",
			     "\"Coulomb\"",0,0);
      write_gauge_info_item(fp,"gauge.fix.tolerance","%g",
			     (char *)&gauge_fix_tol,0,0);
    }
  

}

