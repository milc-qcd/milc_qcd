/*********************** gauge_info.c *************************/
/* MIMD version 6 */

/* For hvy_qpot */

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

#include "string_break_includes.h"

/*---------------------------------------------------------------------------*/
/* This routine writes the ASCII info file.  It is called from one of
   the lattice output routines in io_lat4.c.*/

void write_appl_gauge_info(FILE *fp, gauge_file *gf)
{

  char sums[20];
  Real ape_weight;

  /* Write generic information */
  write_generic_gauge_info(fp, gf);

  /* The rest are optional */

  if(startlat_p != NULL)
    {
      /* To retain some info about the original (or previous)
	 configuration */
      write_gauge_info_item(fp,"gauge.previous.filename","\"%s\"",
                             startlat_p->filename,0,0);
      write_gauge_info_item(fp,"gauge.previous.time_stamp","\"%s\"",
                             startlat_p->header->time_stamp,0,0);
      sprintf(sums,"%x %x",startlat_p->check.sum29,startlat_p->check.sum31);
      write_gauge_info_item(fp,"gauge.previous.checksums","\"%s\"",sums,0,0);
    }

  write_gauge_info_item(fp,"gauge.fix.description","\"%s\"",
			"\"Temporal axial with ploop copies\"",0,0);
  write_gauge_info_item(fp,"gauge.smear.description","\"%s\"",
			"\"Spatial links SU3 projected\"",0,0);
  write_gauge_info_item(fp,"gauge.smear.steps","\"%d\"",
			(char *)&tot_smear,0,0);

  /* For the smearing factor we use Tom and Anna's convention
     which corresponds to the weighting

       U_fat = (1 - a) U_link + a/6 U_staple

       so smear_fac = (1 - a)/(a/6)

       */
     
  ape_weight = 6./(smear_fac + 6.);

  write_gauge_info_item(fp,"gauge.smear.factor","\"%f\"",
			(char *)&ape_weight,0,0);
}
