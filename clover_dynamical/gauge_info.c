/*********************** gauge_info.c *************************/
/* MIMD version 6 */

/* For clover_dynamical */

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

#include "cl_dyn_includes.h"
#include <string.h>

/*---------------------------------------------------------------------------*/
/* This routine writes the ASCII info file.  It is called from one of
   the lattice output routines in io_lat4.c.*/

void write_appl_gauge_info(FILE *fp)
{
  char sums[20];
  Real gauge_fix_tol = GAUGE_FIX_TOL;

  /* Note that the file has already been opened and
     the required magic number, time stamp, and lattice
     dimensions have already been written */

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
  if(fixflag==COULOMB_GAUGE_FIX)
    {
      write_gauge_info_item(fp,"gauge.fix.description","%s",
			     "\"Coulomb\"",0,0);
      write_gauge_info_item(fp,"gauge.fix.tolerance","%g",
			     (char *)&gauge_fix_tol,0,0);
    }
  

}

#define INFOSTRING_MAX 2048
/* For now we simply use the MILC info */
char *create_QCDML(){

  size_t bytes = 0;
  char *info = (char *)malloc(INFOSTRING_MAX);
  size_t max = INFOSTRING_MAX;
  char begin[] = "<?xml version=\"1.0\" encoding=\"UTF-8\"?><info>";
  char end[] = "</info>";
  Real gauge_fix_tol = GAUGE_FIX_TOL;
  char sums[20];

  snprintf(info+bytes, max-bytes,"%s",begin);
  bytes = strlen(info);

  /* The rest are optional */
  if(startlat_p != NULL)
    {
      /* To retain some info about the original (or previous)
	 configuration */
      bytes = strlen(info);
      sprint_gauge_info_item(info+bytes, max-bytes,"gauge.previous.filename",
			     "%s", startlat_p->filename,0,0);
      bytes = strlen(info);
      sprint_gauge_info_item(info+bytes, max-bytes,"gauge.previous.time_stamp",
			     "%s", startlat_p->header->time_stamp,0,0);
      sprintf(sums,"%x %x",startlat_p->check.sum29,startlat_p->check.sum31);
      bytes = strlen(info);
      sprint_gauge_info_item(info+bytes, max-bytes,"gauge.previous.checksums",
			     "%s", sums,0,0);
    }

  if(fixflag==COULOMB_GAUGE_FIX)
    {
      bytes = strlen(info);
      sprint_gauge_info_item(info+bytes, max-bytes,"gauge.fix.description",
			     "%s", "Coulomb",0,0);
      bytes = strlen(info);
      sprint_gauge_info_item(info+bytes, max-bytes,"gauge.fix.tolerance",
			     "%g", (char *)&gauge_fix_tol,0,0);
    }

  bytes = strlen(info);
  snprintf(info+bytes, max-bytes,"%s",end);

  return info;
}

void free_QCDML(char *info){
  if(info != NULL)free(info);
}


