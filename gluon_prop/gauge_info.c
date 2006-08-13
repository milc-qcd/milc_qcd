/*********************** gauge_info.c *************************/
/* MIMD version 7 */

/* For gluon_prop */

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

#include "gluon_prop_includes.h"

/*---------------------------------------------------------------------------*/
/* This routine writes the ASCII info file.  It is called from one of
   the lattice output routines in io_lat4.c.*/

void write_appl_gauge_info(FILE *fp)
{
  char sums[20];
  float gauge_fix_tol = GAUGE_FIX_TOL;

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
  write_gauge_info_item(fp,"gauge.ssplaq","%f",(char *)&g_ssplaq,0,0);
  write_gauge_info_item(fp,"gauge.stplaq","%f",(char *)&g_stplaq,0,0);
  write_gauge_info_item(fp,"gauge.linktr.real","%f",
			(char *)&(linktrsum.real),0,0);
  write_gauge_info_item(fp,"gauge.linktr.imag","%f",
			(char *)&(linktrsum.imag),0,0);
  if(fixflag==COULOMB_GAUGE_FIX)
    {
      write_gauge_info_item(fp,"gauge.fix.description","%s",
			     "\"Coulomb\"",0,0);
      write_gauge_info_item(fp,"gauge.fix.tolerance","%g",
			     (char *)&gauge_fix_tol,0,0);
    }
  if(fixflag==LANDAU_GAUGE_FIX)
    {
      write_gauge_info_item(fp,"gauge.fix.description","%s",
			     "\"Landau\"",0,0);
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
  char sums[20];
  float gauge_fix_tol = GAUGE_FIX_TOL;

  snprintf(info+bytes, max-bytes,"%s",begin);
  bytes = strlen(info);

  if(startlat_p != NULL)
    {
      sprint_gauge_info_item(info+bytes, max-bytes,"gauge.previous.filename","%s",
			     startlat_p->filename,0,0);
      bytes = strlen(info);
      sprint_gauge_info_item(info+bytes, max-bytes,"gauge.previous.time_stamp","%s",
			     startlat_p->header->time_stamp,0,0);
      bytes = strlen(info);
      sprintf(sums,"%x %x",startlat_p->check.sum29,startlat_p->check.sum31);
      sprint_gauge_info_item(info+bytes, max-bytes,"gauge.previous.checksums","%s",
			     sums,0,0);
      bytes = strlen(info);
    }
  sprint_gauge_info_item(info+bytes, max-bytes,"gauge.ssplaq","%f",
			 (char *)&g_ssplaq,0,0);
  sprint_gauge_info_item(info+bytes, max-bytes,"gauge.stplaq","%f",
			 (char *)&g_stplaq,0,0);
  bytes = strlen(info);
  sprint_gauge_info_item(info+bytes, max-bytes,"gauge.linktr.real","%f",
			 (char *)&(linktrsum.real),0,0);
  bytes = strlen(info);
  sprint_gauge_info_item(info+bytes, max-bytes,"gauge.linktr.imag","%f",
			 (char *)&(linktrsum.imag),0,0);
  
  bytes = strlen(info);

  if(fixflag==COULOMB_GAUGE_FIX)
    {
      sprint_gauge_info_item(info+bytes, max-bytes,"gauge.fix.description","%s",
			     "\"Coulomb\"",0,0);
      bytes = strlen(info);
      
      sprint_gauge_info_item(info+bytes, max-bytes,"gauge.fix.tolerance","%g",
			 (char *)&gauge_fix_tol,0,0);
      bytes = strlen(info);
    }
  if(fixflag==LANDAU_GAUGE_FIX)
    {
      sprint_gauge_info_item(info+bytes, max-bytes,"gauge.fix.description","%s",
			     "\"Landau\"",0,0);
      bytes = strlen(info);
      sprint_gauge_info_item(info+bytes, max-bytes,"gauge.fix.tolerance","%g",
			     (char *)&gauge_fix_tol,0,0);
      bytes = strlen(info);
    }
  snprintf(info+bytes, max-bytes,"%s",end);
  return info;
}

void free_QCDML(char *info){
  if(info != NULL)free(info);
}

