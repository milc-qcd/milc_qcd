/*********************** gauge_info.c *************************/
/* MIMD version 6 */

/* For pure_gauge */

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

#include "pure_gauge_includes.h"

/*---------------------------------------------------------------------------*/
/* This routine writes the ASCII info file.  It is called from one of
   the lattice output routines in io_lat4.c.*/

void write_appl_gauge_info(FILE *fp)
{

  /* Note that the file has already been opened and
     the required magic number, time stamp, and lattice
     dimensions have already been written */

  /* The rest are optional */
  write_gauge_info_item(fp,"action.description","%s",
			"\"Pure gauge\"",0,0);
  write_gauge_info_item(fp,"gauge.description","%s",
			"\"One plaquette gauge action.\"",0,0);
  write_gauge_info_item(fp,"gauge.beta11","%f",(char *)&beta,0,0);
}

#define INFOSTRING_MAX 2048

char *create_QCDML(){

  size_t bytes = 0;
  char *info = (char *)malloc(INFOSTRING_MAX);
  size_t max = INFOSTRING_MAX;
  char begin[] = "<?xml version=\"1.0\" encoding=\"UTF-8\"?><info>";
  char end[] = "</info>";

  snprintf(info+bytes, max-bytes,"%s",begin);
  bytes = strlen(info);

  sprint_gauge_info_item(info+bytes, max-bytes,"action.description","%s",
			"\"Pure gauge\"",0,0);
  bytes = strlen(info);
  sprint_gauge_info_item(info+bytes, max-bytes,"gauge.description","%s",
			"\"One plaquette gauge action.\"",0,0);
  bytes = strlen(info);
  sprint_gauge_info_item(info+bytes, max-bytes,"gauge.beta11","%f",
			 (char *)&beta,0,0);

  sprint_gauge_info_item(info+bytes, max-bytes,"gauge.ssplaq","%f",
			 (char *)&g_ssplaq,0,0);
  bytes = strlen(info);
  sprint_gauge_info_item(info+bytes, max-bytes,"gauge.stplaq","%f",
			 (char *)&g_stplaq,0,0);
  bytes = strlen(info);
  sprint_gauge_info_item(info+bytes, max-bytes,"gauge.linktr.real","%f",
			 (char *)&(linktrsum.real),0,0);
  bytes = strlen(info);
  sprint_gauge_info_item(info+bytes, max-bytes,"gauge.linktr.imag","%f",
			 (char *)&(linktrsum.imag),0,0);
  bytes = strlen(info);
  snprintf(info+bytes, max-bytes,"%s",end);
  return info;
}

void free_QCDML(char *info){
  if(info != NULL)free(info);
}
