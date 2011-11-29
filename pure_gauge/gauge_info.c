/*********************** gauge_info.c *************************/
/* MIMD version 7 */

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

void write_appl_gauge_info(FILE *fp, gauge_file *gf)
{
  Real myssplaq = g_ssplaq;  /* Precision conversion */
  Real mystplaq = g_stplaq;  /* Precision conversion */
  Real nersc_linktr = linktrsum.real/3.;  /* Convention and precision */

  /* Write generic information */
  write_generic_gauge_info(fp, gf);

  /* The rest are optional */
  write_gauge_info_item(fp,"action.description","%s",
			"\"Pure gauge\"",0,0);
  write_gauge_info_item(fp,"gauge.description","%s",
			"\"One plaquette gauge action.\"",0,0);
  write_gauge_info_item(fp,"gauge.beta11","%f",(char *)&beta,0,0);
  write_gauge_info_item(fp,"gauge.ssplaq","%f",(char *)&myssplaq,0,0);
  write_gauge_info_item(fp,"gauge.stplaq","%f",(char *)&mystplaq,0,0);
  write_gauge_info_item(fp,"gauge.nersc_linktr","%f",
			(char *)&(nersc_linktr),0,0);
  write_gauge_info_item(fp,"gauge.nersc_checksum","%u",
			(char *)&(nersc_checksum),0,0);
}

#define INFOSTRING_MAX 2048
/* Follow USQCD style for record XML */
char *create_QCDML(){

  size_t bytes = 0;
  char *info = (char *)malloc(INFOSTRING_MAX);
  size_t max = INFOSTRING_MAX;
  char begin[] = "<?xml version=\"1.0\" encoding=\"UTF-8\"?><usqcdInfo><version>1.0</version>";
  char begin_info[] = "<info>";
  char end_info[] = "</info>";
  char end[] = "</usqcdInfo>";
  Real myssplaq = g_ssplaq;  /* Precision conversion */
  Real mystplaq = g_stplaq;  /* Precision conversion */
  Real nersc_linktr = linktrsum.real/3.;  /* Convention and precision */


  snprintf(info+bytes, max-bytes,"%s",begin);
  bytes = strlen(info);

  snprintf(info+bytes, max-bytes,"<plaq>%e</plaq>",(myssplaq+mystplaq)/6.);
  bytes = strlen(info);

  snprintf(info+bytes, max-bytes,"<linktr>%e</linktr>",nersc_linktr);
  bytes = strlen(info);

  snprintf(info+bytes, max-bytes,"%s",begin_info);
  bytes = strlen(info);

  sprint_gauge_info_item(info+bytes, max-bytes,"action.description","%s",
			"\"Pure gauge\"",0,0);
  bytes = strlen(info);
  sprint_gauge_info_item(info+bytes, max-bytes,"gauge.description","%s",
			"\"One plaquette gauge action.\"",0,0);
  bytes = strlen(info);
  sprint_gauge_info_item(info+bytes, max-bytes,"gauge.beta11","%f",
			 (char *)&beta,0,0);

  bytes = strlen(info);
  sprint_gauge_info_item(info+bytes, max-bytes,"gauge.ssplaq","%f",
			 (char *)&myssplaq,0,0);
  bytes = strlen(info);
  sprint_gauge_info_item(info+bytes, max-bytes,"gauge.stplaq","%f",
			 (char *)&mystplaq,0,0);
  bytes = strlen(info);
  sprint_gauge_info_item(info+bytes, max-bytes,"gauge.nersc_linktr","%e",
			 (char *)&nersc_linktr,0,0);
  bytes = strlen(info);
  sprint_gauge_info_item(info+bytes, max-bytes,"gauge.nersc_checksum","%u",
			 (char *)&nersc_checksum,0,0);
  bytes = strlen(info);
  snprintf(info+bytes, max-bytes,"%s",end_info);

  bytes = strlen(info);
  snprintf(info+bytes, max-bytes,"%s",end);
  return info;
}

void free_QCDML(char *info){
  if(info != NULL)free(info);
}
