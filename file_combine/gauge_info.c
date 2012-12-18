/*********************** gauge_info.c *************************/
/* MIMD version 7 */

/* For clover_invert */

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

#include "ext_src_includes.h"
#include <string.h>

/*---------------------------------------------------------------------------*/
/* This routine writes the ASCII info file.  It is called from one of
   the lattice output routines in io_lat4.c.*/

void write_appl_gauge_info(FILE *fp, gauge_file *gf)
{
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

  /* The rest are optional */
  snprintf(info+bytes, max-bytes,"%s",end_info);
  bytes = strlen(info);

  snprintf(info+bytes, max-bytes,"%s",end);
  return info;
}

void free_QCDML(char *info){
  if(info != NULL)free(info);
}


