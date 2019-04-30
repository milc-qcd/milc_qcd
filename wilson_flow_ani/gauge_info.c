/*********************** gauge_info.c *************************/
/* MIMD version 7 */

/* For wilson_flow */

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

/* Definitions, files, and prototypes */
#include "wilson_flow_includes.h"

/* Additional includes */
#include <string.h>
#ifdef HAVE_QIO
#include <qio.h>
#include "../include/io_scidac.h"
#endif

/*---------------------------------------------------------------------------*/
/* This routine writes the ASCII info file.  It is called from one of
   the lattice output routines in io_lat4.c.*/

void write_appl_gauge_info(FILE *fp, gauge_file *gf)
{
  /* Variable to store string combinations (checksum and flow description) */
  char strc[20];

  /* Write generic information */
  write_generic_gauge_info(fp, gf);

  /* The rest are optional */

  /* Retain information about the original configuration (if present) */
  if( startlat_p != NULL ) {
    write_gauge_info_item(fp, "gauge.previous.filename", "\"%s\"",
                          startlat_p->filename, 0, 0);
    write_gauge_info_item(fp, "gauge.previous.time_stamp", "\"%s\"",
                          startlat_p->header->time_stamp, 0, 0);
    sprintf(strc, "%x %x", startlat_p->check.sum29, startlat_p->check.sum31);
    write_gauge_info_item(fp, "gauge.previous.checksums", "\"%s\"",
                          strc, 0, 0);
  }

  /* Flow information:                                       */
  /*  the smear factor is the step size used for integration */
  sprintf(strc, "%s flow", flow_description);
  write_gauge_info_item(fp, "gauge.smear.description", "\"%s\"",
                        strc, 0, 0);
  write_gauge_info_item(fp, "gauge.smear.steps", "\"%d\"",
                        (char *)&total_steps, 0, 0);
  write_gauge_info_item(fp, "gauge.smear.factor", "\"%f\"",
                        (char *)&stepsize, 0, 0);
}

#define INFOSTRING_MAX 2048
/* Follow USQCD style for record XML */
char *create_QCDML(){

  /* For managing info string */
  size_t bytes = 0;
  char *info = (char *)malloc(INFOSTRING_MAX);
  size_t max = INFOSTRING_MAX;

  /* For xml header and USQCD generic info tags */
  char begin[] = ("<?xml version=\"1.0\" encoding=\"UTF-8\"?>"
                  "<usqcdInfo><version>1.0</version>");
  char begin_info[] = "<info>";
  char end_info[] = "</info>";
  char end[] = "</usqcdInfo>";
  
  /* For <linktr/> and <plaq/> tags */
  Real myssplaq = g_ssplaq;
  Real mystplaq = g_stplaq;
  Real nersc_linktr = linktrsum.real/3.;
  
  /* Application specific variables */
  char strc[20];

  /* Combine USQCD generic informations and tags */
  snprintf(info+bytes, max-bytes, "%s", begin);
  bytes = strlen(info);

  snprintf(info+bytes, max-bytes, "<plaq>%e</plaq>", (myssplaq+mystplaq)/6.);
  bytes = strlen(info);

  snprintf(info+bytes, max-bytes, "<linktr>%e</linktr>", nersc_linktr);
  bytes = strlen(info);

  snprintf(info+bytes, max-bytes, "%s", begin_info);
  bytes = strlen(info);

  /* Application specific information starts here */

  /* Retain information about the original configuration (if present) */ 
  if( startlat_p != NULL ) {
      sprint_gauge_info_item(info+bytes, max-bytes, "gauge.previous.filename",
                             "%s", startlat_p->filename, 0, 0);
      bytes = strlen(info);

      sprint_gauge_info_item(info+bytes, max-bytes, "gauge.previous.time_stamp",
                             "%s", startlat_p->header->time_stamp, 0, 0);
      bytes = strlen(info);

      sprintf(strc, "%x %x", startlat_p->check.sum29, startlat_p->check.sum31);
      sprint_gauge_info_item(info+bytes, max-bytes, "gauge.previous.checksums",
                             "%s", strc, 0, 0);
      bytes = strlen(info);
  }

  /* Flow information:                                       */
  /*  the smear factor is the step size used for integration */
  sprintf(strc, "%s flow", flow_description);
  sprint_gauge_info_item(info+bytes, max-bytes, "gauge.smear.description", 
                         "%s", strc, 0, 0);
  bytes = strlen(info);
  sprint_gauge_info_item(info+bytes, max-bytes, "gauge.smear.steps", 
                         "%d", (char *)&total_steps, 0, 0);
  bytes = strlen(info);
  sprint_gauge_info_item(info+bytes, max-bytes, "gauge.smear.factor", 
                         "%f", (char *)&stepsize, 0, 0);
  bytes = strlen(info);

  /* Finish generic USQCD tags */
  snprintf(info+bytes, max-bytes, "%s", end_info);
  bytes = strlen(info);

  snprintf(info+bytes, max-bytes, "%s", end);

  return info;
}

void free_QCDML(char *info)
{
  if(info != NULL) 
    free(info);
}
