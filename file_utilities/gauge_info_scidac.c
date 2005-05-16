/*********************** gauge_info_scidac.c *************************/
/* MIMD version 7 */

/* For utilities */

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

#include <stdio.h>
#include "../include/complex.h"
#include "../include/su3.h"
#include "../include/macros.h"
#include "../include/io_lat.h"
#include "../include/io_wprop.h"
#include <lattice.h>
#include <qio.h>

EXTERN QIO_String *xml_record_in;

/*---------------------------------------------------------------------------*/
/* This routine writes the ASCII info file.  It is called from one of
   the lattice output routines in io_lat4.c.*/

void write_appl_gauge_info(FILE *fp)
{
  /* Copy Scidac record metadata to the info file specified by fp */
  fprintf(fp,"%s",QIO_string_ptr(xml_record_in));
}

char *create_QCDML(){
  char dummy[] = "Dummy QCDML";
  char *qcdml = (char *)malloc(sizeof(dummy)+1);
  strcpy(qcdml, dummy);
  return qcdml;
}

void free_QCDML(char *qcdml){
  if(qcdml != NULL)free(qcdml);
}
