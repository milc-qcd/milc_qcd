/*********************** gauge_info.c *************************/
/* MIMD version 7 */

/* For schroed_cl_inv */

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

#include "schroed_cl_includes.h"

/*---------------------------------------------------------------------------*/
/* This routine writes the ASCII info file.  It is called from one of
   the lattice output routines in io_lat4.c.*/

void write_appl_gauge_info(FILE *fp, gauge_file *gf)
{

  /* Note: this application does not write gauge field configurations! */

}

#define INFOSTRING_MAX 2048
/* For now we simply use the MILC info */
char *create_QCDML(){
  char *info = NULL;
  return info;
}

void free_QCDML(char *info){
  if(info != NULL)free(info);
}
