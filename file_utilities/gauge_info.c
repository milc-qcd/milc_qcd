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

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "../include/complex.h"
#include "../include/su3.h"
#include <lattice.h>    /* global variables for lattice fields */
#include "../include/macros.h"
#include "../include/comdefs.h"
#include "../include/io_lat.h"
#include "../include/generic.h"
#include <string.h>


/*---------------------------------------------------------------------------*/
/* This routine writes the ASCII info file.  It is called from one of
   the lattice output routines in io_lat4.c.*/

void write_appl_gauge_info(FILE *fp, gauge_file *gf)
{

  /* Write generic information */
  write_generic_gauge_info(fp, gf);

  /* The rest are optional */
  write_gauge_info_item(fp,"action.description","%s",
			"\"Pure gauge\"",0,0);
  write_gauge_info_item(fp,"gauge.description","%s",
			"\"One plaquette gauge action.\"",0,0);
  write_gauge_info_item(fp,"gauge.beta11","%f",(char *)&beta,0,0);
}

char *create_QCDML(){
  char dummy[] = "<?xml version=\"1.0\" encoding=\"UTF-8\"?><title>Dummy QCDML</title>";
  char *qcdml = (char *)malloc(sizeof(dummy)+1);
  strcpy(qcdml, dummy);
  return qcdml;
}

void free_QCDML(char *qcdml){
  if(qcdml != NULL)free(qcdml);
}
