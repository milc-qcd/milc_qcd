/*********************** ksprop_info.c *************************/
/* MIMD version 7 */

/* For ks_imp_invert_multi */

/* Application-dependent routine for writing ksprop info file */
/* This file is an ASCII companion to the KS propagator file
   and contains information about the action used to generate it.
   This information is consistently written in the pattern

       keyword  value

   or

       keyword[n] value1 value2 ... valuen

   where n is an integer.

   To maintain a semblance of consistency, the possible keywords are
   listed in io_ksprop.h.

*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#define MAX_XML 2049

/*---------------------------------------------------------------------------*/
/* This routine writes the ASCII info file.  It is called from one of
   the lattice output routines in ../generic_ks/io_prop_ks.c. */

void write_appl_ksprop_info(FILE *fp)
{

  /* Note that the file has already been opened and
     the required magic number, time stamp, and lattice
     dimensions have already been written */

  /* The rest are optional */

}

/* Temporary until we can create XML */
/* TO DO: Proviso for a number of quark action paths different from 6 */
char *create_ks_XML()
{
  char *xml;

  xml = (char *)malloc(MAX_XML);
  strcpy(xml, "KS Propagator");
  return xml;
}


void free_ks_XML(char *xml){
  if(xml != NULL)free(xml);
}
