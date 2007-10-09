/*********************** ksprop_info.c *************************/
/* MIMD version 7 */

/* For ks_spect */

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

#include "w_static_includes.h"
#include <string.h>

/*---------------------------------------------------------------------------*/
/* This routine writes the ASCII info file.  It is called from one of
   the lattice output routines in ../generic_ks/io_prop_ks.c.*/

void write_appl_ksprop_info(FILE *fp)
{

  /* This application does not create KS prop files! */

}

char *create_ks_XML(){
  char dummy[] = "KS Propagator";
  char *filexml = (char *)malloc(sizeof(dummy)+1);
  strcpy(filexml, dummy);
  return filexml;
}

void free_ks_XML(char *filexml){
  if(filexml != NULL)free(filexml);
}
