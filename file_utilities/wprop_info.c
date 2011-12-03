/*********************** wprop_info.c *************************/
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
   listed in io_wprop.h.

*/

#include <stdlib.h>
#include <stdio.h>
#include "lattice.h"
#include "../include/io_wprop.h"
#include <string.h>
#define MAX_XML 2049

/*---------------------------------------------------------------------------*/

/* Fill in the spin table of contents for the propagator header -
   In some projects we may want to write propagators for only
   a couple of source spins, rather than the complete set of 4.
   This table of contents specifies which spin values actually
   appear. */

void build_w_prop_hdr(w_prop_header *wph)
{
  int i;

  /* Note that all other values in the header structure are
     loaded by the io_prop_w.c routines, since they are common
     to all projects */

  /* Copy from values preset in lattice_cl.h */

  wph->n_spins  = n_spins;
  for(i=0;i<n_spins;i++)
    wph->spins[i] = spins[i];

} /* build_w_prop_hdr */

/*---------------------------------------------------------------------------*/
/* This routine writes the ASCII info file.  It is called from one of
   the lattice output routines in ../generic_w/io_prop_w.c. */

void write_appl_w_prop_info(FILE *fp)
{

  /* Note that the file has already been opened and
     the required magic number, time stamp, and lattice
     dimensions have already been written */

  /* The rest are optional */

}

/* Temporary until we can create XML */
/* TO DO: Proviso for a number of quark action paths different from 6 */
char *create_w_QCDML()
{
  char *xml;

  xml = (char *)malloc(MAX_XML);
  strcpy(xml, "Wilson Propagator");
  return xml;
}


void free_w_QCDML(char *xml){
  if(xml != NULL)free(xml);
}
