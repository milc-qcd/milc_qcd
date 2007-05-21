/*********************** wilson_info_dummy.c *************************/
/* MIMD version 7 */

/* For utilities */

/* Application-dependent routine for writing gauge info file
   called from one of the output routines in io_prop_w.c */

/* This file is an ASCII companion to the gauge configuration file
   and contains information about the action used to generate it.
   This information is consistently written in the pattern

       keyword  value

   To maintain a semblance of consistency, the possible keywords are
   listed in io_wprop.h.  Add more as the need arises, but be sure
   to notify the rest of the collaboration.

   */

#include <stdio.h>
#include "../include/complex.h"
#include "../include/su3.h"
#include "../include/macros.h"
#include "../include/io_lat.h"
#include "../include/io_wprop.h"

/* build_w_prop_hdr        Fills in the spin table of contents in the header
                           structure
   write_appl_w_prop_info  Writes supplementary information to the info file */

/*---------------------------------------------------------------------------*/

/* Fill in the spin table of contents for the propagator header -
   In some projects we may want to write propagators for only
   a couple of source spins, rather than the complete set of 4.
   This table of contents specifies which spin values actually
   appear. */

void build_w_prop_hdr(w_prop_header *wph)
{
} /* build_w_prop_hdr */

/*---------------------------------------------------------------------------*/
/* This routine writes the ASCII info file.  It is called from one of
   the lattice output routines in io_prop_w.c.*/

void write_appl_w_prop_info(FILE *fp)
{
}


