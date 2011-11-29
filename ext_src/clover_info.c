/*********************** clover_info.c *************************/
/* MIMD version 7 */

/* For ext_src */

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

/* build_w_prop_hdr        Fills in the spin table of contents in the header
                           structure
   write_appl_w_prop_info  Writes supplementary information to the info file */

#include "../include/io_wprop.h"
#include "../include/generic_quark_types.h"
extern dirac_clover_param dcptmp;
extern quark_source wqstmp;
#include <string.h>


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

  wph->n_spins  = 4;
  for(i=0;i<4;i++)
    wph->spins[i] = i;

} /* build_w_prop_hdr */

/*---------------------------------------------------------------------------*/
/* This routine writes the ASCII info file.  It is called from one of
   the lattice output routines in io_prop_w.c.*/

void write_appl_w_prop_info(FILE *fp)
{
  int n_spins = 4;
  int spins[4] = {0,1,2,3};

  /* Note that the file has already been opened and
     the required magic number, time stamp, and lattice
     dimensions, and checksums have already been written */
  
  /* The rest of these are optional */
  
  write_w_prop_info_item(fp,"quark.description","%s",
			 "\"Clover\"",0,0);
  write_w_prop_info_item(fp,"quark.kappa","%f",(char *)&dcptmp.Kappa,0,0);
  write_w_prop_info_item(fp,"quark.clover.clov_c","%f",
			 (char *)&dcptmp.Clov_c,0,0);
  write_w_prop_info_item(fp,"quark.clover.u0","%f",(char *)&dcptmp.U0,0,0);
  write_w_prop_info_item(fp,"quark.boundary_condition","%s",
			 "\"periodic all directions\"",0,0);
  /* It would be better to code this in w_source, since it
     depends on how the source is defined CD */
  write_w_prop_info_item(fp,"source.description","\"%s\"",
			 wqstmp.descrp,0,0);
  write_w_prop_info_item(fp,"source.size","%f",(char *)&wqstmp.r0,0,0);
  write_w_prop_info_item(fp,"source.x","%d",(char *)&wqstmp.x0,0,0);
  write_w_prop_info_item(fp,"source.y","%d",(char *)&wqstmp.y0,0,0);
  write_w_prop_info_item(fp,"source.z","%d",(char *)&wqstmp.z0,0,0);
  write_w_prop_info_item(fp,"source.t","%d",(char *)&wqstmp.t0,0,0);

  write_w_prop_info_item(fp,"source.n_spins","%d",(char *)&n_spins,0,0);
  write_w_prop_info_item(fp,"source.spins","%d",(char *)&spins[0],n_spins,
			 sizeof(int));
}

#define INFOSTRING_MAX 2048

/* For now we simply use the MILC info */
char *create_w_QCDML(){

  size_t bytes = 0;
  char *info = (char *)malloc(INFOSTRING_MAX);
  size_t max = INFOSTRING_MAX;
  //  char begin[] = "<?xml version=\"1.0\" encoding=\"UTF-8\"?><info>";
  //  char end[] = "</info>";
  int n_spins = 4;
  int spins[4] = {0,1,2,3};
  
  //  snprintf(info+bytes, max-bytes,"%s",begin);
  //  bytes = strlen(info);
  
  sprint_w_prop_info_item(info+bytes, max-bytes,"gauge.fix.description","%s",
			  "No gauge fixing",0,0);
  bytes = strlen(info);
  sprint_w_prop_info_item(info+bytes, max-bytes,"quark.description","%s",
			  "Clover",0,0);
  bytes = strlen(info);
  sprint_w_prop_info_item(info+bytes, max-bytes,"quark.kappa",
			  "%f",(char *)&dcptmp.Kappa,0,0);
  bytes = strlen(info);
  sprint_w_prop_info_item(info+bytes, max-bytes,"quark.clover.clov_c",
			  "%f",(char *)&dcptmp.Clov_c,0,0);
  bytes = strlen(info);
  sprint_w_prop_info_item(info+bytes, max-bytes,"quark.clover.u0",
			  "%f",(char *)&dcptmp.U0,0,0);
  bytes = strlen(info);
  sprint_w_prop_info_item(info+bytes, max-bytes,"quark.boundary_condition",
			  "%s", "periodic all directions",0,0);
  /* It would be better to code this in w_source, since it
     depends on how the source is defined CD */
  bytes = strlen(info);
  sprint_w_prop_info_item(info+bytes, max-bytes,"source.description","%s",
			  wqstmp.descrp,0,0);
  bytes = strlen(info);
  sprint_w_prop_info_item(info+bytes, max-bytes,"source.size",
			  "%f",(char *)&wqstmp.r0,0,0);
  bytes = strlen(info);
  sprint_w_prop_info_item(info+bytes, max-bytes,"source.x",
			  "%d",(char *)&wqstmp.x0,0,0);
  bytes = strlen(info);
  sprint_w_prop_info_item(info+bytes, max-bytes,"source.y",
			  "%d",(char *)&wqstmp.y0,0,0);
  bytes = strlen(info);
  sprint_w_prop_info_item(info+bytes, max-bytes,"source.z",
			  "%d",(char *)&wqstmp.z0,0,0);
  bytes = strlen(info);
  sprint_w_prop_info_item(info+bytes, max-bytes,"source.t",
			  "%d",(char *)&wqstmp.t0,0,0);
  
  bytes = strlen(info);
  sprint_w_prop_info_item(info+bytes, max-bytes,"source.n_spins",
			  "%d",(char *)&n_spins,0,0);
  bytes = strlen(info);
  sprint_w_prop_info_item(info+bytes, max-bytes,"source.spins",
			  "%d",(char *)&spins[0],n_spins,
			  sizeof(int));
  
  // bytes = strlen(info);
  //  snprintf(info+bytes, max-bytes,"%s",end);

  return info;
}

void free_w_QCDML(char *info){
  if(info != NULL)free(info);
}

