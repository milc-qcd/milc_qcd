/*********************** clover_info.c *************************/
/* MIMD version 7 */

/* For clover_invert */

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

#include "cl_inv_includes.h"

#include "../include/io_wprop.h"
#include "../include/generic_quark_types.h"
#include "params.h"
extern dirac_clover_param dcptmp;
extern quark_source wqstmp;
extern gauge_file *startlat_p;
extern gauge_file *savelat_p;
EXTERN  params param;           /* user input parameters */

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
  char sums[20];
  Real gauge_fix_tol = GAUGE_FIX_TOL;
  int n_spins = 4;
  int spins[4] = {0,1,2,3};

  /* Note that the file has already been opened and
     the required magic number, time stamp, and lattice
     dimensions, and checksums have already been written */
  
  /* The rest of these are optional */
  
  if(startlat_p != NULL)
    {
      write_w_prop_info_item(fp,"gauge.filename","\"%s\"",
			     startlat_p->filename,0,0);
      write_w_prop_info_item(fp,"gauge.time_stamp","\"%s\"",
			     startlat_p->header->time_stamp,0,0);
      sprintf(sums,"%x %x",startlat_p->check.sum29,startlat_p->check.sum31);
      write_w_prop_info_item(fp,"gauge.checksums","\"%s\"",sums,0,0);
    }
  if(param.fixflag==COULOMB_GAUGE_FIX)
    {
      write_w_prop_info_item(fp,"gauge.fix.description","%s",
			     "\"Coulomb\"",0,0);
      write_w_prop_info_item(fp,"gauge.fix.tolerance","%g",
			     (char *)&gauge_fix_tol,0,0);
      if(savelat_p != NULL)
	{
	  write_w_prop_info_item(fp,"gauge.fix.filename","\"%s\"",
				 savelat_p->filename,0,0);
	  write_w_prop_info_item(fp,"gauge.fix.time_stamp","\"%s\"",
				 savelat_p->header->time_stamp,0,0);
	  sprintf(sums,"%x %x",savelat_p->check.sum29,savelat_p->check.sum31);
	  write_w_prop_info_item(fp,"gauge.fix.checksums","\"%s\"",sums,0,0);
	}
    }
  else
    write_w_prop_info_item(fp,"gauge.fix.description","%s",
			   "\"No gauge fixing\"",0,0);
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
  char *info;
  size_t max = INFOSTRING_MAX;
  //  char begin[] = "<?xml version=\"1.0\" encoding=\"UTF-8\"?><info>";
  //  char end[] = "</info>";
  char sums[20];
  Real gauge_fix_tol = GAUGE_FIX_TOL;
  int n_spins = 4;
  int spins[4] = {0,1,2,3};

  info = (char *)malloc(INFOSTRING_MAX);  
  info[0] = '\0';
  //  snprintf(info+bytes, max-bytes,"%s",begin);
  //  bytes = strlen(info);
  
  if(startlat_p != NULL)
    {
      bytes = strlen(info);
      sprint_w_prop_info_item(info+bytes, max-bytes,"gauge.filename","%s",
			      startlat_p->filename,0,0);
      bytes = strlen(info);
      sprint_w_prop_info_item(info+bytes, max-bytes,"gauge.time_stamp","%s",
			      startlat_p->header->time_stamp,0,0);
      sprintf(sums,"%x %x",startlat_p->check.sum29,startlat_p->check.sum31);
      bytes = strlen(info);
      sprint_w_prop_info_item(info+bytes, max-bytes,"gauge.checksums",
			      "%s",sums,0,0);
    }
  if(param.fixflag==COULOMB_GAUGE_FIX)
    {
      bytes = strlen(info);
      sprint_w_prop_info_item(info+bytes, max-bytes,"gauge.fix.description",
			      "%s", "Coulomb",0,0);
      bytes = strlen(info);
      sprint_w_prop_info_item(info+bytes, max-bytes,"gauge.fix.tolerance","%g",
			      (char *)&gauge_fix_tol,0,0);
      if(savelat_p != NULL)
	{
	  bytes = strlen(info);
	  sprint_w_prop_info_item(info+bytes, max-bytes,"gauge.fix.filename",
				  "%s", savelat_p->filename,0,0);
	  bytes = strlen(info);
	  sprint_w_prop_info_item(info+bytes, max-bytes,"gauge.fix.time_stamp",
				  "%s", savelat_p->header->time_stamp,0,0);
	  bytes = strlen(info);
	  sprintf(sums,"%x %x",savelat_p->check.sum29,savelat_p->check.sum31);
	  bytes = strlen(info);
	  sprint_w_prop_info_item(info+bytes, max-bytes,"gauge.fix.checksums",
				  "%s",sums,0,0);
	}
    }
  else {
    bytes = strlen(info);
    sprint_w_prop_info_item(info+bytes, max-bytes,"gauge.fix.description","%s",
			    "No gauge fixing",0,0);
  }
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

