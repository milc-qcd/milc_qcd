/***************** accumulate_density.c *****************************************/

/* Read and accumulate the current density over the entire lattice */

/* MIMD version 7 */

/* 04/05/15 C. DeTar */

#include "rcorr_includes.h"

static void
sum_c_field(complex *dest, Real *src, Real wt, int nrand){
  int i, j;

  FORALLFIELDSITES(i){
    for(j = 0; j < nrand; j++){
      dest[nrand*i+j].real += src[nrand*i+j]*wt;
    }    
  }
}

void 
accumulate_current_density(char *filename, complex *qin[], 
			   Real charge, Real *mass, int nrand)
{
  QIO_Reader *infile;
  int status = 0;
  int jrand = 0, k;
  Real *tmp;
  QIO_String *recxml = QIO_string_create();

  QIO_verbose(QIO_VERB_OFF);

  infile = r_open_scidac_file(filename, QIO_SERIAL);
  if(infile == NULL)terminate(1);

  /* Read the lattice fields for all random sources in this file */
  /* Accumulate the result in dest */
  tmp = create_r_array_field(NMU);
  for(k = 0; k < nrand; k++){
    status = read_real_scidac_xml(infile, tmp, NMU, recxml);
    if(qio_status(status) != 0)exit(1);

    /* Parse metadata */
    /* Format is "source index %d mass %g" */
    sscanf(QIO_string_ptr(recxml),"%*s %*s %d %*s %lf", &jrand, mass);

    if(k != jrand){
      fprintf(stderr, "Got random source index %d but wanted %d\n", jrand, k);
      terminate(1);
    }

    sum_c_field(qin[k], tmp, charge, NMU);
  }
  destroy_r_array_field(tmp, NMU);

  r_close_scidac_file(infile);

} /* accumulate_density.c */
