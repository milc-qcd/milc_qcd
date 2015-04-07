/***************** accumulate_density.c *****************************************/

/* Read and accumulate the current density over the entire lattice */

/* MIMD version 7 */

/* 04/05/15 C. DeTar */

#include "rcorr_includes.h"

static void
sum_c_field(complex *dest, complex *src, Real wt, int count){
  int i, j;

  FORALLFIELDSITES(i){
    for(j = 0; j < count; j++){
      dest[count*i+j].real += src[count*i+j].real*wt;
      dest[count*i+j].imag += src[count*i+j].imag*wt;
    }    
  }
}

void 
accumulate_current_density(char *filename, complex *qin[], 
			   double charge, double *mass, int *count)
{
  QIO_Reader *infile;
  int status = 0;
  int jrand = 0, k;
  complex *tmp;
  QIO_String *recxml = QIO_string_create();

  QIO_verbose(QIO_VERB_OFF);

  infile = r_open_complex_scidac_file(filename, QIO_PARALLEL);
  if(infile == NULL)terminate(1);

  /* Read the lattice fields for all random sources in this file */
  /* Accumulate the result in dest */
  tmp = create_c_array_field(NMU);
  for(k = 0; k < MAXRAND; k++){
    status = read_complex_scidac_xml(infile, tmp, NMU, recxml);

    if(qio_status(status) == -1) break;  // EOF

    node0_printf("Record info \n\"%s\"\n",QIO_string_ptr(recxml));
    if(qio_status(status) != 0)exit(1);

    /* Parse metadata */
    /* Format is "source index %d mass %g" */
    sscanf(QIO_string_ptr(recxml),"%*s %*s %d %*s %lf", &jrand, mass);

    if(k != jrand){
      fprintf(stderr, "Got random source index %d but wanted %d\n", jrand, k);
      terminate(1);
    }

    /* Accumulate values in qin, weighted by the charge */
    if(qin[k] == NULL){
      qin[k] = create_c_array_field(NMU);
      if(qin[k] == NULL){
	node0_printf("accumulate_current_density: No room for qin[%d]\n",k);
	terminate(1);
      }
    }

    sum_c_field(qin[k], tmp, charge, NMU);
  }
  destroy_c_array_field(tmp, NMU);

  if(k == MAXRAND && qio_status(status) != -1){
    printf("Too many random sources for dimension %d\n", MAXRAND);
    terminate(1);
  }

  if(*count != 0 && *count != k){
    fprintf(stderr, "record count %d does not match previous count %d in %s\n", 
	    k, *count, filename);
    terminate(1);
  }
  
  *count = k;

  r_close_complex_scidac_file(infile);
}

/* accumulate_density.c */
