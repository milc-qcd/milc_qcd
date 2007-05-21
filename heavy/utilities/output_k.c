/*********************** output_k.c *************************/
/* MIMD version 6 */

/*
 * routines for  summed meson hopping expansion (i.e. fixed kappa_h) ascii 
 * input/output. 
 */

#include "w_sum_includes.h"
#include <errno.h>

/* write ascii meson propagators */
/*
 * format: version_number (int) nx ny nz nt (int) beta   kappa_light (Real) 
 * i1(int) kappa_c   (Real) kappa_h (Real) total_number of iterations (int)
 * iteration number  (int) for(t=...) for(channel=...){prop[channel]}} 
 */

/* this subroutine opens the file for writing */

FILE *(w_ascii_k_i(char *filenam, int i1))
{
  FILE *fp;

  /* node 0 does all the writing */
  if (this_node == 0)
  {
    printf(" meson output file is named %s\n", filenam);
    fp = fopen(filenam, "w");
    if (fp == NULL)
    {
      printf("Can't open file %s, error %d\n", filenam, errno);
      terminate(1);
    }
    if ((fprintf(fp, "%d\n", VERSION_NUMBER)) == EOF)
    {
      printf("Error in writing meson propagator header\n");
      terminate(1);
    }
    if ((fprintf(fp, "%d\t%d\t%d\t%d\n", nx, ny, nz, nt)) == EOF)
    {
      printf("Error in writing meson propagator header\n");
      terminate(1);
    }
    if ((fprintf(fp, "%.7e\t%.7e\t%d\n", (double) beta, (double) kappa, i1)) == EOF)
    {
      printf("Error in writing meson propagator header\n");
      terminate(1);
    }
/**    if ((fprintf(fp, "%.7e\n", (double) kappa_c, i1)) == EOF)   BUG FIXED  **/
    if ((fprintf(fp, "%.7e\n", (double) kappa_c)) == EOF)
    {
      printf("Error in writing meson propagator header\n");
      terminate(1);
    }
    if ((fprintf(fp, "%.7e\n", (double) kappa_h)) == EOF)
    {
      printf("Error in writing meson propagator header\n");
      terminate(1);
    }
    if ((fprintf(fp, "%d\n", nhop)) == EOF)
    {
      printf("Error in writing meson propagator header\n");
      terminate(1);
    }
    fflush(fp);
  }

  else fp = NULL;  /* If not node zero */

  return fp;

}



/* this subroutine writes a meson propagator */
/* only one node does anything in this program */


void w_ascii_k(FILE * fp, int iter, double *prop[])
{
  int t, channel;
  if (this_node == 0)
  {
    /* write out iteration number */
    if ((fprintf(fp, "%d\n", iter)) == EOF)
    {
      printf("Write error in w_ascii_k\n");
      terminate(1);
    }
    fflush(fp);


    for (t = 0; t < nt; t++)
    {
      for (channel = 0; channel < nchannels; channel++)
      {
	if ((fprintf(fp, "%.7E\t", prop[channel][t])) == EOF)
	{
	  printf("Write error in w_ascii_k\n");
	  terminate(1);
	}
      }
      if ((fprintf(fp, "\n")) == EOF)
      {
	printf("Write error in w_ascii_k\n");
	terminate(1);
      }
    }
    fflush(fp);
  }				/* node 0 */
}				/* w_ascii_k  */


/* this subroutine closes the file */

void w_ascii_k_f(FILE * fp, char *filenam)
{
  if (this_node == 0)
  {
      fflush(fp);
    printf("Saved ascii meson propagator in file  %s\n", filenam);
    fclose(fp);
  }
}
