/*********************** output_m.c *************************/
/* MIMD version 7 */

/* routines for  meson hopping expansion ascii  output. */

#include "w_heavy_includes.h"
#include <errno.h>

/* write ascii meson propagators */
/*
 * format: version_number (int) nx ny nz nt (int) beta kappa_light (Real)
 * width^-2 (Real), wallflag(int) wall_cutoff (int), wall_separation (int)
 * kappa_c (Real)    total_number_of_hopping_iters  (int) for spin{ for
 * color{   for(N_iter) spin color N_iter for(t=...)
 * for(channel=...){prop[channel]}} 
 */

/* this subroutine opens the file for writing */
FILE *(w_ascii_m_i(char *filenam, int i1))
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
    if ((fprintf(fp, "%d\n", (int)VERSION_NUMBER)) == EOF)
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
    if ((fprintf(fp, "%.7e\t%d\n", (double) width, wqs.type)) == EOF)
    {
      printf("Error in writing meson propagator header\n");
      terminate(1);
    }
    if ((fprintf(fp, "%d\t%d\n", wqs.wall_cutoff, wall_separation)) == EOF)
    {
      printf("Error in writing meson propagator header\n");
      terminate(1);
    }
    if ((fprintf(fp, "%.7e\n", (double) kappa_c)) == EOF)
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

/* this subroutine opens the file for appending */
FILE *(a_ascii_m_i(char *filenam, int i1))
{
  FILE *fp;

  /* node 0 does all the writing */
  if (this_node == 0)
  {
    printf(" meson output file is named %s\n", filenam);
    fp = fopen(filenam, "a");
    if (fp == NULL)
    {
      printf("Can't open file %s, error %d\n", filenam, errno);
      terminate(1);
    }
  }

  else fp = NULL;  /* If not node zero */

  return fp;
}


/* this subroutine writes a meson propagator */
/*
 * all nodes have the same summed meson propagator at this point so node 0
 * does all the work  
 */

void w_ascii_m(FILE * fp, int spin, int color, int iter, double *prop[])
{
  int t, tp, channel;
  g_sync();
  if (this_node == 0)
  {
    /* first the spin and color and iteration_number  */
    if ((fprintf(fp, "%d %d   %d\n", spin, color, iter)) == EOF)
    {
      printf("Write error in w_ascii_m\n");
      terminate(1);
    }
    fflush(fp);

    /* next the elements */

    for (t = 0; t < nt; t++)
    {
      /*
       * t is time in regular order; but prop stores time-slices in
       * even-first order 
       */

      tp = (t % 2) * (nt / 2) + t / 2;

      for (channel = 0; channel < nchannels; channel++)
      {
	if ((fprintf(fp, "%.7E\t", prop[channel][tp])) == EOF)
	{
	  printf("Write error in w_ascii_m\n");
	  terminate(1);
	}
      }
      if ((fprintf(fp, "\n")) == EOF)
      {
	printf("Write error in w_ascii_m\n");
	terminate(1);
      }
    }
    fflush(fp);
  }				/* node 0 */
}				/* w_ascii_m  */



/* this subroutine closes the file */
void w_ascii_m_f(FILE * fp, char *filenam)
{
    g_sync();
  if (this_node == 0)
  {
    fflush(fp);
    printf("Saved ascii meson propagator in file  %s\n", filenam);
    fclose(fp);
  }
}
