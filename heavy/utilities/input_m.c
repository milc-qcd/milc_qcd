/*********************** input_m.c *************************/
/* MIMD version 6 */

/* routines for  meson hopping expansion ascii  input */

#include "w_sum_includes.h"
#include <errno.h>

/* read ascii meson propagators */
/*
 * format: version_number (int) nx ny nz nt (int) beta kappa_light (Real)
 * width^-2 (Real), wallflag(int) wall_cutoff (int), wall_separation (int)
 * kappa_c (Real)    fnhop (int)   (number of hopping iterations in the
 * file) for spin{ for color{   for(N_iter) spin color N_iter for(t=...)
 * for(channel=...){prop[channel]}} 
 */


/* this subroutine opens the file for reading */
/*
 * opens files, reads header, tells calling routine the number of hopping
 * iterations in file 
 */
FILE *(r_ascii_m_i(char *filenam, int i1, int *file_hops))
{
  FILE *fp;
  int version_number, t;
  Real x1, x2;
  int i_in, dummy1, dummy2;
  /* node 0 does all the reading */
  if (this_node == 0)
  {
    printf(" input file is named %s\n", filenam);
    fp = fopen(filenam, "r");
    if (fp == NULL)
    {
      printf("Can't open file %s, error %d\n", filenam, errno);
      terminate(1);
    }
    if ((fscanf(fp, "%d", &version_number)) != 1)
    {
      printf("Error in reading meson header\n");
      terminate(1);
    }
    if (version_number != VERSION_NUMBER)
    {
      printf("Incorrect version number in meson header\n");
      terminate(1);
    }
    if ((fscanf(fp, "%d%d%d%d", &nx, &ny, &nz, &t)) != 4)
    {
      printf("Error in reading meson header\n");
      terminate(1);
    }
    if (t != nt)
    {
      printf("Incorrect lattice nt %d\n", t);
      terminate(1);
    }
#if PRECISION == 1
    if ((fscanf(fp, "%e %e %d", &x1, &x2, &i_in)) != 3)
#else
    if ((fscanf(fp, "%le %le %d", &x1, &x2, &i_in)) != 3)
#endif
    {
      printf("Error in reading meson header\n");
      terminate(1);
    }
    if (fabs(x1 - beta) > EPS || fabs(x2 - kappa) > EPS)
    {
      printf(
	     "Warning: input couplings, %e, %e, not equal to program's\n",
	     (double) x1, (double) x2);
    }
    if (i_in != i1)
    {
      printf("Warning: file has %d props, prog wants %d\n", i_in, i1);
    }
#if PRECISION == 1
    if ((fscanf(fp, "%e %d", &x1, &i_in)) != 2)
#else
    if ((fscanf(fp, "%le %d", &x1, &i_in)) != 2)
#endif
    {
      printf("Error in reading meson header\n");
      terminate(1);
    }
    if ((fscanf(fp, "%d%d", &dummy1, &dummy2)) != 2)
    {
      printf("Error in reading meson header\n");
      terminate(1);
    }
#if PRECISION == 1
    if ((fscanf(fp, "%e", &x1)) != 1)
#else
    if ((fscanf(fp, "%le", &x1)) != 1)
#endif
    {
      printf("Error in reading meson header\n");
      terminate(1);
    }
    if (fabs(x1 - kappa_c) > EPS)
    {
      printf(
	     "File kappa_c, %e, not equal to program's, %e\n",
	     (double) x1, (double) kappa_c);
      terminate(1);
    }
    if ((fscanf(fp, "%d", file_hops)) != 1)
    {
      printf("Error in reading meson header\n");
      terminate(1);
    }
  }

  else fp = NULL;  /* If not node zero */

  return fp;
}

/* this subroutine reads a meson propagator */
void r_ascii_m(FILE * fp, int spin, int color, int iter, double *prop[])
{
  int t, channel, i, j, k;
  if (this_node == 0)
  {
    if ((fscanf(fp, "%d%d%d", &i, &j, &k)) != 3)
    {
      printf("Error in reading meson header\n");
      terminate(1);
    }
    if (i != spin || j != color || k != iter)
    {
      printf("file spin = %d, file color=%d, file iter=%d\n", i, j, k);
      printf("prog spin = %d, prog color=%d, prog iter=%d \n", spin, color, iter);
      terminate(1);
    }
    /* next the elements */

    for (t = 0; t < nt; t++)
    {
      for (channel = 0; channel < nchannels; channel++)
      {
	if ((fscanf(fp, "%lE", &(prop[channel][t]))) == EOF)
	{
	  printf("Read error in r_ascii_m\n");
	  terminate(1);
	}
      }
    }
  }				/* node 0 */
}				/* r_ascii_m  */



/* this subroutine closes the file */
void r_ascii_m_f(FILE * fp, char *filenam)
{
  if (this_node == 0)
  {
    printf("done reading ascii meson propagator in file  %s\n", filenam);
    fclose(fp);
  }
}
