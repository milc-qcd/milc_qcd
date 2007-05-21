/*********************** output_mb2.c *************************/
/* MIMD version 7 */
/* (Identical to heavy/output_mb2.c) */

/* routines for  meson hopping expansion binary  input/output. */

#include "w_static_includes.h"
#include <fcntl.h>  /** raw UNIX calls **/
#include <errno.h>

/* write binary meson propagators */
/*
  format: 
    version_number (int) 
    nx ny nz nt (int) 
    beta kappa_light (Real)
    width^-2 (Real), wallflag(int) 
    wall_cutoff (int), wall_separation (int)
    kappa_c (Real)    
    total_number_of_hopping_iters (int) 
    for spin{ for color{   
       for(N_iter) 
          spin color N_iter 
             prop 
*/


/* this subroutine opens the file for writing */
int w_binary_m_i(char *filenam, int i1)
{
  int fp;
  int32type dims[4];
  int32type i_out;

  /* node 0 does all the writing */
  if (this_node == 0)
  {
    printf(" meson output file is named %s\n", filenam);
    fp = creat(filenam, 0664);
    if (fp < 0)
    {
      printf("Can't open file %s, error %d\n", filenam, errno);
      terminate(1);
    }
    i_out = (int32type)VERSION_NUMBER;
    if ((write(fp, &i_out, sizeof(int32type))) != sizeof(int32type))
    {
      printf("Error in writing meson prop header\n");
      terminate(1);
    }
    dims[XUP] = (int32type)nx;
    dims[YUP] = (int32type)ny;
    dims[ZUP] = (int32type)nz;
    dims[TUP] = (int32type)nt;
    if ((write(fp, dims, 4 * sizeof(int32type))) != 4 * sizeof(int32type))
    {
      printf("Error in writing meson prop header\n");
      terminate(1);
    }
    if ((write(fp, &beta, sizeof(Real))) != sizeof(Real))
    {
      printf("Error in writing meson header beta\n");
      terminate(1);
    }
    if ((write(fp, &kappa, sizeof(Real))) != sizeof(Real))
    {
      printf("Error in writing meson header kappa\n");
      terminate(1);
    }
    i_out = (int32type)i1;
    if ((write(fp, &i_out, sizeof(int32type))) != sizeof(int32type))
    {
      printf("Error in writing meson header i1\n");
      terminate(1);
    }
    if ((write(fp, &width, sizeof(Real))) != sizeof(Real))
    {
      printf("Error in writing meson header width\n");
      terminate(1);
    }
    i_out = (int32type)wqs.type;
    if ((write(fp, &i_out, sizeof(int32type))) != sizeof(int32type))
    {
      printf("Error in writing meson header wallflag\n");
      terminate(1);
    }
    i_out = (int32type)wqs.wall_cutoff;
    if ((write(fp, &i_out, sizeof(int32type))) != sizeof(int32type))
    {
      printf("Error in writing meson header wall_cutoff\n");
      terminate(1);
    }
    i_out = (int32type)wall_separation;
    if ((write(fp, &i_out, sizeof(int32type))) != sizeof(int32type))
    {
      printf("Error in writing meson header wall_separation\n");
      terminate(1);
    }
    if ((write(fp, &kappa_c, sizeof(Real))) != sizeof(Real))
    {
      printf("Error in writing meson header kappa_c\n");
      terminate(1);
    }
    i_out = (int32type)nhop;
    if ((write(fp, &i_out, sizeof(int32type))) != sizeof(int32type))
    {
      printf("Error in writing meson header nhop\n");
      terminate(1);
    }
  }
  else fp = -1;  /* If not node zero */

  return fp;
}


/* this subroutine opens the file for writing but does not write headers */
int a_binary_m_i(char *filenam, int i1)
{
  int  fd;

  /* node 0 does all the writing */
  if (this_node == 0)
  {
    printf(" meson output file is named %s\n", filenam);
    fd = open(filenam, O_WRONLY, 0);
    if (fd < 0)
    {
      printf("Can't open file %s, error %d\n", filenam, errno);
      terminate(1);
    }
  }

  else fd = -1; /* If not node 0 */
  return fd;

}



/* this subroutine writes a meson propagator */
/*
 * all nodes have the same summed meson propagator at this point so node 0
 * does all the work  
 */


void w_binary_m(int fp, int spin, int color, int iter, double *prop[])
{
  int32type ints[3];
  int channel, t, tp, nw1, nw2;
  long head_size, offset, body_size, meson_prop_size, declare_spin_color_etc;
  long ls1, ls2;
  static Real *temp;
  static int firsttime = 1;
  g_sync();

  head_size = 10 * sizeof(int32type) + 4 * sizeof(Real);
  declare_spin_color_etc = 3 * sizeof(int32type);
  meson_prop_size = nchannels * nt * sizeof(Real);

  body_size = declare_spin_color_etc + meson_prop_size;
  body_size = ((1 + iter) + (1 + nhop) * (color + 3 * spin)) * body_size;

  offset = head_size + body_size;


  if (this_node == 0)
  {
    ls1 = lseek(fp, offset, 0);

    /* first the spin and color and iteration_number  */
    ints[0] = spin;
    ints[1] = color;
    ints[2] = iter;
    if ((nw1 = write(fp, ints, 3 * sizeof(int32type))) != 3 * sizeof(int32type))
    {
      printf("Write error in save_binary\n");
      terminate(1);
    }
    ls2 = lseek(fp, 0L, 1);

    /* next the elements */


    /* write Reals to save space */
    if (firsttime)
    {
      temp = (Real *) malloc(nt * sizeof(Real));
      firsttime = 0;
    }
    for (channel = 0; channel < nchannels; channel++)
    {
      /* undo even-slices-first order */
      for (t = 0; t < nt; t++)
      {
	tp = (t % 2) * (nt / 2) + t / 2;
	temp[t] = prop[channel][tp];
      }
      if ((nw2 = write(fp, temp, nt * sizeof(Real))) != nt * sizeof(Real))
      {
	printf("Write error in w_binary_m\n");
	terminate(1);
      }
    }
  }				/* node 0 */
  g_sync();
}				/* w_binary_m  */


/* this subroutine closes the file */
void w_binary_m_f(int fp, char *filenam)
{
    g_sync();
  if (this_node == 0)
  {
    printf("Saved binary meson propagator in file  %s\n", filenam);
    close(fp);
  }
}
