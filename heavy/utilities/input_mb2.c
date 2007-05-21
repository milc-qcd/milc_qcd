/*********************** input_mb.c *************************/
/* MIMD version 6 */

/* routines for  meson hopping expansion binary  input */

#include "w_sum_includes.h"
#include "../../include/int32type.h"
#include <fcntl.h>  /** raw UNIX calls **/
#include <errno.h>

static int sixtyfourbits;
static int byterevflag;

/*----------------------------------------------------------------------*/

/* read binary meson propagators */
/*
 * format: version_number (int) nx ny nz nt (int) beta kappa_light (Real)
 * width^-2 (Real), wallflag(int) wall_cutoff (int), wall_separation (int)
 * kappa_c (Real)    total_number_of_hopping_iters (int) for spin{ for
 * color{   for(N_iter) spin color N_iter prop 
 */

/*
 * opens files, reads header, tells calling routine the number of hopping
 * iterations in file 
 */

int r_binary_m_i(char *filenam, int i1, int *file_hops)
{
  int t;
  int32type dims[4];
  int32type version_number;
  Real x1, x2;
  int fp;
  int j;
  int32type i_in ;

  /* node 0 does all the reading */
  if (this_node == 0)
  {
    printf(" input file is named %s\n", filenam);
    fp = open(filenam, O_RDONLY, 0);
    if (fp < 0)
    {
      printf("Can't open file %s, error %d\n", filenam, errno);
      terminate(1);
    }
    if ((read(fp, &version_number, sizeof(int32type))) != sizeof(int32type))
    {
      printf("Error in reading meson header\n");
      terminate(1);
    }
    /* For cases in which we made a mistake and created
     a header with 64-bit integers */

    if(version_number == 0)
      {
	sixtyfourbits = 1;
	printf("input_mb: First 4 bytes of meson binary file were zero.\n");
	printf("Trying to interpret with 64 bit integer format.\n");
	if ((read(fp, &version_number, sizeof(int32type))) != sizeof(int32type))
	  {
	    printf("Error in reading meson header\n");
	    terminate(1);
	  }
      }

    else sixtyfourbits = 0;

    if (version_number == VERSION_NUMBER) byterevflag = 0;
    else
      {
	byterevn((int32type *)&version_number,1);
	if(version_number == VERSION_NUMBER)
	  {
	    byterevflag = 1;
	    printf("Reading meson binary file with byte reversal\n");
	  }
	else
	  {
	    printf("Incorrect version number in meson header\n");
	    terminate(1);
	  }

      }

    for(j=0;j<4;j++)
      {
	if ((read(fp, &dims[j], sizeof(int32type))) !=  sizeof(int32type))
	  {
	    printf("Error in reading meson header\n");
	    terminate(1);
	  }
	/* If 64 bit integers, then we have to read 4 more bytes get the
	   correct low-order bits */
	if(sixtyfourbits) 
	  {
	    if ((read(fp, &dims[j], sizeof(int32type))) !=  sizeof(int32type))
	      {
		printf("Error in reading meson header\n");
		terminate(1);
	      }
	  }
      }

    /* Byte reverse if needed */
    if(byterevflag==1) byterevn(&dims[0],4);

    nx = dims[XUP];
    ny = dims[YUP];
    nz = dims[ZUP];
    t = dims[TUP];
    if (t != nt)
    {
      printf("Incorrect nt %d\n", t);
      terminate(1);
    }
    if ((read(fp, &x1, sizeof(Real))) != sizeof(Real))
    {
      printf("Error in reading meson header\n");
      terminate(1);
    }
    if(byterevflag==1) byterevn((int32type *)&x1,1);

    if ((read(fp, &x2, sizeof(Real))) != sizeof(Real))
    {
      printf("Error in reading meson header\n");
      terminate(1);
    }
    if(byterevflag==1) byterevn((int32type *)&x2,1);

    if ((read(fp, &i_in, sizeof(int32type))) != sizeof(int32type))
    {
      printf("Error in reading meson header\n");
      terminate(1);
    }
    if(sixtyfourbits==1)
      {
	if ((read(fp, &i_in, sizeof(int32type))) != sizeof(int32type))
	  {
	    printf("Error in reading meson header\n");
	    terminate(1);
	  }
      }
    if(byterevflag==1) byterevn(&i_in,1);

    if (fabs(x1 - beta) > EPS || fabs(x2 - kappa) > EPS)
    {
      printf(
       "Warning: input couplings, %e, %e, not equal to program's %e, %e\n ",
	     (double) x1, (double) x2, (double) beta, (double) kappa);
    }
    if (i_in != i1)
    {
      printf("Warning: file has %d props, prog wants %d\n", i_in, i1);
    }
    if ((read(fp, &x1, sizeof(Real))) != sizeof(Real))
      /* this is wall width, which we don't need */
    {
      printf("Error in reading meson binary header 1 \n");
      terminate(1);
    }
    if(byterevflag==1) byterevn((int32type *)&x1,1);

    if ((read(fp, &i_in, sizeof(int32type))) != sizeof(int32type))
      /* this is wallflag, which we don't need */
    {
      printf("Error in reading meson header\n");
      terminate(1);
    }
    if(sixtyfourbits==1)
      {
	if ((read(fp, &i_in, sizeof(int32type))) != sizeof(int32type))
	  /* this is wallflag, which we don't need */
	  {
	    printf("Error in reading meson header\n");
	    terminate(1);
	  }
      }
    if(byterevflag==1) byterevn(&i_in,1);

    if ((read(fp, &i_in, sizeof(int32type))) != sizeof(int32type))
      /* this is wall_cutoff, which we don't need */
    {
      printf("Error in reading meson header\n");
      terminate(1);
    }
    if(sixtyfourbits==1)
      {
	if ((read(fp, &i_in, sizeof(int32type))) != sizeof(int32type))
	  /* this is wall_cutoff, which we don't need */
	  {
	    printf("Error in reading meson header\n");
	    terminate(1);
	  }
      }
    if(byterevflag==1) byterevn(&i_in,1);

    if ((read(fp, &i_in, sizeof(int32type))) != sizeof(int32type))
      /* this is wall_separation, which we don't need */
    {
      printf("Error in reading meson header\n");
      terminate(1);
    }
    if(sixtyfourbits==1)
      {
	if ((read(fp, &i_in, sizeof(int32type))) != sizeof(int32type))
	  /* this is wall_separation, which we don't need */
	  {
	    printf("Error in reading meson header\n");
	    terminate(1);
	  }
      }
    if(byterevflag==1) byterevn(&i_in,1);

    if ((read(fp, &x1, sizeof(Real))) != sizeof(Real))
    {
      printf("Error in reading meson binary header 2 \n");
      terminate(1);
    }
    if(byterevflag==1) byterevn((int32type *)&x1,1);
    if (fabs(x1 - kappa_c) > EPS)
    {
      printf(
	     "File kappa_c, %e, not equal to program's, %e\n",
	     (double) x1, (double) kappa_c);
      terminate(1);
    }
    if ((read(fp, &i_in, sizeof(int32type))) != sizeof(int32type))
    {
      printf("Error in reading meson header\n");
      terminate(1);
    }
    if(sixtyfourbits==1)
      {
	if ((read(fp, &i_in, sizeof(int32type))) != sizeof(int32type))
	  {
	    printf("Error in reading meson header\n");
	    terminate(1);
	  }
      }
    if(byterevflag==1) byterevn(&i_in,1);

    *file_hops = i_in;
  }
  else fp = -1;  /* If not node zero */

  return fp;
}

/* this subroutine reads a meson propagator */
void r_binary_m(int fp, int spin, int color, int iter, double *prop[])
{
  int ints[3];
  int32type i_in;
  int channel, t;
  int j;
  static Real *temp;
  static int firsttime = 1;
  if (this_node == 0)
  {
    for(j=0;j<3;j++)
      {
	if ((read(fp, &i_in, sizeof(int32type))) != sizeof(int32type))
	  {
	    printf("Error in reading meson binary header 3 \n");
	    terminate(1);
	  }
	if(sixtyfourbits==1)
	  {
	    if ((read(fp, &i_in, sizeof(int32type))) != sizeof(int32type))
	      {
		printf("Error in reading meson binary header 3 \n");
		terminate(1);
	      }
	  }
	if(byterevflag==1) byterevn(&i_in,1);
	ints[j] = i_in;
      }

    if (ints[0] != spin || ints[1] != color || ints[2] != iter)
    {
      printf("Data mismatch in meson binary:\n");
      printf("file spin = %d, file color=%d, file iter=%d\n",
	     ints[0], ints[1], ints[2]);
      printf("prog spin = %d, prog color=%d, prog iter=%d\n",
	     spin, color, iter);
      terminate(1);
    }
    /* next the elements */
    /* file is written as Reals to save space */
    if (firsttime)
    {
      temp = (Real *) malloc(nt * sizeof(Real));
      firsttime = 0;
    }
    for (channel = 0; channel < nchannels; channel++)
    {
      if ((read(fp, temp, nt * sizeof(Real))) != nt * (int)sizeof(Real))
      {
	printf("Read error in r_binary_m\n");
	terminate(1);
      }
      if(byterevflag==1) byterevn((int32type *)&temp[0],nt);
      for (t = 0; t < nt; t++)
	prop[channel][t] = temp[t];
    }
  }				/* node 0 */
}				/* r_binary_m  */


/* this subroutine closes the file */
void r_binary_m_f(int fp, char *filenam)
{
  /* cmn ? {  * */
    printf("done reading binary meson propagator in file  %s\n", filenam);
  close(fp);
  /**	cmn ? }  **/
}
