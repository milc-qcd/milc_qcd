/******************************************************************************
FILE:   topo_info.c
DESCRIPTION:
   Report descriptive information from a binary topological density file.

AUTHOR:     Mark Stephenson
DATE:       01/29/98
******************************************************************************/

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <errno.h>
#include "../../include/int32type.h"
#define MAXFILENAME  256   /* ASCII string length for all file names */

#define TOPO_VERSION_NUMBER 66051  /* hex 00010203 */

void byterevn(int32type w[], int n);

int main(int argc, char *argv[])
{
   FILE *fp;
   char fname[FILENAME_MAX];
   int32type topo_magic_number[1];
   int32type dim[4];
   char startfile[MAXFILENAME];
   int32type total_sweeps[1];
   Real ape_weight[1];
   Real q[4];
   int byterevflag;

   /* get input parameters */
   if (argc == 2)
   {
      strcpy(fname, argv[1]);
   }
   else
   {
      printf("Enter input density filename       ");
      scanf("%s", fname);
   }

   /* open binary density file and header information */
   if ( (fp = fopen(fname, "rb")) == NULL )
   {
      fprintf(stderr, "Unable to open input file %s, ", fname);
      fprintf(stderr, "error %d.\n", errno);
      fprintf(stderr, "Exiting program.\n");
      exit(1);
   }
   fread(topo_magic_number, sizeof(int32type), 1, fp);
   if(*topo_magic_number == TOPO_VERSION_NUMBER)
   {
      byterevflag = 0;
   }
   else
   {
      byterevn(topo_magic_number, 1);
      if(*topo_magic_number == TOPO_VERSION_NUMBER)
      {
         byterevflag = 1;
      }
      else
      {
         printf("Wrong version number. Terminating program.\n");
         exit(1);
      }
   }
   fread(dim, sizeof(int32type), 4, fp);
   if (byterevflag == 1)
   {
      byterevn(dim, 4);
   }
   fread(startfile, sizeof(char), 80, fp);
   fread(total_sweeps, sizeof(int32type), 1, fp);
   if (byterevflag == 1)
   {
      byterevn(total_sweeps, 1);
   }
   fread(ape_weight, sizeof(Real), 1, fp);
   if (byterevflag == 1)
   {
      byterevn((int32type *)ape_weight, 1);
   }
   fread(q, sizeof(Real), 4, fp);
   fclose(fp);
   if (byterevflag == 1)
   {
      byterevn((int32type *)q, 4);
   }

   /* report topological density description */
   printf("magic_number  %d\n", *topo_magic_number);
   printf("nx ny nz nt   %d %d %d %d\n", dim[0], dim[1], dim[2], dim[3]);
   printf("startfile     %s\n", startfile);
   printf("total_sweeps  %d\n", *total_sweeps);
   printf("ape_weight    %.3f\n", *ape_weight);
   printf("First few density values:\n");
   printf("%f\n", q[0]);
   printf("%f\n", q[1]);
   printf("%f\n", q[2]);
   printf("%f\n", q[3]);

   /* terminate program normally */
   return 0;
}

/*----------------------------------------------------------------------*/

/* For doing byte reversal on 32-bit words */
/* From io_lat4.c */

void byterevn(int32type w[], int n)
{
  register int32type old,new;
  int j;
  
  for(j=0; j<n; j++)
    {
      old = w[j];
      new = old >> 24 & 0x000000ff;
      new |= old >> 8 & 0x0000ff00;
      new |= old << 8 & 0x00ff0000;
      new |= old << 24 & 0xff000000;
      w[j] = new;
    }
} /* byterevn */
