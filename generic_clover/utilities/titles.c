/*
 *  Print out the titles for the creation of smearing functions
 *
 *
 *
 */


#include<stdio.h>

void titles (int nx, int ny, int nz, int seed, int nosmear)
{

  printf("------------------------------------------------------------\n");
  printf("Set up the sources for the static variational code\n");
  printf("------------------------------------------------------------\n");

  printf("\nDimension of the lattice = %d x %d x %d \n",nx,ny,nz);
  printf("The number of smearing functions = %d\n",nosmear);

  printf("Seed for the random number generator = %d\n",seed);


}
