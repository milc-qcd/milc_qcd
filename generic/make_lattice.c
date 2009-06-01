/******** make_lattice.c *********/

/* Common to most applications */

/* 1. Allocates space for the lattice fields, as specified by the
      application site structure.  Fills in coordinates, parity, index.
   2. Allocates gen_pt pointers for gather results 
   3. Initializes site-based random number generator, if specified
      by macro SITERAND */

#include "generic_includes.h"
#include <defines.h>                 /* For SITERAND */

void make_lattice(){
  register int i;		/* scratch */
  short x,y,z,t;		/* coordinates */
  /* allocate space for lattice, fill in parity, coordinates and index */
  node0_printf("Mallocing %.1f MBytes per node for lattice\n",
	       (double)sites_on_node * sizeof(site)/1e6);
  lattice = (site *)malloc( sites_on_node * sizeof(site) );
  if(lattice==NULL){
    printf("NODE %d: no room for lattice\n",this_node);
    terminate(1);
  }

  /* Allocate address vectors */
  for(i=0;i<N_POINTERS;i++){
    gen_pt[i] = (char **)malloc(sites_on_node*sizeof(char *) );
    if(gen_pt[i]==NULL){
      printf("NODE %d: no room for pointer vector\n",this_node);
      terminate(1);
    }
  }
  
  for(t=0;t<nt;t++)for(z=0;z<nz;z++)for(y=0;y<ny;y++)for(x=0;x<nx;x++){
    if(node_number(x,y,z,t)==mynode()){
      i=node_index(x,y,z,t);
      lattice[i].x=x;	lattice[i].y=y;	lattice[i].z=z;	lattice[i].t=t;
      lattice[i].index = x+nx*(y+ny*(z+nz*t));
      if( (x+y+z+t)%2 == 0)lattice[i].parity=EVEN;
      else	         lattice[i].parity=ODD;
#ifdef SITERAND
      initialize_prn( &(lattice[i].site_prn) , iseed, lattice[i].index);
#endif
    }
  }
}

void free_lattice()
{
  int i;

  for(i=0;i<N_POINTERS;i++)
    free(gen_pt[i]);

  free(lattice);
}
