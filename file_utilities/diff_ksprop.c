/*************************** diff_ksprop.c ************************/
/* MIMD version 7 */
/* Compare two KS prop files of any format */
/* C. DeTar 3/26/05 */

/* Usage ...

   diff_ksprop <ksprop_file1> <ksprop_file2>

*/

#define CONTROL

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include "../include/complex.h"
#include "../include/su3.h"
#include <lattice.h>
#include "../include/macros.h"
#include "../include/comdefs.h"
#include <time.h>
#include <string.h>
#include <math.h>
#include "../include/file_types.h"
#include "../include/io_ksprop.h"
#include "../include/generic.h"
#include "../include/generic_ks.h"
#include "../include/io_lat.h"
#ifndef HAVE_QIO
REQUIRES QIO
#else
#include <qio.h>
#endif

/*----------------------------------------------------------------------*/
void make_lattice(){
register int i;               /* scratch */
int x,y,z,t;            /* coordinates */
    /* allocate space for lattice, fill in parity, coordinates and index.  */
    lattice = (site *)malloc( sites_on_node * sizeof(site) );
    if(lattice==NULL){
        printf("NODE %d: no room for lattice\n",this_node);
        terminate(1);
    }

    for(t=0;t<nt;t++)for(z=0;z<nz;z++)for(y=0;y<ny;y++)for(x=0;x<nx;x++){
        if(node_number(x,y,z,t)==mynode()){
            i=node_index(x,y,z,t);
            lattice[i].x=x;     lattice[i].y=y; lattice[i].z=z; lattice[i].t=t;
            lattice[i].index = x+nx*(y+ny*(z+nz*t));
            if( (x+y+z+t)%2 == 0)lattice[i].parity=EVEN;
            else                 lattice[i].parity=ODD;
        }
    }
}

void free_lattice()
{
  free(lattice);
}

/*----------------------------------------------------------------------*/

void setup_refresh() {

  /* Set up lattice */
  broadcast_bytes((char *)&nx,sizeof(int));
  broadcast_bytes((char *)&ny,sizeof(int));
  broadcast_bytes((char *)&nz,sizeof(int));
  broadcast_bytes((char *)&nt,sizeof(int));
  
  setup_layout();
  make_lattice();
}

/*----------------------------------------------------------------------*/

typedef struct {
  int stopflag;
  int startflag;  /* what to do for beginning propagator */
  int saveflag;   /* what to do with propagator at end */
  char startfile[MAXFILENAME],savefile[MAXFILENAME];
} params;

params par_buf;

/*----------------------------------------------------------------------*/

int main(int argc, char *argv[])
{

  int file_type1, file_type2;
  int i,color;
  site *s;
  int dims[4],ndim;
  Real norm2,maxnorm2,avnorm2;
  su3_vector ksdiff;
  su3_vector *ksprop1, *ksprop2;
  char *ksprop_file1, *ksprop_file2;
  quark_source ksqs1, ksqs2;

  if(argc < 3){
    node0_printf("Usage %s <ksprop_file1> <ksprop_file2>\n", argv[0]);
    return 1;
  }

  ksprop_file1 = argv[1];
  ksprop_file2 = argv[2];

  initialize_machine(&argc,&argv);

  this_node = mynode();
  number_of_nodes = numnodes();

  init_qs(&ksqs1);
  init_qs(&ksqs2);

  /* Sniff out the input file types */
  file_type1 = get_file_type(ksprop_file1);
  if(file_type1 < 0){
    node0_printf("Can't determine KS prop file type %s\n", ksprop_file1);
    return 1;
  }
  
  file_type2 = get_file_type(ksprop_file2);
  if(file_type2 < 0){
    node0_printf("Can't determine KS prop file type %s\n", ksprop_file2);
    return 1;
  }
  
  /* Get the lattice dimensions from the first input file */
  if(read_lat_dim_ksprop(ksprop_file1, file_type1, &ndim, dims)!=0){
    terminate(1);
  }
  
  if(this_node == 0)
    {
      nx = dims[0]; ny = dims[1]; nz = dims[2]; nt = dims[3];
      printf("Dimensions %d %d %d %d\n",nx,ny,nz,nt);fflush(stdout);
    }
  
  /* Finish setup - broadcast dimensions */
  setup_refresh();
  
  volume=nx*ny*nz*nt;
  
  /* Allocate space for ksprops */
  ksprop1 = (su3_vector *)malloc(sites_on_node*3*sizeof(su3_vector));
  
  if(ksprop1 == NULL){
    node0_printf("No room for propagator\n");
    terminate(1);
  }
  
  ksprop2 = (su3_vector *)malloc(sites_on_node*3*sizeof(su3_vector));
  
  if(ksprop2 == NULL){
    node0_printf("No room for propagator\n");
    terminate(1);
  }
  
  if(this_node == 0)printf("Comparing file %s with file %s\n",
			   ksprop_file1,ksprop_file2);
  
  /* Read all of both files */
  reload_ksprop_to_field3(RELOAD_SERIAL, ksprop_file1, &ksqs1, ksprop1, 0);
  reload_ksprop_to_field3(RELOAD_SERIAL, ksprop_file2, &ksqs2, ksprop2, 0);
  
  /* Compare data */
  
  maxnorm2 = 0;
  avnorm2 = 0;
  FORALLSITES(i,s){
    for(color = 0; color < 3; color++){
      sub_su3_vector(&ksprop1[3*i+color],&ksprop2[3*i+color],&ksdiff);
      norm2 = magsq_su3vec(&ksdiff);
      avnorm2 += norm2;
      if(norm2 > maxnorm2)maxnorm2 = norm2;
    }
  }
  
  /* Sum over lattice and normalize */
  g_floatsum(&avnorm2);
  avnorm2 /= 3*volume;
  avnorm2 = sqrt(avnorm2);
  maxnorm2 = sqrt(maxnorm2);
  
  fprintf(stderr,"L2 norm difference is mean %e ; max %e\n",
	  avnorm2,maxnorm2);
  
  free(ksprop1);
  free(ksprop2);
  free_lattice();

  normal_exit(0);

  return 0;
}
