/*************************** ksprop_to_scidac.c ************************/
/* MIMD version 6 */
/* Read a KS prop, convert to SciDAC format */
/* C. DeTar 3/22/05 */

/* Usage ...

   ksprop_to_scidac milc_file scidac_file

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
#include "../include/file_types.h"
#include "../include/io_prop_ks.h"
#include "../include/generic.h"
#include "../include/generic_ks.h"
#include "../include/io_lat.h"
#include <qio.h>

#define MAX_RECXML 512

static file_type ksprop_list[N_KSPROP_TYPES] =
  { {FILE_TYPE_KSPROP,       KSPROP_VERSION_NUMBER},
    {FILE_TYPE_KSFMPROP,     KSFMPROP_VERSION_NUMBER},
    {FILE_TYPE_KSQIOPROP,    LIME_MAGIC_NO}
  };

/*----------------------------------------------------------------------*/
void make_lattice(){
register int i,j;               /* scratch */
int x,y,z,t;            /* coordinates */
    /* allocate space for lattice, fill in parity, coordinates and index.  */
    lattice = (site *)malloc( sites_on_node * sizeof(site) );
    if(lattice==NULL){
        printf("NODE %d: no room for lattice\n",this_node);
        terminate(1);
    }
   /* Allocate address vectors */
    for(i=0;i<8;i++){
        gen_pt[i] = (char **)malloc(sites_on_node*sizeof(char *) );
        if(gen_pt[i]==NULL){
            printf("NODE %d: no room for pointer vector\n",this_node);
            terminate(1);
        }
    }

    for(t=0;t<nt;t++)for(z=0;z<nz;z++)for(y=0;y<ny;y++)for(x=0;x<nx;x++){
        if(node_number(x,y,z,t)==mynode()){
            i=node_index(x,y,z,t);
            lattice[i].x=x;     lattice[i].y=y; lattice[i].z=z; lattice[i].t=t;
            lattice[i].index = x+nx*(y+ny*(z+nz*t));
            if( (x+y+z+t)%2 == 0)lattice[i].parity=EVEN;
            else                 lattice[i].parity=ODD;
#ifdef SITERAND
            initialize_prn( &(lattice[i].site_prn) , iseed, lattice[i].index);
#endif
        }
    }
}
/*----------------------------------------------------------------------*/

void setup() {

  /* Set up lattice */
  broadcast_bytes((char *)&nx,sizeof(int));
  broadcast_bytes((char *)&ny,sizeof(int));
  broadcast_bytes((char *)&nz,sizeof(int));
  broadcast_bytes((char *)&nt,sizeof(int));
  
  setup_layout();
  make_lattice();
}

/*----------------------------------------------------------------------*/

int main(int argc, char *argv[])
{

  char *filename_milc,*filename_scidac;
  int file_type;
  char recxml[MAX_RECXML];
  int i;
  int dims[4],ndim;
  field_offset dest = F_OFFSET(prop[0]);
  
  if(argc < 3)
    {
      fprintf(stderr,"Usage %s <MILC file> <SciDAC file>\n",argv[0]);
      exit(1);
    }
  filename_milc   = argv[1];
  filename_scidac = argv[2];

  if(this_node == 0)printf("Converting file %s to SciDAC file %s\n",
			   filename_milc, filename_scidac);

  initialize_machine(argc,argv);

  this_node = mynode();
  number_of_nodes = numnodes();

  /* Sniff out the input file type */
  file_type = io_detect(filename_milc, ksprop_list, N_KSPROP_TYPES);
  if(file_type < 0){
    node0_printf("Can't read file %s\n", filename_milc);
    return 1;
  }

  /* Get the lattice dimensions from the input file */
  read_lat_dim_ksprop(filename_milc, file_type, &ndim, dims);

  if(this_node == 0)
    {
      nx = dims[0]; ny = dims[1]; nz = dims[2]; nt = dims[3];
      printf("Dimensions %d %d %d %d\n",nx,ny,nz,nt);
    }
  
  /* Finish setup - broadcast dimensions */
  setup();

  /* Read the whole file */
  reload_ksprop(RELOAD_SERIAL, filename_milc, dest, 0);

  /* Write file in SciDAC format */
  /* Some arbitrary metadata */
  snprintf(recxml,MAX_RECXML,"Converted from %s",filename_milc);

  save_ks_vector_scidac(filename_scidac, recxml, QIO_SINGLEFILE, dest, 3);

  return 0;
}
