/*************************** lattice_to_scidac.c ************************/
/* MIMD version 7 */
/* Read gauge configuration, convert any gauge cfg lattice to SciDAC format */
/* MPP or single-processor */
/* Reads the entire lattice before conversion, so requires sufficient memory */
/* C. DeTar 1/30/05 */

/* Usage ...

   lattice_to_scidac milc_file scidac_file

*/

#define CONTROL

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <math.h>
#include "../include/complex.h"
#include "../include/su3.h"
#include <lattice.h>
#include "../include/macros.h"
#include "../include/comdefs.h"
#include <fcntl.h>
#include <errno.h>
#include <time.h>
#include <string.h>
#include "../include/io_lat.h"
#include "../include/io_wprop.h"
#include "../include/generic.h"

#include <qio.h>

#define PARALLEL 1
#define SERIAL 0

#define NODE_DUMP_ORDER 1
#define NATURAL_ORDER 0

#define MAX_BUF_LENGTH 4096

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

  gauge_file *gf;
  gauge_header *gh;
  gauge_file *r_serial_i(char *filename);
  void r_serial_f(gauge_file *gf);
  char *filename_milc,*filename_scidac;
  
  if(argc < 3)
    {
      fprintf(stderr,"Usage %s <MILC file> <SciDAC file>\n",argv[0]);
      exit(1);
    }
  filename_milc   = argv[1];
  filename_scidac = argv[2];

  if(this_node == 0)printf("Converting file %s to SciDAC file %s\n",
			   filename_milc, filename_scidac);

  initialize_machine(&argc,&argv);

  this_node = mynode();
  number_of_nodes = numnodes();

  /* Read header */
  nx = ny = nz = nt = -1;  /* To suppress dimension checking */
  gf = r_serial_i(filename_milc);

  /* We don't do a conversion if this is already a SciDAC file */
  if(gf->header->magic_number == LIME_MAGIC_NO){
    printf("No conversion needed.\n");
    return 1;
  }
  gh = gf->header;

  if(this_node == 0)
    {
      nx = gh->dims[0];
      ny = gh->dims[1];
      nz = gh->dims[2];
      nt = gh->dims[3];
      printf("Dimensions %d %d %d %d\n",nx,ny,nz,nt);
      if(gh->magic_number != GAUGE_VERSION_NUMBER_ARCHIVE){
	printf("Time stamp %s\n",gh->time_stamp);
	if(gh->order == NATURAL_ORDER)printf("File in natural order\n");
	if(gh->order == NODE_DUMP_ORDER)printf("File in node dump order\n");
      }
    }

  /* Close file */
  r_serial_f(gf);

  setup();

  /* Read the entire file this time */
  restore_serial(filename_milc);
  
  check_unitarity();

  /* Write file in SciDAC format */
  save_serial_scidac(filename_scidac);

  normal_exit(0);

  return 0;
}
