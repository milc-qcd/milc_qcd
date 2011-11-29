/*************************** diff_colormatrix.c ************************/
/* MIMD version 7 */
/* Compare two color matrix files (SciDAC format only) */
/* C. DeTar 3/27/05 */

/* Usage ...

   diff_colormatrix <colormatrix_file1> <colormatrix_file2>

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

#include <qio.h>
#include "../include/io_scidac.h"

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

  int i,dir;
  site *s;
  int dims[4],ndim;
  Real norm2,maxnorm2,avnorm2;
  su3_matrix cmdiff;
  su3_matrix *colormatrix1, *colormatrix2;
  char *colormatrix_file1, *colormatrix_file2;

  if(argc < 3){
    node0_printf("Usage %s <colormatrix_file1> <colormatrix_file2>\n", argv[0]);
    return 1;
  }

  colormatrix_file1 = argv[1];
  colormatrix_file2 = argv[2];

  initialize_machine(&argc,&argv);

  this_node = mynode();
  number_of_nodes = numnodes();

  /* Get the lattice dimensions from the first input file */
  if(read_lat_dim_scidac(colormatrix_file1, &ndim, dims)!=0){
    terminate(1);
  }
  
  if(this_node == 0)
    {
      nx = dims[0]; ny = dims[1]; nz = dims[2]; nt = dims[3];
      printf("Dimensions %d %d %d %d\n",nx,ny,nz,nt);
    }
  
  /* Finish setup - broadcast dimensions */
  setup_refresh();
  
  volume=nx*ny*nz*nt;
  
  /* Allocate space for 4 color matrices */
  colormatrix1 = (su3_matrix *)malloc(sites_on_node*4*sizeof(su3_matrix));
  
  if(colormatrix1 == NULL){
    node0_printf("No room for color matrix\n");
    terminate(1);
  }
  
  colormatrix2 = (su3_matrix *)malloc(sites_on_node*4*sizeof(su3_matrix));
  
  if(colormatrix2 == NULL){
    node0_printf("No room for color matrix\n");
    terminate(1);
  }
  
  if(this_node == 0)printf("Comparing file %s with file %s\n",
			   colormatrix_file1,colormatrix_file2);
  
  /* Read all of both files */
  restore_color_matrix_scidac_to_field (colormatrix_file1, colormatrix1, 4);
  restore_color_matrix_scidac_to_field (colormatrix_file2, colormatrix2, 4);
  
  /* Compare data */
  
  maxnorm2 = 0;
  avnorm2 = 0;
  FORALLSITES(i,s){
    for(dir = 0; dir < 4; dir++){
      sub_su3_matrix(&colormatrix1[4*i+dir],&colormatrix2[4*i+dir],&cmdiff);
      norm2 = realtrace_su3(&cmdiff,&cmdiff);
      avnorm2 += norm2;
      if(norm2 > maxnorm2)maxnorm2 = norm2;
    }
  }
  
  /* Sum over lattice and normalize */
  g_floatsum(&avnorm2);
  avnorm2 /= 4*volume;
  avnorm2 = sqrt(avnorm2);
  maxnorm2 = sqrt(maxnorm2);
  
  if(this_node==0)fprintf(stderr,"L2 norm difference is mean %e ; max %e\n",
	  avnorm2,maxnorm2);
  
  free(colormatrix1);
  free(colormatrix2);
  free_lattice();

  normal_exit(0);

  return 0;
}
