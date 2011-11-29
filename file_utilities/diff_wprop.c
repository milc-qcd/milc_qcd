/*************************** diff_wprop.c ************************/
/* MIMD version 7 */
/* Compare two Wilson prop files of any format */
/* C. DeTar 3/26/05 */

/* Usage ...

   diff_wprop <wprop_file1> <wprop_file2>

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
#include "../include/io_wprop.h"
#include "../include/generic.h"
#include "../include/generic_wilson.h"
#include "../include/io_lat.h"
#ifdef HAVE_QIO
#include "../include/io_scidac_w.h"
#endif

#include <qio.h>

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
  int i,color,spin;
  site *s;
  int dims[4],ndim;
  Real norm2,maxnorm2,avnorm2;
  wilson_vector wdiff;
  wilson_propagator *wprop1, *wprop2;
  char *wprop_file1, *wprop_file2;
  quark_source wqs1, wqs2;

  if(argc < 3){
    node0_printf("Usage %s <wprop_file1> <wprop_file2>\n", argv[0]);
    return 1;
  }

  wprop_file1 = argv[1];
  wprop_file2 = argv[2];

  initialize_machine(&argc,&argv);

  this_node = mynode();
  number_of_nodes = numnodes();

  /* Sniff out the input file types */
  file_type1 = get_file_type(wprop_file1);
  if(file_type1 < 0){
    node0_printf("Can't determine Wilson prop file type %s\n", wprop_file1);
    return 1;
  }
  
  /* For FNAL types we need to look farther to distinguish Wilson
     prop files from KS prop files  */
  if(file_type1 == FILE_TYPE_FM)
    file_type1 = io_detect_fm(wprop_file1);
  
  /* For QIO(LIME) types, same thing */
  if(file_type1 == FILE_TYPE_LIME){
#ifdef HAVE_QIO
    file_type1 = io_detect_w_usqcd(wprop_file1);
#else
    node0_printf("This looks like a QIO file, but to read it requires QIO compilation\n");
    return 0;
#endif
  }
	
  file_type2 = get_file_type(wprop_file2);
  if(file_type2 < 0){
    node0_printf("Can't determine Wilson prop file type %s\n", wprop_file2);
    return 1;
  }
  
  /* For FNAL types we need to look farther to distinguish Wilson
     prop files from KS prop files  */
  if(file_type2 == FILE_TYPE_FM)
    file_type2 = io_detect_fm(wprop_file2);
  
  /* For QIO(LIME) types, same thing */
  if(file_type2 == FILE_TYPE_LIME){
#ifdef HAVE_QIO
    file_type2 = io_detect_w_usqcd(wprop_file2);
#else
    node0_printf("This looks like a QIO file, but to read it requires QIO compilation\n");
    return 0;
#endif
  }
	
  /* Get the lattice dimensions from the first input file */
  if(read_lat_dim_wprop(wprop_file1, file_type1, &ndim, dims)!=0){
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
  
  /* Allocate space for wprops */
  wprop1 = (wilson_propagator *)malloc(sites_on_node*
				       sizeof(wilson_propagator));
  
  if(wprop1 == NULL){
    node0_printf("No room for propagator\n");
    terminate(1);
  }
  
  wprop2 = (wilson_propagator *)malloc(sites_on_node*
				       sizeof(wilson_propagator));
  
  if(wprop2 == NULL){
    node0_printf("No room for propagator\n");
    terminate(1);
  }
  
  if(this_node == 0)printf("Comparing file %s with file %s\n",
			   wprop_file1,wprop_file2);
  
  /* Read all of both files */
  reload_wprop_to_field(RELOAD_SERIAL, wprop_file1, &wqs1, wprop1, 0);
  reload_wprop_to_field(RELOAD_SERIAL, wprop_file2, &wqs2, wprop2, 0);
  
  /* Compare data */
  
  maxnorm2 = 0;
  avnorm2 = 0;
  FORALLSITES(i,s){
    for(color = 0; color < 3; color++)for(spin = 0; spin < 4; spin++)
      {
	sub_wilson_vector(&wprop1[i].c[color].d[spin],
			  &wprop2[i].c[color].d[spin],
			  &wdiff);
	norm2 = magsq_wvec(&wdiff);
	avnorm2 += norm2;
	if(norm2 > maxnorm2)maxnorm2 = norm2;
    }
  }
  
  /* Sum over lattice and normalize */
  g_floatsum(&avnorm2);
  avnorm2 /= 4*3*volume;
  avnorm2 = sqrt(avnorm2);
  maxnorm2 = sqrt(maxnorm2);
  
  fprintf(stderr,"L2 norm difference is mean %e ; max %e\n",
	  avnorm2,maxnorm2);

  free(wprop1);
  free(wprop2);
  free_lattice();

  normal_exit(0);

  return 0;
}
