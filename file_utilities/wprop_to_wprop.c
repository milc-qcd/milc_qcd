/*************************** wprop_to_scidac.c ************************/
/* MIMD version 6 */
/* Read a Wilson prop, convert to SciDAC format */
/* C. DeTar 3/25/05 */

/* Usage ...

   wprop_to_scidac milc_file scidac_file

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
#include "../include/io_wprop.h"
#include "../include/generic.h"
#include "../include/generic_wilson.h"
#include "../include/io_lat.h"
#include <qio.h>

#define MAX_RECXML 512

static file_type w_prop_list[N_WPROP_TYPES] =
  { {FILE_TYPE_W_PROP,       W_PROP_VERSION_NUMBER},
    {FILE_TYPE_W_PROP_1996,  W_PROP_VERSION_NUMBER_1996},
    {FILE_TYPE_W_FMPROP,     W_FMPROP_VERSION_NUMBER},
    {FILE_TYPE_W_QIOPROP,    LIME_MAGIC_NO}
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

int get_prompt( int *prompt ){
    char initial_prompt[80];
    int status;

    *prompt = -1;
    printf( "type 0 for no prompts  or 1 for prompts\n");
    status = scanf("%s",initial_prompt);
    if(status != 1){
      printf("\nget_prompt: Can't read stdin\n");
      terminate(1);
    }
    if(strcmp(initial_prompt,"prompt") == 0)  {
       scanf("%d",prompt);
    }
    else if(strcmp(initial_prompt,"0") == 0) *prompt=0;
    else if(strcmp(initial_prompt,"1") == 0) *prompt=1;

    if( *prompt==0 || *prompt==1 )return(0);
    else{
        printf("\nget_prompt: ERROR IN INPUT: initial prompt\n");
        return(1);
    }
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

#define IF_OK if(status==0)

int readin(prompt)
{
  int status = 0;

  if(this_node==0) {

    IF_OK status += ask_starting_wprop( prompt, &(par_buf.startflag),
					par_buf.startfile );
    /* find out what to do with lattice at end */
    IF_OK status += ask_ending_wprop( prompt, &(par_buf.saveflag),
				      par_buf.savefile );
    if(status > 0)par_buf.stopflag = 1; else par_buf.stopflag = 0;
  }

  /* Node 0 broadcasts parameter buffer to all other nodes */
  broadcast_bytes((char *)&par_buf,sizeof(par_buf));
  
  return par_buf.stopflag;
}

/*----------------------------------------------------------------------*/

int main(int argc, char *argv[])
{

  int file_type;
  char recxml[MAX_RECXML];
  int dims[4],ndim;
  int prompt;
  wilson_propagator *wprop;

  initialize_machine(argc,argv);

  this_node = mynode();
  number_of_nodes = numnodes();

  if(get_prompt(&prompt) != 0)
    return 0;

  node0_printf("BEGIN\n");

  /* Loop over input requests */
  while(readin(prompt) == 0)
    {
      /* Sniff out the input file type */
      file_type = io_detect(par_buf.startfile, w_prop_list, N_WPROP_TYPES);
      if(file_type < 0){
	node0_printf("Can't determine Wilson prop file type %s\n", par_buf.startfile);
	return 1;
      }
      
      /* Get the lattice dimensions from the input file */
      read_lat_dim_wprop(par_buf.startfile, file_type, &ndim, dims);
      
      if(this_node == 0)
	{
	  nx = dims[0]; ny = dims[1]; nz = dims[2]; nt = dims[3];
	  printf("Dimensions %d %d %d %d\n",nx,ny,nz,nt);
	}

      volume=nx*ny*nz*nt;
      
      /* Finish setup - broadcast dimensions */
      setup_refresh();
      
      /* Allocate space for the full propagator on each node */
      wprop = (wilson_propagator *)malloc(sites_on_node*
					  sizeof(wilson_propagator));
      if(wprop == NULL){
	node0_printf("No room for propagator\n");
	terminate(1);
      }
      
      if(this_node == 0)printf("Converting file %s to file %s\n",
			       par_buf.startfile, par_buf.savefile);
      
      /* Read the whole propagator */
      reload_wprop_to_field(par_buf.startflag, par_buf.startfile, wprop, 1);
      
      /* Write the whole propagator */
      /* Some arbitrary metadata */
      snprintf(recxml,MAX_RECXML,"Converted from %s",par_buf.startfile);
      save_wprop_from_field(par_buf.saveflag, par_buf.savefile, recxml, 
			    wprop, 1);

      free(wprop);
      free_lattice();
    }

  node0_printf("RUNNING COMPLETED\n");

  return 0;
}
