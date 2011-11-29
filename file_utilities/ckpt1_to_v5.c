/**********************  ckpt1_to_v5.c ***************************/
/* MIMD version 7 */
/* Convert gauge configuration from "checkpoint1" format to version 5 format */

/* 01 Nov 1997 C. DeTar */

/* Usage ...

   ckpt1_to_v5 ckpt1_prop v5_prop
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
#include "../include/dirs.h"


#define PARALLEL 1
#define SERIAL 0

#define NODE_DUMP_ORDER 1
#define NATURAL_ORDER 0

#define MAX_BUF_LENGTH 4096

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

#define VERSION_NUMBER 59354

/*----------------------------------------------------------------------*/
/* in binary, a buffered version of parallel    */
/* read; fast version that can only be read with same number   */
/* of nodes as it was written */
void restore_checkpoint1(char *filename,Real *c1,Real *c2)
{
  int fd;
  int version_number,dir,dims[4];
  su3_matrix lbuf[4*MAX_BUF_LENGTH];
  int buf_length,num_buffers,where_in_buf;
  register site *s;
  register int i;
  long head_size,offset;
  head_size = 5*sizeof(int) + 2*sizeof(Real);
  fd = open(filename,O_RDONLY,0);  /*  all nodes open file   */
  
  if(this_node==0){        /* node 0 reads the header */
    if(fd < 0){
      printf("Can't open file %s, error %d\n",filename,errno);
      terminate(1);
    }
    if( (read(fd,&version_number,sizeof(int)))!=sizeof(int) ){
      printf("Error in reading lattice header\n"); terminate(1);
    }
    if(version_number != VERSION_NUMBER){
      printf("Incorrect version number in lattice header\n");terminate(1);
    }
    if( (read(fd,dims,4*sizeof(int)))!=4*sizeof(int) ){
      printf("Error in reading lattice header\n"); terminate(1);
    }
    nx = dims[XUP]; ny=dims[YUP]; nz=dims[ZUP]; nt=dims[TUP];
    printf("Lattice dimensions %d %d %d %d\n",nx,ny,nz,nt);
    if( (read(fd,c1,sizeof(Real)))!=sizeof(Real) ){
      printf("Error in reading lattice header\n"); terminate(1);
    }
    if( (read(fd,c2,sizeof(Real)))!=sizeof(Real) ){
      printf("Error in reading lattice header\n"); terminate(1);
    }
    printf("Header parameters %f %f\n",*c1,*c2);
  }
  
  setup();

  buf_length = 0;
  where_in_buf = 0;
  num_buffers = 0;
  
  FORALLSITES(i,s) {
    
    if(where_in_buf == buf_length){  /* get new buffer */
      /* new buffer length  = remaining sites, but never bigger 
	 than MAX_BUF_LENGTH */
      buf_length = sites_on_node - (num_buffers * MAX_BUF_LENGTH); 
      if(buf_length > MAX_BUF_LENGTH) buf_length = MAX_BUF_LENGTH; 
      /* then do read */
      /* each node reads its sites */
      offset = head_size + (this_node*sites_on_node 
			    +num_buffers*MAX_BUF_LENGTH)*4*sizeof(su3_matrix);  
      /* (all previous buffers are of full size) */
      lseek(fd,offset,0);
      if( (read(fd,lbuf,buf_length*4*sizeof(su3_matrix))) !=
	 buf_length*4*sizeof(su3_matrix)){
	printf("Read error in restore_checkpoint1\n"); terminate(1);
      }
      num_buffers++;
      where_in_buf = 0;  /* reset counter */
    }
    
    
    
    /* put into lattice */
    for(dir=XUP;dir<=TUP;dir++)s->link[dir]=lbuf[4*where_in_buf+dir];
    where_in_buf++;
    
  }
  
  
  close(fd);
  g_sync();
  if(this_node==0){
    printf("Restored binary lattice from file  %s\n",filename);
    fflush(stdout);
  }
}

/*----------------------------------------------------------------------*/

int main(int argc, char *argv[])
{

  gauge_file *gf_v5;
  char *filename_ckpt1,*filename_v5;
  Real c1,c2;
  
  if(argc < 3)
    {
      fprintf(stderr,"Usage %s <ckpt1_prop> <v5_prop>\n",argv[0]);
      exit(1);
    }
  filename_ckpt1 = argv[1];
  filename_v5    = argv[2];

  initialize_machine(&argc,&argv);
  /* Remap standard I/O */
  if(remap_stdio_from_args(argc, argv) == 1)terminate(1);

  this_node = mynode();
  number_of_nodes = numnodes();

  /* Read file in checkpoint format */
  restore_checkpoint1(filename_ckpt1,&c1,&c2);
  
  check_unitarity();

  beta = c1;  /* Only node 0 needs to know */

  /* Write file  in version 5 format */
  gf_v5 = save_parallel(filename_v5);
  
  normal_exit(0);

  return 0;
}
