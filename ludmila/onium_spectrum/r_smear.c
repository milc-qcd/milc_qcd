#include "../generic/generic_includes.h"
#include "../include/io_lat.h"
#include <sys/types.h>
#include <fcntl.h>
#include <errno.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include "../include/io_lat.h" /* for utilities like get_f ,etc */
#include "../include/io_prop_ks.h"

#define IO_UNI_MAGIC 0x71626434

void get_smearings_bi_serial(char *filename)
{
  FILE *fp;
  int destnode;
  int x,y,z,t,i, byterevflag;
  complex w_smear;
  int32type tmp, magic_number,elements_per_site; 
  int32type  size_of_element, order, dim[4]; 
  int32type  t_stamp;

  /*READING FILE HEADER*/

  if(this_node==0){
    fp = fopen(filename,"rb");
    if(fp==NULL){
      printf("Can't open file %s, error %d\n",filename,errno);
      terminate(1);
    }

    if(fread(&magic_number,sizeof(int32type),1,fp) != 1)
      {
	printf("error in reading magic number from file %s\n", filename);
	terminate(1);
      }
    tmp=magic_number;
    if(magic_number == IO_UNI_MAGIC) byterevflag=0;
    else 
      {
	byterevn((int32type *)&magic_number,1);
      if(magic_number == IO_UNI_MAGIC) 
	{
	  byterevflag=1; 
	  printf("Reading with byte reversal\n");
	  if( sizeof(float) != sizeof(int32type)) {
	    printf("%s: Can't byte reverse\n", filename);
	    printf("requires size of int32type(%d) = size of float(%d)\n",
		   (int)sizeof(int32type),(int)sizeof(float));
	    terminate(1);
	  }
	}
      else
	{
	  /* Restore magic number as originally read */
	  magic_number = tmp;
	  
	  /* End of the road. */
	  printf("%s: Unrecognized magic number in prop file header.\n",
		 filename);
	  printf("Expected %x but read %x\n",
		 IO_UNI_MAGIC,tmp);
	  terminate(1);
	}
    };
    if(fread(&t_stamp,sizeof(t_stamp),1,fp) != 1)
     {	
       printf("error in reading time stamp from file %s\n", filename);
       terminate(1);
     }

    if(fread(&size_of_element,sizeof(int32type),1,fp) != 1)
     {	
       printf("error in reading size of element from file %s\n", filename);
       terminate(1);
     }

    if(fread(&elements_per_site,sizeof(int32type),1,fp) != 1)
      {	
	printf("error in reading elements per site from file %s\n", filename);
	terminate(1);
      }
    if(psread_byteorder(byterevflag,0,fp,dim,sizeof(dim),
		       filename,"dimensions")!=0) terminate(1);

    if( dim[0]!=nx || dim[1]!=ny || 
	dim[2]!=nz || dim[3]!=1 )
      {
	printf(" Incorrect lattice size %d,%d,%d,%d\n",
	       dim[0], dim[1], dim[2], dim[3]);
	terminate(1);
      }
    if( size_of_element != sizeof(float) ||
	elements_per_site != 2 /* complex */)
      {	
	printf(" file %s is not a smearing function",
	       filename);
	terminate(1);
      }
    
    if(psread_byteorder(byterevflag,0,fp,&order,sizeof(int32type),
			filename,"order parameter")!=0) terminate(1);
  } /*if this_node==0*/
  
  for(z=0;z<nz;z++)for(y=0;y<ny;y++)for(x=0;x<nx;x++){
    /* Node 0 reads, and sends site to correct node */
    if(this_node==0)
      if(psread_byteorder(byterevflag,0,fp,&w_smear,sizeof(complex),
			  filename,"reading smearing function")!=0) terminate(1);
  
      for(t=0;t<nt;t++){
	destnode=node_number(x,y,z,t);

	if(this_node==0){
	  if(destnode==0){	/* just copy links */
	    i = node_index(x,y,z,t);
	    lattice[i].w = w_smear;
	  }
	  else {		/* send to correct node */
	    send_field((char *)&(w_smear), sizeof(complex),destnode);
	  }
	}
        /* The node which contains this site reads message */
	else {	/* for all nodes other than node 0 */
	  if(this_node==destnode){
	    get_field((char *)&w_smear, sizeof(complex),0);
	    i = node_index(x,y,z,t);
	    lattice[i].w = w_smear;
	  }
	}
      }
  }
  g_sync();
  
  if(this_node==0){
    printf("Restored smearing function from binary file  %s\n",
	   filename);
    fclose(fp);
    fflush(stdout);
  }
}

