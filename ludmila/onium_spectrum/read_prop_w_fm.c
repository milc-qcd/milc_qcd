#include "../generic/generic_includes.h"
#include "../include/io_lat.h"
#include <sys/types.h>
#include <fcntl.h>
#include <errno.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include "../include/io_lat.h" /* for utilities like get_f ,etc */


#define IO_UNI_MAGIC 0x71626434

void r_prop_w_fm(char *filename, field_offset dest)
{

  FILE *fp;
  int destnode;
  int x,y,z,t,i, byterevflag, c0,s0,c1,s1;
  wilson_matrix q;
  int32type tmp, magic_number,elements_per_site; 
  int32type  size_of_element, order, dims[4]; 
  int32type   t_stamp;
  site *s;
  wilson_propagator *qp;

  /*READING FILE HEADER*/

  if(this_node==0){
    fp = fopen(filename,"rb");
    if(fp==NULL){
      printf("Can't open propagator file %s, error %d\n",filename,errno);
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
    if(psread_byteorder(byterevflag,0,fp,dims,sizeof(dims),
		       filename,"dimensions")!=0) terminate(1);

    if( dims[0]!=nx || dims[1]!=ny || 
	dims[2]!=nz || dims[3]!=nt )
      {
	printf(" Incorrect lattice size %d,%d,%d,%d\n",
	       dims[0], dims[1], dims[2], dims[3]);
	terminate(1);
      }
    if( size_of_element != sizeof(float) ||
	elements_per_site != 288 /* wilson_propagator */)
      {	
	printf(" file %s is not a wilson propagator in arch_fm format\n",
	       filename);
	terminate(1);
      }
    
    if(psread_byteorder(byterevflag,0,fp,&order,sizeof(int32type),
			filename,"order parameter")!=0) terminate(1);
  } /*if this_node==0*/

  printf("fm prop header\n magic number %x\n timestamp %x\n size of elem %d\n el per site %d\n dim %d, %d, %d, %d\n order %d\nbyterevflag %d\n\n", magic_number, t_stamp, size_of_element, elements_per_site,
	 dims[0],dims[1],dims[2],dims[3], order,byterevflag);
	 

  
  for(t=0;t<nt;t++)for(z=0;z<nz;z++)for(y=0;y<ny;y++)for(x=0;x<nx;x++){
    destnode=node_number(x,y,z,t);
    
    /* Node 0 reads, and sends site to correct node */
    if(this_node==0){
      if(psread_byteorder(byterevflag,0,fp,&q,sizeof(wilson_matrix),
			  filename,"reading the wilson propagator\n")!=0) terminate(1);
  
      
      
      if(destnode==0){	/* just copy links */
	i = node_index(x,y,z,t);
	for(s0=0;s0<4;s0++)for(c0=0;c0<3;c0++)
	  for(s1=0;s1<4;s1++)for(c1=0;c1<3;c1++)
	    { 
	      s = &lattice[i];
	      qp = (wilson_propagator *)F_PT(s,dest);
	      qp->c[c0].d[s0].d[s1].c[c1].real 
	      = q.d[s0].c[c0].d[s1].c[c1].real;
	      qp->c[c0].d[s0].d[s1].c[c1].imag 
	      = q.d[s0].c[c0].d[s1].c[c1].imag;
	    }
      }
      else {		/* send to correct node */
	send_field((char *)&q, sizeof(wilson_matrix),destnode);
      }
  
    }
    /* The node which contains this site reads message */
    else {	/* for all nodes other than node 0 */
      if(this_node==destnode){
	get_field((char *)&q, sizeof(wilson_matrix),0);
	i = node_index(x,y,z,t);
       	for(s0=0;s0<4;s0++)for(c0=0;c0<3;c0++)
	  for(s1=0;s1<4;s1++)for(c1=0;c1<3;c1++)
	    {
	      s = &lattice[i];
	      qp = (wilson_propagator *)F_PT(s,dest);
	      qp->c[c0].d[s0].d[s1].c[c1].real 
	      = q.d[s0].c[c0].d[s1].c[c1].real;
	      qp->c[c0].d[s0].d[s1].c[c1].imag 
	      = q.d[s0].c[c0].d[s1].c[c1].imag;	      
	    }
      }
    }
  }
  
  g_sync();
 
  if(this_node==0){
    printf("Restored wilson propagator from binary file  %s, magic number %x\n",
	   filename, IO_UNI_MAGIC);
    fclose(fp);
    fflush(stdout);
  }
}
