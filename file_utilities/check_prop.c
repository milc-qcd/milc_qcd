UNSUPPORTED
/**************************** check_prop.c **********************/
/* MIMD version 7 */
/* Read propagator check checksums (used in version 5) */


/* C. DeTar 10/30/97 */

/* Usage ...

   check_prop propfile
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

#define PARALLEL 1
#define SERIAL 0

#define NODE_DUMP_ORDER 1
#define NATURAL_ORDER 0

#undef MAX_BUF_LENGTH
#define MAX_BUF_LENGTH 4096

/*----------------------------------------------------------------------*/

/* Here only node 0 reads the Wilson propagator from a binary file */

int r_check_w(w_prop_file *wpf, int spin, int color)
{
  /* wpf  = propagator file structure */

  FILE *fp;
  w_prop_header *wph;
  char *filename;
  int byterevflag;
  int spinindex;            /* Counts spin records in file   -
			       wph->spins[spinindex] = spin */

  off_t offset ;            /* File stream pointer */
  off_t w_prop_size;        /* Size of propagator blocks for all nodes */
  off_t w_prop_check_size;  /* Size of propagator checksum record */
  off_t coord_list_size;    /* Size of coordinate list in bytes */
  off_t head_size;          /* Size of header plus coordinate list */
  off_t body_size ;         /* Size of propagator blocks for all nodes 
			      plus checksum record */
  int rcv_rank, rcv_coords;
  int destnode;
  int i,k,x,y,z,t;
  int buf_length,where_in_buf;
  w_prop_check test_wpc;
  u_int32type *val;
  int rank29,rank31;
  wilson_vector *lbuf;
  wilson_vector *dest;
  wilson_vector work;
  int status;
  char myname[] = "r_check_w";

  fp = wpf->fp;
  wph = wpf->header;
  filename = wpf->filename;
  byterevflag = wpf->byterevflag;

  status = 0;
  if(this_node == 0)
    {

      if(wpf->parallel == PARALLEL)
	printf("%s: Attempting serial read from parallel file \n",myname);
      
      w_prop_size = volume*sizeof(wilson_vector) ;
      w_prop_check_size =  sizeof(wpf->check.spin) +
	sizeof(wpf->check.color) +  sizeof(wpf->check.sum29);
      /* 1996 format had an unused 32-bit checksum.
	 Version 5 format has two 32-bit checksums */
      if(wph->magic_number == W_PROP_VERSION_NUMBER)
	w_prop_check_size += sizeof(wpf->check.sum31);

      body_size = w_prop_size + w_prop_check_size;
      
      /* Look up requested spin in header spin table of contents */
      
      for(spinindex=0 ; spinindex < wph->n_spins ; spinindex++)
	if(wph->spins[spinindex]==spin)break;
      
      if(spinindex == wph->n_spins)
	{
	  printf("%s: Requested spin %d not in file %s\n",
		 myname,spin,filename);
	  printf("  Table of contents: ");
	  for(spinindex=0; spinindex<wph->n_spins; spinindex++)
	    printf(" %d",wph->spins[spinindex]);
	  printf("\n");fflush(stdout);
	  terminate(1);
	}
      
      if(wph->order == NATURAL_ORDER)coord_list_size = 0;
      else coord_list_size = sizeof(int32type)*volume;
      head_size = wph->header_bytes + coord_list_size;
      
      offset = head_size + body_size*(spinindex*3 + color);
      
      lbuf = (wilson_vector *)malloc(MAX_BUF_LENGTH*sizeof(wilson_vector));
      if(lbuf == NULL)
	{
	  printf("%s: Node %d can't malloc lbuf\n",myname,this_node);
	  fflush(stdout);
	  terminate(1);
	}
      
      /* Position file pointer for reading check record */
      
      if( fseek(fp,offset,SEEK_SET) < 0 ) 
	{
	  printf("%s: Node %d fseek %ld failed error %d file %s\n",
		 myname,this_node,(long)offset,errno,filename);
	  fflush(stdout);
	  status = 1;
	}
      
    }

  broadcast_bytes((char *)&status,sizeof(int));
  if(status != 0)return status;
  
  status = 0;
  if(this_node == 0)
    {
      /* Read check record */

      status += sread_byteorder(byterevflag,fp,&wpf->check.spin,
		      sizeof(wpf->check.spin),myname,"check.spin");
      status += sread_byteorder(byterevflag,fp,&wpf->check.color,
		      sizeof(wpf->check.color),myname,"check.color");
      status += sread_byteorder(byterevflag,fp,&wpf->check.sum29,
		      sizeof(wpf->check.sum29),myname,"check.sum29");
      /* 1996 format had an unused 32-bit checksum.
	 Version 5 format has two 32-bit checksums */
      if(wph->magic_number == W_PROP_VERSION_NUMBER)
	status += sread_byteorder(byterevflag,fp,&wpf->check.sum31,
			sizeof(wpf->check.sum31),myname,"check.sum31");

      /* Verify spin and color - checksums come later */
      
      if(wpf->check.spin != spin || wpf->check.color != color)
	{
	  printf("%s: Spin %d and color %d do not match check record on file %s\n",
		 myname,spin,color,filename);
	  printf("  Check record said %d %d\n",wpf->check.spin,wpf->check.color);
	  fflush(stdout);
	  status++;
	}

      buf_length = 0;
      where_in_buf = 0;

    }

  broadcast_bytes((char *)&status,sizeof(int));
  if(status != 0)return status;
  
  /* all nodes initialize checksums */
  test_wpc.sum31 = 0;
  test_wpc.sum29 = 0;
  /* counts 32-bit words mod 29 and mod 31 in order of appearance
     on file */
  /* Here all nodes see the same sequence because we read serially */
  rank29 = 0;
  rank31 = 0;
  
  g_sync();

  /* Node 0 reads and deals out the values */
  status = 0;
  for(rcv_rank=0; rcv_rank<volume; rcv_rank++)
    {
      /* If file is in coordinate natural order, receiving coordinate
         is given by rank Otherwise, it is found in the table */
      
      if(wpf->header->order == NATURAL_ORDER)
	rcv_coords = rcv_rank;
      else
	rcv_coords = wpf->rank2rcv[rcv_rank];

      x = rcv_coords % nx;   rcv_coords /= nx;
      y = rcv_coords % ny;   rcv_coords /= ny;
      z = rcv_coords % nz;   rcv_coords /= nz;
      t = rcv_coords % nt;
      
      /* The node that gets the next Wilson vector */
      destnode=node_number(x,y,z,t);

      if(this_node==0){
	/* Node 0 fills its buffer, if necessary */
	if(where_in_buf == buf_length)
	  {  /* get new buffer */
	    /* new buffer length  = remaining sites, but never bigger 
	       than MAX_BUF_LENGTH */
	    buf_length = volume - rcv_rank;
	    if(buf_length > MAX_BUF_LENGTH)buf_length = MAX_BUF_LENGTH;
	    /* then do read */
	    
	    if( fread(lbuf,sizeof(wilson_vector),buf_length,fp) != buf_length)
	      {
		if(status == 0)
		  printf("%s: node %d propagator read error %d file %s\n",
			 myname,this_node,errno,filename); 
		fflush(stdout); 
		status = 1;
	      }
	    where_in_buf = 0;  /* reset counter */
	  }  /*** end of the buffer read ****/
	
	if(destnode==0){	/* just copy Wilson vector */
	  i = node_index(x,y,z,t);
	  dest = &work;
	  memcpy((void *)dest, (void *)&lbuf[where_in_buf], 
		 sizeof(wilson_vector));
	}
	else {		/* send to correct node */
	  send_field((char *)&lbuf[where_in_buf],
		     sizeof(wilson_vector),destnode);
	}
	where_in_buf++;
      }
      
      /* The node which contains this site reads message */
      else {	/* for all nodes other than node 0 */
	if(this_node==destnode){
	  i = node_index(x,y,z,t);
	  dest = &work;
	  get_field((char *)dest, sizeof(wilson_vector),0);
	}
      }
      /* The receiving node does the checksum and then byte reversal,
         if needed */
      
      if(this_node==destnode)
	{
	  if(byterevflag==1)
	    byterevn((int32type *)dest,
		     sizeof(wilson_vector)/sizeof(int32type));
	  
	  /* Accumulate checksums */
	  for(k = 0, val = (u_int32type *)dest; 
	      k < sizeof(wilson_vector)/sizeof(int32type); k++, val++)
	    {
	      test_wpc.sum29 ^= (*val)<<rank29 | (*val)>>(32-rank29);
	      test_wpc.sum31 ^= (*val)<<rank31 | (*val)>>(32-rank31);
	      rank29++; if(rank29 >= 29)rank29 = 0;
	      rank31++; if(rank31 >= 31)rank31 = 0;
	    }
	}
      else
	{
	  rank29 += sizeof(wilson_vector)/sizeof(int32type);
	  rank31 += sizeof(wilson_vector)/sizeof(int32type);
	  rank29 %= 29;
	  rank31 %= 31;
	}
    }

  broadcast_bytes((char *)&status,sizeof(int));
  if(status != 0)return status;

  /* Combine node checksum contributions with global exclusive or */
  g_xor32(&test_wpc.sum29);
  g_xor32(&test_wpc.sum31);
  
  status = 0;
  if(this_node==0)
    {
/*      printf("Read prop serially for spin %d color %d from file %s\n",
	   spin,color,filename); */
      
      /* Verify checksum */
      /* Checksums not implemented until version 5 */
      
      if(wph->magic_number == W_PROP_VERSION_NUMBER)
	{
	  if(wpf->check.sum29 != test_wpc.sum29 ||
	     wpf->check.sum31 != test_wpc.sum31)
	    {
	      printf("%s: Checksum violation spin %d color %d file %s\n",
		     myname,wpf->check.spin,wpf->check.color,wpf->filename);
	      printf("Computed %x %x.  Read %x %x.\n",
		     test_wpc.sum29,test_wpc.sum31,
		     wpf->check.sum29,wpf->check.sum31);
	      status = 1;
	    }
	  else
	    printf("Checksums %x %x OK for spin %d color %d file %s\n",
		   wpf->check.sum29,wpf->check.sum31,
		   wpf->check.spin,wpf->check.color,wpf->filename);
	}
      else
	{
	  printf("Checksums not implemented in this format\n");
	}

      fflush(stdout);
      free(lbuf);
    }

  return status;
  
} /* r_check_w */

/*----------------------------------------------------------------------*/

int main(int argc, char *argv[])
{

  w_prop_file *wpf;
  w_prop_header *wph;
  char *filename;
  int color,spin,spinindex;
  int status;
  
  if(argc < 2)
    {
      fprintf(stderr,"Usage %s <propfilename>\n",argv[0]);
      exit(1);
    }
  filename = argv[1];

  initialize_machine(&argc,&argv);

  /* Remap standard I/O */
  if(remap_stdio_from_args(argc, argv) == 1)terminate(1);

  this_node = mynode();
  number_of_nodes = numnodes();

  /* Read header */
  nx = ny = nz = nt = -1;  /* To suppress dimension checking */
  wpf = r_serial_w_i(filename);
  wph = wpf->header;

  if(this_node == 0)
    {
      printf("Checking file %s\n",filename);
      printf("Time stamp %s\n",wph->time_stamp);
      nx = wph->dims[0];
      ny = wph->dims[1];
      nz = wph->dims[2];
      nt = wph->dims[3];
      printf("Dimensions %d %d %d %d\n",nx,ny,nz,nt);
      if(wph->order == NATURAL_ORDER)printf("File in natural order\n");
      if(wph->order == NODE_DUMP_ORDER)printf("File in node dump order\n");
    }

  broadcast_bytes((char *)&nx,sizeof(int));
  broadcast_bytes((char *)&ny,sizeof(int));
  broadcast_bytes((char *)&nz,sizeof(int));
  broadcast_bytes((char *)&nt,sizeof(int));

  setup_layout();

  volume = nx*ny*nz*nt;

  /* Read and check data */
  status = 0;
  for(color = 0; color < 3; color++)
    for(spinindex = 0; spinindex < wph->n_spins; spinindex++)
      {
	spin = wph->spins[spinindex];
	status += r_check_w(wpf,spin,color);
      }
  
  /* Close file */
  r_serial_w_f(wpf);

  normal_exit(0);

  return 0;
}
