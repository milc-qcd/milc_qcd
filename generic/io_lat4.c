/*********************** io_lat4.c *************************/
/* MIMD version 7 */

/* routines for gauge configuration input/output. */
/* This works for most machines.  Wrappers for parallel I/O
   are in io_ansi.c, io_piofs.c, or io_paragon2.c */

/* Modifications */
/* 10/04/01 Removed save_old_binary (but can still read old binary) C.D. */
/* 7/11/01 large file (64 bit addressing) support */
/* 4/16/00 additions to READ ARChive format J.H. */
/*         adapted for version 7 12/21/00 UMH */
/* 4/17/98 r_parallel_w: g_syncs to prevent shmem message pileups C.D. */
/* 9/19/97 version 5 format with checksums C.D. */
/* 9/04/97 parallel files to be written in typewriter order C.D. */
/* 8/30/96 fixed macros for C syntax UMH */
/* 8/27/96 io_lat3.c converted to parallel reads and writes C.D. */
/*         Synchronization done through message passing instead of g_sync */
/*         Attempt at implementing ANSI standard, by UMH */

#include "generic_includes.h"
#include "../include/io_lat.h"
#include <sys/types.h>
#include <fcntl.h>
#include <errno.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#ifdef HAVE_QIO
#include <qio.h>
#endif

#define EPS 1e-6

#define PARALLEL 1   /* Must evaluate to true */
#define SERIAL 0     /* Must evaluate to false */

#define NODE_DUMP_ORDER 1
#define NATURAL_ORDER 0

#undef MAX_BUF_LENGTH
#define MAX_BUF_LENGTH 4096

/* Checksums
    
   The dataset from which each checksum is computed is the full gauge
   configuration for lattice files and for propagator files, the
   propagator for a single source spin-color combination.  Data in these
   files appear as a series of 32-bit floating point numbers.  We treat
   the 32-bit values as unsigned integers v(i) where i = 0,...,N-1 ranges
   over the values in the order of appearance on the file and N is the
   total count of values in the data set.  The checksum is obtained by a
   combination of bit-rotations and exclusive or operations.  It is
   designed to be commutative and associative, unlike the BSD sum
   operation, so that the data set can be read in parallel with checksum
   contributions computed for each portion read, and then combined
   afterwards.  The sum29 checksum does a left bit rotation through i mod
   29 bits and forms an exclusive or with the accumulated checksum.  The
   sum31 checksum does the same thing, but with i mod 31 bits.
   
   In writing the file the bit rotation is done on the number as
   represented on the architecture and consequently as written on the
   file.  In reading and checking file integrity on an architecure with a
   relatively byte-reversed representation, byte reversal of the data
   must be done before doing the bit rotation and the resulting checksum
   must be compared with the checksum recorded on the file after
   byte-reversal. 
*/

#define SUCCESS  0
#define FAILURE -1
#define MAX_LINE_LENGTH 1024
#define MAX_TOKENS 512

/* For NERSC archive format */
typedef float INPUT_TYPE;
typedef float OUTPUT_TYPE;

/* Version: 1.0 */
#define OLDHEADERSIZE 0
#define TOL 0.0000001  
/* tolerance for floating point checks */
/* For checksums we want a 32 bit unsigned int, for which      */
/* we have u_int32type defined in ../include/int32type.h which is    */
/* included in ../include/io_lat.h .                           */
/*=============================================================*/


/*----------------------------------------------------------------------*/

/* Open a binary file for serial writing by node 0 */

static gauge_file *w_serial_i(char *filename)
{
  /* Only node 0 opens the file filename */
  /* Returns a file structure describing the opened file */

  char myname[] = "w_serial_i";
  FILE *fp;
  gauge_file *gf;
  gauge_header *gh;

  /* Set up gauge file and gauge header structures and load header values */
  gf = setup_output_gauge_file();
  gh = gf->header;

  /* Set number of nodes to zero to indicate coordinate natural ordering */

  gh->order = NATURAL_ORDER;

  /* Only node 0 opens the requested file */

  if(this_node == 0)
    {
      fp = fopen(filename, "wb");
      if(fp == NULL)
	{
	  printf("%s: Node %d can't open file %s, error %d\n",
		 myname,this_node,filename,errno);fflush(stdout);
	  terminate(1);
	}

      /* Node 0 writes the header */
      
      swrite_gauge_hdr(fp,gh);

    }
  
  /* Assign values to file structure */

  if(this_node==0)gf->fp = fp; 
  else gf->fp = NULL;                /* Only node 0 knows about this file */

  gf->filename = filename;
  gf->byterevflag    = 0;            /* Not used for writing */
  gf->rank2rcv       = NULL;         /* Not used for writing */
  gf->parallel       = 0;

  return gf;

} /* w_serial_i */


/*---------------------------------------------------------------------------*/

/* Here only node 0 writes gauge configuration to a binary file */

#if 0
static void w_serial_old(gauge_file *gf)
{
  /* gf  = file descriptor as opened by w_serial_i */

  FILE *fp = NULL;
  gauge_header *gh = NULL;
  u_int32type *val;
  int rank29,rank31;
  fsu3_matrix *lbuf = NULL;
  fsu3_matrix tbuf[4];
  int buf_length;
  register int i,j,k;
  off_t offset;             /* File stream pointer */
  off_t coord_list_size;    /* Size of coordinate list in bytes */
  off_t head_size;          /* Size of header plus coordinate list */
  off_t checksum_offset = 0; /* Location of checksum */
  off_t gauge_check_size;   /* Size of checksum record */

  int currentnode,newnode;
  int x,y,z,t;

  if(this_node==0)
    {
      if(gf->parallel)
	printf("w_serial: Attempting serial write to parallel file \n");

      lbuf = (fsu3_matrix *)malloc(MAX_BUF_LENGTH*4*sizeof(fsu3_matrix));
      if(lbuf == NULL)
	{
	  printf("w_serial: Node 0 can't malloc lbuf\n"); 
	  fflush(stdout);terminate(1);
        }

      fp = gf->fp;
      gh = gf->header;
      
      /* No coordinate list was written because fields are to be written
	 in standard coordinate list order */
      
      coord_list_size = 0;
      head_size = gh->header_bytes + coord_list_size;

      checksum_offset = head_size;

      gauge_check_size = sizeof(gf->check.sum29) + sizeof(gf->check.sum31);
      
      offset = head_size + gauge_check_size;

      if( g_seek(fp,offset,SEEK_SET) < 0 ) 
	{
	  printf("w_serial: Node %d g_seek %lld failed error %d file %s\n",
		 this_node,(long long)offset,errno,gf->filename);
	  fflush(stdout);terminate(1);
	}
    }
      
  /* Buffered algorithm for writing fields in serial order */
  
  /* initialize checksums */
  gf->check.sum31 = 0;
  gf->check.sum29 = 0;
  /* counts 32-bit words mod 29 and mod 31 in order of appearance on file */
  /* Here only node 0 uses these values */
  rank29 = 4*sizeof(fsu3_matrix)/sizeof(int32type)*sites_on_node*this_node % 29;
  rank31 = 4*sizeof(fsu3_matrix)/sizeof(int32type)*sites_on_node*this_node % 31;

  g_sync();
  currentnode=0;

  buf_length = 0;

  for(j=0,t=0;t<nt;t++)for(z=0;z<nz;z++)for(y=0;y<ny;y++)for(x=0;x<nx;x++,j++)
    {
      newnode=node_number(x,y,z,t);
      if(newnode != currentnode){	/* switch to another node */
	/* Node 0 sends a few bytes to newnode to say it's OK to
	   send */
	if( this_node==0 && newnode!=0 )send_field((char *)tbuf,32,newnode);
	if( this_node==newnode && newnode!=0 )get_field((char *)tbuf,32,0);
	currentnode=newnode;
      }
      
      /* Node 0 receives the data */
      if(this_node==0)
	{
	  /* Data on node 0 is just copied to tbuf */
	  if(currentnode==0)
	    {
	      i=node_index(x,y,z,t);
	      d2f_4mat(&lattice[i].link[0],tbuf);
	    }
	  else
	    {
	      /* Data on any other node is received in tbuf */
	      get_field((char *)tbuf,4*sizeof(fsu3_matrix),currentnode);
	    }

	  /* Pack tbufs in lbuf */
	  memcpy((void *)&lbuf[4*buf_length], 
		 (void *)tbuf, 4*sizeof(fsu3_matrix));


	  /* Accumulate checksums - contribution from next site */
	  for(k = 0, val = (u_int32type *)&lbuf[4*buf_length]; 
	      k < 4*(int)sizeof(fsu3_matrix)/(int)sizeof(int32type); 
	      k++, val++)
	    {
	      gf->check.sum29 ^= (*val)<<rank29 | (*val)>>(32-rank29);
	      gf->check.sum31 ^= (*val)<<rank31 | (*val)>>(32-rank31);
	      rank29++; if(rank29 >= 29)rank29 = 0;
	      rank31++; if(rank31 >= 31)rank31 = 0;
	    }

	  buf_length++;
	  

	  if( (buf_length == MAX_BUF_LENGTH) || (j == volume-1))
	    {
	      /* write out buffer */
	      
	      if( (int)g_write(lbuf,4*sizeof(fsu3_matrix),buf_length,fp) != buf_length)
		{
		  printf("w_serial: Node %d gauge configuration write error %d file %s\n",
			 this_node,errno,gf->filename); 
		  fflush(stdout);
		  terminate(1);   
		}
	      buf_length = 0;		/* start again after write */
	    }
	} /* if this_node == 0 */
      else  /* for nodes other than 0 */
	{	
	  if(this_node==currentnode){
	    i=node_index(x,y,z,t);
	    /* Convert 4 matrices from generic to single precision in
	       tbuf and send */
	    d2f_4mat(&lattice[i].link[0],tbuf);
	    send_field((char *)tbuf,4*sizeof(fsu3_matrix),0);
	  }
	}
      
    } /*close x,y,z,t loops */
  
  g_sync();
  
  if(this_node==0)
    {
      free(lbuf);
      printf("Saved gauge configuration serially to binary file %s\n",
	     gf->filename);
      printf("Time stamp %s\n",gh->time_stamp);
      
      /* Write checksum */
      /* Position file pointer */
      if( g_seek(fp,checksum_offset,SEEK_SET) < 0 ) 
	{
	  printf("w_serial: Node %d g_seek %lld failed error %d file %s\n",
		 this_node,(long long)checksum_offset,errno,gf->filename);
	  fflush(stdout);terminate(1);
	}
      write_checksum(SERIAL,gf);
    }
  
} /* w_serial_old */
#endif

/* Flush lbuf to output */
/* buf_length is reset */
static void flush_lbuf_to_file(gauge_file *gf, fsu3_matrix *lbuf, 
			       int *buf_length)
{
  FILE *fp = gf->fp;

  if(*buf_length <= 0)return;
  if( (int)g_write(lbuf,4*sizeof(fsu3_matrix),*buf_length,fp) != 
      *buf_length)
    {
      printf("w_serial: Node %d gauge configuration write error %d file %s\n",
	     this_node,errno,gf->filename); 
      fflush(stdout);
      terminate(1);   
    }
  *buf_length = 0;
}

/* Accumulate checksums */
static void accum_cksums(gauge_file *gf, int *rank29, int *rank31,
			 u_int32type *buf, int n){
  int k;
  u_int32type *val;

  for(k = 0, val = buf; k < n; k++, val++)
    {
      gf->check.sum29 ^= (*val)<<(*rank29) | (*val)>>(32-(*rank29));
      gf->check.sum31 ^= (*val)<<(*rank31) | (*val)>>(32-(*rank31));
      (*rank29)++; if(*rank29 >= 29)*rank29 = 0;
      (*rank31)++; if(*rank31 >= 31)*rank31 = 0;
    }
}

/* Flush tbuf to lbuf and accumulate checksums */
/* tbuf_length is not reset here */
static void flush_tbuf_to_lbuf(gauge_file *gf, int *rank29, int *rank31,
			       fsu3_matrix *lbuf, int *buf_length, 
			       fsu3_matrix *tbuf, int tbuf_length){

  int nword;
  u_int32type *buf;

  if(tbuf_length > 0){
    memcpy((void *)&lbuf[4*(*buf_length)], 
	   (void *)tbuf, 4*tbuf_length*sizeof(fsu3_matrix));

    nword= 4*(int)sizeof(fsu3_matrix)/(int)sizeof(int32type)*tbuf_length;
    buf = (u_int32type *)&lbuf[4*(*buf_length)];
    accum_cksums(gf, rank29, rank31, buf, nword);

    *buf_length += tbuf_length;
  }
}

static void send_buf_to_node0(fsu3_matrix *tbuf, int tbuf_length, 
			      int currentnode){
  if(this_node == currentnode){
    send_field((char *)tbuf,4*tbuf_length*sizeof(fsu3_matrix),0);
  }
  else if(this_node == 0){
    get_field((char *)tbuf,4*tbuf_length*sizeof(fsu3_matrix),
	      currentnode);
  }
}

static void w_serial(gauge_file *gf)
{
  /* gf  = file descriptor as opened by w_serial_i */

  FILE *fp = NULL;
  gauge_header *gh = NULL;
  int rank29,rank31;
  fsu3_matrix *lbuf = NULL;
  fsu3_matrix *tbuf = NULL;
  int buf_length, tbuf_length;
  register int i,j;
  off_t offset;             /* File stream pointer */
  off_t coord_list_size;    /* Size of coordinate list in bytes */
  off_t head_size;          /* Size of header plus coordinate list */
  off_t checksum_offset = 0; /* Location of checksum */
  off_t gauge_check_size;   /* Size of checksum record */

  int currentnode,newnode;
  int x,y,z,t;
  char myname[] = "w_serial";

  /* Allocate message buffer space for one x dimension of the local
     hypercube */
  /* The largest possible space we need is nx */

  tbuf = (fsu3_matrix *)malloc(nx*4*sizeof(fsu3_matrix));
  if(tbuf == NULL){
    printf("%s(%d): No room for tbuf\n",myname,this_node);
    terminate(1);
  }

  if(this_node==0)
    {
      if(gf->parallel)
	printf("w_serial: Attempting serial write to parallel file \n");

      lbuf = (fsu3_matrix *)malloc(MAX_BUF_LENGTH*4*sizeof(fsu3_matrix));
      if(lbuf == NULL)
	{
	  printf("w_serial: Node 0 can't malloc lbuf\n"); 
	  fflush(stdout);terminate(1);
        }

      fp = gf->fp;
      gh = gf->header;
      
      /* No coordinate list was written because fields are to be written
	 in standard coordinate list order */
      
      coord_list_size = 0;
      head_size = gh->header_bytes + coord_list_size;

      checksum_offset = head_size;

      gauge_check_size = sizeof(gf->check.sum29) + sizeof(gf->check.sum31);
      
      offset = head_size + gauge_check_size;

      if( g_seek(fp,offset,SEEK_SET) < 0 ) 
	{
	  printf("w_serial: Node %d g_seek %lld failed error %d file %s\n",
		 this_node,(long long)offset,errno,gf->filename);
	  fflush(stdout);terminate(1);
	}
    }
      
  /* Buffered algorithm for writing fields in serial (lexicographic) order */
  
  /* initialize checksums */
  gf->check.sum31 = 0;
  gf->check.sum29 = 0;
  /* counts 32-bit words mod 29 and mod 31 in order of appearance on file */
  /* Here only node 0 uses these values -- both start at 0 */ 
  rank29 = 4*sizeof(fsu3_matrix)/sizeof(int32type)*sites_on_node*this_node % 29;
  rank31 = 4*sizeof(fsu3_matrix)/sizeof(int32type)*sites_on_node*this_node % 31;

  g_sync();
  currentnode=0;  /* The node delivering data */

  buf_length = 0;
  tbuf_length = 0;
  for(j=0,t=0;t<nt;t++)for(z=0;z<nz;z++)for(y=0;y<ny;y++)for(x=0;x<nx;x++,j++)
    {
      newnode=node_number(x,y,z,t);  /* The node providing the next site */
      if(newnode != currentnode || x == 0){
	/* We are switching to a new node or have exhausted a line of nx */
	/* Sweep any data in the retiring node's tbuf to the node0 lbuf*/
	if(tbuf_length > 0){
	  if(currentnode != 0)
	    send_buf_to_node0(tbuf, tbuf_length, currentnode);

	  if(this_node == 0){
	    /* Node 0 flushes tbuf and accumulates checksum */
	    flush_tbuf_to_lbuf(gf, &rank29, &rank31, lbuf, &buf_length, 
			       tbuf, tbuf_length);
	    /* Node 0 writes lbuf if full */
	    if(buf_length > MAX_BUF_LENGTH - nx)
	      flush_lbuf_to_file(gf, lbuf, &buf_length);
	  }
	  tbuf_length = 0;
	}
	/* Node 0 sends a few bytes to newnode as a clear to send signal */
	if(newnode != currentnode){
	  if( this_node==0 && newnode!=0 )send_field((char *)tbuf,32,newnode);
	  if( this_node==newnode && newnode!=0 )get_field((char *)tbuf,32,0);
	  currentnode=newnode;
	}
      } /* currentnode != newnode */

      /* The node with the data just appends to its tbuf */
      if(this_node == currentnode)
	{
	  i=node_index(x,y,z,t);
	  d2f_4mat(&lattice[i].link[0],&tbuf[4*tbuf_length]);
	}

      if(this_node == currentnode || this_node == 0)tbuf_length++;

    } /*close x,y,z,t loops */

  /* Purge any remaining data */

  if(tbuf_length > 0){
    if(currentnode != 0)
      send_buf_to_node0(tbuf, tbuf_length, currentnode);
  }

  if(this_node == 0){
    flush_tbuf_to_lbuf(gf, &rank29, &rank31, lbuf, &buf_length, 
		       tbuf, tbuf_length);
    flush_lbuf_to_file(gf, lbuf, &buf_length);
  }

  g_sync();
  free(tbuf);
  
  if(this_node==0)
    {
      free(lbuf);
      printf("Saved gauge configuration serially to binary file %s\n",
	     gf->filename);
      printf("Time stamp %s\n",gh->time_stamp);
      
      /* Write checksum */
      /* Position file pointer */
      if( g_seek(fp,checksum_offset,SEEK_SET) < 0 ) 
	{
	  printf("w_serial: Node %d g_seek %lld failed error %d file %s\n",
		 this_node,(long long)checksum_offset,errno,gf->filename);
	  fflush(stdout);terminate(1);
	}
      write_checksum(SERIAL,gf);
    }
  
} /* w_serial */

/*----------------------------------------------------------------------*/

/* Here only node 0 reads the gauge configuration from a binary file */

static void r_serial(gauge_file *gf)
{
  /* gf  = gauge configuration file structure */

  FILE *fp;
  gauge_header *gh;
  char *filename;
  int byterevflag;

  off_t offset = 0 ;        /* File stream pointer */
  off_t gauge_check_size;   /* Size of gauge configuration checksum record */
  off_t coord_list_size;    /* Size of coordinate list in bytes */
  off_t head_size;          /* Size of header plus coordinate list */
  off_t checksum_offset = 0; /* Where we put the checksum */
  int rcv_rank, rcv_coords;
  int destnode;
  int k;
  int x,y,z,t;
  int buf_length = 0, where_in_buf = 0;
  gauge_check test_gc;
  u_int32type *val;
  int rank29,rank31;
  fsu3_matrix *lbuf = NULL;
  fsu3_matrix tmpsu3[4];
  char myname[] = "r_serial";
  int idest = 0;

  fp = gf->fp;
  gh = gf->header;
  filename = gf->filename;
  byterevflag = gf->byterevflag;

  if(this_node == 0)
    {
      /* Compute offset for reading gauge configuration */

      /* (1996 gauge configuration files had a 32-bit unused checksum 
	 record before the gauge link data) */
      if(gh->magic_number == GAUGE_VERSION_NUMBER)
	gauge_check_size = sizeof(gf->check.sum29) + 
	  sizeof(gf->check.sum31);
      else
	gauge_check_size = 0;
      
      if(gf->header->order == NATURAL_ORDER)coord_list_size = 0;
      else coord_list_size = sizeof(int32type)*volume;
      checksum_offset = gf->header->header_bytes + coord_list_size;
      head_size = checksum_offset + gauge_check_size;
      
      if(gf->parallel)
	printf("%s: Attempting serial read from parallel file \n",myname);

      /* Allocate single precision read buffer */
      lbuf = (fsu3_matrix *)malloc(MAX_BUF_LENGTH*4*sizeof(fsu3_matrix));
      if(lbuf == NULL)
	{
	  printf("%s: Node %d can't malloc lbuf\n",myname,this_node);
	  fflush(stdout);
	  terminate(1);
	}
  
      /* Position file for reading gauge configuration */
      
      offset = head_size;

      if( g_seek(fp,offset,SEEK_SET) < 0 ) 
	{
	  printf("%s: Node 0 g_seek %lld failed error %d file %s\n",
		 myname,(long long)offset,errno,filename);
	  fflush(stdout);terminate(1);   
	}

      buf_length = 0;
      where_in_buf = 0;
      
    }

  /* all nodes initialize checksums */
  test_gc.sum29 = 0;
  test_gc.sum31 = 0;
  /* counts 32-bit words mod 29 and mod 31 in order of appearance
     on file */
  /* Here all nodes see the same sequence because we read serially */
  rank29 = 0;
  rank31 = 0;

  g_sync();

  /* Node 0 reads and deals out the values */

  for(rcv_rank=0; rcv_rank<volume; rcv_rank++)
    {
      /* If file is in coordinate natural order, receiving coordinate
         is given by rank. Otherwise, it is found in the table */

      if(gf->header->order == NATURAL_ORDER)
	rcv_coords = rcv_rank;
      else
	rcv_coords = gf->rank2rcv[rcv_rank];

      x = rcv_coords % nx;   rcv_coords /= nx;
      y = rcv_coords % ny;   rcv_coords /= ny;
      z = rcv_coords % nz;   rcv_coords /= nz;
      t = rcv_coords % nt;

      /* The node that gets the next set of gauge links */
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
	    
	    if( (int)g_read(lbuf,4*sizeof(fsu3_matrix),buf_length,fp) 
		!= buf_length)
	      {
		printf("%s: node %d gauge configuration read error %d file %s\n",
		       myname,this_node,errno,filename); 
		fflush(stdout); terminate(1);
	      }
	    where_in_buf = 0;  /* reset counter */
	  }  /*** end of the buffer read ****/

	if(destnode==0){	/* just copy links */
	  idest = node_index(x,y,z,t);
	  /* Save 4 matrices in tmpsu3 for further processing */
	  memcpy(tmpsu3,&lbuf[4*where_in_buf],4*sizeof(fsu3_matrix));
	}
	else {		/* send to correct node */
	  send_field((char *)&lbuf[4*where_in_buf],
		     4*sizeof(fsu3_matrix),destnode);
	}
	where_in_buf++;
      }
      
      /* The node that contains this site reads the message */
      else {	/* for all nodes other than node 0 */
	if(this_node==destnode){
	  idest = node_index(x,y,z,t);
	  /* Receive 4 matrices in temporary space for further processing */
	  get_field((char *)tmpsu3,4*sizeof(fsu3_matrix),0);
	}
      }

      /* The receiving node does the byte reversal and then checksum,
         if needed.  At this point tmpsu3 contains the input matrices
         and idest points to the destination site structure. */

      if(this_node==destnode)
	{
	  if(byterevflag==1)
	    byterevn((int32type *)tmpsu3,
		     4*sizeof(fsu3_matrix)/sizeof(int32type));
	  /* Accumulate checksums */
	  for(k = 0, val = (u_int32type *)tmpsu3; 
	      k < 4*(int)sizeof(fsu3_matrix)/(int)sizeof(int32type); 
	      k++, val++)
	    {
	      test_gc.sum29 ^= (*val)<<rank29 | (*val)>>(32-rank29);
	      test_gc.sum31 ^= (*val)<<rank31 | (*val)>>(32-rank31);
	      rank29++; if(rank29 >= 29)rank29 = 0;
	      rank31++; if(rank31 >= 31)rank31 = 0;
	    }
	  /* Copy 4 matrices to lattice[idest], converting to generic
	     precision */
	  f2d_4mat(tmpsu3,&lattice[idest].link[0]);
	}
      else
	{
	  rank29 += 4*sizeof(fsu3_matrix)/sizeof(int32type);
	  rank31 += 4*sizeof(fsu3_matrix)/sizeof(int32type);
	  rank29 %= 29;
	  rank31 %= 31;
	}
    }
  
  /* Combine node checksum contributions with global exclusive or */
  g_xor32(&test_gc.sum29);
  g_xor32(&test_gc.sum31);
  
  if(this_node==0)
    {
      /* Read and verify checksum */
      /* Checksums not implemented until version 5 */
      
      printf("Restored binary gauge configuration serially from file %s\n",
	     filename);
      if(gh->magic_number == GAUGE_VERSION_NUMBER)
	{
	  printf("Time stamp %s\n",gh->time_stamp);
	  if( g_seek(fp,checksum_offset,SEEK_SET) < 0 ) 
	    {
	      printf("%s: Node 0 g_seek %lld failed error %d file %s\n",
		    myname,(long long)offset,errno,filename);
	      fflush(stdout);terminate(1);   
	    }
	  read_checksum(SERIAL,gf,&test_gc);
	}
      else
	{
	  printf("Checksums %x %x\n",test_gc.sum29,test_gc.sum31);
	  printf("Checksums not verified in this format\n");
	}
      fflush(stdout);
      free(lbuf);
    }
  
} /* r_serial */

/*----------------------------------------------------------------------*/

static void r_serial_arch(gauge_file *gf)
{
  /* gf  = gauge configuration file structure */

  FILE *fp;
  char *filename;

  int rcv_rank, rcv_coords;
  int destnode;
  int i,k;
  int x,y,z,t;
  gauge_check test_gc;
  u_int32type *val;
  int rank29,rank31;
  su3_matrix tmpsu3[4];
  int dataformat = gf->dataformat;
  int precision = gf->precision;
  char myname[] = "r_serial_arch";
  int mu,a,b,p;
  float *uin = NULL, *q;
  double *uind = NULL, *qd;
  int big_end = 0;
  float U[4][18];
  double Ud[4][18];
  u_int32type chksum;
  int realspersite;

  if(dataformat == ARCHIVE_3x2)realspersite = 48;
  else realspersite = 72;
  fp = gf->fp;
  filename = gf->filename;

  if(this_node == 0)
    {
      if(gf->parallel)
	printf("%s: Attempting serial read from parallel file \n",myname);

      big_end = big_endian();
      /* printf("big_end is %d\n", big_end); */
      if(precision == 1){
	uin = (float *) malloc(nx*ny*nz*realspersite*sizeof(float));
	if(uin == NULL)
	  {
	    printf("%s: Node %d can't malloc uin buffer to read timeslice\n",
		   myname,this_node);
	    printf("recompile with smaller read buffer: uin\n");
	    fflush(stdout);
	    terminate(1);
	  }
      } else {
	uind = (double *) malloc(nx*ny*nz*realspersite*sizeof(double));
	if(uind == NULL)
	  {
	    printf("%s: Node %d can't malloc uind buffer to read timeslice\n",
		   myname,this_node);
	    printf("recompile with smaller read buffers: uin\n");
	    fflush(stdout);
	    terminate(1);
	  }
      }
    }

  /* Initialize checksums */
  chksum = 0;
  test_gc.sum29 = 0;
  test_gc.sum31 = 0;
  /* counts 32-bit words mod 29 and mod 31 in order of appearance
     on file */
  /* Here all nodes see the same sequence because we read serially */
  rank29 = 0;
  rank31 = 0;

  g_sync();

  /* Node 0 reads and deals out the values */
  for(rcv_rank=0; rcv_rank<volume; rcv_rank++)
    {
      rcv_coords = rcv_rank;

      x = rcv_coords % nx;   rcv_coords /= nx;
      y = rcv_coords % ny;   rcv_coords /= ny;
      z = rcv_coords % nz;   rcv_coords /= nz;
      t = rcv_coords % nt;

      /* The node that gets the next set of gauge links */
      destnode=node_number(x,y,z,t);
      
      if(this_node==0){
	if(precision == 1){
	  if( (int)g_read(uin,realspersite*sizeof(float),1,fp) != 1)
	    {
	      printf("%s: node %d gauge configuration read error %d file %s\n",
		     myname,this_node,errno,filename); 
	      fflush(stdout); terminate(1);
	    }
	  if (!big_end) byterevn((int32type *)uin,realspersite);
	  q = uin;
	  for (mu=0;mu<4;mu++) {
	    for (p=0;p<realspersite/4;p++) {
	      chksum += *(u_int32type *) q;
	      U[mu][p] = (float) *(q++);
	    }
	    if(dataformat == ARCHIVE_3x2)
	      complete_U(U[mu]);
	    /**
	       for (p=0;p<18;p++) printf("p=%d, e=%f\n", p, U[mu][p]);
	    **/
	    
	    /* Copy, converting precision if necessary */
	    for(a=0; a<3; a++) for(b=0; b<3; b++) { 
	      tmpsu3[mu].e[a][b].real = U[mu][2*(3*a+b)];
	      /* printf("real: p=%d, mu=%d, e=%f\n", p,mu,U[mu][2*(3*a+b)]); */
	      tmpsu3[mu].e[a][b].imag = U[mu][2*(3*a+b)+1];
	      /*printf("imag: p=%d, mu=%d, e=%f\n", p,mu,U[mu][2*(3*a+b)+1]); */
	    }
	  }
	} else { /* precision == 2 */
	  
	  if( (int)g_read(uind,realspersite*sizeof(double),1,fp) != 1)
	    {
	      printf("%s: node %d gauge configuration read error %d file %s\n",
		     myname,this_node,errno,filename); 
	      fflush(stdout); terminate(1);
	    }
	  if (!big_end) byterevn64((int32type *)uind,realspersite);
	  qd = uind;
	  for (mu=0;mu<4;mu++) {
	    for (p=0;p<realspersite/4;p++) {
	      chksum += *(u_int32type *) qd;
	      chksum += *((u_int32type *) qd + 1);
	      Ud[mu][p] = (double) *(qd++);
	    }
	    if(dataformat == ARCHIVE_3x2)
	      complete_Ud(Ud[mu]);
	    /**
	       for (p=0;p<18;p++) printf("p=%d, e=%f\n", p, Ud[mu][p]);
	    **/
	    
	    /* Copy, converting precision if necessary */
	    for(a=0; a<3; a++) for(b=0; b<3; b++) { 
	      tmpsu3[mu].e[a][b].real = Ud[mu][2*(3*a+b)];
	      /* printf("real: p=%d, mu=%d, e=%f\n", p,mu,Ud[mu][2*(3*a+b)]); */
	      tmpsu3[mu].e[a][b].imag = Ud[mu][2*(3*a+b)+1];
	      /*printf("imag: p=%d, mu=%d, e=%f\n", p,mu,Ud[mu][2*(3*a+b)+1]); */
	    }
	  } 
	}
	
	if(destnode==0){	/* just copy links */
	  i = node_index(x,y,z,t);
     /*   printf("lattice node_index = %d, mu = %d\n", i, mu); */
	  /* Copy from tmpsu3 to site structure */
	  memcpy((void *)lattice[i].link,(void *)tmpsu3,4*sizeof(su3_matrix));
	} else {		
	  /* send to correct node */
	  send_field((char *)tmpsu3, 4*sizeof(su3_matrix),destnode);
	}
      } 
      /* The node which contains this site reads message */
      else {	
	/* for all nodes other than node 0 */
	if(this_node==destnode){
	  i = node_index(x,y,z,t);
	  get_field((char *)tmpsu3,4*sizeof(su3_matrix),0);
	  /* Store in site structure, converting to generic precision */
	  memcpy((void *)lattice[i].link,(void *)tmpsu3,4*sizeof(su3_matrix));
	}
      }

      /* Any needed byte reversing was already done. Compute MILC
         checksums. At this point tmpsu3 on destnode contains the link
         matrices we just read */

      if(this_node==destnode)
	{
	  /* Accumulate checksums */
	  for(k = 0, val = (u_int32type *)tmpsu3;
	      k < 4*(int)sizeof(su3_matrix)/(int)sizeof(int32type); k++, val++)
   	    {
	      test_gc.sum29 ^= (*val)<<rank29 | (*val)>>(32-rank29);
	      test_gc.sum31 ^= (*val)<<rank31 | (*val)>>(32-rank31);
	      rank29++; if(rank29 >= 29)rank29 = 0;
	      rank31++; if(rank31 >= 31)rank31 = 0;
	    }
	}
      else
	{
	  rank29 += 4*sizeof(su3_matrix)/sizeof(int32type);
	  rank31 += 4*sizeof(su3_matrix)/sizeof(int32type);
	  rank29 %= 29;
	  rank31 %= 31;
	}
    }
  
  /* Combine node checksum contributions with global exclusive or */
  g_xor32(&test_gc.sum29);
  g_xor32(&test_gc.sum31);
  
  if(this_node==0)
    {
      /* Read and verify checksum */
      
      printf("Restored archive gauge configuration serially from file %s\n",
	     filename);
      if (chksum != gf->check.sum31)
	{
	  printf("Archive style checksum violation: computed %x, read %x\n",
		 chksum, gf->check.sum31);
	}
      else
	{
	  printf("Archive style checksum = %x OK\n", chksum);
	}
      fflush(stdout);
      free(uin);

      /* Store MILC style checksums */
      gf->check.sum29 = test_gc.sum29;
      gf->check.sum31 = test_gc.sum31;
    }
  
} /* r_serial_arch */





/*---------------------------------------------------------------------------*/
/* Write parallel gauge configuration in coordinate natural order */

static void w_parallel(gauge_file *gf)
{
  /* gf  = file descriptor as opened by w_parallel_i */

  FILE *fp;
  fsu3_matrix *lbuf;
  int buf_length,where_in_buf;
  u_int32type *val;
  int rank29,rank31;
  off_t checksum_offset;
  register int i;
  int j,k;
  int x,y,z,t;
  struct {
    short x,y,z,t;
    fsu3_matrix link[4];
  } msg;
  int isite,ksite,site_block;
  int rcv_coords,rcv_rank;
  int destnode,sendnode;
  char myname[] = "w_parallel";

  fp = gf->fp;

  lbuf = w_parallel_setup(gf,&checksum_offset);

  /* Collect buffer from other nodes and write when full */

  /* initialize checksums */
  gf->check.sum31 = 0;
  gf->check.sum29 = 0;

  /* Read and deal */

  g_sync();
  buf_length = 0;

  /* Clear buffer as a precaution.  Easier to tell if we botch the
     buffer loading. */
  for(i=0;i<MAX_BUF_LENGTH;i++)
    for(j=0;j<3;j++)for(k=0;k<3;k++)
      { lbuf[i].e[j][k].real = lbuf[i].e[j][k].imag = 0.;}
  
  /* Cycle through nodes, collecting a buffer full of values from the
     appropriate node before proceeding to the next node in sequence.
     We don't know if this pattern is generally optimal.  It is
     possible that messages arrive at a node in an order different
     from the order of sending so we include the site coordinates in
     the message to be sure it goes where it belongs */
  
  /* MUST be a factor of MAX_BUF_LENGTH */
  site_block = MAX_BUF_LENGTH;  
  if(MAX_BUF_LENGTH % site_block != 0)
    {printf("%s: site_block incommensurate with buffer size\n",myname);
     fflush(stdout);terminate(1);}

  for(ksite=0; ksite<sites_on_node; ksite += site_block)
    {
      for(destnode=0; destnode<number_of_nodes; destnode++)
	for(isite=ksite; 
	    isite<sites_on_node && isite<ksite+site_block; isite++)
	  {
	    
	    /* This is the coordinate natural (typewriter) rank
	       of the site the destnode needs next */
	    
	    rcv_rank = rcv_coords = destnode*sites_on_node + isite;
	    
	    /* The coordinate corresponding to this site */
	    
	    x = rcv_coords % nx; rcv_coords /= nx;
	    y = rcv_coords % ny; rcv_coords /= ny;
	    z = rcv_coords % nz; rcv_coords /= nz;
	    t = rcv_coords % nt;
	    
	    /* The node that has this site */
	    sendnode=node_number(x,y,z,t);
	    
	    /* Node sendnode sends site value to destnode */
	    if(this_node==sendnode && destnode!=sendnode){
	      /* Message consists of site coordinates and 4 link matrices */
	      msg.x = x; msg.y = y; msg.z = z; msg.t = t;
	      i = node_index(x,y,z,t);
	      /* Copy 4 matrices and convert to single precision msg
		 structure */
	      d2f_4mat(&lattice[i].link[0],&msg.link[0]);

	      send_field((char *)&msg,sizeof(msg),destnode);
	    }
	    /* Node destnode copies or receives a message */
	    else if(this_node==destnode){
	      if(destnode==sendnode){ 
		/* just copy links to write buffer */
		i = node_index(x,y,z,t);
		where_in_buf = buf_length;
		d2f_4mat(&lattice[i].link[0],&lbuf[4*where_in_buf]);
		rank29 = rank31 = 
		  4*sizeof(fsu3_matrix)/sizeof(int32type)*rcv_rank;
	      }
	      else {
		/* Receive a message */
		/* Note that messages may arrive in any order
		   so we use the x,y,z,t coordinate to tell
		   where it goes in the write buffer */
		get_field((char *)&msg,sizeof(msg),sendnode);
		/* Reconstruct rank from message coordinates */
		i = msg.x+nx*(msg.y+ny*(msg.z+nz*msg.t));
		/* The buffer location is then */
		where_in_buf = (i % sites_on_node) % MAX_BUF_LENGTH;

		/* Move data to buffer */
		memcpy((void *)&lbuf[4*where_in_buf],
		       (void *)msg.link,4*sizeof(fsu3_matrix));
		rank29 = rank31 = 4*sizeof(fsu3_matrix)/sizeof(int32type)*i;
	      }

	      /* Receiving node accumulates checksums as the values
		 are inserted into its buffer */
	      rank29 %= 29; rank31 %= 31;
	      for(k = 0, val = (u_int32type *)&lbuf[4*where_in_buf]; 
		  k < 4*(int)sizeof(fsu3_matrix)/(int)sizeof(int32type); k++, val++)
		{
		  gf->check.sum29 ^= (*val)<<rank29 | (*val)>>(32-rank29);
		  gf->check.sum31 ^= (*val)<<rank31 | (*val)>>(32-rank31);
		  rank29++; if(rank29 >= 29)rank29 = 0;
		  rank31++; if(rank31 >= 31)rank31 = 0;
		}

	      buf_length++;
	      if( (buf_length == MAX_BUF_LENGTH) || 
		 (isite == sites_on_node -1))
		{
		  /* write out buffer */
		  
		  if( (int)g_write(lbuf,4*sizeof(fsu3_matrix),buf_length,fp) 
		     != buf_length)
		    {
		      printf("%s: Node %d gauge configuration write error %d file %s\n",
			     myname,this_node,errno,gf->filename); 
		      fflush(stdout);
		      terminate(1);   
		    }
		  buf_length = 0;		/* start again after write */
		  /* Clear buffer as a precaution */
		  for(i=0;i<MAX_BUF_LENGTH;i++)
		    for(j=0;j<3;j++)for(k=0;k<3;k++)
		      { lbuf[i].e[j][k].real = lbuf[i].e[j][k].imag = 0.;}
		}
	    } /* else if(this_node==destnode) */
	    
	  } /* destnode, isite */
      g_sync();  /* To assure all write buffers are completed before
		   starting on the next buffer */
    } /* ksite */
  
  free(lbuf);

  /* Combine checksums */

  g_xor32(&gf->check.sum29);
  g_xor32(&gf->check.sum31);

  /* Write checksum at end of lattice file */

  /* Position file for writing checksum */
  /* Only node 0 writes checksum data */
      
  if(this_node==0){
    if( g_seek(fp,checksum_offset,SEEK_SET) < 0 ) 
      {
	printf("%s: Node %d g_seek %ld for checksum failed error %d file %s\n",
	       myname,this_node,(long)checksum_offset,errno,gf->filename);
	fflush(stdout);terminate(1);   
      }

    write_checksum(PARALLEL,gf);

    printf("Saved gauge configuration in parallel to binary file %s\n",
	   gf->filename);
    printf("Time stamp %s\n",(gf->header)->time_stamp);
    
  }

} /* w_parallel */

/*-----------------------------------------------------------------------*/

/* Write parallel gauge configuration in node dump order */

static void w_checkpoint(gauge_file *gf)
{
  /* gf  = file descriptor as opened by w_checkpoint_i */

  FILE *fp;
  fsu3_matrix *lbuf;
  u_int32type *val;
  int k;
  int rank29,rank31;
  off_t checksum_offset;
  int buf_length;
  register site *s;
  register int i;
  char myname[] = "w_checkpoint";

  fp = gf->fp;

  lbuf = w_parallel_setup(gf,&checksum_offset);

  /* C. McNeile's algorithm, changed slightly*/

  /* initialize checksums */
  gf->check.sum31 = 0;
  gf->check.sum29 = 0;
  /* counts 32-bit words mod 29 and mod 31 in order of appearance on file */
  /* Here all nodes use these values */
  rank29 = 4*sizeof(fsu3_matrix)/sizeof(int32type)*sites_on_node*this_node % 29;
  rank31 = 4*sizeof(fsu3_matrix)/sizeof(int32type)*sites_on_node*this_node % 31;

  buf_length = 0;

  FORALLSITES(i,s)
  {
        
    /* load the gauge configuration into the buffer */
    /* convert (copy) generic to single precision */
    d2f_4mat(&lattice[i].link[0],&lbuf[4*buf_length]);

    /* Accumulate checksums - contribution from next site moved into buffer*/
    for(k = 0, val = (u_int32type *)&lbuf[4*buf_length]; 
	k < 4*(int)sizeof(fsu3_matrix)/(int)sizeof(int32type); k++, val++)
      {
	gf->check.sum29 ^= (*val)<<rank29 | (*val)>>(32-rank29);
	gf->check.sum31 ^= (*val)<<rank31 | (*val)>>(32-rank31);
	rank29++; if(rank29 >= 29)rank29 = 0;
	rank31++; if(rank31 >= 31)rank31 = 0;
      }

    buf_length++;
    
    if( (buf_length == MAX_BUF_LENGTH) || (i == sites_on_node -1))
      {
	/* write out buffer */
	
	fflush(stdout);
	if( (int)g_write(lbuf,4*sizeof(fsu3_matrix),buf_length,fp) != buf_length)
	  {
	    printf("%s: Node %d gauge configuration write error %d file %s\n",
		   myname,this_node,errno,gf->filename); 
	    fflush(stdout);
	    terminate(1);   
	  }
	buf_length = 0;		/* start again after write */
      }
    
  } 
  
  free(lbuf);

  /* Combine checksums */

  g_xor32(&gf->check.sum29);
  g_xor32(&gf->check.sum31);

  /* Write checksum at end of lattice file */

  /* Position file for writing checksum */
  /* Only node 0 writes checksum data */
      
  if(this_node == 0)
    {
      if( g_seek(fp,checksum_offset,SEEK_SET) < 0 ) 
	{
	  printf("%s: Node %d g_seek %ld for checksum failed error %d file %s\n",
		 myname,this_node,(long)checksum_offset,errno,gf->filename);
	  fflush(stdout);terminate(1);   
	}

      write_checksum(PARALLEL,gf);
      
      printf("Saved gauge configuration checkpoint file %s\n",
	     gf->filename);
      printf("Time stamp %s\n",(gf->header)->time_stamp);
    }

} /* w_checkpoint */

/*---------------------------------------------------------------------------*/

static gauge_file *r_parallel_i(char *filename)
{
  /* Returns file descriptor for opened file */

  gauge_header *gh;
  gauge_file *gf;
  FILE *fp;
  int byterevflag;

  /* All nodes set up a gauge file and guage header structure for reading */

  gf = setup_input_gauge_file(filename);
  gh = gf->header;

  gf->parallel = 1;   /* File was opened for parallel access */

  /* All nodes open a file */

  fp = g_open(filename, "rb");
  if(fp == NULL)
    {
      printf("r_parallel_i: Node %d can't open file %s, error %d\n",
	     this_node,filename,errno);fflush(stdout);terminate(1);
    }

  gf->fp = fp;

  /* Node 0 reads header */

  if(this_node==0)
    byterevflag = read_gauge_hdr(gf,PARALLEL);
  
  /* Broadcast the byterevflag from node 0 to all nodes */

  broadcast_bytes((char *)&byterevflag,sizeof(byterevflag));

  gf->byterevflag = byterevflag;

  /* Broadcasts the header structure from node 0 to all nodes */
  
  broadcast_bytes((char *)gh,sizeof(gauge_header));

  /* Read site list and broadcast to all nodes */

  read_site_list(PARALLEL,gf);

  return gf;

} /* r_parallel_i */

/*----------------------------------------------------------------------*/

/* Read gauge configuration in parallel from a single file */
static void r_parallel(gauge_file *gf)
{
  /* gf  = gauge configuration file structure */

  FILE *fp;
  gauge_header *gh;
  char *filename;
  fsu3_matrix *lbuf;
  struct {
    short x,y,z,t;
    fsu3_matrix link[4];
  } msg;

  int buf_length,where_in_buf;
  gauge_check test_gc;
  u_int32type *val;
  int rank29,rank31;
  int destnode,sendnode,isite,ksite,site_block;
  int x,y,z,t;
  int rcv_rank,rcv_coords;
  register int i,k;

  off_t offset ;            /* File stream pointer */
  off_t gauge_node_size;   /* Size of a gauge configuration block for
                              all sites on one node */
  off_t gauge_check_size;  /* Size of gauge configuration checksum record */
  off_t coord_list_size;    /* Size of coordinate list in bytes */
  off_t head_size;          /* Size of header plus coordinate list */
  off_t checksum_offset;    /* Where we put the checksum */
  char myname[] = "r_parallel";

  fp = gf->fp;
  gh = gf->header;

  filename = gf->filename;

  if(!gf->parallel)
    printf("%s: Attempting parallel read from serial file.\n",myname);

  /* Allocate single precision read buffer */
  lbuf = (fsu3_matrix *)malloc(MAX_BUF_LENGTH*4*sizeof(fsu3_matrix));
  if(lbuf == NULL)
    {
      printf("%s: Node %d can't malloc lbuf\n",myname,this_node); 
      fflush(stdout);terminate(1);
    }

  gauge_node_size = sites_on_node*4*sizeof(fsu3_matrix) ;

  /* (1996 gauge configuration files had a 32-bit unused checksum 
     record before the gauge link data) */
  if(gh->magic_number == GAUGE_VERSION_NUMBER)
    gauge_check_size = sizeof(gf->check.sum29) + 
      sizeof(gf->check.sum31);
  else
    gauge_check_size = 0;

  if(gf->header->order == NATURAL_ORDER)coord_list_size = 0;
  else coord_list_size = sizeof(int32type)*volume;
  checksum_offset = gf->header->header_bytes + coord_list_size;
  head_size = checksum_offset + gauge_check_size;

  offset = head_size;

  /* Position file for reading gauge configuration */
  /* Each node reads */

  offset += gauge_node_size*this_node;
  
  if( g_seek(fp,offset,SEEK_SET) < 0 ) 
    {
      printf("%s: Node %d g_seek %ld failed error %d file %s\n",
	     myname,this_node,(long)offset,errno,filename);
      fflush(stdout);terminate(1);   
    }

  /* initialize checksums */
  test_gc.sum29 = 0;
  test_gc.sum31 = 0;
  /* counts 32-bit words mod 29 and mod 31 in order of appearance on file */
  /* Here all nodes use these values */
  rank29 = 4*sizeof(fsu3_matrix)/sizeof(int32type)*sites_on_node*this_node % 29;
  rank31 = 4*sizeof(fsu3_matrix)/sizeof(int32type)*sites_on_node*this_node % 31;

  /* Read and deal */

  g_sync();
  buf_length = 0;
  where_in_buf = 0;
  
  /* Cycle through nodes, dealing 4 values from each node in sequence.
     (We don't know if this pattern is generally optimal.)

     It is possible that messages arrive at a node in an order
     different from the order of dealing so we include the site
     coordinates in the message to specify where it goes */
  
  site_block = 4;
  for(ksite=0; ksite<sites_on_node; ksite += site_block)
    {
    for(sendnode=0; sendnode<number_of_nodes; sendnode++)
      for(isite=ksite; 
	  isite<sites_on_node && isite<ksite+site_block; isite++)
	{
	  /* Compute destination coordinate for the next field 
	     
	     In coordinate natural order (typewriter order)
	     the rank order of data for site (x,y,z,t) on the 
	     file is given by
	     
	     rcv_coords = x+nx*(y+ny*(z+nz*t))
	     
	     For purposes of reading, the data is divided
	     equally among the nodes with node 0 taking the 1st block,
	     node 1 the second, etc.  */
	  
	  rcv_rank = sendnode*sites_on_node + isite;
	  
	  /* If sites are not in natural order, use the
	     site list */
	  
	  if(gf->header->order == NATURAL_ORDER)
	    rcv_coords = rcv_rank;
	  else
	    rcv_coords = gf->rank2rcv[rcv_rank];
	  
	  x = rcv_coords % nx; rcv_coords /= nx;
	  y = rcv_coords % ny; rcv_coords /= ny;
	  z = rcv_coords % nz; rcv_coords /= nz;
	  t = rcv_coords % nt;
	  
	  
	  /* Destination node for this value */
	  destnode=node_number(x,y,z,t);
	  
	  /* Node sendnode reads, and sends site to correct node */
	  if(this_node==sendnode){
	    
	    if(where_in_buf == buf_length)
	      
	      {  /* get new buffer */
		
		/* new buffer length  = remaining sites, but never bigger 
		   than MAX_BUF_LENGTH */
		buf_length = sites_on_node - isite;
		if(buf_length > MAX_BUF_LENGTH) buf_length = MAX_BUF_LENGTH; 
		/* then do read */
		/* each node reads its sites */
		
		if( g_read(lbuf,buf_length*4*sizeof(fsu3_matrix),1,fp) != 1)
		  {
		    printf("%s: node %d gauge configuration read error %d file %s\n",
			   myname,this_node,errno,filename); 
		    fflush(stdout); terminate(1);
		  }
		where_in_buf = 0;  /* reset counter */
	      }  /*** end of the buffer read ****/
	    
	    /* Do byte reversal if needed */
	    if(gf->byterevflag==1)
	      byterevn((int32type *)&lbuf[4*where_in_buf],
		       4*sizeof(fsu3_matrix)/sizeof(int32type));

	    /* Accumulate checksums - contribution from next site */
	    for(k = 0, val = (u_int32type *)&lbuf[4*where_in_buf]; 
		k < 4*(int)sizeof(fsu3_matrix)/(int)sizeof(int32type); k++, val++)
	      {
		test_gc.sum29 ^= (*val)<<rank29 | (*val)>>(32-rank29);
		test_gc.sum31 ^= (*val)<<rank31 | (*val)>>(32-rank31);
		rank29++; if(rank29 >= 29)rank29 = 0;
		rank31++; if(rank31 >= 31)rank31 = 0;
	      }

	    if(destnode==sendnode){	
	      /* just copy links, converting to generic precision */
	      i = node_index(x,y,z,t);
	      f2d_4mat((fsu3_matrix *)&lbuf[4*where_in_buf],
		       &lattice[i].link[0]);
	    }
	    else {		
	      /* send to correct node */
	      /* Message consists of site coordinates and 4 link matrices */
	      msg.x = x; msg.y = y; msg.z = z; msg.t = t;
	      memcpy((void *)msg.link,
		     (void *)&lbuf[4*where_in_buf],4*sizeof(fsu3_matrix));
	      
	      send_field((char *)&msg,sizeof(msg),destnode);
	    }
	    where_in_buf++;
	  }
	  /* The node which contains this site reads a message */
	  else {	/* for all nodes other than node sendnode */
	    if(this_node==destnode){
	      get_field((char *)&msg,sizeof(msg),sendnode);
	      i = node_index(msg.x,msg.y,msg.z,msg.t);
	      if(this_node!= node_number(msg.x,msg.y,msg.z,msg.t))
		{
		  printf("BOTCH. Node %d received %d %d %d %d\n",
			 this_node,msg.x,msg.y,msg.z,msg.t);
		  fflush(stdout); terminate(1);
		}
	      /* Store in the proper location, converting to generic
		 precision */
	      f2d_4mat(&msg.link[0],&lattice[i].link[0]);
	    }
	  }
	} /** end over the lattice sites in block on all nodes ***/

    g_sync(); /* To prevent incoming message pileups */
  }  /** end over blocks **/

  free(lbuf);

  /* Combine node checksum contributions with global exclusive or */
  g_xor32(&test_gc.sum29);
  g_xor32(&test_gc.sum31);

  /* Read and verify checksum */
  
  if(this_node == 0)
    {
      /* Node 0 positions file for reading checksum */
      
      printf("Restored binary gauge configuration in parallel from file %s\n",
	       filename);
      if(gh->magic_number == GAUGE_VERSION_NUMBER)
	{
	  printf("Time stamp %s\n",gh->time_stamp);
	  if( g_seek(fp,checksum_offset,SEEK_SET) < 0 ) 
	    {
	      printf("%s: Node 0 g_seek %ld for checksum failed error %d file %s\n",
		     myname,(long)offset,errno,filename);
	      fflush(stdout);terminate(1);   
	    }
	  
	  read_checksum(PARALLEL,gf,&test_gc);
	}
      fflush(stdout);
    }  
  
} /* r_parallel */

/*---------------------------------------------------------------------------*/
/* Top level routines */
/*---------------------------------------------------------------------------*/

/* Read a lattice in ASCII format serially (node 0 only) */

/* format
    version_number (int)
    time_stamp (char string enclosed in quotes)
    nx ny nz nt (int)
    for(t=...)for(z=...)for(y=...)for(x=...){
	xlink,ylink,zlink,tlink
    }
        for each link:
            for(i=...)for(j=...){link[i][j].real, link[i][j].imag}

    A separate ASCII info file is also written.
*/
gauge_file *restore_ascii(char *filename) {
  gauge_header *gh;
  gauge_file *gf;
  FILE *fp = NULL;
  int destnode;
  int version_number,i,j,x,y,z,t,dir;
  fsu3_matrix lbuf[4];
  
  /* Set up a gauge file and guage header structure for reading */

  gf = setup_input_gauge_file(filename);
  gh = gf->header;

  /* File opened for serial reading */
  gf->parallel = 0;

  /* Node 0 opens the file and reads the header */

  if(this_node==0){
    fp = fopen(filename,"r");
    if(fp==NULL){
      printf("Can't open file %s, error %d\n",filename,errno);
      terminate(1);
    }

    gf->fp = fp;

    if( (fscanf(fp,"%d",&version_number))!=1 ){
      printf("restore_ascii: Error reading version number\n"); terminate(1);
    }
    gh->magic_number = version_number;
    if(gh->magic_number != GAUGE_VERSION_NUMBER){
      printf("restore_ascii: Incorrect version number in lattice header\n");
      printf("  read %d but expected %d\n",
	     gh->magic_number,GAUGE_VERSION_NUMBER);
      terminate(1);
    }
    /* Time stamp is enclosed in quotes - discard the leading white
       space and the quotes and read the enclosed string */
    if((i = fscanf(fp,"%*[ \f\n\r\t\v]%*[\"]%[^\"]%*[\"]",gh->time_stamp))!=1){
      printf("restore_ascii: Error reading time stamp\n"); 
      printf("count %d time_stamp %s\n",i,gh->time_stamp);
      terminate(1);
    }
    if( (fscanf(fp,"%d%d%d%d",&x,&y,&z,&t))!=4 ){
      printf("restore_ascii: Error in reading dimensions\n"); terminate(1);
    }
    gh->dims[0] = x; gh->dims[1] = y; gh->dims[2] = z; gh->dims[3] = t;
    if( gh->dims[0]!=nx || gh->dims[1]!=ny || 
       gh->dims[2]!=nz || gh->dims[3]!=nt )
      {
	/* So we can use this routine to discover the dimensions,
	   we provide that if nx = ny = nz = nt = -1 initially
	   we don't die */
	if(nx != -1 || ny != -1 || nz != -1 || nt != -1)
	  {
	    printf("restore_ascii: Incorrect lattice size %d,%d,%d,%d\n",
		   gh->dims[0],gh->dims[1],gh->dims[2],gh->dims[3]);
	    terminate(1);
	  }
	else
	  {
	    nx = gh->dims[0];
	    ny = gh->dims[1];
	    nz = gh->dims[2];
	    nt = gh->dims[3];
	    volume = nx*ny*nz*nt;
	  }
      }

    gh->order = NATURAL_ORDER;           /* (Not used) */

  } /* if node 0 */
  
  else gf->fp = NULL;

  gf->byterevflag = 0;    /* (Not used) */

  /* Node 0 broadcasts the header structure to all nodes */
  
  broadcast_bytes((char *)gh,sizeof(gauge_header));

  /* Read gauge field values */  
  g_sync();
  
  for(t=0;t<nt;t++)for(z=0;z<nz;z++)for(y=0;y<ny;y++)for(x=0;x<nx;x++){
    destnode=node_number(x,y,z,t);
    
    /* Node 0 reads, and sends site to correct node */
    if(this_node==0){
      for(dir=XUP;dir<=TUP;dir++){
	for(i=0;i<3;i++)for(j=0;j<3;j++){
	  if( fscanf(fp,"%e%e\n",
		     &(lbuf[dir].e[i][j].real),
		     &(lbuf[dir].e[i][j].imag) )!= 2)
	    {
	    printf("restore_ascii: gauge link read error\n"); 
	    terminate(1);
	  }
	}
      }
      if(destnode==0){	/* just copy links */
	i = node_index(x,y,z,t);
	f2d_4mat(lbuf, lattice[i].link);
      }
      else {		/* send to correct node */
	send_field((char *)lbuf,4*sizeof(fsu3_matrix),destnode);
      }
    }
    
    /* The node which contains this site reads message */
    else {	/* for all nodes other than node 0 */
      if(this_node==destnode){
	get_field((char *)lbuf,4*sizeof(fsu3_matrix),0);
	i = node_index(x,y,z,t);
	f2d_4mat(lbuf, lattice[i].link);
      }
    }
  }
  
  g_sync();
  if(this_node==0){
    printf("Restored gauge configuration from ascii file  %s\n",
	   filename);
    printf("Time stamp %s\n",gh->time_stamp);
    g_close(fp);
    gf->fp = NULL;
    fflush(stdout);
  }

  return gf;
}

/*---------------------------------------------------------------------------*/

/* Save a lattice in ASCII format serially (node 0 only) */

gauge_file *save_ascii(char *filename) {
  FILE *fp = NULL;
  int currentnode,newnode;
  int i,j,x,y,z,t,dir;
  fsu3_matrix lbuf[4];
  gauge_file *gf;
  gauge_header *gh;

  /* Set up gauge file and gauge header structures and load header values */
  gf = setup_output_gauge_file();
  gh = gf->header;

  /* node 0 does all the writing */
  if(this_node==0){

    fp = fopen(filename,"w");
    if(fp==NULL){
      printf("Can't open file %s, error %d\n",filename,errno);terminate(1);
    }

    gf->fp = fp;
    gf->parallel = 0;
    gf->filename        = filename;
    gf->byterevflag    = 0;            /* Not used for writing */

    if( (fprintf(fp,"%d\n",GAUGE_VERSION_NUMBER))==0 ){
      printf("Error in writing version number\n"); terminate(1);
    }
    if( (fprintf(fp,"\"%s\"\n",gh->time_stamp))==0 ){
      printf("Error in writing time stamp\n"); terminate(1);
    }
    
    if( (fprintf(fp,"%d\t%d\t%d\t%d\n",nx,ny,nz,nt))==0 ){
      printf("Error in writing dimensions\n"); terminate(1);
    }

    write_gauge_info_file(gf);
  }

  /* Write gauge field */

  g_sync();
  currentnode=0;
  
  for(t=0;t<nt;t++)for(z=0;z<nz;z++)for(y=0;y<ny;y++)for(x=0;x<nx;x++){
    newnode=node_number(x,y,z,t);
    if(newnode != currentnode){	/* switch to another node */
      /**g_sync();**/
      /* tell newnode it's OK to send */
      if( this_node==0 && newnode!=0 )send_field((char *)lbuf,32,newnode);
      if( this_node==newnode && newnode!=0 )get_field((char *)lbuf,32,0);
      currentnode=newnode;
    }
    
    if(this_node==0){
      if(currentnode==0){
	i=node_index(x,y,z,t);
	d2f_4mat(lattice[i].link, lbuf);
      }
      else{
	get_field((char *)lbuf,4*sizeof(fsu3_matrix),currentnode);
      }
      for(dir=XUP;dir<=TUP;dir++){
	for(i=0;i<3;i++)for(j=0;j<3;j++){
	  if( (fprintf(fp,"%.7e\t%.7e\n",(float)lbuf[dir].e[i][j].real,
		       (float)lbuf[dir].e[i][j].imag))== EOF){
	    printf("Write error in save_ascii\n"); terminate(1);
	  }
	}
      }
    }
    else {	/* for nodes other than 0 */
      if(this_node==currentnode){
	i=node_index(x,y,z,t);
	d2f_4mat(lattice[i].link, lbuf);
	send_field((char *)lbuf,4*sizeof(fsu3_matrix),0);
      }
    }
  }
  g_sync();
  if(this_node==0){
    fflush(fp);
    printf("Saved gauge configuration to ascii file %s\n",
	   gf->filename);
    printf("Time stamp %s\n",gh->time_stamp);
    g_close(fp);
    fflush(stdout);
    }
  return gf;
}

/*---------------------------------------------------------------------*/
/* Restore lattice file by reading serially (node 0 only) */
/* Handles most lattice formats */
  
gauge_file *restore_serial(char *filename)
{
  gauge_file *gf;

  gf = r_serial_i(filename);
  if(gf->header->magic_number == GAUGE_VERSION_NUMBER_ARCHIVE) 
    {
      r_serial_arch(gf);
      r_serial_f(gf);
    } 
  else if(gf->header->magic_number == LIME_MAGIC_NO)
    {
      r_serial_f(gf);
      /* Close this reader and reread to get the header */
      free(gf->header);
      free(gf);
#ifdef HAVE_QIO
      gf = restore_serial_scidac(filename);
#else
      node0_printf("Looks like a SciDAC file.  Recompile with QIO.\n");
      terminate(1);
#endif
    }
  else
    {
      r_serial(gf);
      r_serial_f(gf);
    }

  return gf;
  
} /* restore_serial */

/*---------------------------------------------------------------------------*/
/* Restore lattice file by reading to all nodes simultaneously */
/* Handles most lattice formats */
  
gauge_file *restore_parallel(char *filename)
{
  gauge_file *gf;

  gf = r_parallel_i(filename);
  if(gf->header->magic_number == LIME_MAGIC_NO)
    {
      r_serial_f(gf);
      /* Close this reader and reread to get the header */
      free(gf->header);
      free(gf);
#ifdef HAVE_QIO
      gf = restore_parallel_scidac(filename);
#else
      node0_printf("Looks like a SciDAC file.  Recompile with QIO.\n");
      terminate(1);
#endif
    }
  else
    {
      r_parallel(gf);
      r_parallel_f(gf);
    }

  return gf;
  
} /* restore_parallel */

/*---------------------------------------------------------------------------*/

/* Save lattice in natural order by writing serially (node 0 only) */

gauge_file *save_serial(char *filename)
{
  gauge_file *gf;

  gf = w_serial_i(filename);
  w_serial(gf);
  w_serial_f(gf);

  return gf;

} /* save_serial */

/*---------------------------------------------------------------------------*/

/* Save lattice in natural order by writing from all nodes at once */

gauge_file *save_parallel(char *filename)
{
  gauge_file *gf;

  gf = w_parallel_i(filename);
  w_parallel(gf);
  w_parallel_f(gf);

  return gf;

} /* save_parallel */

/*---------------------------------------------------------------------------*/

/* Save lattice in node-dump order */

/* This is much faster than save_parallel.  Lattices in this format
   can also be read much more quickly to the same number of nodes and
   layout. However, we probably wouldn't share lattices written in this
   order with our friends. */

gauge_file *save_checkpoint(char *filename)
{
  gauge_file *gf;

  gf = w_checkpoint_i(filename);
  w_checkpoint(gf);
  w_parallel_f(gf);

  return gf;

} /* save_checkpoint */

/*---------------------------------------------------------------------------*/
gauge_file *save_serial_archive(char *filename) {
  /* Single node writes in archive file format */

  int currentnode,newnode;
  int i,j,x,y,z,dir;
  su3_matrix lbuf[4];
  gauge_file *gf;
  gauge_header *gh;

  FILE *outfile = NULL;
  site *s;
  u_int32type chksum;
  char sums[30];
  OUTPUT_TYPE *uout = NULL;
  int big_end_p; 
  double ssplaq, stplaq, avgtrace, avgplaq;
  double_complex linktrsum;
  int mu,a,b,vol3=0,tslice;

  /* Check which end is up */
  big_end_p = big_endian();
  /* if(this_node == 0) printf("big_end_p is %d\n", big_end_p); */

  /* Set up gauge file and gauge header structures and load header values */
  gf = setup_output_gauge_file();
  gh = gf->header;
  
  /* Compute plaquette, trace and checksum */
  d_plaquette(&ssplaq, &stplaq);
  avgplaq = (ssplaq+stplaq)/6.0;
  d_linktrsum(&linktrsum);
  avgtrace = linktrsum.real/3.0;
  chksum = nersc_cksum();

  /* node 0 does all the writing */
  if(this_node==0){
    
    printf("trace = %f\n", avgtrace);
    printf("chksum_x = %x\n", chksum);
    printf("chksum_u = %12u\n", chksum);
    printf("plaquette = %f\n", avgplaq);
    
    printf("Writing archive format lattice to %s\n", filename);
    /* Create output file */
    outfile = fopen(filename,"w");
    if (outfile == NULL) {
      printf("error opening output file: %s\n", filename);
      terminate(1);
    }

    fprintf(outfile,"BEGIN_HEADER\n");
    fprintf(outfile,"DATATYPE = 4D_SU3_GAUGE\n");
    fprintf(outfile,"DIMENSION_1 = %d\n",nx);
    fprintf(outfile,"DIMENSION_2 = %d\n",ny);
    fprintf(outfile,"DIMENSION_3 = %d\n",nz);
    fprintf(outfile,"DIMENSION_4 = %d\n",nt);
    fprintf(outfile,"CHECKSUM = %x\n",chksum);
    fprintf(outfile,"LINK_TRACE = %.10f\n",avgtrace);
    fprintf(outfile,"PLAQUETTE = %.10f\n",avgplaq);
    fprintf(outfile,"ENSEMBLE_ID = %s\n", ensemble_id);
    fprintf(outfile,"SEQUENCE_NUMBER = %d\n",sequence_number);
    /* write Milc info section */
    fprintf(outfile,"MILC_INFO = -------BEGIN-------\n");
    write_gauge_info_item(outfile,"time_stamp","\"%s\"",gh->time_stamp,0,0);
    sprintf(sums,"%x %x",gf->check.sum29,gf->check.sum31);
    write_gauge_info_item(outfile,"checksums","\"%s\"",sums,0,0);
    write_gauge_info_item(outfile,"nx","%d",(char *)&nx,0,0);
    write_gauge_info_item(outfile,"ny","%d",(char *)&ny,0,0);
    write_gauge_info_item(outfile,"nz","%d",(char *)&nz,0,0);
    write_gauge_info_item(outfile,"nt","%d",(char *)&nt,0,0);
    write_appl_gauge_info(outfile, gf);
    fprintf(outfile,"MILC_INFO = --------END--------\n");
    fprintf(outfile,"END_HEADER\n");

    vol3 = nx*ny*nz;
    uout = (OUTPUT_TYPE *) malloc(48*vol3*sizeof(OUTPUT_TYPE));
    if(uout == NULL) { 
      printf("can\'t malloc uout timeslice\n"); terminate(1); 
    }
  }

  /* Write gauge field */

  g_sync();
  currentnode=0;

  for(tslice=0; tslice<nt; ++tslice) {
    j = 0;
    for(z=0; z<nz; ++z) for(y=0; y<ny; ++y) for(x=0; x<nx; ++x) {
      newnode=node_number(x,y,z,tslice);
      if(newnode != currentnode){ /* switch to another node */
	/* tell newnode it's OK to send */
	if( this_node==0 && newnode!=0 )send_field((char *)lbuf,32,newnode);
	if( this_node==newnode && newnode!=0 )get_field((char *)lbuf,32,0);
	currentnode=newnode;
      }

      if(this_node==0){
	if(currentnode==0){
	  s = &lattice[node_index(x,y,z,tslice)];
	  for(mu=0; mu<4; ++mu) {
	    for(a=0; a<2; ++a) {
	      for(b=0; b<3; ++b) {
		uout[2*(b+3*a)+12*mu+48*j] 
		    = (OUTPUT_TYPE) s->link[mu].e[a][b].real;
		uout[1+2*(b+3*a)+12*mu+48*j] 
		    = (OUTPUT_TYPE) s->link[mu].e[a][b].imag;
	      }
	    }
	  }
	}
	else{
	  get_field((char *)lbuf,4*sizeof(su3_matrix),currentnode);
	  for(mu=0; mu<4; ++mu) {
	    for(a=0; a<2; ++a) {
	      for(b=0; b<3; ++b) {
		uout[2*(b+3*a)+12*mu+48*j] 
		    = (OUTPUT_TYPE) lbuf[mu].e[a][b].real;
		uout[1+2*(b+3*a)+12*mu+48*j] 
		    = (OUTPUT_TYPE) lbuf[mu].e[a][b].imag;
	      }
	    }
	  }
	}
	++j;
      }
      else {	/* for nodes other than 0 */
	if(this_node==currentnode){
	  i=node_index(x,y,z,tslice);
	  for(dir=XUP;dir<=TUP;dir++)lbuf[dir]=lattice[i].link[dir];
	  send_field((char *)lbuf,4*sizeof(su3_matrix),0);
	}
      }
    }

    if(this_node==0){
      if (!big_end_p) byterevn((int32type *)uout,48*vol3);
      if(g_write(uout,48*vol3*sizeof(OUTPUT_TYPE),1,outfile) != 1)
	printf("g_write bombed...\n");
      fflush(outfile);
    }
  }

  if(this_node==0){
    g_close(outfile);
    printf("Wrote archive gauge file %s\n",filename);
    free(uout);
  }

  g_sync();
  return gf;
}

