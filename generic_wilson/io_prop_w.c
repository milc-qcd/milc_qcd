/*********************** io_prop_w.c *************************/
/* MIMD version 7 */
/* Formerly generic/io_wb.c */
/* This reads/write single precision propagator files */

/* routines for Wilson propagator binary input/output  */
/* The read routines try to accommodate byte reversal */
/* No provision for 64-bit integers here - there should be no such files! */

/* Modifications
   07/11/00 Large file (64-bit addressing) support
   05/10/00 Fixed bit shift for checksums so resuts should be 
            architecture-independent C.D.
   11/20/98 Added w_multidump* and r_multidump* routines for writing 
            and reading separate dump files on each node C.D.
   11/13/98 r_parallel_w: g_syncs to prevent message pileups C.D.
   8/30/96 fixed multiple name conflicts C.D.   
   8/30/96 fixed macros for C syntax UMH   
   08/19/96 Fixed a few misplaced byterev's and corrected some freads and extra closes U.M.H.
   08/12/96 Convert to portable parallel file operations, moving
            system-dependency to com_XXX.c   C.D.
   08/08/96 Upgrade "binary" reads and writes to current format C.D.
   08/05/96 Add parallel reads C.D.
   07/27/96 Add parallel writes C.D.
   07/27/96 Implement new standard for parallel propagator file format. C.D.
*/

/* Includes capability of reading and writing parallel files */

/* The "serial" routines write a "serial" binary file through node 0 alone */
/* The "parallel" routines write a single file from all nodes
   simultaneously in coordinate natural order */
/* The "checkpoint" routines write a single file from all nodes
   simultaneously in node-dump order.  A site list is included.
   This mode is faster if files are to be read to the same
   node number and layout, but it is portable. */
/* The "multidump" routines write a separate file from each node
   simultaneously in node-dump order.  No identifying header is
   provided, but a checksum and color-spin index is recorded to
   provide data checking.  This mode is intended for fast I/O for
   intermediate scratch files on machines that lack a parallel I/O
   feature.  This format is not intended to be portable, so should
   not be used for archive files. */

/* w_serial_w and w_parallel_w writes a file in coordinate 
   natural order with sites ordered as (x,y,z,t) = 
   (0,0,0,0), (1,0,0,0), etc  */

/* w_checkpoint_w writes a file in node_dump order with node 0 coming first */
/* w_multidump_w writes to separate files, one per node */

/* r_serial_w and r_parallel_w read file of all formats */
/* r_parallel_w has all nodes reading simultaneously.  
   r_serial_w - only node 0 */

/* All procedures permit reopening, since many applications can handle
   only one spin and color degree of freedom at a time */

/* w_serial_w_i     Node 0 opens serial file for writing and writes header
   w_serial_w       Node 0 writes propagator to the specified serial file
   w_serial_w_f     Closes the file

   r_serial_w_i     Node 0 Opens serial file for reading and reads header
   r_serial_w       Node 0 reads propagator from specified serial file
   r_serial_w_f     Closes the file

   r_parallel_w_i   Opens binary file for parallel reading and reads header
   r_parallel_w_o   Reopens binary file for parallel reading
   r_parallel_w     All nodes read propagator from binary file
   r_parallel_w_c   Closes the file temporarily
   r_parallel_w_f   Closes the file and frees all structures 

   r_multidump_w_i   Opens private dump file for reading
   r_multidump_w_o   Reopens private dump file for reading
   r_multidump_w     Each node reads dumped propagator from private file
   r_multidump_w_c   Closes the private files temporarily
   r_multidump_w_f   Closes the private files and frees all structures 

   w_parallel_w_i   Opens binary file for parallel writing and writes header
   w_parallel_w_o   Reopens binary file for parallel writing
   w_parallel_w     All nodes write propagator to binary file
   w_parallel_w_c   Closes the file temporarily
   w_parallel_w_f   Closes the file and frees all structures

   w_checkpoint_w_i  Opens binary file for parallel writing and writes header
   w_checkpoint_w_o  Reopens binary file for parallel writing (node-dump order)
   w_checkpoint_w    All nodes write propagator to binary file
   w_checkpoint_w_c  Closes the file temporarily
   w_checkpoint_w_f  Closes the file and frees all structures

   w_multidump_w_i  Opens node-private dump file for writing
   w_multidump_w_o  Reopens node-private dump file for writing
   w_multidump_w    All nodes dump propagator to their private files
   w_multidump_w_c  Closes the dump files temporarily
   w_multidump_w_f  Closes the dump files and frees all structures

   */

/* Requires application-dependent routines:

   build_w_prop_hdr       Fills in table of spins in the header structure
   write_appl_w_prop_info Writes supplemental information to ascii info file

   */


#include "generic_wilson_includes.h"
#include <sys/types.h>
#include <fcntl.h>
#include <errno.h>
#include <time.h>
#include <string.h>
#include <assert.h>
#include "../include/io_lat.h"
#include "../include/io_wprop.h"
#include "../include/file_types.h"

#ifndef HAVE_FSEEKO
#define fseeko fseek
#endif

#define PARALLEL 1
#define SERIAL 0

#define NODE_DUMP_ORDER 1
#define NATURAL_ORDER 0

#define MAXDUMPFILENAME 256

#undef MAX_BUF_LENGTH
#define MAX_BUF_LENGTH 4096

/* The send buffer would normally be the size of a single Wilson
   vector but on we may need to pad for short messages generated
   in serial reads and writes to avoid switch inefficiencies.  In that
   case we define PAD_SEND_BUF on the compilation line to increase
   this.  */
   
#ifndef PAD_SEND_BUF
#define PAD_SEND_BUF 8
#endif

/*---------------------------------------------------------------------------*/
/* Convert a single precision wilson_vector to generic precision */

void f2d_wvec(fwilson_vector *a, wilson_vector *b){
  int c,s;
  
  for(c = 0; c < 3; c++)for(s = 0; s < 4; s++){
    b->d[s].c[c].real = a->d[s].c[c].real;
    b->d[s].c[c].imag = a->d[s].c[c].imag;
  }
}

/* Convert a double precision wilson_vector to single precision */
void d2f_wvec(wilson_vector *a, fwilson_vector *b){
  int c,s;
  
  for(c = 0; c < 3; c++)for(s = 0; s < 4; s++){
    b->d[s].c[c].real = a->d[s].c[c].real;
    b->d[s].c[c].imag = a->d[s].c[c].imag;
  }
}

/*---------------------------------------------------------------------------*/

/* This subroutine writes the propagator header structure */
/* Parallel access version */
/* While the procedures for serial and parallel writing are
   identical, (the header is written only by node 0, no matter what),
   the file which is accessed can be opened either by all
   nodes in w_parallel or one node in w_serial.  We have to
   distinguish between these modes when writing */

void pwrite_w_prop_hdr(FILE *fp, w_prop_header *wph)
{
        /* This would be a slick way to do it, but different
	 compilers align structure members on different
	 word boundaries, so the whole is not the sum of
	 its parts.  For portability, we have to do it piecemeal. */

      /** if( fwrite(&wph,sizeof(wph),1,fp) !=
	 sizeof(wph)) **/

  int i;
  int32type zero32 = 0;
  char myname[] = "pwrite_w_prop_hdr";

  pwrite_data(fp,(void *)&wph->magic_number,sizeof(wph->magic_number),
	      myname,"magic_number");
  pwrite_data(fp,(void *)wph->dims,sizeof(wph->dims),
	      myname,"dimensions");
  pwrite_data(fp,(void *)wph->time_stamp,sizeof(wph->time_stamp),
	      myname,"time_stamp");
  pwrite_data(fp,&wph->order,sizeof(wph->order),
	      myname,"order");
  pwrite_data(fp,&wph->n_spins,sizeof(wph->n_spins),
	      myname,"n_spins");
  for(i=0;i<MAX_SOURCE_SPINS;i++)
    {
      if(i < wph->n_spins)
	pwrite_data(fp,&wph->spins[i],sizeof(wph->spins[i]),
		    myname,"spins");
      else 
	pwrite_data(fp,&zero32,sizeof(zero32),
		    myname,"spins");
    }
    

  /* Header byte length */

  wph->header_bytes = sizeof(wph->magic_number) + sizeof(wph->dims) + 
    sizeof(wph->time_stamp) + sizeof(wph->order) + sizeof(wph->n_spins) +
      sizeof(wph->spins);
  
} /* pwrite_w_prop_hdr */

/*----------------------------------------------------------------------*/
/* This subroutine writes the propagator header structure */
/* Serial access version */

void swrite_w_prop_hdr(FILE *fp, w_prop_header *wph)
{
  int i;
  int32type zero32 = 0;
  char myname[] = "swrite_w_prop_hdr";

  swrite_data(fp,(void *)&wph->magic_number,sizeof(wph->magic_number),
	      myname,"magic_number");
  swrite_data(fp,(void *)wph->dims,sizeof(wph->dims),
	      myname,"dimensions");
  swrite_data(fp,(void *)wph->time_stamp,sizeof(wph->time_stamp),
	      myname,"time_stamp");
  swrite_data(fp,&wph->order,sizeof(wph->order),
	      myname,"order");
  swrite_data(fp,&wph->n_spins,sizeof(wph->n_spins),
	      myname,"n_spins");
  for(i=0;i<MAX_SOURCE_SPINS;i++)
    {
      if(i < wph->n_spins)
	swrite_data(fp,&wph->spins[i],sizeof(wph->spins[i]),
		    myname,"spins");
      else 
	swrite_data(fp,&zero32,sizeof(zero32),
		    myname,"spins");
    }
    

  /* Header byte length */

  wph->header_bytes = sizeof(wph->magic_number) + sizeof(wph->dims) + 
    sizeof(wph->time_stamp) + sizeof(wph->order) + sizeof(wph->n_spins) +
      sizeof(wph->spins);

} /* swrite_w_prop_hdr */

/*------------------------------------------------------------------------*/

/* Write a data item to the propagator info file */
int write_w_prop_info_item( FILE *fpout,    /* ascii file pointer */
		       char *keyword,   /* keyword */
		       char *fmt,       /* output format -
					      must use s, d, e, f, or g */
		       char *src,       /* address of starting data
					   floating point data must be
					   of type (float) */
		       int count,       /* number of data items if > 1 */
		       int stride)      /* byte stride of data if
                                           count > 1 */
{

  int i,k,n;
  char *data;
  Real tt;

  /* Check for valid keyword */

  for(i=0;strlen(w_prop_info_keyword[i])>0 &&
      strcmp(w_prop_info_keyword[i],keyword) != 0; i++);
  if(strlen(w_prop_info_keyword[i])==0)
    printf("write_w_prop_info_item: WARNING: keyword %s not in table\n",
	    keyword);

  /* Write keyword */

  fprintf(fpout,"%s",keyword);

  /* Write count if more than one item */
  if(count > 1)
    fprintf(fpout,"[%d]",count);

  n = count; if(n==0)n = 1;
  
  /* Write data */
  for(k = 0, data = (char *)src; k < n; k++, data += stride)
    {
      fprintf(fpout," ");
      if(strstr(fmt,"s") != NULL)
	fprintf(fpout,fmt,data);
      else if(strstr(fmt,"d") != NULL)
	fprintf(fpout,fmt,*(int *)data);
      else if(strstr(fmt,"e") != NULL || 
	      strstr(fmt,"f") != NULL || 
	      strstr(fmt,"g") != NULL)
	{
	  tt = *(Real *)data;
	  fprintf(fpout,fmt,tt);
	}
      else
	{
	  printf("write_w_prop_info_item: Unrecognized data type %s\n",fmt);
	  return 1;
	}
    }
  fprintf(fpout,"\n");
  return 0;
}

/*------------------------------------------------------------------------*/

/* Write a data item to a character string */
int sprint_w_prop_info_item( 
  char *string,    /* character string */
  size_t nstring,     /* string length */			    
  char *keyword,   /* keyword */
  char *fmt,       /* output format -
		      must use s, d, e, f, or g */
  char *src,       /* address of starting data
		      floating point data must be
		      of type (Real) */
  int count,       /* number of data items if > 1 */
  int stride)      /* byte stride of data if
		      count > 1 */
{

  int i,k,n;
  size_t bytes;
  char *data;
  float tt;

  /* Check for valid keyword */

  for(i=0;strlen(w_prop_info_keyword[i])>0 &&
      strcmp(w_prop_info_keyword[i],keyword) != 0; i++);
  if(strlen(w_prop_info_keyword[i])==0)
    printf("write_w_prop_info_item: WARNING: keyword %s not in table\n",
	    keyword);

  /* Write keyword */
  bytes = 0;

  snprintf(string,nstring-bytes,"%s =",keyword);
  bytes = strlen(string);
  if(bytes >= nstring)return 1;

  /* Write count if more than one item */
  if(count > 1){
    snprintf(string+bytes, nstring-bytes, "[%d]",count);
    bytes = strlen(string);
    if(bytes >= nstring)return 1;
  }
    
  n = count; if(n==0)n = 1;
  
  /* Write data */
  for(k = 0, data = (char *)src; k < n; k++, data += stride)
    {
      snprintf(string+bytes, nstring-bytes," ");
      bytes = strlen(string);
      if(bytes >= nstring)return 1;

      if(strstr(fmt,"s") != NULL){
	snprintf(string+bytes,nstring-bytes, fmt,data);
	bytes = strlen(string);
	if(bytes >= nstring)return 1;
      }
      else if(strstr(fmt,"d") != NULL){
	snprintf(string+bytes,nstring-bytes,fmt,*(int *)data);
	bytes = strlen(string);
	if(bytes >= nstring)return 1;
      }
      else if(strstr(fmt,"e") != NULL || 
	      strstr(fmt,"f") != NULL || 
	      strstr(fmt,"g") != NULL)
	{
	  tt = *(Real *)data;
	  snprintf(string+bytes,nstring-bytes,fmt,tt);
	  bytes = strlen(string);
	  if(bytes >= nstring)return 1;
	}
      else
	{
	  printf("write_w_prop_info_item: Unrecognized data type %s\n",fmt);
	  return 1;
	}
    }
  snprintf(string+bytes,nstring-bytes,"\n");
  bytes = strlen(string);
  if(bytes >= nstring)return 1;

  return 0;
}

/*----------------------------------------------------------------------*/
/* Open, write, and close the ASCII info file */

int write_w_prop_info_file(w_prop_file *wpf)
{
  FILE *info_fp;
  w_prop_header *wph;
  char info_filename[256];
  void write_appl_w_prop_info(FILE *fp);

  wph = wpf->header;

  /* Construct header file name from lattice file name 
   by adding filename extension to lattice file name */

  strcpy(info_filename,wpf->filename);
  strcat(info_filename,ASCII_W_PROP_INFO_EXT);

  /* Open header file */
  
  if((info_fp = fopen(info_filename,"w")) == NULL)
    {
      printf("write_w_prop_info_file: Can't open ascii info file %s\n",info_filename);
      return 1;
    }
  
  /* Write required information */

  write_w_prop_info_item(info_fp,"magic_number","%d",
			 (char *)&wph->magic_number,0,0);
  write_w_prop_info_item(info_fp,"time_stamp","\"%s\"",wph->time_stamp,0,0);
  write_w_prop_info_item(info_fp,"nx","%d",(char *)&nx,0,0);
  write_w_prop_info_item(info_fp,"ny","%d",(char *)&ny,0,0);
  write_w_prop_info_item(info_fp,"nz","%d",(char *)&nz,0,0);
  write_w_prop_info_item(info_fp,"nt","%d",(char *)&nt,0,0);

  /* Write optional information - this routine supplied by application */
  write_appl_w_prop_info(info_fp);

  fclose(info_fp);

/*  printf("Wrote propagator info file %s\n",info_filename); */

  return 0;
} /*write_w_prop_info_file */

/*----------------------------------------------------------------------*/

/* Set up the input propagator file structure and header structure */

w_prop_file *setup_input_w_prop_file(char *filename)
{
  w_prop_file *wpf;
  w_prop_header *wph;

  /* Allocate space for the file structure */

  wpf = (w_prop_file *)malloc(sizeof(w_prop_file));
  if(wpf == NULL)
    {
      printf("setup_input_w_prop_file: Can't malloc wpf\n");
      terminate(1);
    }

  wpf->filename = filename;

  /* Allocate space for the header */

  /* Make sure compilation gave us a 32 bit integer type */
  assert(sizeof(int32type) == 4);

  wph = (w_prop_header *)malloc(sizeof(w_prop_header));
  if(wph == NULL)
    {
      printf("setup_input_w_prop_file: Can't malloc wph\n");
      terminate(1);
    }

  wpf->header = wph;

  return wpf;
} /* setup_input_w_prop_file */

/*----------------------------------------------------------------------*/

/* Set up the output w_prop file and w_prop header structure */

void build_w_prop_hdr(w_prop_header *);

w_prop_file *setup_output_w_prop_file()
{
  w_prop_file *wpf;
  w_prop_header *wph;
  time_t time_stamp;
  int i;

  /* Allocate space for a new file structure */

  wpf = (w_prop_file *)malloc(sizeof(w_prop_file));
  if(wpf == NULL)
    {
      printf("setup_w_prop_header: Can't malloc wpf\n");
      terminate(1);
    }

  /* Allocate space for a new header structure */

  /* Make sure compilation gave us a 32 bit integer type */
  assert(sizeof(int32type) == 4);

  wph = (w_prop_header *)malloc(sizeof(w_prop_header));
  if(wph == NULL)
    {
      printf("setup_w_prop_header: Can't malloc wph\n");
      terminate(1);
    }

  /* Load header pointer and file name */
  wpf->header = wph;

  /* Load header values */

  wph->magic_number = W_PROP_VERSION_NUMBER;

  wph->dims[0] = nx;
  wph->dims[1] = ny;
  wph->dims[2] = nz;
  wph->dims[3] = nt;

  /* Get date and time stamp. (We use local time on node 0) */

  if(this_node==0)
    {
      time(&time_stamp);
      strcpy(wph->time_stamp,ctime(&time_stamp));

      /* For aesthetic reasons, don't leave trailing junk bytes here to be
	 written to the file */
      for(i = strlen(wph->time_stamp) + 1; 
	  i < (int)sizeof(wph->time_stamp); i++)
	wph->time_stamp[i] = '\0';
      
      /* Remove trailing end-of-line character */
      if(wph->time_stamp[strlen(wph->time_stamp) - 1] == '\n')
	wph->time_stamp[strlen(wph->time_stamp) - 1] = '\0';
    }

  /* Broadcast to all nodes */
  broadcast_bytes(wph->time_stamp,sizeof(wph->time_stamp));

  /* Add spin index information to header structure */
  
  build_w_prop_hdr(wph);
  
  return wpf;
} /* setup_output_w_prop_file */

/*----------------------------------------------------------------------*/

/* Open a binary file for serial writing by node 0 */

w_prop_file *w_serial_w_i(char *filename)
{
  /* Only node 0 opens the file filename */
  /* Returns a file structure describing the opened file */

  FILE *fp;
  w_prop_file *wpf;
  w_prop_header *wph;

  /* Set up w_prop file and w_prop header structures and load header values */
  wpf = setup_output_w_prop_file();
  wph = wpf->header;

  /* Indicate coordinate natural ordering */

  wph->order = NATURAL_ORDER;

  /* Only node 0 opens the requested file */

  if(this_node == 0)
    {
      fp = fopen(filename, "wb");
      if(fp == NULL)
	{
	  printf("w_serial_w_i: Node %d can't open file %s, error %d\n",
		 this_node,filename,errno);fflush(stdout);
	  terminate(1);
	}

/*      printf("Opened prop file %s for serial writing\n",filename); */
      
      /* Node 0 writes the header */
      
      swrite_w_prop_hdr(fp,wph);

    }
  
  /* Assign values to file structure */

  if(this_node==0)wpf->fp = fp; 
  else wpf->fp = NULL;                /* Only node 0 knows about this file */

  wpf->filename       = filename;
  wpf->byterevflag    = 0;            /* Not used for writing */
  wpf->rank2rcv       = NULL;         /* Not used for writing */
  wpf->parallel       = SERIAL;

  /* Node 0 writes ascii info file */

  if(this_node == 0)write_w_prop_info_file(wpf);

  return wpf;

} /* w_serial_w_i */

/*---------------------------------------------------------------------------*/
/* Write checksum to lattice file.  It is assumed that the file
   is already correctly positioned.

   Should be called only by one node */

void write_checksum_w(int parallel, w_prop_file *wpf)
{

  char myname[] = "write_checksum";

  pswrite_data(parallel,wpf->fp,
	       &wpf->check.spin,sizeof(wpf->check.spin),myname,"checksum");
  pswrite_data(parallel,wpf->fp,
	       &wpf->check.color,sizeof(wpf->check.color),myname,"checksum");
  pswrite_data(parallel,wpf->fp,
	       &wpf->check.sum29,sizeof(wpf->check.sum29),myname,"checksum");
  pswrite_data(parallel,wpf->fp,
	       &wpf->check.sum31,sizeof(wpf->check.sum31),myname,"checksum");
  /* printf("Prop checksums %x %x for spin %d color %d\n",
	 wpf->check.sum29,wpf->check.sum31,wpf->check.spin,wpf->check.color);*/
} /* write_checksum_w */

/*---------------------------------------------------------------------------*/
/* Here only node 0 writes propagator to a serial file for this spin and color */

void w_serial_w(w_prop_file *wpf, int spin, int color, field_offset src_site,
		wilson_vector *src_field)
{
  /* wpf  = file descriptor as opened by w_serial_w_i 
     src  = field offset for propagator Wilson vector (type wilson_vector)  */

  FILE *fp = NULL;
  w_prop_header *wph;
  u_int32type *val;
  int rank29,rank31;
  fwilson_vector *lbuf = NULL;
  wilson_vector *src;
  int fseek_return;  /* added by S.G. for large file debugging */
  struct {
    fwilson_vector wv;
    char pad[PAD_SEND_BUF];    /* Introduced because some switches
				  perform better if message lengths are longer */
  } msg;
  int buf_length;
  register int i,j,k;
  off_t offset;             /* File stream pointer */
  off_t w_prop_size;        /* Size of propagator blocks for all nodes */
  off_t w_prop_check_size;  /* Size of propagator checksum record */
  off_t coord_list_size;    /* Size of coordinate list in bytes */
  off_t head_size = 0;      /* Size of header plus coordinate list */
  off_t body_size = 0;      /* Size of propagator blocks for all nodes 
			      plus checksum record */
  int currentnode,newnode;
  int x,y,z,t;
  int spinindex = 0;

  if(this_node==0)
    {
      if(wpf->parallel == PARALLEL)
	printf("w_serial_w: Attempting serial write to file opened in parallel \n");

      lbuf = (fwilson_vector *)malloc(MAX_BUF_LENGTH*sizeof(fwilson_vector));
      if(lbuf == NULL)
	{
	  printf("w_serial_w: Node 0 can't malloc lbuf\n"); 
	  fflush(stdout);terminate(1);
        }

      fp = wpf->fp;
      wph = wpf->header;
      
      w_prop_size = volume*sizeof(fwilson_vector) ;
      w_prop_check_size =  sizeof(wpf->check.spin) +
	sizeof(wpf->check.color) +  sizeof(wpf->check.sum29) + 
	  sizeof(wpf->check.sum31);
      body_size = w_prop_size + w_prop_check_size;
      
      /* No coordinate list was written because fields are to be written
	 in standard coordinate list order */
      
      coord_list_size = 0;
      head_size = wph->header_bytes + coord_list_size;
      
      /* Find spin in table of contents */
      
      for(spinindex=0 ; spinindex < wph->n_spins ; spinindex++)
	if(wph->spins[spinindex]==spin)break;
      
      if(spinindex == wph->n_spins)
	{
	  printf("w_serial_w: Requested spin %d not planned for file %s\n",spin,wpf->filename);
	  printf("  Table of contents: ");
	  for(spinindex=0; spinindex<wph->n_spins; spinindex++)
	    printf(" %d",wph->spins[spinindex]);
	  printf("\n");fflush(stdout);
	  terminate(1);
	}
      
      offset = head_size + body_size*(spinindex*3 + color)
           + w_prop_check_size;

      fseek_return=fseeko(fp,offset,SEEK_SET);
      /* printf("w_serial_w: Node %d fseek_return = %d\n",this_node,fseek_return); */
      if( fseek_return < 0 ) 
	{
	  printf("w_serial_w: Node %d fseeko %lld failed error %d file %s\n",
		 this_node,(long long)offset,errno,wpf->filename);
	  fflush(stdout);terminate(1);
	}
    }
      
  /* Buffered algorithm for writing fields in serial order */
  
  /* initialize checksums */
  wpf->check.sum31 = 0;
  wpf->check.sum29 = 0;
  /* counts 32-bit words mod 29 and mod 31 in order of appearance on file */
  /* Here only node 0 uses these values */
  rank29 = sizeof(fwilson_vector)/sizeof(int32type)*sites_on_node*this_node % 29;
  rank31 = sizeof(fwilson_vector)/sizeof(int32type)*sites_on_node*this_node % 31;

  g_sync();
  currentnode=0;

  buf_length = 0;

  for(j=0,t=0;t<nt;t++)for(z=0;z<nz;z++)for(y=0;y<ny;y++)for(x=0;x<nx;x++,j++)
    {
      newnode=node_number(x,y,z,t);
      if(newnode != currentnode){	/* switch to another node */
	/* Send a few bytes of garbage to tell newnode it's OK to send */
	if( this_node==0 && newnode!=0 )
	  send_field((char *)&msg,sizeof(msg),newnode);
	if( this_node==newnode && newnode!=0 )
	  get_field((char *)&msg,sizeof(msg),0);
	currentnode=newnode;
      }
      
      if(this_node==0)
	{
	  if(currentnode==0)
	    {
	      i=node_index(x,y,z,t);
	      /* Load msg structure with single precision wilson vector */

 	      if(src_site == (field_offset)(-1))
		src = src_field + i;
	      else
		src = (wilson_vector *)F_PT( &(lattice[i]), src_site );

	      d2f_wvec(src, &msg.wv);
	    }
	  else
	    {
	      get_field((char *)&msg, sizeof(msg),currentnode);
	    }

	  lbuf[buf_length] = msg.wv;

	  /* Accumulate checksums - contribution from next site */
	  for(k = 0, val = (u_int32type *)&lbuf[buf_length]; 
	      k < (int)sizeof(fwilson_vector)/(int)sizeof(int32type); k++, val++)
	    {
	      wpf->check.sum29 ^= (*val)<<rank29 | (*val)>>(32-rank29);
	      wpf->check.sum31 ^= (*val)<<rank31 | (*val)>>(32-rank31);
	      rank29++; if(rank29 >= 29)rank29 = 0;
	      rank31++; if(rank31 >= 31)rank31 = 0;
	    }

	  buf_length++;
	  
	  if( (buf_length == MAX_BUF_LENGTH) || (j == volume-1))
	    {
	      /* write out buffer */
	      
	      if( (int)fwrite(lbuf,sizeof(fwilson_vector),buf_length,fp) 
		  != buf_length)
		{
		  printf("w_serial_w: Node %d propagator write error %d file %s\n",
			 this_node,errno,wpf->filename); 
		  fflush(stdout);
		  terminate(1);   
		}
	      buf_length = 0;		/* start again after write */
	    }
	}
      else  /* for nodes other than 0 */
	{	
	  if(this_node==currentnode){
	    i=node_index(x,y,z,t);
	    /* Copy or convert data into send buffer and send to node
	       0 with padding */

	    if(src_site == (field_offset)(-1))
	      src = src_field + i;
	    else
	      src = (wilson_vector *)F_PT( &(lattice[i]), src_site );

	    d2f_wvec(src, &msg.wv);
	    send_field((char *)&msg, sizeof(msg),0);
	  }
	}
      
    } /*close x,y,z,t loops */
  
  g_sync();
  
  if(this_node==0)
    {
/*      printf("Wrote prop serially for spin %d color %d to file %s\n",
	     spin,color,wpf->filename);fflush(stdout); */
      free(lbuf);

      /* Construct check record */

      /* Convert to 32 bit integer */
      wpf->check.spin = spin;
      wpf->check.color = color;
      
      /* Position file pointer for writing check record */
      /* Record PRECEDES the propagator data for a given spin and color */
      
      offset = head_size + body_size*(spinindex*3 + color);
      if( fseeko(fp,offset,SEEK_SET) < 0 ) 
	{
	  printf("w_serial_w: Node %d fseeko %lld failed error %d file %s\n",
		 this_node,(long long)offset,errno,wpf->filename);
	  fflush(stdout);terminate(1);   
	}
      
      write_checksum_w(SERIAL,wpf);
    }
      
} /* w_serial_w */

/*---------------------------------------------------------------------------*/
/* Here only node 0 writes propagator to a serial file for this spin and color */

void w_serial_w_from_site(w_prop_file *wpf, int spin, int color, 
			  field_offset src_site)
{
  /* wpf  = file descriptor as opened by w_serial_w_i 
     src  = field offset for propagator Wilson vector (type wilson_vector)  */

  w_serial_w(wpf, spin, color, src_site, NULL);
}

/*---------------------------------------------------------------------------*/
/* Here only node 0 writes propagator to a serial file for this spin and color */

void w_serial_w_from_field(w_prop_file *wpf, int spin, int color, 
		wilson_vector *src_field)
{
  /* wpf  = file descriptor as opened by w_serial_w_i 
     src  = field offset for propagator Wilson vector (type wilson_vector)  */

  w_serial_w(wpf, spin, color, (field_offset)(-1), src_field);
}

/*---------------------------------------------------------------------------*/

void w_serial_w_f(w_prop_file *wpf)

/* this subroutine closes the file and frees associated structures */
{
  g_sync();
  if(this_node==0)
    {
      if(wpf->parallel == PARALLEL)
	printf("w_serial_w_f: Attempting serial close on file opened in parallel \n");

      printf("Wrote prop file %s time stamp %s\n",wpf->filename,
	     (wpf->header)->time_stamp);

      if(wpf->fp != NULL)fclose(wpf->fp);
    }

  /* Free header and file structures */
  free(wpf->header);
  free(wpf);

} /* w_serial_w_f */

/*---------------------------------------------------------------------------*/

/* Subroutine for reading site list from propagator file */
/* Only node 0 reads this list, so same for parallel and serial reading */

void read_site_list_w(int parallel, w_prop_file *wpf)
{

  /* All nodes allocate space for site list table, 
     but not needed if file is in natural order */

  if(wpf->header->order != NATURAL_ORDER)
    {
      wpf->rank2rcv = (int32type *)malloc(volume*sizeof(int32type));
      if(wpf->rank2rcv == NULL)
	{
	  printf("read_site_list_w: Can't malloc rank2rcv table\n");
	  terminate(1);
	}

      /* Only node 0 reads the site list */
      
      if(this_node==0)
	{
	  
	  /* Reads receiving site coordinate if file is not in natural order */
          if(parallel)
            {
	      if((int)g_read(wpf->rank2rcv,sizeof(int32type),volume,wpf->fp) 
		 != volume )
		{
		  printf("read_site_list_w: Node %d site list read error %d\n",
			 this_node,errno);
		  terminate(1);	
		}
            }
	  else
 	    {
	      if((int)fread(wpf->rank2rcv,sizeof(int32type),volume,wpf->fp) 
		 != volume )
		{
		  printf("read_site_list_w: Node %d site list read error %d\n",
			 this_node,errno);
		  terminate(1);	
		}
	    }
 	  if(wpf->byterevflag==1)byterevn(wpf->rank2rcv,volume);
	}

      /* Broadcast result to all nodes */

      broadcast_bytes((char *)wpf->rank2rcv,volume*sizeof(int32type));
    }
      
  else wpf->rank2rcv = NULL;

} /* read_site_list_w */

/*----------------------------------------------------------------------*/

int read_1996_w_prop_hdr(w_prop_file *wpf, int parallel, int *byterevflag)
{
  /* parallel = 1 (true) if all nodes are accessing the file */
  /*            0 for access from node 0 only */

  /* Return value  0 if successful
		   1 if not 1996 format */

  FILE *fp;
  w_prop_header *wph;
  int32type tmp;
  int j;
  w_prop_header_1996 wph96;
  char myname[] = "read_1996_w_prop_hdr";
  
  fp = wpf->fp;
  wph = wpf->header;
  
  /* Assumes the magic number has already been read */
  
  tmp = wph->magic_number;
  
  if(wph->magic_number == W_PROP_VERSION_NUMBER_1996) 
    {
      printf("Reading as 1996-style propagator field configuration.\n");
      *byterevflag=0;
    }
  else 
    {
      byterevn((int32type *)&wph->magic_number,1);
      if(wph->magic_number == W_PROP_VERSION_NUMBER_1996) 
	{
	  *byterevflag=1;
	  printf("Reading as 1996-style propagator file with byte reversal\n");
	  if( sizeof(float) != sizeof(int32type)) {
	    printf("%s: Can't byte reverse\n",myname);
	    printf("requires size of int32type(%d) = size of float(%d)\n",
		   (int)sizeof(int32type),(int)sizeof(float));
	    terminate(1);
	  }
	}
      /* Not recognized as 1996 format */
      else
	{
	  /* Not recognized as 1996 format */
	  /* Restore header to entry state */
	  wph->magic_number = tmp;
	  return 1;  /* error signal */
	}
    }
  
  /* So we have a 1996-format propagator file
     Read header into wph96 structure, do byte reversal, 
     if necessary, and check consistency */
  
  /* Lattice dimensions */
  
  if(psread_byteorder(*byterevflag,parallel,fp,wph96.dims,sizeof(wph96.dims),
		   myname,"dimensions")!=0)terminate(1);
  
  if(wph96.dims[0] != nx || 
     wph96.dims[1] != ny ||
     wph96.dims[2] != nz ||
     wph96.dims[3] != nt)
    {
      /* So we can use this routine to discover the dimensions,
	 we provide that if nx = ny = nz = nt = -1 initially
	 we don't die */
      if(nx != -1 || ny != -1 || nz != -1 || nt != -1)
	{
	  printf("%s: Incorrect lattice dimensions ",myname);
	  for(j=0;j<4;j++)
	    printf("%d ",wph96.dims[j]); 
	  printf("\n");fflush(stdout);terminate(1);
	}
      else
	{
	  nx = wph->dims[0];
	  ny = wph->dims[1];
	  nz = wph->dims[2];
	  nt = wph->dims[3];
	  volume = nx*ny*nz*nt;
	}
    }
  
  /* Header byte length */
  
  if(psread_byteorder(*byterevflag,parallel,fp,
		   &wph96.header_bytes,sizeof(wph96.header_bytes),
		   myname,"header size")!=0)terminate(1);
  
  /* Data order */
  
  if(psread_byteorder(*byterevflag,parallel,fp,
		   &wph96.order,sizeof(wph96.order),
		   myname,"order")!=0)terminate(1);
  
  /* Length of gauge field descriptor */
  
  if(psread_byteorder(*byterevflag,parallel,fp,
		   &wph96.gauge_field.n_descript,
		   sizeof(wph96.gauge_field.n_descript),
		   myname,"n_descript")!=0)terminate(1);
  
  if(wph96.gauge_field.n_descript > MAX_GAUGE_FIELD_DESCRIPT)
    {
      printf("%s: gauge field descriptor length %d\n",myname,
	     wph96.gauge_field.n_descript);
      printf(" exceeds allocated space %d\n",
	     MAX_GAUGE_FIELD_DESCRIPT);
      terminate(1);
    }
  
  /* Gauge field descriptor */
  
  /* We read the specified length, rather than the allocated length */
  /* Read without bytereversal */
  
  if(psread_data(parallel,fp,wph96.gauge_field.descript,
	      sizeof(wph96.gauge_field.descript),
	      myname,"descrip")!=0)terminate(1);
  
  /* Assures termination of string */
  wph96.gauge_field.descript
    [wph96.gauge_field.n_descript-1] = '\0';
  
  /* Number of gauge field parameters */
  
  if(psread_byteorder(*byterevflag,parallel,fp,
		   &wph96.gauge_field.n_param,sizeof(wph96.gauge_field.n_param),
		   myname,"n_param")!=0)terminate(1);
  
  if(wph96.gauge_field.n_param > MAX_GAUGE_FIELD_PARAM )
    {
      printf("%s: gauge field parameter vector length %d\n",myname,
	     wph96.gauge_field.n_param);
      printf(" exceeds allocated space %d\n",
	     MAX_GAUGE_FIELD_PARAM);
      terminate(1);
    }
  
  /* Gauge field parameters */
  
  for(j=0;j<wph96.gauge_field.n_param;j++)
    if(psread_byteorder(*byterevflag,parallel,fp,
		     &wph96.gauge_field.param[j],
		     sizeof(wph96.gauge_field.param[j]),
		     myname,"gauge param")!=0)terminate(1);
  
  /* Length of Dirac matrix descriptor */
  
  if(psread_byteorder(*byterevflag,parallel,fp,
		   &wph96.dirac.n_descript,sizeof(wph96.dirac.n_descript),
		   myname,"Dirac n_descript")!=0)terminate(1);
  
  if(wph96.dirac.n_descript > MAX_DIRAC_DESCRIPT)
    {
      printf("%s: Dirac matrix descriptor length %d\n",myname,
	     wph96.dirac.n_descript);
      printf(" exceeds allocated space %d\n",
	     MAX_DIRAC_DESCRIPT);
      terminate(1);
    }
  
  /* Dirac matrix descriptor */
  
  /* Read without bytereversal */

  if(psread_data(parallel,fp,
	      wph96.dirac.descript,sizeof(wph96.dirac.descript),
	      myname,"Dirac descript")!=0)terminate(1);
  
  /* Assures termination of string */
  wph96.dirac.descript
    [wph96.dirac.n_descript-1] = '\0';
  
  /* Number of Dirac matrix parameters */
  
  if(psread_byteorder(*byterevflag,parallel,fp,
		     &wph96.dirac.n_param,sizeof(wph96.dirac.n_param),
		     myname,"Dirac n_param")!=0)terminate(1);

  if(wph96.dirac.n_param > MAX_DIRAC_PARAM )
    {
      printf("%s: Dirac matrix parameter vector length %d\n",myname,
	     wph96.dirac.n_param);
      printf(" exceeds allocated space %d\n",
	     MAX_DIRAC_PARAM);
      terminate(1);
    }
  /* Dirac matrix parameters */
  
  for(j=0;j<wph96.dirac.n_param;j++)
    if(psread_byteorder(*byterevflag,parallel,fp,
		     &wph96.dirac.param[j],sizeof(wph96.dirac.param[j]),
		     myname,"Dirac param")!=0)terminate(1);
  
  /* Length of source descriptor */
  
  if(psread_byteorder(*byterevflag,parallel,fp,
		   &wph96.source.n_descript,sizeof(wph96.source.n_descript),
		   myname,"source n_descript")!=0)terminate(1);

  if(wph96.source.n_descript > MAX_SOURCE_DESCRIPT)
    {
      printf("%s: Source descriptor length %d\n",myname,
	     wph96.source.n_descript);
      printf(" exceeds allocated space %d\n",
	     MAX_SOURCE_DESCRIPT);
      terminate(1);
    }
  
  /* Source descriptor */
  
  /* We read the specified length, rather than the allocated length */
  /* Read without bytereversal */
  
  if(psread_data(parallel,fp,
	      wph96.source.descript,sizeof(wph96.source.descript),
	      myname,"source descript")!=0)terminate(1);
  
  /* Assures termination of string */
  wph96.source.descript
    [wph96.source.n_descript-1] = '\0';
  
  /* Number of source parameters */
  
  if(psread_byteorder(*byterevflag,parallel,fp,
		   &wph96.source.n_param,sizeof(wph96.source.n_param),
		   myname,"source n_param")!=0)terminate(1);
  
  if(wph96.source.n_param > MAX_SOURCE_PARAM )
    {
      printf("%s: Source parameter vector length %d\n",myname,
	     wph96.source.n_param);
      printf(" greater than allocated length %d\n",
	     MAX_SOURCE_PARAM);
      terminate(1);
    }
  
  
  /* Source parameters */
  
  if(psread_byteorder(*byterevflag,parallel,fp,
		   &wph96.source.param.i1,sizeof(wph96.source.param.i1),
		   myname,"source param i1")!=0)terminate(1);
  if(psread_byteorder(*byterevflag,parallel,fp,
		   &wph96.source.param.c1,sizeof(wph96.source.param.c1),
		   myname,"source param c1")!=0)terminate(1);
  
  /* Number of source spin values recorded */
  
  if(psread_byteorder(*byterevflag,parallel,fp,
		   &wph96.source.n_spins,sizeof(wph96.source.n_spins),
		   myname,"source n_spins")!=0)terminate(1);
  
  if(wph96.source.n_spins > MAX_SOURCE_SPINS)
    {
      printf("%s: Source spins vector length %d\n",myname,
	     wph96.source.n_spins);
      printf(" exceeds allocated space %d\n",
	     MAX_SOURCE_SPINS);
      terminate(1);
    }
  
  
  /* List of source spin values recorded */
  
  for(j=0;j<wph96.source.n_spins;j++)
    {
      if(psread_byteorder(*byterevflag,parallel,fp,
		       &wph96.source.spins[j],sizeof(wph96.source.spins[j]),
		       myname,"source spins")!=0)terminate(1);
      if(wph96.source.spins[j] > 
	 wph96.source.n_spins)
	{
	  printf("%s: spin value %d exceeds %d\n",myname,
		 wph96.source.spins[j],
		 wph96.source.n_spins);
	  terminate(1);
	}
    }

  /* Load relevant quantities into working header structure */

  for(j=0;j<4;j++)wph->dims[j] = wph96.dims[j];
  wph->header_bytes = wph96.header_bytes;
  wph->order = wph96.order;
  wph->n_spins = wph96.source.n_spins;
  for(j=0;j<wph96.source.n_spins;j++)
    wph->spins[j] = wph96.source.spins[j];
  
  return 0;
  
} /* read_1996_w_prop_hdr */

/*----------------------------------------------------------------------*/

int read_w_prop_hdr(w_prop_file *wpf, int parallel)
{
  /* parallel = 1 (TRUE) if all nodes are accessing the file */
  /*            0        for access from node 0 only */

  /* Returns byterevflag  = 0 or 1 */

  FILE *fp;
  w_prop_header *wph;
  int32type tmp;
  int j;
  int byterevflag;
  char myname[] = "read_w_prop_hdr";

  fp = wpf->fp;
  wph = wpf->header;

  /* Read and verify magic number */

  if(psread_data(parallel, fp,&wph->magic_number,sizeof(wph->magic_number),
	      myname,"magic number")!=0)terminate(1);

  tmp = wph->magic_number;
  
  if(wph->magic_number == W_PROP_VERSION_NUMBER) 
    byterevflag=0;
  else 
    {
      byterevn((int32type *)&wph->magic_number,1);
      if(wph->magic_number == W_PROP_VERSION_NUMBER) 
	{
	  byterevflag=1;
	  /** printf("Reading with byte reversal\n"); **/
	  if( sizeof(float) != sizeof(int32type)) {
	    printf("%s: Can't byte reverse\n",myname);
	    printf("requires size of int32type(%d) = size of float(%d)\n",
		   (int)sizeof(int32type),(int)sizeof(float));
	    terminate(1);
	  }
	}
      else
	{
	  /* Restore magic number as originally read */
	  wph->magic_number = tmp;
	  
	  /* Try old-style configurations */
	  if(read_1996_w_prop_hdr(wpf,parallel,&byterevflag) != 0)
	    {
	      /* End of the road. */
	      printf("%s: Unrecognized magic number in prop file header.\n",
		     myname);
	      printf("Expected %x but read %x\n",
		     W_PROP_VERSION_NUMBER,tmp);
	      terminate(1);
	    }
	  return byterevflag;
	}
    }
  
  /* Read header, do byte reversal, 
     if necessary, and check consistency */
  
  /* Lattice dimensions */
  
  if(psread_byteorder(byterevflag,parallel,fp,wph->dims,sizeof(wph->dims),
		   myname,"dimensions")!=0)terminate(1);

  if(wph->dims[0] != nx || 
     wph->dims[1] != ny ||
     wph->dims[2] != nz ||
     wph->dims[3] != nt)
    {
      /* So we can use this routine to discover the dimensions,
	 we provide that if nx = ny = nz = nt = -1 initially
	 we don't die */
      if(nx != -1 || ny != -1 || nz != -1 || nt != -1)
	{
	  printf("%s: Incorrect lattice dimensions ",myname);
	  for(j=0;j<4;j++)
	    printf("%d ",wph->dims[j]); 
	  printf("\n");fflush(stdout);terminate(1);
	}
      else
	{
	  nx = wph->dims[0];
	  ny = wph->dims[1];
	  nz = wph->dims[2];
	  nt = wph->dims[3];
	  volume = nx*ny*nz*nt;
	}
    }
  
  /* Date and time stamp */

  if(psread_data(parallel,fp,wph->time_stamp,sizeof(wph->time_stamp),
	      myname,"time stamp")!=0)terminate(1);

  /* Header byte length */

  wph->header_bytes = sizeof(wph->magic_number) + sizeof(wph->dims) + 
    sizeof(wph->time_stamp) + sizeof(wph->order) + sizeof(wph->n_spins) +
      sizeof(wph->spins);
  
  /* Data order */
  
  if(psread_byteorder(byterevflag,parallel,fp,&wph->order,sizeof(wph->order),
		   myname,"order parameter")!=0)terminate(1);
  
  /* Number of source spins */
  
  if(psread_byteorder(byterevflag,parallel,fp,&wph->n_spins,sizeof(wph->n_spins),
		   myname,"n_spins")!=0)terminate(1);
  
  /* List of spin values */

  for(j=0;j<wph->n_spins;j++)
    if(psread_byteorder(byterevflag,parallel,fp,&wph->spins[j],sizeof(wph->spins[j]),
		     myname,"spins[j")!=0)terminate(1);
  
  return byterevflag;
  
} /* read_w_prop_hdr */

/*---------------------------------------------------------------------------*/

w_prop_file *r_serial_w_i(char *filename)
{
  /* Returns file descriptor for opened file */

  w_prop_header *wph;
  w_prop_file *wpf;
  FILE *fp;
  int byterevflag;

  /* All nodes set up a propagator file and propagator header
     structure for reading */

  wpf = setup_input_w_prop_file(filename);
  wph = wpf->header;

  /* File opened for serial reading */
  wpf->parallel = 0;

  /* Node 0 alone opens a file and reads the header */

  if(this_node==0)
    {
      fp = fopen(filename, "rb");
      if(fp == NULL)
	{
	  printf("r_serial_w_i: Node %d can't open file %s, error %d\n",
		 this_node,filename,errno);fflush(stdout);terminate(1);
	}
      
/*      printf("Opened prop file %s for serial reading\n",filename); */
      
      wpf->fp = fp;

      byterevflag = read_w_prop_hdr(wpf,SERIAL);

    }

  else wpf->fp = NULL;

  /* Broadcast the byterevflag from node 0 to all nodes */
      
  broadcast_bytes((char *)&byterevflag,sizeof(byterevflag));
  wpf->byterevflag = byterevflag;
  
  /* Node 0 broadcasts the header structure to all nodes */
  
  broadcast_bytes((char *)wph,sizeof(w_prop_header));

  /* Node 0 reads site list and assigns wpf->rank2rcv */

  read_site_list_w(SERIAL,wpf);

  return wpf;

}/* r_serial_w_i */

/*----------------------------------------------------------------------*/

/* Here only node 0 reads the Wilson propagator from a binary file */

int r_serial_w(w_prop_file *wpf, int spin, int color, field_offset dest_site,
	       wilson_vector *dest_field)
{
  /* wpf  = propagator file structure 
     dest_site  = field offset for propagator Wilson vector
     dest_field = pointer to field
     use only one of the dest parameters!  */

  /* 0 is normal exit code
     1 for seek, read error, or missing data error */

  FILE *fp;
  w_prop_header *wph;
  char *filename;
  int byterevflag;
  int spinindex = 0;        /* Counts spin records in file   -
			       wph->spins[spinindex] = spin */

  off_t offset ;            /* File stream pointer */
  off_t w_prop_size;        /* Size of propagator blocks for all nodes */
  off_t w_prop_check_size;  /* Size of propagator checksum record */
  off_t coord_list_size;    /* Size of coordinate list in bytes */
  off_t head_size;          /* Size of header plus coordinate list */
  off_t body_size = 0 ;     /* Size of propagator blocks for all nodes 
			      plus checksum record */
  int rcv_rank, rcv_coords;
  int destnode;
  int i,k,x,y,z,t;
  int status;
  int buf_length = 0,where_in_buf = 0;
  w_prop_check test_wpc;
  u_int32type *val;
  int rank29,rank31;
  fwilson_vector *lbuf = NULL;
  wilson_vector *dest = NULL;

  struct {
    fwilson_vector wv;
    char pad[PAD_SEND_BUF];    /* Introduced because some switches
				  perform better if message lengths are longer */
  } msg;

  char myname[] = "r_serial_w";

  fp = wpf->fp;
  wph = wpf->header;
  filename = wpf->filename;
  byterevflag = wpf->byterevflag;

  status = 0;
  if(this_node == 0)
    {

      if(wpf->parallel == PARALLEL)
	printf("%s: Attempting serial read from parallel file \n",myname);
      
      w_prop_size = volume*sizeof(fwilson_vector) ;
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
	  status = 1;
	}
    }
  broadcast_bytes((char *)&status,sizeof(int));
  if(status != 0)return status;

  status = 0;
  if(this_node == 0)
    {
      if(wph->order == NATURAL_ORDER)coord_list_size = 0;
      else coord_list_size = sizeof(int32type)*volume;
      head_size = wph->header_bytes + coord_list_size;
      
      offset = head_size + body_size*(spinindex*3 + color);
      
      lbuf = (fwilson_vector *)malloc(MAX_BUF_LENGTH*sizeof(fwilson_vector));
      if(lbuf == NULL)
	{
	  printf("%s: Node %d can't malloc lbuf\n",myname,this_node);
	  fflush(stdout);
	  terminate(1);
	}
      
      /* Position file pointer for reading check record */
      
      if( fseeko(fp,offset,SEEK_SET) < 0 ) 
	{
	  printf("%s: Node %d fseeko %lld failed error %d file %s\n",
		 myname,this_node,(long long)offset,errno,filename);
	  fflush(stdout);
	  status = 1;
	}
      
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
	  status = 1;
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
	    
	    if( (int)fread(lbuf,sizeof(fwilson_vector),buf_length,fp) 
		!= buf_length)
	      {
		if(status == 0)
		  printf("%s: node %d propagator read error %d file %s\n",
			 myname,this_node,errno,filename); 
		fflush(stdout); 
		status = 1;
	      }
	    where_in_buf = 0;  /* reset counter */
	  }  /*** end of the buffer read ****/
	
	if(destnode==0){	
	  /* no need to send the Wilson vector - will just copy */
	  i = node_index(x,y,z,t);

	  if(dest_site == (field_offset)(-1))
	    dest = dest_field + i;
	  else
	    dest = (wilson_vector *)F_PT( &(lattice[i]), dest_site );

	  /* Copy to msg.wv for further processing */
	  msg.wv = lbuf[where_in_buf];
	}
	else {		
	  /* send to correct node */
	  msg.wv = lbuf[where_in_buf];
	  send_field((char *)&msg, sizeof(msg), destnode);
	}
	where_in_buf++;
      }
      
      /* The node which contains this site reads message */
      else {	
	/* for all nodes other than node 0, receive the message */
	if(this_node==destnode){
	  i = node_index(x,y,z,t);

	  if(dest_site == (field_offset)(-1))
	    dest = dest_field + i;
	  else
	    dest = (wilson_vector *)F_PT( &(lattice[i]), dest_site );
	  
	  /* Receive padded message in msg */
	  get_field((char *)&msg, sizeof(msg), 0);
	}
      }

      /* The receiving node does the byte reversal and then checksum,
         if needed.  At this point the wilson vector on the destnode
         is in msg.wv and dest points to its ultimate destination. */
      
      if(this_node==destnode)
	{
	  if(byterevflag==1)
	    byterevn((int32type *)&msg.wv,
		     sizeof(fwilson_vector)/sizeof(int32type));
	  /* Accumulate checksums */
	  for(k = 0, val = (u_int32type *)(&msg.wv); 
	      k < (int)sizeof(fwilson_vector)/(int)sizeof(int32type); 
	      k++, val++)
	    {
	      test_wpc.sum29 ^= (*val)<<rank29 | (*val)>>(32-rank29);
	      test_wpc.sum31 ^= (*val)<<rank31 | (*val)>>(32-rank31);
	      rank29++; if(rank29 >= 29)rank29 = 0;
	      rank31++; if(rank31 >= 31)rank31 = 0;
	    }

	  /* Convert to generic precision and put in destination */
	  f2d_wvec(&msg.wv,dest);
	}
      else
	{
	  rank29 += sizeof(fwilson_vector)/sizeof(int32type);
	  rank31 += sizeof(fwilson_vector)/sizeof(int32type);
	  rank29 %= 29;
	  rank31 %= 31;
	}
    }

  broadcast_bytes((char *)&status,sizeof(int));
  if(status != 0)return status;

  /* Combine node checksum contributions with global exclusive or */
  g_xor32(&test_wpc.sum29);
  g_xor32(&test_wpc.sum31);
  
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
	    }
/*	  else
	    printf("Checksums %x %x OK for spin %d color %d file %s\n",
		   wpf->check.sum29,wpf->check.sum31,
		   wpf->check.spin,wpf->check.color,wpf->filename); */
	}
      fflush(stdout);
      free(lbuf);
    }

  return 0;
  
} /* r_serial_w */

/*----------------------------------------------------------------------*/

/* Here only node 0 reads the Wilson propagator from a binary file */

int r_serial_w_to_site(w_prop_file *wpf, int spin, int color, 
		       field_offset dest_site)
{
  /* wpf  = propagator file structure 
     dest  = field offset for propagator Wilson vector  */

  /* 0 is normal exit code
     1 for seek, read error, or missing data error */

  return r_serial_w(wpf, spin, color, dest_site, NULL);
}

/*----------------------------------------------------------------------*/

/* Here only node 0 reads the Wilson propagator from a binary file */

int r_serial_w_to_field(w_prop_file *wpf, int spin, int color, 
			wilson_vector *dest_field)
{
  /* wpf  = propagator file structure 
     dest  = field offset for propagator Wilson vector  */

  /* 0 is normal exit code
     1 for seek, read error, or missing data error */

  return r_serial_w(wpf, spin, color, (field_offset)(-1), dest_field);
}

/*----------------------------------------------------------------------*/

void r_serial_w_f(w_prop_file *wpf)

/* Close the file and free associated structures */
{
  g_sync();
  if(this_node==0)
    {
      if(wpf->parallel == PARALLEL)
	printf("r_serial_w_f: Attempting serial close on parallel file \n");
      
      if(wpf->fp != NULL)fclose(wpf->fp);
/*      printf("Closed prop file %s\n",wpf->filename);*/
      fflush(stdout);
    }
  
  if(wpf->rank2rcv != NULL)free(wpf->rank2rcv);
  free(wpf->header);
  free(wpf);
  
} /* r_serial_w_f */

/*---------------------------------------------------------------------------*/

/* Write site list - only for checkpoint files */

void write_site_list_w(FILE *fp, w_prop_header *wph)
{
  off_t offset;
  int i;
  int buf_length;
  register site *s;
  int32type coords, *cbuf;

  /* All nodes write their site coordinate list in sequential
     blocks after the header.  The list is in the order of appearance
     in the lattice array.  Node 0 writes to the first block
     followed by node 1, etc.  The result is a contiguous table
     that can be used to remap the data to the corresponding space-time
     coordinate */

  /* Location of site list for this node */
  
  offset = wph->header_bytes + 
    sizeof(int32type)*sites_on_node*this_node;

  cbuf = (int32type *)malloc(sites_on_node*sizeof(int32type));
  if(cbuf == NULL)
    {
      printf("write_site_list_w: node %d can't malloc cbuf\n",this_node);
      fflush(stdout);terminate(1);
    }

  if( g_seek(fp,offset,SEEK_SET) < 0 ) 
    {
      printf("write_site_list_w: Node %d fseek %ld failed errno %d\n",
	  this_node,(long)offset,errno);
      fflush(stdout);terminate(1);   
    }
  
  buf_length = 0;

  FORALLSITES(i,s)
    {
      /* Encode the space-time coordinate vector as a 32-bit integer */
      coords = nx*(ny*(nz*s->t + s->z) + s->y) + s->x;
      cbuf[buf_length] = coords;
      buf_length++;
    }

  if( (int)g_write(cbuf,sizeof(int32type),sites_on_node,fp) != sites_on_node)
    {
      printf("write_site_list_w: Node %d coords write error %d\n",
	     this_node,errno);fflush(stdout);terminate(1);
    }

  free(cbuf);
} /* write_site_list_w */

/*---------------------------------------------------------------------------*/

/* Open a file for parallel writing */
w_prop_file *parallel_open_w(int order, char *filename)
{
  /* All nodes open the same filename */
  /* Returns a file structure describing the opened file */

  /* order = NATURAL_ORDER for coordinate natural order 
           = NODE_DUMP_ORDER for node-dump order */

  FILE *fp;
  w_prop_file *wpf;
  w_prop_header *wph;

  /* Set up propagator file and header structures and load header values */
  wpf = setup_output_w_prop_file();
  wph = wpf->header;

  wph->order = order;

  /* All nodes open the requested file */

  fp = g_open(filename, "wb");
  if(fp == NULL)
    {
      printf("parallel_open_w: Node %d can't open file %s, error %d\n",
	     this_node,filename,errno);fflush(stdout);terminate(1);
    }

/*  if(this_node==0)printf("Opened prop file %s for parallel writing\n",filename); */
      


  /* Node 0 writes the header */

  if(this_node==0)
    pwrite_w_prop_hdr(fp,wph);

  broadcast_bytes((char *)&wph->header_bytes,sizeof(wph->header_bytes));

  /* All nodes write site list to file */

  write_site_list_w(fp,wph);
  
  /* Assign values to file structure */

  wpf->fp             = fp;
  wpf->filename        = filename;
  wpf->byterevflag    = 0;            /* Not used for writing */
  wpf->parallel       = 1;            /* File opened in parallel */

  /* Node 0 writes ascii info file */

  if(this_node == 0)write_w_prop_info_file(wpf);

  return wpf;

} /* parallel_open_w */

/*-----------------------------------------------------------------------*/

/* Open a file for parallel writing in natural order */
w_prop_file *w_parallel_w_i(char *filename)
{
  /* All nodes open the same filename */
  /* Returns a file structure describing the opened file */

  return parallel_open_w(NATURAL_ORDER,filename);

} /* w_parallel_w_i */


/*---------------------------------------------------------------------------*/
/* Open a file for parallel writing in node-dump order */
w_prop_file *w_checkpoint_w_i(char *filename)
{
  /* All nodes open the same filename */
  /* Returns a file structure describing the opened file */

  return parallel_open_w(NODE_DUMP_ORDER,filename);

} /* w_checkpoint_w_i */

/*---------------------------------------------------------------------------*/

void w_parallel_w_o(w_prop_file *wpf)
{
  /* Reopen previously opened file, specified by wpf, for parallel writing */

  if(!wpf->parallel == PARALLEL)
    printf("w_parallel_w_o: Attempting parallel reopen on file previously opened as serial.\n");
  
  wpf->fp = g_open(wpf->filename, "rb+");  /* Binary write with update */
  
  if(wpf->fp == NULL)
    {
      printf("w_parallel_w_o: Node %d can't reopen file %s, error %d\n",
	     this_node,wpf->filename,errno);fflush(stdout);
      terminate(1);
    }

/*  if(this_node==0)printf("Reopened prop file %s for parallel writing.\n",
			 wpf->filename); */
  wpf->parallel = PARALLEL;

} /* w_parallel_w_o */

/*---------------------------------------------------------------------------*/

void w_checkpoint_w_o(w_prop_file *wpf)
{
  w_parallel_w_o(wpf);
}
/*---------------------------------------------------------------------------*/

/* Position propagator file for writing in parallel */

fwilson_vector *w_parallel_setup_w(w_prop_file *wpf, off_t *checksum_offset,
			       int spin, int color)
{
  /* wpf  = file descriptor as opened by w_parallel_w_i or w_checkpoint_w_i */

  FILE *fp;
  w_prop_header *wph;
  fwilson_vector *lbuf;
  off_t offset ;            /* File stream pointer */
  off_t w_prop_node_size;   /* Size of a w_prop configuration block for
                              all sites on one node */
  off_t coord_list_size;    /* Size of coordinate list in bytes */
  off_t head_size;          /* Size of header plus coordinate list */
  off_t w_prop_check_size;  /* Size of propagator checksum record */
  off_t w_prop_size;        /* Size of propagator blocks for all nodes */
  off_t body_size ;         /* Size of propagator blocks for all nodes 
			      plus checksum record */
  int spinindex;
  char myname[] = "w_parallel_setup_w";

  if(!wpf->parallel == PARALLEL)
    printf("%s: Attempting parallel write to serial file.\n",myname);

  lbuf = (fwilson_vector *)malloc(MAX_BUF_LENGTH*sizeof(fwilson_vector));
  if(lbuf == NULL)
    {
      printf("%s: Node %d can't malloc lbuf\n",myname,this_node);
      fflush(stdout);
      terminate(1);
    }

  fp = wpf->fp;
  wph = wpf->header;

  w_prop_node_size = sites_on_node*sizeof(fwilson_vector) ;
  w_prop_size = volume*sizeof(fwilson_vector) ;
  w_prop_check_size =  sizeof(wpf->check.spin) +
     sizeof(wpf->check.color) +
     sizeof(wpf->check.sum29) + sizeof(wpf->check.sum31);
  body_size = w_prop_size + w_prop_check_size;

  if(wph->order == NATURAL_ORDER)coord_list_size = 0;
  else coord_list_size = sizeof(int32type)*volume;
  head_size = wph->header_bytes + coord_list_size;

  /* Find spin in table of contents */
  for(spinindex=0 ; spinindex < wph->n_spins ; spinindex++)
    if(wph->spins[spinindex]==spin)break;
  
  if(spinindex == wph->n_spins)
    {
      printf("w_parallel_w: Requested spin %d not planned for file %s\n",
	     spin,wpf->filename);
      printf("  Table of contents: ");
      for(spinindex=0; spinindex<wph->n_spins; spinindex++)
	printf(" %d",wph->spins[spinindex]);
      printf("\n");fflush(stdout);
      terminate(1);
    }

  /* Offset for propagator check data */

  *checksum_offset = head_size + body_size*(spinindex*3 + color);

  /* Offset for Wilson vector data for this node */
  
  offset = *checksum_offset + w_prop_check_size + w_prop_node_size*this_node;

  if( g_seek(fp,offset,SEEK_SET) < 0 ) 
    {
      printf("%s: Node %d fseek %ld failed error %d file %s\n",
	     myname,this_node,(long)offset,errno,wpf->filename);
      fflush(stdout);terminate(1);
    }

  return lbuf;
} /* w_parallel_setup_w */

/*---------------------------------------------------------------------------*/

/* Write parallel propagator for this spin and color */

void w_parallel_w(w_prop_file *wpf, int spin, int color, field_offset src_site,
		  wilson_vector *src_field)
{
  /* wpf  = file descriptor as opened by w_parallel_w_i 
     src_site  = field offset for propagator Wilson vector
     src_field  = pointer to Wilson vector */

  FILE *fp;
  fwilson_vector *lbuf;
  wilson_vector *src;
  int buf_length,where_in_buf;
  u_int32type *val;
  int rank29,rank31;
  off_t checksum_offset;
  register int i;
  int k;
  int x,y,z,t;
  struct {
    short x,y,z,t;
    fwilson_vector wv;
    char pad[PAD_SEND_BUF];    /* Introduced because some switches
				  perform better if message lengths are longer */
  } msg;
  int isite,ksite,site_block;
  int rcv_coords,rcv_rank;
  int destnode,sendnode;
  int c,s;
  char myname[] = "w_parallel_w";

  fp = wpf->fp;

  lbuf = w_parallel_setup_w(wpf,&checksum_offset,spin,color);

  /* initialize checksums */
  wpf->check.sum31 = 0;
  wpf->check.sum29 = 0;

  /* Deal and write */

  g_sync();
  buf_length = 0;

  /* Clear buffer as a precaution.  Easier to tell if we botch the
     buffer loading. */
  for(i=0;i<MAX_BUF_LENGTH;i++)
    for(s=0;s<4;s++)for(c=0;c<3;c++){
      lbuf[i].d[s].c[c].real = 0.; 
      lbuf[i].d[s].c[c].imag = 0.;
    }

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
	      /* Message consists of site coordinates and Wilson vector */
	      msg.x = x; msg.y = y; msg.z = z; msg.t = t;
	      i = node_index(x,y,z,t);

	      if(src_site == (field_offset)(-1))
		src = src_field + i;
	      else
		src = (wilson_vector *)F_PT( &(lattice[i]), src_site );

	      /* Copy from site structure to msg structure, converting
		 to single precision */

	      d2f_wvec(src,&msg.wv);
	      /* Then send the message */
	      send_field((char *)&msg,sizeof(msg),destnode);
	    }
	    /* Node destnode copies or receives a message */
	    else if(this_node==destnode){
	      if(destnode==sendnode){ 
		/* just copy Wilson vector */
		i = node_index(x,y,z,t);
		if(src_site == (field_offset)(-1))
		  src = src_field + i;
		else
		  src = (wilson_vector *)F_PT( &(lattice[i]), src_site );
		where_in_buf = buf_length;
		/* Copy from site structure to write buffer,
		   converting to single precision */
		d2f_wvec(src,&lbuf[where_in_buf]);

		rank29 = rank31 = 
		  sizeof(fwilson_vector)/sizeof(int32type)*rcv_rank;
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

		/* Move data from msg structure to write buffer */
		lbuf[where_in_buf] = msg.wv;
		rank29 = rank31 = sizeof(fwilson_vector)/sizeof(int32type)*i;
	      }

	      /* Receiving node accumulates checksums as the values
		 are inserted into its buffer */
	      rank29 %= 29; rank31 %= 31;
	      for(k = 0, val = (u_int32type *)&lbuf[where_in_buf]; 
		  k < (int)sizeof(fwilson_vector)/(int)sizeof(int32type); 
		  k++, val++)
		{
		  wpf->check.sum29 ^= (*val)<<rank29 | (*val)>>(32-rank29);
		  wpf->check.sum31 ^= (*val)<<rank31 | (*val)>>(32-rank31);
		  rank29++; if(rank29 >= 29)rank29 = 0;
		  rank31++; if(rank31 >= 31)rank31 = 0;
		}

	      buf_length++;
	      if( (buf_length == MAX_BUF_LENGTH) || 
		 (isite == sites_on_node -1))
		{
		  /* write out buffer */
		  
		  if( (int)g_write(lbuf,sizeof(fwilson_vector),buf_length,fp) 
		     != buf_length)
		    {
		      printf("%s: Node %d Wilson propagator write error %d file %s\n",
			     myname,this_node,errno,wpf->filename); 
		      fflush(stdout);
		      terminate(1);   
		    }
		  buf_length = 0;		/* start again after write */
		  /* Clear buffer as a precaution */
		  for(i=0;i<MAX_BUF_LENGTH;i++)
		    for(s=0;s<4;s++)for(c=0;c<3;c++){
		      lbuf[i].d[s].c[c].real = 0.; 
		      lbuf[i].d[s].c[c].imag = 0.;
		    }
		}
	    } /* else if(this_node==destnode) */
	    
	  } /* destnode, isite */
      g_sync();  /* To assure all write buffers are completed before
		   starting on the next buffer */
    } /* ksite */

  
  free(lbuf);

  /* Combine checksums */

  g_xor32(&wpf->check.sum29);
  g_xor32(&wpf->check.sum31);
  wpf->check.spin = spin;
  wpf->check.color = color;

  /* Write checksum at end of lattice file */

  if(this_node == 0)
    {
/*      printf("Wrote prop in parallel for spin %d color %d to file %s\n",
	     spin,color,wpf->filename); */

      /* Node 0 positions file for writing checksum */
      
      if( g_seek(fp,checksum_offset,SEEK_SET) < 0 ) 
	{
	  printf("%s: Node 0 g_seek %ld for checksum failed error %d file %s\n",
		 myname,(long)checksum_offset,errno,wpf->filename);
	  fflush(stdout);terminate(1);   
	}
      
      write_checksum_w(PARALLEL,wpf);
      fflush(stdout);
    }

} /* w_parallel_w */
        
/*---------------------------------------------------------------------------*/

/* Write parallel propagator for this spin and color */

void w_parallel_w_from_site(w_prop_file *wpf, int spin, int color, 
			    field_offset src_site)
{
  /* wpf  = file descriptor as opened by w_parallel_w_i 
     src  = field offset for propagator Wilson vector  */

  w_parallel_w(wpf, spin, color, src_site, NULL);
}

/*---------------------------------------------------------------------------*/

/* Write parallel propagator for this spin and color */

void w_parallel_w_from_field(w_prop_file *wpf, int spin, int color, 
			     wilson_vector *src_field)
{
  /* wpf  = file descriptor as opened by w_parallel_w_i 
     src  = field offset for propagator Wilson vector  */

  w_parallel_w(wpf, spin,color, (field_offset)(-1), src_field);
}

/*---------------------------------------------------------------------------*/

/* Write checkpoint propagator for this spin and color */

void w_checkpoint_w(w_prop_file *wpf, int spin, int color, 
		    field_offset src_site, wilson_vector *src_field)
{
  /* wpf  = file descriptor as opened by w_checkpoint_w_i 
     src  = field offset for propagator Wilson vector  */

  FILE *fp;
  fwilson_vector *lbuf;
  wilson_vector *src;
  u_int32type *val;
  int k;
  int rank29,rank31;
  off_t checksum_offset;
  int buf_length;
  register site *s;
  register int i;
  char myname[] = "w_checkpoint_w";

  fp = wpf->fp;

  lbuf = w_parallel_setup_w(wpf,&checksum_offset,spin,color);

  /* C. McNeile's algorithm, changed slightly*/

  /* initialize checksums */
  wpf->check.sum31 = 0;
  wpf->check.sum29 = 0;
  /* counts 32-bit words mod 29 and mod 31 in order of appearance on file */
  /* Here all nodes use these values */
  rank29 = sizeof(fwilson_vector)/sizeof(int32type)*sites_on_node*this_node % 29;
  rank31 = sizeof(fwilson_vector)/sizeof(int32type)*sites_on_node*this_node % 31;

  buf_length = 0;

  FORALLSITES(i,s)
  {
        
    /* load the quark propagator into the write buffer, converting to
       single precision */

    if(src_site == (field_offset)(-1))
      src = src_field + i;
    else
      src = (wilson_vector *)F_PT( s, src_site );

    d2f_wvec(src,&lbuf[buf_length]);

    /* Accumulate checksums - contribution from next site moved into buffer*/
    for(k = 0, val = (u_int32type *)&lbuf[buf_length]; 
	k < (int)sizeof(fwilson_vector)/(int)sizeof(int32type); k++, val++)
      {
	wpf->check.sum29 ^= (*val)<<rank29 | (*val)>>(32-rank29);
	wpf->check.sum31 ^= (*val)<<rank31 | (*val)>>(32-rank31);
	rank29++; if(rank29 >= 29)rank29 = 0;
	rank31++; if(rank31 >= 31)rank31 = 0;
      }

    buf_length++;
    
    if( (buf_length == MAX_BUF_LENGTH) || (i == sites_on_node -1))
      {
	/* write out buffer */
	
	if( (int)g_write(lbuf,sizeof(fwilson_vector),buf_length,fp) != buf_length)
	  {
	    printf("%s: Node %d propagator write error %d file %s\n",
		   myname,this_node,errno,wpf->filename); 
	    fflush(stdout);
	    terminate(1);   
	  }
	buf_length = 0;		/* start again after write */
      }
    
  } 
  
  free(lbuf);

  /* Combine checksums */

  g_xor32(&wpf->check.sum29);
  g_xor32(&wpf->check.sum31);
  wpf->check.spin = spin;
  wpf->check.color = color;

  /* Write checksum at end of lattice file */

  if(this_node==0)
    {
/*      printf("Wrote prop for spin %d color %d to checkpoint file %s\n",
	     spin,color,wpf->filename);fflush(stdout); */

      /* Node 0 positions file for writing checksum */
      
      if( g_seek(fp,checksum_offset,SEEK_SET) < 0 ) 
	{
	  printf("%s: Node 0 g_seek %ld for checksum failed error %d file %s\n",
		 myname,(long)checksum_offset,errno,wpf->filename);
	  fflush(stdout);terminate(1);   
	}
      
      write_checksum_w(PARALLEL,wpf);
    }
      
} /* w_checkpoint_w */

/*---------------------------------------------------------------------------*/

/* Write checkpoint propagator for this spin and color */

void w_checkpoint_w_from_site(w_prop_file *wpf, int spin, int color, 
			      field_offset src_site)
{
  /* wpf  = file descriptor as opened by w_checkpoint_w_i 
     src  = field offset for propagator Wilson vector  */

  w_checkpoint_w(wpf, spin, color, src_site, NULL);
}

/*---------------------------------------------------------------------------*/

/* Write checkpoint propagator for this spin and color */

void w_checkpoint_w_from_field(w_prop_file *wpf, int spin, int color, 
			       wilson_vector *src_field)
{
  /* wpf  = file descriptor as opened by w_checkpoint_w_i 
     src  = field offset for propagator Wilson vector  */

  w_checkpoint_w(wpf, spin, color, (field_offset)(-1), src_field);
}

/*---------------------------------------------------------------------------*/

void w_parallel_w_c(w_prop_file *wpf)
{
  /* Close file temporarily -- releases only the system FILE structure */
  g_sync();
  if(!wpf->parallel == PARALLEL)
    printf("w_parallel_w_c: Attempting parallel close on serial file.\n");

  if(wpf->fp != NULL)g_close(wpf->fp);
  wpf->fp = NULL;
/*  if(this_node==0)printf("Closed prop file %s temporarily\n",wpf->filename);*/
} /* w_parallel_w_c */

/*---------------------------------------------------------------------------*/

void w_checkpoint_w_c(w_prop_file *wpf)
{
  w_parallel_w_c(wpf);
}
/*---------------------------------------------------------------------------*/

void w_parallel_w_f(w_prop_file *wpf)
{
  /* Close file (if still active) and release header and file structures */
  if(wpf->fp != NULL)w_parallel_w_c(wpf);
  if(this_node==0)
    printf("Wrote prop file %s time stamp %s\n",wpf->filename,
	   (wpf->header)->time_stamp);
  free(wpf->header);
  free(wpf);
} /* w_parallel_w_f */


/*---------------------------------------------------------------------------*/

void w_checkpoint_w_f(w_prop_file *wpf)
{
  w_parallel_w_f(wpf);
}

/*---------------------------------------------------------------------------*/

w_prop_file *r_parallel_w_i(char *filename)
{
  /* Returns file descriptor for opened file */

  w_prop_header *wph;
  w_prop_file *wpf;
  FILE *fp;
  int byterevflag;

  /* All nodes set up a propagator file and header structure for reading */

  wpf = setup_input_w_prop_file(filename);
  wph = wpf->header;

  wpf->parallel = 1;   /* File was opened for parallel access */

  /* All nodes open a file */

  fp = g_open(filename, "rb");
  if(fp == NULL)
    {
      printf("r_parallel_w_i: Node %d can't open file %s, error %d\n",
	     this_node,filename,errno);fflush(stdout);terminate(1);
    }

/*  if(this_node==0)printf("Opened prop file %s for parallel reading\n",
			 filename); */
  
  wpf->fp = fp;

  /* Node 0 reads header */

  if(this_node==0)
    byterevflag = read_w_prop_hdr(wpf,PARALLEL);
  
  /* Broadcast the byterevflag from node 0 to all nodes */

  broadcast_bytes((char *)&byterevflag,sizeof(byterevflag));

  wpf->byterevflag = byterevflag;

  /* Broadcasts the header structure from node 0 to all nodes */
  
  broadcast_bytes((char *)wph,sizeof(w_prop_header));

  /* Read site list and broadcast to all nodes. */

  read_site_list_w(PARALLEL,wpf);

  return wpf;

} /* r_parallel_w_i */

/*----------------------------------------------------------------------*/

void r_parallel_w_o(w_prop_file *wpf)
{
  /* Reopen previously opened file, specified by wpf, for parallel reading */

  if(!wpf->parallel == PARALLEL)
    printf("r_parallel_w_o: Attempting parallel open on file previously opened as serial.\n");

  wpf->fp = g_open(wpf->filename, "rb");  /* Binary read */

  if(wpf->fp == NULL)
    {
      printf("r_parallel_w_o: Node %d can't reopen file %s, error %d\n",
	     this_node,wpf->filename,errno);fflush(stdout);
      terminate(1);
    }
  
/*  if(this_node==0)printf("Reopened prop file %s for parallel reading\n",
			 wpf->filename); */
  wpf->parallel = PARALLEL;
}

/*----------------------------------------------------------------------*/

/* Read Wilson propagator in parallel from a single file */
int r_parallel_w(w_prop_file *wpf, int spin, int color, field_offset dest_site,
		 wilson_vector *dest_field)
{
  /* wpf  = propagator file structure 
     dest_site  = field offset for propagator Wilson vector
     dest_field = pointer to field
     use only one of the dest parameters!  */

  /* 0 is normal exit code
     1 for seek, read error, or missing data error */

  FILE *fp;
  w_prop_header *wph;
  char *filename;
  int byterevflag;
  int spinindex;            /* Counts spin records in file   -
			       wph->spins[spinindex] = spin */
  fwilson_vector *lbuf;
  wilson_vector *dest;
  struct {
    short x,y,z,t;
    fwilson_vector wv;
  } msg;

  int buf_length,where_in_buf;
  w_prop_check test_wpc;
  u_int32type *val;
  int rank29,rank31;
  int destnode,sendnode,isite,ksite,site_block;
  int x,y,z,t;
  int rcv_rank,rcv_coords;
  register int i,k;
  int status;
  Real xstatus;

  off_t offset ;            /* File stream pointer */
  off_t w_prop_node_size;   /* Size of a propagator block for all sites on
                              one node */
  off_t w_prop_size;        /* Size of propagator blocks for all nodes */
  off_t w_prop_check_size;  /* Size of propagator checksum record */
  off_t coord_list_size;    /* Size of coordinate list in bytes */
  off_t head_size;          /* Size of header plus coordinate list */
  off_t body_size ;         /* Size of propagator blocks for all nodes 
			      plus checksum record */
  char myname[] = "r_parallel_w";

  fp = wpf->fp;
  wph = wpf->header;
  filename = wpf->filename;
  byterevflag = wpf->byterevflag;

  if(!wpf->parallel == PARALLEL)
    printf("%s: Attempting parallel read from serial file.\n",myname);

  /* Allocate read buffer */
  lbuf = (fwilson_vector *)malloc(MAX_BUF_LENGTH*sizeof(fwilson_vector));
  if(lbuf == NULL)
    {
      printf("%s: Node %d can't malloc lbuf\n",myname,this_node);
      fflush(stdout);
      terminate(1);
    }

  w_prop_node_size = sites_on_node*sizeof(fwilson_vector) ;
  w_prop_size = volume*sizeof(fwilson_vector) ;
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

  status = 0;
  if(spinindex == wph->n_spins && this_node == 0)
    {
      printf("%s: Requested spin %d not in file %s\n",myname,spin,filename);
      printf("  Table of contents: ");
      for(spinindex=0; spinindex<wph->n_spins; spinindex++)
	printf(" %d",wph->spins[spinindex]);
      printf("\n");fflush(stdout);
      status = 1;
    }

  broadcast_bytes((char *)&status,sizeof(int));
  if(status != 0)return status;

  if(wpf->header->order == NATURAL_ORDER)coord_list_size = 0;
  else coord_list_size = sizeof(int32type)*volume;
  head_size = wpf->header->header_bytes + coord_list_size;

  offset = head_size + body_size*(spinindex*3 + color);

  /* Only node 0 reads and verifies check record */

  status = 0;
  if(this_node == 0)
    {
      /* Position file pointer for reading check record */

      if( g_seek(fp,offset,SEEK_SET) < 0 ) 
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

      status += pread_byteorder(byterevflag,fp,&wpf->check.spin,
		      sizeof(wpf->check.spin),myname,"check.spin");
      status += pread_byteorder(byterevflag,fp,&wpf->check.color,
		      sizeof(wpf->check.color),myname,"check.color");
      status += pread_byteorder(byterevflag,fp,&wpf->check.sum29,
		      sizeof(wpf->check.sum29),myname,"check.sum29");
      /* 1996 format had an unused 32-bit checksum.
	 Version 5 format has two 32-bit checksums */
      if(wph->magic_number == W_PROP_VERSION_NUMBER)
	status += pread_byteorder(byterevflag,fp,&wpf->check.sum31,
			sizeof(wpf->check.sum31),myname,"check.sum31");

      /* Verify spin and color - checksums come later */
      
      if(wpf->check.spin != spin || wpf->check.color != color)
	{
	  printf("%s: Spin %d and color %d do not match check record on file %s\n",
		 myname,spin,color,filename);
	  printf("  Check record said %d %d\n",wpf->check.spin,wpf->check.color);
	  fflush(stdout);
	  status = 1;
	}

    }

  broadcast_bytes((char *)&status,sizeof(int));
  if(status != 0)return status;
  
  /* Position file for reading propagator */
  
  offset += w_prop_check_size;

  /* Each node reads its propagator values */

  offset += w_prop_node_size*this_node;
  
  xstatus = 0.;
  if( g_seek(fp,offset,SEEK_SET) < 0 ) 
    {
      printf("%s: Node %d fseek %ld failed error %d file %s\n",
	     myname,this_node,(long)offset,errno,filename);
      fflush(stdout);
      xstatus = 1.;
    }
  g_floatsum(&xstatus);
  if(xstatus != 0.)return 1;

  /* initialize checksums */
  test_wpc.sum31 = 0;
  test_wpc.sum29 = 0;
  /* counts 32-bit words mod 29 and mod 31 in order of appearance on file */
  /* Here all nodes use these values */
  rank29 = sizeof(fwilson_vector)/sizeof(int32type)*sites_on_node*this_node % 29;
  rank31 = sizeof(fwilson_vector)/sizeof(int32type)*sites_on_node*this_node % 31;

  /* Read and deal */

  g_sync();
  buf_length = 0;
  where_in_buf = 0;
  
  /* All nodes participate in reading and dealing.
     Cycle through nodes, dealing "site_block" values from each node
     in sequence.  (We haven't worked hard at optimizing this value.)
     We sync at the end of each block to prevent message pileups
     on slow nodes.  This is especially necessary on time-sharing
     machines like the SGI Cray Origin 2000.

     It is possible that messages arrive at a node in an order
     different from the order of dealing so we include the destination
     site coordinates in each message to specify where it goes */

  xstatus = 0.;
  site_block = 16;
  for(ksite=0; ksite<sites_on_node; ksite += site_block)
    {
      for(sendnode=0; sendnode<number_of_nodes; sendnode++)
	for(isite=ksite; 
	    isite<sites_on_node && isite<ksite+site_block; isite++)
	  {
	    /* Compute destination coordinate for the next field 
	       that this node reads and sends
	       
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
	    
	    if(wpf->header->order == NATURAL_ORDER)
	      rcv_coords = rcv_rank;
	    else
	      rcv_coords = wpf->rank2rcv[rcv_rank];
	    
	    /* Decode site coordinates */
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
		  
		  if( (int)g_read(lbuf,buf_length*sizeof(fwilson_vector),1,fp) != 1)
		    {
		      if(xstatus == 0.)
			printf("%s: node %d propagator read error %d file %s\n",
			       myname,this_node,errno,filename); 
		      fflush(stdout); 
		      xstatus = 1.;
		    }
		  where_in_buf = 0;  /* reset counter */
		}  /*** end of the buffer read ****/
	      
	      /* Do byte reversal if needed */
	      if(wpf->byterevflag==1)
		byterevn((int32type *)&lbuf[where_in_buf],
			 sizeof(fwilson_vector)/sizeof(int32type));
	      
	      /* Accumulate checksums - contribution from next site */
	      for(k = 0, val = (u_int32type *)&lbuf[where_in_buf]; 
		  k < (int)sizeof(fwilson_vector)/(int)sizeof(int32type); 
		  k++, val++)
		{
		  test_wpc.sum29 ^= (*val)<<rank29 | (*val)>>(32-rank29);
		  test_wpc.sum31 ^= (*val)<<rank31 | (*val)>>(32-rank31);
		  rank29++; if(rank29 >= 29)rank29 = 0;
		  rank31++; if(rank31 >= 31)rank31 = 0;
		}
	      
	      if(destnode==sendnode){	
		/* just copy Wilson vector */
		i = node_index(x,y,z,t);

		if(dest_site == (field_offset)(-1))
		  dest = dest_field + i;
		else
		  dest = (wilson_vector *)F_PT( &(lattice[i]), dest_site );

		/* copy from read buffer to site structure, converting
		   to double precision */
		f2d_wvec(&lbuf[where_in_buf],dest);
	      }
	      else {		/* send to correct node */
		/* Message consists of site coordinates and Wilson vector */
		msg.x = x; msg.y = y; msg.z = z; msg.t = t;
		msg.wv = lbuf[where_in_buf];
		
		send_field((char *)&msg,sizeof(msg),destnode);
	      }
	      where_in_buf++;
	    }
	    /* The node which contains this site reads a message */
	    else {	/* for all nodes other than node sendnode */
	      if(this_node==destnode){
		get_field((char *)&msg,sizeof(msg),sendnode);
		i = node_index(msg.x,msg.y,msg.z,msg.t);

		if(dest_site == (field_offset)(-1))
		  dest = dest_field + i;
		else
		  dest = (wilson_vector *)F_PT( &(lattice[i]), dest_site );

		if(this_node!= node_number(msg.x,msg.y,msg.z,msg.t))
		  {
		    printf("BOTCH. Node %d received %d %d %d %d\n",
			   this_node,msg.x,msg.y,msg.z,msg.t);
		    fflush(stdout); terminate(1);
		  }
		/* Copy from message structure to site structure,
		   converting to double precision */
		f2d_wvec(&msg.wv,dest);
	      }
	    }
	  }  /** sendnode, isite **/
      g_sync();  /* To prevent message pileups on slow nodes */
    } /** ksite **/
  
  g_floatsum(&xstatus);
  if(xstatus != 0)return 1;

  free(lbuf);

  /* Read and verify checksum */
  /* Checksums implemented with version 5 */
  
  /* Combine node checksum contributions with global exclusive or */
  g_xor32(&test_wpc.sum29);
  g_xor32(&test_wpc.sum31);

  if(this_node == 0)
    {
/*      printf("Read prop in parallel for spin %d color %d in from file %s\n",
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
	    }
/*	  else
	    printf("Checksums %x %x OK for spin %d color %d file %s\n",
		   wpf->check.sum29,wpf->check.sum31,
		   wpf->check.spin,wpf->check.color,wpf->filename); */
	}
      fflush(stdout);
    }

  return 0;
  
} /* r_parallel_w */

/*----------------------------------------------------------------------*/

/* Read Wilson propagator in parallel from a single file */
int r_parallel_w_to_site(w_prop_file *wpf, int spin, int color, 
			 field_offset dest_site)
{
  /* wpf  = propagator file structure 
     dest_site  = field offset for propagator Wilson vector  */

  /* 0 is normal exit code
     1 for seek, read error, or missing data error */
  return r_parallel_w(wpf, spin, color, dest_site, NULL);
}

/*----------------------------------------------------------------------*/

/* Read Wilson propagator in parallel from a single file */
int r_parallel_w_to_field(w_prop_file *wpf, int spin, int color, 
			  wilson_vector *dest_field)
{
  /* wpf  = propagator file structure 
     dest_field  = pointer to Wilson vector  */

  /* 0 is normal exit code
     1 for seek, read error, or missing data error */

  return r_parallel_w(wpf, spin, color, (field_offset)(-1), dest_field);
}

/*---------------------------------------------------------------------------*/

void r_parallel_w_c(w_prop_file *wpf)
{
  /* Close file temporarily -- releases only the system FILE structure */
  g_sync();
  if(!wpf->parallel == PARALLEL)
    printf("r_parallel_w_c: Attempting parallel close on serial file.\n");
  if(wpf->fp != NULL)g_close(wpf->fp);
  wpf->fp = NULL;
/*  if(this_node==0)printf("Closed prop file %s temporarily\n",wpf->filename);*/
} /* r_parallel_w_c */

/*---------------------------------------------------------------------------*/

void r_parallel_w_f(w_prop_file *wpf)
{
  /* Close file (if active) and release header and file structures */

  if(wpf->fp != NULL)r_parallel_w_c(wpf);
  if(wpf->rank2rcv != NULL)free(wpf->rank2rcv);
  free(wpf->header);
  free(wpf);
/*  if(this_node==0)printf("Released prop file %s\n",wpf->filename); */
 } /* r_parallel_w_f */

/*---------------------------------------------------------------------------*/

/* ASCII file format

   format:
    version_number (int)
    time_stamp (char string enclosed in quotes)
    nx ny nz nt (int)
  The rest comes in the order of writing and must be reread in the
  same order:
    source_spin, source_color (int) (for values to follow)
    for(t=...)for(z=...)for(y=...)for(x=...){
        for each wilson_vector:
            for(i=...)for(j=...){p.d[i].c[j].real, p.d[i].c[j].imag}}}
*/

/*---------------------------------------------------------------------------*/
/* Open and write header info for ascii propagator file */

w_prop_file *w_ascii_w_i(char *filename)
{
  w_prop_header *wph;
  w_prop_file *wpf;
  FILE *fp;

  wpf = setup_output_w_prop_file();
  wph = wpf->header;

  /* node 0 does all the writing */
  if(this_node==0){
    /* Set up w_prop file and w_prop header structures and load header values */

    fp = fopen(filename,"w");
    if(fp==NULL){
      printf("Can't open file %s, error %d\n",filename,errno);terminate(1);
    }

    wpf->fp = fp;

    if( (fprintf(fp,"%d\n",W_PROP_VERSION_NUMBER))==0 ){
      printf("Error in writing version number\n"); terminate(1);
    }
    if( (fprintf(fp,"\"%s\"\n",wph->time_stamp))==0 ){
      printf("Error in writing time stamp\n"); terminate(1);
    }
    
    if( (fprintf(fp,"%d\t%d\t%d\t%d\n",nx,ny,nz,nt))==0 ){
      printf("Error in writing dimensions\n"); terminate(1);
    }

  }
  else wpf->fp = NULL;

  /* Assign remaining values to propagator file structure */
  wpf->parallel = 0;
  wpf->filename       = filename;
  wpf->rank2rcv       = NULL;         /* Not used for writing */
  wpf->byterevflag    = 0;            /* Not used for writing */

  /* Node 0 writes info file */
  if(this_node==0)write_w_prop_info_file(wpf);

  return wpf;
} /* w_ascii_w_i */
  
/*---------------------------------------------------------------------------*/
/* Write ASCII propagator */

void w_ascii_w(w_prop_file *wpf,int spin,int color,field_offset src)
{
  FILE *fp;
  int currentnode,newnode;
  int i,j,l,x,y,z,t;
  wilson_vector lbuf;
  int node0=0;
  g_sync();
  currentnode=0;

  fp = wpf->fp;
  if(this_node==0)
    /* first the color and spin */
    {
      if( (fprintf(fp,"%d %d\n",spin,color))== EOF)
	{
	  printf("w_ascii_w: error writing spin,color\n"); 
	  terminate(1);
	} 
      fflush(fp);
    }
  /* next the elements */
  
  for(t=0;t<nt;t++)for(z=0;z<nz;z++)for(y=0;y<ny;y++)for(x=0;x<nx;x++)
    {
      newnode=node_number(x,y,z,t);
      if(newnode != currentnode)
	{	/* switch to another node */
	  g_sync();
	  currentnode=newnode;
	}
      
      if(this_node==0)
	{
	  if(currentnode==0)
	    {
	      l=node_index(x,y,z,t);
	      lbuf = *( (wilson_vector *)F_PT( &(lattice[l]), src ) );
	    }
	  else
	    {
	      get_field((char *)&lbuf,sizeof(wilson_vector),currentnode);
	    }
	  for(i=0;i<4;i++)for(j=0;j<3;j++)
	    {
	      if( (fprintf(fp,"%.7e\t%.7e\n",(double)lbuf.d[i].c[j].real,
			   (double)lbuf.d[i].c[j].imag))== EOF)
		{
		  printf("w_ascii_w: error writing prop\n"); 
		  terminate(1);
		} 
	    }
	}
      else
	{	/* for nodes other than 0 */
	  if(this_node==currentnode)
	    {
	      l=node_index(x,y,z,t);
	      lbuf = *( (wilson_vector *)F_PT( &(lattice[l]), src ) );
	      send_field((char *)&lbuf,sizeof(wilson_vector),node0);
	    }
	}
    }
  g_sync();
  if(this_node==0)
    {
      fflush(fp);
/*      printf("Wrote prop for spin %d color %d to ASCII file  %s\n",
	     spin,color,wpf->filename); */
    }

} /* w_ascii_w */
/*---------------------------------------------------------------------------*/
/* Close ASCII propagator file */
void w_ascii_w_f(w_prop_file *wpf)
{
  g_sync();
  if(this_node==0)
    {
      fflush(wpf->fp);
      if(wpf->fp != NULL)fclose(wpf->fp);
      printf("Wrote prop file %s time stamp %s\n",wpf->filename,
	     (wpf->header)->time_stamp);
    }

  /* Free header and file structures */
  free(wpf->header);
  free(wpf);
}



/*---------------------------------------------------------------------------*/
/* Open ASCII propagator file and read header information */

w_prop_file *r_ascii_w_i(char *filename)
{
  w_prop_file *wpf;
  w_prop_header *wph;
  FILE *fp;

  /* All nodes set up a propagator file and propagator header
     structure for reading */

  wpf = setup_input_w_prop_file(filename);
  wph = wpf->header;

  /* File opened for serial reading */
  wpf->parallel = 0;
  wpf->byterevflag = 0;  /* Unused for ASCII */
  wpf->rank2rcv = NULL;  /* Unused for ASCII */

  /* Indicate coordinate natural ordering */
  wph->order = NATURAL_ORDER;

  /* Node 0 alone opens a file and reads the header */

  if(this_node==0)
    {
      fp = fopen(filename,"r");
      if(fp==NULL)
	{
	  printf("r_ascii_w_i: Node %d can't open file %s, error %d\n",
		 this_node,filename,errno);fflush(stdout);terminate(1);
        }
      wpf->fp = fp;

      if( (fscanf(fp,"%d",&wph->magic_number))!=1 )
	{
	  printf("r_ascii_w_i: Error in reading version number\n"); terminate(1);
	}
      if(wph->magic_number != W_PROP_VERSION_NUMBER)
	{
	  printf("r_ascii_w_i: Unrecognized magic number in propagator file header.\n");
	  printf("Expected %d but read %d\n",
		     W_PROP_VERSION_NUMBER,wph->magic_number);
	  terminate(1);
	}
      if(fscanf(fp,"%*[ \f\n\r\t\v]%*[\"]%[^\"]%*[\"]",
		     wph->time_stamp)!=1)
	{
	  printf("r_ascii_w_i: Error reading time stamp\n"); terminate(1);
	}
      if( (fscanf(fp,"%d%d%d%d",&wph->dims[0],&wph->dims[1],
		  &wph->dims[2],&wph->dims[3]))!=4 )
	{
	  printf("r_ascii_w_i: Error reading lattice dimensions\n"); terminate(1);
	}
      if( wph->dims[0]!=nx || wph->dims[1]!=ny 
	 || wph->dims[2]!=nz || wph->dims[3]!=nt )
	{
	  /* So we can use this routine to discover the dimensions,
	     we provide that if nx = ny = nz = nt = -1 initially
	     we don't die */
	  if(nx != -1 || ny != -1 || nz != -1 || nt != -1)
	    {
	      printf("r_ascii_w_i: Incorrect lattice size %d,%d,%d,%d\n",wph->dims[0],
		     wph->dims[1],wph->dims[2],wph->dims[3]);terminate(1);
	    }
	  else
	    {
	      nx = wph->dims[0];
	      ny = wph->dims[1];
	      nz = wph->dims[2];
	      nt = wph->dims[3];
	      volume = nx*ny*nz*nt;
	    }
	}
      wph->header_bytes = 0;    /* Unused for ASCII */
    }

  else wpf->fp = NULL;  /* Other nodes don't know about this file */

  /* Broadcasts the header structure from node 0 to all nodes */
  
  broadcast_bytes((char *)wph,sizeof(w_prop_header));

  return wpf;
} /* r_ascii_w_i */

/*---------------------------------------------------------------------------*/
/* Read a propagator for specified source spin and color from an ASCII file */
/* Actually spin and color are ignored.  The next Wilson vector in the
   file is read regardless.  The file must be written and read in the
   same order */

int r_ascii_w(w_prop_file *wpf,int spin,int color,field_offset src)
{
  /* 0 normal exit code
     1 read error */

  w_prop_header *wph;
  FILE *fp;
  int destnode;
  int i,j,x,y,z,t;
  wilson_vector lbuf;
  int status;

  wph = wpf->header;
  fp = wpf->fp;

  g_sync();

  status = 0;
  if(this_node == 0)
    {
      if( (fscanf(fp,"%d%d",&i,&j))!=2 )
	{
	  printf("r_ascii_w: Error reading spin-color\n");
	  status = 1;
	}
      if(status == 0 && (i != spin || j != color))
	{
	  printf("r_ascii_w: file spin = %d, file color=%d, prog spin=%d, prog color=%d\n",
		 i,j,spin,color);
	  status = 1;
	}
    }
  broadcast_bytes((char *)&status,sizeof(int));
  if(status != 0)return status;

  status = 0;
  for(t=0;t<nt;t++)for(z=0;z<nz;z++)for(y=0;y<ny;y++)for(x=0;x<nx;x++)
    {
      destnode=node_number(x,y,z,t);

      /* Node 0 reads, and sends site to correct node */
      if(this_node==0)
	{
	  for(i=0;i<4;i++)for(j=0;j<3;j++)
	    {
	      if( 
#if PRECISION == 1
		 fscanf(fp,"%e%e\n",
#else
		 fscanf(fp,"%le%le\n",
#endif
			&(lbuf.d[i].c[j].real),
			&(lbuf.d[i].c[j].imag) )!= 2)
		{
		  if(status == 0)
		    printf("r_ascii_w: Error reading wilson_vector\n"); 
		  status = 1;
		}
	    }
	  if(destnode==0)
	    { /* just copy Wilson vector */
	      i = node_index(x,y,z,t);
	      *( (wilson_vector *)F_PT( &(lattice[i]), src ) )=lbuf;
	    }
	  else 
	    {              /* send to correct node */
	      send_field((char *)&lbuf,sizeof(wilson_vector),destnode);
	    }
	}
      
      /* The node which contains this site reads message */
      else
	{ 
	  /* for all nodes other than node 0 */
	  if(this_node==destnode)
	    {
	      get_field((char *)&lbuf,sizeof(wilson_vector),0);
	      i = node_index(x,y,z,t);
	      *( (wilson_vector *)F_PT( &(lattice[i]), src ) )=lbuf;
	    }
	}
    }

  broadcast_bytes((char *)&status,sizeof(int));

  return status;
} /* r_ascii_w */

/*---------------------------------------------------------------------------*/
/* Close propagator file */
void r_ascii_w_f(w_prop_file *wpf)
{
  FILE *fp;

  fp = wpf->fp;

  g_sync();
  if(this_node==0)
    {
/*      printf("Closed ASCII prop file  %s\n",wpf->filename);*/
      fclose(fp);
      wpf->fp = NULL;
      fflush(stdout);
    }
} /* r_ascii_w_f */
/*---------------------------------------------------------------------------*/
/* Construct a unique dump filename for this node from a given stem */
char *make_multidump_filename(char *filename)
{
  char *multidump_filename;

  multidump_filename = (char *)malloc(MAXDUMPFILENAME*sizeof(char));

  sprintf(multidump_filename,"%s.n%04d",filename,this_node);

  return multidump_filename;

} /* make_multidump_filename */

/*---------------------------------------------------------------------------*/
/* Opens separate file for each node for writing a no-frills dump file */

w_prop_file *w_multidump_w_i(char *filename)
{
  w_prop_file *wpf;
  FILE *fp;
  char *multidump_filename;

  /* Allocate space for a new file structure */

  wpf = (w_prop_file *)malloc(sizeof(w_prop_file));
  if(wpf == NULL)
    {
      printf("w_multidump_w_i: Can't malloc wpf\n");
      terminate(1);
    }

  /* Create file name */

  multidump_filename = make_multidump_filename(filename);
  
  /* All nodes open a separate serial file */

  fp = fopen(multidump_filename, "wb");

  if(fp == NULL)
    {
      printf("w_multidump_w_i: Node %d can't open file %s, error %d\n",
	     this_node,filename,errno);fflush(stdout);terminate(1);
    }

  /* Assign values to file structure */

  wpf->fp             = fp;
  wpf->filename       = multidump_filename;
  wpf->byterevflag    = 0;            /* Not used for writing */
  wpf->parallel       = 0;            /* File opened in serial */

  /** if(this_node==0)
    printf("Opened prop file %s for multidump writing\n",wpf->filename);**/

  return wpf;

} /* w_multidump_w_i */

/*---------------------------------------------------------------------------*/
/* Dumps wilson vector to a no-frills dump file for this spin and color */
/* Since these files are certainly not archived, but more likely used
   as temporary extended storage, the file is in generic precision */

void w_multidump_w(w_prop_file *wpf, int spin, int color, 
		   field_offset src_site, wilson_vector *src_field)
{

  FILE *fp;
  wilson_vector *lbuf;
  wilson_vector *src;
  u_int32type *val;
  int k;
  int rank29,rank31;
  int buf_length;
  off_t w_prop_size,w_prop_check_size,body_size,offset;
  register site *s;
  register int i;
  char myname[] = "w_multidump_w";

  lbuf = (wilson_vector *)malloc(MAX_BUF_LENGTH*sizeof(wilson_vector));
  if(lbuf == NULL)
    {
      printf("%s: Node %d can't malloc lbuf\n",myname,this_node);
      fflush(stdout);
      terminate(1);
    }

  fp = wpf->fp;

  /* initialize checksums */
  wpf->check.sum31 = 0;
  wpf->check.sum29 = 0;
  /* counts 32-bit words mod 29 and mod 31 in order of appearance on file */
  /* Here all nodes use these values */
  rank29 = sizeof(wilson_vector)/sizeof(int32type)*sites_on_node*this_node % 29;
  rank31 = sizeof(wilson_vector)/sizeof(int32type)*sites_on_node*this_node % 31;

  /* Position file for writing */

  w_prop_size = sites_on_node*sizeof(wilson_vector) ;
  w_prop_check_size =  sizeof(wpf->check.spin) +
    sizeof(wpf->check.color) +  sizeof(wpf->check.sum29) +
      sizeof(wpf->check.sum31);
  body_size = w_prop_size + w_prop_check_size;
  offset = body_size*(spin*3 + color);

  if( fseeko(fp,offset,SEEK_SET) < 0 ) 
    {
      printf("%s: Node %d fseeko %lld failed error %d file %s\n",
	     myname,this_node,(long long)offset,errno,wpf->filename);
      fflush(stdout);
      terminate(1);
    }

  /* Dump node data */

  buf_length = 0;

  FORALLSITES(i,s)
  {
        
    /* load the quark propagator into the buffer */
    if(src_site == (field_offset)(-1))
      src = src_field + i;
    else
      src = (wilson_vector *)F_PT( s, src_site );

    /* Just write as is.  No precision conversion. */
    lbuf[buf_length]= *src;

    /* Accumulate checksums - contribution from next site moved into buffer*/
    for(k = 0, val = (u_int32type *)&lbuf[buf_length]; 
	k < (int)sizeof(wilson_vector)/(int)sizeof(int32type); k++, val++)
      {
	wpf->check.sum29 ^= (*val)<<rank29 | (*val)>>(32-rank29);
	wpf->check.sum31 ^= (*val)<<rank31 | (*val)>>(32-rank31);
	rank29++; if(rank29 >= 29)rank29 = 0;
	rank31++; if(rank31 >= 31)rank31 = 0;
      }

    buf_length++;
    
    if( (buf_length == MAX_BUF_LENGTH) || (i == sites_on_node -1))
      {
	/* write out buffer */
	
	if( (int)fwrite(lbuf,sizeof(wilson_vector),buf_length,fp) 
	    != buf_length)
	  {
	    printf("%s: Node %d propagator write error %d file %s\n",
		   myname,this_node,errno,wpf->filename); 
	    fflush(stdout);
	    terminate(1);   
	  }
	buf_length = 0;		/* start again after write */
      }
    
  } 
  
  free(lbuf);

  /* Combine checksums */

  g_xor32(&wpf->check.sum29);
  g_xor32(&wpf->check.sum31);
  wpf->check.spin = spin;
  wpf->check.color = color;

  /* Write checksum at end of lattice file */

  /* Serial, because file is private to this node */
  write_checksum_w(SERIAL,wpf);
      
  /* printf("Dumped prop for spin %d color %d to checkpoint file %s\n",
	     spin,color,wpf->filename);fflush(stdout);*/

} /* w_multidump_w */

/*---------------------------------------------------------------------------*/
/* Dumps wilson vector to a no-frills dump file for this spin and color */
/* Since these files are certainly not archived, but more likely used
   as temporary extended storage, the file is in generic precision */

void w_multidump_w_from_site(w_prop_file *wpf, int spin, int color, 
		   field_offset src_site)
{

  w_multidump_w(wpf, spin, color, src_site, NULL);
}

/*---------------------------------------------------------------------------*/
/* Dumps wilson vector to a no-frills dump file for this spin and color */
/* Since these files are certainly not archived, but more likely used
   as temporary extended storage, the file is in generic precision */

void w_multidump_w_from_field(w_prop_file *wpf, int spin, int color, 
			      wilson_vector *src_field)
{

  w_multidump_w(wpf, spin, color, (field_offset)(-1), src_field);
}


/*---------------------------------------------------------------------------*/

void w_multidump_w_o(w_prop_file *wpf)
{
  /* Reopen previously opened file, specified by wpf, for multidump writing */

  wpf->fp = fopen(wpf->filename, "rb+");  /* Binary write and update */

  if(wpf->fp == NULL)
    {
      printf("w_multidump_w_o: Node %d can't reopen file %s, error %d\n",
	     this_node,wpf->filename,errno);fflush(stdout);
      terminate(1);
    }

  /*if(this_node==0)printf("Reopened prop file %s for multidump writing.\n",
			 wpf->filename);*/

} /* w_multidump_w_o */

/*---------------------------------------------------------------------------*/

void w_multidump_w_c(w_prop_file *wpf)
{
  /* Close file temporarily -- releases only the system FILE structure */
  g_sync();
  if(wpf->fp != NULL)fclose(wpf->fp);
  wpf->fp = NULL;
  /*if(this_node==0)printf("Closed dump file %s temporarily\n",wpf->filename);*/
} /* w_multidump_w_c */

/*---------------------------------------------------------------------------*/

void w_multidump_w_f(w_prop_file *wpf)
{
  /* Close file (if still active) and release header and file structures */
  if(wpf->fp != NULL)w_multidump_w_c(wpf);
  if(this_node==0)
    printf("Wrote multidump file %s\n",wpf->filename);
  free(wpf->filename);
  free(wpf);
} /* w_multidump_w_f */

/*---------------------------------------------------------------------------*/

w_prop_file *r_multidump_w_i(char *filename)
{
  /* Returns file descriptor for opened file */

  w_prop_file *wpf;
  FILE *fp;
  char *multidump_filename;

  /* Allocate space for a new file structure */

  wpf = (w_prop_file *)malloc(sizeof(w_prop_file));
  if(wpf == NULL)
    {
      printf("r_multidump_w_i: Can't malloc wpf\n");
      terminate(1);
    }

  /* Create file name */

  multidump_filename = make_multidump_filename(filename);
  
  /* All nodes open a separate serial file */

  fp = fopen(multidump_filename, "rb");

  if(fp == NULL)
    {
      printf("r_multidump_w_i: Node %d can't open file %s, error %d\n",
	     this_node,filename,errno);fflush(stdout);terminate(1);
    }

  /* Assign values to file structure */

  wpf->fp             = fp;
  wpf->filename       = multidump_filename;
  wpf->byterevflag    = 0;            /* Not used for dump format */
  wpf->parallel       = 0;            /* File opened in serial */

  return wpf;

} /* r_multidump_w_i */

/*----------------------------------------------------------------------*/

/* Read Wilson propagator from a separate multidump file for each node */
/* Spin, color combinations must be read in the same order as written */
/* WARNING: The file is assumed to have been written in the same
   precision as read.  */

int r_multidump_w(w_prop_file *wpf, int spin, int color, 
		  field_offset dest_site, wilson_vector *dest_field)
{
  /* wpf  = propagator file structure 
     dest_site  = field offset for propagator Wilson vector
     dest_field = pointer to field
     use only one of the dest parameters!  */

  /* 0 is normal exit code
     1 for seek, read error, or missing data error */

  FILE *fp;
  char *filename;
  int byterevflag;
  wilson_vector *lbuf;
  wilson_vector *dest;

  int buf_length,where_in_buf;
  int w_prop_size,w_prop_check_size,body_size,offset;
  w_prop_check test_wpc;
  u_int32type *val;
  int rank29,rank31;
  register int i,k;
  register site *s;
  int status;
  Real xstatus;

  char myname[] = "r_multidump_w";

  fp = wpf->fp;
  filename = wpf->filename;
  byterevflag = wpf->byterevflag;

  /* Allocate read buffer */
  lbuf = (wilson_vector *)malloc(MAX_BUF_LENGTH*sizeof(wilson_vector));
  if(lbuf == NULL)
    {
      printf("%s: Node %d can't malloc lbuf\n",myname,this_node);
      fflush(stdout);
      terminate(1);
    }

  /* Position the file for reading */

  w_prop_size = sites_on_node*sizeof(wilson_vector) ;
  w_prop_check_size =  sizeof(wpf->check.spin) +
    sizeof(wpf->check.color) +  sizeof(wpf->check.sum29) +
      sizeof(wpf->check.sum31);
  body_size = w_prop_size + w_prop_check_size;
  offset = body_size*(spin*3 + color);

  xstatus = 0.;
  if( fseeko(fp,offset,SEEK_SET) < 0 ) 
    {
      printf("%s: Node %d fseeko %lld failed error %d file %s\n",
	     myname,this_node,(long long)offset,errno,filename);
      fflush(stdout);
      xstatus = 1.;
    }
  g_floatsum(&xstatus);
  if(xstatus != 0.)return 1;

  /* Each node reads its propagator values */

  /* initialize checksums */
  test_wpc.sum31 = 0;
  test_wpc.sum29 = 0;
  /* counts 32-bit words mod 29 and mod 31 in order of appearance on file */
  /* Here all nodes use these values */
  rank29 = sizeof(wilson_vector)/sizeof(int32type)*sites_on_node*this_node % 29;
  rank31 = sizeof(wilson_vector)/sizeof(int32type)*sites_on_node*this_node % 31;

  buf_length = 0;
  where_in_buf = 0;
  
  FORALLSITES(i,s)
  {
    if(where_in_buf == buf_length)
      
      {  /* get new buffer */
	
	/* new buffer length  = remaining sites, but never bigger 
	   than MAX_BUF_LENGTH */
	buf_length = sites_on_node - i;
	if(buf_length > MAX_BUF_LENGTH) buf_length = MAX_BUF_LENGTH; 
	/* then do read */
	/* each node reads its sites */
	
	if( (int)fread(lbuf,buf_length*sizeof(wilson_vector),1,fp) != 1)
	  {
	    if(xstatus == 0.)
	      printf("%s: node %d propagator read error %d file %s\n",
		     myname,this_node,errno,filename); 
	    fflush(stdout); 
	    xstatus = 1.;
	  }
	where_in_buf = 0;  /* reset counter */
      }  /*** end of the buffer read ****/
    
    /* Accumulate checksums - contribution from next site */
    for(k = 0, val = (u_int32type *)&lbuf[where_in_buf]; 
	k < (int)sizeof(wilson_vector)/(int)sizeof(int32type); k++, val++)
      {
	test_wpc.sum29 ^= (*val)<<rank29 | (*val)>>(32-rank29);
	test_wpc.sum31 ^= (*val)<<rank31 | (*val)>>(32-rank31);
	rank29++; if(rank29 >= 29)rank29 = 0;
	rank31++; if(rank31 >= 31)rank31 = 0;
      }

    if(dest_site == (field_offset)(-1))
      dest = dest_field + i;
    else
      dest = (wilson_vector *)F_PT( &(lattice[i]), dest_site );

    /* No precision conversion.  Data is written as is */
    *dest = lbuf[where_in_buf];

    where_in_buf++;
  } /** i,s **/
  
  g_floatsum(&xstatus);
  if(xstatus != 0)return 1;

  free(lbuf);

  /* Combine node checksum contributions with global exclusive or */
  g_xor32(&test_wpc.sum29);
  g_xor32(&test_wpc.sum31);

  /* Read spin, color, and checksum record */

  status = 0;

  status += sread_data(fp,&wpf->check.spin,
			    sizeof(wpf->check.spin),myname,"check.spin");
  status += sread_data(fp,&wpf->check.color,
			    sizeof(wpf->check.color),myname,"check.color");
  status += sread_data(fp,&wpf->check.sum29,
			    sizeof(wpf->check.sum29),myname,"check.sum29");
  status += sread_data(fp,&wpf->check.sum31,
			    sizeof(wpf->check.sum31),myname,"check.sum31");

  /* Verify spin and color */
  
  if(wpf->check.spin != spin || wpf->check.color != color)
    {
      printf("%s: Spin %d and color %d do not match check record on file %s\n",
	     myname,spin,color,filename);
      printf("  Check record said %d %d\n",wpf->check.spin,wpf->check.color);
      fflush(stdout);
      status = 1;
    }
  
  broadcast_bytes((char *)&status,sizeof(int));
  if(status != 0)return status;
  
  /* Verify checksum */
  
  if(wpf->check.sum29 != test_wpc.sum29 ||
     wpf->check.sum31 != test_wpc.sum31)
    {
      printf("%s: Checksum violation spin %d color %d file %s\n",
	     myname,wpf->check.spin,wpf->check.color,wpf->filename);
      printf("Computed %x %x.  Read %x %x.\n",
	     test_wpc.sum29,test_wpc.sum31,
	     wpf->check.sum29,wpf->check.sum31);
    }
  else
    printf("Checksums %x %x OK for spin %d color %d file %s\n",
	   wpf->check.sum29,wpf->check.sum31,
	   wpf->check.spin,wpf->check.color,wpf->filename);
  
  /*if(this_node == 0)
    {
      printf("Read prop from multidump for spin %d color %d in from file %s\n",
	     spin,color,filename);
    }*/
  
  fflush(stdout);
  
  return 0;
  
} /* r_multidump_w */


/*----------------------------------------------------------------------*/

/* Read Wilson propagator from a separate multidump file for each node */
/* Spin, color combinations must be read in the same order as written */
/* WARNING: The file is assumed to have been written in the same
   precision as read.  */

int r_multidump_w_to_site(w_prop_file *wpf, int spin, int color, 
			  field_offset dest_site)
{
  /* wpf  = propagator file structure 
     dest_site  = field offset for propagator Wilson vector */

  /* 0 is normal exit code
     1 for seek, read error, or missing data error */

  return r_multidump_w(wpf, spin, color, dest_site, NULL);
}

/*----------------------------------------------------------------------*/

/* Read Wilson propagator from a separate multidump file for each node */
/* Spin, color combinations must be read in the same order as written */
/* WARNING: The file is assumed to have been written in the same
   precision as read.  */

int r_multidump_w_to_field(w_prop_file *wpf, int spin, int color, 
			   wilson_vector *dest_field)
{
  /* wpf  = propagator file structure 
     dest_field = pointer to field */

  /* 0 is normal exit code
     1 for seek, read error, or missing data error */

  return r_multidump_w(wpf, spin, color, (field_offset)(-1), dest_field);
}

/*----------------------------------------------------------------------*/

void r_multidump_w_o(w_prop_file *wpf)
{
  /* Reopen previously opened file, specified by wpf, for multidump reading */
  
  wpf->fp = fopen(wpf->filename, "rb");  /* Binary read */
  
  if(wpf->fp == NULL)
    {
      printf("r_multidump_w_o: Node %d can't reopen file %s, error %d\n",
	     this_node,wpf->filename,errno);fflush(stdout);
      terminate(1);
    }
  
  /*if(this_node==0)printf("Reopened prop file %s for multidump reading\n",
			 wpf->filename);*/
  wpf->parallel = SERIAL;
  
} /* r_multidump_w_o */

/*---------------------------------------------------------------------------*/

void r_multidump_w_c(w_prop_file *wpf)
{
  /* Close dump file temporarily -- releases only the system FILE structure */
  g_sync();
  if(wpf->fp != NULL)fclose(wpf->fp);
  wpf->fp = NULL;
  /*if(this_node==0)printf("Closed dump prop file %s temporarily\n",wpf->filename);*/
} /* r_multidump_w_c */

/*---------------------------------------------------------------------------*/

void r_multidump_w_f(w_prop_file *wpf)
{
  /* Close dump file (if active) and release header and file structures */

  if(wpf->fp != NULL)r_multidump_w_c(wpf);
  free(wpf->filename);
  free(wpf);
  /*if(this_node==0)printf("Released dump prop file %s\n",wpf->filename);*/
} /* r_multidump_w_f */

/*---------------------------------------------------------------------------*/
