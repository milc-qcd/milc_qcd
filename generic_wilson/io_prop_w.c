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
#include "../include/io_scidac_w.h"
#ifdef HAVE_QIO
#include <qio.h>
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
  /* Make exception for FNAL format timeslice checksum 
     since there are too many to list */
  if(strstr(keyword,"quark.t") == NULL &&
     strlen(w_prop_info_keyword[i])==0)
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
      else if(strstr(fmt,"lu") != NULL){
	snprintf(string+bytes,nstring-bytes,fmt,*(unsigned long *)data);
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
  strcat(info_filename,ASCII_INFO_EXT);

  /* Open header file */
  
  if((info_fp = g_open(info_filename,"w")) == NULL)
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

  g_close(info_fp);

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
  wpf->fp = NULL;

  /* Allocate space for the header */

  /* Make sure compilation gave us a 32 bit integer type */
  assert(sizeof(int32type) == 4);

  wph = (w_prop_header *)malloc(sizeof(w_prop_header));
  if(wph == NULL)
    {
      printf("setup_input_w_prop_file: Can't malloc wph\n");
      terminate(1);
    }

  wpf->header         = wph;
  wpf->prop           = NULL;
  wpf->byterevflag    = 0;
  wpf->parallel       = 0;
  wpf->info_fp        = NULL;
  wpf->file_type      = FILE_TYPE_UNKNOWN;
#ifdef HAVE_QIO
  wpf->infile         = NULL;
#endif

  return wpf;
} /* setup_input_w_prop_file */

/*----------------------------------------------------------------------*/

void clear_input_w_prop_file(w_prop_file *wpf){
  /* Free structures associated with w_prop_file */
      
  if(wpf->fp       != NULL)g_close(wpf->fp);
  if(wpf->header   != NULL)free(wpf->header);
  if(wpf->prop     != NULL)free(wpf->prop);
#ifdef HAVE_QIO
  if(wpf->infile   != NULL)QIO_close_read(wpf->infile);
#endif
  free(wpf);
}

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

  wph->magic_number = 0;  /* We don't support the old MILC formats */

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
      for(i = strlen(wph->time_stamp) + 1; i < (int)sizeof(wph->time_stamp); i++)
	wph->time_stamp[i] = '\0';
      
      /* Remove trailing end-of-line character */
      if(wph->time_stamp[strlen(wph->time_stamp) - 1] == '\n')
	wph->time_stamp[strlen(wph->time_stamp) - 1] = '\0';
    }
  
  /* Broadcast to all nodes */
  broadcast_bytes(wph->time_stamp,sizeof(wph->time_stamp));

  build_w_prop_hdr(wph);
  
  wpf->fp             = NULL;
  wpf->prop           = NULL;
  wpf->file_type      = FILE_TYPE_UNKNOWN;

  return wpf;
} /* setup_output_w_prop_file */

/*----------------------------------------------------------------------*/

void clear_output_w_prop_file(w_prop_file *wpf){
  /* Free structures associated with w_prop_file */
      
  if(wpf->fp       != NULL)g_close(wpf->fp);
  if(wpf->header   != NULL)free(wpf->header);
  if(wpf->prop     != NULL)free(wpf->prop);
#ifdef HAVE_QIO
  if(wpf->outfile != NULL)
    w_close_usqcd_wprop_file(wpf->outfile);
#endif
  free(wpf);
}

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

    fp = g_open(filename,"w");
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
  wpf->byterevflag    = 0;            /* Not used for writing */

  /* Node 0 writes info file */
  if(this_node==0)write_w_prop_info_file(wpf);

  return wpf;
} /* w_ascii_w_i */
  
/*---------------------------------------------------------------------------*/
/* Write ASCII propagator */

void w_ascii_w(w_prop_file *wpf,int spin,int color, wilson_vector *src)
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
	      lbuf = *( src + l );
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
	      lbuf = *( src + l );
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
      if(wpf->fp != NULL)g_close(wpf->fp);
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

  /* Indicate coordinate natural ordering */
  wph->order = NATURAL_ORDER;

  /* Node 0 alone opens a file and reads the header */

  if(this_node==0)
    {
      fp = g_open(filename,"r");
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

  FILE *fp;
  int destnode;
  int i,j,x,y,z,t;
  wilson_vector lbuf;
  int status;

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
      g_close(fp);
      wpf->fp = NULL;
      fflush(stdout);
    }
} /* r_ascii_w_f */
