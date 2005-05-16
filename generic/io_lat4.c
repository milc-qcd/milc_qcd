/*********************** io_lat4.c *************************/
/* MIMD version 6 */

/* routines for gauge configuration input/output. */
/* This works for most machines.  Wrappers for parallel I/O
   are in io_ansi.c, io_piofs.c, or io_paragon2.c */

/* Modifications */
/* 10/04/01 Removed save_old_binary (but can still read old binary) C.D. */
/* 7/11/01 large file (64 bit addressing) support */
/* 4/16/00 additions to READ ARChive format J.H. */
/*         adapted for version 6 12/21/00 UMH */
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

#ifndef HAVE_FSEEKO
#define fseeko fseek
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


/*----------------------------------------------------------------------
   Routines for archive I/O
   -----------------------------------------------------------------------*/

int qcdhdr_get_str(char *s, QCDheader *hdr, char **q) {     
  /* find a token and return the value */
  int i;
  for (i=0; i<(char)(*hdr).ntoken; i++) {
    if (strcmp(s,(char *)(*hdr).token[i])==0) {
      *q = (*hdr).value[i];
      return(SUCCESS);
    }
  }
  *q = NULL;
  return (FAILURE);
}
  
int qcdhdr_get_int(char *s,QCDheader *hdr,int *q) {
  char *p;
  qcdhdr_get_str(s,hdr,&p);
  if (p==NULL) return (FAILURE);
  sscanf(p,"%d",q);
  return (SUCCESS);
}
int qcdhdr_get_int32x(char *s,QCDheader *hdr,u_int32type *q) {
  char *p;
  int r;
  qcdhdr_get_str(s,hdr,&p);
  if (p==NULL) return (FAILURE);
  sscanf(p,"%x",&r);
  *q = r;
  return (SUCCESS);
}
int qcdhdr_get_float(char *s, QCDheader *hdr, Real *q) {
  char *p;
  qcdhdr_get_str(s,hdr,&p);
  if (p==NULL) return (FAILURE);
#if PRECISION == 1
  sscanf(p,"%f",q);
#else
  sscanf(p,"%lf",q);
#endif
  return (SUCCESS);
}

void error_exit(char *s) { printf("%s\n",s); terminate(1);}

void complete_U(float *u) {
  u[12] = u[ 2]*u[10] - u[ 4]*u[ 8] - u[ 3]*u[11] + u[ 5]*u[ 9];
  u[13] = u[ 4]*u[ 9] - u[ 2]*u[11] + u[ 5]*u[ 8] - u[ 3]*u[10];
  u[14] = u[ 4]*u[ 6] - u[ 0]*u[10] - u[ 5]*u[ 7] + u[ 1]*u[11];
  u[15] = u[ 0]*u[11] - u[ 4]*u[ 7] + u[ 1]*u[10] - u[ 5]*u[ 6];
  u[16] = u[ 0]*u[ 8] - u[ 2]*u[ 6] - u[ 1]*u[ 9] + u[ 3]*u[ 7];
  u[17] = u[ 2]*u[ 7] - u[ 0]*u[ 9] + u[ 3]*u[ 6] - u[ 1]*u[ 8];
}


int big_endian() {
  union  {
    long l;
    char c[sizeof (long)];
  } u;
  u.l = 1;
  return (u.c[sizeof (long) - 1] == 1);
}

QCDheader * qcdhdr_get_hdr(FILE *in)
{
#define MAX_LINE_LENGTH 1024
#define MAX_TOKENS 512
  char line[MAX_LINE_LENGTH];
  int n,len;
  QCDheader *hdr;
  char **tokens, **values;
  char *p, *q;

  /* Begin reading, and check for "BEGIN_HEADER" token */
  fgets(line,MAX_LINE_LENGTH,in);
  /*
  if (strcmp(line,"BEGIN_HEADER\n")!=0)
    error_exit("qcdhdr_get_hdr: Missing \"BEGIN_HEADER\"; punting \n");
  */
  /* Allocate space for QCDheader and its pointers */
  tokens = (char **) malloc(MAX_TOKENS*sizeof(char *));
  values = (char **) malloc(MAX_TOKENS*sizeof(char *));
  hdr = (QCDheader *) malloc(sizeof(QCDheader));
  (*hdr).token = tokens;
  (*hdr).value = values;

  /* Begin loop on tokens */
  n = 0;
  printf("reading Archive header:\n");
  while (1) {
    fgets(line,MAX_LINE_LENGTH,in);
    printf("%s", line);

    if (strcmp(line,"END_HEADER\n")==0) break;

    /* Tokens are terminated by a space */
    q = strchr(line, (int)' ');

    /* Overwrite space with a terminating null */
    *q = '\0';
    len = (int) q - (int) line;

    /* allocate space and copy the token in to it */
    p = (char *)malloc(len+1);
    (*hdr).token[n] = p;
    strcpy(p,line);

    q = strchr(++q, (int)'='); q++;
    len = strlen(q);
    q[len-1] = 0;
    p = (char *)malloc(len);
    (*hdr).value[n] = p;
    strcpy(p,q);
    n++;
  }
  (*hdr).ntoken = n;
  return (hdr);
}

/* Destroy header - for freeing up storage */
void qcdhdr_destroy_hdr(QCDheader *hdr){
  int i;
  
  if(hdr == NULL)return;

  for(i = 0; i < hdr->ntoken; i++){
    free(hdr->value[i]);
    free(hdr->token[i]);
  }

  free(hdr->token);
  free(hdr->value);
  free(hdr);
}

/*---------------------------------------------------------------------------*/
/* Convert (or copy) four single precision su3_matrices to generic precision */

void f2d_4mat(fsu3_matrix *a, su3_matrix *b){
  int dir,i,j;
  
  for(dir = 0; dir < 4; dir++){
    for(i = 0; i < 3; i++)for(j = 0; j < 3; j++){
      b[dir].e[i][j].real = a[dir].e[i][j].real;
      b[dir].e[i][j].imag = a[dir].e[i][j].imag;
    }
  }
}

/* Convert (or copy) four generic precision su3_matrices to single precision */
void d2f_4mat(su3_matrix *a, fsu3_matrix *b){
  int dir,i,j;
  
  for(dir = 0; dir < 4; dir++){
    for(i = 0; i < 3; i++)for(j = 0; j < 3; j++){
      b[dir].e[i][j].real = a[dir].e[i][j].real;
      b[dir].e[i][j].imag = a[dir].e[i][j].imag;
    }
  }
}

/*---------------------------------------------------------------------------*/
void swrite_data(FILE* fp, void *src, size_t size, char *myname, char *descrip)
{
  if(fwrite(src,size,1,fp) != 1)
    {
      printf("%s: Node %d %s write error %d\n",
	    myname,this_node,descrip,errno);
      fflush(stdout);
      terminate(1);
    }
}
/*---------------------------------------------------------------------------*/
void pwrite_data(FILE* fp, void *src, size_t size, char *myname, char *descrip)
{
  if(g_write(src,size,1,fp) != 1)
    {
      printf("%s: Node %d %s descrip,write error %d\n",
	    myname,this_node,descrip,errno);
      fflush(stdout);
      terminate(1);
    }
}
/*---------------------------------------------------------------------------*/
void pswrite_data(int parallel, FILE* fp, void *src, size_t size, 
		 char *myname, char *descrip)
{
  if(parallel)pwrite_data(fp,src,size,myname,descrip);
  else        swrite_data(fp,src,size,myname,descrip);
}
/*---------------------------------------------------------------------------*/
int sread_data(FILE* fp, void *src, size_t size, char *myname, char *descrip)
{
  if(fread(src,size,1,fp) != 1)
    {
      printf("%s: Node %d %s read error %d\n",
	    myname,this_node,descrip,errno);
      fflush(stdout);
      return 1;
    }
  return 0;
}
/*---------------------------------------------------------------------------*/
int pread_data(FILE* fp, void *src, size_t size, char *myname, char *descrip)
{
  if(g_read(src,size,1,fp) != 1)
    {
      printf("%s: Node %d %s read error %d\n",
	    myname,this_node,descrip,errno);
      fflush(stdout);
      return 1;
    }
  return 0;
}
/*---------------------------------------------------------------------------*/
int pread_byteorder(int byterevflag, FILE* fp, void *src, size_t size, char *myname, char *descrip)
{
  int status;

  status = pread_data(fp,src,size,myname,descrip);
  if(byterevflag==1)
    byterevn((int32type *)src,size/sizeof(int32type));
  return status;
}
/*---------------------------------------------------------------------------*/
int sread_byteorder(int byterevflag, FILE* fp, void *src, size_t size, char *myname, char *descrip)
{
  int status;

  status = sread_data(fp,src,size,myname,descrip);
  if(byterevflag==1)
    byterevn((int32type *)src,size/sizeof(int32type));
  return status;
}
/*---------------------------------------------------------------------------*/
int psread_data(int parallel, FILE* fp, void *src, size_t size, 
		 char *myname, char *descrip)
{
  if(parallel)return pread_data(fp,src,size,myname,descrip);
  else        return sread_data(fp,src,size,myname,descrip);
}
/*---------------------------------------------------------------------------*/
int psread_byteorder(int byterevflag, int parallel, FILE* fp, 
		      void *src, size_t size, 
		      char *myname, char *descrip)
{
  if(parallel)return pread_byteorder(byterevflag,fp,src,size,myname,descrip);
  else        return sread_byteorder(byterevflag,fp,src,size,myname,descrip);
}
/*---------------------------------------------------------------------------*/

/* This subroutine writes the gauge configuration header structure */
/* Parallel access version */
/* While the procedures for serial and parallel writing are
   identical, (the header is written only by node 0, no matter what),
   the file which is accessed can be opened either by all
   nodes in w_parallel or one node in w_serial.  We have to
   distinguish between these modes when writing */

void pwrite_gauge_hdr(FILE *fp, gauge_header *gh)
{

  char myname[] = "pwrite_gauge_hdr";

  pwrite_data(fp,(void *)&gh->magic_number,sizeof(gh->magic_number),
	      myname,"magic_number");
  pwrite_data(fp,(void *)gh->dims,sizeof(gh->dims),
	      myname,"dimensions");
  pwrite_data(fp,(void *)gh->time_stamp,sizeof(gh->time_stamp),
	      myname,"time_stamp");
  pwrite_data(fp,&gh->order,sizeof(gh->order),
	      myname,"order");

  /* Header byte length */

  gh->header_bytes = sizeof(gh->magic_number) + sizeof(gh->dims) + 
    sizeof(gh->time_stamp) + sizeof(gh->order);

} /* pwrite_gauge_hdr */

/*----------------------------------------------------------------------*/
/* This subroutine writes the gauge configuration header structure */
/* Serial access version */

void swrite_gauge_hdr(FILE *fp, gauge_header *gh)
{

  char myname[] = "swrite_gauge_hdr";

  swrite_data(fp,(void *)&gh->magic_number,sizeof(gh->magic_number),
	      myname,"magic_number");
  swrite_data(fp,(void *)gh->dims,sizeof(gh->dims),
	      myname,"dimensions");
  swrite_data(fp,(void *)gh->time_stamp,sizeof(gh->time_stamp),
	      myname,"time_stamp");
  swrite_data(fp,&gh->order,sizeof(gh->order),
	      myname,"order");

  /* Header byte length */

  gh->header_bytes = sizeof(gh->magic_number) + sizeof(gh->dims) + 
    sizeof(gh->time_stamp) + sizeof(gh->order);
  
} /* swrite_gauge_hdr */

/*------------------------------------------------------------------------*/

/* Write a data item to the gauge info file */
int write_gauge_info_item( FILE *fpout,    /* ascii file pointer */
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
  char *data;
  float tt;

  /* Check for valid keyword */

  for(i=0;strlen(gauge_info_keyword[i])>0 &&
      strcmp(gauge_info_keyword[i],keyword) != 0; i++);
  if(strlen(gauge_info_keyword[i])==0)
    printf("write_gauge_info_item: WARNING: keyword %s not in table\n",
	    keyword);

  /* Write keyword */

  fprintf(fpout,"%s =",keyword);

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
	  printf("write_gauge_info_item: Unrecognized data type %s\n",fmt);
	  return 1;
	}
    }
  fprintf(fpout,"\n");
  return 0;
}

/*------------------------------------------------------------------------*/

/* Write a data item to a character string */
int sprint_gauge_info_item( 
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
  int bytes;
  char *data;
  float tt;

  /* Check for valid keyword */

  for(i=0;strlen(gauge_info_keyword[i])>0 &&
      strcmp(gauge_info_keyword[i],keyword) != 0; i++);
  if(strlen(gauge_info_keyword[i])==0)
    printf("write_gauge_info_item: WARNING: keyword %s not in table\n",
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
	  printf("write_gauge_info_item: Unrecognized data type %s\n",fmt);
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

void write_gauge_info_file(gauge_file *gf)
{
  FILE *info_fp;
  gauge_header *gh;
  char info_filename[256];
  char sums[20];

  gh = gf->header;

  /* Construct header file name from lattice file name 
   by adding filename extension to lattice file name */

  strcpy(info_filename,gf->filename);
  strcat(info_filename,ASCII_GAUGE_INFO_EXT);

  /* Open header file */
  
  if((info_fp = fopen(info_filename,"w")) == NULL)
    {
      printf("write_gauge_info_file: Can't open ascii info file %s\n",info_filename);
      return;
    }
  
  /* Write required information */

  write_gauge_info_item(info_fp,"magic_number","%d",(char *)&gh->magic_number,0,0);
  write_gauge_info_item(info_fp,"time_stamp","\"%s\"",gh->time_stamp,0,0);
  sprintf(sums,"%x %x",gf->check.sum29,gf->check.sum31);
  write_gauge_info_item(info_fp,"checksums","\"%s\"",sums,0,0);
  write_gauge_info_item(info_fp,"nx","%d",(char *)&nx,0,0);
  write_gauge_info_item(info_fp,"ny","%d",(char *)&ny,0,0);
  write_gauge_info_item(info_fp,"nz","%d",(char *)&nz,0,0);
  write_gauge_info_item(info_fp,"nt","%d",(char *)&nt,0,0);

  write_appl_gauge_info(info_fp);

  fclose(info_fp);

  printf("Wrote info file %s\n",info_filename);

} /*write_gauge_info_file */

/*----------------------------------------------------------------------*/

/* Set up the input gauge file and gauge header structures */

gauge_file *setup_input_gauge_file(char *filename)
{
  char myname[] = "setup_input_gauge_file";
  gauge_file *gf;
  gauge_header *gh;

  /* Allocate space for the file structure */

  gf = (gauge_file *)malloc(sizeof(gauge_file));
  if(gf == NULL)
    {
      printf("%s: Can't malloc gf\n",myname);
      terminate(1);
    }

  gf->filename = filename;

  /* Allocate space for the header */

  /* Make sure compilation gave us a 32 bit integer type */
  assert(sizeof(int32type) == 4);

  gh = (gauge_header *)malloc(sizeof(gauge_header));
  if(gh == NULL)
    {
      printf("%s: Can't malloc gh\n",myname);
      terminate(1);
    }

  gf->header = gh;
  gf->rank2rcv = NULL;
  gf->check.sum29 = 0;
  gf->check.sum31 = 0;

  return gf;
}

/*----------------------------------------------------------------------*/

/* Set up the output gauge file an gauge header structure */

gauge_file *setup_output_gauge_file()
{
  char myname[] = "setup_output_gauge_file";
  gauge_file *gf;
  gauge_header *gh;
  time_t time_stamp;
  int i;

  /* Allocate space for a new file structure */

  gf = (gauge_file *)malloc(sizeof(gauge_file));
  if(gf == NULL)
    {
      printf("%s: Can't malloc gf\n",myname);
      terminate(1);
    }

  /* Allocate space for a new header structure */

  /* Make sure compilation gave us a 32 bit integer type */
  assert(sizeof(int32type) == 4);

  gh = (gauge_header *)malloc(sizeof(gauge_header));
  if(gh == NULL)
    {
      printf("%s: Can't malloc gh\n",myname);
      terminate(1);
    }

  /* Load header pointer and file name */
  gf->header = gh;

  /* Initialize */
  gf->check.sum29 = 0;
  gf->check.sum31 = 0;

  /* Load header values */

  gh->magic_number = GAUGE_VERSION_NUMBER;

  gh->dims[0] = nx;
  gh->dims[1] = ny;
  gh->dims[2] = nz;
  gh->dims[3] = nt;

  /* Get date and time stamp. (We use local time on node 0) */

  if(this_node==0)
    {
      time(&time_stamp);
      strcpy(gh->time_stamp,ctime(&time_stamp));
      /* For aesthetic reasons, don't leave trailing junk bytes here to be
	 written to the file */
      for(i = strlen(gh->time_stamp) + 1; i < (int)sizeof(gh->time_stamp); i++)
	gh->time_stamp[i] = '\0';
      
      /* Remove trailing end-of-line character */
      if(gh->time_stamp[strlen(gh->time_stamp) - 1] == '\n')
	gh->time_stamp[strlen(gh->time_stamp) - 1] = '\0';
    }
  
  /* Broadcast to all nodes */
  broadcast_bytes(gh->time_stamp,sizeof(gh->time_stamp));

  return gf;
} /* setup_output_gauge_file */
/*----------------------------------------------------------------------*/

/* Open a binary file for serial writing by node 0 */

gauge_file *w_serial_i(char *filename)
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
/* Read checksum and compare.  It is assumed that the file is already
   correctly positioned.

   Should be called only by one node */

void read_checksum(int parallel, gauge_file *gf, gauge_check *test_gc)
{

  char myname[] = "read_checksum";
  
  /* Read checksums with byte reversal */
  
  if(psread_byteorder(gf->byterevflag,parallel,gf->fp,
	 &gf->check.sum29,sizeof(gf->check.sum29), myname,"checksum")!=0)
    terminate(1);
  if(psread_byteorder(gf->byterevflag,parallel,gf->fp,
	 &gf->check.sum31,sizeof(gf->check.sum31), myname,"checksum")!=0)
    terminate(1);

  if(gf->check.sum29 != test_gc->sum29 ||
     gf->check.sum31 != test_gc->sum31)
    printf("%s: Checksum violation. Computed %x %x.  Read %x %x.\n",
	    myname,test_gc->sum29,test_gc->sum31,
	   gf->check.sum29,gf->check.sum31);
  else
    {
      printf("Checksums %x %x OK\n",gf->check.sum29,gf->check.sum31);
      fflush(stdout);
    }
} /* read_checksum */

/*---------------------------------------------------------------------------*/
/* Write checksum to lattice file.  It is assumed that the file
   is already correctly positioned.

   Should be called only by one node */

void write_checksum(int parallel, gauge_file *gf)
{

  char myname[] = "write_checksum";

  pswrite_data(parallel,gf->fp,
	       &gf->check.sum29,sizeof(gf->check.sum29),myname,"checksum");
  pswrite_data(parallel,gf->fp,
	       &gf->check.sum31,sizeof(gf->check.sum31),myname,"checksum");
  printf("Checksums %x %x\n",gf->check.sum29,gf->check.sum31);
}

/*---------------------------------------------------------------------------*/

/* Here only node 0 writes gauge configuration to a binary file */

void w_serial(gauge_file *gf)
{
  /* gf  = file descriptor as opened by w_serial_i */

  FILE *fp;
  gauge_header *gh;
  u_int32type *val;
  int rank29,rank31;
  fsu3_matrix *lbuf;
  fsu3_matrix tbuf[4];
  int buf_length;
  register int i,j,k;
  off_t offset;             /* File stream pointer */
  off_t coord_list_size;    /* Size of coordinate list in bytes */
  off_t head_size;          /* Size of header plus coordinate list */
  off_t checksum_offset;    /* Location of checksum */
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

      if( fseeko(fp,offset,SEEK_SET) < 0 ) 
	{
	  printf("w_serial: Node %d fseeko %lld failed error %d file %s\n",
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
	if( this_node==0 && newnode!=0 )send_field((char *)tbuf,4,newnode);
	if( this_node==newnode && newnode!=0 )get_field((char *)tbuf,4,0);
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
	      
	      if( (int)fwrite(lbuf,4*sizeof(fsu3_matrix),buf_length,fp) != buf_length)
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
      if( fseeko(fp,checksum_offset,SEEK_SET) < 0 ) 
	{
	  printf("w_serial: Node %d fseeko %lld failed error %d file %s\n",
		 this_node,(long long)checksum_offset,errno,gf->filename);
	  fflush(stdout);terminate(1);
	}
      write_checksum(SERIAL,gf);
    }
  
} /* w_serial */

/*---------------------------------------------------------------------------*/

void w_serial_f(gauge_file *gf)

/* Close the file and free associated structures */
{
  g_sync();
  if(this_node==0)
    {
      if(gf->parallel)
	printf("w_serial_f: Attempting serial close on parallel file \n");

      fclose(gf->fp);
    }

  /* Node 0 writes ascii info file */

  if(this_node == 0)write_gauge_info_file(gf);

  /* Do not free gf and gf->header so calling program can use them */

} /* w_serial_f */

/*---------------------------------------------------------------------------*/

/* Subroutine for reading site list from gauge configuration file */
/* Only node 0 reads this list, so same for parallel and serial reading */

void read_site_list(int parallel,gauge_file *gf)
{

  /* All nodes allocate space for site list table, if file is not in
     natural order */

  if(gf->header->order != NATURAL_ORDER)
    {
      gf->rank2rcv = (int32type *)malloc(volume*sizeof(int32type));
      if(gf->rank2rcv == NULL)
	{
	  printf("read_site_list: Can't malloc rank2rcv table\n");
	  terminate(1);
	}

      /* Only node 0 reads the site list */
      
      if(this_node==0)
	{
	  
	  /* Reads receiving site coordinate if file is not in natural order */
	  if(parallel)
	    {
	      if((int)g_read(gf->rank2rcv,sizeof(int32type),volume,gf->fp) != volume )
		{
		  printf("read_site_list: Node %d site list read error %d\n",
			 this_node,errno);
		  terminate(1);	
		}
	    }
	  else
	    {
	      if((int)fread(gf->rank2rcv,sizeof(int32type),volume,gf->fp) != volume )
		{
		  printf("read_site_list: Node %d site list read error %d\n",
			 this_node,errno);
		  terminate(1);	
		}
	    }
	  
	  if(gf->byterevflag==1)byterevn(gf->rank2rcv,volume);
	}

      /* Broadcast result to all nodes */

      broadcast_bytes((char *)gf->rank2rcv,volume*sizeof(int32type));
    }
      
  else gf->rank2rcv = NULL;  /* If no site list */

} /* read_site_list */

/*----------------------------------------------------------------------*/
/* Kept for compatibility */

int read_v3_gauge_hdr(gauge_file *gf, int parallel, int *byterevflag)
{
  /* Provides compatibility with old-style gauge field configurations */

  /* parallel = 1 (true) if all nodes are accessing the file */
  /*            0 for access from node 0 only */

  FILE *fp;
  gauge_header *gh;
  int32type tmp;
  int j;
  int sixtyfourbits;
  Real c1,c2;
  float fc1,fc2;
  char myname[] = "read_v3_gauge_hdr";

  fp = gf->fp;
  gh = gf->header;

  /* Assumes the magic number has already been read */

  /* For cases in which we made a mistake on the T3D and created
     a header with 64-bit integers */

  if(gh->magic_number == 0)
    {
      sixtyfourbits = 1;
      printf("First 4 bytes were zero. Trying to interpret with 64 bit integer format.\n");

      /* Read next 32 bits (without byte reversal) and hope we find it now */
      if(psread_data(parallel,fp,&gh->magic_number,sizeof(gh->magic_number),
	     myname,"magic number")!=0)terminate(1);
    }

  else sixtyfourbits = 0;

  tmp = gh->magic_number;

  if(gh->magic_number == GAUGE_VERSION_NUMBER_V1) 
    {
      printf("Reading as old-style gauge field configuration.\n");
      *byterevflag=0;
    }
  else 
    {
      byterevn((int32type *)&gh->magic_number,1);
      if(gh->magic_number == GAUGE_VERSION_NUMBER_V1) 
	{
	  *byterevflag=1;
	  printf("Reading as old-style gauge field configuration with byte reversal\n");
	  if( sizeof(float) != sizeof(int32type)) {
	    printf("read_v3_gauge_hdr: Can't byte reverse\n");
	    printf("requires size of int32type(%d) = size of float(%d)\n",
		   (int)sizeof(int32type),(int)sizeof(float));
	    terminate(1);
	  }
	}
      else 
	{
	  /* Not recognized as V3 format */
	  /* Restore header to entry state */
	  gh->magic_number = tmp;
	  return 1;  /* error signal */
	}
    }

  /* Read header, do byte reversal, 
     if necessary, and check consistency */
  
  /* Lattice dimensions */
  
  for(j=0;j<4;j++)
    {
      if(psread_byteorder(*byterevflag,parallel,
		       fp,&gh->dims[j],sizeof(gh->dims[j]),
		       myname,"dimensions")!=0)terminate(1);
      /* If 64 bit integers, then we have to read 4 more bytes get the
	 correct low-order bits */
      if(sixtyfourbits)
	if(psread_byteorder(*byterevflag,parallel,
			 fp,&gh->dims[j],sizeof(gh->dims[j]),
			 myname,"dimensions")!=0)terminate(1);
    }

  if(gh->dims[0] != nx || 
     gh->dims[1] != ny ||
     gh->dims[2] != nz ||
     gh->dims[3] != nt)
    {
      /* So we can use this routine to discover the dimensions,
	 we provide that if nx = ny = nz = nt = -1 initially
	 we don't die */
      if(nx != -1 || ny != -1 || nz != -1 || nt != -1)
	{
	  printf("read_v3_gauge_hdr: Incorrect lattice dimensions ");
	  for(j=0;j<4;j++)
	    printf("%d ",gh->dims[j]); 
	  printf("\n");fflush(stdout);terminate(1);
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
  /* Header byte length for this file */
  /* This value is used later in g_seek for locating the gauge link
     matrices */

  if( sixtyfourbits == 0 )
    gh->header_bytes = 2*4 + 5*4;
  else
    gh->header_bytes = 2*4 + 5*8;

  /* Data order - old configuration files have no coordinate list */
  
  gh->order = NATURAL_ORDER;
  
  /* Gauge field parameters */
  
  if(psread_byteorder(*byterevflag,parallel,fp,&fc1,sizeof(float),
		   myname,"c1")!=0)terminate(1);
  if(psread_byteorder(*byterevflag,parallel,fp,&fc2,sizeof(float),
		   myname,"c2")!=0)terminate(1);
  c1=(double)fc1;
  c2=(double)fc2;

  printf("Old format header parameters are %f %f\n",fc1,fc2);
  
  return 0;
} /* read_v3_gauge_hdr */
/*----------------------------------------------------------------------*/
/* Kept for compatibility. */

int read_1996_gauge_hdr(gauge_file *gf, int parallel, int *byterevflag)
{
  /* parallel = 1 (true) if all nodes are accessing the file */
  /*            0 for access from node 0 only */

  FILE *fp;
  gauge_header *gh;
  int32type tmp;
  int j;
  /* We keep this part of the old gauge header, but
     we ignore all but the two parameters */

  struct {                      /* Gauge field parameters */
    int32type n_descript;          /* Number of bytes in character string */
    char   descript[MAX_GAUGE_FIELD_DESCRIPT];  /* Describes gauge field */
    int32type n_param;             /* Number of gauge field parameters */
    float  param[MAX_GAUGE_FIELD_PARAM];        /* GF parameters */
  } gauge_field;
  char myname[] = "read_1996_gauge_hdr";

  fp = gf->fp;
  gh = gf->header;
  
  /* Assumes the magic number has already been read */
  
  tmp = gh->magic_number;
  
  if(gh->magic_number == GAUGE_VERSION_NUMBER_1996) 
    {
      printf("Reading as 1996-style gauge field configuration.\n");
      *byterevflag=0;
    }
  else 
    {
      byterevn((int32type *)&gh->magic_number,1);
      if(gh->magic_number == GAUGE_VERSION_NUMBER_1996) 
	{
	  *byterevflag=1;
	  printf("Reading as 1996-style gauge field configuration with byte reversal\n");
	  if( sizeof(float) != sizeof(int32type)) {
	    printf("read_1996_gauge_hdr: Can't byte reverse\n");
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
	gh->magic_number = tmp;
	return 1;  /* error signal */
      }
    }
  
  /* Read header, do byte reversal, 
     if necessary, and check consistency */
  
  /* Lattice dimensions */
  
  if(psread_byteorder(*byterevflag,parallel,fp,gh->dims,sizeof(gh->dims),
		   myname,"dimensions")!=0)terminate(1);

  if(gh->dims[0] != nx || 
     gh->dims[1] != ny ||
     gh->dims[2] != nz ||
     gh->dims[3] != nt)
    {
      /* So we can use this routine to discover the dimensions,
	 we provide that if nx = ny = nz = nt = -1 initially
	 we don't die */
      if(nx != -1 || ny != -1 || nz != -1 || nt != -1)
	{
	  printf("read_1996_gauge_hdr: Incorrect lattice dimensions ");
	  for(j=0;j<4;j++)
	    printf("%d ",gh->dims[j]); 
	  printf("\n");fflush(stdout);terminate(1);
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
  
  /* Header byte length */
  
  if(psread_byteorder(*byterevflag,parallel,fp,
		   &gh->header_bytes,sizeof(gh->header_bytes),
		   myname,"header size")!=0)terminate(1);
  
  /* Data order */
  
  if(psread_byteorder(*byterevflag,parallel,fp,
		   &gh->order,sizeof(gh->order),
		   myname,"order")!=0)terminate(1);
  
  /* Length of gauge field descriptor */
  
  if(psread_byteorder(*byterevflag,parallel,fp,
		   &gauge_field.n_descript,sizeof(gauge_field.n_descript),
		   myname,"n_descript")!=0)terminate(1);

  if(gauge_field.n_descript > MAX_GAUGE_FIELD_DESCRIPT)
    {
      printf("read_1996_gauge_hdr: gauge field descriptor length %d\n",
	     gauge_field.n_descript);
      printf(" exceeds allocated space %d\n",
	     MAX_GAUGE_FIELD_DESCRIPT);
      terminate(1);
    }
  
  /* Gauge field descriptor */
  
  /* We read the specified length, rather than the allocated length */
  /* Read without byte reversal */

  if(psread_data(parallel,fp,gauge_field.descript,sizeof(gauge_field.descript),
	      myname,"descrip")!=0)terminate(1);

  /* Assures termination of string */
  gauge_field.descript
    [gauge_field.n_descript-1] = '\0';

  printf("gauge_field.descript: %s\n", gauge_field.descript);

  /* Number of gauge field parameters */
  
  if(psread_byteorder(*byterevflag,parallel,fp,
		   &gauge_field.n_param,sizeof(gauge_field.n_param),
		   myname,"n_param")!=0)terminate(1);

  if(gauge_field.n_param > MAX_GAUGE_FIELD_PARAM )
    {
      printf("read_1996_gauge_hdr: gauge field parameter vector length %d\n",
	     gauge_field.n_param);
      printf(" exceeds allocated space %d\n",
	     MAX_GAUGE_FIELD_PARAM);
      terminate(1);
    }
  
  /* Gauge field parameters */
  
  for(j=0;j<gauge_field.n_param;j++)
    {
      if(psread_byteorder(*byterevflag,parallel,fp,
		     &gauge_field.param[j],sizeof(gauge_field.param[j]),
		     myname,"gauge param")!=0)terminate(1);
      printf("gauge_field.param[%d] = %f\n", j, gauge_field.param[j]);
    }
  
  /* Since there aren't many of these lattices in circulation, 
     we simply ignore the information in this header */
  
  return 0;
  
} /* read_1996_gauge_hdr */

/*----------------------------------------------------------------------*/
/* Kept for compatibility. */

int read_fnal_gauge_hdr(gauge_file *gf, int parallel, int *byterevflag)
{
  /* parallel = 1 (true) if all nodes are accessing the file */
  /*            0 for access from node 0 only */

  FILE *fp;
  gauge_header *gh;
  int32type tmp;
  int j;
  char myname[] = "read_fnal_gauge_hdr";
  int32type size_of_element, elements_per_site, gmtime_stamp;


  fp = gf->fp;
  gh = gf->header;
  
  /* Assumes the magic number has already been read */
  
  tmp = gh->magic_number;
  
  if(gh->magic_number == GAUGE_VERSION_NUMBER_FNAL) 
    {
      printf("Reading as FNAL-style gauge field configuration.\n");
      *byterevflag=0;
    }
  else 
    {
      byterevn((int32type *)&gh->magic_number,1);
      if(gh->magic_number == GAUGE_VERSION_NUMBER_FNAL) 
	{
	  *byterevflag=1;
	  printf("Reading as FNAL-style gauge field configuration with byte reversal\n");
	  if( sizeof(float) != sizeof(int32type)) {
	    printf("read_fnal_gauge_hdr: Can't byte reverse\n");
	    printf("requires size of int32type(%d) = size of float(%d)\n",
		   (int)sizeof(int32type),(int)sizeof(float));
	    terminate(1);
	  }
	}
      /* Not recognized as fnal format */
      else
      {
	/* Not recognized as fnal format */
	/* Restore header to entry state */
	gh->magic_number = tmp;
	return 1;  /* error signal */
      }
    }
  
  /* Read header, do byte reversal, 
     if necessary, and check consistency */

  if(psread_byteorder(*byterevflag,parallel,fp,&gmtime_stamp,
	sizeof(int32type), myname,"gmtime_stamp")!=0)terminate(1);

  if(psread_byteorder(*byterevflag,parallel,fp,&size_of_element,
	sizeof(int32type), myname,"size_of_element")!=0)terminate(1);
  
  if(psread_byteorder(*byterevflag,parallel,fp,&elements_per_site,
	sizeof(int32type), myname,"elements_per_site")!=0)terminate(1);
  
  if( (size_of_element != 4 )|| (elements_per_site != 72) ) 
	node0_printf("This does not look like a single precision gauge field\nsize-of-element= %d\t elements-per-site= %d\n",size_of_element,elements_per_site);

  /* Lattice dimensions */
  
  if(psread_byteorder(*byterevflag,parallel,fp,gh->dims,sizeof(gh->dims),
		   myname,"dimensions")!=0)terminate(1);

  if(gh->dims[0] != nx || 
     gh->dims[1] != ny ||
     gh->dims[2] != nz ||
     gh->dims[3] != nt)
    {
      /* So we can use this routine to discover the dimensions,
	 we provide that if nx = ny = nz = nt = -1 initially
	 we don't die */
      if(nx != -1 || ny != -1 || nz != -1 || nt != -1)
	{
	  printf("read_fnal_gauge_hdr: Incorrect lattice dimensions ");
	  for(j=0;j<4;j++)
	    printf("%d ",gh->dims[j]); 
	  printf("\n");fflush(stdout);terminate(1);
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
  
  /* Data order */
  
  if(psread_byteorder(*byterevflag,parallel,fp,
		   &gh->order,sizeof(gh->order),
		   myname,"order")!=0)terminate(1);
  
  /* This is the end of the Fermilab style header */
  /* set the header_bytes head field to 36 since it is not an FNAL field */
	gh->header_bytes=36;
  
  return 0;
  
} /* read_fnal_gauge_hdr */
/*----------------------------------------------------------------------*/

int read_gauge_hdr(gauge_file *gf, int parallel)
{
  /* parallel = 1 (TRUE) if all nodes are accessing the file */
  /*            0        for access from node 0 only */

  FILE *fp;
  gauge_header *gh;
  int32type tmp, btmp;
  int j;
  int byterevflag;
  char myname[] = "read_gauge_hdr";
  int i;
  QCDheader *hdr;
  int dims[4];
  int ARCHYES=0;
  u_int32type chksum;

  fp = gf->fp;
  gh = gf->header;

  /* Read header, do byte reversal, 
     if necessary, and check consistency */
  
  /* Read and verify magic number */

  if(psread_data(parallel, fp,&gh->magic_number,sizeof(gh->magic_number),
			 myname,"magic number")!=0)terminate(1);

  tmp = gh->magic_number;
  btmp = gh->magic_number;
  byterevn((int32type *)&btmp,1);

  /** See if header chunk is BEGI = 1111836489 for big endian
      or the byte reverse 1229407554 for little endian **/

  if(tmp == GAUGE_VERSION_NUMBER_ARCHIVE) 
    {
      printf("reading as archive format\n"); 
      ARCHYES=1;
      byterevflag=0;
    }
  else if(btmp == GAUGE_VERSION_NUMBER_ARCHIVE) 
	{
	  printf("reading as archive format with byte reversal\n"); 
	  ARCHYES=1;
	  byterevflag=1;	/* not really needed */
	  gh->magic_number = btmp;
	  if( sizeof(float) != sizeof(int32type)) {
	    printf("%s: Can't byte reverse\n",myname);
	    printf("requires size of int32type(%d) = size of float(%d)\n",
		   (int)sizeof(int32type),(int)sizeof(float));
	    terminate(1);
	  }
	}
  else if(tmp == GAUGE_VERSION_NUMBER) 
    {
      byterevflag=0;
    }
  else if(btmp == GAUGE_VERSION_NUMBER) 
    {
      byterevflag=1;
      gh->magic_number = btmp;
      /**      printf("Reading with byte reversal\n"); **/
      if( sizeof(float) != sizeof(int32type)) {
	printf("%s: Can't byte reverse\n",myname);
	printf("requires size of int32type(%d) = size of float(%d)\n",
	       (int)sizeof(int32type),(int)sizeof(float));
	terminate(1);
      }
    }
  else if(tmp == LIME_MAGIC_NO || btmp == LIME_MAGIC_NO)
    {
      /* LIME format suggests a SciDAC file */
      /* We do not read any further here:  Set flag and return */
      printf("%s: Reading as a SciDAC formatted file\n",myname);
      gh->magic_number = LIME_MAGIC_NO;
      return 0;
    }
  else
    {
      /* Try old-style configurations */
      if( (read_fnal_gauge_hdr(gf,parallel,&byterevflag) != 0) &&
	  (read_v3_gauge_hdr(gf,parallel,&byterevflag) != 0) &&
	  (read_1996_gauge_hdr(gf,parallel,&byterevflag) != 0) )
	{
	  /* End of the road. */
	  printf("%s: Unrecognized magic number in gauge configuration file header.\n",myname);
	  printf("Expected %x but read %x\n",
		 GAUGE_VERSION_NUMBER,tmp);
	  printf("Expected %s but read %s\n",
		 (char *)GAUGE_VERSION_NUMBER,(char *)tmp);
	  terminate(1);
	}
      return byterevflag;
    }
  
  /* Read and process header information */
  /* Get lattice dimensions */
  
  /* Special processing for NERSC archive format files */
  if(ARCHYES == 1) 
    {
      gf->header->order = NATURAL_ORDER;
      
      if(parallel) {
	fprintf(stderr,
		"%s: Must use reload_serial with archive files for now.\n",
		myname);
	terminate(1);
      }

      /* Reads the entire header of the archive file */
      hdr = qcdhdr_get_hdr(fp);

      /* Get dimensions */
      if (qcdhdr_get_int("DIMENSION_1",hdr,dims+0)==FAILURE)
	error_exit("DIMENSION_1 not present");
      if (qcdhdr_get_int("DIMENSION_2",hdr,dims+1)==FAILURE)
	error_exit("DIMENSION_2 not present");
      if (qcdhdr_get_int("DIMENSION_3",hdr,dims+2)==FAILURE)
	error_exit("DIMENSION_3 not present");
      if (qcdhdr_get_int("DIMENSION_4",hdr,dims+3)==FAILURE)
	error_exit("DIMENSION_4 not present");

      for(i=0; i<4; i++) gh->dims[i] = dims[i];

      /* Get archive checksum */
      if (qcdhdr_get_int32x("CHECKSUM",hdr,&chksum)==FAILURE)
	error_exit("CHECKSUM not present");
      gf->check.sum31 = chksum;
    }

  /* not a archive lattice - read lattice dimensions */
  else
    {
      if(psread_byteorder(byterevflag,parallel,fp,gh->dims,sizeof(gh->dims),
			  myname,"dimensions")!=0)terminate(1);
    }

  /* Check lattice dimensions for consistency */

  if(gh->dims[0] != nx || 
     gh->dims[1] != ny ||
     gh->dims[2] != nz ||
     gh->dims[3] != nt)
    {
      /* So we can use this routine to discover the dimensions,
	 we provide that if nx = ny = nz = nt = -1 initially
	 we don't die */
      if(nx != -1 || ny != -1 || nz != -1 || nt != -1)
	{
	  printf("%s: Incorrect lattice dimensions ",myname);
	  for(j=0;j<4;j++)
	    printf("%d ",gh->dims[j]); 
	  printf("\n");fflush(stdout);terminate(1);
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

  if(ARCHYES) {

  /* After we are done processing the archive header information, we
     discard it */

    qcdhdr_destroy_hdr(hdr);
  }

  else {
    
    /* Read date and time stamp */
    
    if(psread_data(parallel,fp,gh->time_stamp,sizeof(gh->time_stamp),
		   myname,"time stamp")!=0)terminate(1);
    
    /* Read header byte length */
    
    gh->header_bytes = sizeof(gh->magic_number) + sizeof(gh->dims) + 
      sizeof(gh->time_stamp) + sizeof(gh->order);
    
    /* Read data order */
    
    if(psread_byteorder(byterevflag,parallel,fp,&gh->order,sizeof(gh->order),
			myname,"order parameter")!=0)terminate(1);
  }  

  return byterevflag;
  
} /* read_gauge_hdr */

/*---------------------------------------------------------------------------*/

gauge_file *r_serial_i(char *filename)
{
  /* Returns file descriptor for opened file */

  gauge_header *gh;
  gauge_file *gf;
  FILE *fp;
  int byterevflag;
  char editfilename[513];

  /* All nodes set up a gauge file and gauge header structure for reading */

  gf = setup_input_gauge_file(filename);
  gh = gf->header;

  /* File opened for serial reading */
  gf->parallel = 0;

  /* Node 0 alone opens the file and reads the header */

  g_sync();

  if(this_node==0)
    {
      fp = fopen(filename, "rb");
      if(fp == NULL)
	{
	  /* If this is a partition format SciDAC file the node 0 name
	     has an extension ".vol0000".  So try again. */
	  printf("r_serial_i: Node %d can't open file %s, error %d\n",
		 this_node,filename,errno);fflush(stdout);
	  strncpy(editfilename,filename,504);
	  editfilename[504] = '\0';  /* Just in case of truncation */
	  strcat(editfilename,".vol0000");
	  printf("r_serial_i: Trying SciDAC partition volume %s\n",editfilename);
	  fp = fopen(editfilename, "rb");
	  if(fp == NULL)
	    {
	      printf("r_serial_i: Node %d can't open file %s, error %d\n",
		     this_node,editfilename,errno);fflush(stdout);terminate(1);
	    }
	  printf("r_serial_i: Open succeeded\n");
	}
      
      gf->fp = fp;

      byterevflag = read_gauge_hdr(gf,SERIAL);

    }

  else gf->fp = NULL;  /* The other nodes don't know about this file */

  /* Broadcast the byterevflag from node 0 to all nodes */

  broadcast_bytes((char *)&byterevflag,sizeof(byterevflag));
  gf->byterevflag = byterevflag;
  
  /* Node 0 broadcasts the header structure to all nodes */
  
  broadcast_bytes((char *)gh,sizeof(gauge_header));

  /* No further processing here if this is a SciDAC file */
  if(gh->magic_number == LIME_MAGIC_NO)
    return gf;

  /* Read site list and broadcast to all nodes */

  read_site_list(SERIAL,gf);

  return gf;

}/* r_serial_i */

/*----------------------------------------------------------------------*/

/* Here only node 0 reads the gauge configuration from a binary file */

void r_serial(gauge_file *gf)
{
  /* gf  = gauge configuration file structure */

  FILE *fp;
  gauge_header *gh;
  char *filename;
  int byterevflag;

  off_t offset ;            /* File stream pointer */
  off_t gauge_check_size;   /* Size of gauge configuration checksum record */
  off_t coord_list_size;    /* Size of coordinate list in bytes */
  off_t head_size;          /* Size of header plus coordinate list */
  off_t checksum_offset;    /* Where we put the checksum */
  int rcv_rank, rcv_coords;
  int destnode;
  int i,k;
  int x,y,z,t;
  int buf_length,where_in_buf;
  gauge_check test_gc;
  u_int32type *val;
  int rank29,rank31;
  fsu3_matrix *lbuf;
  fsu3_matrix tmpsu3[4];
  char myname[] = "r_serial";
  int ii,dir,row,col;
  int idest;

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
      else if(gh->magic_number == GAUGE_VERSION_NUMBER_1996)
	gauge_check_size =  4;
      else
	gauge_check_size = 0;
      
      if(gf->header->order == NATURAL_ORDER)coord_list_size = 0;
      else coord_list_size = sizeof(int32type)*volume;
      checksum_offset = gf->header->header_bytes + coord_list_size;
      head_size = checksum_offset + gauge_check_size;
      
      /* Allocate space for read buffer */

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

      if( fseeko(fp,offset,SEEK_SET) < 0 ) 
	{
	  printf("%s: Node 0 fseeko %lld failed error %d file %s\n",
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
	    
	    if( (int)fread(lbuf,4*sizeof(fsu3_matrix),buf_length,fp) 
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
	  if( fseeko(fp,checksum_offset,SEEK_SET) < 0 ) 
	    {
	      printf("%s: Node 0 fseeko %lld failed error %d file %s\n",
		    myname,(long long)offset,errno,filename);
	      fflush(stdout);terminate(1);   
	    }
	  read_checksum(SERIAL,gf,&test_gc);
	}
      fflush(stdout);
      free(lbuf);
    }
  
} /* r_serial */

/*----------------------------------------------------------------------*/

void r_serial_arch(gauge_file *gf)
{
  /* gf  = gauge configuration file structure */

  FILE *fp;
  gauge_header *gh;
  char *filename;
  int byterevflag;

  off_t gauge_check_size;   /* Size of gauge configuration checksum record */
  int rcv_rank, rcv_coords;
  int destnode;
  int i,k;
  int x,y,z,t;
  gauge_check test_gc;
  u_int32type *val;
  int rank29,rank31;
  fsu3_matrix tmpsu3[4];
  char myname[] = "r_serial_arch";

  int mu,a,b,p;
  float *uin, *q;
  int big_end;
  float U[4][18];
  u_int32type chksum;
  
  fp = gf->fp;
  gh = gf->header;
  filename = gf->filename;
  byterevflag = gf->byterevflag;

  if(this_node == 0)
    {
      gauge_check_size = 0;
      
      if(gf->parallel)
	printf("%s: Attempting serial read from parallel file \n",myname);

      big_end = big_endian();
      /* printf("big_end is %d\n", big_end); */
      uin = (float *) malloc(nx*ny*nz*48*sizeof(float));
      if(uin == NULL)
	{
	  printf("%s: Node %d can't malloc uin buffer to read timeslice\n",
		 myname,this_node);
	  printf("recompile with smaller read buffers: uin\n");
	  fflush(stdout);
	  terminate(1);
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
	if( (int)fread(uin,48*sizeof(float),1,fp) != 1)
	  {
	    printf("%s: node %d gauge configuration read error %d file %s\n",
		   myname,this_node,errno,filename); 
	    fflush(stdout); terminate(1);
	  }

	if (!big_end) byterevn((int32type *)uin,48);
	q = uin;
	for (mu=0;mu<4;mu++) {
	  for (p=0;p<12;p++) {
	    chksum += *(u_int32type *) q;
	    U[mu][p] = (float) *(q++);
	  }
	  complete_U(U[mu]);
	  /**
	  for (p=0;p<18;p++) printf("p=%d, e=%f\n", p, U[mu][p]);
	  **/
		 
          for(a=0; a<3; a++) for(b=0; b<3; b++) { 
	    tmpsu3[mu].e[a][b].real = U[mu][2*(3*a+b)];
     /*	    printf("real: p=%d, mu=%d, e=%f\n", p,mu,U[mu][2*(3*a+b)]); */
	    tmpsu3[mu].e[a][b].imag = U[mu][2*(3*a+b)+1];
     /*	    printf("imag: p=%d, mu=%d, e=%f\n", p,mu,U[mu][2*(3*a+b)+1]); */
	  } 
	}

	if(destnode==0){	
	  /* just copy links */
	  i = node_index(x,y,z,t);
     /*   printf("lattice node_index = %d, mu = %d\n", i, mu); */
	  /* Copy from tmpsu3 to site structure, converting to double */
	  f2d_4mat(tmpsu3,&lattice[i].link[0]);
	} else {		
	  /* send to correct node */
	  send_field((char *)tmpsu3, 4*sizeof(fsu3_matrix),destnode);
	}
      } 
      /* The node which contains this site reads message */
      else {	
	/* for all nodes other than node 0 */
	if(this_node==destnode){
	  i = node_index(x,y,z,t);
	  get_field((char *)tmpsu3,4*sizeof(fsu3_matrix),0);
	  /* Store in site structure, converting to generic precision */
	  f2d_4mat(tmpsu3,&lattice[i].link[0]);
	}
      }

      /* Any needed byte reversing was already done. Compute MILC
         checksums. At this point tmpsu3 on destnode contains the link
         matrices we just read */

      if(this_node==destnode)
	{
	  /* Accumulate checksums */
	  for(k = 0, val = (u_int32type *)tmpsu3;
	      k < 4*(int)sizeof(fsu3_matrix)/(int)sizeof(int32type); k++, val++)
   	    {
	      test_gc.sum29 ^= (*val)<<rank29 | (*val)>>(32-rank29);
	      test_gc.sum31 ^= (*val)<<rank31 | (*val)>>(32-rank31);
	      rank29++; if(rank29 >= 29)rank29 = 0;
	      rank31++; if(rank31 >= 31)rank31 = 0;
	    }
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





/*----------------------------------------------------------------------*/

void r_serial_f(gauge_file *gf)

/* Close the file and free associated structures */
{
  g_sync();
  if(this_node==0)
    {
      if(gf->parallel)
	printf("r_serial_f: Attempting serial close on parallel file \n");

      fclose(gf->fp);
    }
  
  if(gf->rank2rcv != NULL)free(gf->rank2rcv);
  
  /* Do not free gf and gf->header so calling program can use them */

} /* r_serial_f */

/*---------------------------------------------------------------------------*/

/* Write site list - only for checkpoint files */

void write_site_list(FILE *fp, gauge_header *gh)
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
  
  offset = gh->header_bytes + 
    sizeof(int32type)*sites_on_node*this_node;

  cbuf = (int32type *)malloc(sites_on_node*sizeof(int32type));
  if(cbuf == NULL)
    {
      printf("write_site_list: node %d can't malloc cbuf\n",this_node);
      fflush(stdout);terminate(1);   
    }

  if( g_seek(fp,offset,SEEK_SET) < 0 ) 
    {
      printf("write_site_list: node %d g_seek %ld failed errno %d\n",
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
	printf("write_site_list: Node %d coords write error %d\n",
	       this_node,errno);fflush(stdout);terminate(1);   
      }

  free(cbuf);
} /* write_site_list */

/*---------------------------------------------------------------------------*/

/* Open a file for parallel writing */
gauge_file *parallel_open(int order, char *filename)
{
  /* All nodes open the same filename */
  /* Returns a file structure describing the opened file */

  /* order = NATURAL_ORDER for coordinate natural order 
           = NODE_DUMP_ORDER for node-dump order */

  FILE *fp;
  gauge_file *gf;
  gauge_header *gh;

  /* Set up gauge file and gauge header structures and load header values */
  gf = setup_output_gauge_file();
  gh = gf->header;

  gh->order = order;

  /* All nodes open the requested file */

  fp = g_open(filename, "wb");
  if(fp == NULL)
    {
      printf("parallel_open: Node %d can't open file %s, error %d\n",
	     this_node,filename,errno);fflush(stdout);terminate(1);
    }

  /* Node 0 writes the header */

  if(this_node==0)
    pwrite_gauge_hdr(fp,gh);

  broadcast_bytes((char *)&gh->header_bytes,sizeof(gh->header_bytes));
  
  /* All nodes write site list to file if order is not natural */

  if(order != NATURAL_ORDER)write_site_list(fp,gh);
  
  /* Assign values to file structure */

  gf->fp             = fp;
  gf->filename        = filename;
  gf->byterevflag    = 0;            /* Not used for writing */
  gf->parallel       = 1;            /* File opened in parallel */

  return gf;
} /* parallel_open */

/*-----------------------------------------------------------------------*/

/* Open a file for parallel writing in natural order */
gauge_file *w_parallel_i(char *filename)
{
  /* All nodes open the same filename */
  /* Returns a file structure describing the opened file */

  return parallel_open(NATURAL_ORDER,filename);

} /* w_parallel_i */


/*---------------------------------------------------------------------------*/
/* Open a file for parallel writing in node-dump order */
gauge_file *w_checkpoint_i(char *filename)
{
  /* All nodes open the same filename */
  /* Returns a file structure describing the opened file */

  return parallel_open(NODE_DUMP_ORDER,filename);

} /* w_checkpoint_i */

/*---------------------------------------------------------------------------*/

/* Position gauge configuration file for writing in parallel */
/* Returns pointer to malloc'ed write buffer */

fsu3_matrix *w_parallel_setup(gauge_file *gf, off_t *checksum_offset)
{
  /* gf  = file descriptor as opened by w_checkpoint_i */

  FILE *fp;
  gauge_header *gh;
  fsu3_matrix *lbuf;

  off_t offset ;           /* File stream pointer */
  off_t gauge_node_size;   /* Size of a gauge configuration block for
                              all sites on one node */
  off_t coord_list_size;   /* Size of coordinate list in bytes */
  off_t head_size;         /* Size of header plus coordinate list */
  off_t gauge_check_size;  /* Size of checksum */
  char myname[] = "w_parallel_setup";

  if(!gf->parallel)
    printf("%s: Attempting parallel write to serial file.\n",myname);

  lbuf = (fsu3_matrix *)malloc(MAX_BUF_LENGTH*4*sizeof(fsu3_matrix));
  if(lbuf == NULL)
    {
      printf("%s: Node %d can't malloc lbuf\n",myname,this_node);
      fflush(stdout);
      terminate(1);
    }

  fp = gf->fp;
  gh = gf->header;

  gauge_node_size = sites_on_node*4*sizeof(fsu3_matrix) ;

  if(gf->header->order == NATURAL_ORDER)coord_list_size = 0;
  else coord_list_size = sizeof(int32type)*volume;
  head_size = gf->header->header_bytes + coord_list_size;
  *checksum_offset = head_size;
  gauge_check_size = sizeof(gf->check.sum29) + sizeof(gf->check.sum31);

  offset = head_size + gauge_check_size;

  /* Each node writes its gauge configuration values */

  offset += gauge_node_size*this_node;
  
  if( g_seek(fp,offset,SEEK_SET) < 0 ) 
    {
      printf("%s: Node %d g_seek %ld failed error %d file %s\n",
	     myname,this_node,(long)offset,errno,gf->filename);
      fflush(stdout);terminate(1);
    }

  return lbuf;
} /* w_parallel_setup */

/*---------------------------------------------------------------------------*/
/* Write parallel gauge configuration in coordinate natural order */

void w_parallel(gauge_file *gf)
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

void w_checkpoint(gauge_file *gf)
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

void w_parallel_f(gauge_file *gf)
{
  /* Close file (if still active) and release header and file structures */

  g_sync();
  if(gf->fp != NULL)
    {
      if(!gf->parallel)
	printf("w_parallel_f: Attempting parallel close on serial file.\n");
      
      g_close(gf->fp);
      gf->fp = NULL;
    }

  /* Node 0 writes ascii info file */

  if(this_node == 0)write_gauge_info_file(gf);

  /* Do not free gf and gf->header so calling program can use them */

} /* w_parallel_f */


/*---------------------------------------------------------------------------*/

gauge_file *r_parallel_i(char *filename)
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
void r_parallel(gauge_file *gf)
{
  /* gf  = gauge configuration file structure */

  FILE *fp;
  gauge_header *gh;
  char *filename;
  int byterevflag;
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
  byterevflag = gf->byterevflag;

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
  else if(gh->magic_number == GAUGE_VERSION_NUMBER_1996)
    gauge_check_size =  4;
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

void r_parallel_f(gauge_file *gf)
{
  /* Close file (if active) and release header and file structures */

  g_sync();
  if(gf->fp != NULL)
    {
      if(!gf->parallel)
	printf("r_parallel_f: Attempting parallel close on serial file.\n");
      g_close(gf->fp);
      gf->fp = NULL;
    }

  /* Do not free gf and gf->header so calling program can use them */

 } /* r_parallel_f */

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
  FILE *fp;
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
    fclose(fp);
    gf->fp = NULL;
    fflush(stdout);
  }

  return gf;
}

/*---------------------------------------------------------------------------*/

/* Save a lattice in ASCII format serially (node 0 only) */

gauge_file *save_ascii(char *filename) {
  FILE *fp;
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
      if( this_node==0 && newnode!=0 )send_field((char *)lbuf,4,newnode);
      if( this_node==newnode && newnode!=0 )get_field((char *)lbuf,4,0);
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
    fclose(fp);
    fflush(stdout);
    }
  return gf;
}

/* Read lattice dimensions from a binary file and close the file */
void read_lat_dim_gf(char *filename, int *ndim, int dims[]){
  gauge_file *gf;
  int i;
  
  /* Only four dimensions here */
  *ndim = 4;

  /* Open the file */
  nx = -1; ny = -1; nz = -1; nt = -1;
  gf = r_serial_i(filename);

  for(i = 0; i < *ndim; i++)
    dims[i] = gf->header->dims[i];

  r_serial_f(gf);
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
  r_parallel(gf);
  r_parallel_f(gf);

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

  FILE *outfile;
  site *s;
  u_int32type chksum, utmp, *p32;
  char sums[30];
  OUTPUT_TYPE *uout;
  int big_end_p; 
  double ssplaq, stplaq, avgtrace, avgplaq;
  float tmpflt;
  double trace;
  int mu,a,b,vol3,tslice;

  /* Check which end is up */
  big_end_p = big_endian();
  /* if(this_node == 0) printf("big_end_p is %d\n", big_end_p); */

  /* Set up gauge file and gauge header structures and load header values */
  gf = setup_output_gauge_file();
  gh = gf->header;
  
  /* Compute plaquette, trace and checksum */
  d_plaquette(&ssplaq, &stplaq);
  avgplaq = (ssplaq+stplaq)/6.0;
  trace = 0.0;
  chksum = 0;
  FORALLSITES(i,s) {
    for(mu=0; mu<4; ++mu) {
      trace += (trace_su3(&(s->link[mu]))).real;
      for(a=0; a<2; a++) for(b=0; b<3; b++) {
	tmpflt = s->link[mu].e[a][b].real;
	p32 = (u_int32type *) &tmpflt;
	chksum += *p32;
	tmpflt = s->link[mu].e[a][b].imag;
	p32 = (u_int32type *) &tmpflt;
	chksum += *p32;
      }
    }
  }
  g_doublesum( &trace);
  avgtrace = trace/(double)(volume*12);

  /* All nodes send their checksum to node 0 */
  for(j=1; j<numnodes(); j++){
    if(this_node == 0)send_field((char *)lbuf,4,j);
    if(this_node == j){
      get_field((char *)lbuf,4,0);
      send_integer(0, (int *)&chksum);
    }
    if(this_node == 0){
      receive_integer(j, (int *)&utmp);
      chksum += utmp;
    }
  }

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
    write_appl_gauge_info(outfile);
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
	if( this_node==0 && newnode!=0 )send_field((char *)lbuf,4,newnode);
	if( this_node==newnode && newnode!=0 )get_field((char *)lbuf,4,0);
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
      if(fwrite(uout,48*vol3*sizeof(OUTPUT_TYPE),1,outfile) != 1)
	printf("fwrite bombed...\n");
      fflush(outfile);
    }
  }

  if(this_node==0){
    fclose(outfile);
    printf("Wrote archive gauge file %s\n",filename);
    free(uout);
  }

  g_sync();
  return gf;
}

/*---------------------------------------------------------------------------*/
gauge_file *save_parallel_archive(char *filename) {
  /* All nodes write in archive file format */
  printf("Parallel archive saves are not implemented, yet\n");
  return NULL;
}
