/*********************** io_lat_utils.c *************************/
/* MIMD version 7 */

/* routines for gauge configuration input/output. */
/* This works for most machines.  Wrappers for parallel I/O
   are in io_ansi.c, io_piofs.c, or io_paragon2.c */

/* Modifications */
/* 7/19/05 Separated from io_lat4.c C.D. */
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


void complete_Ud(double *u) {
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
    len = (int)( q -  line );

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
  if(g_write(src,size,1,fp) != 1)
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
  if(g_read(src,size,1,fp) != 1)
    {
      int errsv = errno;
      printf("%s: Node %d %s read error %d %s\n",
	    myname,this_node,descrip,errsv,strerror(errsv));
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
      int errsv = errno;
      printf("%s: Node %d %s read error %d %s\n",
	    myname,this_node,descrip,errsv,strerror(errsv));
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
					      must use s, d, e, f, lu, or g */
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
      else if(strstr(fmt,"u") != NULL)
	fprintf(fpout,fmt,*(unsigned int *)data);
      else if(strstr(fmt,"lu") != NULL)
	fprintf(fpout,fmt,*(unsigned long *)data);
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
  size_t bytes;
  char *data;
  float tt;

  /* Check for valid keyword */

  for(i=0;strlen(gauge_info_keyword[i])>0 &&
      strcmp(gauge_info_keyword[i],keyword) != 0; i++);
  if(strlen(gauge_info_keyword[i])==0)
    printf("sprint_gauge_info_item: WARNING: keyword %s not in table\n",
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
	  printf("sprint_gauge_info_item: Unrecognized data type %s\n",fmt);
	  return 1;
	}
    }
  snprintf(string+bytes,nstring-bytes,"\n");
  bytes = strlen(string);
  if(bytes >= nstring)return 1;

  return 0;
}

/*----------------------------------------------------------------------*/
/* Write generic information to info file */

void write_generic_gauge_info(FILE *info_fp, gauge_file *gf)
{
  char sums[20];
  gauge_header *gh = gf->header;

  write_gauge_info_item(info_fp,"magic_number","%d",(char *)&gh->magic_number,0,0);
  write_gauge_info_item(info_fp,"time_stamp","\"%s\"",gh->time_stamp,0,0);
  sprintf(sums,"%x %x",gf->check.sum29,gf->check.sum31);
  write_gauge_info_item(info_fp,"checksums","\"%s\"",sums,0,0);
  write_gauge_info_item(info_fp,"nx","%d",(char *)&nx,0,0);
  write_gauge_info_item(info_fp,"ny","%d",(char *)&ny,0,0);
  write_gauge_info_item(info_fp,"nz","%d",(char *)&nz,0,0);
  write_gauge_info_item(info_fp,"nt","%d",(char *)&nt,0,0);
}

/*----------------------------------------------------------------------*/
/* Open, write, and close the ASCII info file */

void write_gauge_info_file(gauge_file *gf)
{
  FILE *info_fp;
  char info_filename[256];

  /* Construct header file name from lattice file name 
   by adding filename extension to lattice file name */

  strcpy(info_filename,gf->filename);
  strcat(info_filename,ASCII_INFO_EXT);

  /* Open header file */
  
  if((info_fp = g_open(info_filename,"w")) == NULL)
    {
      printf("write_gauge_info_file: Can't open ascii info file %s\n",info_filename);
      return;
    }
  
  /* Write application information to info file */
  write_appl_gauge_info(info_fp, gf);

  g_close(info_fp);

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

  memset(gh, 0, sizeof(gauge_header));

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
		  int errsv = errno;
		  printf("read_site_list: Node %d site list read error %d %s\n",
			 this_node,errsv,strerror(errsv));
		  terminate(1);	
		}
	    }
	  else
	    {
	      if((int)g_read(gf->rank2rcv,sizeof(int32type),volume,gf->fp) != volume )
		{
		  int errsv = errno;
		  printf("read_site_list: Node %d site list read error %d %s\n",
			 this_node,errno,strerror(errsv));
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
  
  if(gh->magic_number == IO_UNI_MAGIC) 
    {
      printf("Reading as FNAL-style gauge field configuration.\n");
      *byterevflag=0;
    }
  else 
    {
      byterevn((int32type *)&gh->magic_number,1);
      if(gh->magic_number == IO_UNI_MAGIC) 
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
  QCDheader *hdr = NULL;
  int dims[4];
  int ARCHYES=0;
  u_int32type chksum = 0;
  char *datatype;
  char *floatpt;

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
	  printf("Expected %lu but read %lu\n",
		 (unsigned long)GAUGE_VERSION_NUMBER,(unsigned long)tmp);
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

      /* Get archive datatype */
      if (qcdhdr_get_str("DATATYPE",hdr,&datatype)==FAILURE)
	error_exit("DATATYPE not present");
      /* Two choices currently */
      gf->dataformat = ARCHIVE_3x2;
      if(strcmp(" 4D_SU3_GAUGE_3x3",datatype) == 0)
	gf->dataformat = ARCHIVE_3x3;

      /* Get archive floating point format */
      gf->precision = 1;
      if (qcdhdr_get_str("FLOATING_POINT",hdr,&floatpt)==FAILURE)
	fprintf(stderr,"FLOATING_POINT not present.  Assuming IEEE32BIG.\n");
      else if(strcmp(" IEEE64BIG",floatpt) == 0)
	gf->precision = 2;
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
  g_sync(); /* Make sure everyon has opened before attempting to write */
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

/*---------------------------------------------------------------------------*/

/* Position gauge configuration file for writing in parallel */
/* Returns pointer to malloc'ed write buffer */

fsu3_matrix *w_parallel_setup(gauge_file *gf, off_t *checksum_offset)
{
  /* gf  = file descriptor as opened by w_checkpoint_i */

  FILE *fp;
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

void w_serial_f(gauge_file *gf)

/* Close the file and free associated structures */
{
  g_sync();
  if(this_node==0)
    {
      if(gf->parallel)
	printf("w_serial_f: Attempting serial close on parallel file \n");

      g_close(gf->fp);
    }

  /* Node 0 writes ascii info file */

  if(this_node == 0)write_gauge_info_file(gf);

  /* Do not free gf and gf->header so calling program can use them */

} /* w_serial_f */

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
      fp = g_open(filename, "rb");
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
	  fp = g_open(editfilename, "rb");
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

void r_serial_f(gauge_file *gf)

/* Close the file and free associated structures */
{
  g_sync();
  if(this_node==0)
    {
      if(gf->parallel)
	printf("r_serial_f: Attempting serial close on parallel file \n");

      g_close(gf->fp);
    }
  
  if(gf->rank2rcv != NULL)free(gf->rank2rcv);
  
  /* Do not free gf and gf->header so calling program can use them */

} /* r_serial_f */

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

/*---------------------------------------------------------------------------*/
#define INFO_DATA_SIZE 4096

char *read_info_file(char base_filename[]){
  char info_filename[512];
  char *info, *buf;
  FILE *info_fp;
  int n;

  /* Allocate space for metadata on all nodes */

  info = (char *)malloc(INFO_DATA_SIZE*sizeof(char));
  if(info == NULL){
    printf("read_info_file: No room for info\n");
    terminate(1);
  }
  info[0] = '\0';

  /* Node 0 reads the file */
  if(this_node==0)
    {
      /* Construct metadata file name from propagator file name 
	 by adding filename extension to propagator file name */
      strcpy(info_filename,base_filename);
      strcat(info_filename,ASCII_INFO_EXT);

      if((info_fp = fopen(info_filename,"r")) != NULL){
	buf = info;
	n = 0;
	while(fgets(buf, INFO_DATA_SIZE-n, info_fp)){
	  n = strlen(buf);
	  buf = info + n;
	}
      }
    }
      
  /* Node 0 broadcasts the result to all nodes */
  broadcast_bytes(info, strlen(info));

  return info;
}
