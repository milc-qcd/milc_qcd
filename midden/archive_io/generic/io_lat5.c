/*********************** io_lat4.c *************************/
/* MIMD version 5 */

/* routines for gauge configuration input/output. */
/* This works for most machines.  Wrappers for parallel I/O
   are in io_ansi.c, io_piofs.c, or io_paragon2.c */

/* Modifications */
/* 4/xx/00 additions to read archive format J.H. */
/* 4/17/98 r_parallel_w: g_syncs to prevent shmem message pileups C.D. */
/* 9/19/97 version 5 format with checksums C.D. */
/* 9/04/97 parallel files to be written in typewriter order C.D. */
/* 8/30/96 fixed macros for C syntax UMH */
/* 8/27/96 io_lat3.c converted to parallel reads and writes C.D. */
/*         Synchronization done through message passing instead of g_sync */
/*         Attempt at implementing ANSI standard, by UMH */

#include "generic_includes.h"
#include <io_lat.h>
#include <sys/types.h>
#include <fcntl.h>
#include <errno.h>
#include <string.h>
#include <time.h>

#define EPS 1e-6

#define PARALLEL 1   /* Must evaluate to true */
#define SERIAL 0     /* Must evaluate to false */

#define NODE_DUMP_ORDER 1
#define NATURAL_ORDER 0

#undef MAX_BUF_LENGTH
#define MAX_BUF_LENGTH 4096

#define SUCCESS  0
#define FAILURE -1
#define MAX_LINE_LENGTH 1024
#define MAX_TOKENS 512

typedef Real INPUT_TYPE;
typedef Real OUTPUT_TYPE;

/* Version: 1.0 */
#define OLDHEADERSIZE 0
#define TOL 0.0000001  
/* tolerance for Realing point checks */
/* For checksums we want a 32 bit unsigned int, for which      */
/* which we define a type uint32_t, which may eventually       */
/* become a C semi-standard.  Here we just check if INT_MAX    */
/* is 2^31-1; if not, we assume we are on a T3E or equivalent  */
/*=============================================================*/
#if (INT_MAX==2147483647)
#define INTS_ARE_32BIT
/*
  typedef unsigned int uint32_t;
*/
  typedef type32 uint32_t;
#else
#undef INTS_ARE_32BIT
  typedef unsigned short uint32_t;
#endif



/*----------------------------------------------------------------------
   Routines for archive I/O
   -----------------------------------------------------------------------*/

qcdhdr_get_str(char *s, QCDheader *hdr, char **q) {     
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
  
qcdhdr_get_int(char *s,QCDheader *hdr,int *q) {
  char *p;
  qcdhdr_get_str(s,hdr,&p);
  if (p==NULL) return (FAILURE);
  sscanf(p,"%d",q);
  return (SUCCESS);
}
qcdhdr_get_Real(char *s, QCDheader *hdr, Real *q) {
  char *p;
  qcdhdr_get_str(s,hdr,&p);
  if (p==NULL) return (FAILURE);
  sscanf(p,"%lfHELP",q);
  return (SUCCESS);
}

error_exit(char *s) { printf("%s\n",s); terminate(1);}

complete_U(double *u) {
  u[12] = u[ 2]*u[10] - u[ 4]*u[ 8] - u[ 3]*u[11] + u[ 5]*u[ 9];
  u[13] = u[ 4]*u[ 9] - u[ 2]*u[11] + u[ 5]*u[ 8] - u[ 3]*u[10];
  u[14] = u[ 4]*u[ 6] - u[ 0]*u[10] - u[ 5]*u[ 7] + u[ 1]*u[11];
  u[15] = u[ 0]*u[11] - u[ 4]*u[ 7] + u[ 1]*u[10] - u[ 5]*u[ 6];
  u[16] = u[ 0]*u[ 8] - u[ 2]*u[ 6] - u[ 1]*u[ 9] + u[ 3]*u[ 7];
  u[17] = u[ 2]*u[ 7] - u[ 0]*u[ 9] + u[ 3]*u[ 6] - u[ 1]*u[ 8];
}


big_endian() {
  union  {
    long l;
    char c[sizeof (long)];
  } u;
  u.l = 1;
  return (u.c[sizeof (long) - 1] == 1);
}

void byte_swap(int n, type32 *u);
/*
byte_swap(int n, type32 *u) {
  union {
    type32 uint_wrd;
    char c[sizeof(type32)];
  } wrd;
  register char chr;
  int s;
  for (s=0;s<n;s++) {
    wrd.uint_wrd = u[s];
    chr = wrd.c[0];
    wrd.c[0] = wrd.c[3];
    wrd.c[3] = chr;
    chr = wrd.c[2];
    wrd.c[2] = wrd.c[1];
    wrd.c[1] = chr;
    u[s] = wrd.uint_wrd;
  }
}
*/

void byte_swap(int n, type32 w[])
{
  register type32 old,new;
  int j;

  for(j=0; j<n; j++)
    {
      old = w[j];
      new = old >> 24 & 0x000000ff;
      new |= old >> 8 & 0x0000ff00;
      new |= old << 8 & 0x00ff0000;
      new |= old << 24 & 0xff000000;
      w[j] = new;
    }
} /* byterevn */



QCDheader * qcdhdr_get_hdr(in)
     FILE *in;
{
#define MAX_LINE_LENGTH 1024
#define MAX_TOKENS 512
  char line[MAX_LINE_LENGTH];
  int i,n,len;
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
  while (1) {
    fgets(line,MAX_LINE_LENGTH,in);
    if (strcmp(line,"END_HEADER\n")==0) break;

    printf("%s\n", line);

    /* Tokens are terminated by a space */
    q = index(line,' ');

    /* Overwrite space with a terminating null */
    *q = '\0';
    len = (int) q - (int) line;

    /* allocate space and copy the token in to it */
    p = malloc(len+1);
    (*hdr).token[n] = p;
    strcpy(p,line);

    q = index(++q,'='); q++;
    len = strlen(q);
    q[len-1] = 0;
    p = malloc(len);
    (*hdr).value[n] = p;
    strcpy(p,q);
    n++;
  }
  (*hdr).ntoken = n;
  return (hdr);
}





/*----------------------------------------------------------------------*/

/* For doing byte reversal on 32-bit words */

void byterevn(type32 w[], int n)
{
  register type32 old,new;
  int j;
  
  for(j=0; j<n; j++)
    {
      old = w[j];
      new = old >> 24 & 0x000000ff;
      new |= old >> 8 & 0x0000ff00;
      new |= old << 8 & 0x00ff0000;
      new |= old << 24 & 0xff000000;
      w[j] = new;
    }
} /* byterevn */

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
    byterevn((type32 *)src,size/sizeof(type32));
  return status;
}
/*---------------------------------------------------------------------------*/
int sread_byteorder(int byterevflag, FILE* fp, void *src, size_t size, char *myname, char *descrip)
{
  int status;

  status = sread_data(fp,src,size,myname,descrip);
  if(byterevflag==1)
    byterevn((type32 *)src,size/sizeof(type32));
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
		       void *src,       /* address of starting data
					   Realing point data must be
					   of type (Real) */
		       int count,       /* number of data items if > 1 */
		       int stride)      /* byte stride of data if
                                           count > 1 */
{

  int i,k,n;
  char *data;

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
      else if(strstr(fmt,"e") != NULL)
	fprintf(fpout,fmt,(double)(*(Real *)data));
      else if(strstr(fmt,"f") != NULL)
	fprintf(fpout,fmt,(double)(*(Real *)data));
      else if(strstr(fmt,"g") != NULL)
	fprintf(fpout,fmt,(double)(*(Real *)data));
      else
	{
	  printf("write_gauge_info_item: Unrecognized data type %s\n",fmt);
	  return 1;
	}
    }
  fprintf(fpout,"\n");
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

  write_gauge_info_item(info_fp,"magic_number","%d",&gh->magic_number,0,0);
  write_gauge_info_item(info_fp,"time_stamp","\"%s\"",gh->time_stamp,0,0);
  sprintf(sums,"%x %x",gf->check.sum29,gf->check.sum31);
  write_gauge_info_item(info_fp,"checksums","\"%s\"",sums,0,0);
  write_gauge_info_item(info_fp,"nx","%d",&nx,0,0);
  write_gauge_info_item(info_fp,"ny","%d",&ny,0,0);
  write_gauge_info_item(info_fp,"nz","%d",&nz,0,0);
  write_gauge_info_item(info_fp,"nt","%d",&nt,0,0);

  write_appl_gauge_info(info_fp);

  fclose(info_fp);

  printf("Wrote info file %s\n",info_filename);

} /*write_gauge_info_file */

/*----------------------------------------------------------------------*/

/* Set up the input gauge file and gauge header structures */

gauge_file *setup_input_gauge_file(char *filename)
{
  gauge_file *gf;
  gauge_header *gh;

  /* Allocate space for the file structure */

  gf = (gauge_file *)malloc(sizeof(gauge_file));
  if(gf == NULL)
    {
      printf("setup_input_gauge_file: Can't malloc gf\n");
      terminate(1);
    }

  gf->filename = filename;

  /* Allocate space for the header */

  gh = (gauge_header *)malloc(sizeof(gauge_header));
  if(gh == NULL)
    {
      printf("setup_input_gauge_file: Can't malloc gh\n");
      terminate(1);
    }

  gf->header = gh;
  gf->check.sum29 = 0;
  gf->check.sum31 = 0;

  return gf;
}

/*----------------------------------------------------------------------*/

/* Set up the output gauge file an gauge header structure */

gauge_file *setup_output_gauge_file()
{
  gauge_file *gf;
  gauge_header *gh;
  time_t time_stamp;
  int i;

  /* Allocate space for a new file structure */

  gf = (gauge_file *)malloc(sizeof(gauge_file));
  if(gf == NULL)
    {
      printf("setup_gauge_header: Can't malloc gf\n");
      terminate(1);
    }

  /* Allocate space for a new header structure */

  gh = (gauge_header *)malloc(sizeof(gauge_header));
  if(gh == NULL)
    {
      printf("setup_gauge_header: Can't malloc gh\n");
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
      for(i = strlen(gh->time_stamp) + 1; i < sizeof(gh->time_stamp); i++)
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
	  printf("w_serial_i: Node %d can't open file %s, error %d\n",
		 this_node,filename,errno);fflush(stdout);
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
  type32 *val;
  int rank29,rank31;
  su3_matrix *lbuf;
  su3_matrix tbuf[4];
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

      lbuf = (su3_matrix *)malloc(MAX_BUF_LENGTH*4*sizeof(su3_matrix));
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

      if( fseek(fp,offset,SEEK_SET) < 0 ) 
	{
	  printf("w_serial: Node %d fseek %ld failed error %d file %s\n",
		 this_node,(long)offset,errno,gf->filename);
	  fflush(stdout);terminate(1);
	}
    }
      
  /* Buffered algorithm for writing fields in serial order */
  
  /* initialize checksums */
  gf->check.sum31 = 0;
  gf->check.sum29 = 0;
  /* counts 32-bit words mod 29 and mod 31 in order of appearance on file */
  /* Here only node 0 uses these values */
  rank29 = 4*sizeof(su3_matrix)/sizeof(type32)*sites_on_node*this_node % 29;
  rank31 = 4*sizeof(su3_matrix)/sizeof(type32)*sites_on_node*this_node % 31;

  g_sync();
  currentnode=0;

  buf_length = 0;

  for(j=0,t=0;t<nt;t++)for(z=0;z<nz;z++)for(y=0;y<ny;y++)for(x=0;x<nx;x++,j++)
    {
      newnode=node_number(x,y,z,t);
      if(newnode != currentnode){	/* switch to another node */
	/* Send a few bytes to tell newnode it's OK to send */
	if( this_node==0 && newnode!=0 )send_field((char *)tbuf,4,newnode);
	if( this_node==newnode && newnode!=0 )get_field((char *)tbuf,4);
	currentnode=newnode;
      }
      
      if(this_node==0)
	{
	  if(currentnode==0)
	    {
	      i=node_index(x,y,z,t);
	      for(k=0;k<4;k++)tbuf[k] = lattice[i].link[k];
	    }
	  else
	    {
	      get_field((char *)tbuf,4*sizeof(su3_matrix));
	    }

	  memcpy((void *)&lbuf[4*buf_length], 
		 (void *)tbuf, 4*sizeof(su3_matrix));


	  /* Accumulate checksums - contribution from next site */
	  for(k = 0, val = (type32 *)&lbuf[4*buf_length]; 
	      k < 4*sizeof(su3_matrix)/sizeof(type32); k++, val++)
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
	      
	      if( fwrite(lbuf,4*sizeof(su3_matrix),buf_length,fp) != buf_length)
		{
		  printf("w_serial: Node %d gauge configuration write error %d file %s\n",
			 this_node,errno,gf->filename); 
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
	    send_field((char *)lattice[i].link,4*sizeof(su3_matrix),0);
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
      if( fseek(fp,checksum_offset,SEEK_SET) < 0 ) 
	{
	  printf("w_serial: Node %d fseek %ld failed error %d file %s\n",
		 this_node,(long)offset,errno,gf->filename);
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
      gf->rank2rcv = (type32 *)malloc(volume*sizeof(type32));
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
	      if(g_read(gf->rank2rcv,sizeof(type32),volume,gf->fp) != volume )
		{
		  printf("read_site_list: Node %d site list read error %d\n",
			 this_node,errno);
		  terminate(1);	
		}
	    }
	  else
	    {
	      if(fread(gf->rank2rcv,sizeof(type32),volume,gf->fp) != volume )
		{
		  printf("read_site_list: Node %d site list read error %d\n",
			 this_node,errno);
		  terminate(1);	
		}
	    }
	  
	  if(gf->byterevflag==1)byterevn(gf->rank2rcv,volume);
	}

      /* Broadcast result to all nodes */

      broadcast_bytes((char *)gf->rank2rcv,volume*sizeof(type32));
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
  type32 tmp;
  int j;
  int sixtyfourbits;
  Real c1,c2;
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
      byterevn((type32 *)&gh->magic_number,1);
      if(gh->magic_number == GAUGE_VERSION_NUMBER_V1) 
	{
	  *byterevflag=1;
	  printf("Reading as old-style gauge field configuration with byte reversal\n");
	  if( sizeof(Real) != sizeof(type32)) {
	    printf("read_v3_gauge_hdr: Can't byte reverse\n");
	    printf("requires size of type32(%d) = size of Real(%d)\n",
		   (int)sizeof(type32),(int)sizeof(Real));
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
  
  if(psread_byteorder(*byterevflag,parallel,fp,&c1,sizeof(Real),
		   myname,"c1")!=0)terminate(1);
  if(psread_byteorder(*byterevflag,parallel,fp,&c2,sizeof(Real),
		   myname,"c2")!=0)terminate(1);

  printf("Old format header parameters are %f %f\n",c1,c2);
  
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
  type32 tmp;
  int j;
  /* We keep this part of the old gauge header, but
     we ignore all but the two parameters */

  struct {                      /* Gauge field parameters */
    type32 n_descript;          /* Number of bytes in character string */
    char   descript[MAX_GAUGE_FIELD_DESCRIPT];  /* Describes gauge field */
    type32 n_param;             /* Number of gauge field parameters */
    Real  param[MAX_GAUGE_FIELD_PARAM];        /* GF parameters */
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
      byterevn((type32 *)&gh->magic_number,1);
      if(gh->magic_number == GAUGE_VERSION_NUMBER_1996) 
	{
	  *byterevflag=1;
	  printf("Reading as 1996-style gauge field configuration with byte reversal\n");
	  if( sizeof(Real) != sizeof(type32)) {
	    printf("read_1996_gauge_hdr: Can't byte reverse\n");
	    printf("requires size of type32(%d) = size of Real(%d)\n",
		   (int)sizeof(type32),(int)sizeof(Real));
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

int read_gauge_hdr(gauge_file *gf, int parallel) 
{
  /* parallel = 1 (TRUE) if all nodes are accessing the file */
  /*            0        for access from node 0 only */

  FILE *fp;
  gauge_header *gh;
  type32 tmp, btmp;
  int j;
  int byterevflag;
  char myname[] = "read_gauge_hdr";
  /* for Archive format */
  char line[MAX_LINE_LENGTH];
  int i,n,len;
  QCDheader *hdr;
  char **tokens, **values;
  char *p, *q;
  char *str;
  int dims[4];
  int ARCHYES=0;

  fp = gf->fp;
  gh = gf->header;

  /* Read and verify magic number */

  if(psread_data(parallel, fp,&gh->magic_number,sizeof(gh->magic_number),
			 myname,"magic number")!=0)terminate(1);

  tmp = gh->magic_number;
  byterevn((type32 *)&gh->magic_number,1);
  btmp = gh->magic_number;
  byterevn((type32 *)&gh->magic_number,1);


  /** See if header chunk is BEGI = 1229407554 **/

  if(tmp == GAUGE_VERSION_NUMBER_ARCHIVE) 
    {
      printf("reading Archive format\n"); 
      ARCHYES=1;
      byterevflag=0;
    }
  else if(btmp == GAUGE_VERSION_NUMBER_ARCHIVE) 
	{
	  printf("reading Archive format with byte reversal\n"); 
	  ARCHYES=1;
	  byterevflag=1;
	  if( sizeof(Real) != sizeof(type32)) {
	    printf("read_gauge_hdr: Can't byte reverse\n");
	    printf("requires size of type32(%d) = size of Real(%d)\n",
		   (int)sizeof(type32),(int)sizeof(Real));
	    terminate(1);
	  }
	}
  else if(gh->magic_number == GAUGE_VERSION_NUMBER) 
    {
      byterevflag=0;
    }
  else 
    {
      byterevn((type32 *)&gh->magic_number,1);
      if(gh->magic_number == GAUGE_VERSION_NUMBER) 
	{
	  byterevflag=1;
	  printf("Reading with byte reversal\n");
	  if( sizeof(Real) != sizeof(type32)) {
	    printf("read_gauge_hdr: Can't byte reverse\n");
	    printf("requires size of type32(%d) = size of Real(%d)\n",
		   (int)sizeof(type32),(int)sizeof(Real));
	    terminate(1);
	  }
	}
      else
	{
	  /* Restore magic number as originally read */
	  gh->magic_number = tmp;

	  /* Try old-style configurations */
	  if((read_v3_gauge_hdr(gf,parallel,&byterevflag) != 0) &&
	     (read_1996_gauge_hdr(gf,parallel,&byterevflag) != 0))
	    {
	      /* End of the road. */
	      printf("read_gauge_hdr: Unrecognized magic number in gauge configuration file header.\n");
	      printf("Expected %x but read %x\n",
		     GAUGE_VERSION_NUMBER,tmp);
	      printf("Expected %s but read %s\n",
		     GAUGE_VERSION_NUMBER,tmp);
	      terminate(1);
	    }
	  return byterevflag;
	}
    }
  
  /* Read header, do byte reversal, 
     if necessary, and check consistency */
  
  /* Lattice dimensions */


  if(ARCHYES == 1) 
    {

      gf->header->order == NATURAL_ORDER;

  /* MUST BE FIXED FOR bytereversal !!!!! */

    fgets(line,MAX_LINE_LENGTH,fp);

    /*
    printf("sizeof(char) = %d\n", sizeof(char));
    printf("sizeof(type32) = %d\n", sizeof(type32));
    */

    /*
    if (strcmp(line,"N_HEADER\n")!=0) {
      printf("reading Archive format: Missing \"BEGIN_HEADER\"; punting...\n");
      terminate(1);
    }
    */
    /* Allocate space for QCDheader and its pointers */
    tokens = (char **) malloc(MAX_TOKENS*sizeof(char *));
    values = (char **) malloc(MAX_TOKENS*sizeof(char *));
    hdr = (QCDheader *)malloc(sizeof(QCDheader));
    (*hdr).token = tokens;
    (*hdr).value = values;
    
    /* Begin loop on tokens */
    n = 0;
    printf("reading Archive header:\n");
    while (1) {
      fgets(line,MAX_LINE_LENGTH,fp);
      
      printf("%s", line);
      
      if (strcmp(line,"END_HEADER\n")==0) break;
      
      /* Tokens are terminated by a space */
      q = index(line,' ');
      
      /* Overwrite space with a terminating null */
      *q = '\0';
      len = (int) q - (int) line;
      
      /* allocate space and copy the token in to it */
      p = malloc(len+1);
      (*hdr).token[n] = p;
      strcpy(p,line);
      
      q = index(++q,'='); q++;
      len = strlen(q);
      q[len-1] = 0;
      p = malloc(len);
      (*hdr).value[n] = p;
      strcpy(p,q);
      n++;
    }
    (*hdr).ntoken = n;
    
    /* Get dimensions */
    if (qcdhdr_get_int("DIMENSION_1",hdr,dims+0)==FAILURE) error_exit("DIMENSION_1 not pr
esent");
    if (qcdhdr_get_int("DIMENSION_2",hdr,dims+1)==FAILURE) error_exit("DIMENSION_2 not pr
esent");
    if (qcdhdr_get_int("DIMENSION_3",hdr,dims+2)==FAILURE) error_exit("DIMENSION_3 not pr
esent");
    if (qcdhdr_get_int("DIMENSION_4",hdr,dims+3)==FAILURE) error_exit("DIMENSION_4 not pr
esent");

    for(i=0; i<4; i++) gh->dims[i] = dims[i];

    }
/* not a Archive lattice */
  else
    {

    if(psread_byteorder(byterevflag,parallel,fp,gh->dims,sizeof(gh->dims),
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
	  printf("read_gauge_hdr: Incorrect lattice dimensions ");
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


  if(!ARCHYES) {

    /* Date and time stamp */
    
    if(psread_data(parallel,fp,gh->time_stamp,sizeof(gh->time_stamp),
		   myname,"time stamp")!=0)terminate(1);
    
    /* Header byte length */
    
    gh->header_bytes = sizeof(gh->magic_number) + sizeof(gh->dims) + 
      sizeof(gh->time_stamp) + sizeof(gh->order);
    
    /* Data order */
    
    if(psread_byteorder(byterevflag,parallel,fp,&gh->order,sizeof(gh->order),
			myname,"order parameter")!=0)terminate(1);
    
    return byterevflag;
  }  
  
} /* read_gauge_hdr */

/*---------------------------------------------------------------------------*/

gauge_file *r_serial_i(char *filename)
{
  /* Returns file descriptor for opened file */

  gauge_header *gh;
  gauge_file *gf;
  FILE *fp;
  int byterevflag;

  /* All nodes set up a gauge file and gauge header structure for reading */

  gf = setup_input_gauge_file(filename);
  gh = gf->header;

  /* File opened for serial reading */
  gf->parallel = 0;

  /* Node 0 alone opens the file and reads the header */

  if(this_node==0)
    {
      fp = fopen(filename, "rb");
      if(fp == NULL)
	{
	  printf("r_serial_i: Node %d can't open file %s, error %d\n",
		 this_node,filename,errno);fflush(stdout);terminate(1);
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
  type32 *val;
  int rank29,rank31;
  su3_matrix *lbuf;
  char myname[] = "r_serial";

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
      else coord_list_size = sizeof(type32)*volume;
      checksum_offset = gf->header->header_bytes + coord_list_size;
      head_size = checksum_offset + gauge_check_size;
      
      /* Allocate space for read buffer */

      if(gf->parallel)
	printf("%s: Attempting serial read from parallel file \n",myname);

      lbuf = (su3_matrix *)malloc(MAX_BUF_LENGTH*4*sizeof(su3_matrix));
      if(lbuf == NULL)
	{
	  printf("%s: Node %d can't malloc lbuf\n",myname,this_node);
	  fflush(stdout);
	  terminate(1);
	}
  
      /* Position file for reading gauge configuration */
      
      offset = head_size;

      if( fseek(fp,offset,SEEK_SET) < 0 ) 
	{
	  printf("%s: Node 0 fseek %ld failed error %d file %s\n",
		 myname,(long)offset,errno,filename);
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
	    
	    if( fread(lbuf,4*sizeof(su3_matrix),buf_length,fp) != buf_length)
	      {
		printf("%s: node %d gauge configuration read error %d file %s\n",
		       myname,this_node,errno,filename); 
		fflush(stdout); terminate(1);
	      }
	    where_in_buf = 0;  /* reset counter */
	  }  /*** end of the buffer read ****/

	if(destnode==0){	/* just copy links */
	  i = node_index(x,y,z,t);
	  memcpy((void *)&lattice[i].link[0], 
		 (void *)&lbuf[4*where_in_buf], 4*sizeof(su3_matrix));
	}
	else {		/* send to correct node */
	  send_field((char *)&lbuf[4*where_in_buf],
		     4*sizeof(su3_matrix),destnode);
	}
	where_in_buf++;
      }
      
      /* The node which contains this site reads message */
      else {	/* for all nodes other than node 0 */
	if(this_node==destnode){
	  i = node_index(x,y,z,t);
	  get_field((char *)&lattice[i].link[0],4*sizeof(su3_matrix));
	}
      }

      /* The receiving node does the byte reversal and then checksum,
         if needed */

      if(this_node==destnode)
	{
	  if(byterevflag==1)
	    byterevn((type32 *)&lattice[i].link[0],
		     4*sizeof(su3_matrix)/sizeof(type32));
	  /* Accumulate checksums */
	  for(k = 0, val = (type32 *)&lattice[i].link[0]; 
	      k < 4*sizeof(su3_matrix)/sizeof(type32); k++, val++)
	    {
	      test_gc.sum29 ^= (*val)<<rank29 | (*val)>>(32-rank29);
	      test_gc.sum31 ^= (*val)<<rank31 | (*val)>>(32-rank31);
	      rank29++; if(rank29 >= 29)rank29 = 0;
	      rank31++; if(rank31 >= 31)rank31 = 0;
	    }
	}
      else
	{
	  rank29 += 4*sizeof(su3_matrix)/sizeof(type32);
	  rank31 += 4*sizeof(su3_matrix)/sizeof(type32);
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
	  if( fseek(fp,checksum_offset,SEEK_SET) < 0 ) 
	    {
	      printf("%s: Node 0 fseek %ld failed error %d file %s\n",
		    myname,(long)offset,errno,filename);
	      fflush(stdout);terminate(1);   
	    }
	  read_checksum(SERIAL,gf,&test_gc);
	}
      fflush(stdout);
      free(lbuf);
    }
  
} /* r_serial */

void r_serial_arch(gauge_file *gf)
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
  gauge_check test_gc;
  type32 *val;
  int rank29,rank31;
  su3_matrix tmpsu3[4];
  char myname[] = "r_serial";

  int vol3,s,mu,j,a,b,p;
  Real *uin, *q;
  int big_end;
  double U[4][18];
  type32 *chk;
  
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
      uin = (Real *) malloc(nx*ny*nz*48*sizeof(Real));
      if(uin == NULL)
	{
	  printf("%s: Node %d can't malloc uin buffer to read timeslice\n",
		 myname,this_node);
	  printf("recompile with smaller read buffer: uin\n");
	  fflush(stdout);
	  terminate(1);
	}
    }

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
	fread(uin,48*sizeof(Real),1,fp);
	if (!big_end) byte_swap(48,uin);
	q = uin;
	for (mu=0;mu<4;mu++) {
	  for (p=0;p<12;p++) {
	    j += *(type32 *) q;
	    U[mu][p] = (double) *q++;
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

	if(destnode==0){	/* just copy links */
	  i = node_index(x,y,z,t);
     /*   printf("lattice node_index = %d, mu = %d\n", i, mu); */
	  memcpy((void *)&lattice[i].link[0], 
		 (void *)&tmpsu3, 4*sizeof(su3_matrix));
	} else {		/* send to correct node */
	  send_field((char *)&tmpsu3, 4*sizeof(su3_matrix),destnode);
	}
      } 
      /* The node which contains this site reads message */
      else {	/* for all nodes other than node 0 */
	if(this_node==destnode){
	  i = node_index(x,y,z,t);
	  get_field((char *)&lattice[i].link[0],4*sizeof(su3_matrix));
	}
      }

      /* The receiving node does the byte reversal and then checksum,
         if needed */

      if(this_node==destnode)
	{
	  if(byterevflag==1)
	    byterevn((type32 *)&lattice[i].link[0],
		     4*sizeof(su3_matrix)/sizeof(type32));
	  /* Accumulate checksums */
	  for(k = 0, val = (type32 *)&lattice[i].link[0]; 
	      k < 4*sizeof(su3_matrix)/sizeof(type32); k++, val++)
   	    {
	      test_gc.sum29 ^= (*val)<<rank29 | (*val)>>(32-rank29);
	      test_gc.sum31 ^= (*val)<<rank31 | (*val)>>(32-rank31);
	      rank29++; if(rank29 >= 29)rank29 = 0;
	      rank31++; if(rank31 >= 31)rank31 = 0;
	    }
	}
      else
	{
	  rank29 += 4*sizeof(su3_matrix)/sizeof(type32);
	  rank31 += 4*sizeof(su3_matrix)/sizeof(type32);
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
	  if( fseek(fp,checksum_offset,SEEK_SET) < 0 ) 
	    {
	      printf("%s: Node 0 fseek %ld failed error %d file %s\n",
		    myname,(long)offset,errno,filename);
	      fflush(stdout);terminate(1);   
	    }
	  read_checksum(SERIAL,gf,&test_gc);
	}
      fflush(stdout);
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
  type32 coords, *cbuf;

  /* All nodes write their site coordinate list in sequential
     blocks after the header.  The list is in the order of appearance
     in the lattice array.  Node 0 writes to the first block
     followed by node 1, etc.  The result is a contiguous table
     that can be used to remap the data to the corresponding space-time
     coordinate */

  /* Location of site list for this node */
  
  offset = gh->header_bytes + 
    sizeof(type32)*sites_on_node*this_node;

  cbuf = (type32 *)malloc(sites_on_node*sizeof(type32));
  if(cbuf == NULL)
    {
      printf("write_site_list: node %d can't malloc cbuf\n",this_node);
      fflush(stdout);terminate(1);   
    }

  if( g_seek(fp,offset,SEEK_SET) < 0 ) 
    {
      printf("write_site_list: node %d fseek %ld failed errno %d\n",
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

    if( g_write(cbuf,sizeof(type32),sites_on_node,fp) != sites_on_node)
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
  
  /* All nodes write site list to file */

  write_site_list(fp,gh);
  
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

su3_matrix *w_parallel_setup(gauge_file *gf, off_t *checksum_offset)
{
  /* gf  = file descriptor as opened by w_checkpoint_i */

  FILE *fp;
  gauge_header *gh;
  su3_matrix *lbuf;

  off_t offset ;           /* File stream pointer */
  off_t gauge_node_size;   /* Size of a gauge configuration block for
                              all sites on one node */
  off_t coord_list_size;   /* Size of coordinate list in bytes */
  off_t head_size;         /* Size of header plus coordinate list */
  off_t gauge_check_size;  /* Size of checksum */
  char myname[] = "w_parallel_setup";

  if(!gf->parallel)
    printf("%s: Attempting parallel write to serial file.\n",myname);

  lbuf = (su3_matrix *)malloc(MAX_BUF_LENGTH*4*sizeof(su3_matrix));
  if(lbuf == NULL)
    {
      printf("%s: Node %d can't malloc lbuf\n",myname,this_node);
      fflush(stdout);
      terminate(1);
    }

  fp = gf->fp;
  gh = gf->header;

  gauge_node_size = sites_on_node*4*sizeof(su3_matrix) ;

  if(gf->header->order == NATURAL_ORDER)coord_list_size = 0;
  else coord_list_size = sizeof(type32)*volume;
  head_size = gf->header->header_bytes + coord_list_size;
  *checksum_offset = head_size;
  gauge_check_size = sizeof(gf->check.sum29) + sizeof(gf->check.sum31);

  offset = head_size + gauge_check_size;

  /* Each node writes its gauge configuration values */

  offset += gauge_node_size*this_node;
  
  if( g_seek(fp,offset,SEEK_SET) < 0 ) 
    {
      printf("%s: Node %d fseek %ld failed error %d file %s\n",
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
  su3_matrix *lbuf;
  int buf_length,where_in_buf;
  type32 *val;
  int rank29,rank31;
  off_t checksum_offset;
  register int i;
  int j,k;
  int x,y,z,t;
  struct {
    short x,y,z,t;
    su3_matrix link[4];
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
	      memcpy((void *)msg.link,
		     (void *)lattice[i].link, 4*sizeof(su3_matrix));
	      
	      send_field((char *)&msg,sizeof(msg),destnode);
	    }
	    /* Node destnode receives a message */
	    else if(this_node==destnode){
	      if(destnode==sendnode){ /* just copy links */
		i = node_index(x,y,z,t);
		where_in_buf = buf_length;
		memcpy((void *)&lbuf[4*where_in_buf],
		       (void *)lattice[i].link,4*sizeof(su3_matrix));
		rank29 = rank31 = 
		  4*sizeof(su3_matrix)/sizeof(type32)*rcv_rank;
	      }
	      else {
		/* Receive a message */
		/* Note that messages may arrive in any order
		   so we use the x,y,z,t coordinate to tell
		   where it goes in the write buffer */
		get_field((char *)&msg,sizeof(msg));
		/* Reconstruct rank from message coordinates */
		i = msg.x+nx*(msg.y+ny*(msg.z+nz*msg.t));
		/* The buffer location is then */
		where_in_buf = (i % sites_on_node) % MAX_BUF_LENGTH;

		/* Move data to buffer */
		memcpy((void *)&lbuf[4*where_in_buf],
		       (void *)msg.link,4*sizeof(su3_matrix));
		rank29 = rank31 = 4*sizeof(su3_matrix)/sizeof(type32)*i;
	      }

	      /* Receiving node accumulates checksums as the values
		 are inserted into its buffer */
	      rank29 %= 29; rank31 %= 31;
	      for(k = 0, val = (type32 *)&lbuf[4*where_in_buf]; 
		  k < 4*sizeof(su3_matrix)/sizeof(type32); k++, val++)
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
		  
		  if( g_write(lbuf,4*sizeof(su3_matrix),buf_length,fp) 
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

  if(this_node == 0)
    {
      printf("Saved gauge configuration in parallel to binary file %s\n",
	     gf->filename);
      printf("Time stamp %s\n",(gf->header)->time_stamp);

      /* Node 0 positions file for writing checksum */
      
      if( g_seek(fp,checksum_offset,SEEK_SET) < 0 ) 
	{
	  printf("%s: Node 0 g_seek %ld for checksum failed error %d file %s\n",
		 myname,(long)checksum_offset,errno,gf->filename);
	  fflush(stdout);terminate(1);   
	}
      
      write_checksum(PARALLEL,gf);
    }

} /* w_parallel */

/*-----------------------------------------------------------------------*/

/* Write parallel gauge configuration in node dump order */

void w_checkpoint(gauge_file *gf)
{
  /* gf  = file descriptor as opened by w_checkpoint_i */

  FILE *fp;
  su3_matrix *lbuf;
  type32 *val;
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
  rank29 = 4*sizeof(su3_matrix)/sizeof(type32)*sites_on_node*this_node % 29;
  rank31 = 4*sizeof(su3_matrix)/sizeof(type32)*sites_on_node*this_node % 31;

  buf_length = 0;

  FORALLSITES(i,s)
  {
        
    /* load the gauge configuration into the buffer */
    memcpy((void *)&lbuf[4*buf_length], 
	   (void *)lattice[i].link, 4*sizeof(su3_matrix));

    /* Accumulate checksums - contribution from next site moved into buffer*/
    for(k = 0, val = (type32 *)&lbuf[4*buf_length]; 
	k < 4*sizeof(su3_matrix)/sizeof(type32); k++, val++)
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
	if( g_write(lbuf,4*sizeof(su3_matrix),buf_length,fp) != buf_length)
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

  if(this_node == 0)
    {
      printf("Saved gauge configuration checkpoint file %s\n",
	     gf->filename);
      printf("Time stamp %s\n",(gf->header)->time_stamp);

      /* Node 0 positions file for writing checksum */
      
      if( g_seek(fp,checksum_offset,SEEK_SET) < 0 ) 
	{
	  printf("%s: Node 0 g_seek %ld for checksum failed error %d file %s\n",
		 myname,(long)checksum_offset,errno,gf->filename);
	  fflush(stdout);terminate(1);   
	}
      
      write_checksum(PARALLEL,gf);
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
  su3_matrix *lbuf;
  struct {
    short x,y,z,t;
    su3_matrix link[4];
  } msg;

  int buf_length,where_in_buf;
  gauge_check test_gc;
  type32 *val;
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

  /* Allocate read buffer */
  lbuf = (su3_matrix *)malloc(MAX_BUF_LENGTH*4*sizeof(su3_matrix));
  if(lbuf == NULL)
    {
      printf("%s: Node %d can't malloc lbuf\n",myname,this_node); 
      fflush(stdout);terminate(1);
    }

  gauge_node_size = sites_on_node*4*sizeof(su3_matrix) ;

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
  else coord_list_size = sizeof(type32)*volume;
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
  rank29 = 4*sizeof(su3_matrix)/sizeof(type32)*sites_on_node*this_node % 29;
  rank31 = 4*sizeof(su3_matrix)/sizeof(type32)*sites_on_node*this_node % 31;

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
		
		if( g_read(lbuf,buf_length*4*sizeof(su3_matrix),1,fp) != 1)
		  {
		    printf("%s: node %d gauge configuration read error %d file %s\n",
			   myname,this_node,errno,filename); 
		    fflush(stdout); terminate(1);
		  }
		where_in_buf = 0;  /* reset counter */
	      }  /*** end of the buffer read ****/
	    
	    /* Do byte reversal if needed */
	    if(gf->byterevflag==1)
	      byterevn((type32 *)&lbuf[4*where_in_buf],
		       4*sizeof(su3_matrix)/sizeof(type32));

	    /* Accumulate checksums - contribution from next site */
	    for(k = 0, val = (type32 *)&lbuf[4*where_in_buf]; 
		k < 4*sizeof(su3_matrix)/sizeof(type32); k++, val++)
	      {
		test_gc.sum29 ^= (*val)<<rank29 | (*val)>>(32-rank29);
		test_gc.sum31 ^= (*val)<<rank31 | (*val)>>(32-rank31);
		rank29++; if(rank29 >= 29)rank29 = 0;
		rank31++; if(rank31 >= 31)rank31 = 0;
	      }

	    if(destnode==sendnode){	/* just copy links */
	      i = node_index(x,y,z,t);
	      memcpy((void *)lattice[i].link,
		     (void *)&lbuf[4*where_in_buf],4*sizeof(su3_matrix));
	    }
	    else {		/* send to correct node */
	      /* Message consists of site coordinates and 4 link matrices */
	      msg.x = x; msg.y = y; msg.z = z; msg.t = t;
	      memcpy((void *)msg.link,
		     (void *)&lbuf[4*where_in_buf],4*sizeof(su3_matrix));
	      
	      send_field((char *)&msg,sizeof(msg),destnode);
	    }
	    where_in_buf++;
	  }
	  /* The node which contains this site reads a message */
	  else {	/* for all nodes other than node sendnode */
	    if(this_node==destnode){
	      get_field((char *)&msg,sizeof(msg));
	      i = node_index(msg.x,msg.y,msg.z,msg.t);
	      if(this_node!= node_number(msg.x,msg.y,msg.z,msg.t))
		{
		  printf("BOTCH. Node %d received %d %d %d %d\n",
			 this_node,msg.x,msg.y,msg.z,msg.t);
		  fflush(stdout); terminate(1);
		}
	      memcpy((void *)lattice[i].link,
		     (void *)msg.link,4*sizeof(su3_matrix));
	    }
	  }
	} /** end over the lattice sites in block on all nodes ***/

    g_sync(); /* To prevent incoming message pileups */
  }  /** end over blocks **/

  free(lbuf);

  /* Read and verify checksum */
  /* Checksums implemented with version 5 */
  
  /* Combine node checksum contributions with global exclusive or */
  g_xor32(&test_gc.sum29);
  g_xor32(&test_gc.sum31);

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
  su3_matrix lbuf[4];
  
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
	  if( (fscanf(fp,"%leHELP%leHELP\n",&(lbuf[dir].e[i][j].real),
		      &(lbuf[dir].e[i][j].imag)) )!= 2){
	    printf("restore_ascii: gauge link read error\n"); 
	    terminate(1);
	  }
	}
      }
      if(destnode==0){	/* just copy links */
	i = node_index(x,y,z,t);
	for(dir=XUP;dir<=TUP;dir++)lattice[i].link[dir]=lbuf[dir];
      }
      else {		/* send to correct node */
	send_field((char *)lbuf,4*sizeof(su3_matrix),destnode);
      }
    }
    
    /* The node which contains this site reads message */
    else {	/* for all nodes other than node 0 */
      if(this_node==destnode){
	get_field((char *)lbuf,4*sizeof(su3_matrix));
	i = node_index(x,y,z,t);
	for(dir=XUP;dir<=TUP;dir++)lattice[i].link[dir]=lbuf[dir];
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
  su3_matrix lbuf[4];
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
      if( this_node==newnode && newnode!=0 )get_field((char *)lbuf,4);
      currentnode=newnode;
    }
    
    if(this_node==0){
      if(currentnode==0){
	i=node_index(x,y,z,t);
	for(dir=XUP;dir<=TUP;dir++)lbuf[dir]=lattice[i].link[dir];
      }
      else{
	get_field((char *)lbuf,4*sizeof(su3_matrix));
      }
      for(dir=XUP;dir<=TUP;dir++){
	for(i=0;i<3;i++)for(j=0;j<3;j++){
	  if( (fprintf(fp,"%.7e\t%.7e\n",(double)lbuf[dir].e[i][j].real,
		       (double)lbuf[dir].e[i][j].imag))== EOF){
	    printf("Write error in save_ascii\n"); terminate(1);
	  }
	}
      }
    }
    else {	/* for nodes other than 0 */
      if(this_node==currentnode){
	i=node_index(x,y,z,t);
	for(dir=XUP;dir<=TUP;dir++)lbuf[dir]=lattice[i].link[dir];
	send_field((char *)lbuf,4*sizeof(su3_matrix),0);
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


/*---------------------------------------------------------------------*/
/* Restore lattice file by reading serially (node 0 only) */
/* Handles most lattice formats */
  
gauge_file *restore_serial(char *filename)
{
  gauge_file *gf;

  gf = r_serial_i(filename);
  if(gf->header->magic_number == GAUGE_VERSION_NUMBER_ARCHIVE) {
    r_serial_arch(gf);
  } else {
    r_serial(gf);
  }
  r_serial_f(gf);

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
/* Kept for compatibility */
/* Save lattice in natural order (version 3-4 format) by writing serially */

#if defined(T3D) || defined(T3E)
#define VERSION_NUMBER 59354
/* T3D and T3E version copied from version 4 io_t3d2.c */
gauge_file *save_old_binary(char *filename, Real c1, Real c2) {
int fd;
int currentnode,newnode;
int i,j,x,y,z,t,dir;
double dtime,dclock();
type32 stemp,dims[4];
su3_matrix lbuf[4];
char *writebuf;
Real x1,x2;
gauge_file *gf;

    /* node 0 does all the writing */
    dtime = -dclock();
    if(this_node==0){
        writebuf = (char *)malloc(nx*4*sizeof(su3_matrix));
        if(writebuf==NULL){
            printf("No room for write buffer\n"); terminate(1);
        }
        fd = creat(filename,0644);
        if(fd < 0){
	    printf("Can't open file %s, error %d\n",filename,errno);
	    fflush(stdout);terminate(1);
        }
	stemp=VERSION_NUMBER;
        if( (write(fd,&stemp,sizeof(type32))) != sizeof(type32) ){
	    printf("Error in writing lattice header\n"); terminate(1);
        }
	dims[XUP]=nx; dims[YUP]=ny; dims[ZUP]=nz; dims[TUP]=nt;
        if( (write(fd,dims,4*sizeof(type32))) != 4*sizeof(type32) ){
	    printf("Error in writing lattice header\n"); terminate(1);
        }
	x1 = c1; x2 = c2; /* can't take address of argument */
        if( (write(fd,&x1,sizeof(Real))) != sizeof(Real) ){
	    printf("Error in writing lattice header\n"); terminate(1);
        }
        if( (write(fd,&x2,sizeof(Real))) != sizeof(Real) ){
	    printf("Error in writing lattice header\n"); terminate(1);
        }
    }
    g_sync();
    currentnode=0;

    for(t=0;t<nt;t++)for(z=0;z<nz;z++)for(y=0;y<ny;y++){
      for(x=0;x<nx;x++){
	newnode=node_number(x,y,z,t);
	if(newnode != currentnode){	/* switch to another node */
	    /**g_sync();**/
	    /* tell newnode it's OK to send */
	    if( this_node==0 && newnode!=0 )send_field((char *)lbuf,4,newnode);
	    if( this_node==newnode && newnode!=0 )get_field((char *)lbuf,4);
	    currentnode=newnode;
	}

	if(this_node==0){
	    if(currentnode==0){
		i=node_index(x,y,z,t);
		for(dir=XUP;dir<=TUP;dir++)lbuf[dir]=lattice[i].link[dir];
	    }
	    else{
/*TEMP*//**pvm_initsend(PvmDataRaw );**/
/*TEMP*//**pvm_send(currentnode,OK_TYPE);**/
		get_field((char *)lbuf,4*sizeof(su3_matrix));
	    }
            memcpy(writebuf+4*x*sizeof(su3_matrix), lbuf,
                4*sizeof(su3_matrix));
	}
	else {	/* for nodes other than 0 */
	    if(this_node==currentnode){
		i=node_index(x,y,z,t);
		for(dir=XUP;dir<=TUP;dir++)lbuf[dir]=lattice[i].link[dir];
/*TEMP*//**pvm_recv(0,OK_TYPE);**/
		send_field((char *)lbuf,4*sizeof(su3_matrix),0);
	    }
	}
/**if(this_node==0){printf("Site %d %d %d %d done\n",x,y,z,t);fflush(stdout);}*/
      } /*close x, y and z loops */
      if(this_node==0){
        if( (write(fd,writebuf,nx*4*sizeof(su3_matrix)))
            != nx*4*sizeof(su3_matrix) ){
            printf("Write error in save_old_binary\n"); terminate(1);
        }
      }
    } /*close t loop */
    g_sync();
    if(this_node==0){
	printf("Saved gauge configuration in binary old-format file  %s\n",
	       filename);
        dtime += dclock();
    	printf("Time to save lattice %e\n",dtime);
        free(writebuf);
	close(fd);
	fflush(stdout);
    }

    /* Set up gauge file and gauge header structures and load header values */
    gf = setup_output_gauge_file();
    gf->fp = NULL;                     /* Not used */
    gf->filename = filename;
    gf->byterevflag    = 0;            /* Not used */
    gf->rank2rcv       = NULL;         /* Not used */
    gf->parallel       = 0;
    return gf;
}

/*---------------------------------------------------------------------------*/
#else
#define VERSION_NUMBER 59354
/* Non-T3D version copied from io_lat2.c */
gauge_file *save_old_binary(char *filename, Real c1, Real c2) 
{
  /* For compatibility.  Single node writes in old gauge file format */

int fd;
int currentnode,newnode;
int i,x,y,z,t,dir,dims[4];
su3_matrix lbuf[4];
Real x1,x2;
gauge_file *gf;

    /* node 0 does all the writing */
    if(this_node==0){
        fd = creat(filename,0644);
        if(fd < 0){
	    printf("Can't open file %s, error %d\n",filename,errno);terminate(1);
        }
	i=VERSION_NUMBER;
        if( (write(fd,&i,sizeof(int))) != sizeof(int) ){
	    printf("Error in writing lattice header\n"); terminate(1);
        }
	dims[XUP]=nx; dims[YUP]=ny; dims[ZUP]=nz; dims[TUP]=nt;
        if( (write(fd,dims,4*sizeof(int))) != 4*sizeof(int) ){
	    printf("Error in writing lattice header\n"); terminate(1);
        }
	x1 = c1; x2 = c2; /* can't take address of argument */
        if( (write(fd,&x1,sizeof(Real))) != sizeof(Real) ){
	    printf("Error in writing lattice header\n"); terminate(1);
        }
        if( (write(fd,&x2,sizeof(Real))) != sizeof(Real) ){
	    printf("Error in writing lattice header\n"); terminate(1);
        }
    }
    g_sync();
    currentnode=0;

    for(t=0;t<nt;t++)for(z=0;z<nz;z++)for(y=0;y<ny;y++)for(x=0;x<nx;x++){
	newnode=node_number(x,y,z,t);
	if(newnode != currentnode){	/* switch to another node */
	    /**g_sync();**/
	    /* tell newnode it's OK to send */
	    if( this_node==0 && newnode!=0 )send_field((char *)lbuf,4,newnode);
	    if( this_node==newnode && newnode!=0 )get_field((char *)lbuf,4);
	    currentnode=newnode;
	}

	if(this_node==0){
	    if(currentnode==0){
		i=node_index(x,y,z,t);
		for(dir=XUP;dir<=TUP;dir++)lbuf[dir]=lattice[i].link[dir];
	    }
	    else{
		get_field((char *)lbuf,4*sizeof(su3_matrix));
	    }
	    if( (write(fd,lbuf,4*sizeof(su3_matrix))) != 4*sizeof(su3_matrix) ){
		printf("save_old_binary: Write error %d file %s\n",
		       errno,filename); terminate(1);
	    }
	}
	else {	/* for nodes other than 0 */
	    if(this_node==currentnode){
		i=node_index(x,y,z,t);
		for(dir=XUP;dir<=TUP;dir++)lbuf[dir]=lattice[i].link[dir];
		send_field((char *)lbuf,4*sizeof(su3_matrix),0);
	    }
	}
    }
    g_sync();
    if(this_node==0){
	printf("Saved gauge configuration in binary old-format file  %s\n",
	       filename);
	close(fd);
	fflush(stdout);
    }

    /* Set up gauge file and gauge header structures and load header values */
    gf = setup_output_gauge_file();
    gf->fp = NULL;                     /* Not used */
    gf->filename = filename;
    gf->byterevflag    = 0;            /* Not used */
    gf->rank2rcv       = NULL;         /* Not used */
    gf->parallel       = 0;
    return gf;
}
#endif


gauge_file *save_archive(char *filename) {
  /* For compatibility.  Single node writes in old gauge file format */

  int currentnode,newnode;
  int i,j,x,y,z,t,dir,dims[4];
  su3_matrix lbuf[4];
  Real x1,x2;
  gauge_file *gf;
  gauge_header *gh;

  FILE *outfile, *headerfile;
  site *s;
  unsigned int chksum, *p32;
  /*  uint32_t chksum, *p32; */
  char *str, sums[30];
  time_t timestamp;
  INPUT_TYPE *uin;
  INPUT_TYPE buf;
  OUTPUT_TYPE *uout;
  int DIM1, DIM2, DIM3, DIM4;
  int big_end_p; 
  Real ssplaq,stplaq,avgtrace, avgplaq, tmpflt;
  double dtime, trace;
  int mu,a,b,vol3,nvol4, tslice;

  /* Check which end is up */
  big_end_p = big_endian();
  printf("big_end_p is %d\n", big_end_p);

  /* Set up gauge file and gauge header structures and load header values */
  gf = setup_output_gauge_file();
  gh = gf->header;
  
  /* node 0 does all the writing */
  if(this_node==0){
    
    trace = 0.0;
    chksum = 0;
    FORALLSITES(i,s) {
      for(mu=0; mu<4; ++mu) {
	trace += (trace_su3(&(s->link[mu]))).real;
	for(a=0; a<2; a++) for(b=0; b<3; b++) {
	  tmpflt = s->link[mu].e[a][b].real;
	  p32 = (unsigned int *) &tmpflt;
	  chksum += *p32;
	  tmpflt = s->link[mu].e[a][b].imag;
	  p32 = (unsigned int *) &tmpflt;
	  chksum += *p32;
	}
      }
    }
    avgtrace = trace/volume/12.0;
    printf("trace = %f\n", avgtrace);
    printf("chksum_x = %x\n", chksum);
    printf("chksum_u = %12u\n", chksum);
    plaquette(&ssplaq, &stplaq);
    avgplaq = (ssplaq+stplaq)/6.0;
    printf("plaquette = %f\n", avgplaq);
    
    printf("Writing archive format lattice to %s\n", savefile);
    /* Create output file */
    outfile = fopen(savefile,"w");
    if (outfile < 0) {
      printf("error opening output file: %s\n", savefile);
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
    fprintf(outfile,"ENSEMBLE_LABEL = %s\n", ensemble_label);
    fprintf(outfile,"ENSEMBLE_ID = %s\n", ensemble_id);
    fprintf(outfile,"SEQUENCE_NUMBER = %d\n",sequence_number);
    /* write Milc info section */
    fprintf(outfile,"MILC_INFO = -------BEGIN-------\n");
    write_gauge_info_item(outfile,"magic_number","%d",&gh->magic_number,0,0);
    write_gauge_info_item(outfile,"time_stamp","\"%s\"",gh->time_stamp,0,0);
    sprintf(sums,"%x %x",gf->check.sum29,gf->check.sum31);
    write_gauge_info_item(outfile,"checksums","\"%s\"",sums,0,0);
    write_gauge_info_item(outfile,"nx","%d",&nx,0,0);
    write_gauge_info_item(outfile,"ny","%d",&ny,0,0);
    write_gauge_info_item(outfile,"nz","%d",&nz,0,0);
    write_gauge_info_item(outfile,"nt","%d",&nt,0,0);
    write_appl_gauge_info(outfile);
    fprintf(outfile,"MILC_INFO = --------END--------\n");
    fprintf(outfile,"END_HEADER\n");

    /*  
    printf("BEGIN_HEADER\n");
    printf("HDR_VERSION = 1.0\n");
    printf("DATATYPE = 4D_SU3_GAUGE\n");
    printf("DIMENSION_1 = %d\n",nx);
    printf("DIMENSION_2 = %d\n",ny);
    printf("DIMENSION_3 = %d\n",nz);
    printf("DIMENSION_4 = %d\n",nt);
    printf("CHECKSUM = %x\n",chksum);
    printf("LINK_TRACE = %.10f\n",avgtrace);
    printf("PLAQUETTE = %.10f\n",avgplaq);
    printf("CREATOR = %s",
	   CREATOR
	   );
    printf("CREATOR_HARDWARE = %s",
	   CREATOR_HARDWARE
	   );
    timestamp = time(NULL);
    printf("ARCHIVE_DATE = %s",ctime(&timestamp));
    printf("ENSEMBLE_LABEL = %s",ensemble_label);
    printf("MILC_REALING_POINT = IEEE32\n");
    printf("ENSEMBLE_ID = %s",ensemble_id);
    printf("SEQUENCE_NUMBER = %d\n",sequence_number);
    printf("BETA = %f\n", beta);
    printf("MASS = %f\n", mass);
    printf("END_HEADER\n");
*/  
    
    vol3 = nx*ny*nz;
    uout = (OUTPUT_TYPE *) malloc(48*vol3*sizeof(OUTPUT_TYPE));
    if(uout == NULL) { 
      printf("can\'t malloc uout timeslice\n"); terminate(1); 
    }
    
    for(tslice=0; tslice<nt; ++tslice) {
      j = 0;
      for(z=0; z<nz; ++z) for(y=0; y<ny; ++y) for(x=0; x<nx; ++x) {
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
	++j;
      }
      if (!big_end_p) byte_swap(48*vol3,uout);
      if(fwrite(uout,48*vol3*sizeof(OUTPUT_TYPE),1,outfile) != 1)
	printf("fwrite bombed...\n");
      fflush(outfile);
    }
    
  }
  fclose(outfile);
  printf("Wrote archive gauge file %s\n",filename);

  ++sequence_number;
}







