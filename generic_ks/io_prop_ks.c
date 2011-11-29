/* io_prop_ks.c -- reads and writes KS quark propagators
   MIMD version 7
   MBW: 21 Feb 2001 -- just ASCII format for now, modified from
                      routines in ../generic/io_lat4.c
   MBW: Apr 2002 -- adding binary I/O
*/
/* This version assumes internal storage is at the prevailing
   precision, but the files are always 32 bit.  This code
   converts to and from the prevailing precision.  CD 11/29/04 */

#include "generic_ks_includes.h"
#include <sys/types.h>
#include <fcntl.h>
#include <errno.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include "../include/io_lat.h" /* for utilities like get_f ,etc */
#include "../include/io_ksprop.h"

#define PARALLEL 1
#define SERIAL 0

#define NATURAL_ORDER 0

#undef MAX_BUF_LENGTH
#define MAX_BUF_LENGTH 4096

/* The send buffer would normally be the size of a single SU(3) vector
   but we may need to pad for short messages generated in serial reads
   and writes to avoid switch inefficiencies.  In that case we define
   PAD_SEND_BUF on the compilation line to increase this.  */
   
#ifndef PAD_SEND_BUF
#define PAD_SEND_BUF 8
#endif

/*---------------------------------------------------------------------------*/
/* Convert (or copy) one single precision su3_vector to generic precision */

void f2d_vec(fsu3_vector *a, su3_vector *b){
  int i;
  
  for(i = 0; i < 3; i++){
      b->c[i].real = a->c[i].real;
      b->c[i].imag = a->c[i].imag;
  }
}

/* Convert (or copy) one generic precision su3_vector to single precision */
void d2f_vec(su3_vector *a, fsu3_vector *b){
  int i;
  
  for(i = 0; i < 3; i++){
    b->c[i].real = a->c[i].real;
    b->c[i].imag = a->c[i].imag;
  }
}

/*------------------------------------------------------------------------*/
/* This stuff is pretty identical to both gauge and wilson prop routines.
   Oh well, don't want to mess with them, so here are the new ones for me.
*/

/*----------------------------------------------------------------------*/
/* This subroutine writes the propagator header structure */
/* Serial access version */

void swrite_ks_prop_hdr(FILE *fp, ks_prop_header *ksph)
{
  char myname[] = "swrite_ks_prop_hdr";

  swrite_data(fp,(void *)&ksph->magic_number,sizeof(ksph->magic_number),
	      myname,"magic_number");
  swrite_data(fp,(void *)ksph->dims,sizeof(ksph->dims),
	      myname,"dimensions");
  swrite_data(fp,(void *)ksph->time_stamp,sizeof(ksph->time_stamp),
	      myname,"time_stamp");
  swrite_data(fp,&ksph->order,sizeof(ksph->order),
	      myname,"order");    

  /* Header byte length */

  ksph->header_bytes = sizeof(ksph->magic_number) + sizeof(ksph->dims) + 
    sizeof(ksph->time_stamp) + sizeof(ksph->order);

} /* swrite_ks_prop_hdr */

/*------------------------------------------------------------------------*/
/* Write a data item to the prop info file */
int write_ksprop_info_item( FILE *fpout,    /* ascii file pointer */
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
#define MAXITEM 128
  char item[MAXITEM];

  /* Check for valid keyword */

  for(i=0;strlen(ks_prop_info_keyword[i])>0 &&
      strcmp(ks_prop_info_keyword[i],keyword) != 0; i++);
  /**  if(strlen(ks_prop_info_keyword[i])==0)
    printf("write_ksprop_info_item: WARNING: keyword %s not in table\n",
    keyword); **/

  /* Write keyword */
  /* All fprintf's converted to fwrite's to humor dcache utilities at FNAL */

  fwrite(keyword,strlen(keyword),1,fpout);
  /*  fprintf(fpout,"%s",keyword); */

  /* Write count if more than one item */
  if(count > 1){
    sprintf(item,"[%d]",count);
    fwrite(item,strlen(item),1,fpout);
    /*fprintf(fpout,"[%d]",count);  */
  }

  n = count; if(n==0)n = 1;
  
  /* Write data */
  for(k = 0, data = (char *)src; k < n; k++, data += stride)
    {
      fwrite(" ",1,1,fpout);
      /*fprintf(fpout," ");*/
      if(strstr(fmt,"s") != NULL){
	sprintf(item,fmt,data);
	fwrite(item,strlen(item),1,fpout);
	/*fprintf(fpout,fmt,data);*/
      }
      else if(strstr(fmt,"d") != NULL){
	sprintf(item,fmt,*(int *)data);
	fwrite(item,strlen(item),1,fpout);
	/*fprintf(fpout,fmt,*(int *)data);*/
      }
      else if(strstr(fmt,"e") != NULL){
	sprintf(item,fmt,(double)(*(Real *)data));
	fwrite(item,strlen(item),1,fpout);
	/*fprintf(fpout,fmt,(double)(*(Real *)data));*/
      }
      else if(strstr(fmt,"f") != NULL){
	sprintf(item,fmt,(double)(*(Real *)data));
	fwrite(item,strlen(item),1,fpout);
	/*fprintf(fpout,fmt,(double)(*(Real *)data));*/
      }
      else if(strstr(fmt,"g") != NULL){
	sprintf(item,fmt,(double)(*(Real *)data));
	fwrite(item,strlen(item),1,fpout);
	/*fprintf(fpout,fmt,(double)(*(Real *)data));*/
      }
      else
	{
	  printf("write_ksprop_info_item: Unrecognized data type %s\n",fmt);
	  return 1;
	}
    }
  /*fprintf(fpout,"\n");*/
  fwrite("\n",1,1,fpout);
  return 0;

} /* end write_ksprop_info_item() */

/*----------------------------------------------------------------------*/
/* Open, write, and close the ASCII info file */

void write_ksprop_info_file(ks_prop_file *pf)
{
  FILE *info_fp;
  ks_prop_header *ph;
  char info_filename[256];
  char sums[20];

  ph = pf->header;

  /* Construct metadata file name from propagator file name 
   by adding filename extension to propagator file name */

  strcpy(info_filename,pf->filename);
  strcat(info_filename,ASCII_INFO_EXT);

  /* Open metadata file */
  
  if((info_fp = fopen(info_filename,"w")) == NULL)
    {
      printf("write_ksprop_info_file: Can't open ascii info file %s\n",
	     info_filename);
      return;
    }
  
  /* Write required information */

  write_ksprop_info_item(info_fp,"magic_number","%d",(char *)&ph->magic_number,0,0);
  write_ksprop_info_item(info_fp,"time_stamp","\"%s\"",ph->time_stamp,0,0);
  sprintf(sums,"%x %x",pf->check.sum29,pf->check.sum31);
  write_ksprop_info_item(info_fp,"checksums","\"%s\"",sums,0,0);
  write_ksprop_info_item(info_fp,"nx","%d",(char *)&nx,0,0);
  write_ksprop_info_item(info_fp,"ny","%d",(char *)&ny,0,0);
  write_ksprop_info_item(info_fp,"nz","%d",(char *)&nz,0,0);
  write_ksprop_info_item(info_fp,"nt","%d",(char *)&nt,0,0);

  write_appl_ksprop_info(info_fp);

  fclose(info_fp);

  printf("Wrote info file %s\n",info_filename); 
  fflush(stdout);

} /*write_ksprop_info_file */

/*----------------------------------------------------------------------*/

/* Set up the input prop file and prop header structures */

ks_prop_file *create_input_ksprop_file_handle(char *filename)
{
  ks_prop_file *pf;
  ks_prop_header *ph;
  char myname[] = "create_input_ksprop_file_handle";

  /* Allocate space for the file structure */

  pf = (ks_prop_file *)malloc(sizeof(ks_prop_file));
  if(pf == NULL)
    {
      printf("%s: Can't malloc pf\n", myname);
      terminate(1);
    }

  pf->filename = filename;
  pf->fp = NULL;

  /* Allocate space for the header */

  /* Make sure compilation gave us a 32 bit integer type */
  assert(sizeof(int32type) == 4);

  ph = (ks_prop_header *)malloc(sizeof(ks_prop_header));
  if(ph == NULL)
    {
      printf("%s: Can't malloc ph\n", myname);
      terminate(1);
    }

  pf->header       = ph;
  pf->prop         = NULL;
  pf->info         = NULL;
  pf->byterevflag  = 0;
  pf->parallel     = 0;
  pf->info_fp      = NULL;
  pf->file_type    = FILE_TYPE_UNKNOWN;
#ifdef HAVE_QIO
  pf->infile       = NULL;
  pf->outfile      = NULL;
#endif

  return pf;
}

/*----------------------------------------------------------------------*/

/* Set up the output prop file and prop header structure */

ks_prop_file *create_output_ksprop_file_handle(void)
{
  ks_prop_file *pf;
  ks_prop_header *ph;
  time_t time_stamp;
  int i;

  /* Allocate space for a new file structure */

  pf = (ks_prop_file *)malloc(sizeof(ks_prop_file));
  if(pf == NULL)
    {
      printf("setup_ksprop_header: Can't malloc pf\n");
      terminate(1);
    }

  /* Allocate space for a new header structure */

  ph = (ks_prop_header *)malloc(sizeof(ks_prop_header));
  if(ph == NULL)
    {
      printf("setup_ksprop_header: Can't malloc ph\n");
      terminate(1);
    }

  /* Load header pointer and file name */
  pf->header = ph;

  /* Initialize */
  pf->check.sum29 = 0;
  pf->check.sum31 = 0;
  pf->prop        = NULL;
  pf->info        = NULL;
  pf->file_type   = FILE_TYPE_UNKNOWN;

  /* Load header values */

  ph->magic_number = KSPROP_VERSION_NUMBER;

  ph->dims[0] = nx;
  ph->dims[1] = ny;
  ph->dims[2] = nz;
  ph->dims[3] = nt;

  /* Get date and time stamp. (We use local time on node 0) */

  if(this_node==0)
    {
      time(&time_stamp);
      strcpy(ph->time_stamp,ctime(&time_stamp));
      /* For aesthetic reasons, don't leave trailing junk bytes here to be
	 written to the file */
      for(i = strlen(ph->time_stamp) + 1; i < (int)sizeof(ph->time_stamp); i++)
	ph->time_stamp[i] = '\0';
      
      /* Remove trailing end-of-line character */
      if(ph->time_stamp[strlen(ph->time_stamp) - 1] == '\n')
	ph->time_stamp[strlen(ph->time_stamp) - 1] = '\0';
    }
  
  /* Broadcast to all nodes */
  broadcast_bytes(ph->time_stamp,sizeof(ph->time_stamp));

  return pf;

} /* setup_output_prop_file */


void destroy_ksprop_file_handle(ks_prop_file *kspf){
  if(kspf == NULL)return;

  if(kspf->header   != NULL)free(kspf->header);
  if(kspf->prop     != NULL)free(kspf->prop);
  if(kspf->info     != NULL)free(kspf->info);
  free(kspf);
}

/*----------------------------------------------------------------------*/

/* Open a binary file for serial writing by node 0 */

ks_prop_file *w_serial_ks_i(char *filename)
{
  /* Only node 0 opens the file filename */
  /* Returns a file structure describing the opened file */

  FILE *fp;
  ks_prop_file *kspf;
  ks_prop_header *ksph;

  /* Set up ks_prop file and ks_prop header structs and load header values */
  kspf = create_output_ksprop_file_handle();
  ksph = kspf->header;

  /* Indicate coordinate natural ordering */

  ksph->order = NATURAL_ORDER;

  /* Only node 0 opens the requested file */

  if(this_node == 0)
    {
      fp = fopen(filename, "wb");
      if(fp == NULL)
	{
	  printf("w_serial_ks_i: Node %d can't open file %s, error %d\n",
		 this_node,filename,errno);fflush(stdout);
	  terminate(1);
	}

/*      printf("Opened prop file %s for serial writing\n",filename); */
      
      /* Node 0 writes the header */
      
      swrite_ks_prop_hdr(fp,ksph);

    }
  
  /* Assign values to file structure */

  if(this_node==0) kspf->fp = fp; 
  else kspf->fp = NULL;                /* Only node 0 knows about this file */

  kspf->filename       = filename;
  kspf->byterevflag    = 0;            /* Not used for writing */
  kspf->parallel       = SERIAL;

  /* Node 0 writes ascii info file */

  if(this_node == 0) write_ksprop_info_file(kspf);

  return kspf;

} /* w_serial_ks_i */

/*---------------------------------------------------------------------------*/
/* Write checksum to lattice file.  It is assumed that the file
   is already correctly positioned.

   Should be called only by one node */

void write_checksum_ks(int parallel, ks_prop_file *kspf)
{

  char myname[] = "write_checksum";

  pswrite_data(parallel,kspf->fp,
	       &kspf->check.color,sizeof(kspf->check.color),myname,"checksum");
  pswrite_data(parallel,kspf->fp,
	       &kspf->check.sum29,sizeof(kspf->check.sum29),myname,"checksum");
  pswrite_data(parallel,kspf->fp,
	       &kspf->check.sum31,sizeof(kspf->check.sum31),myname,"checksum");
  /* printf("Prop checksums %x %x for color %d\n", 
     kspf->check.sum29,kspf->check.sum31,kspf->check.color);*/

} /* write_checksum_ks */

/*---------------------------------------------------------------------------*/
/* Here only node 0 writes propagator to a serial file 
   The propagator is passed as 3 su3_vectors, one for each source color 
*/

void w_serial_ks(ks_prop_file *kspf, int color, field_offset src_site,
		 su3_vector *src_field)
{
  /* kspf  = file descriptor as opened by ks_serial_w_i 
     src   = field offset for propagator su3 vector (type su3_vector)  */

  FILE *fp = NULL;
  ks_prop_header *ksph;
  u_int32type *val;
  int rank29,rank31;
  fsu3_vector *pbuf = NULL;
  int fseek_return;  /* added by S.G. for large file debugging */
  struct {
    fsu3_vector ksv;
    char pad[PAD_SEND_BUF]; /* Introduced because some switches
			       perform better if message lengths are longer */
  } msg;
  int buf_length;

  register int i,j,k;
  off_t offset;             /* File stream pointer */
  off_t ks_prop_size;       /* Size of propagator blocks for all nodes */
  off_t ks_prop_check_size; /* Size of propagator checksum record */
  off_t coord_list_size;    /* Size of coordinate list in bytes */
  off_t head_size = 0;      /* Size of header plus coordinate list */
  off_t body_size = 0;      /* Size of propagator blocks for all nodes 
			      plus checksum record */
  int currentnode,newnode;
  int x,y,z,t;

  if(this_node==0)
    {
      if(kspf->parallel == PARALLEL)
	printf("w_serial_ks: Attempting serial write to file opened in parallel \n");

      pbuf = (fsu3_vector *)malloc(MAX_BUF_LENGTH*sizeof(fsu3_vector));
      if(pbuf == NULL)
	{
	  printf("w_serial_ks: Node 0 can't malloc pbuf\n"); 
	  fflush(stdout); terminate(1);
        }

      fp = kspf->fp;
      ksph = kspf->header;

      ks_prop_size = volume*sizeof(fsu3_vector);
      ks_prop_check_size = sizeof(kspf->check.color) +
	sizeof(kspf->check.sum29) + sizeof(kspf->check.sum31);
      body_size = ks_prop_size + ks_prop_check_size;

      /* No coordinate list was written because fields are to be written
	 in standard coordinate list order */
      
      coord_list_size = 0;
      head_size = ksph->header_bytes + coord_list_size;
      
      offset = head_size + body_size*color + ks_prop_check_size;

      fseek_return=g_seek(fp,offset,SEEK_SET);
      /* printf("w_serial_ks: Node %d fseek_return = %d\n",this_node,fseek_return); */
      if( fseek_return < 0 ) 
	{
	  printf("w_serial_ks: Node %d g_seek %lld failed error %d file %s\n",
		 this_node, (long long)offset, errno, kspf->filename);
	  fflush(stdout); terminate(1);
	}

    } /* end if(this_node==0) */

  /* Buffered algorithm for writing fields in serial order */
  
  /* initialize checksums */
  kspf->check.sum31 = 0;
  kspf->check.sum29 = 0;
  /* counts 32-bit words mod 29 and mod 31 in order of appearance on file */
  /* Here only node 0 uses these values */
  rank29 = sizeof(fsu3_vector)/sizeof(int32type)*sites_on_node*this_node % 29;
  rank31 = sizeof(fsu3_vector)/sizeof(int32type)*sites_on_node*this_node % 31;

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
	      /* Copy or convert from structure to msg */
	      if(src_site == (field_offset)(-1))
		d2f_vec(src_field+3*i+color, &msg.ksv);
	      else
		d2f_vec((su3_vector *)F_PT(&lattice[i],src_site), &msg.ksv);
	    }
	  else
	    {
	      get_field((char *)&msg, sizeof(msg),currentnode);
	    }

	  pbuf[buf_length] = msg.ksv;

	  /* Accumulate checksums - contribution from next site */
	  for(k = 0, val = (u_int32type *)&pbuf[buf_length]; 
	      k < (int)sizeof(fsu3_vector)/(int)sizeof(int32type); k++, val++)
	    {
	      kspf->check.sum29 ^= (*val)<<rank29 | (*val)>>(32-rank29);
	      kspf->check.sum31 ^= (*val)<<rank31 | (*val)>>(32-rank31);
	      rank29++; if(rank29 >= 29)rank29 = 0;
	      rank31++; if(rank31 >= 31)rank31 = 0;
	    }

	  buf_length++;
	  
	  if( (buf_length == MAX_BUF_LENGTH) || (j == volume-1))
	    {
	      /* write out buffer */
	      
	      if( (int)fwrite(pbuf,sizeof(fsu3_vector),buf_length,fp) 
		  != buf_length)
		{
		  printf("w_serial_ks: Node %d prop write error %d file %s\n",
			 this_node,errno,kspf->filename); 
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
	    /* Copy data into send buffer and send to node 0 with padding */
	    if(src_site == (field_offset)(-1))
	      d2f_vec(src_field+3*i+color, &msg.ksv);
	    else
	      d2f_vec((su3_vector *)F_PT(&lattice[i],src_site), &msg.ksv);
	    send_field((char *)&msg, sizeof(msg),0);
	  }
	}
      
    } /*close x,y,z,t loops */
  
  g_sync();

  if(this_node==0)
    {
/*      printf("Wrote prop serially to file %s\n", kspf->filename); 
	fflush(stdout); */
      free(pbuf);

      /* Construct check record */

      /* Convert to 32 bit integer */
      kspf->check.color = color;

      /* Position file pointer for writing check record */

      offset = head_size + body_size*color;
      if( g_seek(fp,offset,SEEK_SET) < 0 ) 
	{
	  printf("w_serial_ks: Node %d g_seek %lld failed error %d file %s\n",
		 this_node,(long long)offset,errno,kspf->filename);
	  fflush(stdout); terminate(1);   
	}
      
      write_checksum_ks(SERIAL,kspf);
    }

} /* w_serial_ks */

/*---------------------------------------------------------------------------*/
/* Here only node 0 writes propagator to a serial file 
   The propagator is passed as a single su3_vector in the site structure
*/

void w_serial_ks_from_site(ks_prop_file *kspf, int color, 
			   field_offset src_site)
{
  w_serial_ks(kspf, color, src_site, NULL);
}
/*---------------------------------------------------------------------------*/
/* Here only node 0 writes propagator to a serial file 
   The propagator is passed as a single su3_vector in a field.
   We assume the field is organized with three contiguous su3_vectors
   per site.  But this call reads the vector for one color at a time.
   So the values for the ith site are taken from src_field[color + 3*i].
*/

void w_serial_ks_from_field(ks_prop_file *kspf, int color, 
			    su3_vector *src_field)
{
  w_serial_ks(kspf, color, (field_offset)(-1), src_field);
}
/*---------------------------------------------------------------------------*/

void w_serial_ks_f(ks_prop_file *kspf)

/* this subroutine closes the file and frees associated structures */
{
  g_sync();
  if(this_node==0)
    {
      if(kspf->parallel == PARALLEL)
	printf("w_serial_ks_f: Attempting serial close on file opened in parallel \n");

      printf("Wrote prop file %s time stamp %s\n", kspf->filename,
	     (kspf->header)->time_stamp);

      if(kspf->fp != NULL) fclose(kspf->fp);
    }

  /* Free header and file structures */
  destroy_ksprop_file_handle(kspf);

} /* w_serial_ks_f */

/*----------------------------------------------------------------------*/

int read_ks_prop_hdr(ks_prop_file *kspf, int parallel)
{
  /* parallel = 1 (TRUE) if all nodes are accessing the file */
  /*            0        for access from node 0 only */

  /* Returns byterevflag  = 0 or 1 */

  FILE *fp;
  ks_prop_header *ksph;
  int32type tmp;
  int j;
  int byterevflag = 0;
  char myname[] = "read_ks_prop_hdr";

  fp = kspf->fp;
  ksph = kspf->header;

  /* Read and verify magic number */

  if(psread_data(parallel, fp,&ksph->magic_number,sizeof(ksph->magic_number),
	      myname,"magic number")!=0)terminate(1);

  tmp = ksph->magic_number;
  
  if(ksph->magic_number == KSPROP_VERSION_NUMBER) 
    byterevflag=0;
  else 
    {
      byterevn((int32type *)&ksph->magic_number,1);
      if(ksph->magic_number == KSPROP_VERSION_NUMBER) 
	{
	  byterevflag=1;
	  /** printf("Reading with byte reversal\n"); **/
	  if( sizeof(Real) != sizeof(int32type)) {
	    printf("%s: Can't byte reverse\n",myname);
	    printf("requires size of int32type(%d) = size of Real(%d)\n",
		   (int)sizeof(int32type),(int)sizeof(Real));
	    terminate(1);
	  }
	}
      else
	{
	  /* Restore magic number as originally read */
	  ksph->magic_number = tmp;
	  
	  /* End of the road. */
	  printf("%s: Unrecognized magic number in prop file header.\n",
		 myname);
	  printf("Expected %x but read %x\n",
		 KSPROP_VERSION_NUMBER,tmp);
	  terminate(1);
	}
    }
  
  /* Read header, do byte reversal, 
     if necessary, and check consistency */
  
  /* Lattice dimensions */
  
  if(psread_byteorder(byterevflag,parallel,fp,ksph->dims,sizeof(ksph->dims),
		   myname,"dimensions")!=0) terminate(1);

  if(ksph->dims[0] != nx || 
     ksph->dims[1] != ny ||
     ksph->dims[2] != nz ||
     ksph->dims[3] != nt)
    {
      /* So we can use this routine to discover the dimensions,
	 we provide that if nx = ny = nz = nt = -1 initially
	 we don't die */
      if(nx != -1 || ny != -1 || nz != -1 || nt != -1)
	{
	  printf("%s: Incorrect lattice dimensions ",myname);
	  for(j=0;j<4;j++)
	    printf("%d ",ksph->dims[j]); 
	  printf("\n"); fflush(stdout); terminate(1);
	}
      else
	{
	  nx = ksph->dims[0];
	  ny = ksph->dims[1];
	  nz = ksph->dims[2];
	  nt = ksph->dims[3];
	  volume = nx*ny*nz*nt;
	}
    }
  
  /* Date and time stamp */

  if(psread_data(parallel,fp,ksph->time_stamp,sizeof(ksph->time_stamp),
	      myname,"time stamp")!=0) terminate(1);

  /* Header byte length */

  ksph->header_bytes = sizeof(ksph->magic_number) + sizeof(ksph->dims) + 
    sizeof(ksph->time_stamp) + sizeof(ksph->order);
  
  /* Data order */
  
  if(psread_byteorder(byterevflag,parallel,fp,&ksph->order,sizeof(ksph->order),
		   myname,"order parameter")!=0) terminate(1);
  
  return byterevflag;
  
} /* read_ks_prop_hdr */

/*---------------------------------------------------------------------------*/

ks_prop_file *r_serial_ks_i(char *filename)
{
  /* Returns file descriptor for opened file */

  ks_prop_header *ksph;
  ks_prop_file *kspf;
  FILE *fp;
  int byterevflag;

  /* All nodes set up a propagator file and propagator header
     structure for reading */

  kspf = create_input_ksprop_file_handle(filename);
  ksph = kspf->header;

  /* File opened for serial reading */
  kspf->parallel = 0;

  /* Node 0 alone opens a file and reads the header */

  if(this_node==0)
    {
      fp = fopen(filename, "rb");
      if(fp == NULL)
	{
	  printf("r_serial_ks_i: Node %d can't open file %s, error %d\n",
		 this_node,filename,errno);fflush(stdout);terminate(1);
	}
      
/*      printf("Opened prop file %s for serial reading\n",filename); */
      
      kspf->fp = fp;

      byterevflag = read_ks_prop_hdr(kspf,SERIAL);
    }

  /* Read the metadata file */

  kspf->info = read_info_file(filename);

  /* Node 0 broadcasts the byterevflag from node 0 to all nodes */
      
  broadcast_bytes((char *)&byterevflag,sizeof(byterevflag));
  kspf->byterevflag = byterevflag;
  
  /* Node 0 broadcasts the header structure to all nodes */
  
  broadcast_bytes((char *)ksph,sizeof(ks_prop_header));
  return kspf;

}/* r_serial_ks_i */

/*----------------------------------------------------------------------*/

/* Here only node 0 reads the KS propagator from a binary file */

int r_serial_ks(ks_prop_file *kspf, int color, field_offset dest_site,
		su3_vector *dest_field)
{
  /* 0 is normal exit code
     1 for seek, read error, or missing data error */

  FILE *fp;
  ks_prop_header *ksph;
  char *filename;
  int byterevflag;

  off_t offset ;            /* File stream pointer */
  off_t ks_prop_size;       /* Size of propagator blocks for all nodes */
  off_t ks_prop_check_size; /* Size of propagator checksum record */
  off_t coord_list_size;    /* Size of coordinate list in bytes */
  off_t head_size = 0;      /* Size of header plus coordinate list */
  off_t body_size = 0;      /* Size of propagator blocks for all nodes 
			      plus checksum record */
  int rcv_rank, rcv_coords;
  int destnode;
  int k,x,y,z,t;
  int status;
  int buf_length = 0, where_in_buf = 0;
  ks_prop_check test_kspc;
  u_int32type *val;
  int rank29,rank31;
  fsu3_vector *pbuf = NULL;
  su3_vector *dest;
  int idest = 0;

  struct {
    fsu3_vector ksv;
    char pad[PAD_SEND_BUF];    /* Introduced because some switches
				  perform better if message lengths are longer */
  } msg;

  char myname[] = "r_serial_ks";

  fp = kspf->fp;
  ksph = kspf->header;
  filename = kspf->filename;
  byterevflag = kspf->byterevflag;

  status = 0;
  if(this_node == 0)
    {
      if(kspf->parallel == PARALLEL)
	printf("%s: Attempting serial read from parallel file \n",myname);
      
      ks_prop_size = volume*sizeof(fsu3_vector) ;
      ks_prop_check_size = sizeof(kspf->check.color) +
	sizeof(kspf->check.sum29) + sizeof(kspf->check.sum31);

      body_size = ks_prop_size + ks_prop_check_size;
    }
  broadcast_bytes((char *)&status,sizeof(int));
  if(status != 0) return status;

  status = 0;
  if(this_node == 0)
    {
      if(ksph->order == NATURAL_ORDER) coord_list_size = 0;
      else coord_list_size = sizeof(int32type)*volume;
      head_size = ksph->header_bytes + coord_list_size;
     
      offset = head_size + body_size*color;

      pbuf = (fsu3_vector *)malloc(MAX_BUF_LENGTH*sizeof(fsu3_vector));
      if(pbuf == NULL)
	{
	  printf("%s: Node %d can't malloc pbuf\n",myname,this_node);
	  fflush(stdout);
	  terminate(1);
	}
      
      /* Position file pointer for reading check record */

      if( g_seek(fp,offset,SEEK_SET) < 0 ) 
	{
	  printf("%s: Node %d g_seek %lld failed error %d file %s\n",
		 myname,this_node,(long long)offset,errno,filename);
	  fflush(stdout);
	  status = 1;
	}
      
      /* Read check record */

      status += sread_byteorder(byterevflag,fp,&kspf->check.color,
		      sizeof(kspf->check.color),myname,"check.color");
      status += sread_byteorder(byterevflag,fp,&kspf->check.sum29,
		      sizeof(kspf->check.sum29),myname,"check.sum29");
      status += sread_byteorder(byterevflag,fp,&kspf->check.sum31,
		      sizeof(kspf->check.sum31),myname,"check.sum31");

      /* Verify spin and color - checksums come later */
      if(kspf->check.color != color)
	{
	  printf("%s: color %d does not match check record on file %s\n",
		 myname,color,filename);
	  printf("  Check record said %d\n",kspf->check.color);
	  fflush(stdout); 
	  status = 1;
	}

      buf_length = 0;
      where_in_buf = 0;

    }
  
  broadcast_bytes((char *)&status,sizeof(int));
  if(status != 0) return status;

  /* all nodes initialize checksums */
  test_kspc.sum31 = 0;
  test_kspc.sum29 = 0;
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
      /* If file is always in coordinate natural order */
      
      rcv_coords = rcv_rank;

      x = rcv_coords % nx;   rcv_coords /= nx;
      y = rcv_coords % ny;   rcv_coords /= ny;
      z = rcv_coords % nz;   rcv_coords /= nz;
      t = rcv_coords % nt;
      
      /* The node that gets the next su3_vector */
      destnode=node_number(x,y,z,t);

      if(this_node==0){
	/* Node 0 fills its buffer, if necessary */
	if(where_in_buf == buf_length)
	  {  /* get new buffer */
	    /* new buffer length  = remaining sites, but never bigger 
	       than MAX_BUF_LENGTH */
	    buf_length = volume - rcv_rank;
	    if(buf_length > MAX_BUF_LENGTH) buf_length = MAX_BUF_LENGTH;
	    /* then do read */
	    
	    if( (int)fread(pbuf,sizeof(fsu3_vector),buf_length,fp) 
		!= buf_length)
	      {
		if(status == 0)
		  printf("%s: node %d propagator read error %d file %s\n",
			 myname, this_node, errno, filename); 
		fflush(stdout); 
		status = 1;
	      }
	    where_in_buf = 0;  /* reset counter */
	  }  /*** end of the buffer read ****/

	/* Save vector in msg structure for further processing */
	msg.ksv = pbuf[where_in_buf];
	if(destnode==0){	/* just copy su3_vector */
	  idest = node_index(x,y,z,t);
	}
	else {		        /* send to correct node */
	  send_field((char *)&msg, sizeof(msg), destnode);
	}
	where_in_buf++;
      }

      /* The node that contains this site reads the message */
      else {	/* for all nodes other than node 0 */
	if(this_node==destnode){
	  idest = node_index(x,y,z,t);
	  /* Receive padded message in msg */
	  get_field((char *)&msg, sizeof(msg),0);
	}
      }

      /* The receiving node does the byte reversal and then checksum,
         if needed.  At this point msg.ksv contains the input vector
         and idest points to the destination site structure. */
      
      if(this_node==destnode)
	{
	  if(byterevflag==1)
	    byterevn((int32type *)(&msg.ksv),
		     sizeof(fsu3_vector)/sizeof(int32type));
	  /* Accumulate checksums */
	  for(k = 0, val = (u_int32type *)(&msg.ksv); 
	      k < (int)sizeof(fsu3_vector)/(int)sizeof(int32type); k++, val++)
	    {
	      test_kspc.sum29 ^= (*val)<<rank29 | (*val)>>(32-rank29);
	      test_kspc.sum31 ^= (*val)<<rank31 | (*val)>>(32-rank31);
	      rank29++; if(rank29 >= 29)rank29 = 0;
	      rank31++; if(rank31 >= 31)rank31 = 0;
	    }
	  /* Copy or convert vector from msg to lattice[idest] */
	  if(dest_site == (field_offset)(-1))
	    dest = dest_field + 3*idest + color;
	  else
	    dest = (su3_vector *)F_PT( &(lattice[idest]), dest_site );
	  f2d_vec(&msg.ksv, dest);
	}
      else
	{
	  rank29 += sizeof(fsu3_vector)/sizeof(int32type);
	  rank31 += sizeof(fsu3_vector)/sizeof(int32type);
	  rank29 %= 29;
	  rank31 %= 31;
	}
    }

  broadcast_bytes((char *)&status,sizeof(int));
  if(status != 0) return status;

  /* Combine node checksum contributions with global exclusive or */
  g_xor32(&test_kspc.sum29);
  g_xor32(&test_kspc.sum31);
  
  if(this_node==0)
    {
      printf("Read prop serially for color %d from file %s\n", 
	     color, filename);
      
      /* Verify checksum */
      /* Checksums not implemented until version 5 */

      if(ksph->magic_number == KSPROP_VERSION_NUMBER)
	{
	  if(kspf->check.sum29 != test_kspc.sum29 ||
	     kspf->check.sum31 != test_kspc.sum31)
	    {
	      printf("%s: Checksum violation color %d file %s\n",
		     myname, kspf->check.color, kspf->filename);
	      printf("Computed checksum %x %x.  Read %x %x.\n",
		     test_kspc.sum29, test_kspc.sum31,
		     kspf->check.sum29, kspf->check.sum31);
	    }
/*	  else
	    printf("Checksums %x %x OK for file %s\n",
		   kspf->check.sum29,kspf->check.sum31,kspf->filename); */
	}
      fflush(stdout);
      free(pbuf);
    }

  return 0;

} /* r_serial_ks */

/*----------------------------------------------------------------------*/

/* Node 0 reads the KS propagator from a binary file with result in the
   site structure*/

int r_serial_ks_to_site(ks_prop_file *kspf, int color, field_offset dest_site)
{
  return r_serial_ks(kspf, color, dest_site, NULL);
}
/*----------------------------------------------------------------------*/

/* Node 0 reads the KS propagator from a binary file with result in 
   a field.
   The propagator is read as a single su3_vector in a field.
   We assume the field is organized with three contiguous su3_vectors
   per site.  But this call reads the vector for one color at a time.
   So the values for the ith site are read to dest_field[color + 3*i]
 */

int r_serial_ks_to_field(ks_prop_file *kspf, int color, su3_vector *dest_field)
{
  return r_serial_ks(kspf, color, (field_offset)(-1), dest_field);
}
/*----------------------------------------------------------------------*/

void r_serial_ks_f(ks_prop_file *kspf)

/* Close the file and free associated structures */
{
  g_sync();

  if(kspf == NULL)return;

  if(this_node==0)
    {
      if(kspf->parallel == PARALLEL)
	printf("r_serial_w_f: Attempting serial close on parallel file \n");
      
      if(kspf->fp != NULL) fclose(kspf->fp);
/*      printf("Closed prop file %s\n",kspf->filename);*/
      fflush(stdout);
    }
  
#ifdef HAVE_QIO
  if(kspf->infile   != NULL){
    QIO_close_read(kspf->infile);
    kspf->infile = NULL;
  }
#endif
  destroy_ksprop_file_handle(kspf);

} /* r_serial_ks_f */

/*---------------------------------------------------------------------------*/

/* ASCII file format

   format:
    version_number (int)
    time_stamp (char string enclosed in quotes)
    nx ny nz nt (int)
    for(t=...)for(z=...)for(y=...)for(x=...){
       for(i=...)for(j=...){prop[i][j].real, prop[i][j].imag}
    }
*/

/*---------------------------------------------------------------------------*/
/* Open and write header info for ascii propagator file */

ks_prop_file *w_ascii_ks_i(char *filename)
{
  ks_prop_header *ksph;
  ks_prop_file *kspf;
  FILE *fp;

  kspf = create_output_ksprop_file_handle();
  ksph = kspf->header;

  /* node 0 does all the writing */
  if(this_node==0){

    /* Set up ksprop file and ksprop header structures & load header values */

    fp = fopen(filename,"w");
    if(fp==NULL){
      printf("Can't open file %s, error %d\n",filename,errno); terminate(1);
    }

    kspf->fp = fp;

    if( (fprintf(fp,"%d\n", KSPROP_VERSION_NUMBER))==0 ){
      printf("Error in writing version number\n"); terminate(1);
    }
    if( (fprintf(fp,"\"%s\"\n",ksph->time_stamp))==0 ){
      printf("Error in writing time stamp\n"); terminate(1);
    }
    
    if( (fprintf(fp,"%d\t%d\t%d\t%d\n",nx,ny,nz,nt))==0 ){
      printf("Error in writing dimensions\n"); terminate(1);
    }

  }
  else kspf->fp = NULL;

  /* Assign remaining values to propagator file structure */
  kspf->parallel = 0;
  kspf->filename       = filename;
  kspf->byterevflag    = 0;            /* Not used for writing */

  /* Node 0 writes info file */
  if(this_node==0) write_ksprop_info_file(kspf);

  return kspf;

} /* w_ascii_ks_i */
  
/*---------------------------------------------------------------------------*/
/* Write ASCII propagator from field */

void w_ascii_ks(ks_prop_file *kspf, int color, su3_vector *src_field)
{
  FILE *fp;
  int currentnode,newnode;
  int b,l,x,y,z,t;
  fsu3_vector pbuf;
  int node0=0;

  g_sync();
  currentnode=0;

  fp = kspf->fp;
  
  /* first the color and spin */
  if(this_node==0)
    {
      if( (fprintf(fp,"%d\n",color))== EOF)
	{
	  printf("w_ascii_ks: error writing color\n"); 
	  terminate(1);
	} 
      fflush(fp);
    }

  /* next the elements */
  for(t=0;t<nt;t++)for(z=0;z<nz;z++)for(y=0;y<ny;y++)for(x=0;x<nx;x++)
    {
      newnode = node_number(x,y,z,t);
      if(newnode != currentnode)
	{	/* switch to another node */
	  g_sync();
	  /* Send a few bytes of garbage to tell newnode it's OK to send */
	  if( this_node==0 && newnode!=0 )send_field((char *)&pbuf,1,newnode);
	  if( this_node==newnode && newnode!=0 )get_field((char *)&pbuf,1,0);
	  currentnode=newnode;
	}
      
      if(this_node==0)
	{
	  if(currentnode==0)
	    {
	      l=node_index(x,y,z,t);
	      /* Copy, converting precision if necessary */
	      d2f_vec(src_field+3*l+color, &pbuf);
	      //  d2f_vec((su3_vector *)F_PT( &(lattice[l]), src ), &pbuf);
	    }
	  else
	    {
	      get_field((char *)&pbuf,sizeof(fsu3_vector),currentnode);
	    }
	  for(b=0;b<3;b++)
	    {
	      if( (fprintf(fp,"%.7e\t%.7e\n",(double)pbuf.c[b].real,
			   (double)pbuf.c[b].imag)) == EOF)
		{
		  printf("w_ascii_ks: error writing prop\n"); 
		  terminate(1);
		} 
	    }
	}
      else
	{	/* for nodes other than 0 */
	  if(this_node==currentnode)
	    {
	      l=node_index(x,y,z,t);
	      /* Copy, converting precision if necessary */
	      d2f_vec(src_field+3*l+color, &pbuf);
	      // d2f_vec((su3_vector *)F_PT( &(lattice[l]), src ), &pbuf);
	      send_field((char *)&pbuf,sizeof(fsu3_vector),node0);
	    }
	}
    }
  g_sync();
  if(this_node==0)
    {
      fflush(fp);
      printf("Wrote prop to ASCII file  %s\n", kspf->filename);
    }

} /* w_ascii_ks */

/*---------------------------------------------------------------------------*/
/* Close ASCII propagator file */

void w_ascii_ks_f(ks_prop_file *kspf)
{
  g_sync();
  if(this_node==0)
    {
      fflush(kspf->fp);
      if(kspf->fp != NULL)fclose(kspf->fp);
      printf("Wrote ksprop file %s time stamp %s\n",kspf->filename,
	     (kspf->header)->time_stamp);
    }

  /* Free header and file structures */
  destroy_ksprop_file_handle(kspf);

}

/*---------------------------------------------------------------------------*/
/* Open ASCII propagator file and read header information */

ks_prop_file *r_ascii_ks_i(char *filename)
{
  ks_prop_file *kspf;
  ks_prop_header *ksph;
  FILE *fp;

  /* All nodes set up a propagator file and propagator header
     structure for reading */

  kspf = create_input_ksprop_file_handle(filename);
  ksph = kspf->header;

  /* File opened for serial reading */
  kspf->parallel = 0;
  kspf->byterevflag = 0;  /* Unused for ASCII */

  /* Indicate coordinate natural ordering */
  ksph->order = NATURAL_ORDER;

  /* Node 0 alone opens a file and reads the header */

  if(this_node==0)
    {
      fp = fopen(filename,"r");
      if(fp==NULL)
	{
	  printf("r_ascii_ks_i: Node %d can't open file %s, error %d\n",
		 this_node,filename,errno); fflush(stdout); terminate(1);
        }
      kspf->fp = fp;

      if( (fscanf(fp,"%d",&ksph->magic_number))!=1 )
	{
	  printf("r_ascii_ks_i: Error in reading version number\n"); 
	  terminate(1);
	}
      if(ksph->magic_number != KSPROP_VERSION_NUMBER)
	{
	  printf("r_ascii_ks_i: Unrecognized magic number in propagator file header.\n");
	  printf("Expected %d but read %d\n",
		     KSPROP_VERSION_NUMBER, ksph->magic_number);
	  terminate(1);
	}
      if(fscanf(fp,"%*[ \f\n\r\t\v]%*[\"]%[^\"]%*[\"]",
		     ksph->time_stamp)!=1)
	{
	  printf("r_ascii_ks_i: Error reading time stamp\n"); 
	  terminate(1);
	}
      if( (fscanf(fp,"%d%d%d%d",&ksph->dims[0],&ksph->dims[1],
		  &ksph->dims[2],&ksph->dims[3]))!=4 )
	{
	  printf("r_ascii_ks_i: Error reading lattice dimensions\n"); 
	  terminate(1);
	}
      if( ksph->dims[0]!=nx || ksph->dims[1]!=ny 
	 || ksph->dims[2]!=nz || ksph->dims[3]!=nt )
	{
	  /* So we can use this routine to discover the dimensions,
	     we provide that if nx = ny = nz = nt = -1 initially
	     we don't die */
	  if(nx != -1 || ny != -1 || nz != -1 || nt != -1)
	    {
	      printf("r_ascii_ks_i: Incorrect lattice size %d,%d,%d,%d\n",
		     ksph->dims[0],ksph->dims[1],ksph->dims[2],ksph->dims[3]);
	      terminate(1);
	    }
	  else
	    {
	      nx = ksph->dims[0];
	      ny = ksph->dims[1];
	      nz = ksph->dims[2];
	      nt = ksph->dims[3];
	      volume = nx*ny*nz*nt;
	    }
	}
      ksph->header_bytes = 0;    /* Unused for ASCII */
    }

  else kspf->fp = NULL;  /* Other nodes don't know about this file */

  /* Broadcasts the header structure from node 0 to all nodes */
  
  broadcast_bytes((char *)ksph, sizeof(ks_prop_header));

  return kspf;

} /* r_ascii_ks_i */

/*---------------------------------------------------------------------------*/
/* Read a propagator */

int r_ascii_ks(ks_prop_file *kspf, int color, su3_vector *dest_field)
{
  /* 0 normal exit code
     1 read error */

  FILE *fp;
  int destnode;
  int i,j,x,y,z,t;
  fsu3_vector pbuf;
  int status;

  fp = kspf->fp;

  g_sync();

  status = 0;
  if(this_node == 0)
    {
      if( (fscanf(fp,"%d",&j)) != 1 )
	{
	  printf("r_ascii_ks: Error reading color\n");
	  printf("r_ascii_ks: color=%d, j=%d\n",color,j);
	  status = 1;
	}
      if(status == 0 && (j != color))
	{
	  printf("r_ascii_ks: file file color=%d, prog color=%d\n", j, color);
	  status = 1;
	}
    }
  broadcast_bytes((char *)&status,sizeof(int));
  if(status != 0) return status;

  status = 0;
  for(t=0;t<nt;t++)for(z=0;z<nz;z++)for(y=0;y<ny;y++)for(x=0;x<nx;x++)
    {
      destnode=node_number(x,y,z,t);

      /* Node 0 reads, and sends site to correct node */
      if(this_node==0)
	{
	  for(j=0;j<3;j++)
	    {
	      if( (fscanf(fp,"%e%e\n",&(pbuf.c[j].real),
			  &(pbuf.c[j].imag)) )!= 2)
		{
		  if(status == 0)
		    printf("r_ascii_ks: Error reading su3_vector\n"); 
		  status = 1;
		}
	    }
	  if(destnode==0)
	    {              /* just copy su3_vector */
	      i = node_index(x,y,z,t);
	      /* Copy, converting precision if necessary */
	      f2d_vec(&pbuf, dest_field + 3*i + color);
	      // f2d_vec(&pbuf, (su3_vector *)F_PT( &(lattice[i]), src ));
	    }
	  else 
	    {              /* send to correct node */
	      send_field((char *)&pbuf, sizeof(fsu3_vector), destnode);
	    }
	}
      
      /* The node which contains this site reads message */
      else
	{ 
	  /* for all nodes other than node 0 */
	  if(this_node==destnode)
	    {
	      get_field((char *)&pbuf, sizeof(fsu3_vector),0);
	      i = node_index(x,y,z,t);
	      /* Copy, converting precision if necessary */
	      f2d_vec(&pbuf, dest_field + 3*i + color);
	      // f2d_vec(&pbuf, (su3_vector *)F_PT( &(lattice[i]), src ));
	    }
	}
    }

  broadcast_bytes((char *)&status, sizeof(int));

  return status;

} /* r_ascii_ks */

/*---------------------------------------------------------------------------*/
/* Close propagator file */

void r_ascii_ks_f(ks_prop_file *kspf)
{
  FILE *fp;

  fp = kspf->fp;

  g_sync();
  if(this_node==0)
    {
/*      printf("Closed ASCII prop file  %s\n", kspf->filename);*/
      fclose(fp);
      kspf->fp = NULL;
      fflush(stdout);
    }

} /* r_ascii_ks_f */

/*---------------------------------------------------------------------------*/
/* Quick and dirty code for binary output of propagator, separated
   into one file per timeslice */

void w_serial_ksprop_tt( char *filename, field_offset prop)
{

  char myname[] = "w_serial_ksprop_tt";
  FILE *fp = NULL;
  ks_prop_file *kspf;
  ks_prop_header *ksph;

  char tfilename[256];
  char *tag=".t";

  off_t offset;             /* File stream pointer */
  off_t ks_prop_size;        /* Size of propagator blocks for all nodes */
  off_t coord_list_size;    /* Size of coordinate list in bytes */
  int fseek_return;  /* added by S.G. for large file debugging */
  int currentnode, newnode;
  int x,y,z,t;
  register int i,a,b;
  fsu3_matrix pbuf;
  su3_vector *proppt;
  site *s;
    
  /* Set up ks_prop file and ks_prop header structs and load header values */
  kspf = create_output_ksprop_file_handle();
  ksph = kspf->header;
  ksph->order = NATURAL_ORDER;

  ks_prop_size = volume*sizeof(fsu3_matrix)/nt;

  /* No coordinate list was written because fields are to be written
     in standard coordinate list order */
  
  coord_list_size = 0;
  
  /* OLD, forget the header now:  offset = head_size; */
  offset = 0;

  for(t=0; t<nt; t++) {
    
    /* Only node 0 opens the requested file */
    if(this_node == 0) {

      sprintf(tfilename, "%s%s%d", filename, tag, t);
      fp = fopen(tfilename, "wb");
      if(fp == NULL)
	{
	  printf("%s: Node %d can't open file %s, error %d\n",
		 myname,this_node,filename,errno); fflush(stdout);
	  terminate(1);
	}

      /* Node 0 writes the header */
      /* forget the header swrite_ks_prop_hdr(fp,ksph); */

    }

    /* Assign values to file structure */
    
    if(this_node==0) kspf->fp = fp; 
    else kspf->fp = NULL;             /* Only node 0 knows about this file */
    
    kspf->filename       = filename;
    kspf->byterevflag    = 0;            /* Not used for writing */
    kspf->parallel       = SERIAL;
    

    if(this_node==0) {

      fseek_return=g_seek(fp,offset,SEEK_SET);
      /* printf("w_serial_ksprop_t: Node %d fseek_return = %d\n",this_node,fseek_return); */
      if( fseek_return < 0 ) 
	{
	  printf("%s: Node %d g_seek %lld failed error %d file %s\n",
		 myname,this_node, (long long)offset, errno, kspf->filename);
	  fflush(stdout); terminate(1);
	}

    } /* end if(this_node==0) */

    /* Synchronize */  
    g_sync();
    
    /* Write propagator */
    currentnode = 0;

    for(z=0;z<nz;z++)for(y=0;y<ny;y++)for(x=0;x<nx;x++){
      newnode = node_number(x,y,z,t);
      if(newnode != currentnode){	/* switch to another node */
	/* Send a few bytes of garbage to tell newnode it's OK to send */
	if( this_node==0 && newnode!=0 )send_field((char *)&pbuf,1,newnode);
	if( this_node==newnode && newnode!=0 )get_field((char *)&pbuf,1,0);
	currentnode=newnode;
      }

      if(this_node==0){
	
	if(currentnode==0){ 
	  /* just copy */
	  i = node_index(x,y,z,t);
	  s = &(lattice[i]);
	  proppt = (su3_vector *)F_PT(s,prop);
	  /* Copy, converting precision if necessary */
	  for(a=0; a<3; a++){
	    for(b=0; b<3; b++){
	      pbuf.e[a][b].real = proppt[a].c[b].real;
	      pbuf.e[a][b].imag = proppt[a].c[b].imag;
	    }
	  }
	}
	else {
	  get_field((char *)&pbuf,sizeof(fsu3_matrix),currentnode);
	}
	
	/* write out */
	for(a=0; a<3; a++){
	  for(b=0; b<3; b++){
	    if( (int)fwrite(&(pbuf.e[a][b].real), sizeof(Real), 1, fp)
		!= 1 ) {
	      printf("w_serial_ksprop_tt: Node %d prop write error %d file %s\n",
		     this_node,errno,kspf->filename); 
	      fflush(stdout);
	      terminate(1);   
	    }
	    if( (int)fwrite(&(pbuf.e[a][b].imag), sizeof(Real), 1, fp)
		!= 1 ) {
	      printf("w_serial_ksprop_tt: Node %d prop write error %d file %s\n",
		     this_node,errno,kspf->filename); 
	      fflush(stdout);
	      terminate(1);   
	    }
	  }
	}

      }  else {    /* for nodes other than 0 */

	  if(this_node==currentnode){
	    i=node_index(x,y,z,t);
	    /* Copy data into send buffer and send to node 0 with padding */
	    s = &(lattice[i]);
	    proppt = (su3_vector *)F_PT(s,prop);
	    /* Copy, converting precision if necessary */
	    for(a=0; a<3; a++){
	      for(b=0; b<3; b++){
		pbuf.e[a][b].real = proppt[a].c[b].real;
		pbuf.e[a][b].imag = proppt[a].c[b].imag;
	      }
	    }
	    send_field((char *)&pbuf,sizeof(fsu3_matrix),0);
	  }
      }

    } /*close x,y,z loops */
  
    g_sync();
    if(this_node==0) fclose(fp);

  } /* end loop over t */

} /* end write_serial_ksprop_t */

/*---------------------------------------------------------------------------*/
/* Quick and dirty code for ascii output of propagator, separated
   into one file per timeslice */

void w_ascii_ksprop_tt( char *filename, field_offset prop) 
{

  char myname[] = "w_ascii_ksprop_tt";
  FILE *fp = NULL;
  ks_prop_file *kspf;
  ks_prop_header *ksph;

  char tfilename[256];
  char *tag=".t";
  int currentnode, newnode;
  int x,y,z,t;
  register int i,a,b;
  fsu3_matrix pbuf;
  su3_vector *proppt;
  site *s;

  /* Set up ks_prop file and ks_prop header structs and load header values */
  kspf = create_output_ksprop_file_handle();
  ksph = kspf->header;
  ksph->order = NATURAL_ORDER;

  for(t=0; t<nt; t++) {

    if(this_node == 0) {
      sprintf(tfilename, "%s%s%d", filename, tag, t);
      fp = fopen(tfilename, "w");
      if(fp == NULL){
	printf("%s: Node %d can't open file %s, error %d\n",
	       myname,this_node,filename,errno); fflush(stdout);
	terminate(1);
      }
    }

    /* Assign values to file structure */

    if(this_node==0) kspf->fp = fp; 
    else kspf->fp = NULL;             /* Only node 0 knows about this file */

    kspf->filename       = filename;
    kspf->byterevflag    = 0;            /* Not used for writing */
    kspf->parallel       = SERIAL;

    /* Synchronize */  
    g_sync();
    
    /* Write propagator */
    currentnode = 0;

    for(z=0;z<nz;z++)for(y=0;y<ny;y++)for(x=0;x<nx;x++){
      newnode = node_number(x,y,z,t);
      if(newnode != currentnode){	/* switch to another node */
	/**g_sync();**/
	/* tell newnode it's OK to send */
	if( this_node==0 && newnode!=0 )send_field((char *)&pbuf,1,newnode);
	if( this_node==newnode && newnode!=0 )get_field((char *)&pbuf,1,0);
	currentnode=newnode;
      }

      if(this_node==0){

	if(currentnode==0){ 
	  i = node_index(x,y,z,t);
	  s = &(lattice[i]);
	  proppt = (su3_vector *)F_PT(s,prop);

	  /* Copy, converting precision if necessary */
	  for(a=0; a<3; a++){
	    for(b=0; b<3; b++){
	      pbuf.e[a][b].real = proppt[a].c[b].real;
	      pbuf.e[a][b].imag = proppt[a].c[b].imag;
	    }
	  }
	}
	else {
	  get_field((char *)&pbuf,sizeof(fsu3_matrix),currentnode);
	}
	
	for(a=0;a<3;a++)for(b=0;b<3;b++){
	  if( (fprintf(fp,"%.7e\t%.7e\n",(double)pbuf.e[a][b].real,
		       (double)pbuf.e[a][b].imag))== EOF){
	    printf("Write error in save_ksprop_ascii\n"); terminate(1);
	  }
	}

      }

      else {	/* for all nodes other than node 0 */
	if(this_node==currentnode){
	  i = node_index(x,y,z,t);
	  s = &(lattice[i]);
	  proppt = (su3_vector *)F_PT(s,prop);
	  /* Copy, converting precision if necessary */
	  for(a=0; a<3; a++){
	    for(b=0; b<3; b++){
	      pbuf.e[a][b].real = proppt[a].c[b].real;
	      pbuf.e[a][b].imag = proppt[a].c[b].imag;
	    }
	  }
	  send_field((char *)&pbuf,sizeof(fsu3_matrix),0);
	} /* if */
      } /* else */

    } /* end loop over x,y,z */
  
    g_sync();
    if(this_node==0){
      fflush(fp);
      fclose(fp);
      fflush(stdout);
    }

  } /* end loop over t */

} /* end w_ascii_ksprop_tt */

/*********************************************************************/
/*********************************************************************/
/*  The following: restore_ksprop_ascii and save_ksprop_ascii are
    deprecated in favor of [w/r]_ascii_ks_i, [w/r]_ascii_ks, and
    [w/r]_ascii_ks_f
*/

/* Read a KS propagator in ASCII format serially (node 0 only) */

/* format:
    version_number (int)
    time_stamp (char string enclosed in quotes)
    nx ny nz nt (int)
    for(t=...)for(z=...)for(y=...)for(x=...){
       for(i=...)for(j=...){prop[i][j].real, prop[i][j].imag}
    }
*/

/* one su3_vector for each source color */
ks_prop_file *restore_ksprop_ascii( char *filename, field_offset prop )
{

  ks_prop_header *ph;
  ks_prop_file *pf;
  FILE *fp = NULL;
  int destnode;
  int version_number,i,a,b,x,y,z,t;
  fsu3_matrix pbuf;
  int src_clr;
  su3_vector *proppt;
  site *s;

  /* Set up a prop file and prop header structure for reading */

  pf = create_input_ksprop_file_handle(filename);
  ph = pf->header;

  /* File opened for serial reading */
  pf->parallel = 0;

  /* Node 0 opens the file and reads the header */

  if(this_node==0){
    fp = fopen(filename,"r");
    if(fp==NULL){
      printf("Can't open file %s, error %d\n",filename,errno);
      terminate(1);
    }

    pf->fp = fp;

    if( (fscanf(fp,"%d",&version_number))!=1 ){
      printf("restore_ksprop_ascii: Error reading version number\n"); 
      terminate(1);
    }
    ph->magic_number = version_number;
    if(ph->magic_number != KSPROP_VERSION_NUMBER_V0){
      printf("restore_ksprop_ascii: Incorrect version number in lattice header\n");
      printf("  read %d but expected %d\n",
	     ph->magic_number, KSPROP_VERSION_NUMBER_V0);
      terminate(1);
    }
    /* Time stamp is enclosed in quotes - discard the leading white
       space and the quotes and read the enclosed string */
    if((i = fscanf(fp,"%*[ \f\n\r\t\v]%*[\"]%[^\"]%*[\"]",ph->time_stamp))!=1){
      printf("restore_ksprop_ascii: Error reading time stamp\n"); 
      printf("count %d time_stamp %s\n",i,ph->time_stamp);
      terminate(1);
    }
    if( (fscanf(fp,"%d%d%d%d",&x,&y,&z,&t))!=4 ){
      printf("restore_ksprop_ascii: Error in reading dimensions\n"); 
      terminate(1);
    }

    ph->dims[0] = x; ph->dims[1] = y; ph->dims[2] = z; ph->dims[3] = t;
    if( ph->dims[0]!=nx || ph->dims[1]!=ny || 
       ph->dims[2]!=nz || ph->dims[3]!=nt )
      {
	/* So we can use this routine to discover the dimensions,
	   we provide that if nx = ny = nz = nt = -1 initially
	   we don't die */
	if(nx != -1 || ny != -1 || nz != -1 || nt != -1)
	  {
	    printf("restore_ksprop_ascii: Incorrect lattice size %d,%d,%d,%d\n",
		   ph->dims[0],ph->dims[1],ph->dims[2],ph->dims[3]);
	    terminate(1);
	  }
	else
	  {
	    nx = ph->dims[0];
	    ny = ph->dims[1];
	    nz = ph->dims[2];
	    nt = ph->dims[3];
	    volume = nx*ny*nz*nt;
	  }
      }

  } /* if node 0 */

  else pf->fp = NULL;

  /* Node 0 broadcasts the header structure to all nodes */
  
  broadcast_bytes((char *)ph,sizeof(ks_prop_header));

  /* Synchronize */  
  g_sync();
  
  for(t=0;t<nt;t++)for(z=0;z<nz;z++)for(y=0;y<ny;y++)for(x=0;x<nx;x++){
    destnode=node_number(x,y,z,t);
    
    /* Node 0 reads, and sends site to correct node */
    if(this_node==0){
      for(a=0;a<3;a++)for(b=0;b<3;b++){
	if( (fscanf(fp,"%e%e\n",&(pbuf.e[a][b].real),
		      &(pbuf.e[a][b].imag)) )!= 2){
	    printf("restore_ksprop_ascii: propagator read error\n"); 
	    terminate(1);
	}
      }

      if(destnode==0){	/* just copy */
	i = node_index(x,y,z,t);
	s = &(lattice[i]);
	proppt = (su3_vector *)F_PT(s,prop);
	/* Copy, converting precision if necessary */
	for(src_clr=0; src_clr<3; src_clr++){
	  for(b=0; b<3; b++){
	    proppt[src_clr].c[b].real = pbuf.e[src_clr][b].real;
	    proppt[src_clr].c[b].imag = pbuf.e[src_clr][b].imag;
	  }
	}
      }
      else {		/* send to correct node */
	send_field((char *)&pbuf,sizeof(fsu3_matrix),destnode);
      }
    }

    /* The node which contains this site reads message */
    else {	/* for all nodes other than node 0 */
      if(this_node==destnode){
	get_field((char *)&pbuf,sizeof(fsu3_matrix),0);
	i = node_index(x,y,z,t);
	s = &(lattice[i]);
	proppt = (su3_vector *)F_PT(s,prop);
	/* Copy, converting precision if necessary */
	for(src_clr=0; src_clr<3; src_clr++){
	  for(b=0; b<3; b++){
	    proppt[src_clr].c[b].real = pbuf.e[src_clr][b].real;
	    proppt[src_clr].c[b].imag = pbuf.e[src_clr][b].imag;
	  }
	}
      } /* if */
    } /* else */

  } /* end loop over x,y,z,t */
  
  g_sync();
  if(this_node==0){
    printf("Restored propagator from ascii file  %s\n",
	   filename);
    printf("Time stamp %s\n",ph->time_stamp);
    fclose(fp);
    pf->fp = NULL;
    fflush(stdout);
  }

  return pf;

} /* end restore_ksprop_ascii() */

/*---------------------------------------------------------------------*/

/* Save a KS propagator in ASCII format serially (node 0 only) */

/* one su3_vector for each source color */
ks_prop_file *save_ksprop_ascii(char *filename, field_offset prop)
{

  ks_prop_header *ph;
  ks_prop_file *pf;
  FILE *fp = NULL;
  int currentnode, newnode;
  int i,a,b,x,y,z,t;
  fsu3_matrix pbuf;
  int src_clr;
  su3_vector *proppt;
  site *s;

  /* Set up a prop file and prop header structure for reading */

  pf = create_output_ksprop_file_handle();
  ph = pf->header;

  /* node 0 does all the writing */

  if(this_node==0){
    fp = fopen(filename,"w");
    if(fp==NULL){
      printf("Can't open file %s, error %d\n",filename,errno);
      terminate(1);
    }

    pf->fp = fp;
    pf->parallel = 0;
    pf->filename = filename;
    
    if( (fprintf(fp,"%d\n",KSPROP_VERSION_NUMBER_V0))==0 ){
      printf("Error in writing version number\n"); terminate(1);
    }
    if( (fprintf(fp,"\"%s\"\n",ph->time_stamp))==0 ){
      printf("Error in writing time stamp\n"); terminate(1);
    }
    
    if( (fprintf(fp,"%d\t%d\t%d\t%d\n",nx,ny,nz,nt))==0 ){
      printf("Error in writing dimensions\n"); terminate(1);
    }

    write_ksprop_info_file(pf);

  } /* if node 0 */

  /* Synchronize */  
  g_sync();
  
  /* Write propagator */
  currentnode = 0;

  for(t=0;t<nt;t++)for(z=0;z<nz;z++)for(y=0;y<ny;y++)for(x=0;x<nx;x++){
    newnode = node_number(x,y,z,t);
    if(newnode != currentnode){	/* switch to another node */
      /**g_sync();**/
      /* tell newnode it's OK to send */
      if( this_node==0 && newnode!=0 )send_field((char *)&pbuf,1,newnode);
      if( this_node==newnode && newnode!=0 )get_field((char *)&pbuf,1,0);
      currentnode=newnode;
    }

    if(this_node==0){

      if(currentnode==0){ 
	i = node_index(x,y,z,t);
	s = &(lattice[i]);
	proppt = (su3_vector *)F_PT(s,prop);
	/* Copy, converting precision if necessary */
	for(src_clr=0; src_clr<3; src_clr++){
	  for(b=0; b<3; b++){
	    pbuf.e[src_clr][b].real = proppt[src_clr].c[b].real;
	    pbuf.e[src_clr][b].imag = proppt[src_clr].c[b].imag;
	  }
	}
      }
      else {
	get_field((char *)&pbuf,sizeof(fsu3_matrix),currentnode);
      }

      for(a=0;a<3;a++)for(b=0;b<3;b++){
	if( (fprintf(fp,"%.7e\t%.7e\n",(double)pbuf.e[a][b].real,
		     (double)pbuf.e[a][b].imag))== EOF){
	  printf("Write error in save_ksprop_ascii\n"); terminate(1);
	}
      }

    }

    else {	/* for all nodes other than node 0 */
      if(this_node==currentnode){
	i = node_index(x,y,z,t);
	s = &(lattice[i]);
	proppt = (su3_vector *)F_PT(s,prop);
	/* Copy, converting precision if necessary */
	for(src_clr=0; src_clr<3; src_clr++){
	  for(b=0; b<3; b++){
	    pbuf.e[src_clr][b].real = proppt[src_clr].c[b].real;
	    pbuf.e[src_clr][b].imag = proppt[src_clr].c[b].imag;
	  }
	}
	send_field((char *)&pbuf,sizeof(fsu3_matrix),0);
      } /* if */
    } /* else */

  } /* end loop over x,y,z,t */
  
  g_sync();
  if(this_node==0){
    fflush(fp);
    printf("Saved propagator to ascii file  %s\n",
	   filename);
    printf("Time stamp %s\n",ph->time_stamp);
    fclose(fp);
    fflush(stdout);
  }

  return pf;

} /* end save_ksprop_ascii() */

