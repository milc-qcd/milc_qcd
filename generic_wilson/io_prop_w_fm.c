/********************** io_prop_w_fm.c ***************************/
/* MIMD version 7 */
/* CD 4/07
   For reading and writing Fermilab formatted propagator files
*/

#include "generic_wilson_includes.h"
#include "../include/file_types.h"
#include "../include/io_lat.h"
#include "../include/io_wprop.h"
#include <sys/types.h>
#include <fcntl.h>
#include <errno.h>
#include <string.h>
#include <time.h>
#include <assert.h>

#define PARALLEL 1
#define SERIAL 0

#undef MAX_BUF_LENGTH
#define MAX_BUF_LENGTH 4096

#ifndef PAD_SEND_BUF
#define PAD_SEND_BUF 8
#endif
#define NATURAL_ORDER 0

/*----------------------------------------------------------------------*/
int read_w_fm_prop_hdr(w_prop_file *wpf)
{
  w_prop_header *wph = wpf->header;
  int *dims = wph->dims;

  int32type   tmp; 
  int32type   size_of_element, order; 
  int32type   t_stamp;
  int byterevflag = 0;
  char myname[] = "read_w_fm_prop_hdr";

  if( sizeof(float) != sizeof(int32type)) {
    printf("%s: Can't byte reverse\n", wpf->filename);
    printf("requires size of int32type(%d) = size of float(%d)\n",
	   (int)sizeof(int32type),(int)sizeof(float));
    terminate(1);
  }

  if(sread_data(wpf->fp,&wph->magic_number,sizeof(int32type),
		myname,"magic_no")) terminate(1);

  tmp = wph->magic_number;
  if(wph->magic_number == IO_UNI_MAGIC) byterevflag = 0;
  else 
    {
      byterevn((int32type *)&wph->magic_number,1);
      if(wph->magic_number == IO_UNI_MAGIC) 
	{
	  byterevflag = 1; 
	  /** printf("Reading with byte reversal\n"); **/
	}
      else
	{
	  /* Restore magic number as originally read */
	  wph->magic_number = tmp;
	  
	  /* End of the road. */
	  printf("%s: Unrecognized magic number in prop file header.\n",
		 wpf->filename);
	  printf("Expected %x but read %x\n", IO_UNI_MAGIC,tmp);
	  terminate(1);
	}
    };

  if(sread_byteorder(byterevflag,wpf->fp,&t_stamp,
	     sizeof(t_stamp),myname,"t_stamp"))terminate(1);
  if(sread_byteorder(byterevflag,wpf->fp,&size_of_element,
	     sizeof(int32type),myname,"size_of_element"))terminate(1);
  if(sread_byteorder(byterevflag,wpf->fp,&wph->elements_per_site,
	     sizeof(int32type),myname,"elements_per_site"))terminate(1);
  if(sread_byteorder(byterevflag,wpf->fp,dims,
	     4*sizeof(int32type),myname,"dimensions"))terminate(1);

  if( dims[0]!=nx || dims[1]!=ny || dims[2]!=nz || dims[3]!=nt )
    {
      /* So we can use this routine to discover the dimensions,
	 we provide that if nx = ny = nz = nt = -1 initially
	 we don't die */
      
      if(nx != -1 || ny != -1 || nz != -1 || nt != -1)
	{
	  printf(" Incorrect lattice size %d,%d,%d,%d\n",
		 dims[0], dims[1], dims[2], dims[3]);
	  terminate(1);
	}
    }
  else
    {
      nx = dims[0]; ny = dims[1]; nz = dims[2]; nt = dims[3];
      volume = nx*ny*nz*nt;
    }
  
  /* Consistency checks */
  
  if( size_of_element != sizeof(float) ||
      ( wph->elements_per_site != 288 /* wilson_propagator */ &&
	wph->elements_per_site != 24  /* wilson_vector */ ))
    {	
      printf(" File %s is not a wilson propagator or wilson vector file",
	     wpf->filename);
      printf(" got size_of_element %d and elements_per_site %d\n",
	     size_of_element, wph->elements_per_site);
      terminate(1);
    }
  
  /* The site order parameter is ignored */
  
  if(sread_byteorder(byterevflag,wpf->fp,&order,sizeof(int32type),
		     myname,"order parameter"))terminate(1);
  
  return byterevflag;
}
 
/*----------------------------------------------------------------------*/
/* This subroutine writes the propagator header structure */
/* Serial access version */

void swrite_w_fm_prop_hdr(FILE *fp, w_prop_header *wph)
{

  int32type   elements_per_site = wph->elements_per_site;
  int32type   size_of_element   = sizeof(float);
  int32type   zero32 = 0;
  int32type   t_stamp = wph->t_stamp;

  char myname[] = "swrite_w_fm_prop_hdr";

  swrite_data(fp,(void *)&wph->magic_number,sizeof(wph->magic_number),
	      myname,"magic_number");
  swrite_data(fp,(void *)&t_stamp,sizeof(t_stamp),
	      myname,"t_stamp");
  swrite_data(fp,(void *)&size_of_element,sizeof(size_of_element),
	      myname,"size_of_element");
  swrite_data(fp,(void *)&elements_per_site,sizeof(elements_per_site),
	      myname,"elements_per_site");
  swrite_data(fp,(void *)wph->dims,sizeof(wph->dims),
	      myname,"dimensions");
  swrite_data(fp,&zero32,sizeof(zero32),
	      myname,"order");

  /* Header byte length */

  wph->header_bytes = sizeof(wph->magic_number) + sizeof(wph->time_stamp) 
    + sizeof(size_of_element) + sizeof(elements_per_site) +
    + sizeof(wph->dims) + sizeof(zero32);

} /* swrite_w_fm_prop_hdr */

/*----------------------------------------------------------------------*/
/* Open and  write to the ASCII info file */

static char info_filename[256];

int write_w_fm_prop_info_file(w_prop_file *wpf)
{
  FILE *info_fp;
  w_prop_header *wph = wpf->header;
  int size_of_element = sizeof(float);
  int elements_per_site = wph->elements_per_site;
  int site_order = 0;
  void write_appl_w_prop_info(FILE *fp);

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

  /* Save the file pointer */
  wpf->info_fp = info_fp;
  
  /* Write required information */

  write_w_prop_info_item(info_fp,"header.time_stamp","\"%s\"",
			 wph->time_stamp,0,0);
  write_w_prop_info_item(info_fp,"header.size_of_element","%d",
			 (char *)&size_of_element,0,0);
  write_w_prop_info_item(info_fp,"header.elements_per_site","%d",
			 (char *)&elements_per_site,0,0);
  write_w_prop_info_item(info_fp,"header.dim[0]","%d",(char *)&nx,0,0);
  write_w_prop_info_item(info_fp,"header.dim[1]","%d",(char *)&ny,0,0);
  write_w_prop_info_item(info_fp,"header.dim[2]","%d",(char *)&nz,0,0);
  write_w_prop_info_item(info_fp,"header.dim[3]","%d",(char *)&nt,0,0);
  write_w_prop_info_item(info_fp,"header.site_order","%d",(char *)&site_order,0,0);

/*  printf("Wrote propagator info file %s\n",info_filename); */

  return 0;
} /*write_w_fm_prop_info_file */

/*----------------------------------------------------------------------*/

/* Append one string item to info file */

void append_w_fm_prop_info_file(w_prop_file *wpf, char *keyword, char *src)
{
  FILE *info_fp = wpf->info_fp;
  write_w_prop_info_item(info_fp,keyword,"%s",src,0,0);

} /*append_w_fm_prop_info_file */

/*----------------------------------------------------------------------*/

/* Open a binary file for serial writing by node 0 */
/* Two formats are supported: sc = 0 full propagator field
   sc = 1 twelve dirac vector fields */

w_prop_file *w_serial_w_fm_generic_i(char *filename, int elements_per_site)
{
  /* Only node 0 opens the file filename */
  /* Returns a file structure describing the opened file */

  FILE *fp;
  w_prop_file *wpf;
  w_prop_header *wph;

  /* Set up w_prop file and w_prop header structures and load header values */
  /* First reuse the MILC procedure */
  wpf = setup_output_w_prop_file();
  wpf->file_type = FILE_TYPE_W_FMPROP;
  wph = wpf->header;
  
  /* We want the FNAL magic number */

  wph->magic_number = IO_UNI_MAGIC;

  /* Indicate coordinate natural ordering */

  wph->order = NATURAL_ORDER;

  /* Set elements per site */

  wph->elements_per_site = elements_per_site;

  /* Only node 0 opens the requested file */

  if(this_node == 0)
    {
      fp = g_open(filename, "wb");
      if(fp == NULL)
	{
	  printf("w_serial_w_fm_i: Node %d can't open file %s, error %d\n",
		 this_node,filename,errno);fflush(stdout);
	  terminate(1);
	}

/*      printf("Opened prop file %s for serial writing\n",filename); */
      
      /* Node 0 writes the header */
      
      swrite_w_fm_prop_hdr(fp,wph);
    }
  
  /* Assign values to file structure */

  if(this_node==0)wpf->fp = fp; 
  else wpf->fp = NULL;                /* Only node 0 knows about this file */

  wpf->filename       = filename;
  wpf->byterevflag    = 0;            /* Not used for writing */
  wpf->parallel       = SERIAL;
  wpf->prop           = NULL;
  wpf->info_fp        = NULL;

  /* Node 0 writes ascii info file */

  if(this_node == 0)write_w_fm_prop_info_file(wpf);

  return wpf;

} /* w_serial_w_fm_generic_i */

/*---------------------------------------------------------------------------*/
w_prop_file *w_serial_w_fm_i(char *filename)
{
  return w_serial_w_fm_generic_i(filename, 288);
}

/*---------------------------------------------------------------------------*/
w_prop_file *w_serial_w_fm_sc_i(char *filename)
{
  return w_serial_w_fm_generic_i(filename, 24);
}

/*---------------------------------------------------------------------------*/
/* Copy source propagator to message structure
   converting word order from wilson_propagator to
   wilson_matrix and possibly converting to single
   precision in the process */

void d2f_wmat(wilson_propagator *src, fwilson_matrix *dest){
  int s0,c0,s1,c1;
  for(s0=0;s0<4;s0++)
    for(c0=0;c0<3;c0++)
      for(s1=0;s1<4;s1++)
	for(c1=0;c1<3;c1++)
	  {
	    dest->d[s0].c[c0].d[s1].c[c1].real =
	      src->c[c0].d[s0].d[s1].c[c1].real; 
	    dest->d[s0].c[c0].d[s1].c[c1].imag =
	      src->c[c0].d[s0].d[s1].c[c1].imag; 
	  }
}

/*---------------------------------------------------------------------------*/
/* Here only node 0 writes the full propagator to a serial file in the
   Fermilab format */

void w_serial_w_fm(w_prop_file *wpf, field_offset src_site,
		wilson_propagator *src_field)
{
  /* wpf  = file descriptor as opened by w_serial_w_fm_i 
     src  = field offset for propagator (type wilson_propagator)  */

  FILE *fp = NULL;
  u_int32type *val;
  int rank29=0,rank31=0;
  fwilson_matrix *lbuf = NULL;
  wilson_propagator *src;
  struct {
    fwilson_matrix wm;
    char pad[PAD_SEND_BUF];    /* Introduced because some switches
				  perform better if message lengths are longer */
  } msg;
  int buf_length;
  register int i,j,k;
  int currentnode,newnode;
  int x,y,z,t;

  if(this_node==0)
    {
      lbuf = (fwilson_matrix *)malloc(MAX_BUF_LENGTH*sizeof(fwilson_matrix));
      if(lbuf == NULL)
	{
	  printf("w_serial_w_fm: Node 0 can't malloc lbuf\n"); 
	  fflush(stdout);terminate(1);
        }

      fp = wpf->fp;
    }
      
  /* Buffered algorithm for writing fields in serial order */
  
  g_sync();
  currentnode=0;

  buf_length = 0;

  for(j=0,t=0;t<nt;t++)for(z=0;z<nz;z++)for(y=0;y<ny;y++)for(x=0;x<nx;x++,j++)
    {
      newnode=node_number(x,y,z,t);
      if(newnode != currentnode){	/* switch to another node */
	/* Send a few bytes of garbage to tell newnode it's OK to send */
	if( this_node==0 && newnode!=0 )
	  send_field((char *)&msg,16,newnode);
	if( this_node==newnode && newnode!=0 )
	  get_field((char *)&msg,16,0);
	currentnode=newnode;
      }
      
      /* All nodes initialize timeslice checksums at the beginning of
	 a time slice */
      if(x == 0 && y == 0 && z == 0)
	{
	  wpf->check.sum31 = 0;
	  wpf->check.sum29 = 0;
	  /* counts 32-bit words mod 29 and mod 31 in order of appearance
	     on file */
	  /* Here all nodes see the same sequence because we read serially */
	  rank29 = 0;
	  rank31 = 0;
	}

      if(this_node==0)
	{
	  if(currentnode == 0)
	    {
	      i=node_index(x,y,z,t);
	      /* Load msg structure with single precision wilson propagator */
	      
	      if(src_site == (field_offset)(-1))
		src = src_field + i;
	      else
		src = (wilson_propagator *)F_PT( &(lattice[i]), src_site );
	      
	      d2f_wmat(src, &msg.wm);
	    }
	  else
	    {
	      get_field((char *)&msg, sizeof(msg),currentnode);
	    }
	  
	  lbuf[buf_length] = msg.wm;

	  /* Accumulate checksums - contribution from next site */
	  for(k = 0, val = (u_int32type *)&lbuf[buf_length]; 
	      k < (int)sizeof(fwilson_matrix)/(int)sizeof(int32type); k++, val++)
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
	      
	      if( (int)g_write(lbuf,sizeof(fwilson_matrix),buf_length,fp) 
		  != buf_length)
		{
		  printf("w_serial_w_fm: Node %d propagator write error %d file %s\n",
			 this_node,errno,wpf->filename); 
		  fflush(stdout);
		  terminate(1);   
		}
	      buf_length = 0;		/* start again after write */
	    }
	} /* this_node == 0 */
      else  /* for nodes other than 0 */
	{
	  if(this_node == currentnode){
	    i=node_index(x,y,z,t);
	    /* Copy or convert data into send buffer and send to node
	       0 with padding */
	    
	    if(src_site == (field_offset)(-1))
	      src = src_field + i;
	    else
	      src = (wilson_propagator *)F_PT( &(lattice[i]), src_site );
	    
	    d2f_wmat(src, &msg.wm);
	    send_field((char *)&msg, sizeof(msg),0);
	  }
	}
      
      /* Accumulate and print checksum at the end of each time slice */
      if(x == nx - 1 && y == ny - 1 && z == nz - 1)
	{
	  /* Combine node checksum contributions with global exclusive or */
	  g_xor32(&wpf->check.sum29);
	  g_xor32(&wpf->check.sum31);
	  
	  if(this_node == 0){
	    char keyword[30];
	    char sums[30];
	    snprintf(keyword,30,"quark.t[%d].checksum",t);
	    snprintf(sums,30,"\"0x%0x 0x%0x\"",
		     wpf->check.sum29,wpf->check.sum31);
	    append_w_fm_prop_info_file(wpf,keyword,sums);
	    printf("quark.t[%d].checksum  \"%0x %0x\"\n",t,
		   wpf->check.sum29, wpf->check.sum31);
	    wpf->check.sum31 = 0;
	    wpf->check.sum29 = 0;
	  }
	}
    } /*close x,y,z,t loops */
  
  g_sync();
  
  if(this_node==0)
    free(lbuf);
      
} /* w_serial_w_fm */

/*---------------------------------------------------------------------------*/
/* Here only node 0 writes the full propagator to a serial file in the
   Fermilab sc format */

void w_serial_w_fm_sc(w_prop_file *wpf, field_offset src_site,
		      wilson_propagator *src_field)
{
  /* wpf  = file descriptor as opened by w_serial_w_fm_i 
     src  = field offset for propagator (type wilson_propagator)  */

  FILE *fp = NULL;
  u_int32type *val;
  int rank29,rank31;
  fwilson_vector *lbuf = NULL;
  wilson_propagator *src;
  struct {
    fwilson_vector wm;
    char pad[PAD_SEND_BUF];    /* Introduced because some switches
				  perform better if message lengths are longer */
  } msg;
  int buf_length;
  register int i,j,k;
  int currentnode,newnode;
  int x,y,z,t;
  int color,spin,c1,s1;
  int irecord;
  if(this_node==0)
    {
      lbuf = (fwilson_vector *)malloc(MAX_BUF_LENGTH*sizeof(fwilson_vector));
      if(lbuf == NULL)
	{
	  printf("w_serial_w_fm: Node 0 can't malloc lbuf\n"); 
	  fflush(stdout);terminate(1);
        }

      fp = wpf->fp;
    }
      
  /* Buffered algorithm for writing fields in serial order */

  g_sync();
  currentnode=0;


  for(irecord=0,spin=0;spin<4;spin++)for(color=0;color<3;color++,irecord++){
      /* initialize checksums for this color and spin */

      buf_length = 0;

      rank29 = 0;
      rank31 = 0;

      wpf->check.sum31 = 0;
      wpf->check.sum29 = 0;

      for(j=0,t=0;t<nt;t++)for(z=0;z<nz;z++)for(y=0;y<ny;y++)for(x=0;x<nx;x++,j++)
     {
       newnode=node_number(x,y,z,t);
       if(newnode != currentnode){	/* switch to another node */
	 /* Send a few bytes of garbage to tell newnode it's OK to send */
	 if( this_node==0 && newnode!=0 )
	   send_field((char *)&msg,16,newnode);
	 if( this_node==newnode && newnode!=0 )
	   get_field((char *)&msg,16,0);
	 currentnode=newnode;
       }
       
       if(this_node == currentnode){
	 
	 i=node_index(x,y,z,t);
	 /* Load msg structure with single precision wilson propagator */
	 
	 if(src_site == (field_offset)(-1))
	   src = src_field + i;
	 else
	   src = (wilson_propagator *)F_PT( &(lattice[i]), src_site );
	 
	 /* Copy source propagator to message structure extracting the
	    wilson vector and possibly converting to single precision
	    in the process */
	 
	 for(s1=0;s1<4;s1++)
	   for(c1=0;c1<3;c1++)
	     {
	       msg.wm.d[s1].c[c1].real =
		 src->c[color].d[spin].d[s1].c[c1].real; 
	       msg.wm.d[s1].c[c1].imag =
		 src->c[color].d[spin].d[s1].c[c1].imag; 
	     }
	 /* Send the converted message to node 0 */
	 if(this_node != 0)
	   send_field((char *)&msg, sizeof(msg),0);
       }
       
       if(this_node==0)
	 {
	   if(currentnode!=0)
	     {
	       get_field((char *)&msg, sizeof(msg),currentnode);
	     }
	   
	   lbuf[buf_length] = msg.wm;
	   
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
	       
	       if( (int)g_write(lbuf,sizeof(fwilson_vector),buf_length,fp) 
		   != buf_length)
		 {
		   printf("w_serial_w_fm: Node %d propagator write error %d file %s\n",
			  this_node,errno,wpf->filename); 
		   fflush(stdout);
		   terminate(1);   
		 }
	       buf_length = 0;		/* start again after write */
	     }
	 }
       
     } /*close x,y,z,t loops */

      /* Print checksum at the end of each source spin/color */
      
      if(this_node == 0){
	char keyword[30];
	char sums[30];
	snprintf(keyword,30,"record[%d].checksum",irecord);
	snprintf(sums,30,"\"0x%0x 0x%0x\"",wpf->check.sum29,wpf->check.sum31);
	append_w_fm_prop_info_file(wpf,keyword,sums);
	printf("record[%d].checksum  \"0x%0x 0x%0x\" \n",irecord,
	       wpf->check.sum29,wpf->check.sum31);
      }
    } /* source spin, color */
  
  g_sync();
  
  if(this_node==0)
    free(lbuf);
      
} /* w_serial_w_fm_sc */

/*---------------------------------------------------------------------------*/
/* Write propagator, according to format  */

void w_prop_w_fm_generic(w_prop_file *wpf, field_offset src_site,
			 wilson_propagator *src_field)
{
  /* wpf  = file descriptor as opened by w_serial_w_fm_i 
     src  = field offset for propagator (type wilson_propagator)  */

  if(wpf->header->elements_per_site == 288){
    /* Full propagator per site */
    w_serial_w_fm(wpf, src_site, src_field);
  }
  else if(wpf->header->elements_per_site == 24){
    /* Twelve records. One Wilson vector per site */
    w_serial_w_fm_sc(wpf, src_site, src_field);
  }
  else {
    printf("w_prop_w_fm: Bad elements_per_site %d in file header struct.\n",
	   wpf->header->elements_per_site);
    terminate(1);
  }
}

/*--------------------------------------------------------------------*/
/* Write according to format */
void w_serial_w_fm_from_site(w_prop_file *wpf, field_offset src_site)
{
  w_prop_w_fm_generic(wpf, src_site, NULL);
}

/*--------------------------------------------------------------------*/
/* Write, according to format */
void w_serial_w_fm_from_field(w_prop_file *wpf, wilson_propagator *src_field)
{
  w_prop_w_fm_generic(wpf, (field_offset)(-1), src_field);
}

/*---------------------------------------------------------------------------*/

void w_serial_w_fm_f(w_prop_file *wpf)

/* this subroutine closes the file and frees associated structures */
{
  g_sync();
  if(this_node==0)
    {
      if(wpf->parallel == PARALLEL)
	printf("w_serial_w_fm_f: Attempting serial close on file opened in parallel \n");

      printf("Wrote prop file %s time stamp %s\n",wpf->filename,
	     wpf->header->time_stamp);

      if(wpf->fp != NULL)g_close(wpf->fp);
    }

  /* Close the info file */
  if(wpf->info_fp != NULL)
    g_close(wpf->info_fp);

  free(wpf->header);
  free(wpf);

} /* w_serial_w_fm_f */

/*---------------------------------------------------------------------------*/

w_prop_file *r_serial_w_fm_i(char *filename)
{
  /* Returns file descriptor for opened file */

  w_prop_file *wpf;
  w_prop_header *wph;
  int byterevflag;

  wpf = setup_input_w_prop_file(filename);
  wph = wpf->header;

  if(this_node==0){
    wpf->fp = g_open(filename,"rb");
    if(wpf->fp==NULL){
      printf("Can't open propagator file %s, error %d\n",filename,errno);
      terminate(1);
    }
    byterevflag = read_w_fm_prop_hdr(wpf);

  }
  else wpf->fp = NULL;

  /* Broadcast the byterevflag from node 0 to all nodes */
  broadcast_bytes((char *)&byterevflag,sizeof(byterevflag));
  wpf->byterevflag = byterevflag;
  
  /* Assign other default values to file structure */
  wpf->prop           = NULL;
  wpf->info_fp        = NULL;

  /* Node 0 broadcasts the header structure to all nodes */
  broadcast_bytes((char *)wph,sizeof(w_prop_header));

  return wpf;

} /* r_serial_w_fm_i */


/* Read the entire file */

void r_serial_w_fm(w_prop_file *wpf, field_offset dest_site, 
		   wilson_propagator *dest_field)
{
  int rcv_rank, rcv_coords, status;
  int destnode;
  int x,y,z,t,i=0, byterevflag, c0,s0,c1,s1,a;
  struct {
    fwilson_matrix q;
    char pad[PAD_SEND_BUF];    /* Introduced because some switches
				  perform better if message lengths
				  are longer */
  } msg;

  int buf_length=0, where_in_buf=0;
  fwilson_matrix *pbuff=NULL;
  wilson_propagator *qp;
  w_prop_check test_wpc;
  u_int32type *val;
  int rank29=0,rank31=0;
  int k;

  byterevflag = wpf->byterevflag;

  if(this_node == 0)
    {
      pbuff = (fwilson_matrix *)malloc(MAX_BUF_LENGTH*sizeof(fwilson_matrix));
      if(pbuff == NULL)
	{
	  printf("Node %d can't malloc pbuff\n",this_node);
	  fflush(stdout);
	  terminate(1);
	}
      
      buf_length = 0;
      where_in_buf = 0;

    } /* end of if(this_node == 0)*/

  g_sync();

  /* Node 0 reads and deals out the values */
  status = 0;
  for(rcv_rank=0; rcv_rank<volume; rcv_rank++)
    {
      /* We do only natural order here */
      rcv_coords = rcv_rank;
      
      x = rcv_coords % nx;   rcv_coords /= nx;
      y = rcv_coords % ny;   rcv_coords /= ny;
      z = rcv_coords % nz;   rcv_coords /= nz;
      t = rcv_coords % nt;

      /* All nodes initialize timeslice checksums at the beginning of
	 a time slice */
      if(x == 0 && y == 0 && z == 0)
	{
	  test_wpc.sum31 = 0;
	  test_wpc.sum29 = 0;
	  /* counts 32-bit words mod 29 and mod 31 in order of appearance
	     on file */
	  /* Here all nodes see the same sequence because we read serially */
	  rank29 = 0;
	  rank31 = 0;
	}

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
	    a=(int)g_read(pbuff,sizeof(fwilson_matrix),buf_length,wpf->fp);
	    
	    if( a  != buf_length)
	      {
		
		if(status == 0)
		  printf(" node %d propagator read error %d file %s\n",
			 this_node, errno, wpf->filename); 
		fflush(stdout); 
		status = 1;
	      }
	    where_in_buf = 0;  /* reset counter */
	  }  /*** end of the buffer read ****/
	
	/* Save Wilson matrix in msg.q for further processing */
	msg.q = pbuff[where_in_buf];

	if(destnode==0){	/* just copy su3_matrix */
	  i = node_index(x,y,z,t);
	}
	else {		        /* send to correct node */
	  send_field((char *)&msg, sizeof(msg), destnode);
	}
	where_in_buf++;
      }
       /* if(this_node==0) */
      else {	/* for all nodes other than node 0 */
	if(this_node==destnode){
	  i = node_index(x,y,z,t);
	  get_field((char *)&msg, sizeof(msg),0);
	}
      }

      /* The receiving node does the byte reversal.  At this point msg
	 contains the input vectors and i points to the destination
	 site structure */
      
      if(this_node==destnode)
	{
	  if(byterevflag==1)
	    byterevn((int32type *)&msg.q, 
		     sizeof(fwilson_matrix)/sizeof(int32type));
	  
	  /* Accumulate checksums */
	  for(k = 0, val = (u_int32type *)(&msg.q); 
	      k < (int)sizeof(fwilson_matrix)/(int)sizeof(int32type); 
	      k++, val++)
	    {
	      test_wpc.sum29 ^= (*val)<<rank29 | (*val)>>(32-rank29);
	      test_wpc.sum31 ^= (*val)<<rank31 | (*val)>>(32-rank31);
	      rank29++; if(rank29 >= 29)rank29 = 0;
	      rank31++; if(rank31 >= 31)rank31 = 0;
	    }

	  /* Now copy the site data into the site structure or field
	     switching to propagator storage and converting to generic
	     precision if needed */

	  if(dest_site == (field_offset)(-1))
	    qp = dest_field + i;
	  else
	    qp = (wilson_propagator *)F_PT(&lattice[i],dest_site);

	  for(s0=0;s0<4;s0++)for(c0=0;c0<3;c0++)
	    for(s1=0;s1<4;s1++)for(c1=0;c1<3;c1++)
	      {
		qp->c[c0].d[s0].d[s1].c[c1].real 
		  = msg.q.d[s0].c[c0].d[s1].c[c1].real;
		qp->c[c0].d[s0].d[s1].c[c1].imag 
		  = msg.q.d[s0].c[c0].d[s1].c[c1].imag;
	      }
	}
      else
	{
	  rank29 += sizeof(fwilson_matrix)/sizeof(int32type);
	  rank31 += sizeof(fwilson_matrix)/sizeof(int32type);
	  rank29 %= 29;
	  rank31 %= 31;
	}

      /* Accumulate and print checksum at the end of each time slice */
      if(x == nx - 1 && y == ny - 1 && z == nz - 1)
	{
	  /* Combine node checksum contributions with global exclusive or */
	  g_xor32(&test_wpc.sum29);
	  g_xor32(&test_wpc.sum31);
	  
	  node0_printf("quark.t[%d].checksum  \"%0x %0x\"\n",t,
		       test_wpc.sum29, test_wpc.sum31);
	}
    } /* rcv_rank */

  if(pbuff != NULL)free(pbuff);  pbuff = NULL;

  /**  if(this_node==0)
    {
      printf("Read Wilson prop serially from file %s\n", wpf->filename);
      }**/
}

/* We read the entire propagator file, even though it is sorted
   according to source spin-color */

void r_serial_w_fm_sc(w_prop_file *wpf, field_offset dest_site, 
		      wilson_propagator *dest_field)
{
  int rcv_rank, rcv_coords, status;
  int destnode;
  int x,y,z,t,i=0, byterevflag, c1,s1,a;
  struct {
    fwilson_vector q;
    char pad[PAD_SEND_BUF];    /* Introduced because some switches
				  perform better if message lengths
				  are longer */
  } msg;

  int buf_length=0, where_in_buf=0;
  fwilson_vector *pbuff=NULL;
  wilson_propagator *qp;
  w_prop_check test_wpc;
  u_int32type *val;
  int rank29=0,rank31=0;
  int k;
  int spin, color;
  int irecord;

  byterevflag = wpf->byterevflag;

  if(this_node == 0)
    {
      pbuff = (fwilson_vector *)malloc(MAX_BUF_LENGTH*sizeof(fwilson_vector));
      if(pbuff == NULL)
	{
	  printf("Node %d can't malloc pbuff\n",this_node);
	  fflush(stdout);
	  terminate(1);
	}
    } /* end of if(this_node == 0)*/

  g_sync();

  /* Node 0 reads and deals out the values */
  status = 0;
  for(irecord=0,spin=0;spin<4;spin++)
    for(color=0;color<3;color++,irecord++){
      
      buf_length = 0;
      where_in_buf = 0;

      test_wpc.sum31 = 0;
      test_wpc.sum29 = 0;
      
      rank29 = 0;
      rank31 = 0;
	  
      for(rcv_rank=0; rcv_rank<volume; rcv_rank++)
	{
	  /* We do only natural order here */
	  rcv_coords = rcv_rank;
	  
	  x = rcv_coords % nx;   rcv_coords /= nx;
	  y = rcv_coords % ny;   rcv_coords /= ny;
	  z = rcv_coords % nz;   rcv_coords /= nz;
	  t = rcv_coords % nt;
	  
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
		a=(int)g_read(pbuff,sizeof(fwilson_vector),buf_length,wpf->fp);
		
		if( a  != buf_length)
		  {
		    
		    if(status == 0)
		      printf(" node %d propagator read error %d file %s\n",
			     this_node, errno, wpf->filename); 
		    fflush(stdout); 
		    status = 1;
		  }
		where_in_buf = 0;  /* reset counter */
	      }  /*** end of the buffer read ****/
	    
	    /* Save Wilson vector in msg.q for further processing */
	    msg.q = pbuff[where_in_buf];
	    
	    if(destnode==0){	/* just copy su3_vector */
	      i = node_index(x,y,z,t);
	    }
	    else {		        /* send to correct node */
	      send_field((char *)&msg, sizeof(msg), destnode);
	    }
	    where_in_buf++;
	  }
	  /* if(this_node==0) */
	  else {	/* for all nodes other than node 0 */
	    if(this_node==destnode){
	      i = node_index(x,y,z,t);
	      get_field((char *)&msg, sizeof(msg),0);
	    }
	  }
	  
	  /* The receiving node does the byte reversal.  At this point msg
	     contains the input vectors and i points to the destination
	     site structure */
	  
	  if(this_node==destnode)
	    {
	      if(byterevflag==1)
		byterevn((int32type *)&msg.q, 
			 sizeof(fwilson_vector)/sizeof(int32type));
	      
	      /* Accumulate checksums */
	      for(k = 0, val = (u_int32type *)(&msg.q); 
		  k < (int)sizeof(fwilson_vector)/(int)sizeof(int32type); 
		  k++, val++)
		{
		  test_wpc.sum29 ^= (*val)<<rank29 | (*val)>>(32-rank29);
		  test_wpc.sum31 ^= (*val)<<rank31 | (*val)>>(32-rank31);
		  rank29++; if(rank29 >= 29)rank29 = 0;
		  rank31++; if(rank31 >= 31)rank31 = 0;
		}
	      
	      /* Now copy the site data into the site structure or field
		 switching to propagator storage and converting to generic
		 precision if needed */
	      
	      if(dest_site == (field_offset)(-1))
		qp = dest_field + i;
	      else
		qp = (wilson_propagator *)F_PT(&lattice[i],dest_site);
	      
	      for(s1=0;s1<4;s1++)
		for(c1=0;c1<3;c1++)
		  {
		    qp->c[color].d[spin].d[s1].c[c1].real 
		      = msg.q.d[s1].c[c1].real;
		    qp->c[color].d[spin].d[s1].c[c1].imag 
		      = msg.q.d[s1].c[c1].imag;
		  }
	    }
	  else
	    {
	      rank29 += sizeof(fwilson_vector)/sizeof(int32type);
	      rank31 += sizeof(fwilson_vector)/sizeof(int32type);
	      rank29 %= 29;
	      rank31 %= 31;
	    }
	} /* rcv_rank */

      /* Combine node checksum contributions with global exclusive or */
      g_xor32(&test_wpc.sum29);
      g_xor32(&test_wpc.sum31);

      node0_printf("record[%d].checksum  \"%0x %0x\"\n",irecord,
		   test_wpc.sum29, test_wpc.sum31);
    } /* source color, spin */
  
  if(pbuff != NULL)free(pbuff);  pbuff = NULL;

  /**  if(this_node==0)
    {
      printf("Read Wilson prop serially from file %s\n", wpf->filename);
      }**/
} /* r_serial_w_fm_sc */

/*----------------------------------------------------------------------*/

void r_serial_w_fm_f(w_prop_file *wpf)

/* Close the file and free associated structures */
{
  g_sync();
  if(this_node==0)
    {
      if(wpf->fp != NULL)g_close(wpf->fp);
      fflush(stdout);
    }
  
  free(wpf->header);
  free(wpf);
  
} /* r_serial_w_fm_f */

/*--------------------------------------------------------------------*/
/* Read according to format */
void r_serial_w_fm_generic(w_prop_file *wpf, field_offset dest_site, 
			   wilson_propagator *dest_field)
{
  /* Read according to file format */
  if(wpf->header->elements_per_site == 288){
    /* One record: Each site has a full propagator matrix */
    r_serial_w_fm(wpf, dest_site, dest_field);
  }
  else if (wpf->header->elements_per_site == 24){
    /* Twelve records, one for each source spin and color */
    r_serial_w_fm_sc(wpf, dest_site, dest_field);
  }
  else {
    printf("r_serial_w_fm_generic: Bad elements per site in file header struct.");
    terminate(1);
  }
}

/*--------------------------------------------------------------------*/
/* Read according to format */
void r_serial_w_fm_to_site(w_prop_file *wpf, field_offset dest_site)
{
  r_serial_w_fm_generic(wpf, dest_site, NULL);
}

/*--------------------------------------------------------------------*/
/* Read, according to format */
void r_serial_w_fm_to_field(w_prop_file *wpf, wilson_propagator *dest_field)
{
  r_serial_w_fm_generic(wpf, (field_offset)(-1), dest_field);
}

/*--------------------------------------------------------------------*/
/* Open, read, close */
void r_prop_w_fm(char *filename, field_offset dest_site, 
		 wilson_propagator *dest_field)
{
  w_prop_file *wpf;
  wpf = r_serial_w_fm_i(filename);
  r_serial_w_fm_generic(wpf, dest_site, dest_field);
  r_serial_w_fm_f(wpf);
  
}

/*--------------------------------------------------------------------*/
/* Open, read, close */
void r_prop_w_fm_to_site(char *filename, field_offset dest_site)
{
  r_prop_w_fm(filename, dest_site, NULL);
}

/*--------------------------------------------------------------------*/
/* Open, read, close */
void r_prop_w_fm_to_field(char *filename, wilson_propagator *dest_field)
{
  r_prop_w_fm(filename, (field_offset)(-1), dest_field);
}

/*--------------------------------------------------------------------*/
