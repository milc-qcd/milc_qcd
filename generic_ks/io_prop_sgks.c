/* io_prop_sgks.c -- reads (maybe) and writes KS quark propagators
   in a site major format compatible with Fermilab running
   MIMD version 6
   S. Gottlieb: Nov 2002 -- binary I/O (based on work by Wingate
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
#include "../include/io_prop_ks.h"

#ifndef HAVE_FSEEKO
#define fseeko fseek
#endif

#define PARALLEL 1
#define SERIAL 0

#define NATURAL_ORDER 0

#undef MAX_BUF_LENGTH
#define MAX_BUF_LENGTH 512

/* The send buffer would normally be the size of a single su3 matrix
   but on we may need to pad for short messages generated
   in serial reads and writes to avoid switch inefficiencies.  In that
   case we define PAD_SEND_BUF on the compilation line to increase
   this.  */
   
#ifndef PAD_SEND_BUF
#define PAD_SEND_BUF 8
#endif

/*----------------------------------------------------------------------*/
/* This subroutine writes the Fermilab compatible propagator header structure */
/* Serial access version */

void swrite_ks_fm_prop_hdr(FILE *fp, ks_prop_header *ksph)
{
  int i;
  int32type size_of_element = sizeof(Real);  /* Real */
  int32type elements_per_site = sizeof(fsu3_matrix)/size_of_element; /* propagator is 3 X 3 complex */
  char myname[] = "swrite_ks_fm_prop_hdr";

  swrite_data(fp,(void *)&ksph->magic_number,sizeof(ksph->magic_number),
	      myname,"magic_number");
/*NEED GMTIME*/
  swrite_data(fp,&ksph->gmtime_stamp,sizeof(ksph->gmtime_stamp), myname,"gmtime_stamp");    
  swrite_data(fp,&size_of_element,sizeof(int32type), myname,"size_of_element"); 
  swrite_data(fp,&elements_per_site,sizeof(int32type), myname,"elements_per_site"); 
  swrite_data(fp,(void *)ksph->dims,sizeof(ksph->dims),
	      myname,"dimensions");
  swrite_data(fp,&ksph->order,sizeof(ksph->order), myname,"site_order");    

  /* Header byte length */

/* includes space for size_of_element and elements_per_site */
  ksph->header_bytes = sizeof(ksph->magic_number) + 
    sizeof(ksph->gmtime_stamp) + 2*sizeof(int32type) +
    sizeof(ksph->dims) + sizeof(ksph->order);

} /* swrite_ks_fmprop_hdr */

/*------------------------------------------------------------------------*/
/* Write a data item to the prop info file THIS DID NOT NEED TO CHANGE*/

/*----------------------------------------------------------------------*/
/* Open, write, and close the ASCII info file */

void write_ks_fmprop_info_file(ks_prop_file *pf)
{
  FILE *info_fp;
  ks_prop_header *ph;
  char info_filename[256];
  char sums[20];
  int32type natural_order = NATURAL_ORDER;
  int32type size_of_element = sizeof(Real);  /* Real */
  int32type elements_per_site = sizeof(fsu3_matrix)/size_of_element; /* propagator is 3 X 3 complex */

  ph = pf->header;

  /* Construct header file name from propagator file name 
   by adding filename extension to propagator file name */

  strcpy(info_filename,pf->filename);
  strcat(info_filename,ASCII_PROP_INFO_EXT);

  /* Open header file */
  
  if((info_fp = fopen(info_filename,"w")) == NULL)
    {
      printf("write_ksprop_info_file: Can't open ascii info file %s\n",
	     info_filename);
      return;
    }
  
  /* Write required information */

  write_ksprop_info_item(info_fp,"header.time_stamp","\"%s\"",ph->time_stamp,0,0);
  /*sprintf(sums,"%x %x",pf->check.sum29,pf->check.sum31);
  write_ksprop_info_item(info_fp,"checksums","\"%s\"",sums,0,0); */
  write_ksprop_info_item(info_fp,"header.dim","%d",(char *)&ph->dims,4,
	sizeof(int32type));
  write_ksprop_info_item(info_fp,"header.magic_number","%d",(char *)&ph->magic_number,0,0);
  write_ksprop_info_item(info_fp,"header.size-of-elment","%d",
		(char *)&size_of_element,0,0);
  write_ksprop_info_item(info_fp,"header.elements-per-site","%d",
		(char *)&elements_per_site,0,0);
/* only natural order is supported now */
  write_ksprop_info_item(info_fp,"header.site-order","%d",
		(char *)&natural_order,0,0);

  write_appl_ksprop_info(info_fp);

  fclose(info_fp);

  printf("Wrote info file %s\n",info_filename); 
  fflush(stdout);

} /*write_ks_fmprop_info_file */

/*----------------------------------------------------------------------*/

/* Set up the output prop file and prop header structure */

ks_prop_file *setup_output_ks_fmprop_file()
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

  /* Load header values */

  ph->magic_number = KSFMPROP_VERSION_NUMBER;

  ph->dims[0] = nx;
  ph->dims[1] = ny;
  ph->dims[2] = nz;
  ph->dims[3] = nt;

  /* Get date and time stamp. (We use local time on node 0) */

  if(this_node==0)
    {
      ph->gmtime_stamp = time(&time_stamp);
/*printf("time_stamp= %d \t gm_time_stamp= %d\n",time_stamp,ph->gmtime_stamp);*/
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

} /* setup_output_ks_fmprop_file */

/*----------------------------------------------------------------------*/

/* Open a binary file for serial writing by node 0. Version for FNAL header */

ks_prop_file *w_serial_ks_fm_i(char *filename)
{
  /* Only node 0 opens the file filename */
  /* Returns a file structure describing the opened file */

  FILE *fp;
  ks_prop_file *kspf;
  ks_prop_header *ksph;

  /* Set up ks_prop file and ks_prop header structs and load header values */
  kspf = setup_output_ks_fmprop_file();
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
      
      swrite_ks_fm_prop_hdr(fp,ksph);

    }
  
  /* Assign values to file structure */

  if(this_node==0) kspf->fp = fp; 
  else kspf->fp = NULL;                /* Only node 0 knows about this file */

  kspf->filename       = filename;
  kspf->byterevflag    = 0;            /* Not used for writing */
  kspf->rank2rcv       = NULL;         /* Not used for writing */
  kspf->parallel       = SERIAL;

  /* Node 0 writes ascii info file */

  if(this_node == 0) write_ks_fmprop_info_file(kspf);

  return kspf;

} /* w_serial_ks_fm_i */

/*---------------------------------------------------------------------------*/

void w_serial_ks_fm(ks_prop_file *kspf, field_offset prop0,
                       field_offset prop1, field_offset prop2)
{
  /* kspf  = file descriptor as opened by w_serial_w_i 
     prop[012]   = field offsets 3 propagators su3 vector (type su3_vector)  */

  FILE *fp;
  ks_prop_header *ksph;
  u_int32type *val;
  int rank29,rank31;
  fsu3_matrix *pbuf;
  int fseek_return;  /* added by S.G. for large file debugging */
  struct {
    fsu3_matrix ksv;
    char pad[PAD_SEND_BUF]; /* Introduced because some switches
			       perform better if message lengths are longer */
  } msg;
  int buf_length;

  register int i,j,k,a,b;
  off_t offset;             /* File stream pointer */
  off_t ks_prop_size;        /* Size of propagator blocks for all nodes */
  off_t ks_prop_check_size;  /* Size of propagator checksum record */
  off_t coord_list_size;    /* Size of coordinate list in bytes */
  off_t head_size;          /* Size of header plus coordinate list */
  off_t body_size ;         /* Size of propagator blocks for all nodes 
			      plus checksum record */
  int currentnode,newnode;
  int x,y,z,t;
  su3_vector *proppt[3];
  site *s;

  if(this_node==0)
    {
      if(kspf->parallel == PARALLEL)
	printf("w_serial_ks: Attempting serial write to file opened in parallel \n");

      pbuf = (fsu3_matrix *)malloc(MAX_BUF_LENGTH*sizeof(su3_matrix));
      if(pbuf == NULL)
	{
	  printf("w_serial_ks: Node 0 can't malloc pbuf\n"); 
	  fflush(stdout); terminate(1);
        }

      fp = kspf->fp;
      ksph = kspf->header;

/* these should no longer be needed
      ks_prop_size = volume*sizeof(fsu3_vector);
      ks_prop_check_size = sizeof(kspf->check.color) +
	sizeof(kspf->check.sum29) + sizeof(kspf->check.sum31);
      body_size = ks_prop_size + ks_prop_check_size;
*/

      /* No coordinate list was written because fields are to be written
	 in standard coordinate list order */
      
      coord_list_size = 0;
      head_size = ksph->header_bytes + coord_list_size;
      
      offset = head_size /*+ body_size*color + ks_prop_check_size*/ ;

      fseek_return=fseeko(fp,offset,SEEK_SET);
      /* printf("w_serial_ks: Node %d fseek_return = %d\n",this_node,fseek_return); */
      if( fseek_return < 0 ) 
	{
	  printf("w_serial_ks: Node %d fseeko %lld failed error %d file %s\n",
		 this_node, (long long)offset, errno, kspf->filename);
	  fflush(stdout); terminate(1);
	}

    } /* end if(this_node==0) */

  /* Buffered algorithm for writing fields in serial order */
  
  /* initialize checksums */
/* NEED TO MODIFY IF CHECKSUM FOR EACH TIME SLICE*/
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
          	s = &(lattice[i]);
          	proppt[0] = (su3_vector *)F_PT(s,prop0);
          	proppt[1] = (su3_vector *)F_PT(s,prop1);
          	proppt[2] = (su3_vector *)F_PT(s,prop2);
		/* Convert precision if necessary */
          	for(a=0; a<3; a++){
            		for(b=0; b<3; b++){
              			msg.ksv.e[a][b].real = proppt[a]->c[b].real;
              			msg.ksv.e[a][b].imag = proppt[a]->c[b].imag;
            		}
	 	}
	    }
	  else
	    {
	      get_field((char *)&msg, sizeof(msg),currentnode);
	    }

	  pbuf[buf_length] = msg.ksv;

	  /* Accumulate checksums - contribution from next site */
	  for(k = 0, val = (u_int32type *)&pbuf[buf_length]; 
	      k < (int)sizeof(fsu3_matrix)/(int)sizeof(int32type); k++, val++)
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
	      
	      if( (int)fwrite(pbuf,sizeof(fsu3_matrix),buf_length,fp) 
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
          	s = &(lattice[i]);
          	proppt[0] = (su3_vector *)F_PT(s,prop0);
          	proppt[1] = (su3_vector *)F_PT(s,prop1);
          	proppt[2] = (su3_vector *)F_PT(s,prop2);
		/* Convert precision here if necessary */
          	for(a=0; a<3; a++){
            		for(b=0; b<3; b++){
              			msg.ksv.e[a][b].real = proppt[a]->c[b].real;
              			msg.ksv.e[a][b].imag = proppt[a]->c[b].imag;
            		}
	 	}
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

      /* NOT DONE NOW NEED FOR EACH TIME SLICE Construct check record */

    }

} /* w_serial_ks_fm */
