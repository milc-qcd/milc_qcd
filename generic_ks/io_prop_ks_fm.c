/********************** io_prop_ks_fm.c *********************************/
/* MIMD version 7 */
/* Reads and write KS quark propagators
   in a site major format compatible with Fermilab running
   (previous versions variously called io_prop_sgks.c io_prop_ks_write.c
    io_prop_sgks_double.c)
   S. Gottlieb: Nov 2002 -- binary I/O (based on work by Wingate)
   C. DeTar: Mar 2005 Support for internal double precision with
                      single precision binary data.
   C. DeTar: Mar 2005 Support for byte reversal.
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
#define MAX_BUF_LENGTH 512

/* The send buffer would normally be the size of a single SU(3) vector
   but we may need to pad for short messages generated in serial reads
   and writes to avoid switch inefficiencies.  In that case we define
   PAD_SEND_BUF on the compilation line to increase this.  */
   
#ifndef PAD_SEND_BUF
#define PAD_SEND_BUF 8
#endif


/*----------------------------------------------------------------------*/

ks_prop_file *create_input_ks_fmprop_file_handle(char *filename)
{
  ks_prop_file *pf;
  ks_prop_header *ph;
  char myname[] = "create_input_ks_fmprop_file_handle";

  /* Allocate space for the file structure */

  pf = (ks_prop_file *)malloc(sizeof(ks_prop_file));
  if(pf == NULL)
    {
      printf("%s: Can't malloc pf\n", myname);
      terminate(1);
    }

  pf->filename = filename;

  /* Allocate space for the header */

  /* Make sure compilation gave us a 32 bit integer type */
  assert(sizeof(int32type) == 4);

  ph = (ks_prop_header *)malloc(sizeof(ks_prop_header));
  if(ph == NULL)
    {
      printf("%s: Can't malloc ph\n", myname);
      terminate(1);
    }

  pf->header = ph;
  pf->check.sum29 = 0;
  pf->check.sum31 = 0;
  pf->file_type   = FILE_TYPE_KS_FMPROP;
#ifdef HAVE_QIO
  pf->infile       = NULL;
  pf->outfile      = NULL;
#endif

  return pf;
} /* create_input_ks_fmprop_file_handle*/


void destroy_ks_fmprop_file_handle(ks_prop_file *kspf){
  if(kspf == NULL)return;

  if(kspf->header != NULL)
    free(kspf->header);
  if(kspf->info != NULL)
    free(kspf->info);
  free(kspf);

}

/*----------------------------------------------------------------------*/
/* This subroutine writes the Fermilab compatible propagator header structure */
/* Serial access version */

void swrite_ks_fm_prop_hdr(FILE *fp, ks_prop_header *ksph)
{
  int32type size_of_element = sizeof(float);  /* float */
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

  /*  printf ( "swrite_ks_fmprop_hdr: sizeof(header) = %d\n", ksph->header_bytes ); */
} /* swrite_ks_fmprop_hdr */

/*------------------------------------------------------------------------*/
/* Write a data item to the prop info file THIS DID NOT NEED TO CHANGE*/

/*----------------------------------------------------------------------*/
/* Open, write, and close the ASCII info file */

void open_write_ks_fmprop_info_file(ks_prop_file *pf)
{
  FILE *info_fp;
  char info_filename[FILENAME_MAX];

  /* If file is already open, do nothing */
  if(pf->info_fp != NULL)return;

  /* Construct header file name from propagator file name 
     by adding filename extension to propagator file name */

  strcpy(info_filename,pf->filename);
  strcat(info_filename,ASCII_INFO_EXT);

  /* Open header file */
  
  info_fp = fopen(info_filename, "w");

  /* Save info file pointer in ks_prop_file structure */
  pf->info_fp = info_fp;

  if(info_fp == NULL)
    {
      printf("open_write_ksprop_fminfo_file: Can't open ascii info file %s\n",
	     info_filename);
    }

  /*  node0_printf("Writing info file %s\n",info_filename);  */
  fflush(stdout);
}
  
void close_write_ks_fmprop_info_file(ks_prop_file *pf)
{
  FILE *info_fp = pf->info_fp;

  if(info_fp == NULL)return;
  fclose(info_fp);

  pf->info_fp = NULL;
}

/* Open, write, and close the ASCII info file */

void write_ks_fmprop_info_file(ks_prop_file *pf)
{
  FILE *info_fp;
  ks_prop_header *ph;
  int32type natural_order = NATURAL_ORDER;
  int32type size_of_element = sizeof(float);  /* float */
  int32type elements_per_site = sizeof(fsu3_matrix)/size_of_element; /* propagator is 3 X 3 complex */

  ph = pf->header;

  /* Open header file */
  
  open_write_ks_fmprop_info_file(pf);
  info_fp = pf->info_fp;
  if(info_fp == NULL)
    return;

  /* Write required information */

  write_ksprop_info_item(info_fp,"ksprop.time_stamp","\"%s\"",
			 ph->time_stamp,0,0);
  /*sprintf(sums,"%x %x",pf->check.sum29,pf->check.sum31);
  write_ksprop_info_item(info_fp,"checksums","\"%s\"",sums,0,0); */
  write_ksprop_info_item(info_fp,"ksprop.dim","%d",(char *)&ph->dims,4,
	sizeof(int32type));
  write_ksprop_info_item(info_fp,"ksprop.magic_number","%d",(char *)&ph->magic_number,0,0);
  write_ksprop_info_item(info_fp,"ksprop.size-of-element","%d",
		(char *)&size_of_element,0,0);
  write_ksprop_info_item(info_fp,"ksprop.elements-per-site","%d",
		(char *)&elements_per_site,0,0);
/* only natural order is supported now */
  write_ksprop_info_item(info_fp,"ksprop.site-order","%d",
		(char *)&natural_order,0,0);

  write_appl_ksprop_info(info_fp);

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
  pf->check.color = 0;
  pf->check.sum29 = 0;
  pf->check.sum31 = 0;
  pf->file_type   = FILE_TYPE_KS_FMPROP;

  /* Initialize pointer to info file */
  pf->info_fp = NULL;
  pf->info    = NULL;

  /* Load header values */

  ph->magic_number = IO_UNI_MAGIC;

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

  FILE *fp = NULL;
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
  kspf->parallel       = SERIAL;

  /* Node 0 writes ascii info file */

  if(this_node == 0) write_ks_fmprop_info_file(kspf);

  return kspf;

} /* w_serial_ks_fm_i */

/*---------------------------------------------------------------------------*/

void w_serial_ks_fm(ks_prop_file *kspf, field_offset src_site, 
		    su3_vector *src_field)
{
  /* kspf  = file descriptor as opened by ks_serial_w_i 
     src_site[3]   = field offset of an array of three su3_vector types  */

  FILE *fp = NULL;
  ks_prop_header *ksph;
  u_int32type *val;
  int rank29 = 0,rank31 = 0;
  fsu3_matrix *pbuf = NULL;
  int fseek_return;  /* added by S.G. for large file debugging */
  struct {
    fsu3_matrix ksv;
    char pad[PAD_SEND_BUF]; /* Introduced because some switches
			       perform better if message lengths are longer */
  } msg;
  int buf_length;

  register int i,j,k,a,b;
  off_t offset;             /* File stream pointer */
  off_t coord_list_size;    /* Size of coordinate list in bytes */
  off_t head_size;          /* Size of header plus coordinate list */

  int currentnode,newnode;
  int x,y,z,t;
  su3_vector *proppt;
  FILE *info_fp = kspf->info_fp;

  if(this_node==0)
    {
      if(kspf->parallel == PARALLEL)
	printf("w_serial_ks: Attempting serial write to file opened in parallel \n");

      pbuf = (fsu3_matrix *)malloc(MAX_BUF_LENGTH*sizeof(fsu3_matrix));
      if(pbuf == NULL)
	{
	  printf("w_serial_ks: Node 0 can't malloc pbuf\n"); 
	  fflush(stdout); terminate(1);
        }

      fp = kspf->fp;
      ksph = kspf->header;

      /* No coordinate list was written because fields are to be written
	 in standard coordinate list order */
      
      coord_list_size = 0;
      head_size = ksph->header_bytes + coord_list_size;
      
      offset = head_size;

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
  
  g_sync();
  currentnode=0;

  buf_length = 0;

  for(j=0,t=0;t<nt;t++)for(z=0;z<nz;z++)for(y=0;y<ny;y++)for(x=0;x<nx;x++,j++)
    {

      /* All nodes initialize timeslice checksums at the beginning of
	 a time slice */
      if(x == 0 && y == 0 && z == 0)
	{
	  kspf->check.sum31 = 0;
	  kspf->check.sum29 = 0;
	  /* counts 32-bit words mod 29 and mod 31 in order of appearance
	     on file */
	  /* Here all nodes see the same sequence because we read serially */
	  rank29 = 0;
	  rank31 = 0;
	}

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

	      if(src_site == (field_offset)(-1))
		proppt = src_field + 3*i;
	      else
		proppt = (su3_vector *)F_PT(&lattice[i],src_site);

	      /* Convert precision if necessary */
	      for(a=0; a<3; a++){
		for(b=0; b<3; b++){
		  msg.ksv.e[a][b].real = proppt[a].c[b].real;
		  msg.ksv.e[a][b].imag = proppt[a].c[b].imag;
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

		if(src_site == (field_offset)(-1))
		  proppt = src_field + 3*i;
		else
		  proppt = (su3_vector *)F_PT(&lattice[i],src_site);

		/* Convert precision here if necessary */
          	for(a=0; a<3; a++){
            		for(b=0; b<3; b++){
              			msg.ksv.e[a][b].real = proppt[a].c[b].real;
              			msg.ksv.e[a][b].imag = proppt[a].c[b].imag;
            		}
	 	}
	    send_field((char *)&msg, sizeof(msg),0);
	  }
	}
      
      /* Accumulate and print checksum at the end of each time slice */
      if(x == nx - 1 && y == ny - 1 && z == nz - 1)
	{
	  /* Combine node checksum contributions with global exclusive or */
	  g_xor32(&kspf->check.sum29);
	  g_xor32(&kspf->check.sum31);
	  
	  if(this_node == 0){
	    char sums[30]; char key[20];
	    sprintf(key,"quark.t[%d].checksum ",t);
	    sprintf(sums,"%0x %0x",kspf->check.sum29,kspf->check.sum31);
	    write_ksprop_info_item(info_fp,key,"\"%s\"",sums,0,0);
	  }
	}

    } /*close x,y,z,t loops */
  
  g_sync();
  
  if(this_node==0)
    {
      close_write_ks_fmprop_info_file(kspf);
      /*      printf("Wrote KS prop serially to file %s\n", kspf->filename);  */
      fflush(stdout);
      free(pbuf);
    }
  
} /* w_serial_ks_fm */

/*---------------------------------------------------------------------------*/

void w_serial_ks_fm_from_site(ks_prop_file *kspf, field_offset src_site)
{
  w_serial_ks_fm(kspf, src_site, NULL);
}
/*---------------------------------------------------------------------------*/

void w_serial_ks_fm_from_field(ks_prop_file *kspf, su3_vector *src_field)
{
  w_serial_ks_fm(kspf, (field_offset)(-1), src_field);
}
/*---------------------------------------------------------------------------*/

void w_serial_ks_fm_f(ks_prop_file *kspf)

/* this subroutine closes the file and frees associated structures */
{
  g_sync();
  if(this_node==0)
    {
      if(kspf->parallel == PARALLEL)
	printf("w_serial_ks_f: Attempting serial close on file opened in parallel \n");

      printf("Wrote KS prop file %s time stamp %s\n", kspf->filename,
	     (kspf->header)->time_stamp);

      if(kspf->fp != NULL) fclose(kspf->fp);
    }

  /* Free header and file structures */
  destroy_ks_fmprop_file_handle(kspf);

} /* w_serial_ks_fm_f */


/* !!!!!!!!!!!!!!!!!!!!!!! My code for reading a binary fm staggered propagator*/

int read_ks_fmprop_hdr(ks_prop_file *kspf, int parallel)
{
  /* parallel = 1 (TRUE) if all nodes are accessing the file */
  /*            0        for access from node 0 only */

  /* Returns byterevflag  = 0 or 1 */

  FILE *fp = NULL;
  ks_prop_header *ksph;
  int32type tmp;
  int32type elements_per_site, size_of_element;
  int j;
  int byterevflag = 0;
  char myname[] = "read_ks_fmprop_hdr";

  fp = kspf->fp;
  ksph = kspf->header;

  /* Read and verify magic number */

  if(psread_data(parallel, fp,&ksph->magic_number,sizeof(ksph->magic_number),
	      myname,"magic number")!=0)terminate(1);

  tmp = ksph->magic_number;
  
  if(ksph->magic_number == IO_UNI_MAGIC) 
    byterevflag=0;
  else 
    {
      byterevn((int32type *)&ksph->magic_number,1);
      if(ksph->magic_number == IO_UNI_MAGIC) 
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
	  ksph->magic_number = tmp;
	  
	  /* End of the road. */
	  printf("%s: Unrecognized magic number in prop file header.\n",
		 myname);
	  printf("Expected %x but read %x\n",
		 IO_UNI_MAGIC,tmp);
	  terminate(1);
	}
    }
  /**printf("byterevflag %d, parallel %d\n, magic number %0x\n", byterevflag, parallel, ksph->magic_number);  **/
  /* Read header, do byte reversal, 
     if necessary, and check consistency */

    /* Date and time stamp */

  if(psread_data(parallel,fp, &(ksph->gmtime_stamp), sizeof(ksph->gmtime_stamp),
		 myname,"time stamp")!=0) terminate(1);

  if(byterevflag)  byterevn(&(ksph->gmtime_stamp), 1);
  //printf("hey! time: %s\n", ksph->gmtime_stamp);
  //elements per site, size of element
  if(psread_data(parallel,fp, &size_of_element, sizeof(int32type),
	      myname,"size of element")!=0) terminate(1);
  if(byterevflag)  byterevn(&size_of_element, 1);
  if(size_of_element != sizeof(float)){printf("wrong size of element. read %d\n", size_of_element);terminate(1);}
  if(psread_data(parallel,fp, &elements_per_site, sizeof(int32type),
	      myname,"elements per site")!=0) terminate(1);
  if(byterevflag)  byterevn(&elements_per_site, 1);
  if(elements_per_site != 18){printf("wrong number of elements. read %d\n", elements_per_site);terminate(1);}
  /* Lattice dimensions */
  
  if(psread_byteorder(byterevflag,parallel,fp,(void *)(ksph->dims),sizeof(ksph->dims),
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
  

  /* Header byte length */

  ksph->header_bytes = sizeof(ksph->magic_number) + sizeof(ksph->dims) + 
    sizeof(ksph->gmtime_stamp) + sizeof(ksph->order)+ 2*sizeof(int32type);
  
  /* Data order */
  
  if(psread_byteorder(byterevflag,parallel,fp,&ksph->order,sizeof(ksph->order),
		   myname,"order parameter")!=0) terminate(1);
  if(byterevflag)  byterevn(&ksph->order, 1);
  
  return byterevflag;
  
} /* read_ks_fmprop_hdr */


ks_prop_file *r_serial_ks_fm_i(char *filename)
{
  /* Returns file descriptor for opened file */

  ks_prop_header *ksph;
  ks_prop_file *kspf;
  FILE *fp = NULL;
  int byterevflag;

  /* All nodes set up a propagator file and propagator header
     structure for reading */

  kspf = create_input_ks_fmprop_file_handle(filename);
  ksph = kspf->header;

  /* File opened for serial reading */
  kspf->parallel = 0;

  /* Node 0 alone opens a file and reads the header */

  if(this_node==0)
    {
      fp = fopen(filename, "rb");
      if(fp == NULL)
	{
	  printf("r_serial_ks_fm_i: Node %d can't open file %s, error %d\n",
		 this_node,filename,errno);fflush(stdout);terminate(1);
	}
      
/*      printf("Opened prop file %s for serial reading\n",filename); */
      
      kspf->fp = fp;

      byterevflag = read_ks_fmprop_hdr(kspf,SERIAL);

    }

  else kspf->fp = NULL;

  /* Broadcast the byterevflag from node 0 to all nodes */
      
  broadcast_bytes((char *)&byterevflag,sizeof(byterevflag));
  kspf->byterevflag = byterevflag;
  
  /* Node 0 broadcasts the header structure to all nodes */
  
  broadcast_bytes((char *)ksph,sizeof(ks_prop_header));

  return kspf;

}/* r_serial_ks_fm_i */



int r_serial_ks_fm(ks_prop_file *kspf, field_offset dest_site, 
		   su3_vector *dest_field)
{
  /* 0 is normal exit code
     1 for seek, read error, or missing data error */

  FILE *fp = NULL;
  ks_prop_header *ksph;
  char *filename;
  int byterevflag,a,b;

  off_t offset ;            /* File stream pointer */
  off_t coord_list_size;    /* Size of coordinate list in bytes */
  off_t head_size;          /* Size of header plus coordinate list */

  int rcv_rank, rcv_coords;
  int destnode;
  int i = 0,k,x,y,z,t;
  int status;
  int buf_length = 0, where_in_buf = 0;
  ks_prop_check test_kspc;
  u_int32type *val;
  int rank29 = 0,rank31 = 0;
  fsu3_vector *pbuf = NULL;
  su3_vector *dest;

  struct {
    fsu3_vector ksv[3];
    char pad[PAD_SEND_BUF];    /* Introduced because some switches
				  perform better if message lengths are longer */
  } msg;

  char myname[] = "r_serial_ks_fm";

  fp = kspf->fp;
  ksph = kspf->header;
  filename = kspf->filename;
  byterevflag = kspf->byterevflag;

  status = 0;
  if(this_node == 0)
    {
      if(kspf->parallel == PARALLEL)
	printf("%s: Attempting serial read from parallel file \n",myname);
    }
  broadcast_bytes((char *)&status,sizeof(int));
  if(status != 0) return status;

  status = 0;
  if(this_node == 0)
    {
      if(ksph->order == NATURAL_ORDER) coord_list_size = 0;
      else coord_list_size = sizeof(int32type)*volume;
      head_size = ksph->header_bytes + coord_list_size;
     
      offset = head_size;
      /**      printf("OFFSET %d\n", offset);**/
      pbuf = (fsu3_vector *)malloc(MAX_BUF_LENGTH*3*sizeof(fsu3_vector));
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
      /* Check it for compatibility!* ??????????????????????????????????????????*/

      /*status += sread_byteorder(byterevflag,fp,&kspf->check.sum29,
		      sizeof(kspf->check.sum29),myname,"check.sum29");
		      status += sread_byteorder(byterevflag,fp,&kspf->check.sum31,
		      sizeof(kspf->check.sum31),myname,"check.sum31");*/

      buf_length = 0;
      where_in_buf = 0;

    } /* end of if(this_node == 0)*/

  broadcast_bytes((char *)&status,sizeof(int));
  if(status != 0) return status;

  g_sync();

  /**  printf("   Start the reading from stag prop file, buf_length %d   \n", buf_length); **/
  /* Node 0 reads and deals out the values */
  status = 0;
  for(rcv_rank=0; rcv_rank<volume; rcv_rank++)
    {
      /* File is always in coordinate natural order, so receiving coordinate
         is given by rank */
      
      rcv_coords = rcv_rank;

      x = rcv_coords % nx;   rcv_coords /= nx;
      y = rcv_coords % ny;   rcv_coords /= ny;
      z = rcv_coords % nz;   rcv_coords /= nz;
      t = rcv_coords % nt;
     
      /* All nodes initialize timeslice checksums at the beginning of
	 a time slice */
      if(x == 0 && y == 0 && z == 0)
	{
	  test_kspc.sum31 = 0;
	  test_kspc.sum29 = 0;
	  /* counts 32-bit words mod 29 and mod 31 in order of appearance
	     on file */
	  /* Here all nodes see the same sequence because we read serially */
	  rank29 = 0;
	  rank31 = 0;
	}
  
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
	    a=(int)fread(pbuf,3*sizeof(fsu3_vector),buf_length,fp);
	    //if( (int)fread(pbuf,3*sizeof(fsu3_vector),buf_length,fp) 
	    if( a	!= buf_length)
	      {
		printf("InSIDE: buf length %d, max_buff_length %d, read elements %d \n",
		       buf_length, MAX_BUF_LENGTH, a );
		if(status == 0)
		  printf("%s: node %d propagator read error %d file %s\n",
			 myname, this_node, errno, filename); 
		fflush(stdout); 
		status = 1;
	      }
	    where_in_buf = 0;  /* reset counter */
	  }  /*** end of the buffer read ****/

	/* Save matrix in msg.ksv for further processing */
	for(a = 0; a < 3; a++)
	  msg.ksv[a] = pbuf[3*where_in_buf + a];
	
	if(destnode==0) {	/* Our value is already in msg */
	  i = node_index(x,y,z,t);
	}
	else {		        /* send to correct node */
	  send_field((char *)&msg, sizeof(msg), destnode);
	}
	where_in_buf++;
      }
      /* if(this_node==0) */



      /* The node which contains this site reads message */
      else {	/* for all nodes other than node 0 */
	if(this_node==destnode){
	  i = node_index(x,y,z,t);
	  /* Receive padded message in msg and copy to temporary location */
	  get_field((char *)&msg, sizeof(msg),0);
	}
      }

      /* The receiving node does the byte reversal and then checksum,
         if needed.  At this point msg contains the input vectors
	 and i points to the destination site structure */
      
      if(this_node==destnode)
	{
	  if(byterevflag==1)
	    byterevn((int32type *)&msg.ksv, 
		     3*sizeof(fsu3_vector)/sizeof(int32type));
	  /* Accumulate checksums */
	  for(k = 0, val = (u_int32type *)(&msg.ksv); 
	      k < 3*(int)sizeof(fsu3_vector)/(int)sizeof(int32type); 
	      k++, val++)
	    {
	      test_kspc.sum29 ^= (*val)<<rank29 | (*val)>>(32-rank29);
	      test_kspc.sum31 ^= (*val)<<rank31 | (*val)>>(32-rank31);
	      rank29++; if(rank29 >= 29)rank29 = 0;
	      rank31++; if(rank31 >= 31)rank31 = 0;
	    }
	  /* Now copy to dest, converting to generic precision if needed */

	  if(dest_site == (field_offset)(-1))
	    dest = dest_field + 3*i;
	  else
	    dest = (su3_vector *)F_PT(lattice+i,dest_site);

	  for(a=0; a<3; a++)for(b=0;b<3;b++){
	    dest[a].c[b].real = msg.ksv[a].c[b].real;
	    dest[a].c[b].imag = msg.ksv[a].c[b].imag;
	  }
	}
      else
	{
	  rank29 += 3*sizeof(fsu3_vector)/sizeof(int32type);
	  rank31 += 3*sizeof(fsu3_vector)/sizeof(int32type);
	  rank29 %= 29;
	  rank31 %= 31;
	}

      /* Accumulate and print checksum at the end of each time slice */
      if(x == nx - 1 && y == ny - 1 && z == nz - 1)
	{
	  /* Combine node checksum contributions with global exclusive or */
	  g_xor32(&test_kspc.sum29);
	  g_xor32(&test_kspc.sum31);

	  node0_printf("quark.t[%d].checksum  \"%0x %0x\"\n",t,
		       test_kspc.sum29, test_kspc.sum31);
	}

    } /* rcv_rank */

  broadcast_bytes((char *)&status,sizeof(int));
  if(status != 0) return status;

  if(this_node==0)
    {
      printf("Read KS prop serially from file %s\n", filename);
      
      /* Verify checksum */
      /* Checksums not implemented until version 5 */

      if(ksph->magic_number == IO_UNI_MAGIC)
	{
	  /*if(kspf->check.sum29 != test_kspc.sum29 ||
	    kspf->check.sum31 != test_kspc.sum31)
	    {
	    printf("%s: Checksum violation file %s\n",
	    myname, kspf->filename);*/

	  printf("Computed checksum %x %x.  Read %x %x.\n",
		 test_kspc.sum29, test_kspc.sum31,
		 kspf->check.sum29, kspf->check.sum31);
	}
/*	  else
	    printf("Checksums %x %x OK for file %s\n",
		   kspf->check.sum29,kspf->check.sum31,kspf->filename); */
      
      fflush(stdout);
      free(pbuf);
    }

  return 0;

} /* r_serial_ks_fm */

int r_serial_ks_fm_to_site(ks_prop_file *kspf, field_offset dest_site)
{
  int status = r_serial_ks_fm(kspf, dest_site, NULL);
  kspf->info = read_info_file(kspf->filename);
  return status;
}

int r_serial_ks_fm_to_field(ks_prop_file *kspf, su3_vector *dest_field)
{
  int status = r_serial_ks_fm(kspf, (field_offset)(-1), dest_field);
  kspf->info = read_info_file(kspf->filename);
  return status;
}

void r_serial_ks_fm_f(ks_prop_file *kspf)

/* Close the file and free associated structures */
{
  g_sync();
  if(this_node==0)
    {
      if(kspf->parallel == PARALLEL)
	printf("r_serial_ks_fm_f: Attempting serial close on parallel file \n");
      
      if(kspf->fp != NULL) fclose(kspf->fp);
/*      printf("Closed prop file %s\n",kspf->filename);*/
      fflush(stdout);
    }
  
  destroy_ks_fmprop_file_handle(kspf);
  
} /* r_serial_ks_f */
