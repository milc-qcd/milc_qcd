/************************** io_source_cmplx_fm.c ****************************/
/* MIMD version 7 */
/*
 *  Load the Fermilab format file for a 3D complex field.
 *
 */


#include "generic_includes.h"
#include <string.h>
#include <assert.h>
#include <errno.h>
#define MAX_BUF_LENGTH 4096
#define ALL_T_SLICES -1

#define IO_UNI_MAGIC 0x71626434 /* "qbd4" */

typedef struct {
  int32type magic_number;
  int32type gmtime;
  int32type size_of_element;
  int32type elements_per_site;
  int32type dims[4];
  int32type site_order;
} cmplx_source_header;

typedef struct {
  cmplx_source_header *header;
  FILE *fp;
  char filename[MAXFILENAME];
  int byterevflag;
} cmplx_source_file;

#ifndef PAD_SEND_BUF
#define PAD_SEND_BUF 8
#endif
#define NATURAL_ORDER 0

static int 
read_cmplx_fm_source_hdr(cmplx_source_file *csf)
{
  cmplx_source_header *csh = csf->header;
  int *dims = csh->dims;

  int32type   t_stamp;
  int32type   size_of_element;
  int32type   tmp, elements_per_site; 
  int32type   order; 
  int byterevflag = 0;
  char myname[] = "read_cmplx_fm_source_hdr";

  if( sizeof(float) != sizeof(int32type)) {
    printf("%s: Can't byte reverse\n", csf->filename);
    printf("requires size of int32type(%d) = size of float(%d)\n",
	   (int)sizeof(int32type),(int)sizeof(float));
    terminate(1);
  }

  if(sread_data(csf->fp,&csh->magic_number,sizeof(int32type),
		myname,"magic_no")) terminate(1);

  tmp = csh->magic_number;
  if(csh->magic_number == IO_UNI_MAGIC) byterevflag = 0;
  else 
    {
      byterevn((int32type *)&csh->magic_number,1);
      if(csh->magic_number == IO_UNI_MAGIC) 
	{
	  byterevflag = 1; 
	  /** printf("Reading with byte reversal\n"); **/
	}
      else
	{
	  /* Restore magic number as originally read */
	  csh->magic_number = tmp;
	  
	  /* End of the road. */
	  printf("%s: Unrecognized magic number in source file header.\n",
		 csf->filename);
	  printf("Expected %x but read %x\n", IO_UNI_MAGIC,tmp);
	  terminate(1);
	}
    };

  if(sread_byteorder(byterevflag,csf->fp,&t_stamp,
	     sizeof(t_stamp),myname,"t_stamp"))terminate(1);
  if(sread_byteorder(byterevflag,csf->fp,&size_of_element,
	     sizeof(int32type),myname,"size_of_element"))terminate(1);
  if(sread_byteorder(byterevflag,csf->fp,&elements_per_site,
	     sizeof(int32type),myname,"elements_per_site"))terminate(1);
  if(sread_byteorder(byterevflag,csf->fp,dims,
	     4*sizeof(int32type),myname,"dimensions"))terminate(1);

  /* Consistency checks */
  
  if ( nx != dims[0] ||
       ny != dims[1] ||
       nz != dims[2]  ||
       1  != dims[3] ){
    printf("Source field dimensions %d %d %d %d are incompatible with this lattice %d %d %d %d\n", dims[0], dims[1], dims[2], dims[3],  nx, ny, nz, nt);
    terminate(1);
  }

  if( size_of_element != sizeof(float) ||
      elements_per_site != 2 /* complex field */)
    {	
      printf(" File %s is not a complex field",  csf->filename);
      printf(" got size_of_element %d and elements_per_site %d\n",
	     size_of_element, elements_per_site);
      terminate(1);
    }
  
  /* The site order parameter is ignored */
  
  if(sread_byteorder(byterevflag,csf->fp,&order,sizeof(int32type),
		     myname,"order parameter"))terminate(1);
  
  return byterevflag;
}
 
static cmplx_source_file *
setup_input_cmplx_source_file(char *filename)
{
  cmplx_source_file *csf;
  cmplx_source_header *cph;

  /* Allocate space for the file structure */

  csf = (cmplx_source_file *)malloc(sizeof(cmplx_source_file));
  if(csf == NULL)
    {
      printf("setup_input_cmplx_source_file: Can't malloc csf\n");
      terminate(1);
    }

  strncpy(csf->filename,filename,MAXFILENAME);
  csf->filename[MAXFILENAME-1] = '\0';

  /* Allocate space for the header */

  /* Make sure compilation gave us a 32 bit integer type */
  assert(sizeof(int32type) == 4);

  cph = (cmplx_source_header *)malloc(sizeof(cmplx_source_header));
  if(cph == NULL)
    {
      printf("setup_input_cmplx_source_file: Can't malloc cph\n");
      terminate(1);
    }

  csf->header = cph;

  return csf;
} /* setup_input_cmplx_source_file */

static cmplx_source_file *
r_source_cmplx_fm_i(char *filename)
{
  /* Returns file descriptor for opened file */

  cmplx_source_file   *csf;
  cmplx_source_header *csh;
  int byterevflag;

  csf = setup_input_cmplx_source_file(filename);
  csh = csf->header;

  if(this_node==0){
    csf->fp = g_open(filename,"rb");
    if(csf->fp==NULL){
      printf("Can't open source file %s, error %d\n",filename,errno);
      terminate(1);
    }
    byterevflag = read_cmplx_fm_source_hdr(csf);

  }
  else csf->fp = NULL;

  /* Broadcast the byterevflag from node 0 to all nodes */
  broadcast_bytes((char *)&byterevflag,sizeof(byterevflag));
  csf->byterevflag = byterevflag;
  
  /* Node 0 broadcasts the header structure to all nodes */
  broadcast_bytes((char *)csh,sizeof(cmplx_source_header));

  return csf;

} /* r_source_cmplx_fm_i */


static void 
r_source_cmplx_fm(cmplx_source_file *csf, field_offset dest_site, 
		  complex *dest_field, int stride,
		  int x0, int y0, int z0, int t0)
{
  int rcv_rank, rcv_coords, status;
  int destnode;
  int x,y,z,t,i=0, byterevflag,a;
  int tmin, tmax;

  struct {
    fcomplex q;
    char pad[PAD_SEND_BUF];    /* Introduced because some switches
				  perform better if message lengths
				  are longer */
  } cmsg;

  int buf_length=0, where_in_buf=0;
  fcomplex *cbuff=NULL;
  complex *c;
  fcomplex cfix;
  int vol3 = nx*ny*nz;

  byterevflag = csf->byterevflag;

  if(this_node == 0)
    {
      cbuff = (fcomplex *)malloc(MAX_BUF_LENGTH*sizeof(fcomplex));
      buf_length = 0;
      where_in_buf = 0;

    } /* end of if(this_node == 0)*/

  g_sync();

  /* If requested, we replicate the source function on ALL time slices */
  /* Otherwise we read only to time slice t0 */

  if(t0 == ALL_T_SLICES){ tmin = 0; tmax = nt-1; }
  else      { tmin = t0; tmax = t0; }

  /* Node 0 reads and deals out the values */
  status = 0;

  /* Iterate only over timeslice t0 */
  for(rcv_rank=0; rcv_rank<vol3; rcv_rank++)
    {
      /* We do only natural (lexicographic) order here */
      rcv_coords = rcv_rank;
      
      /* Include the requested translation in the conversion from
	 lexicographic to Cartesian coordinates */
      x = (rcv_coords + x0) % nx;   rcv_coords /= nx;
      y = (rcv_coords + y0) % ny;   rcv_coords /= ny;
      z = (rcv_coords + z0) % nz;
      
      if(this_node==0){
	/* Node 0 fills its buffer, if necessary */
	if(where_in_buf == buf_length)
	  {  /* get new buffer */
	    /* new buffer length  = remaining sites, but never bigger 
	       than MAX_BUF_LENGTH */
	    buf_length = vol3 - rcv_rank;
	    if(buf_length > MAX_BUF_LENGTH) buf_length = MAX_BUF_LENGTH;
	    /* then do read */
	    a=(int)g_read(cbuff,sizeof(fcomplex),buf_length,csf->fp);
	    
	    if( a != buf_length)
	      {
		if(status == 0)
		  printf(" node %d source file read error %d file %s\n",
			 this_node, errno, csf->filename); 
		fflush(stdout); 
		status = 1;
	      }
	    where_in_buf = 0;  /* reset counter */
	  }  /*** end of the buffer read ****/
	
	/* Save data in msg.q for further processing */
	cmsg.q = cbuff[where_in_buf];
      }
      
      /* Loop either over all time slices or just one time slice */
      for(t = tmin; t <= tmax; t++){
	
	destnode=node_number(x,y,z,t);
	
	if(this_node == 0){
	  /* node 0 doesn't send to itself */
	  if(destnode != 0){
	    /* send to correct node */
	    send_field((char *)&cmsg, sizeof(cmsg), destnode);
	  }
	} /* if(this_node==0) */
	else {	/* for all nodes other than node 0 */
	  if(this_node==destnode){
	    get_field((char *)&cmsg, sizeof(cmsg),0);
	  }
	}
	
	/* The receiving node does the byte reversal.  At this point msg
	   contains the input vectors and i points to the destination
	   site structure */
	
	if(this_node==destnode)
	{
	    /* Byte reverse a copy, since we may need to reuse cmsg.q */
	    cfix = cmsg.q;
	    if(byterevflag==1){
		byterevn((int32type *)&cfix, 
			 sizeof(fcomplex)/sizeof(int32type));
	    }
	    
	    /* Now copy the site data into the destination converting
	       to generic precision if needed */
	    
	    i = node_index(x,y,z,t);
	    
	    if(dest_site == (field_offset)(-1))
		c = dest_field + i*stride;
	    else
		c = (complex *)F_PT(&lattice[i],dest_site);
	    
	    c->real = cfix.real;
	    c->imag = cfix.imag;
	    //printf("C %d %d %d %g %g\n", 
	    //	       x, y, z, cmsg.q.real, cmsg.q.imag);
	} /* if this_node == destnode */
      } /* t */
      
      where_in_buf++;
      
    }  /* rcv_rank, irecord */
  
  if(cbuff != NULL)free(cbuff); cbuff = NULL;
}

/*----------------------------------------------------------------------*/

static void 
r_source_cmplx_fm_f(cmplx_source_file *csf)

/* Close the file and free associated structures */
{
  g_sync();
  if(this_node==0)
    {
      if(csf->fp != NULL)g_close(csf->fp);
      fflush(stdout);
    }
  free(csf->header);
  free(csf);
  
} /* r_source_cmplx_fm_f */

/*--------------------------------------------------------------------*/
static void 
r_source_cmplx_fm_generic(char *filename, 
			  field_offset dest_site,
			  complex *dest_field, int stride,
			  int x0, int y0, int z0, int t0)
{
  cmplx_source_file *csf;
  csf = r_source_cmplx_fm_i(filename);
  r_source_cmplx_fm(csf, dest_site, dest_field, stride, x0, y0, z0, t0);
  r_source_cmplx_fm_f(csf);
  
}

/*--------------------------------------------------------------------*/
void r_source_cmplx_fm_to_site(char *filename, field_offset dest_site,
			       int x0, int y0, int z0, int t0)
{
  r_source_cmplx_fm_generic(filename, dest_site, (complex *)NULL, 0, x0, y0, z0, t0);
}

/*--------------------------------------------------------------------*/
void r_source_cmplx_fm_to_field(char *filename, complex *dest_field, int stride,
				int x0, int y0, int z0, int t0)
{
  r_source_cmplx_fm_generic(filename, (field_offset)(-1), dest_field, stride,
			    x0, y0, z0, t0);
}

