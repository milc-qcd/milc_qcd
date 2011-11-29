/**************************** io_source_ks_fm.c ****************************/
/* MIMD version 7 */
/*
 *  Read a Fermilab format source file for staggered fermions
 *
 */


#include "generic_ks_includes.h"
#include <string.h>
#include <assert.h>
#include <errno.h>
#define MAX_BUF_LENGTH 4096

#define IO_UNI_MAGIC 0x71626434 /* "qbd4" */

#ifndef PAD_SEND_BUF
#define PAD_SEND_BUF 8
#endif
#define NATURAL_ORDER 0

static int 
read_ks_fm_source_hdr(ks_fm_source_file *kssf)
{
  ks_source_header *kssh = kssf->header;
  int *dims = kssh->dims;

  int32type   t_stamp;
  int32type   size_of_element;
  int32type   tmp, elements_per_site; 
  int32type   order; 
  int byterevflag = 0;
  char myname[] = "read_ks_fm_source_hdr";

  if( sizeof(float) != sizeof(int32type)) {
    printf("%s: Can't byte reverse\n", kssf->filename);
    printf("requires size of int32type(%d) = size of float(%d)\n",
	   (int)sizeof(int32type),(int)sizeof(float));
    terminate(1);
  }

  if(sread_data(kssf->fp,&kssh->magic_number,sizeof(int32type),
		myname,"magic_no")) terminate(1);

  tmp = kssh->magic_number;
  if(kssh->magic_number == IO_UNI_MAGIC) byterevflag = 0;
  else 
    {
      byterevn((int32type *)&kssh->magic_number,1);
      if(kssh->magic_number == IO_UNI_MAGIC) 
	{
	  byterevflag = 1; 
	  /** printf("Reading with byte reversal\n"); **/
	}
      else
	{
	  /* Restore magic number as originally read */
	  kssh->magic_number = tmp;
	  
	  /* End of the road. */
	  printf("%s: Unrecognized magic number in source file header.\n",
		 kssf->filename);
	  printf("Expected %x but read %x\n", IO_UNI_MAGIC,tmp);
	  terminate(1);
	}
    };

  if(sread_byteorder(byterevflag,kssf->fp,&t_stamp,
	     sizeof(t_stamp),myname,"t_stamp"))terminate(1);
  if(sread_byteorder(byterevflag,kssf->fp,&size_of_element,
	     sizeof(int32type),myname,"size_of_element"))terminate(1);
  if(sread_byteorder(byterevflag,kssf->fp,&elements_per_site,
	     sizeof(int32type),myname,"elements_per_site"))terminate(1);
  if(sread_byteorder(byterevflag,kssf->fp,dims,
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
      elements_per_site != 6 /* su3_vector */)
    {	
      printf(" File %s is not a vector field",  kssf->filename);
      printf(" got size_of_element %d and elements_per_site %d\n",
	     size_of_element, elements_per_site);
      terminate(1);
    }
  
  /* The site order parameter is ignored */
  
  if(sread_byteorder(byterevflag,kssf->fp,&order,sizeof(int32type),
		     myname,"order parameter"))terminate(1);
  
  return byterevflag;
}
 
static ks_fm_source_file *
setup_input_ks_fm_source_file(char *filename)
{
  ks_fm_source_file *kssf;
  ks_source_header *ksph;

  /* Allocate space for the file structure */

  kssf = (ks_fm_source_file *)malloc(sizeof(ks_fm_source_file));
  if(kssf == NULL)
    {
      printf("setup_input_ks_fm_source_file: Can't malloc kssf\n");
      terminate(1);
    }

  strncpy(kssf->filename,filename,MAXFILENAME);
  kssf->filename[MAXFILENAME-1] = '\0';

  /* Allocate space for the header */

  /* Make sure compilation gave us a 32 bit integer type */
  assert(sizeof(int32type) == 4);

  ksph = (ks_source_header *)malloc(sizeof(ks_source_header));
  if(ksph == NULL)
    {
      printf("setup_input_ks_fm_source_file: Can't malloc ksph\n");
      terminate(1);
    }

  kssf->header = ksph;

  return kssf;
} /* setup_input_ks_fm_source_file */

/*--------------------------------------------------------------------*/
ks_fm_source_file *
r_source_ks_fm_i(char *filename)
{
  /* Returns file descriptor for opened file */

  ks_fm_source_file   *kssf;
  ks_source_header *kssh;
  int byterevflag;

  kssf = setup_input_ks_fm_source_file(filename);
  kssh = kssf->header;

  if(this_node==0){
    kssf->fp = g_open(filename,"rb");
    if(kssf->fp==NULL){
      printf("Can't open source file %s, error %d\n",filename,errno);
      terminate(1);
    }
    byterevflag = read_ks_fm_source_hdr(kssf);

  }
  else kssf->fp = NULL;

  /* Broadcast the byterevflag from node 0 to all nodes */
  broadcast_bytes((char *)&byterevflag,sizeof(byterevflag));
  kssf->byterevflag = byterevflag;
  
  /* Node 0 broadcasts the header structure to all nodes */
  broadcast_bytes((char *)kssh,sizeof(ks_source_header));

  return kssf;

} /* r_source_ks_fm_i */

/*--------------------------------------------------------------------*/
void 
r_source_ks_fm(ks_fm_source_file *kssf, 
	      su3_vector *dest_field,
	      int x0, int y0, int z0, int t0)
{
  int rcv_rank, rcv_coords, status;
  int destnode;
  int x,y,z,t,i=0, byterevflag,c,a;
  int tmin, tmax;

  struct {
    fsu3_vector q;
    char pad[PAD_SEND_BUF];    /* Introduced because some switches
				  perform better if message lengths
				  are longer */
  } ksmsg;

  int buf_length=0, where_in_buf=0;
  fsu3_vector *vbuff=NULL;
  fsu3_vector vfix;
  su3_vector *v;
  int vol3 = nx*ny*nz;

  byterevflag = kssf->byterevflag;

  if(this_node == 0)
    {
      vbuff = (fsu3_vector *)malloc(MAX_BUF_LENGTH*sizeof(fsu3_vector));

      if(vbuff == NULL)
	{
	  printf("Node %d can't malloc vbuff\n",this_node);
	  fflush(stdout);
	  terminate(1);
	}
      
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
      
      /* Build in the requested translation in converting
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
	    a=(int)g_read(vbuff,sizeof(fsu3_vector),buf_length,kssf->fp);
	    
	    if( a != buf_length)
	      {
		if(status == 0)
		  printf(" node %d source file read error %d file %s\n",
			 this_node, errno, kssf->filename); 
		fflush(stdout); 
		status = 1;
	      }
	    where_in_buf = 0;  /* reset counter */
	  }  /*** end of the buffer read ****/
	
	/* Save data in msg.q for further processing */
	ksmsg.q = vbuff[where_in_buf];
      }
      
      /* Loop either over all time slices or just one time slice */
      for(t = tmin; t <= tmax; t++){
	
	destnode=node_number(x,y,z,t);
	
	if(this_node == 0){
	  /* node 0 doesn't send to itself */
	  if(destnode != 0){
	    /* send to correct node */
	    send_field((char *)&ksmsg, sizeof(ksmsg), destnode);
	  }
	} /* if(this_node==0) */
	else {	/* for all nodes other than node 0 */
	  if(this_node==destnode){
	    get_field((char *)&ksmsg, sizeof(ksmsg),0);
	  }
	}
	
	/* The receiving node does the byte reversal.  At this point msg
	   contains the input vectors and i points to the destination
	   site structure */
	
	if(this_node==destnode)
	  {
	    vfix = ksmsg.q;
	    if(byterevflag==1)
	      byterevn((int32type *)&vfix, 
		       sizeof(fsu3_vector)/sizeof(int32type));
	    
	    /* Now copy the site data into the destination converting
	       to generic precision if needed */
	    
	    i = node_index(x,y,z,t);
	    
	    v = dest_field + i;
	    
	    for(c = 0; c < 3; c++){
	      v->c[c].real = vfix.c[c].real;
	      v->c[c].imag = vfix.c[c].imag;
	    }
	  } /* if this_node == destnode */
      } /* t */
      
      where_in_buf++;
      
    } /* rcv_rank */
  
  
  if(vbuff != NULL)free(vbuff); vbuff = NULL;
}

/*----------------------------------------------------------------------*/
/* Close the file and free associated structures */

void 
r_source_ks_fm_f(ks_fm_source_file *kssf)
{
  g_sync();
  if(this_node==0)
    {
      if(kssf->fp != NULL)g_close(kssf->fp);
      fflush(stdout);
    }
  free(kssf->header);
  free(kssf);
  
} /* r_source_ks_fm_f */

