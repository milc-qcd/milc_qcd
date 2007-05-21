/**************************** io_source_w_fm.c ****************************/
/* MIMD version 7 */
/*
 *  Load the Fermilab format file for Wilson propagator sources
 *  Supports the "wavefunction" format - a 3D complex field
 *  and the "timeslice" format - one 3D Dirac field for each
 *  of 12 source colors and spins.
 *
 */


#include "generic_wilson_includes.h"
#include <string.h>
#include <assert.h>
#include <errno.h>
#define MAX_BUF_LENGTH 4096

#define IO_UNI_MAGIC 0x71626434 /* "qbd4" */

typedef struct {
  int32type magic_number;
  int32type gmtime;
  int32type size_of_element;
  int32type elements_per_site;
  int32type dims[4];
  int32type site_order;
} w_source_header;

typedef struct {
  w_source_header *header;
  FILE *fp;
  char filename[MAXFILENAME];
  int byterevflag;
  int type;
  wilson_propagator *source;
} w_source_file;

#ifndef PAD_SEND_BUF
#define PAD_SEND_BUF 8
#endif
#define NATURAL_ORDER 0

static int 
read_w_fm_source_hdr(w_source_file *wsf, int source_type)
{
  w_source_header *wsh = wsf->header;
  int *dims = wsh->dims;

  int32type   t_stamp;
  int32type   size_of_element;
  int32type   tmp, elements_per_site; 
  int32type   order; 
  int byterevflag = 0;
  char myname[] = "read_w_fm_source_hdr";

  if( sizeof(float) != sizeof(int32type)) {
    printf("%s: Can't byte reverse\n", wsf->filename);
    printf("requires size of int32type(%d) = size of float(%d)\n",
	   (int)sizeof(int32type),(int)sizeof(float));
    terminate(1);
  }

  if(sread_data(wsf->fp,&wsh->magic_number,sizeof(int32type),
		myname,"magic_no")) terminate(1);

  tmp = wsh->magic_number;
  if(wsh->magic_number == IO_UNI_MAGIC) byterevflag = 0;
  else 
    {
      byterevn((int32type *)&wsh->magic_number,1);
      if(wsh->magic_number == IO_UNI_MAGIC) 
	{
	  byterevflag = 1; 
	  /** printf("Reading with byte reversal\n"); **/
	}
      else
	{
	  /* Restore magic number as originally read */
	  wsh->magic_number = tmp;
	  
	  /* End of the road. */
	  printf("%s: Unrecognized magic number in source file header.\n",
		 wsf->filename);
	  printf("Expected %x but read %x\n", IO_UNI_MAGIC,tmp);
	  terminate(1);
	}
    };

  if(sread_byteorder(byterevflag,wsf->fp,&t_stamp,
	     sizeof(t_stamp),myname,"t_stamp"))terminate(1);
  if(sread_byteorder(byterevflag,wsf->fp,&size_of_element,
	     sizeof(int32type),myname,"size_of_element"))terminate(1);
  if(sread_byteorder(byterevflag,wsf->fp,&elements_per_site,
	     sizeof(int32type),myname,"elements_per_site"))terminate(1);
  if(sread_byteorder(byterevflag,wsf->fp,dims,
	     4*sizeof(int32type),myname,"dimensions"))terminate(1);

  /* Consistency checks */
  
  if ( nx != dims[0] ||
       ny != dims[1] ||
       nz != dims[2]  ||
       1  != dims[3] ){
    printf("Source field dimensions %d %d %d %d are incompatible with this lattice %d %d %d %d\n", dims[0], dims[1], dims[2], dims[3],  nx, ny, nz, nt);
    terminate(1);
  }

  if( source_type == COMPLEX_FIELD ){
    if( size_of_element != sizeof(float) ||
	elements_per_site != 2 /* complex field */)
      {	
	printf(" File %s is not a complex field",  wsf->filename);
	printf(" got size_of_element %d and elements_per_site %d\n",
	       size_of_element, elements_per_site);
	terminate(1);
      }
  }
  else if( source_type == DIRAC_FIELD ){
    if( size_of_element != sizeof(float) ||
	elements_per_site != 24 /* wilson_vector */)
      {	
	printf(" File %s is not a Dirac field",  wsf->filename);
	printf(" got size_of_element %d and elements_per_site %d\n",
	       size_of_element, elements_per_site);
	terminate(1);
      }
  }
  else {
    printf("Unknown source type %d\n",source_type);
    terminate(1);
  }
  
  /* The site order parameter is ignored */
  
  if(sread_byteorder(byterevflag,wsf->fp,&order,sizeof(int32type),
		     myname,"order parameter"))terminate(1);
  
  return byterevflag;
}
 
static w_source_file *
setup_input_w_source_file(char *filename, int source_type)
{
  w_source_file *wsf;
  w_source_header *wph;

  /* Allocate space for the file structure */

  wsf = (w_source_file *)malloc(sizeof(w_source_file));
  if(wsf == NULL)
    {
      printf("setup_input_w_source_file: Can't malloc wsf\n");
      terminate(1);
    }

  strncpy(wsf->filename,filename,MAXFILENAME);
  wsf->filename[MAXFILENAME-1] = '\0';
  wsf->type = source_type;

  /* Allocate space for the header */

  /* Make sure compilation gave us a 32 bit integer type */
  assert(sizeof(int32type) == 4);

  wph = (w_source_header *)malloc(sizeof(w_source_header));
  if(wph == NULL)
    {
      printf("setup_input_w_source_file: Can't malloc wph\n");
      terminate(1);
    }

  wsf->header = wph;

  return wsf;
} /* setup_input_w_source_file */

static w_source_file *
r_source_w_fm_i(char *filename, int source_type)
{
  /* Returns file descriptor for opened file */

  w_source_file   *wsf;
  w_source_header *wsh;
  int byterevflag;

  wsf = setup_input_w_source_file(filename, source_type);
  wsh = wsf->header;

  if(this_node==0){
    wsf->fp = g_open(filename,"rb");
    if(wsf->fp==NULL){
      printf("Can't open source file %s, error %d\n",filename,errno);
      terminate(1);
    }
    byterevflag = read_w_fm_source_hdr(wsf, source_type);

  }
  else wsf->fp = NULL;

  /* Make room for the whole Dirac source */
  if(source_type == DIRAC_FIELD){
    wsf->source = 
      (wilson_propagator *)malloc(sizeof(wilson_propagator)*sites_on_node);
    if(wsf->source == NULL){
      printf("r_source_w_fm_i: No room for Dirac source\n");
      terminate(1);
    }
    memset(wsf->source, '\0', sizeof(wilson_propagator)*sites_on_node);
  }
  else
    wsf->source = NULL;

  /* Broadcast the byterevflag from node 0 to all nodes */
  broadcast_bytes((char *)&byterevflag,sizeof(byterevflag));
  wsf->byterevflag = byterevflag;
  
  /* Node 0 broadcasts the header structure to all nodes */
  broadcast_bytes((char *)wsh,sizeof(w_source_header));

  return wsf;

} /* r_source_w_fm_i */


static void 
r_source_w_fm(w_source_file *wsf, 
	      field_offset dest_site, 
	      wilson_vector *dest_field,
	      int spin, int color, int t0)
{
  int rcv_rank, rcv_coords, status;
  int destnode;
  int x,y,z,i=0, byterevflag, c0,s0,c1,s1,a;
  site *s;

  struct {
    fwilson_vector q;
    char pad[PAD_SEND_BUF];    /* Introduced because some switches
				  perform better if message lengths
				  are longer */
  } wmsg;

  struct {
    fcomplex q;
    char pad[PAD_SEND_BUF];    /* Introduced because some switches
				  perform better if message lengths
				  are longer */
  } cmsg;

  int buf_length=0, where_in_buf=0;
  int irecord, nrecord;
  fwilson_vector *wbuff=NULL;
  fcomplex *cbuff=NULL;
  wilson_vector *wv;
  wilson_propagator *wp;
  int vol3 = nx*ny*nz;
  int source_type;

  byterevflag = wsf->byterevflag;
  source_type = wsf->type;

  if(this_node == 0)
    {
      if(source_type == COMPLEX_FIELD)
	cbuff = (fcomplex *)malloc(MAX_BUF_LENGTH*sizeof(fcomplex));
      else if(source_type == DIRAC_FIELD)
	wbuff = (fwilson_vector *)malloc(MAX_BUF_LENGTH*sizeof(fwilson_vector));
      else {
	printf("r_source_w_fm: Unknown source type %d\n", source_type);
	fflush(stdout);
	terminate(1);
      }
      if(wbuff == NULL && cbuff == NULL)
	{
	  printf("Node %d can't malloc wbuff or cbuff\n",this_node);
	  fflush(stdout);
	  terminate(1);
	}
      
      buf_length = 0;
      where_in_buf = 0;

    } /* end of if(this_node == 0)*/

  if(source_type == COMPLEX_FIELD)
    nrecord = 1;
  else
    nrecord = 12;

  g_sync();

  /* Node 0 reads and deals out the values */
  status = 0;
  /* Iterate only over timeslice t0 */
  for(irecord = 0; irecord < nrecord; irecord++)
    for(rcv_rank=vol3*t0; rcv_rank<vol3*(t0+1); rcv_rank++)
      {
	/* We do only natural order here */
	rcv_coords = rcv_rank;
	
	x = rcv_coords % nx;   rcv_coords /= nx;
	y = rcv_coords % ny;   rcv_coords /= ny;
	z = rcv_coords % nz;

	destnode=node_number(x,y,z,t0);
	
	if(this_node==0){
	  /* Node 0 fills its buffer, if necessary */
	  if(where_in_buf == buf_length)
	    {  /* get new buffer */
	      /* new buffer length  = remaining sites, but never bigger 
		 than MAX_BUF_LENGTH */
	      buf_length = volume - rcv_rank;
	      if(buf_length > MAX_BUF_LENGTH) buf_length = MAX_BUF_LENGTH;
	      /* then do read */
	      if(source_type == COMPLEX_FIELD)
		a=(int)g_read(cbuff,sizeof(fcomplex),buf_length,wsf->fp);
	      else
		a=(int)g_read(wbuff,sizeof(fwilson_vector),buf_length,wsf->fp);
	      
	      if( a  != buf_length)
		{
		  
		  if(status == 0)
		    printf(" node %d source file read error %d file %s\n",
			   this_node, errno, wsf->filename); 
		  fflush(stdout); 
		  status = 1;
		}
	      where_in_buf = 0;  /* reset counter */
	    }  /*** end of the buffer read ****/
	  
	  /* Save data in msg.q for further processing */
	  if(source_type == COMPLEX_FIELD)
	    cmsg.q = cbuff[where_in_buf];
	  else
	    wmsg.q = wbuff[where_in_buf];
	  
	  if(destnode==0){	/* just copy su3_vector */
	    i = node_index(x,y,z,t0);
	  }
	  else {		        /* send to correct node */
	    if(source_type == COMPLEX_FIELD)
	      send_field((char *)&cmsg, sizeof(cmsg), destnode);
	    else
	      send_field((char *)&wmsg, sizeof(wmsg), destnode);
	  }
	  where_in_buf++;
	}
	/* if(this_node==0) */
	else {	/* for all nodes other than node 0 */
	  if(this_node==destnode){
	    if(source_type == COMPLEX_FIELD)
	      get_field((char *)&cmsg, sizeof(cmsg),0);
	    else
	      get_field((char *)&wmsg, sizeof(wmsg),0);
	  }
	}
	
	/* The receiving node does the byte reversal.  At this point msg
	   contains the input vectors and i points to the destination
	   site structure */
	
	if(this_node==destnode)
	  {
	    if(byterevflag==1){
	      if(source_type == COMPLEX_FIELD)
		byterevn((int32type *)&cmsg.q, 
			 sizeof(fcomplex)/sizeof(int32type));
	      else
		byterevn((int32type *)&wmsg.q, 
			 sizeof(fwilson_vector)/sizeof(int32type));
	    }
	    
	    /* Now copy the site data into the destination converting
	       to generic precision if needed */
	    
	    i = node_index(x,y,z,t0);

	    /* For complex source, we copy directly to the user
	       destination, filling only the designated spin and color
	       with the field and zero the rest */
	    
	    if(source_type == COMPLEX_FIELD){

	      if(dest_site == (field_offset)(-1))
		wv = dest_field + i;
	      else
		wv = (wilson_vector *)F_PT(&lattice[i],dest_site);

	      clear_wvec(wv);
	      wv->d[spin].c[color].real = cmsg.q.real;
	      wv->d[spin].c[color].imag = cmsg.q.imag;
	    }
	    
	    /* For Dirac field source, we copy first to our buffer,
	       filling all components */
	    else{
	      wp = wsf->source + i;
	      for(s1=0;s1<4;s1++)for(c1=0;c1<3;c1++)
		{
		  c0 = irecord % 3;
		  s0 = irecord / 3;
		  wp->c[c0].d[s0].d[s1].c[c1].real = wmsg.q.d[s1].c[c1].real;
		  wp->c[c0].d[s0].d[s1].c[c1].imag = wmsg.q.d[s1].c[c1].imag;
		}
	    }
	  }
      }  /* rcv_rank, irecord */

  if(cbuff != NULL)free(cbuff); cbuff = NULL;
  if(wbuff != NULL)free(wbuff); wbuff = NULL;

  /* Finally for the Dirac source type, convert spin basis from FNAL
     to MILC and then copy the requested MILC spin, color to the user
     Wilson vector */
  if(source_type == DIRAC_FIELD){
    convert_wprop_fnal_to_milc_field(wsf->source);
    /* Copy converted source to source field for this spin and color */
    FORALLSITES(i,s){

      if(dest_site == (field_offset)(-1))
	wv = dest_field + i;
      else
	wv = (wilson_vector *)F_PT(&lattice[i],dest_site);

      *wv = wsf->source[i].c[color].d[spin];
    }
  }
}

/*----------------------------------------------------------------------*/

static void 
r_source_w_fm_f(w_source_file *wsf)

/* Close the file and free associated structures */
{
  g_sync();
  if(this_node==0)
    {
      if(wsf->fp != NULL)g_close(wsf->fp);
      fflush(stdout);
    }
  if(wsf->source != NULL)
    free(wsf->source);
  free(wsf->header);
  free(wsf);
  
} /* r_source_w_fm_f */

/*--------------------------------------------------------------------*/
static void 
r_source_w_fm_generic(char *filename, 
			   field_offset dest_site,
			   wilson_vector *dest_field, 
			   int spin, int color, int t0, int source_type)
{
  w_source_file *wsf;
  wsf = r_source_w_fm_i(filename, source_type);
  r_source_w_fm(wsf, dest_site, dest_field, spin, color, t0);
  r_source_w_fm_f(wsf);
  
}

/*--------------------------------------------------------------------*/
void r_source_w_fm_to_site(char *filename, field_offset dest_site,
			   int spin, int color, int t0, int source_type)
{
  r_source_w_fm_generic(filename, dest_site, (wilson_vector *)NULL,
			spin, color, t0, source_type);
}

/*--------------------------------------------------------------------*/
void r_source_w_fm_to_field(char *filename, wilson_vector *dest_field,
			    int spin, int color, int t0, int source_type)
{
  r_source_w_fm_generic(filename, (field_offset)(-1), dest_field, 
			spin, color, t0, source_type);
}

