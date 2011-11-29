/**************************** io_source_w_fm.c ****************************/
/* MIMD version 7 */
/*
 *  Load the Fermilab format file for Wilson propagator sources
 *  This version reads 12 source vectors for each site in
 *  the so-called FNAL "timeslice" format, i.e.
 *  one 3D Dirac field for each of 12 source colors and spins.
 *  It converts the canopy Dirac basis to MILC Dirac basis.
 *  Then it returns one Dirac vector of the requested spin and color.
 *  This is inefficient, since the source must be reread and retransformed
 *  for each spin and color.  But we are trying to get folks to use
 *  the SciDAC/USQCD source format files instead.
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
  wilson_propagator *source;
} w_source_file;

#ifndef PAD_SEND_BUF
#define PAD_SEND_BUF 8
#endif
#define NATURAL_ORDER 0

static int 
read_w_fm_source_hdr(w_source_file *wsf)
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

  if( size_of_element != sizeof(float) ||
      elements_per_site != 24 /* wilson_vector */)
    {	
      printf(" File %s is not a Dirac field",  wsf->filename);
      printf(" got size_of_element %d and elements_per_site %d\n",
	     size_of_element, elements_per_site);
      terminate(1);
    }
  
  /* The site order parameter is ignored */
  
  if(sread_byteorder(byterevflag,wsf->fp,&order,sizeof(int32type),
		     myname,"order parameter"))terminate(1);
  
  return byterevflag;
}
 
static w_source_file *
setup_input_w_source_file(char *filename)
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
r_source_w_fm_i(char *filename)
{
  /* Returns file descriptor for opened file */

  w_source_file   *wsf;
  w_source_header *wsh;
  int byterevflag;

  wsf = setup_input_w_source_file(filename);
  wsh = wsf->header;

  if(this_node==0){
    wsf->fp = g_open(filename,"rb");
    if(wsf->fp==NULL){
      printf("Can't open source file %s, error %d\n",filename,errno);
      terminate(1);
    }
    byterevflag = read_w_fm_source_hdr(wsf);

  }
  else wsf->fp = NULL;

  /* Make room for the whole Dirac source */
  wsf->source = 
    (wilson_propagator *)malloc(sizeof(wilson_propagator)*sites_on_node);
  memset(wsf->source, '\0', sizeof(wilson_propagator)*sites_on_node);

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
	      int spin, int color, int x0, int y0, int z0, 
	      int t0)
{
  int rcv_rank, rcv_coords, status;
  int destnode;
  int x,y,z,t,i=0, byterevflag, c0,s0,c1,s1,a;
  int tmin, tmax;
  site *s;

  struct {
    fwilson_vector q;
    char pad[PAD_SEND_BUF];    /* Introduced because some switches
				  perform better if message lengths
				  are longer */
  } wmsg;

  int buf_length=0, where_in_buf=0;
  int irecord, nrecord;
  fwilson_vector *wbuff=NULL;
  fwilson_vector wfix;
  wilson_vector *wv;
  wilson_propagator *wp;
  int vol3 = nx*ny*nz;

  byterevflag = wsf->byterevflag;

  if(this_node == 0)
    {
      wbuff = (fwilson_vector *)malloc(MAX_BUF_LENGTH*sizeof(fwilson_vector));
      if(wbuff == NULL)
	{
	  printf("Node %d can't malloc wbuff\n",this_node);
	  fflush(stdout);
	  terminate(1);
	}
      
      buf_length = 0;
      where_in_buf = 0;

    } /* end of if(this_node == 0)*/

  nrecord = 12;

  g_sync();

  /* If requested, we replicate the source function on ALL time slices */
  /* Otherwise we read only to time slice t0 */

  if(t0 == ALL_T_SLICES){ tmin = 0; tmax = nt-1; }
  else      { tmin = t0; tmax = t0; }

  /* Node 0 reads and deals out the values */
  status = 0;
  /* Iterate only over timeslice t0 */
  for(irecord = 0; irecord < nrecord; irecord++)
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
	      a=(int)g_read(wbuff,sizeof(fwilson_vector),buf_length,wsf->fp);
	      
	      if( a != buf_length)
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
	  wmsg.q = wbuff[where_in_buf];
	}
	  
	/* Loop either over all time slices or just one time slice */
	for(t = tmin; t <= tmax; t++){

	  destnode=node_number(x,y,z,t);
	  
	  if(this_node == 0){
	    /* node 0 doesn't send to itself */
	    if(destnode != 0){
	      /* send to correct node */
	      send_field((char *)&wmsg, sizeof(wmsg), destnode);
	    }
	  } /* if(this_node==0) */
	  else {	/* for all nodes other than node 0 */
	    if(this_node==destnode){
	      get_field((char *)&wmsg, sizeof(wmsg),0);
	    }
	  }
	  
	  /* The receiving node does the byte reversal.  At this point msg
	     contains the input vectors and i points to the destination
	     site structure */
	  
	  if(this_node==destnode)
	  {
	    wfix = wmsg.q;
	    if(byterevflag==1)
	      byterevn((int32type *)&wfix, 
		       sizeof(fwilson_vector)/sizeof(int32type));
	      
	      /* Now copy the site data into the destination converting
		 to generic precision if needed */
	      
	      i = node_index(x,y,z,t);
	      
	      /* We copy first to our buffer, filling all components */
	      wp = wsf->source + i;
	      for(s1=0;s1<4;s1++)for(c1=0;c1<3;c1++)
  	         {
		    c0 = irecord % 3;
		    s0 = irecord / 3;
		    wp->c[c0].d[s0].d[s1].c[c1].real = wfix.d[s1].c[c1].real;
		    wp->c[c0].d[s0].d[s1].c[c1].imag = wfix.d[s1].c[c1].imag;
		 }
	  } /* if this_node == destnode */
	} /* t */
	
	where_in_buf++;
	
      }  /* rcv_rank, irecord */
  
  if(wbuff != NULL)free(wbuff); wbuff = NULL;
  
  /* Finally convert spin basis from FNAL to MILC and then copy the
     requested MILC spin, color to the user Wilson vector */
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
		      field_offset dest_site, wilson_vector *dest_field, 
		      int spin, int color, int x0, int y0, int z0, int t0)
{
  w_source_file *wsf;
  wsf = r_source_w_fm_i(filename);
  r_source_w_fm(wsf, dest_site, dest_field, spin, color, x0, y0, z0, t0);
  r_source_w_fm_f(wsf);
  
}

/*--------------------------------------------------------------------*/
void r_source_w_fm_to_site(char *filename, field_offset dest_site,
			   int spin, int color, int x0, int y0, int z0, int t0)
{
  r_source_w_fm_generic(filename, dest_site, (wilson_vector *)NULL,
			spin, color, x0, y0, z0, t0);
}

/*--------------------------------------------------------------------*/
void r_source_w_fm_to_field(char *filename, wilson_vector *dest_field,
			    int spin, int color, int x0, int y0, int z0, int t0)
{
  r_source_w_fm_generic(filename, (field_offset)(-1), dest_field, 
			spin, color, x0, y0, z0, t0);
}

