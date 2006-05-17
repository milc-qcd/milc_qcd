/*************************** v5_to_scidac.c ************************/
/* MIMD version 7 */
/* Read gauge configuration, convert MILD v5 to SciDAC format */
/* Optionally include ILDG records in the conversion */
/* C. DeTar 1/30/05 */
/* C. DeTar 5/4/06 support for ILDG format */

/* Usage ...

   lattice_to_scidac [--ildg] milc_file scidac_file
   [LFN string]

   The LFN string is used as the logical file name for ILDG usage.
   If the MILC info file is also present, it is copied into the
   user record XML.
*/

#define CONTROL

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <math.h>
#include "../include/complex.h"
#include "../include/su3.h"
#include <lattice.h>
#include "../include/macros.h"
#include "../include/comdefs.h"
#include <fcntl.h>
#include <errno.h>
#include <time.h>
#include <string.h>
#include "../include/io_lat.h"
#include "../include/io_wprop.h"
#include "../include/generic.h"
#include "../include/io_scidac.h"
#include <qio.h>
#include <dml.h>

#define LATDIM 4
#define PARALLEL 1
#define SERIAL 0
#define MAX_ILDGLFN 513
#define MAX_INFO 4097
#define MAX_INFO_LINE 129
#define MAX_INFO_FILENAME 513
#define MAX_RECXML_LENGTH 4097

#define NODE_DUMP_ORDER 1
#define NATURAL_ORDER 0

#define MAX_BUF_LENGTH 4096

typedef struct {
  gauge_file *gf;
  int buf_length;
  int where_in_buf;
  int siterank;
  int rank29,rank31;
  fsu3_matrix * lbuf;
  off_t checksum_offset;
  gauge_check test_gc;
  uint32_t crc;
} r_serial_site_reader;

void read_checksum(int parallel, gauge_file *gf, gauge_check *test_gc);
gauge_file *r_serial_i(char *filename);
void r_serial_f(gauge_file *gf);

/*----------------------------------------------------------------------*/
/* Read the lattice info file if it exists and put its contents in "buf"
   with tags suitable for the user record XML */
int read_info(char *buf, char *latfilename, int maxlength){
  char infofilename[MAX_INFO_FILENAME];
  FILE *infofp;
  char begin[] = "<?xml version=\"1.0\" encoding=\"UTF-8\"?><info>";
  char end[] = "</info>";
  char missing[] = "missing";
  char infoline[MAX_INFO_LINE];
  int length = 0;

  /* Construct file name */
  if(strlen(latfilename) + 5 >= MAX_INFO_FILENAME){
    fprintf(stderr,"Not enough room for info filename\n");
    return 1;
  }
  strcpy(infofilename,latfilename);
  strcat(infofilename,".info");

  /* Start the buffer */
  length += strlen(begin);
  if(length < maxlength)
    strcpy(buf, begin);

  /* Look for info file in directory */
  infofp = fopen(infofilename,"r");
  if(infofp == NULL){
    fprintf(stderr,"Can't find the info file %s\nIgnoring it.\n",
	    infofilename);
  }

  /* Read and copy if we have a file */
  if(infofp != NULL){
    while(feof(infofp) == 0){
      infoline[0] = '\0';
      fgets(infoline, MAX_INFO_LINE, infofp);
      length += strlen(infoline);
      if(length < maxlength)
	strcat(buf, infoline);
    }
  }
  else{
      length += strlen(missing);
      if(length < maxlength)
	strcat(buf, missing);
  }

  /* Finish the buffer */
  length += strlen(end);
  if(length < maxlength)
    strcat(buf, end);

  if(length >= maxlength){
    fprintf(stderr,"Not enough room for the info file data\n");
    return 1;
  }
  return 0;
}

/*----------------------------------------------------------------------*/
void make_lattice(){
register int i;               /* scratch */
int x,y,z,t;            /* coordinates */
    /* allocate space for lattice, fill in parity, coordinates and index.  */
    lattice = (site *)malloc( sites_on_node * sizeof(site) );
    if(lattice==NULL){
        printf("NODE %d: no room for lattice\n",this_node);
        terminate(1);
    }
   /* Allocate address vectors */
    for(i=0;i<8;i++){
      /** gen_pt[i] = (char **)malloc(sites_on_node*sizeof(char *) );
        if(gen_pt[i]==NULL){
            printf("NODE %d: no room for pointer vector\n",this_node);
            terminate(1);
	    }**/
      gen_pt[i] = NULL;
    }

    for(t=0;t<nt;t++)for(z=0;z<nz;z++)for(y=0;y<ny;y++)for(x=0;x<nx;x++){
        if(node_number(x,y,z,t)==mynode()){
            i=node_index(x,y,z,t);
            lattice[i].x=x;     lattice[i].y=y; lattice[i].z=z; lattice[i].t=t;
            lattice[i].index = x+nx*(y+ny*(z+nz*t));
            if( (x+y+z+t)%2 == 0)lattice[i].parity=EVEN;
            else                 lattice[i].parity=ODD;
#ifdef SITERAND
            initialize_prn( &(lattice[i].site_prn) , iseed, lattice[i].index);
#endif
        }
    }
}
/*----------------------------------------------------------------------*/

void setup() {

  /* Set up lattice */
  broadcast_bytes((char *)&nx,sizeof(int));
  broadcast_bytes((char *)&ny,sizeof(int));
  broadcast_bytes((char *)&nz,sizeof(int));
  broadcast_bytes((char *)&nt,sizeof(int));
  
  setup_layout();
  make_lattice();
}

/*----------------------------------------------------------------------*/

void r_serial_start_lattice(gauge_file *gf, r_serial_site_reader *state)
{
  /* gf  = gauge configuration file structure */
  /* state of the writer for a single site */

  FILE *fp;
  gauge_header *gh;
  fsu3_matrix *lbuf;
  off_t offset;             /* File stream pointer */
  off_t coord_list_size;    /* Size of coordinate list in bytes */
  off_t head_size;          /* Size of header plus coordinate list */
  off_t checksum_offset;    /* Location of checksum */
  off_t gauge_check_size;   /* Size of checksum record */

  char myname[] = "r_serial_start_lattice";
  
  fp = gf->fp;
  gh = gf->header;

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
  
  coord_list_size = 0;
  if(gf->header->order != NATURAL_ORDER)
    {
      printf("%s: Can't convert a file that is not in natural order\n",
	     myname);
      terminate(1);
    }
  checksum_offset = gf->header->header_bytes + coord_list_size;
  head_size = checksum_offset + gauge_check_size;
  
  /* Allocate space for read buffer */
  
  if(gf->parallel)
    printf("%s: Attempting serial read from parallel file \n",myname);
  
  /* Allocate single precision read buffer.  We will read only one site
     at a time */
  lbuf = (fsu3_matrix *)malloc(MAX_BUF_LENGTH*4*sizeof(fsu3_matrix));
  if(lbuf == NULL)
    {
      printf("%s: Node %d can't malloc lbuf\n",myname,this_node);
      fflush(stdout);
      terminate(1);
    }
  
  /* Position file for reading gauge configuration */
  
  offset = head_size;
  
  if( fseeko(fp,offset,SEEK_SET) < 0 ) 
    {
      printf("%s: Node 0 fseeko %lld failed error %d file %s\n",
	     myname,(long long)offset,errno,gf->filename);
      fflush(stdout);terminate(1);   
    }

  /* Save state */
  state->gf              = gf;
  state->buf_length      = 0;
  state->where_in_buf    = 0;
  state->siterank        = 0;
  state->rank29          = 0;
  state->rank31          = 0;
  state->lbuf            = lbuf;
  state->checksum_offset = checksum_offset;
  state->test_gc.sum29   = 0;
  state->test_gc.sum31   = 0;
  state->crc             = 0;

}
/*----------------------------------------------------------------------*/
/* Factory function for the SciDAC reader.  Reads the data from the
   MILC v5 file in the order this function is called regardless of "index". */

void r_serial_site_links(char *buf, size_t index, int count, void *arg)
{
  r_serial_site_reader *state = (r_serial_site_reader *)arg;
  gauge_file *gf    = state->gf;
  int buf_length    = state->buf_length;
  int where_in_buf  = state->where_in_buf;
  fsu3_matrix *lbuf = state->lbuf;
  int x,y,z,t;
  int rcv_coords    = state->siterank;
  char myname[]     = "r_serial_site_links";

  FILE *fp = gf->fp;
  u_int32type *val;
  int k;

  x = rcv_coords % nx;   rcv_coords /= nx;
  y = rcv_coords % ny;   rcv_coords /= ny;
  z = rcv_coords % nz;   rcv_coords /= nz;
  t = rcv_coords % nt;

  if(node_index(x,y,z,t) != index){
    printf("%s: expected index %d but got index %d\n",
	   myname,state->siterank,index);
  }

  /* Fill buffer, if necessary */
  if(where_in_buf == buf_length)
    {  /* get new buffer */
      /* new buffer length  = remaining sites, but never bigger 
	 than MAX_BUF_LENGTH */
      buf_length = volume - state->siterank;
      if(buf_length > MAX_BUF_LENGTH)buf_length = MAX_BUF_LENGTH;
      /* then do read */
      
      if( (int)fread(lbuf,4*sizeof(fsu3_matrix),buf_length,fp) 
	  != buf_length)
	{
	  printf("%s: gauge configuration read error %d file %s\n",
		 myname,errno,gf->filename); 
	  fflush(stdout); terminate(1);
	}
      where_in_buf = 0;  /* reset counter */
    }  /*** end of the buffer read ****/
  
  /* Save 4 matrices in buf for further processing */
  /* We always work in single precision */
  memcpy(buf,&lbuf[4*where_in_buf],4*sizeof(fsu3_matrix));
  where_in_buf++;
      
  /* Do the byte reversal and then checksum,
     if needed.  At this point tmpsu3 contains the input matrices
     and buf is where we want to put them. */
  
  if(gf->byterevflag==1)
    byterevn((int32type *)buf,
	     4*sizeof(fsu3_matrix)/sizeof(int32type));
  /* Accumulate checksums */
  for(k = 0, val = (u_int32type *)buf; 
      k < 4*(int)sizeof(fsu3_matrix)/(int)sizeof(int32type); 
      k++, val++)
    {
      state->test_gc.sum29 ^= (*val)<<state->rank29 | 
	(*val)>>(32-state->rank29);
      state->test_gc.sum31 ^= (*val)<<state->rank31 |
	(*val)>>(32-state->rank31);
      state->rank29++; if(state->rank29 >= 29)state->rank29 = 0;
      state->rank31++; if(state->rank31 >= 31)state->rank31 = 0;
    }
  state->crc = 
    DML_crc32(state->crc, (char *)buf, 4*(int)sizeof(fsu3_matrix));
  state->buf_length  = buf_length;
  state->where_in_buf = where_in_buf;
  state->siterank++;

} /* r_serial_site_links */

/*----------------------------------------------------------------------*/
void r_serial_finish_lattice(r_serial_site_reader *state)
{
  /* gf  = file descriptor as opened by r_serial_i */

  gauge_file *gf = state->gf;
  fsu3_matrix *lbuf = state->lbuf;
  FILE *fp = gf->fp;
  gauge_header *gh = gf->header;
  off_t checksum_offset = state->checksum_offset;    /* Location of checksum */
  char myname[] = "r_serial_finish_lattice";

  /* Read and verify checksum */
  /* Checksums not implemented until version 5 */
  
  printf("Restored binary gauge configuration serially from file %s\n",
	 gf->filename);
  if(gh->magic_number == GAUGE_VERSION_NUMBER)
    {
      printf("Time stamp %s\n",gh->time_stamp);
      if( fseeko(fp,checksum_offset,SEEK_SET) < 0 ) 
	{
	  printf("%s: Node 0 fseeko %lld failed error %d file %s\n",
		 myname,(long long)checksum_offset,errno,gf->filename);
	  fflush(stdout);terminate(1);   
	}
      read_checksum(SERIAL,gf,&(state->test_gc));
    }
  fflush(stdout);
  free(lbuf);
}

/*----------------------------------------------------------------------*/

int main(int argc, char *argv[])
{

  int ndim,dims[4];
  gauge_file *gf;
  FILE *fp;
  char *filename_milc,*filename_scidac;
  QIO_Layout layout;
  QIO_Writer *outfile;
  QIO_RecordInfo *rec_info;
  int status;
  int ildgstyle;
  /* We assume input precision is single */
  int datum_size = sizeof(fsu3_matrix);
  int count = 4;
  int word_size = sizeof(float);
  int length;
  r_serial_site_reader state;
  QIO_String *xml_record_out;
  char ildg_lfn[MAX_ILDGLFN];
  char recxml[MAX_RECXML_LENGTH];
  char filexml[] = "<?xml version=\"1.0\" encoding=\"UTF-8\"?><title>MILC ILDG archival gauge configuration</title>";

  if(argc < 3 || 
     (argc < 4 && argv[1][0] == '-') || 
     (argc >= 4 && (strcmp(argv[1],"--ildg") != 0)))
    {
      fprintf(stderr,"Usage %s [--ildg] <MILC file> <SciDAC file>\n",argv[0]);
      return 1;
    }

  if(argc < 4){
    filename_milc   = argv[1];
    filename_scidac = argv[2];
    ildgstyle = QIO_ILDGNO;
  }
  else{
    filename_milc   = argv[2];
    filename_scidac = argv[3];
    ildgstyle = QIO_ILDGLAT;
    if(fgets(ildg_lfn, MAX_ILDGLFN, stdin) == NULL){
      fprintf(stderr,"Couldn't read the LFN\n");
      return 1;
    }
    else{
      /* Chop end-of-line character */
      length = strlen(ildg_lfn);
      if(ildg_lfn[length-1] == '\n'){
	ildg_lfn[length-1] = '\0';
      }
    }
  }

  if(read_info(recxml, filename_milc, MAX_RECXML_LENGTH) != 0)return 1;


  if(this_node == 0)printf("Converting MILC v5 file %s to SciDAC file %s\n",
			   filename_milc, filename_scidac);
  if(ildgstyle == QIO_ILDGLAT)
    printf("in ILDG compatible format with LFN\n%s\n",ildg_lfn);

  initialize_machine(argc,argv);
#ifdef HAVE_QDP
  QDP_initialize(&argc, &argv);
#endif

  this_node = mynode();
  number_of_nodes = numnodes();

  if(number_of_nodes != 1){
    printf("This is single-processor code. Please rebuild as such.\n");
    terminate(1);
  }

  /* Open the MILC file and discover the lattice dimensions.  Then
     close. */

  read_lat_dim_gf(filename_milc, &ndim, dims);
  
  if(ndim != 4){
    printf("Wanted ndims = 4 in %s but got %d\n",filename_milc,ndim);
    terminate(1);
  }

  nx = dims[0]; ny = dims[1]; nz = dims[2]; nt = dims[3];
  volume = nx*ny*nz*nt;

  /* Finish setting up, now we know the dimensions */
  setup();

  /* Build the QIO layout */
  build_qio_layout(&layout);

  /* Open the MILC v5 file for reading */
  fp = fopen(filename_milc, "r");
  if(fp == NULL)
    {
      printf("Can't open file %s, error %d\n",
	     filename_milc,errno);fflush(stdout);
      terminate(1);
    }

  /* Read MILC header */
  gf = r_serial_i(filename_milc);

  /* Open the SciDAC file for writing */
  outfile = open_scidac_output(filename_scidac, QIO_SINGLEFILE,
			       QIO_SERIAL, ildgstyle, ildg_lfn, &layout,
			       filexml);
  if(outfile == NULL)terminate(1);

  /* Initialize reading the MILC lattice data */
  r_serial_start_lattice(gf, &state);

  /* Write the SciDAC record. The factory function reads the
     site links from the MILC file */

  xml_record_out = QIO_string_create();
  QIO_string_set(xml_record_out, recxml);
  rec_info = QIO_create_record_info(QIO_FIELD, "QDP_F3_ColorMatrix", "F", 
				    3, 0, datum_size, 4);
  status = QIO_write(outfile, rec_info, xml_record_out,
		     r_serial_site_links, datum_size*count, word_size, 
		     (void *)&state);
  if(status != QIO_SUCCESS)terminate(1);

  node0_printf("SciDAC checksums %x %x\n",
	       QIO_get_writer_last_checksuma(outfile),
	       QIO_get_writer_last_checksumb(outfile));
  printf("crc32 checksum node %d %lu\n",this_node,state.crc);

  /* Close the SciDAC file */
  QIO_close_write(outfile);

  /* Finish the MILC v5 file */
  r_serial_finish_lattice(&state);

  r_serial_f(gf);

  QIO_string_destroy(xml_record_out);

#ifdef HAVE_QDP
  QDP_finalize();
#endif  
  normal_exit(0);

  return 0;
}
