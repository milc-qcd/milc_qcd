/*************************** scidac_to_v5.c ************************/
/* MIMD version 7 */
/* Read SciDAC format gauge configuration, convert to MILC v5 format */
/* This is single-processor code. */
/* It converts on the fly, so the memory requirement is low */
/* C. DeTar 5/16/05 */

/* Usage ...

   scidac_to_v5 scidac_file milc_file

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

#define LATDIM 4
#define PARALLEL 1
#define SERIAL 0

#define NODE_DUMP_ORDER 1
#define NATURAL_ORDER 0

#define MAX_BUF_LENGTH 4096

typedef struct {
  gauge_file *gf;
  int buf_length;
  int siterank;
  int rank29,rank31;
  fsu3_matrix * lbuf;
  off_t checksum_offset;
  int prec;
} w_serial_site_writer;

void d2f_4mat(su3_matrix *a, fsu3_matrix *b);
void write_checksum(int parallel, gauge_file *gf);
void swrite_gauge_hdr(FILE *fp, gauge_header *gh);
void w_serial_f(gauge_file *gf);

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
      /**gen_pt[i] = (char **)malloc(sites_on_node*sizeof(char *) );
        if(gen_pt[i]==NULL){
            printf("NODE %d: no room for pointer vector\n",this_node);
            terminate(1);
	    } **/
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

static fsu3_matrix *outputbuf;


/*----------------------------------------------------------------------*/
void w_serial_start_lattice(gauge_file *gf, w_serial_site_writer *state,
			    int input_prec)
{
  /* gf  = file descriptor as opened by w_serial_i */
  /* state of the writer for a single site */

  FILE *fp;
  gauge_header *gh;
  off_t offset;             /* File stream pointer */
  off_t coord_list_size;    /* Size of coordinate list in bytes */
  off_t head_size;          /* Size of header plus coordinate list */
  off_t checksum_offset;    /* Location of checksum */
  off_t gauge_check_size;   /* Size of checksum record */

  if(gf->parallel)
    printf("w_serial_start_lattice: Attempting serial write to parallel file \n");
  
  outputbuf = (fsu3_matrix *)malloc(MAX_BUF_LENGTH*4*sizeof(fsu3_matrix));
  if(outputbuf == NULL)
    {
      printf("w_serial: Node 0 can't malloc outputbuf\n"); 
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
  
  if( g_seek(fp,offset,SEEK_SET) < 0 ) 
    {
      printf("w_serial: Node %d g_seek %lld failed error %d file %s\n",
	     this_node,(long long)offset,errno,gf->filename);
      fflush(stdout);terminate(1);
    }

  /* initialize checksums */
  gf->check.sum31 = 0;
  gf->check.sum29 = 0;

  state->gf              = gf;
  state->buf_length      = 0;
  state->siterank        = 0;
  state->rank29          = 0;
  state->rank31          = 0;
  state->lbuf            = outputbuf;
  state->checksum_offset = checksum_offset;
  state->prec            = input_prec;
} /* w_serial_start_lattice */

/*----------------------------------------------------------------------*/
/* Factory function for the SciDAC reader.  Writes the data to the
   MILC v5 file in the order this function is called regardless of "index". */
void w_serial_site_links(char *buf, size_t index, int count, void *arg)
{
  w_serial_site_writer *state = (w_serial_site_writer *)arg;
  gauge_file *gf    = state->gf;
  int buf_length    = state->buf_length;
  int siterank      = state->siterank;
  int rank29        = state->rank29;
  int rank31        = state->rank31;
  fsu3_matrix *lbuf = state->lbuf;
  fsu3_matrix *fbuf;
  dsu3_matrix *dbuf;
  int input_prec    = state->prec;
  int dir,i,j;

  FILE *fp = gf->fp;
  u_int32type *val;
  register int k;

  if(count != 4){
    printf("w_serial_site_links: expecting 4 color matrices but got %d\n",
	   count);
    terminate(1);
  }

  /* Copy buf to lbuf, converting precision if necessary */
  
  for(dir = 0; dir < 4; dir++){
    for(i = 0; i < 3; i++)for(j = 0; j < 3; j++)
      {
	if(input_prec == 1){
	  fbuf = (fsu3_matrix *)buf;
	  lbuf[4*buf_length+dir].e[i][j].real = fbuf[dir].e[i][j].real;
	  lbuf[4*buf_length+dir].e[i][j].imag = fbuf[dir].e[i][j].imag;
	} else { /* input_prec == 2 */
	  dbuf = (dsu3_matrix *)buf;
	  lbuf[4*buf_length+dir].e[i][j].real = dbuf[dir].e[i][j].real;
	  lbuf[4*buf_length+dir].e[i][j].imag = dbuf[dir].e[i][j].imag;
	}
      }
  }


  /* Accumulate checksums - contribution from next site */
  for(k = 0, val = (u_int32type *)&lbuf[4*buf_length]; 
      k < 4*(int)sizeof(fsu3_matrix)/(int)sizeof(int32type); 
      k++, val++)
    {
      gf->check.sum29 ^= (*val)<<rank29 | (*val)>>(32-rank29);
      gf->check.sum31 ^= (*val)<<rank31 | (*val)>>(32-rank31);
      rank29++; if(rank29 >= 29)rank29 = 0;
      rank31++; if(rank31 >= 31)rank31 = 0;
    }
  
  buf_length++;

  if( (buf_length == MAX_BUF_LENGTH) || (siterank == volume-1))
    {
      /* write out buffer */
      
      if( (int)fwrite(lbuf,4*sizeof(fsu3_matrix),buf_length,fp) != buf_length)
	{
	  printf("w_serial_site_links: gauge configuration write error %d file %s\n",
		 errno,gf->filename); 
	  fflush(stdout);
	  terminate(1);   
	}
      buf_length = 0;		/* start again after write */
    }

  state->buf_length  = buf_length;
  state->siterank++;
  state->rank29      = rank29;
  state->rank31      = rank31;
  
} /* w_serial_site_links */

/*----------------------------------------------------------------------*/
void w_serial_finish_lattice(w_serial_site_writer *state)
{
  /* gf  = file descriptor as opened by w_serial_i */

  gauge_file *gf = state->gf;
  fsu3_matrix *lbuf = state->lbuf;
  FILE *fp = gf->fp;
  gauge_header *gh = gf->header;
  off_t checksum_offset = state->checksum_offset;    /* Location of checksum */
  
  free(lbuf);
  printf("Saved gauge configuration serially to binary file %s\n",
	 gf->filename);
  printf("Time stamp %s\n",gh->time_stamp);
  
  /* Write checksum */
  /* Position file pointer */
  if( g_seek(fp,checksum_offset,SEEK_SET) < 0 ) 
    {
      printf("w_serial: Node %d g_seek %lld failed error %d file %s\n",
	     this_node,(long long)checksum_offset,errno,gf->filename);
      fflush(stdout);terminate(1);
    }
  write_checksum(SERIAL,gf);

} /* w_serial_finish_lattice */

/*----------------------------------------------------------------------*/

EXTERN QIO_String *xml_record_in;

int main(int argc, char *argv[])
{

  int ndim,dims[4];
  gauge_file *gf;
  gauge_header *gh;
  FILE *fp;
  char *filename_milc,*filename_scidac;
  QIO_Layout layout;
  QIO_Reader *infile;
  QIO_RecordInfo rec_info;
  char *datatype;
  int status;
  int datum_size;
  int input_prec;
  int count = 4;
  int word_size;
  int typesize;
  w_serial_site_writer state;
  
  if(argc < 3)
    {
      fprintf(stderr,"Usage %s <SciDAC file> <MILC file>\n",argv[0]);
      exit(1);
    }
  filename_scidac = argv[1];
  filename_milc   = argv[2];

  if(this_node == 0)printf("Converting file %s to MILC v5 file %s\n",
			   filename_scidac, filename_milc);

  initialize_machine(&argc,&argv);

  this_node = mynode();
  number_of_nodes = numnodes();

  if(number_of_nodes != 1){
    printf("This is single-processor code. Please rebuild as such.\n");
    terminate(1);
  }

  /* Open the SciDAC file and discover the lattice dimensions.  Then
     close. */

  status = read_lat_dim_scidac(filename_scidac, &ndim, dims);
  if(status)terminate(1);
  
  if(ndim != 4){
    printf("Wanted ndims = 4 in %s but got %d\n",filename_scidac,ndim);
    terminate(1);
  }

  nx = dims[0]; ny = dims[1]; nz = dims[2]; nt = dims[3];
  volume = nx*ny*nz*nt;

  /* Finish setting up, now we know the dimensions */
  setup();

  /* Build the QIO layout */
  build_qio_layout(&layout);

  /* Open the SciDAC file for reading */
  infile = open_scidac_input(filename_scidac, &layout, 0, QIO_SERIAL);
  if(infile == NULL)terminate(1);

  /* Open the MILC v5 file for writing */
  fp = fopen(filename_milc, "wb");
  if(fp == NULL)
    {
      printf("Can't open file %s, error %d\n",
	     filename_milc,errno);fflush(stdout);
      terminate(1);
    }
  gf = setup_output_gauge_file();
  gh = gf->header;

  /* Read the SciDAC record header. */
  xml_record_in = QIO_string_create();
  status = QIO_read_record_info(infile, &rec_info, xml_record_in);
  if(status != QIO_SUCCESS)terminate(1);
  node0_printf("Record info \n\"%s\"\n",QIO_string_ptr(xml_record_in));

  /* Make sure this is a lattice field */
  datatype = QIO_get_datatype(&rec_info);
  typesize = QIO_get_typesize(&rec_info);
  if(strcmp(datatype, "QDP_F3_ColorMatrix") == 0 ||
     strcmp(datatype, "USQCD_F3_ColorMatrix") == 0 ||
     typesize == 72){
    datum_size = sizeof(fsu3_matrix);  
    input_prec = 1;
    word_size = sizeof(float);
  }
  else if(strcmp(datatype, "QDP_D3_ColorMatrix") == 0 ||
	  strcmp(datatype, "USQCD_F3_ColorMatrix") == 0 ||
	  typesize == 144){
    datum_size = sizeof(dsu3_matrix);  
    input_prec = 2;
    word_size = sizeof(double);
  }
  else {
    printf("Unrecognized datatype %s\n",datatype);
    terminate(1);
  }

  /* Copy the time stamp from the SciDAC file */
  strncpy(gh->time_stamp, QIO_get_record_date(&rec_info), 
	  MAX_TIME_STAMP);
  gh->time_stamp[MAX_TIME_STAMP-1] = '\0';

  /* Write the MILC v5 header */
  gh->order = NATURAL_ORDER;

  /* Node 0 writes the header */
  
  swrite_gauge_hdr(fp,gh);
  
  /* Assign values to file structure */

  gf->fp = fp; 
  gf->filename = filename_milc;
  gf->byterevflag    = 0;            /* Not used for writing */
  gf->rank2rcv       = NULL;         /* Not used for writing */
  gf->parallel       = 0;

  /* Initialize writing the lattice data */
  w_serial_start_lattice(gf, &state, input_prec);

  /* Read the SciDAC record data.  The factory function writes the
     site links to a file. */

  status = QIO_read_record_data(infile, w_serial_site_links, 
				datum_size*count, word_size, 
				(void *)&state);
  if(status != QIO_SUCCESS)terminate(1);

  node0_printf("SciDAC checksums %x %x\n",
	       QIO_get_reader_last_checksuma(infile),
	       QIO_get_reader_last_checksumb(infile));

  /* Close the SciDAC file */
  QIO_close_read(infile);

  /* Finish the MILC v5 file */
  w_serial_finish_lattice(&state);

  w_serial_f(gf);

  QIO_string_destroy(xml_record_in);

  normal_exit(0);

  return 0;
}
