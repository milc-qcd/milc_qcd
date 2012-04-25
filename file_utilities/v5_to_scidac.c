/*************************** v5_to_scidac.c ************************/
/* MIMD version 7 */
/* Read gauge configuration, convert MILD v5 to SciDAC format */
/* This is single processor code. */
/* It converts on the fly, so the memory requirement is low */
/* Optionally include ILDG records in the conversion */
/* C. DeTar 1/30/05 */
/* C. DeTar 5/4/06 support for ILDG format */

/* Usage ...

   v5_to_scidac [--ildg] milc_file scidac_file
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
#define MAX_STRING 32

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

/*--------------------------------------------------------------------*/
/* Parse MILC info line for a desired value */
void read_info_val(char tag[], float *val, int *found, char *line)
{

  if(strstr(line,tag) == NULL)return;
  /* Support two line formats */
  if(strstr(line,"=") == NULL)
    {
      /* TAG value */
      if(sscanf(line,"%*s %f",val) == 1)*found = 1;
    }
  else
    {
      /* TAG = value */
      if(sscanf(line,"%*s %*s %f",val) == 1)*found = 1;
    }
}

/*--------------------------------------------------------------------*/
/* Parse MILC info line for a desired value */
void read_info_valstring(char tag[], char *valstring, int *found, char *line)
{

  if(strstr(line,tag) == NULL)return;
  /* Support two line formats */
  if(strstr(line,"=") == NULL)
    {
      /* TAG value */
      if(sscanf(line,"%*s %s",valstring) == 1)*found = 1;
    }
  else
    {
      /* TAG = value */
      if(sscanf(line,"%*s %*s %s",valstring) == 1)*found = 1;
    }
}

/*----------------------------------------------------------------------*/
/* Print string formatting with digits based on intended precision */
void print_prec(char string[], size_t n, Real value, int prec){
  if(prec == 1){
    snprintf(string,n,"%.6e",value);  /* single precision */
  }
  else
    snprintf(string,n,"%.15e",value); /* double precision */
}

static QIO_String *xml_record;

/*----------------------------------------------------------------------*/
/* Read the lattice info file if it exists and put its contents in "buf"
   with tags suitable for the user record XML */
QIO_String *create_recxml(char *latfilename){
  char buf[MAX_RECXML_LENGTH];
  int maxlength = MAX_RECXML_LENGTH;
  char infofilename[MAX_INFO_FILENAME];
  FILE *infofp;
  QIO_USQCDLatticeInfo *record_info;
  char missing[] = "missing";
  char infoline[MAX_INFO_LINE];
  int length = 0;
  int foundssplaq = 0, foundstplaq = 0, foundlinktr = 0;
  float ssplaq = 0, stplaq = 0;
  char plaqstring[MAX_STRING];
  char linktrstring[MAX_STRING];

  /* Construct file name */
  if(strlen(latfilename) + 5 >= MAX_INFO_FILENAME){
    fprintf(stderr,"Not enough room for info filename\n");
    return NULL;
  }
  strcpy(infofilename,latfilename);
  strcat(infofilename,".info");

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
      /* Look for needed metadata */
      read_info_val("gauge.ssplaq", &ssplaq, &foundssplaq, infoline);
      read_info_val("gauge.stplaq", &stplaq, &foundstplaq, infoline);
      read_info_valstring("gauge.nersc_linktr", linktrstring, 
			  &foundlinktr, infoline);
    }
  }
  else{
      length += strlen(missing);
      if(length < maxlength)
	strcat(buf, missing);
  }

  if(length >= maxlength){
    fprintf(stderr,"Not enough room for the info file data\n");
    return NULL;
  }

  printf("found ss %d st %d tr %d\n",foundssplaq,foundstplaq,foundlinktr);

  /* Append the USQCD tags and values if we have them */
  if(foundssplaq && foundstplaq)
    print_prec(plaqstring, 32, (ssplaq+stplaq)/6., 1);

  record_info = QIO_create_usqcd_lattice_info(plaqstring, linktrstring, buf);

  xml_record = QIO_string_create();
  QIO_encode_usqcd_lattice_info(xml_record, record_info);
  QIO_destroy_usqcd_lattice_info(record_info);
  
  return xml_record;
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
  
  if( g_seek(fp,offset,SEEK_SET) < 0 ) 
    {
      printf("%s: Node 0 g_seek %lld failed error %d file %s\n",
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
  state->crc             = 0xffffffffL;

}

/*------------------------------------------------------------------*/
/* Posix 1003.2 cksum algorithm taken from 
   http://fxr.googlebit.com/source/usr.bin/cksum/crc.c?v=8-CURRENT#L116
   modified for ILDG use  C. DeTar 10/25/2009 */


typedef unsigned long uLong;            /* At least 32 bits */

/*
  Generate a table for a byte-wise 32-bit CRC calculation on the polynomial:
  x^32+x^26+x^23+x^22+x^16+x^12+x^11+x^10+x^8+x^7+x^5+x^4+x^2+x+1.

  Polynomials over GF(2) are represented in binary, one bit per coefficient,
  with the lowest powers in the most significant bit.  Then adding polynomials
  is just exclusive-or, and multiplying a polynomial by x is a right shift by
  one.  If we call the above polynomial p, and represent a byte as the
  polynomial q, also with the lowest power in the most significant bit (so the
  byte 0xb1 is the polynomial x^7+x^3+x+1), then the CRC is (q*x^32) mod p,
  where a mod b means the remainder after dividing a by b.

  This calculation is done using the shift-register method of multiplying and
  taking the remainder.  The register is initialized to zero, and for each
  incoming bit, x^32 is added mod p to the register if the bit is a one (where
  x^32 mod p is p+x^32 = x^26+...+1), and the register is multiplied mod p by
  x (which is shifting right by one and adding x^32 mod p if the bit shifted
  out is a one).  We start with the highest power (least significant bit) of
  q and repeat for all eight bits of q.

  The table is simply the CRC of all possible eight bit values.  This is all
  the information needed to generate CRC's on data a byte at a time for all
  combinations of CRC register values and incoming bytes.
*/

/* ========================================================================
 * Table of CRC-32's of all single-byte values 
 */
static uLong crc_table[256] = {
  0x0L,
  0x04c11db7L, 0x09823b6eL, 0x0d4326d9L, 0x130476dcL, 0x17c56b6bL,
  0x1a864db2L, 0x1e475005L, 0x2608edb8L, 0x22c9f00fL, 0x2f8ad6d6L,
  0x2b4bcb61L, 0x350c9b64L, 0x31cd86d3L, 0x3c8ea00aL, 0x384fbdbdL,
  0x4c11db70L, 0x48d0c6c7L, 0x4593e01eL, 0x4152fda9L, 0x5f15adacL,
  0x5bd4b01bL, 0x569796c2L, 0x52568b75L, 0x6a1936c8L, 0x6ed82b7fL,
  0x639b0da6L, 0x675a1011L, 0x791d4014L, 0x7ddc5da3L, 0x709f7b7aL,
  0x745e66cdL, 0x9823b6e0L, 0x9ce2ab57L, 0x91a18d8eL, 0x95609039L,
  0x8b27c03cL, 0x8fe6dd8bL, 0x82a5fb52L, 0x8664e6e5L, 0xbe2b5b58L,
  0xbaea46efL, 0xb7a96036L, 0xb3687d81L, 0xad2f2d84L, 0xa9ee3033L,
  0xa4ad16eaL, 0xa06c0b5dL, 0xd4326d90L, 0xd0f37027L, 0xddb056feL,
  0xd9714b49L, 0xc7361b4cL, 0xc3f706fbL, 0xceb42022L, 0xca753d95L,
  0xf23a8028L, 0xf6fb9d9fL, 0xfbb8bb46L, 0xff79a6f1L, 0xe13ef6f4L,
  0xe5ffeb43L, 0xe8bccd9aL, 0xec7dd02dL, 0x34867077L, 0x30476dc0L,
  0x3d044b19L, 0x39c556aeL, 0x278206abL, 0x23431b1cL, 0x2e003dc5L,
  0x2ac12072L, 0x128e9dcfL, 0x164f8078L, 0x1b0ca6a1L, 0x1fcdbb16L,
  0x018aeb13L, 0x054bf6a4L, 0x0808d07dL, 0x0cc9cdcaL, 0x7897ab07L,
  0x7c56b6b0L, 0x71159069L, 0x75d48ddeL, 0x6b93dddbL, 0x6f52c06cL,
  0x6211e6b5L, 0x66d0fb02L, 0x5e9f46bfL, 0x5a5e5b08L, 0x571d7dd1L,
  0x53dc6066L, 0x4d9b3063L, 0x495a2dd4L, 0x44190b0dL, 0x40d816baL,
  0xaca5c697L, 0xa864db20L, 0xa527fdf9L, 0xa1e6e04eL, 0xbfa1b04bL,
  0xbb60adfcL, 0xb6238b25L, 0xb2e29692L, 0x8aad2b2fL, 0x8e6c3698L,
  0x832f1041L, 0x87ee0df6L, 0x99a95df3L, 0x9d684044L, 0x902b669dL,
  0x94ea7b2aL, 0xe0b41de7L, 0xe4750050L, 0xe9362689L, 0xedf73b3eL,
  0xf3b06b3bL, 0xf771768cL, 0xfa325055L, 0xfef34de2L, 0xc6bcf05fL,
  0xc27dede8L, 0xcf3ecb31L, 0xcbffd686L, 0xd5b88683L, 0xd1799b34L,
  0xdc3abdedL, 0xd8fba05aL, 0x690ce0eeL, 0x6dcdfd59L, 0x608edb80L,
  0x644fc637L, 0x7a089632L, 0x7ec98b85L, 0x738aad5cL, 0x774bb0ebL,
  0x4f040d56L, 0x4bc510e1L, 0x46863638L, 0x42472b8fL, 0x5c007b8aL,
  0x58c1663dL, 0x558240e4L, 0x51435d53L, 0x251d3b9eL, 0x21dc2629L,
  0x2c9f00f0L, 0x285e1d47L, 0x36194d42L, 0x32d850f5L, 0x3f9b762cL,
  0x3b5a6b9bL, 0x0315d626L, 0x07d4cb91L, 0x0a97ed48L, 0x0e56f0ffL,
  0x1011a0faL, 0x14d0bd4dL, 0x19939b94L, 0x1d528623L, 0xf12f560eL,
  0xf5ee4bb9L, 0xf8ad6d60L, 0xfc6c70d7L, 0xe22b20d2L, 0xe6ea3d65L,
  0xeba91bbcL, 0xef68060bL, 0xd727bbb6L, 0xd3e6a601L, 0xdea580d8L,
  0xda649d6fL, 0xc423cd6aL, 0xc0e2d0ddL, 0xcda1f604L, 0xc960ebb3L,
  0xbd3e8d7eL, 0xb9ff90c9L, 0xb4bcb610L, 0xb07daba7L, 0xae3afba2L,
  0xaafbe615L, 0xa7b8c0ccL, 0xa379dd7bL, 0x9b3660c6L, 0x9ff77d71L,
  0x92b45ba8L, 0x9675461fL, 0x8832161aL, 0x8cf30badL, 0x81b02d74L,
  0x857130c3L, 0x5d8a9099L, 0x594b8d2eL, 0x5408abf7L, 0x50c9b640L,
  0x4e8ee645L, 0x4a4ffbf2L, 0x470cdd2bL, 0x43cdc09cL, 0x7b827d21L,
  0x7f436096L, 0x7200464fL, 0x76c15bf8L, 0x68860bfdL, 0x6c47164aL,
  0x61043093L, 0x65c52d24L, 0x119b4be9L, 0x155a565eL, 0x18197087L,
  0x1cd86d30L, 0x029f3d35L, 0x065e2082L, 0x0b1d065bL, 0x0fdc1becL,
  0x3793a651L, 0x3352bbe6L, 0x3e119d3fL, 0x3ad08088L, 0x2497d08dL,
  0x2056cd3aL, 0x2d15ebe3L, 0x29d4f654L, 0xc5a92679L, 0xc1683bceL,
  0xcc2b1d17L, 0xc8ea00a0L, 0xd6ad50a5L, 0xd26c4d12L, 0xdf2f6bcbL,
  0xdbee767cL, 0xe3a1cbc1L, 0xe760d676L, 0xea23f0afL, 0xeee2ed18L,
  0xf0a5bd1dL, 0xf464a0aaL, 0xf9278673L, 0xfde69bc4L, 0x89b8fd09L,
  0x8d79e0beL, 0x803ac667L, 0x84fbdbd0L, 0x9abc8bd5L, 0x9e7d9662L,
  0x933eb0bbL, 0x97ffad0cL, 0xafb010b1L, 0xab710d06L, 0xa6322bdfL,
  0xa2f33668L, 0xbcb4666dL, 0xb8757bdaL, 0xb5365d03L, 0xb1f740b4L
};


/* ========================================================================= */
#define DO1(buf) crc = crc << 8 ^ crc_table[(int)(crc >> 24 ^ (*buf++)) & 0xff];
#define DO2(buf)  DO1(buf); DO1(buf);
#define DO4(buf)  DO2(buf); DO2(buf);
#define DO8(buf)  DO4(buf); DO4(buf);

/* ========================================================================= */
uLong crc32(uLong crc, const unsigned char *buf, size_t len)
{
  if (buf == NULL) return 0L;
  crc = crc ^ 0xffffffffL;
  while (len >= 8)
    {
      DO8(buf);
      len -= 8;
    }
  if (len) do {
      DO1(buf);
    } while (--len);
  return crc ^ 0xffffffffL;
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
  unsigned char cbuf[4*sizeof(fsu3_matrix)];
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
    printf("%s: expected index %d but got index %lu\n",
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
      
  /* Do the input byte reversal and then checksum, if needed. */
  
  if(gf->byterevflag==1)
    byterevn((int32type *)buf,
	     4*sizeof(fsu3_matrix)/sizeof(int32type));
  /* Accumulate MILC v5 checksums */
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

  /* For the crc checksum we have to get the output byte ordering right */
  memcpy(cbuf,buf,4*sizeof(fsu3_matrix));
  if(! big_endian())
    byterevn((int32type *)cbuf,4*sizeof(fsu3_matrix)/sizeof(int32type));
  state->crc = 
    crc32(state->crc, cbuf, 4*(int)sizeof(fsu3_matrix));

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
      if( g_seek(fp,checksum_offset,SEEK_SET) < 0 ) 
	{
	  printf("%s: Node 0 g_seek %lld failed error %d file %s\n",
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
  QIO_Filesystem fs;
  QIO_RecordInfo *rec_info;
  int status;
  int ildgstyle;
  /* We assume input precision is single */
  int datum_size = sizeof(fsu3_matrix);
  int count = 4;
  int word_size = sizeof(float);
  int length;
  off_t payload_bytes;
  r_serial_site_reader state;
  QIO_String *xml_record_out;
  char ildg_lfn[MAX_ILDGLFN];
  QIO_String *filexml;
  char default_file_xml[] = "<?xml version=\"1.0\" encoding=\"UTF-8\"?><title>MILC ILDG archival gauge configuration</title>";

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

  xml_record_out = create_recxml(filename_milc);
  if(xml_record_out == NULL)return 1;

  if(this_node == 0)printf("Converting MILC v5 file %s to SciDAC file %s\n",
			   filename_milc, filename_scidac);
  if(ildgstyle == QIO_ILDGLAT)
    printf("in ILDG compatible format with LFN\n%s\n",ildg_lfn);

  initialize_machine(&argc,&argv);

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
  build_qio_filesystem(&fs);
  filexml = QIO_string_create();
  QIO_string_set(filexml, default_file_xml);
  outfile = open_scidac_output(filename_scidac, QIO_SINGLEFILE,
			       QIO_SERIAL, ildgstyle, ildg_lfn, &layout,
			       &fs, filexml);
  if(outfile == NULL)terminate(1);
  QIO_string_destroy(filexml);

  /* Initialize reading the MILC lattice data */
  r_serial_start_lattice(gf, &state);

  /* Write the SciDAC record. The factory function reads the
     site links from the MILC file */

  rec_info = QIO_create_record_info(QIO_FIELD, NULL, NULL, 0,
				    "USQCD_F3_ColorMatrix", "F", 
				    3, 0, datum_size, 4);
  status = QIO_write(outfile, rec_info, xml_record_out,
		     r_serial_site_links, datum_size*count, word_size, 
		     (void *)&state);
  if(status != QIO_SUCCESS)terminate(1);

  node0_printf("SciDAC checksums %x %x\n",
	       QIO_get_writer_last_checksuma(outfile),
	       QIO_get_writer_last_checksumb(outfile));

  /* cksum folds in the data length in bytes in the checksum */

  payload_bytes = (off_t)volume * count * datum_size;

  /* Take care of endianness before computing crc32 on length */
  if(big_endian()){
    if(sizeof(off_t) == 4)
      byterevn((int32type *)&payload_bytes, 1);
    else if(sizeof(off_t) == 8)
      byterevn64((int32type *)&payload_bytes, 1);
    else{
      printf("UNEXPECTED sizeof(off_t) = %d. Don't trust cksum!\n", 
	     sizeof(off_t));
    }
  }

  for(; payload_bytes != 0; payload_bytes >>= 8){
    unsigned char a = payload_bytes & 0xff;
    // printf("%02x", a);
    state.crc = crc32(state.crc,&a,1);
  }
  
  printf("cksum (crc32) checksum node %d %lu\n",this_node,state.crc);

  /* Close the SciDAC file */
  QIO_close_write(outfile);

  /* Finish the MILC v5 file */
  r_serial_finish_lattice(&state);

  r_serial_f(gf);

  QIO_string_destroy(xml_record_out);

  normal_exit(0);

  return 0;
}
