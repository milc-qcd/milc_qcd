/*************************** check_gauge.c ************************/
/* MIMD version 7 */
/* Read gauge configuration, check checksums (in version 5),
   and unitarity */

/* C. DeTar 10/30/97 */

/* Usage ...

   check_gauge gaugefile
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
#include "../include/generic.h"
#ifdef HAVE_QIO
#include <qio.h>
#endif

gauge_file *r_serial_i(char *filename);
void r_serial_f(gauge_file *gf);

#define PARALLEL 1
#define SERIAL 0

#define NODE_DUMP_ORDER 1
#define NATURAL_ORDER 0

#define MAX_BUF_LENGTH 4096

#define MAXERRCOUNT 10
#define TOLERANCE (0.0001)

/*----------------------------------------------------------------------*/
float ck_unitarity(su3_matrix *work,int x, int y, int z, int t)
{
  register int dir;
  register su3_matrix *mat;
  static int errcount = 0;
  int ii,jj;
  float deviation,max_deviation;
  union {
    float fval;
    int ival;
  } ifval;
  float check_su3(su3_matrix *c);
  
  mat = work;
  
  max_deviation=0;
  FORALLUPDIR(dir)
    {
      deviation = check_su3(&work[dir]);
      if (deviation>TOLERANCE){
	printf("Unitarity problem on node %d, site (%d,%d,%d,%d) dir %d, deviation=%f\n",
	       mynode(),x,y,z,t,dir,deviation);
	printf("SU3 matrix:\n");
	for(ii=0;ii<=2;ii++){
	  for(jj=0;jj<=2;jj++){
	    printf("%f ",mat[dir].e[ii][jj].real); 
	    printf("%f ",mat[dir].e[ii][jj].imag); 
	  }
	  printf("\n");
	}
	printf("repeat in hex:\n");
	for(ii=0;ii<=2;ii++){
	  for(jj=0;jj<=2;jj++){
	    ifval.fval = mat[dir].e[ii][jj].real; 
	    printf("%08x ", ifval.ival); 
	    ifval.fval = mat[dir].e[ii][jj].imag; 
	    printf("%08x ", ifval.ival); 
	  }
	  printf("\n");
	}
	printf("  \n \n");
	if(errcount++ >= MAXERRCOUNT)
	  {
	    fflush(stdout); terminate(1);
	  }
      }
      if(max_deviation<deviation) max_deviation=deviation;
    }
  return max_deviation;
}
/*----------------------------------------------------------------------*/

/* Here only node 0 reads the gauge configuration from a binary file */

void byterevn(int32type w[], int n);
void read_checksum(int parallel, gauge_file *gf, gauge_check *test_gc);

void r_check(gauge_file *gf, float *max_deviation)
{
  /* gf  = gauge configuration file structure */

  FILE *fp;
  gauge_header *gh;
  char *filename;
  int byterevflag;

  off_t offset ;            /* File stream pointer */
  off_t gauge_check_size;   /* Size of gauge configuration checksum record */
  off_t coord_list_size;    /* Size of coordinate list in bytes */
  off_t head_size;          /* Size of header plus coordinate list */
  off_t checksum_offset;    /* Where we put the checksum */
  int rcv_rank, rcv_coords;
  int destnode;
  int i,k;
  int x,y,z,t;
  int buf_length,where_in_buf;
  gauge_check test_gc;
  u_int32type *val;
  int rank29,rank31;
  su3_matrix *lbuf;
  su3_matrix work[4];
  float deviation;

  char myname[] = "r_check";

  fp = gf->fp;
  gh = gf->header;
  filename = gf->filename;
  byterevflag = gf->byterevflag;

  if(this_node == 0)
    {
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
      
      if(gf->header->order == NATURAL_ORDER)coord_list_size = 0;
      else coord_list_size = sizeof(int32type)*volume;
      checksum_offset = gf->header->header_bytes + coord_list_size;
      head_size = checksum_offset + gauge_check_size;
      
      /* Allocate space for read buffer */

      if(gf->parallel)
	printf("%s: Attempting serial read from parallel file \n",myname);

      lbuf = (su3_matrix *)malloc(MAX_BUF_LENGTH*4*sizeof(su3_matrix));
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
	  printf("%s: Node 0 g_seek %ld failed error %d file %s\n",
		 myname,(long)offset,errno,filename);
	  fflush(stdout);terminate(1);   
	}

      buf_length = 0;
      where_in_buf = 0;
      
    }

  /* all nodes initialize checksums */
  test_gc.sum29 = 0;
  test_gc.sum31 = 0;
  /* counts 32-bit words mod 29 and mod 31 in order of appearance
     on file */
  /* Here all nodes see the same sequence because we read serially */
  rank29 = 0;
  rank31 = 0;
  *max_deviation = 0;

  g_sync();

  /* Node 0 reads and deals out the values */

  for(rcv_rank=0; rcv_rank<volume; rcv_rank++)
    {
      /* If file is in coordinate natural order, receiving coordinate
         is given by rank. Otherwise, it is found in the table */

      if(gf->header->order == NATURAL_ORDER)
	rcv_coords = rcv_rank;
      else
	rcv_coords = gf->rank2rcv[rcv_rank];

      x = rcv_coords % nx;   rcv_coords /= nx;
      y = rcv_coords % ny;   rcv_coords /= ny;
      z = rcv_coords % nz;   rcv_coords /= nz;
      t = rcv_coords % nt;

      /* The node that gets the next set of gauge links */
      destnode=node_number(x,y,z,t);
      
      if(this_node==0){
	/* Node 0 fills its buffer, if necessary */
	if(where_in_buf == buf_length)
	  {  /* get new buffer */
	    /* new buffer length  = remaining sites, but never bigger 
	       than MAX_BUF_LENGTH */
	    buf_length = volume - rcv_rank;
	    if(buf_length > MAX_BUF_LENGTH)buf_length = MAX_BUF_LENGTH;
	    /* then do read */
	    
	    if( (int)g_read(lbuf,4*sizeof(su3_matrix),buf_length,fp) != buf_length)
	      {
		printf("%s: node %d gauge configuration read error %d file %s\n",
		       myname,this_node,errno,filename); 
		fflush(stdout); terminate(1);
	      }
	    where_in_buf = 0;  /* reset counter */
	  }  /*** end of the buffer read ****/

	if(destnode==0){	/* just copy links */
	  i = node_index(x,y,z,t);
	  memcpy((void *)&work[0],
		 (void *)&lbuf[4*where_in_buf], 4*sizeof(su3_matrix));
	}
	else {		/* send to correct node */
	  send_field((char *)&lbuf[4*where_in_buf],
		     4*sizeof(su3_matrix),destnode);
	}
	where_in_buf++;
      }
      
      /* The node which contains this site reads message */
      else {	/* for all nodes other than node 0 */
	if(this_node==destnode){
	  i = node_index(x,y,z,t);
	  get_field((char *)&work[0],4*sizeof(su3_matrix),0);
	}
      }

      /* The receiving node does the byte reversal and then checksum,
         if needed */

      if(this_node==destnode)
	{
	  if(byterevflag==1)
	    byterevn((int32type *)&work[0],
		     4*sizeof(su3_matrix)/sizeof(int32type));
	  /* Accumulate checksums */
	  for(k = 0, val = (u_int32type *)&work[0]; 
	      k < 4*(int)sizeof(su3_matrix)/(int)sizeof(int32type); k++, val++)
	    {
	      test_gc.sum29 ^= (*val)<<rank29 | (*val)>>(32-rank29);
	      test_gc.sum31 ^= (*val)<<rank31 | (*val)>>(32-rank31);
	      rank29++; if(rank29 >= 29)rank29 = 0;
	      rank31++; if(rank31 >= 31)rank31 = 0;
	    }
	  deviation = ck_unitarity(work,x,y,z,t);
	  if(deviation > *max_deviation)*max_deviation = deviation;
	}
      else
	{
	  rank29 += 4*sizeof(su3_matrix)/sizeof(int32type);
	  rank31 += 4*sizeof(su3_matrix)/sizeof(int32type);
	  rank29 %= 29;
	  rank31 %= 31;
	}
    }
  
  /* Combine node checksum contributions with global exclusive or */
  g_xor32(&test_gc.sum29);
  g_xor32(&test_gc.sum31);
  
  if(this_node==0)
    {
      /* Read and verify checksum */
      /* Checksums not implemented until version 5 */
      
      printf("Restored binary gauge configuration serially from file %s\n",
	     filename);
      if(gh->magic_number == GAUGE_VERSION_NUMBER)
	{
	  if( g_seek(fp,checksum_offset,SEEK_SET) < 0 ) 
	    {
	      printf("%s: Node 0 g_seek %ld failed error %d file %s\n",
		    myname,(long)offset,errno,filename);
	      fflush(stdout);terminate(1);   
	    }
	  read_checksum(SERIAL,gf,&test_gc);
	}
      else
	{
	  printf("Checksums %x %x\n",test_gc.sum29,test_gc.sum31);
	  printf("Checksums not verified in this format\n");
	}
      fflush(stdout);
      free(lbuf);
    }
  
} /* r_check */

/*--------------------------------------------------------------*/

void r_check_arch(gauge_file *gf, float *max_deviation)
{
  /* gf  = gauge configuration file structure */

  FILE *fp;
  gauge_header *gh;
  char *filename;
  int byterevflag;

  off_t gauge_check_size;   /* Size of gauge configuration checksum record */
  int rcv_rank, rcv_coords;
  int destnode;
  int i,k;
  int x,y,z,t;
  gauge_check test_gc;
  u_int32type *val;
  int rank29,rank31;
  su3_matrix work[4];
  float deviation;
  su3_matrix tmpsu3[4];
  int dataformat = gf->dataformat;
  int precision = gf->precision;
  char myname[] = "r_check_arch";
  int mu,a,b,p;
  float *uin = NULL, *q;
  double *uind = NULL, *qd;
  int big_end = 0;
  float U[4][18];
  double Ud[4][18];
  u_int32type chksum;
  int realspersite;


  if(dataformat == ARCHIVE_3x2)realspersite = 48;
  else realspersite = 72;
  fp = gf->fp;
  gh = gf->header;
  filename = gf->filename;
  byterevflag = gf->byterevflag;

  if(this_node == 0)
    {
      gauge_check_size = 0;
      
      if(gf->parallel)
	printf("%s: Attempting serial read from parallel file \n",myname);

      big_end = big_endian();
      /* printf("big_end is %d\n", big_end); */
      if(precision == 1){
	uin = (float *) malloc(nx*ny*nz*realspersite*sizeof(float));
	if(uin == NULL)
	  {
	    printf("%s: Node %d can't malloc uin buffer to read timeslice\n",
		   myname,this_node);
	    printf("recompile with smaller read buffer: uin\n");
	    fflush(stdout);
	    terminate(1);
	  }
      } else {
	uind = (double *) malloc(nx*ny*nz*realspersite*sizeof(double));
	if(uind == NULL)
	  {
	    printf("%s: Node %d can't malloc uind buffer to read timeslice\n",
		   myname,this_node);
	    printf("recompile with smaller read buffers: uin\n");
	    fflush(stdout);
	    terminate(1);
	  }
      }
    }
      
  /* Initialize checksums */
  chksum = 0;
  test_gc.sum29 = 0;
  test_gc.sum31 = 0;
  /* counts 32-bit words mod 29 and mod 31 in order of appearance
     on file */
  /* Here all nodes see the same sequence because we read serially */
  rank29 = 0;
  rank31 = 0;
  *max_deviation = 0;

  g_sync();

  /* Node 0 reads and deals out the values */
  for(rcv_rank=0; rcv_rank<volume; rcv_rank++)
    {
      rcv_coords = rcv_rank;

      x = rcv_coords % nx;   rcv_coords /= nx;
      y = rcv_coords % ny;   rcv_coords /= ny;
      z = rcv_coords % nz;   rcv_coords /= nz;
      t = rcv_coords % nt;

      /* The node that gets the next set of gauge links */
      destnode=node_number(x,y,z,t);
      
      if(this_node==0){
	if(precision == 1){
	  if( (int)g_read(uin,48*sizeof(float),1,fp) != 1)
	    {
	      printf("%s: node %d gauge configuration read error %d file %s\n",
		     myname,this_node,errno,filename); 
	      fflush(stdout); terminate(1);
	    }
	  if (!big_end) byterevn((int32type *)uin,48);
	  q = uin;
	  for (mu=0;mu<4;mu++) {
	    for (p=0;p<realspersite/4;p++) {
	      chksum += *(u_int32type *) q;
	      U[mu][p] = (float) *(q++);
	    }
	    if(dataformat == ARCHIVE_3x2)
	      complete_U(U[mu]);
	    /**
	       for (p=0;p<18;p++) printf("p=%d, e=%f\n", p, U[mu][p]);
	    **/
	    
	    /* Copy, converting precision if necessary */
	    for(a=0; a<3; a++) for(b=0; b<3; b++) { 
	      tmpsu3[mu].e[a][b].real = U[mu][2*(3*a+b)];
	      /* printf("real: p=%d, mu=%d, e=%f\n", p,mu,U[mu][2*(3*a+b)]); */
	      tmpsu3[mu].e[a][b].imag = U[mu][2*(3*a+b)+1];
	      /*printf("imag: p=%d, mu=%d, e=%f\n", p,mu,U[mu][2*(3*a+b)+1]); */
	    } 
	  }
	} else { /* precision == 2 */
	  
	  if( (int)g_read(uind,realspersite*sizeof(double),1,fp) != 1)
	    {
	      printf("%s: node %d gauge configuration read error %d file %s\n",
		     myname,this_node,errno,filename); 
	      fflush(stdout); terminate(1);
	    }
	  if (!big_end) byterevn64((int32type *)uind,realspersite);
	  qd = uind;
	  for (mu=0;mu<4;mu++) {
	    for (p=0;p<realspersite/4;p++) {
	      chksum += *(u_int32type *) qd;
	      chksum += *((u_int32type *) qd + 1);
	      Ud[mu][p] = (double) *(qd++);
	    }
	    if(dataformat == ARCHIVE_3x2)
	      complete_Ud(Ud[mu]);
	    /**
	       for (p=0;p<18;p++) printf("p=%d, e=%f\n", p, Ud[mu][p]);
	    **/
	    
	    /* Copy, converting precision if necessary */
	    for(a=0; a<3; a++) for(b=0; b<3; b++) { 
	      tmpsu3[mu].e[a][b].real = Ud[mu][2*(3*a+b)];
	      /* printf("real: p=%d, mu=%d, e=%f\n", p,mu,Ud[mu][2*(3*a+b)]); */
	      tmpsu3[mu].e[a][b].imag = Ud[mu][2*(3*a+b)+1];
	      /*printf("imag: p=%d, mu=%d, e=%f\n", p,mu,Ud[mu][2*(3*a+b)+1]); */
	    }
	  } 
	}
	
	if(destnode==0){	/* just copy links */
	  i = node_index(x,y,z,t);
     /*   printf("lattice node_index = %d, mu = %d\n", i, mu); */
	  memcpy((void *)&work[0],
		 (void *)&tmpsu3, 4*sizeof(su3_matrix));
	} else {		/* send to correct node */
	  send_field((char *)tmpsu3, 4*sizeof(su3_matrix),destnode);
	}
      } 
      /* The node which contains this site reads message */
      else {	/* for all nodes other than node 0 */
	if(this_node==destnode){
	  i = node_index(x,y,z,t);
	  get_field((char *)&work[0],4*sizeof(su3_matrix),0);
	}
      }

      /* Any needed byte reversing was already done. Compute MILC checksums */

      if(this_node==destnode)
	{
	  /* Accumulate checksums */
	  for(k = 0, val = (u_int32type *)&work[0]; 
	      k < 4*(int)sizeof(su3_matrix)/(int)sizeof(int32type); k++, val++)
   	    {
	      test_gc.sum29 ^= (*val)<<rank29 | (*val)>>(32-rank29);
	      test_gc.sum31 ^= (*val)<<rank31 | (*val)>>(32-rank31);
	      rank29++; if(rank29 >= 29)rank29 = 0;
	      rank31++; if(rank31 >= 31)rank31 = 0;
	    }
	  deviation = ck_unitarity(work,x,y,z,t);
	  if(deviation > *max_deviation)*max_deviation = deviation;
	}
      else
	{
	  rank29 += 4*sizeof(su3_matrix)/sizeof(int32type);
	  rank31 += 4*sizeof(su3_matrix)/sizeof(int32type);
	  rank29 %= 29;
	  rank31 %= 31;
	}
    }
  
  /* Combine node checksum contributions with global exclusive or */
  g_xor32(&test_gc.sum29);
  g_xor32(&test_gc.sum31);
  
  if(this_node==0)
    {
      /* Read and verify checksum */
      
      printf("Restored archive gauge configuration serially from file %s\n",
	     filename);
      if (chksum != gf->check.sum31)
	{
	  printf("Archive style checksum violation: computed %x, read %x\n",
		 chksum, gf->check.sum31);
	}
      else
	{
	  printf("Archive style checksum = %x OK\n", chksum);
	}
      fflush(stdout);
      free(uin);

      /* Store MILC style checksums */
      gf->check.sum29 = test_gc.sum29;
      gf->check.sum31 = test_gc.sum31;
    }
  
} /* r_check_arch */

/*----------------------------------------------------------------------*/

int main(int argc, char *argv[])
{

  gauge_file *gf;
  gauge_header *gh;
  char *filename;
  float max_deviation;
  
  if(argc < 2)
    {
      fprintf(stderr,"Usage %s <gaugefilename>\n",argv[0]);
      exit(1);
    }
  filename = argv[1];

  initialize_machine(&argc,&argv);

  this_node = mynode();
  number_of_nodes = numnodes();

  if(this_node == 0)printf("Checking file %s\n",filename);

  /* Read header */
  nx = ny = nz = nt = -1;  /* To suppress dimension checking */
  gf = r_serial_i(filename);
  gh = gf->header;

  if(gh->magic_number == LIME_MAGIC_NO){
    printf("Can't check a SciDAC file yet\n");
    return 1;
  }

  if(this_node == 0)
    {
      nx = gh->dims[0];
      ny = gh->dims[1];
      nz = gh->dims[2];
      nt = gh->dims[3];
      printf("Dimensions %d %d %d %d\n",nx,ny,nz,nt);
      if(gh->magic_number != GAUGE_VERSION_NUMBER_ARCHIVE){
	printf("Time stamp %s\n",gh->time_stamp);
	if(gh->order == NATURAL_ORDER)printf("File in natural order\n");
	if(gh->order == NODE_DUMP_ORDER)printf("File in node dump order\n");
      }
    }

  broadcast_bytes((char *)&nx,sizeof(int));
  broadcast_bytes((char *)&ny,sizeof(int));
  broadcast_bytes((char *)&nz,sizeof(int));
  broadcast_bytes((char *)&nt,sizeof(int));

  setup_layout();

  volume = nx*ny*nz*nt;

  if(gh->magic_number == GAUGE_VERSION_NUMBER_ARCHIVE){
    r_check_arch(gf,&max_deviation);
    g_floatmax(&max_deviation);
    if(this_node==0)
      printf("Max unitarity deviation = %0.2g\n",max_deviation);
  }
  else{
    r_check(gf,&max_deviation);
    g_floatmax(&max_deviation);
    if(this_node==0)
      printf("Max unitarity deviation = %0.2g\n",max_deviation);
  }


  
  /* Close file */
  r_serial_f(gf);

  normal_exit(0);

  return 0;
}
