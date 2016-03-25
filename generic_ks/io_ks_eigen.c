/* io_ks_eigen.c -- reads and writes KS eigenvectors
   MIMD version 7

   H. Ohno: 11/20/2014 -- derived from io_prop_ks.c
*/
/* This version assumes internal storage is at the prevailing
   precision, but the files are always 32 bit.  This code
   converts to and from the prevailing precision.  CD 11/29/04 */

#include "generic_ks_includes.h"
#include "../include/io_ks_eigen.h"
#include <sys/types.h>
#include <fcntl.h>
#include <errno.h>
#include <string.h>
#include <time.h>
#include <assert.h>

#define PARALLEL 1
#define SERIAL 0

#define NATURAL_ORDER 0

#undef MAX_BUF_LENGTH
#define MAX_BUF_LENGTH 4096

/* The send buffer would normally be the size of a single SU(3) vector
   but we may need to pad for short messages generated in serial reads
   and writes to avoid switch inefficiencies.  In that case we define
   PAD_SEND_BUF on the compilation line to increase this.  */
#ifndef PAD_SEND_BUF
#define PAD_SEND_BUF 8
#endif

/*----------------------------------------------------------------------*/

/* This subroutine writes the KS eigenvector header structure */
/* Serial access version */
void swrite_ks_eigen_hdr(FILE *fp, ks_eigen_header *kseigh){

  char myname[] = "swrite_ks_eigen_hdr";

  //  swrite_data(fp, (void *)&kseigh->magic_number, sizeof(kseigh->magic_number),
  //              myname, "magic_number");
  swrite_data(fp, kseigh->dims, sizeof(kseigh->dims),
              myname, "dimensions");
  swrite_data(fp, kseigh->time_stamp, sizeof(kseigh->time_stamp),
              myname, "time_stamp");
  swrite_data(fp, &kseigh->order, sizeof(kseigh->order),
              myname, "order");

  /* Header byte length */
  kseigh->header_bytes = //sizeof(ksegih->magic_number) +
    sizeof(kseigh->dims) + sizeof(kseigh->time_stamp) + sizeof(kseigh->order);

} /* swrite_ks_eigen_hdr */

/*------------------------------------------------------------------------*/

/* parallel = 1 (TRUE) if all nodes are accessing the file */
/*            0        for access from node 0 only */
/* Returns byterevflag  = 0 or 1 */
int read_ks_eigen_hdr(ks_eigen_file *kseigf, int parallel){

  FILE *fp;
  ks_eigen_header *kseigh;
  //int32type tmp;
  int j;
  int byterevflag = 0;
  char myname[] = "read_ks_eigen_hdr";

  fp = kseigf->fp;
  kseigh = kseigf->header;

  /* Read and verify magic number */
  //if(psread_data(parallel, fp, &kseigh->magic_number,
  //		 sizeof(kseigh->magic_number), myname, "magic number") != 0)
  //terminate(1);

  //tmp = kseigh->magic_number;
  //if(kseigh->magic_number == KS_EIGEN_VERSION_NUMBER) 
  //  byterevflag = 0;
  //else{
  //  byterevn((int32type *)&kseigh->magic_number, 1);
  //  if(kseigh->magic_number == KS_EIGEN_VERSION_NUMBER){
  //    byterevflag = 1;
      /** printf("Reading with byte reversal\n"); **/
  //    if(sizeof(Real) != sizeof(int32type)){
  //	  printf("%s: Can't byte reverse\n", myname);
  //	  printf("requires size of int32type(%d) = size of Real(%d)\n",
  //	         (int)sizeof(int32type), (int)sizeof(Real));
  //  	  terminate(1);
  //    }
  //  }
  //  else{
      /* Restore magic number as originally read */
  //    kseigh->magic_number = tmp;
	  
      /* End of the road. */
  //    printf("%s: Unrecognized magic number in KS eigenvalue file header.\n",
  //	     myname);
  //    printf("Expected %x but read %x\n", KS_EIGEN_VERSION_NUMBER, tmp);
  //    terminate(1);
  //  }
  //}
  
  /* Read header, do byte reversal, 
     if necessary, and check consistency */
  
  /* Lattice dimensions */
  if(psread_byteorder(byterevflag, parallel, fp, kseigh->dims,
		      sizeof(kseigh->dims), myname, "dimensions") != 0)
    terminate(1);

  if(kseigh->dims[0] != nx || kseigh->dims[1] != ny ||
     kseigh->dims[2] != nz || kseigh->dims[3] != nt){
    /* So we can use this routine to discover the dimensions,
       we provide that if nx = ny = nz = nt = -1 initially
       we don't die */
    if(nx != -1 || ny != -1 || nz != -1 || nt != -1){
      printf("%s: Incorrect lattice dimensions ", myname);
      for(j = 0; j < 4; j++) printf("%d ", kseigh->dims[j]); 
      printf("\n");
      fflush(stdout);
      terminate(1);
    }
    else{
      nx = kseigh->dims[0];
      ny = kseigh->dims[1];
      nz = kseigh->dims[2];
      nt = kseigh->dims[3];
      volume = nx*ny*nz*nt;
    }
  }
  
  /* Date and time stamp */
  if(psread_data(parallel, fp, kseigh->time_stamp,
		 sizeof(kseigh->time_stamp), myname, "time stamp") != 0)
    terminate(1);

  /* Header byte length */
  kseigh->header_bytes = //sizeof(kseigh->magic_number) +
    sizeof(kseigh->dims) + sizeof(kseigh->time_stamp) + sizeof(kseigh->order);
  
  /* Data order */  
  if(psread_byteorder(byterevflag, parallel, fp, &kseigh->order,
		      sizeof(kseigh->order), myname, "order parameter") != 0)
    terminate(1);
  
  return byterevflag;
  
} /* read_ks_eigen_hdr */

/*------------------------------------------------------------------------*/

/* Write a data item to the KS eigenvector info file */
int write_ks_eigen_info_item(FILE *fpout,    /* ascii file pointer */
			     char *keyword,  /* keyword */
			     char *fmt,      /* output format -
						must use s, d, e, f, or g */
			     char *src,      /* address of starting data
						floating point data must be
						of type (Real) */
			     int count,      /* number of data items if > 1 */
			     int stride){    /* byte stride of data if
						count > 1 */

  int i, n;
  char *data;
#define MAXITEM 128
  char myname[] = "write_ks_eigen_info_item";

  /* Check for valid keyword */
  for(i=0; strlen(ks_eigen_info_keyword[i])>0 &&
	strcmp(ks_eigen_info_keyword[i], keyword) != 0; i++);

  /* Write keyword */
  fprintf(fpout, "%s", keyword);

  /* Write count if more than one item */
  if(count > 1) fprintf(fpout, "[%d]", count);

  n = count;
  if(n==0) n = 1;
  
  /* Write data */
  for(i = 0, data = (char *)src; i < n; i++, data += stride){
    fprintf(fpout, " ");
    if(strstr(fmt, "s") != NULL)
      fprintf(fpout, fmt, data);
    else if(strstr(fmt, "d") != NULL)
      fprintf(fpout, fmt, *(int *)data);
    else if(strstr(fmt, "e") != NULL)
      fprintf(fpout, fmt, (double)(*(Real *)data));
    else if(strstr(fmt, "f") != NULL)
      fprintf(fpout, fmt, (double)(*(Real *)data));
    else if(strstr(fmt, "g") != NULL)
      fprintf(fpout, fmt, (double)(*(Real *)data));
    else{
      printf("%s: Unrecognized data type %s\n", myname, fmt);
      return 1;
    }
  }
  fprintf(fpout, "\n");

  return 0;
} /* write_ks_eigen_info_item */

/*----------------------------------------------------------------------*/
/* Open, write, and close the ASCII info file */
void write_ks_eigen_info_file(ks_eigen_file *kseigf){

  int i;
  FILE *info_fp;
  ks_eigen_header *kseigh;
  char info_filename[256];
  char sbuffer[64];
  char myname[] = "write_ks_eigen_info_file";

  kseigh = kseigf->header;

  /* Construct metadata file name from KS eigenvector file name 
     by adding filename extension to KS eigenvector file name */
  strcpy(info_filename, kseigf->filename);
  strcat(info_filename, ASCII_INFO_EXT);

  /* Open metadata file */
  if((info_fp = fopen(info_filename,"w")) == NULL){
    printf("%s: Can't open ascii info file %s\n", myname, info_filename);
    return;
  }
  
  /* Write required information */
  //  write_ks_eigen_info_item(info_fp, "magic_number", "%d",
  //			   (char *)&kseigh->magic_number, 0, 0);
  write_ks_eigen_info_item(info_fp, "time_stamp", "\"%s\"", kseigh->time_stamp, 0, 0);
  sprintf(sbuffer, "%x %x", kseigf->check.sum29, kseigf->check.sum31);
  write_ks_eigen_info_item(info_fp, "checksums", "\"%s\"", sbuffer, 0, 0);
  write_ks_eigen_info_item(info_fp, "nx","%d", (char *)&nx, 0, 0);
  write_ks_eigen_info_item(info_fp, "ny","%d", (char *)&ny, 0, 0);
  write_ks_eigen_info_item(info_fp, "nz","%d", (char *)&nz, 0, 0);
  write_ks_eigen_info_item(info_fp, "nt", "%d", (char *)&nt, 0 ,0);
  write_ks_eigen_info_item(info_fp, "parity", "%d", (char *)&(kseigf->parity), 0 ,0);
  write_ks_eigen_info_item(info_fp, "Nvecs", "%d", (char *)&(kseigf->Nvecs), 0 ,0);
  for(i = 0; i < kseigf->Nvecs; i++){
    sprintf(sbuffer,"[%d] %e  (resid = %e)", i, kseigf->eigVal[i], kseigf->resid[i]);
    write_ks_eigen_info_item(info_fp, "eigVal", "%s", sbuffer, 0, 0);
  }

  fclose(info_fp);

  printf("Wrote info file %s\n", info_filename); 
  fflush(stdout);

} /*write_ks_eigen_info_file */

/*----------------------------------------------------------------------*/

/* Set up the input KS eigenvector file and header structures */
ks_eigen_file *create_input_ks_eigen_file_handle(char *filename){

  ks_eigen_file *kseigf;
  ks_eigen_header *kseigh;
  char myname[] = "create_input_ks_eigen_file_handle";

  /* Allocate space for the file structure */
  kseigf = (ks_eigen_file *)malloc(sizeof(ks_eigen_file));
  if(kseigf == NULL){
    printf("%s: Can't malloc KS eigenvector file\n", myname);
    terminate(1);
  }

  kseigf->filename = filename;
  kseigf->fp = NULL;

  /* Allocate space for the header */
  /* Make sure compilation gave us a 32 bit integer type */
  assert(sizeof(int32type) == 4);

  kseigh = (ks_eigen_header *)malloc(sizeof(ks_eigen_header));
  if(kseigh == NULL){
    printf("%s: Can't malloc KS eigenvector header\n", myname);
    terminate(1);
  }

  kseigf->header = kseigh;
  kseigf->Nvecs = 0;
  kseigf->eigVal = NULL;
  kseigf->resid = NULL;
  //kseigf->info = NULL;
  kseigf->byterevflag = 0;
  kseigf->parallel = 0;
  kseigf->info_fp = NULL;

  return kseigf;
} /* create_input_ks_eigen_file_handle */

/*----------------------------------------------------------------------*/

/* Set up the output KS eigenvector file and  header structure */
ks_eigen_file *create_output_ks_eigen_file_handle(void){

  ks_eigen_file *kseigf;
  ks_eigen_header *kseigh;
  time_t time_stamp;
  int i;
  char myname[] = "create_output_ks_eigen_file_handle";

  /* Allocate space for a new file structure */
  kseigf = (ks_eigen_file *)malloc(sizeof(ks_eigen_file));
  if(kseigf == NULL){
    printf("%s: Can't malloc KS eigenvector file\n", myname);
    terminate(1);
  }
  
  /* Allocate space for a new header structure */
  /* Make sure compilation gave us a 32 bit integer type */
  assert(sizeof(int32type) == 4);

  kseigh = (ks_eigen_header *)malloc(sizeof(ks_eigen_header));
  if(kseigh == NULL){
    printf("%s: Can't malloc KS eigenvector header\n", myname);
    terminate(1);
  }

  /* Load header pointer and file name */
  kseigf->header = kseigh;

  /* Initialize */
  kseigf->check.sum29 = 0;
  kseigf->check.sum31 = 0;
  kseigf->Nvecs = 0;
  kseigf->eigVal = NULL;
  kseigf->resid = NULL;
  //kseigf->info = NULL;

  /* Load header values */
  //kseigh->magic_number = KS_EIGEN_VERSION_NUMBER;
  kseigh->dims[0] = nx;
  kseigh->dims[1] = ny;
  kseigh->dims[2] = nz;
  kseigh->dims[3] = nt;

  /* Get date and time stamp. (We use local time on node 0) */
  if(this_node==0){
    time(&time_stamp);
    strcpy(kseigh->time_stamp, ctime(&time_stamp));
    /* For aesthetic reasons, don't leave trailing junk bytes here to be
       written to the file */
    for(i = strlen(kseigh->time_stamp) + 1;
	i < (int)sizeof(kseigh->time_stamp); i++)
      kseigh->time_stamp[i] = '\0';
      
    /* Remove trailing end-of-line character */
    if(kseigh->time_stamp[strlen(kseigh->time_stamp) - 1] == '\n')
      kseigh->time_stamp[strlen(kseigh->time_stamp) - 1] = '\0';
  }
  
  /* Broadcast to all nodes */
  broadcast_bytes(kseigh->time_stamp, sizeof(kseigh->time_stamp));

  return kseigf;
} /* setup_output_ks_eigen_file_handle */

/* Destroy the KS eigenvector file and  header structure */
void destroy_ks_eigen_file_handle(ks_eigen_file *kseigf){

  if(kseigf == NULL) return;

  if(kseigf->header != NULL) free(kseigf->header);
  //if(kseigf->info != NULL) free(kseigf->info);

  free(kseigf);
}

/*----------------------------------------------------------------------*/

/* Open a binary file for serial writing by node 0 */
/* Only node 0 opens the file filename */
/* Returns a file structure describing the opened file */
ks_eigen_file *w_serial_ks_eigen_i(char *filename, int parity){

  FILE *fp;
  ks_eigen_file *kseigf;
  ks_eigen_header *kseigh;
  char myname[] = "w_serial_ks_eigen_i";

  /* Set up KS eigenvector file and header structs and load header values */
  kseigf = create_output_ks_eigen_file_handle();
  kseigh = kseigf->header;

  /* Indicate coordinate natural ordering */
  kseigh->order = NATURAL_ORDER;

  /* Only node 0 opens the requested file */
  if(this_node == 0){
    fp = fopen(filename, "wb");
    if(fp == NULL){
      printf("%s: Node %d can't open file %s, error %d\n",
	     myname, this_node, filename, errno);
      fflush(stdout);
      terminate(1);
    }

    printf("Opened KS eigenvector file %s for serial writing\n", filename);

    /* Node 0 writes the header */
    swrite_ks_eigen_hdr(fp, kseigh);
  }
  
  /* Assign values to file structure */
  if(this_node==0) kseigf->fp = fp; 
  else kseigf->fp = NULL;       /* Only node 0 knows about this file */

  kseigf->filename = filename;
  kseigf->byterevflag = 0;      /* Not used for writing */
  kseigf->parallel = SERIAL;
  kseigf->parity = parity;

  return kseigf;
} /* w_serial_ks_eigen_i */

/*---------------------------------------------------------------------------*/

/* Here only node 0 writes eigenvectors to a serial file */
void w_serial_ks_eigen(ks_eigen_file *kseigf, int Nvecs, double *eigVal, su3_vector **eigVec,
		       double *resid){

  FILE *fp = NULL;
  u_int32type *val;
  int rank29, rank31;
  su3_vector *eigbuf = NULL;
  struct {
    su3_vector ksv;
    char pad[PAD_SEND_BUF]; /* Introduced because some switches
			       perform better if message lengths are longer */
  } msg;
  int buf_length = 0;

  int i, j, k = 0;
  int parity, ndata = 0;
  int currentnode, newnode;
  int x, y, z, t, ivecs;
  char myname[] = "w_serial_ks_eigen";

  if(this_node == 0){
    if(kseigf->parallel == PARALLEL)
      printf("%s: Attempting serial write to file opened in parallel \n", myname);

    kseigf->Nvecs = Nvecs;
    kseigf->eigVal = eigVal;
    kseigf->resid = resid;
    
    /* Node 0 writes ascii info file */
    write_ks_eigen_info_file(kseigf);
    
    eigbuf = (su3_vector *)malloc(MAX_BUF_LENGTH*sizeof(su3_vector));
    if(eigbuf == NULL){
      printf("%s: Node 0 can't malloc eigbuf\n", myname); 
      fflush(stdout);
      terminate(1);
    }

    fp = kseigf->fp;

    /* Write eigenvalues */
    swrite_data(fp, &(kseigf->parity), sizeof(int), myname, "parity");
    swrite_data(fp, &Nvecs, sizeof(int), myname, "Nvecs");
    if((int)fwrite(eigVal, sizeof(double), Nvecs, fp) != Nvecs){
      printf("%s: Node %d eigenvalue write error %d file %s\n",
	     myname, this_node, errno, kseigf->filename); 
      fflush(stdout);
      terminate(1);
    }

    ndata = (kseigf->parity == EVEN || kseigf->parity == ODD) ? volume/2 : volume;
    ndata *= Nvecs;
  } /* end if(this_node==0) */
  
  /* initialize checksums */
  kseigf->check.sum31 = 0;
  kseigf->check.sum29 = 0;
  /* counts 32-bit words mod 29 and mod 31 in order of appearance on file */
  /* Here only node 0 uses these values */
  rank29 = sizeof(su3_vector)/sizeof(int32type)*sites_on_node*this_node % 29;
  rank31 = sizeof(su3_vector)/sizeof(int32type)*sites_on_node*this_node % 31;

  currentnode = 0;

  /* Write eigenvectors */   
  for(ivecs = 0; ivecs < Nvecs; ivecs++)
  for(t = 0; t < nt; t++)
  for(z = 0; z < nz; z++)
  for(y = 0; y < ny; y++)
  for(x = 0; x < nx; x++){
    parity = 2 - (x + y + z + t)%2;
    if(parity == kseigf->parity || kseigf->parity == EVENANDODD){
      newnode = node_number(x, y, z, t);
      if(newnode != currentnode){	/* switch to another node */
	/* Send a few bytes of garbage to tell newnode it's OK to send */
	if(this_node == 0 && newnode != 0)
	  send_field((char *)&msg, sizeof(msg), newnode);
	if(this_node == newnode && newnode != 0)
	  get_field((char *)&msg, sizeof(msg), 0);
	currentnode = newnode;
      }
      
      if(this_node == 0){
	if(currentnode == 0){
	  i = node_index(x, y, z, t);
	  msg.ksv = eigVec[ivecs][i];
	}
	else
	  get_field((char *)&msg, sizeof(msg), currentnode);

	eigbuf[buf_length] = msg.ksv;

	/* Accumulate checksums - contribution from next site */
	for(j = 0, val = (u_int32type *)&eigbuf[buf_length]; 
	    j < (int)sizeof(su3_vector)/(int)sizeof(int32type); j++, val++){
	  kseigf->check.sum29 ^= (*val)<<rank29 | (*val)>>(32-rank29);
	  kseigf->check.sum31 ^= (*val)<<rank31 | (*val)>>(32-rank31);
	  rank29++; if(rank29 >= 29) rank29 = 0;
	  rank31++; if(rank31 >= 31) rank31 = 0;
	}

	k++;
	buf_length++;
	  
	if((buf_length == MAX_BUF_LENGTH) || (k == ndata)){
	  /* write out buffer */
	  if((int)fwrite(eigbuf, sizeof(su3_vector), buf_length, fp) != buf_length){
	    printf("%s: Node %d eigenvector write error %d file %s\n",
		   myname, this_node, errno, kseigf->filename); 
	    fflush(stdout);
	    terminate(1);
	  }
	  buf_length = 0;		/* start again after write */
	}
      }
      else{  /* for nodes other than 0 */
	if(this_node == currentnode){
	  i = node_index(x, y, z, t);
	  msg.ksv = eigVec[ivecs][i];
	  send_field((char *)&msg, sizeof(msg), 0);
	}
      }
    } /* if(parity == kseigf->parity || kseigf->parity == EVENANDODD) */
  } /* close x, y, z, t, ivecs loops */
  
  g_sync();

  if(this_node==0){
    printf("Wrote eigenvectors serially to file %s\n", kseigf->filename); 
    fflush(stdout);

    free(eigbuf);

    /* Construct check record */
    pswrite_data(kseigf->parallel, fp, &kseigf->check.sum29,
		 sizeof(kseigf->check.sum29), myname, "checksum");
    pswrite_data(kseigf->parallel, fp, &kseigf->check.sum31,
		 sizeof(kseigf->check.sum31), myname, "checksum");
  }
} /* w_serial_ks_eigen */

/*---------------------------------------------------------------------------*/

/* This subroutine closes the file and frees associated structures */
void w_serial_ks_eigen_f(ks_eigen_file *kseigf){

  char myname[] = "w_serial_ks_eigen_f";

  g_sync();
  
  if(this_node==0){
    if(kseigf->parallel == PARALLEL)
      printf("%s: Attempting serial close on file opened in parallel\n",
	     myname);
    printf("Closed KS eigenvector file %s time stamp %s\n", kseigf->filename,
	   (kseigf->header)->time_stamp);

    if(kseigf->fp != NULL) fclose(kseigf->fp);
  }

  /* Free header and file structures */
  destroy_ks_eigen_file_handle(kseigf);

} /* w_serial_ks_eigen_f */

/*---------------------------------------------------------------------------*/

/* Returns file descriptor for opened file */
ks_eigen_file *r_serial_ks_eigen_i(char *filename){

  FILE *fp;
  ks_eigen_file *kseigf;
  ks_eigen_header *kseigh;
  int byterevflag;
  char myname[] = "r_serial_ks_eigen_i";

  /* All nodes set up an eigenvector file and header structure for reading */
  kseigf = create_input_ks_eigen_file_handle(filename);
  kseigh = kseigf->header;

  /* File opened for serial reading */
  kseigf->parallel = 0;

  /* Node 0 alone opens a file and reads the header */
  if(this_node==0){
    fp = fopen(filename, "rb");
    if(fp == NULL){
      printf("%s: Node %d can't open file %s, error %d\n",
	     myname, this_node, filename, errno);
      fflush(stdout);
      terminate(1);
    }
      
    printf("Opened KS eigenvector file %s for serial reading\n", filename);
      
    kseigf->fp = fp;

    byterevflag = read_ks_eigen_hdr(kseigf, SERIAL);
  }

  /* Read the metadata file */
  //kseigf->info = read_info_file(filename);

  /* Node 0 broadcasts the byterevflag from node 0 to all nodes */
  broadcast_bytes((char *)&byterevflag, sizeof(byterevflag));
  kseigf->byterevflag = byterevflag;
  
  /* Node 0 broadcasts the header structure to all nodes */
  broadcast_bytes((char *)kseigh, sizeof(ks_eigen_header));

  return kseigf;
}/* r_serial_ks_eigen_i */

/*----------------------------------------------------------------------*/

/* Here only node 0 reads the KS eigenvectors from a binary file */
/* 0 is normal exit code
   >1 for seek, read error, or missing data error */
int r_serial_ks_eigen(ks_eigen_file *kseigf, int Nvecs, double *eigVal, su3_vector **eigVec){

  FILE *fp;
  char *filename;
  int byterevflag;

  int rcv_rank, rcv_coords;
  int parity, ndata = 0, data_count = 0;
  int destnode;
  int k, x, y, z, t, ivecs;
  int status;
  int buf_length = 0, where_in_buf = 0;
  ks_eigen_check test_kseigc;
  u_int32type *val;
  int rank29, rank31;
  su3_vector *eigbuf = NULL;
  int idest = 0;
  double tmp;

  struct {
    su3_vector ksv;
    char pad[PAD_SEND_BUF]; /* Introduced because some switches
			       perform better if message lengths are longer */
  } msg;

  char myname[] = "r_serial_ks_eigen";

  fp = kseigf->fp;
  filename = kseigf->filename;
  byterevflag = kseigf->byterevflag;
  status = 0;

  if(this_node == 0){
    if(kseigf->parallel == PARALLEL)
      printf("%s: Attempting serial read from parallel file \n", myname);

    /* Read eigenvalues */
    status += sread_byteorder(byterevflag, fp, &(kseigf->parity), sizeof(int),
			      myname, "parity");
    status += sread_byteorder(byterevflag, fp, &(kseigf->Nvecs), sizeof(int),
			      myname, "Nvecs");
    for(ivecs = 0; ivecs < kseigf->Nvecs; ivecs++){
      status += sread_byteorder(byterevflag, fp, &tmp, sizeof(double),
				myname, "eigVal");
      if(ivecs < Nvecs) eigVal[ivecs] = tmp;
    }

    eigbuf = (su3_vector *)malloc(MAX_BUF_LENGTH*sizeof(su3_vector));
    if(eigbuf == NULL){
      printf("%s: Node %d can't malloc eigbuf\n", myname, this_node);
      fflush(stdout);
      terminate(1);
    }

    ndata = (kseigf->parity == EVEN || kseigf->parity == ODD) ? volume/2 : volume;
    ndata *= kseigf->Nvecs;
  }
  broadcast_bytes((char *)&status, sizeof(int));
  if(status != 0) return status;
  broadcast_bytes((char *)&(kseigf->parity), sizeof(int));
  broadcast_bytes((char *)&(kseigf->Nvecs), sizeof(int));
  broadcast_bytes((char *)eigVal, Nvecs*sizeof(double));

  /* all nodes initialize checksums */
  test_kseigc.sum31 = 0;
  test_kseigc.sum29 = 0;
  /* counts 32-bit words mod 29 and mod 31 in order of appearance on file */
  /* Here all nodes see the same sequence because we read serially */
  rank29 = 0;
  rank31 = 0;
  
  g_sync();

  /* Node 0 reads and deals out the values */
  for(ivecs = 0; ivecs < kseigf->Nvecs; ivecs++)
  for(rcv_rank = 0; rcv_rank < volume; rcv_rank++){
    /* If file is always in coordinate natural order */  
    rcv_coords = rcv_rank;

    x = rcv_coords % nx; rcv_coords /= nx;
    y = rcv_coords % ny; rcv_coords /= ny;
    z = rcv_coords % nz; rcv_coords /= nz;
    t = rcv_coords % nt;

    parity = 2 - (x + y + z + t)%2;
    if(parity == kseigf->parity || kseigf->parity == EVENANDODD){
      /* The node that gets the next su3_vector */
      destnode = node_number(x, y, z, t);

      /* Read eigenvectors */
      if(this_node==0){ /* Node 0 fills its buffer, if necessary */
	if(where_in_buf == buf_length){ /* get new buffer*/
	  /* new buffer length  = remaining data length, but never bigger 
	     than MAX_BUF_LENGTH */
	  buf_length = ndata - data_count;
	  if(buf_length > MAX_BUF_LENGTH) buf_length = MAX_BUF_LENGTH; 
	  /* then do read */
	  //	  if((int)fread(eigbuf, sizeof(su3_vector), buf_length, fp) != buf_length){
	  idest = (int)fread(eigbuf, sizeof(su3_vector), buf_length, fp);
	  if(idest != buf_length){
	    if(status == 0){
	      printf("%s: node %d eigenvector read error %d file %s\n",
		     myname, this_node, errno, filename); 
	      fflush(stdout);
	    } 
	    status++;
	  }
	  where_in_buf = 0;  /* reset counter */
	}  /*** end of the buffer read ****/

	/* Save vector in msg structure for further processing */
	msg.ksv = eigbuf[where_in_buf];
	if(destnode == 0)
	  /* just copy su3_vector */
	  idest = node_index(x, y, z, t); 
	else
	  /* send to correct node */
	  send_field((char *)&msg, sizeof(msg), destnode); 
	where_in_buf++;
	data_count++;
      }
      /* The node that contains this site reads the message */
      else {	/* for all nodes other than node 0 */
	if(this_node == destnode){
	  idest = node_index(x, y, z, t);
	  /* Receive padded message in msg */
	  get_field((char *)&msg, sizeof(msg),0);
	}
      } /* if(parity == kseigf->parity || kseigf->parity == EVENANDODD) */

      /* The receiving node does the byte reversal and then checksum,
	 if needed.  At this point msg.ksv contains the input vector
	 and idest points to the destination site structure. */
      if(this_node == destnode){
	if(byterevflag == 1)
	  byterevn((int32type *)(&msg.ksv), sizeof(su3_vector)/sizeof(int32type));

	/* Accumulate checksums */
	for(k = 0, val = (u_int32type *)(&msg.ksv); 
	    k < (int)sizeof(su3_vector)/(int)sizeof(int32type); k++, val++){
	  test_kseigc.sum29 ^= (*val)<<rank29 | (*val)>>(32-rank29);
	  test_kseigc.sum31 ^= (*val)<<rank31 | (*val)>>(32-rank31);
	  rank29++; if(rank29 >= 29)rank29 = 0;
	  rank31++; if(rank31 >= 31)rank31 = 0;
	}
      
	if(ivecs < Nvecs) eigVec[ivecs][idest] = msg.ksv;
      }
      else{
	rank29 += sizeof(su3_vector)/sizeof(int32type);
	rank31 += sizeof(su3_vector)/sizeof(int32type);
	rank29 %= 29;
	rank31 %= 31;
      }
    } /* if(parity == kseigf->parity || kseigf->parity == EVENANDODD) */
  }
  broadcast_bytes((char *)&status, sizeof(int));
  if(status != 0) return status;

  /* Combine node checksum contributions with global exclusive or */
  g_xor32(&test_kseigc.sum29);
  g_xor32(&test_kseigc.sum31);
  
  if(this_node == 0){
    printf("Read eigenvectors serially from file %s\n", filename);

    /* Read check record */
    status += sread_byteorder(byterevflag, fp, &kseigf->check.sum29,
			      sizeof(kseigf->check.sum29), myname, "check.sum29");
    status += sread_byteorder(byterevflag, fp, &kseigf->check.sum31,
			      sizeof(kseigf->check.sum31), myname, "check.sum31");
      
    /* Verify checksum */
    /* Checksums not implemented until version 5 */
    // if(kseigh->magic_number == KS_EIGEN_VERSION_NUMBER){
    if(kseigf->check.sum29 != test_kseigc.sum29 ||
       kseigf->check.sum31 != test_kseigc.sum31){
      printf("%s: Checksum violation file %s\n", myname, kseigf->filename);
      printf("Computed checksum %x %x.  Read %x %x.\n",
	     test_kseigc.sum29, test_kseigc.sum31,
	     kseigf->check.sum29, kseigf->check.sum31);
    }
    /* else
       printf("Checksums %x %x OK for file %s\n",
       kseigf->check.sum29, kseigf->check.sum31, kseigf->filename); */
    fflush(stdout);
    // }

    free(eigbuf);
  }
  broadcast_bytes((char *)&status, sizeof(int));

  return status;
} /* r_serial_ks_eigen */

/*----------------------------------------------------------------------*/

/* This subroutine closes the file and frees associated structures */
void r_serial_ks_eigen_f(ks_eigen_file *kseigf){

  char myname[] = "r_serial_ks_eigen_f";

  g_sync();

  if(kseigf == NULL) return;

  if(this_node==0){
    if(kseigf->parallel == PARALLEL)
      printf("%s: Attempting serial close on parallel file\n", myname);
      
    if(kseigf->fp != NULL) fclose(kseigf->fp);
    printf("Closed KS eigenvector file %s\n", kseigf->filename);
    fflush(stdout);
  }

  destroy_ks_eigen_file_handle(kseigf);
} /* r_serial_ks_eigen_f */

/*---------------------------------------------------------------------------*/

/* ASCII file format
   format:
   //version_number (int)
   time_stamp (char string enclosed in quotes)
   nx ny nz nt (int)
   Nvecs (int)
   for(ivecs=...) eigVal
   for(ivecs=...)for(t=...)for(z=...)for(y=...)for(x=...){
     for(i=...){eigVec[i].real, eigVec[i].imag}
   }
*/

/*---------------------------------------------------------------------------*/

/* Open and write header info for ASCII KS eigenvector file */
ks_eigen_file *w_ascii_ks_eigen_i(char *filename, int parity){

  FILE *fp;
  ks_eigen_file *kseigf;
  ks_eigen_header *kseigh;
  char myname[] = "w_ascii_ks_eigen_i";

  kseigf = create_output_ks_eigen_file_handle();
  kseigh = kseigf->header;

  /* node 0 does all the writing */
  if(this_node==0){

    /* Set up KS eigenvector file and header structures & load header values */

    fp = fopen(filename,"w");
    if(fp == NULL){
      printf("%s; Can't open file %s, error %d\n", myname, filename, errno);
      terminate(1);
    }

    printf("Opened KS eigenvector file %s for serial writing\n", filename);

    kseigf->fp = fp;

    /* if((fprintf(fp,"%d\n", KS_EIGEN_VERSION_NUMBER)) != 1){
         printf("%s: Error in writing version number\n", myname);
         terminate(1);
       }*/
    if((fprintf(fp,"\"%s\"\n", kseigh->time_stamp)) == 0){
      printf("%s: Error in writing time stamp\n", myname);
      terminate(1);
    }    
    if((fprintf(fp,"%d\t%d\t%d\t%d\n", nx, ny, nz, nt)) == 0){
      printf("%s: Error in writing dimensions\n", myname);
      terminate(1);
    }
  }
  else kseigf->fp = NULL;

  /* Assign remaining values to KS eigenvector file structure */
  kseigf->parallel = 0;
  kseigf->filename = filename;
  kseigf->byterevflag = 0; /* Not used for writing */
  kseigf->parity = parity;

  return kseigf;
} /* w_ascii_ks_eigen_i */
  
/*---------------------------------------------------------------------------*/

/* Write ASCII KS eigenvector from field */
void w_ascii_ks_eigen(ks_eigen_file *kseigf, int  Nvecs, double *eigVal, su3_vector **eigVec,
		      double *resid){

  FILE *fp;
  int currentnode, newnode;
  int i, j, x, y, z, t, ivecs;
  su3_vector eigbuf;
  int node0 = 0;
  int parity;
  char myname[] = "w_ascii_ks_eigen";

  currentnode = 0;
  fp = kseigf->fp;

  /* first Nvecs */
  if(this_node == 0){
    kseigf->Nvecs = Nvecs;
    kseigf->eigVal = eigVal;
    kseigf->resid = resid;

    /* Node 0 writes info file */
    write_ks_eigen_info_file(kseigf);

    /* Write eigenvalues */
    if((fprintf(fp, "%d\n", kseigf->parity)) == EOF){
      printf("%s: Node %d parity write error %d file %s\n",
	     myname, this_node, errno, kseigf->filename); 
      fflush(stdout);
      terminate(1);
    }
    if((fprintf(fp, "%d\n", Nvecs)) == EOF){
      printf("%s: Node %d Nvecs write error %d file %s\n",
	     myname, this_node, errno, kseigf->filename); 
      fflush(stdout);
      terminate(1);
    }
    for(ivecs = 0; ivecs < Nvecs; ivecs++){
      if((fprintf(fp, "%.7e\n", eigVal[ivecs])) == EOF){
	printf("%s: Node %d eigenvalue write error %d file %s\n",
	       myname, this_node, errno, kseigf->filename); 
	fflush(stdout);
	terminate(1);
      }
    }

    fflush(fp);
  }

  /* Write eigenvectors */
  for(ivecs = 0; ivecs < Nvecs; ivecs++)
  for(t = 0; t < nt; t++)
  for(z = 0; z < nz; z++)
  for(y = 0; y < ny; y++)
  for(x = 0; x < nx; x++){
    parity = 2 - (x + y + z + t)%2;
    if(parity == kseigf->parity || kseigf->parity == EVENANDODD){
      newnode = node_number(x, y, z, t);
      if(newnode != currentnode){ /* switch to another node */
	/* Send a few bytes of garbage to tell newnode it's OK to send */
	if(this_node == 0 && newnode != 0) send_field((char *)&eigbuf, 1, newnode);
	if(this_node == newnode && newnode != 0) get_field((char *)&eigbuf, 1, 0);
	currentnode = newnode;
      }
      
      if(this_node==0){
	if(currentnode==0){
	  i = node_index(x, y, z, t);
	  eigbuf = eigVec[ivecs][i];
	}
	else
	  get_field((char *)&eigbuf, sizeof(su3_vector), currentnode);
	for(j = 0;  j < 3; j++){
	  if((fprintf(fp, "%.7e\t%.7e\n", (double)eigbuf.c[j].real,
		      (double)eigbuf.c[j].imag)) == EOF){
	    printf("%s: error writing eigenvector\n", myname); 
	    terminate(1);
	  }
	}
      }
      else{ /* for nodes other than 0 */
	if(this_node == currentnode){
	  i = node_index(x, y, z, t);
	  eigbuf = eigVec[ivecs][i];
	  send_field((char *)&eigbuf, sizeof(su3_vector), node0);
	}
      }
    } /* if(parity == kseigf->parity || kseigf->parity == EVENANDODD) */
  }

  g_sync();

  if(this_node == 0){
    fflush(fp);
    printf("Wrote eigenvectors to ASCII file  %s\n", kseigf->filename);
  }
} /* w_ascii_ks_eigen */

/*---------------------------------------------------------------------------*/

/* Close ASCII KS eigenvector file */
void w_ascii_ks_eigen_f(ks_eigen_file *kseigf){

  //char myname[] = "w_ascii_ks_eigen_f";

  g_sync();

  if(this_node==0){
    fflush(kseigf->fp);
    if(kseigf->fp != NULL) fclose(kseigf->fp);
    printf("Closed ASCII KS eigenvector file %s time stamp %s\n", kseigf->filename,
	   (kseigf->header)->time_stamp);
  }

  /* Free header and file structures */
  destroy_ks_eigen_file_handle(kseigf);
}

/*---------------------------------------------------------------------------*/

/* Open ASCII KS eigenvector file and read header information */
ks_eigen_file *r_ascii_ks_eigen_i(char *filename){

  FILE *fp;
  ks_eigen_file *kseigf;
  ks_eigen_header *kseigh;
  char myname[] = "r_ascii_ks_eigen_i";

  /* All nodes set up a KS eigenvector file and header structure for reading */
  kseigf = create_input_ks_eigen_file_handle(filename);
  kseigh = kseigf->header;

  /* File opened for serial reading */
  kseigf->parallel = 0;
  kseigf->byterevflag = 0;  /* Unused for ASCII */

  /* Indicate coordinate natural ordering */
  kseigh->order = NATURAL_ORDER;

  /* Node 0 alone opens a file and reads the header */
  if(this_node == 0){
    fp = fopen(filename, "r");
    if(fp == NULL){
      printf("%s: Node %d can't open file %s, error %d\n",
	     myname, this_node,filename,errno);
      fflush(stdout);
      terminate(1);
    }

    printf("Opened KS eigenvector file %s for serial reading\n", filename);

    kseigf->fp = fp;

    /*
    if((fscanf(fp, "%d", &kseigh->magic_number)) ! = 1){
      printf("%s: Error in reading version number\n", myname); 
      terminate(1);
    }
    if(kseigh->magic_number != KS_EIGEN_VERSION_NUMBER){
      printf("%s: Unrecognized magic number in KS eigenvector file header.\n", myname);
      printf("Expected %d but read %d\n",
	     KS_EIGEN_VERSION_NUMBER, kseigh->magic_number);
      terminate(1);
    }
    */
    if((fscanf(fp, "%*[\"]%[^\"]%*[\"]", kseigh->time_stamp)) != 1){
      printf("%s: Error reading time stamp\n", myname); 
      terminate(1);
    }
    if((fscanf(fp, "%d%d%d%d", &kseigh->dims[0], &kseigh->dims[1],
	       &kseigh->dims[2], &kseigh->dims[3])) != 4){
      printf("%s: Error reading lattice dimensions\n", myname); 
      terminate(1);
    }
    if(kseigh->dims[0] != nx || kseigh->dims[1] != ny ||
       kseigh->dims[2] != nz || kseigh->dims[3] != nt ){
      /* So we can use this routine to discover the dimensions,
	 we provide that if nx = ny = nz = nt = -1 initially
	 we don't die */
      if(nx != -1 || ny != -1 || nz != -1 || nt != -1){
	printf("%s: Incorrect lattice size %d,%d,%d,%d\n",
	       myname, kseigh->dims[0], kseigh->dims[1], kseigh->dims[2], kseigh->dims[3]);
	terminate(1);
      }
      else{
	nx = kseigh->dims[0];
	ny = kseigh->dims[1];
	nz = kseigh->dims[2];
	nt = kseigh->dims[3];
	volume = nx*ny*nz*nt;
      }
    }
    kseigh->header_bytes = 0; /* Unused for ASCII */
  }

  else kseigf->fp = NULL;  /* Other nodes don't know about this file */

  /* Broadcasts the header structure from node 0 to all nodes */
  broadcast_bytes((char *)kseigh, sizeof(ks_eigen_header));

  return kseigf;
} /* r_ascii_ks_eigen_i */

/*---------------------------------------------------------------------------*/

/* Read eigenvectors */
/* 0 is normal exit code
   >1 for seek, read error, or missing data error */
int r_ascii_ks_eigen(ks_eigen_file *kseigf, int Nvecs, double *eigVal, su3_vector **eigVec){

  FILE *fp;
  int destnode;
  int i, j, x, y, z, t, ivecs;
  su3_vector eigbuf;
  int status;
  int parity;
  double tmp;
  char myname[] = "r_ascii_ks_eigen";

  fp = kseigf->fp;
  status = 0;

  if(this_node == 0){
    /* Read eigenvalues */
    if((fscanf(fp,"%d", &(kseigf->parity))) != 1){
      printf("%s: Error reading parity\n", myname);
      status++;
    }
    if((fscanf(fp,"%d", &(kseigf->Nvecs))) != 1){
      printf("%s: Error reading Nvecs\n", myname);
      status++;
    }
    for(ivecs = 0; ivecs < kseigf->Nvecs; ivecs++){
      if((fscanf(fp, "%lf\n", &tmp)) != 1){
	printf("%s: Error reading eigenvalues\n", myname);
	status++;
      }
      if(ivecs < Nvecs) eigVal[ivecs] = tmp;
    }
  }
  broadcast_bytes((char *)&status, sizeof(int));
  if(status != 0) return status;
  broadcast_bytes((char *)&(kseigf->parity), sizeof(int));
  broadcast_bytes((char *)&(kseigf->Nvecs), sizeof(int));
  broadcast_bytes((char *)eigVal, Nvecs*sizeof(double));

  g_sync();

  /* Read eigenvectors */
  for(ivecs = 0; ivecs < Nvecs; ivecs++)
  for(t = 0; t < nt; t++)
  for(z = 0; z < nz; z++)
  for(y = 0; y < ny; y++)
  for(x = 0; x < nx; x++){
    parity = 2 - (x + y + z + t)%2;
    if(parity == kseigf->parity || kseigf->parity == EVENANDODD){
      destnode = node_number(x, y, z, t);
      /* Node 0 reads, and sends site to correct node */
      if(this_node == 0){
	for(j = 0; j < 3; j++){
	  if((fscanf(fp, "%lf%lf\n",&(eigbuf.c[j].real), &(eigbuf.c[j].imag))) != 2){
	    if(status == 0) printf("%s: Error reading su3_vector\n", myname); 
	    status++;
	  }
	}
	if(destnode == 0){ /* just copy su3_vector */
	  i = node_index(x, y, z, t);
	  eigVec[ivecs][i] = eigbuf;
	}
	else /* send to correct node */
	  send_field((char *)&eigbuf, sizeof(su3_vector), destnode);
      }
      /* The node which contains this site reads message */
      else{ /* for all nodes other than node 0 */
	if(this_node == destnode){
	  get_field((char *)&eigbuf, sizeof(su3_vector), 0);
	  i = node_index(x, y, z, t);
	  eigVec[ivecs][i] = eigbuf;
	}
      }
    } /* if(parity == kseigf->parity || kseigf->parity == EVENANDODD) */
  }
  broadcast_bytes((char *)&status, sizeof(int));

  if(this_node == 0){
    fflush(fp);
    printf("Read eigenvectors to ASCII file  %s\n", kseigf->filename);
  }

  return status;
} /* r_ascii_ks_eigen */

/*---------------------------------------------------------------------------*/

/* Close KS eigenvector file */
void r_ascii_ks_eigen_f(ks_eigen_file *kseigf){

  FILE *fp;
  //char myname[] = "r_ascii_ks_eigen_f";

  fp = kseigf->fp;

  g_sync();

  if(this_node==0){
    printf("Closed ASCII KS eigenvector file  %s\n", kseigf->filename);
    fclose(fp);
    kseigf->fp = NULL;
    fflush(stdout);
  }
} /* r_ascii_ks_eigen_f */
