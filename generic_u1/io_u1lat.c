/***************** io_u1lat.c************************************/

/* I/O utilities for the U(1) field                             */

/* MIMD version 7 */

/* ************************************************************ */
/*								*/
/*  1. ask_ending_u1_lattice()					*/
/*  2. save_u1_lattice()					*/
/*  3. save_u1_ascii()						*/
/*  4. save_u1_serial()						*/
/*  5. save_u1_parallel()					*/
/*  6. setup_output_u1gauge_file()				*/
/*  7. write_u1gauge_info_file()				*/
/*  8. ask_starting_u1_lattice()				*/
/*  9. reload_u1_lattice()					*/
/* 10. restore_u1_ascii()					*/
/* 11. restore_u1_serial()					*/
/* 12. restore_u1_parallel()					*/
/* 13. setup_input_u1gauge_file()				*/
/* 14. f2d_4cmplx()						*/
/* 15. d2f_4cmplx()						*/
/* 16. f2d_4real_cmplx()					*/
/* 17. d2f_4cmplx_real()					*/
/* 18. cold_u1lat()						*/
/* 19. w_u1_serial_i()						*/
/* 20. w_u1_serial()						*/
/* 21. w_u1_serial_f()						*/
/* 22. setup_u1_output_gauge_file()				*/
/* 23. flush_u1_lbuf_to_file()					*/
/* 24. flush_u1_tbuf_to_lbuf()					*/
/* 25. send_u1_buf_to_node0()					*/
/* 26. r_u1_serial()						*/
/* 27. r_u1_serial_i()						*/
/*								*/
/* Author S. Basak						*/
/* Last updated on 07.24.07					*/
/* Updated 1/04/2012 by SG to include binary saved lattice	*/
/*  binary format saves only (single precsion) real part of field u1gf	*/
/* CD 5/24/12 added SciDAC file support				*/
/*								*/
/*								*/
/* ************************************************************ */

#include <assert.h>
#include "generic_u1_includes.h"
#include "../include/io_u1lat.h"
#include "../include/generic_u1.h"
#include <string.h>
#include <errno.h>
#include <time.h>
#ifdef HAVE_QIO
#include "../include/io_scidac.h"
#endif
#define PARALLEL 1   /* Must evaluate to true */
#define SERIAL 0     /* Must evaluate to false */

#define NODE_DUMP_ORDER 1
#define NATURAL_ORDER 0

#undef MAX_BUF_LENGTH
#define MAX_BUF_LENGTH 4096

/* Forward declarations */
static gauge_file *w_u1_serial_i(char *filename);
static void w_u1_serial(gauge_file *gf);
static void w_u1_serial_f(gauge_file *gf);
static void flush_u1_tbuf_to_lbuf(gauge_file *gf, int *rank29, int *rank31,
                               float *lbuf, int *buf_length,
                               float *tbuf, int tbuf_length);
static void send_u1_buf_to_node0(float *tbuf, int tbuf_length, int currentnode);
static void r_u1_serial(gauge_file *gf);
static void accum_cksums(gauge_file *gf, int *rank29, int *rank31,
                         u_int32type *buf, int n);
static void flush_u1_lbuf_to_file(gauge_file *gf, float *lbuf,
				  int *buf_length);

/* find out what do to with lattice at
   end, and lattice name if necessary. */

int ask_ending_u1_lattice(FILE *fp,
			  int prompt, int *flag, char *filename)
{
  
  int status;
  char savebuf[256];
  char myname[] = "ask_ending_u1_lattice";
  
  if(prompt!=0)
    {
      printf("enter: 'forget_u1', 'save_u1_ascii', \n");
      printf("'save_u1_serial' or 'save_u1_parallel: \n");
      printf("'save_u1_serial_scidac' or 'save_u1_parallel_scidac: \n");
    }
  status=fscanf(fp,"%s",savebuf);
  if(status!=1)
    {
      node0_printf("%s: wrong ending lattice commandr --\n", myname);
      return(1);
    }
  printf("%s ",savebuf);
  
  if(strcmp("save_u1_ascii",savebuf)==0)
    *flag=SAVE_ASCII;
  else if(strcmp("save_u1_serial",savebuf)==0)
    *flag=SAVE_SERIAL;
  else if(strcmp("save_u1_parallel",savebuf)==0)
    *flag=SAVE_PARALLEL;
  else if(strcmp("save_u1_serial_scidac",savebuf)==0)
    *flag=SAVE_SERIAL_SCIDAC;
  else if(strcmp("save_u1_parallel_scidac",savebuf)==0)
    *flag=SAVE_PARALLEL_SCIDAC;
  else if(strcmp("forget_u1",savebuf)==0)
    *flag=FORGET;
  else
    {
      node0_printf("%s: \"%s\" is not a ", myname, savebuf);
      node0_printf("save u(1) lattice command\n");
      return(1);
    }
  
  if(*flag!=FORGET)
    {
      if(prompt!=0) printf("enter filename: \n");
      status=fscanf(fp,"%s",filename);
      if(status !=1)
	{
	  printf("%s: error reading filename\n", myname);
      return(1);
	}
      printf("%s\n",filename);
    }
  
  return(0);
  
} /* end of ask_ending_u1_lattice() */

/* save_lattice: writes lattice in ascii to a file */
gauge_file *save_u1_lattice(int flag,char *filename)
{

  double dtime;
  gauge_file *gf=NULL;

  u1plaq(&g_splaq,&g_tplaq);

  dtime=-dclock();
  switch(flag){
	case FORGET:
		gf=NULL;
		break;
	case SAVE_ASCII:
		gf=save_u1_ascii(filename);
		break;
	case SAVE_SERIAL:
		gf=save_u1_serial(filename);
		break;
	case SAVE_PARALLEL:
		gf=save_u1_parallel(filename);
		break;
#ifdef HAVE_QIO
	case SAVE_SERIAL_SCIDAC:
		gf=save_u1_serial_scidac(filename);
		break;
	case SAVE_PARALLEL_SCIDAC:
		gf=save_u1_parallel_scidac(filename);
		break;
#endif
	default:
	node0_printf("save_u1_lattice(): unknown type of saving!\n");
	terminate(1);
  }
  dtime+=dclock();

  if(flag!=FORGET){
    node0_printf("Saved lattice to ascii file %s!\n",gf->filename);
    node0_printf("Time to save = %e\n",dtime);
  }
#if (PRECISION==1)
  node0_printf("\nCHECK U(1) PLAQ: %e %e\n\n",g_splaq,g_tplaq);
#else
  node0_printf("\nCHECK U(1) PLAQ: %.16e %.16e\n\n",g_splaq,g_tplaq);
#endif

  return(gf);

} /* end of save_u1_lattice() */

/* saves u(1) lattice in ascii format serially */
gauge_file *save_u1_ascii(char *filename)
{

  FILE *fp=NULL;
  int currentnode,newnode;
  int i,x,y,z,t,dir;
  fcomplex lbuf[4];
  gauge_file *gf;
  gauge_header *gh;

  /* Setup gauge file and header structures and load values */
  gf=setup_output_u1gauge_file();
  gh=gf->header;

  /* Open gauge file */
  if(this_node==0)
    {
    fp=fopen(filename,"w");
    if(fp==NULL)
      {
      printf("Can't open file %s for gauge field!\n",filename);
      terminate(1);
      }
    gf->fp=fp;
    gf->parallel=0;
    gf->filename=filename;
    gf->byterevflag=0;            /* not used for writing */

    if((fprintf(fp,"%d\n",U1GAUGE_VERSION_NUMBER))==0)
      {
      printf("Error in writing version number\n"); terminate(1);
      }
    if((fprintf(fp,"\"%s\"\n",gh->time_stamp))==0)
      {
      printf("Error in writing time stamp\n");terminate(1);
      }
    if((fprintf(fp,"%d\t%d\t%d\t%d\n",nx,ny,nz,nt))==0)
      {
      printf("Error in writing dimensions\n");terminate(1);
      }
    write_u1gauge_info_file(gf);
    }

  /* Begin writing gauge field */
  g_sync();
  currentnode=0;

  for(t=0;t<nt;t++)for(z=0;z<nz;z++)
  for(y=0;y<ny;y++)for(x=0;x<nx;x++){
      newnode=node_number(x,y,z,t);
      if(newnode!=currentnode)
        {
        if(this_node==0 && newnode!=0) send_field((char *)lbuf,4,newnode);
        if(this_node==newnode && newnode!=0) get_field((char *)lbuf,4,0);
        }
      if(this_node==0)
        {
        if(currentnode==0)
          {
	    complex c[4];
	    i=node_index(x,y,z,t);
	    FORALLUPDIR(dir){
	      c[dir].real = u1_A[4*i+dir];
	      c[dir].imag = 0.;
	    }
	    d2f_4cmplx(c,lbuf);
          }
        else
          {
          get_field((char *)lbuf,4*sizeof(fcomplex),currentnode);
          }
        FORALLUPDIR(dir){
            if((fprintf(fp,"%.7e\t%.7e\n",
               (float)lbuf[dir].real,(float)lbuf[dir].imag))==EOF)
              {
              printf("Write error in save_u1_lattice()!\n");terminate(1);
              }
           }
        }
      else      /* for nodes other than 0 */
        {
	  if(this_node==currentnode)
	    {
	      complex c[4];
	      i=node_index(x,y,z,t);
	      FORALLUPDIR(dir){
		c[dir].real = u1_A[4*i+dir];
		c[dir].imag = 0.;
	      }
	      d2f_4cmplx(c,lbuf);
	      send_field((char *)lbuf,4*sizeof(fcomplex),0);
	    }
        }
     } /* x,y,z,t-loops end */

  g_sync();
  if(this_node==0)
    {
    fflush(fp);
    printf("Saved U(1) configuration to ascii file %s ...\n",gf->filename);
    printf("Time stamp %s\n",gh->time_stamp);
    fclose(fp);
    fflush(stdout);
    }

  return(gf);

} /* end of save_u1_ascii() */

/* suppose to save u(1) lattice serially */
gauge_file *save_u1_serial(char *filename)
{
  gauge_file *gf;

  gf = w_u1_serial_i(filename);
  w_u1_serial(gf);
  w_u1_serial_f(gf);

  return gf;


} /* end of save_u1_serial() */

/* suppose to save u(1) lattice parallel-ly */
gauge_file *save_u1_parallel(char *filename)
{

  printf("Can't handle parallel saving presently!!\n");
  terminate(1);

  return(NULL);

} /* end of save_u1_parallel() */

#ifdef HAVE_QIO
/* suppose to save u(1) lattice serially */
gauge_file *save_u1_serial_scidac(char *filename)
{
  save_real_scidac_from_field(filename, "", "", QIO_SINGLEFILE, u1_A, 4);
  return NULL;
} /* end of save_u1_serial() */

/* suppose to save u(1) lattice in parallel */
gauge_file *save_u1_parallel_scidac(char *filename)
{
  printf("Can't handle parallel saving presently!!\n");
  terminate(1);

  return(NULL);
} /* end of save_u1_parallel() */
#endif

/* setup_output_u1gauge_file: gauge header in gauge file */
gauge_file *setup_output_u1gauge_file(void)
{

  int i;
  gauge_file *gf;
  gauge_header *gh;
  time_t time_stamp;
  char myname[]="setup_output_u1gauge_file";

  /* Allocate space for a new file structure */
  assert(sizeof(int32type)==4);
  gf=(gauge_file *)malloc(sizeof(gauge_file));
  if(gf==NULL)
    {
    printf("%s: Can't malloc gf\n",myname);terminate(1);
    }
  gh=(gauge_header *)malloc(sizeof(gauge_header));
  if(gh==NULL)
    {
    printf("%s: Can't malloc gh\n",myname);terminate(1);
    }

  /* Load header values */
  gf->header=gh;
  gf->check.sum29=0;
  gf->check.sum31=0;
  gh->dims[0]=nx;
  gh->dims[1]=ny;
  gh->dims[2]=nz;
  gh->dims[3]=nt;
  gh->magic_number=U1GAUGE_VERSION_NUMBER;

  /* Get date and time stamp */
  if(this_node==0)
    {
    time(&time_stamp);
    strcpy(gh->time_stamp,ctime(&time_stamp));
    for(i=strlen(gh->time_stamp)+1;i<(int)sizeof(gh->time_stamp);i++)
        gh->time_stamp[i]='\0';
    if(gh->time_stamp[strlen(gh->time_stamp)-1]=='\n')
      gh->time_stamp[strlen(gh->time_stamp)-1]='\0';
    }
  /* Broadcast to all nodes */
  broadcast_bytes(gh->time_stamp,sizeof(gh->time_stamp));

  return(gf);

} /* end of setup_output_u1gauge_file() */

/* Open, write, and close the ASCII info file */
void write_u1gauge_info_file(gauge_file *gf)
{

  /* empty */

} /* end of write_u1gauge_info_file() */

/* ************************************************************ */

/* find out what kind of starting lattice
   to use, and lattice name if necessary. */
int ask_starting_u1_lattice(FILE *fp,
    int prompt,int *flag,char *filename)
{
  
  int status;
  char savebuf[256];
  char myname[] = "ask_starting_u1_lattice";
  
  if(prompt!=0)
    {
      printf("enter: 'fresh_u1', 'continue_u1' or 'reload_u1_ascii'\n");
      printf(" or 'reload_u1_serial' or 'reload_u1_parallel\n");
    }
  status=fscanf(fp,"%s",savebuf);
  if(status==EOF)
    {
      node0_printf("%s: EOF on STDIN!\n", myname);
      return(1);
    }
  if(status!=1)
    {
      node0_printf("%s: starting lattice", myname);
      node0_printf("command \"%s\" is invalid!\n",savebuf);
      return(1);
    }
  printf("%s \n",savebuf);
  
  if(strcmp("fresh_u1",savebuf)==0)
    *flag=FRESH;
  else if(strcmp("continue_u1",savebuf)==0)
    *flag=CONTINUE;
  else if(strcmp("reload_u1_ascii",savebuf)==0)
    *flag=RELOAD_ASCII;
  else if(strcmp("reload_u1_serial",savebuf)==0)
    *flag=RELOAD_SERIAL;
  else if(strcmp("reload_u1_parallel",savebuf)==0)
    *flag=RELOAD_PARALLEL;
  else
    {
      node0_printf("%s: lattice command", myname);
      node0_printf("\"%s\" is invalid!\n",savebuf);
      return(1);
    }
  
  /*read name of file and load it */
  if(*flag!=FRESH && *flag!=CONTINUE)
    {
      if(prompt!=0) printf("enter file containing lattice:\n");
      status=fscanf(fp,"%s",filename);
      if(status!=1)
	{
	  printf("ask_starting_u1_lattice(): error reading file name!\n");
	  return(1);
	}
      printf("%s\n",filename);
    }
  
  return(0);
  
} /* end of ask_starting_u1_lattice() */

/* read ascii lattice serially */
gauge_file *reload_u1_lattice(int flag, char *filename)
{
  
  double dtime;
  gauge_file *gf=NULL;
  
  dtime=-dclock();
  switch(flag){
  case FRESH:
    cold_u1lat();
    gf=NULL;
    break;
  case CONTINUE:
    gf=NULL;
    break;
  case RELOAD_ASCII:
    gf=restore_u1_ascii(filename);
    break;
  case RELOAD_SERIAL:
    gf=restore_u1_serial(filename);
    break;
  case RELOAD_PARALLEL:
    gf=restore_u1_parallel(filename);
    break;
  default:
    printf("reload_u1_lattice(): bad startflag %d!\n",flag);
    terminate(1);
  }
  dtime+=dclock();
  
  if(flag!=FRESH && flag!=CONTINUE)
    node0_printf("Time to reload gauge configuration = %e\n",dtime);
  u1plaq(&g_splaq,&g_tplaq);

#if (PRECISION==1)
  node0_printf("\nCHECK U(1) PLAQ: %e %e\n\n",g_splaq,g_tplaq);
  fflush(stdout);
#else
  node0_printf("\nCHECK U(1) PLAQ: %.16e %.16e\n\n",g_splaq,g_tplaq);
  fflush(stdout);
#endif
  
  return(gf);
  
} /* end of reload_u1_lattice() */

/* read u(1) lattice in ascii format serially */
gauge_file *restore_u1_ascii(char *filename)
{

  FILE *fp = NULL;
  int destnode;
  int version_number;
  int i,x,y,z,t,dir;
  fcomplex lbuf[4];
  gauge_header *gh;
  gauge_file *gf;

  /* Setup gauge file and header structure for reading */
  gf=setup_input_u1gauge_file(filename);
  gh=gf->header;

  /* File opened for serial reading */
  gf->parallel=0;
  if(this_node==0)
    {
    fp=fopen(filename,"r");
    if(fp==NULL)
      {
      printf("restore_u1_ascii(): Can't open file %s for reading!\n",
		filename);
      terminate(1);
      }
    gf->fp=fp;

    /* Read gauge file header structure */
    if((fscanf(fp,"%d",&version_number))!=1)
      {
      printf("restore_u1_ascii: Error reading version number\n");
      terminate(1);
      }
    gh->magic_number=version_number;
    if(gh->magic_number!=U1GAUGE_VERSION_NUMBER)
      {
      printf("restore_u1_ascii: Incorrect version number!\n");
      printf("\t read %d but expected %d\n",
		gh->magic_number,U1GAUGE_VERSION_NUMBER);
      terminate(1);
      }
    if((i=fscanf(fp,"%*[\f\n\r\t\v]%*[\"]%[^\"]%*[\"]",gh->time_stamp))!=1)
      {
      printf("restore_u1_ascii(): Error reading time stamp\n");
      printf("count %d time_stamp %s\n",i,gh->time_stamp);
      terminate(1);
      }
    if((fscanf(fp,"%d%d%d%d",&x,&y,&z,&t))!=4)
      {
      printf("restore_u1_ascii(): Error reading dimensions\n");
      terminate(1);
      }
    gh->dims[0]=x; gh->dims[1]=y;
    gh->dims[2]=z; gh->dims[3]=t;
    if(gh->dims[0]!=nx||gh->dims[1]!=ny||
       gh->dims[2]!=nz||gh->dims[3]!=nt)
      {
      if(nx!=-1||ny!=-1||nz!=-1||nt!=-1)
        {
        printf("restore_u1_ascii(): Incorrect lattice size %d,%d,%d,%d\n",
                    gh->dims[0],gh->dims[1],gh->dims[2],gh->dims[3]);
        terminate(1);
        }
      else
        {
        nx=gh->dims[0]; ny=gh->dims[1];
        nz=gh->dims[2]; nt=gh->dims[3];
        volume=nx*ny*nz*nt;
        }
      }
    } /* if node=0 */
  else gf->fp=NULL;

  /* Broadcast to all nodes */
  broadcast_bytes((char *)gh,sizeof(gauge_header));
  g_sync();

  /* Read gauge fields */
  for(t=0;t<nt;t++)for(z=0;z<nz;z++)
  for(y=0;y<ny;y++)for(x=0;x<nx;x++){
      destnode=node_number(x,y,z,t);
      if(this_node==0)
	{
        FORALLUPDIR(dir){
            if(fscanf(fp,"%e %e\n",
                    &(lbuf[dir].real),&(lbuf[dir].imag))!=2)
              {
              printf("restore_u1_ascii(): gauge field read error!\n");
              terminate(1);
              }
           }
	if(destnode==0)
	  {
	    complex c[4];
	    i=node_index(x,y,z,t);
	    f2d_4cmplx(lbuf,c);
	    FORALLUPDIR(dir){
	      u1_A[4*i+dir] = c[dir].real;
	    }
	  }
	else
	  {
	  send_field((char *)lbuf,4*sizeof(fcomplex),destnode);
	  }
	} /* if node=0 */
      else
	{
	if(this_node==destnode)
	  {
	    complex c[4];
	    get_field((char *)lbuf,4*sizeof(fcomplex),0);
	    i=node_index(x,y,z,t);
	    f2d_4cmplx(lbuf,c);
	    FORALLUPDIR(dir){
	      u1_A[4*i+dir] = c[dir].real;
	    }
	  }
	}

     } /* x,y,z,t-loops end */
  g_sync();

  if(this_node==0)
    {
    printf("Restored lattice from ascii file %s\n",filename);
    printf("Time stamp %s\n",gh->time_stamp);
    fclose(fp);
    gf->fp=NULL;
    fflush(stdout);
    }

  return(gf);

} /* end of restore_u1_ascii() */

/* suppose to read u(1) lattice serially */
gauge_file *restore_u1_serial(char *filename)
{
  
  gauge_file *gf = NULL;

  gf = r_u1_serial_i(filename);
  if(gf->header->magic_number == LIME_MAGIC_NO)
    {
      r_serial_f(gf);
      /* Close this reader and reread to get the header */
      free(gf->header);
      free(gf);
#ifdef HAVE_QIO
      restore_real_scidac_to_field(filename, QIO_SERIAL, u1_A, 4);
#else
      node0_printf("Looks like a SciDAC file.  Recompile with QIO.\n");
      terminate(1);
#endif
    }
  else if(gf->header->magic_number == U1GAUGE_VERSION_NUMBER_BINARY)
    {
      r_u1_serial(gf);
      r_serial_f(gf);
    }
  else
    {
      printf("Unrecognized U(1) binary format.\n");
      terminate(1);
    }
  
  return gf;
  
} /* end of restore_u1_serial() */

/* suppose to read u(1) lattice n parallel */
gauge_file *restore_u1_parallel(char *filename)
{

  printf("Can't handle parallel restoring presently!!\n");
  terminate(1);

  return(NULL);

} /* end of restore_u1_parallel() */

/* input gauge file and header */
gauge_file *setup_input_u1gauge_file(char *filename)
{

  gauge_file *gf;
  gauge_header *gh;
  char myname[]="setup_input_u1gauge_file";

  gf=(gauge_file *)malloc(sizeof(gauge_file));
  if(gf==NULL)
    {
    printf("%s: Can't malloc gf!\n",myname); terminate(1);
    }
  gf->filename = filename;

  assert(sizeof(int32type)==4);
  gh=(gauge_header *)malloc(sizeof(gauge_header));
  if(gh==NULL)
    {
    printf("%s: Can't malloc gh!\n",myname); terminate(1);
    }
  gf->header=gh;
  gf->rank2rcv = NULL;
  gf->check.sum29 = 0;
  gf->check.sum31 = 0;

  return(gf);

} /* end of setup_input_u1gauge_file() */

/* ************************************************************ */

/* Convert four single precision complex to generic precision */
void f2d_4cmplx(fcomplex *a,complex *b)
{

  int dir;

  for(dir=0;dir<4;dir++){
      b[dir].real=a[dir].real;
      b[dir].imag=a[dir].imag;
    }

} /* end of f2d_4cmplx() */

/* Convert four generic precision complex to single precision */
void d2f_4cmplx(complex *a,fcomplex *b)
{

  int dir;

  for(dir=0;dir<4;dir++){
      b[dir].real=a[dir].real;
      b[dir].imag=a[dir].imag;
    }

} /* end of d2f_4cmplx() */

/* Convert four single precision real to generic precision complex */
void f2d_4real_cmplx(float *a,complex *b)
{

  int dir;

  for(dir=0;dir<4;dir++){
      b[dir].real=a[dir];
      b[dir].imag=0.0;
    }

} /* end of f2d_4real_cmplx() */

/* Convert four generic precision complex to single precision real */
void d2f_4cmplx_real(complex *a,float *b)
{

  int dir;

  for(dir=0;dir<4;dir++){
      b[dir]=a[dir].real;
    }

} /* end of d2f_4cmplx_real() */

/* ************************************************************ */

/* sets u1-link variables to unity */
void cold_u1lat(void)
{

  register int i,dir;

  FORALLFIELDSITES(i)for(dir=XUP;dir<=TUP;dir++){
      u1_A[4*i+dir] = 0.0;
  }

} /* end of cold_u1lat() */

/* ************************************************************ */

/* Open a binary u1 file for serial writing by node 0 */

static gauge_file *w_u1_serial_i(char *filename)
{
  /* Only node 0 opens the file filename */
  /* Returns a file structure describing the opened file */

  char myname[] = "w_u1_serial_i";
  FILE *fp;
  gauge_file *gf;
  gauge_header *gh;

  /* Set up gauge file and gauge header structures and load header values */
  gf = setup_u1_output_gauge_file();
  gh = gf->header;

  /* Set number of nodes to zero to indicate coordinate natural ordering */

  gh->order = NATURAL_ORDER;

  /* Only node 0 opens the requested file */

  if(this_node == 0)
    {
      fp = fopen(filename, "wb");
      if(fp == NULL)
        {
          printf("%s: Node %d can't open file %s, error %d\n",
                 myname,this_node,filename,errno);fflush(stdout);
          terminate(1);
        }

      /* Node 0 writes the header */

      swrite_gauge_hdr(fp,gh);

    }

  /* Assign values to file structure */

  if(this_node==0)gf->fp = fp;
  else gf->fp = NULL;                /* Only node 0 knows about this file */

  gf->filename = filename;
  gf->byterevflag    = 0;            /* Not used for writing */
  gf->rank2rcv       = NULL;         /* Not used for writing */
  gf->parallel       = 0;

  return gf;
} /* w_u1_serial_i */

static void w_u1_serial(gauge_file *gf)
{
  /* gf  = file descriptor as opened by w_u1_serial_i */

  FILE *fp = NULL;
  gauge_header *gh = NULL;
  int rank29,rank31;
  float *lbuf = NULL;
  float *tbuf = NULL;
  int buf_length, tbuf_length;
  register int i,j;
  off_t offset;             /* File stream pointer */
  off_t coord_list_size;    /* Size of coordinate list in bytes */
  off_t head_size;          /* Size of header plus coordinate list */
  off_t checksum_offset = 0; /* Location of checksum */
  off_t gauge_check_size;   /* Size of checksum record */

  int currentnode,newnode;
  int x,y,z,t;
  int dir;
  char myname[] = "w_u1_serial";

  /* Allocate message buffer space for one x dimension of the local
     hypercube */
  /* The largest possible space we need is nx */

  tbuf = (float *)malloc(nx*4*sizeof(float));
  if(tbuf == NULL){
    printf("%s(%d): No room for tbuf\n",myname,this_node);
    terminate(1);
  }

  if(this_node==0)
    {
      if(gf->parallel)
        printf("w_serial: Attempting serial write to parallel file \n");

      lbuf = (float *)malloc(MAX_BUF_LENGTH*4*sizeof(float));
      if(lbuf == NULL)
        {
          printf("w_serial: Node 0 can't malloc lbuf\n");
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
    }

  /* Buffered algorithm for writing fields in serial (lexicographic) order */

  /* initialize checksums */
  gf->check.sum31 = 0;
  gf->check.sum29 = 0;
  /* counts 32-bit words mod 29 and mod 31 in order of appearance on file */
  /* Here only node 0 uses these values -- both start at 0 */
  rank29 = 4*sizeof(float)/sizeof(int32type)*sites_on_node*this_node % 29;
  rank31 = 4*sizeof(float)/sizeof(int32type)*sites_on_node*this_node % 31;

  g_sync();
  currentnode=0;  /* The node delivering data */

  buf_length = 0;
  tbuf_length = 0;
  for(j=0,t=0;t<nt;t++)for(z=0;z<nz;z++)for(y=0;y<ny;y++)for(x=0;x<nx;x++,j++)
    {
      newnode=node_number(x,y,z,t);  /* The node providing the next site */
      if(newnode != currentnode || x == 0){
        /* We are switching to a new node or have exhausted a line of nx */
        /* Sweep any data in the retiring node's tbuf to the node0 lbuf*/
        if(tbuf_length > 0){
          if(currentnode != 0)
            send_u1_buf_to_node0(tbuf, tbuf_length, currentnode);

          if(this_node == 0){
            /* Node 0 flushes tbuf and accumulates checksum */
            flush_u1_tbuf_to_lbuf(gf, &rank29, &rank31, lbuf, &buf_length,
                               tbuf, tbuf_length);
            /* Node 0 writes lbuf if full */
            if(buf_length > MAX_BUF_LENGTH - nx)
              flush_u1_lbuf_to_file(gf, lbuf, &buf_length);
          }
          tbuf_length = 0;

        }
        /* Node 0 sends a few bytes to newnode as a clear to send signal */
        if(newnode != currentnode){
          if( this_node==0 && newnode!=0 )send_field((char *)tbuf,4,newnode);
          if( this_node==newnode && newnode!=0 )get_field((char *)tbuf,4,0);
          currentnode=newnode;
        }
      } /* currentnode != newnode */ 

      /* The node with the data just appends to its tbuf */
      if(this_node == currentnode)
        {
          i=node_index(x,y,z,t);
	  FORALLUPDIR(dir){
	    tbuf[4*tbuf_length+dir] = u1_A[4*i+dir];
	  }
        }
  
      if(this_node == currentnode || this_node == 0)tbuf_length++;
  
    } /*close x,y,z,t loops */
  
  /* Purge any remaining data */
  
  if(tbuf_length > 0){
    if(currentnode != 0)
      send_u1_buf_to_node0(tbuf, tbuf_length, currentnode);
  }
  
  if(this_node == 0){
    flush_u1_tbuf_to_lbuf(gf, &rank29, &rank31, lbuf, &buf_length,
                       tbuf, tbuf_length);
    flush_u1_lbuf_to_file(gf, lbuf, &buf_length);
  } 

  g_sync();
  free(tbuf); 

  if(this_node==0)
    {
      free(lbuf);
      printf("Saved u1 gauge configuration serially to binary file %s\n",
             gf->filename);
      printf("Time stamp %s\n",gh->time_stamp);

      /* Write checksum */
      /* Position file pointer */ 
      if( g_seek(fp,checksum_offset,SEEK_SET) < 0 )
        {
          printf("w_u1_serial: Node %d g_seek %lld failed error %d file %s\n",
                 this_node,(long long)checksum_offset,errno,gf->filename);
          fflush(stdout);terminate(1);
        }
      write_checksum(SERIAL,gf);
    }

} /* w_u1_serial */

/* Close the file and free associated structures */
void w_u1_serial_f(gauge_file *gf)
{
  g_sync();
  if(this_node==0)
    {
      if(gf->parallel)
        printf("w_u1_serial_f: Attempting serial close on parallel file \n");

      g_close(gf->fp);
    }

  /* Node 0 writes ascii info file. Hope don't need special u1 routine */

  if(this_node == 0)write_gauge_info_file(gf);

  /* Do not free gf and gf->header so calling program can use them */

} /* w_u1_serial_f */

/* Set up the output gauge file an gauge header structure */

gauge_file *setup_u1_output_gauge_file()
{
  char myname[] = "setup_u1_output_gauge_file";
  gauge_file *gf;
  gauge_header *gh;
  time_t time_stamp;
  int i;

  /* Allocate space for a new file structure */

  gf = (gauge_file *)malloc(sizeof(gauge_file));
  if(gf == NULL)
    {
      printf("%s: Can't malloc gf\n",myname);
      terminate(1);
    }

  /* Allocate space for a new header structure */

  /* Make sure compilation gave us a 32 bit integer type */
  assert(sizeof(int32type) == 4);

  gh = (gauge_header *)malloc(sizeof(gauge_header));
  if(gh == NULL)
    {
      printf("%s: Can't malloc gh\n",myname);
      terminate(1);
    }

  /* Load header pointer and file name */
  gf->header = gh;

  /* Initialize */
  gf->check.sum29 = 0;
  gf->check.sum31 = 0;

  /* Load header values */

  gh->magic_number = U1GAUGE_VERSION_NUMBER_BINARY;

  gh->dims[0] = nx;
  gh->dims[1] = ny;
  gh->dims[2] = nz;
  gh->dims[3] = nt;

  /* Get date and time stamp. (We use local time on node 0) */

  if(this_node==0)
    {
      time(&time_stamp);
      strcpy(gh->time_stamp,ctime(&time_stamp));
      /* For aesthetic reasons, don't leave trailing junk bytes here to be
         written to the file */
      for(i = strlen(gh->time_stamp) + 1; i < (int)sizeof(gh->time_stamp); i++)
        gh->time_stamp[i] = '\0';

      /* Remove trailing end-of-line character */
      if(gh->time_stamp[strlen(gh->time_stamp) - 1] == '\n')
        gh->time_stamp[strlen(gh->time_stamp) - 1] = '\0';
    }

  /* Broadcast to all nodes */
  broadcast_bytes(gh->time_stamp,sizeof(gh->time_stamp));

  return gf;
} /* setup_u1_output_gauge_file */


/* Flush lbuf to output */
/* buf_length is reset */
static void flush_u1_lbuf_to_file(gauge_file *gf, float *lbuf,
                               int *buf_length)
{
  FILE *fp = gf->fp;

  if(*buf_length <= 0)return;
  if( (int)g_write(lbuf,4*sizeof(float),*buf_length,fp) !=
      *buf_length)
    {
      printf("w_u1_serial: Node %d gauge configuration write error %d file %s\n",
             this_node,errno,gf->filename);
      fflush(stdout);
      terminate(1);
    }
  *buf_length = 0;
}

/* Flush tbuf to lbuf and accumulate checksums */
/* tbuf_length is not reset here */
static void flush_u1_tbuf_to_lbuf(gauge_file *gf, int *rank29, int *rank31,
                               float *lbuf, int *buf_length,
                               float *tbuf, int tbuf_length)
{

  int nword;
  u_int32type *buf; 

  if(tbuf_length > 0){
    memcpy((void *)&lbuf[4*(*buf_length)],
           (void *)tbuf, 4*tbuf_length*sizeof(float));

    nword= 4*(int)sizeof(float)/(int)sizeof(int32type)*tbuf_length;
    buf = (u_int32type *)&lbuf[4*(*buf_length)];
    accum_cksums(gf, rank29, rank31, buf, nword);

    *buf_length += tbuf_length;
  }
}

static void send_u1_buf_to_node0(float *tbuf, int tbuf_length,
                              int currentnode)
{
  if(this_node == currentnode){
    send_field((char *)tbuf,4*tbuf_length*sizeof(float),0);
  }
  else if(this_node == 0){
    get_field((char *)tbuf,4*tbuf_length*sizeof(float),
              currentnode);
  }
}

/* Here only node 0 reads the gauge configuration from a binary file */

static void r_u1_serial(gauge_file *gf)
{
  /* gf  = gauge configuration file structure */

  FILE *fp;
  gauge_header *gh;
  char *filename;
  int byterevflag;

  off_t offset = 0 ;        /* File stream pointer */
  off_t gauge_check_size;   /* Size of gauge configuration checksum record */
  off_t coord_list_size;    /* Size of coordinate list in bytes */
  off_t head_size;          /* Size of header plus coordinate list */
  off_t checksum_offset = 0; /* Where we put the checksum */
  int rcv_rank, rcv_coords;
  int destnode;
  int k;
  int x,y,z,t;
  int dir;
  int buf_length = 0, where_in_buf = 0;
  gauge_check test_gc;
  u_int32type *val;
  int rank29,rank31;
  float *lbuf = NULL;
  float tmpu1[4];
  char myname[] = "r_u1_serial";
  int idest = 0;

  fp = gf->fp;
  gh = gf->header;
  filename = gf->filename;
  byterevflag = gf->byterevflag;

  if(this_node == 0)
    {
      /* Compute offset for reading gauge configuration */

	/* binary U1 files have a check sum */
        gauge_check_size = sizeof(gf->check.sum29) +
          sizeof(gf->check.sum31);

      if(gf->header->order == NATURAL_ORDER)coord_list_size = 0;
      else coord_list_size = sizeof(int32type)*volume;
      checksum_offset = gf->header->header_bytes + coord_list_size;
      head_size = checksum_offset + gauge_check_size;
      /* Allocate space for read buffer */

      if(gf->parallel)
        printf("%s: Attempting serial read from parallel file \n",myname);

      /* Allocate single precision read buffer */
      lbuf = (float *)malloc(MAX_BUF_LENGTH*4*sizeof(float));
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
                 myname,(long long)offset,errno,filename);
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

            if( (int)g_read(lbuf,4*sizeof(float),buf_length,fp)
                != buf_length)
              {
                printf("%s: node %d gauge configuration read error %d file %s\n",
                       myname,this_node,errno,filename);
                fflush(stdout); terminate(1);
              }
            where_in_buf = 0;  /* reset counter */
          }  /*** end of the buffer read ****/
  
        if(destnode==0){        /* just copy links */
          idest = node_index(x,y,z,t);
          /* Save 4 u1 phases in tmpu1 for further processing */
          memcpy(tmpu1,&lbuf[4*where_in_buf],4*sizeof(float));
        }
        else {          /* send to correct node */
          send_field((char *)&lbuf[4*where_in_buf],
                     4*sizeof(float),destnode);
        } 
        where_in_buf++;
      }
     
      /* The node that contains this site reads the message */
      else {    /* for all nodes other than node 0 */
        if(this_node==destnode){
          idest = node_index(x,y,z,t);

          /* Receive 4 u1 in temporary space for further processing */
          get_field((char *)tmpu1,4*sizeof(float),0);
        }
      }

      /* The receiving node does the byte reversal and then checksum,
         if needed.  At this point tmpu1 contains the input u1
         and idest points to the destination site structure. */

      if(this_node==destnode)
        {
          if(byterevflag==1)
            byterevn((int32type *)tmpu1,
                     4*sizeof(float)/sizeof(int32type));
          /* Accumulate checksums */
          for(k = 0, val = (u_int32type *)tmpu1;
              k < 4*(int)sizeof(float)/(int)sizeof(int32type);
              k++, val++)
            {
              test_gc.sum29 ^= (*val)<<rank29 | (*val)>>(32-rank29);
              test_gc.sum31 ^= (*val)<<rank31 | (*val)>>(32-rank31);
              rank29++; if(rank29 >= 29)rank29 = 0;
              rank31++; if(rank31 >= 31)rank31 = 0;
            }
          /* Copy 4 u1 to lattice[idest], converting to generic
             precision */
	  FORALLUPDIR(dir){
	    u1_A[4*idest+dir] = tmpu1[dir];
	  }
        }
      else
        {
          rank29 += 4*sizeof(float)/sizeof(int32type);
          rank31 += 4*sizeof(float)/sizeof(int32type);
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

      printf("Restored binary u1 gauge configuration serially from file %s\n",
             filename);
     if(gh->magic_number == U1GAUGE_VERSION_NUMBER_BINARY)
        {
          printf("Time stamp %s\n",gh->time_stamp);
          if( g_seek(fp,checksum_offset,SEEK_SET) < 0 )
            {
              printf("%s: Node 0 g_seek %lld failed error %d file %s\n",
                    myname,(long long)offset,errno,filename);
              fflush(stdout);terminate(1);   
            }
          read_checksum(SERIAL,gf,&test_gc);
        }
      fflush(stdout);
      free(lbuf);
    }

} /* r_u1_serial */


gauge_file *r_u1_serial_i(char *filename)
{
  /* Returns file descriptor for opened file */

  gauge_header *gh;
  gauge_file *gf;
  FILE *fp;
  int byterevflag = 0;

  /* All nodes set up a gauge file and gauge header structure for reading */

  gf = setup_input_gauge_file(filename);
  gh = gf->header;

  /* File opened for serial reading */
  gf->parallel = 0;

  /* Node 0 alone opens the file and reads the header */

  g_sync();

  if(this_node==0)
    {
      fp = g_open(filename, "rb");
      if(fp == NULL)
        {
          printf("r_u1_serial_i: Node %d can't open file %s, error %d\n",
                 this_node,filename,errno);fflush(stdout);terminate(1);
        }

      gf->fp = fp;

      byterevflag = read_u1_gauge_hdr(gf,SERIAL);
    }

  else gf->fp = NULL;  /* The other nodes don't know about this file */

  /* Broadcast the byterevflag from node 0 to all nodes */

  broadcast_bytes((char *)&byterevflag,sizeof(byterevflag));
  gf->byterevflag = byterevflag;

  /* Node 0 broadcasts the header structure to all nodes */

  broadcast_bytes((char *)gh,sizeof(gauge_header));

  /* Read site list and broadcast to all nodes */

  read_site_list(SERIAL,gf);

  return gf;

}/* r_u1_serial_i */

/* Accumulate checksums (copied from io_lat4.c) */
static void accum_cksums(gauge_file *gf, int *rank29, int *rank31,
                         u_int32type *buf, int n){
  int k;
  u_int32type *val;

  for(k = 0, val = buf; k < n; k++, val++)
    {
      gf->check.sum29 ^= (*val)<<(*rank29) | (*val)>>(32-(*rank29));
      gf->check.sum31 ^= (*val)<<(*rank31) | (*val)>>(32-(*rank31));
      (*rank29)++; if(*rank29 >= 29)*rank29 = 0;
      (*rank31)++; if(*rank31 >= 31)*rank31 = 0;
    }
}

int read_u1_gauge_hdr(gauge_file *gf, int parallel)
{
  /* parallel = 1 (TRUE) if all nodes are accessing the file */
  /*            0        for access from node 0 only */

  FILE *fp;
  gauge_header *gh;
  int32type tmp, btmp;
  int j;
  int byterevflag = 0;
  char myname[] = "read_gauge_hdr";

  fp = gf->fp;
  gh = gf->header;

  /* Read header, do byte reversal, 
     if necessary, and check consistency */
  
  /* Read and verify magic number */

  if(psread_data(parallel, fp,&gh->magic_number,sizeof(gh->magic_number),
			 myname,"magic number")!=0)terminate(1);

  tmp = gh->magic_number;
  btmp = gh->magic_number;
  byterevn((int32type *)&btmp,1);

  if(tmp == U1GAUGE_VERSION_NUMBER_BINARY) 
    {
      byterevflag=0;
      node0_printf("Reloading U1 gauge field\n");
    }
  else if(btmp == U1GAUGE_VERSION_NUMBER_BINARY) 
    {
      byterevflag=1;
      gh->magic_number = btmp;
      node0_printf("Reloading U1 gauge field\n");
      /**      printf("Reading with byte reversal\n"); **/
      if( sizeof(float) != sizeof(int32type)) {
	printf("%s: Can't byte reverse\n",myname);
	printf("requires size of int32type(%d) = size of float(%d)\n",
	       (int)sizeof(int32type),(int)sizeof(float));
	terminate(1);
      }
    }
  else
    {
	  /* End of the road. */
	  printf("%s: Unrecognized magic number in gauge configuration file header.\n",myname);
	  printf("Expected %x but read %x\n",
		 GAUGE_VERSION_NUMBER,tmp);
	  printf("Expected %lu but read %lu\n",
		 (unsigned long)GAUGE_VERSION_NUMBER,(unsigned long)tmp);
	  terminate(1);
      return byterevflag;
    }
  
  /* Read and process header information */
  /* Get lattice dimensions */

      if(psread_byteorder(byterevflag,parallel,fp,gh->dims,sizeof(gh->dims),
			  myname,"dimensions")!=0)terminate(1);

  /* Check lattice dimensions for consistency */

  if(gh->dims[0] != nx || 
     gh->dims[1] != ny ||
     gh->dims[2] != nz ||
     gh->dims[3] != nt)
    {
      /* So we can use this routine to discover the dimensions,
	 we provide that if nx = ny = nz = nt = -1 initially
	 we don't die */
      if(nx != -1 || ny != -1 || nz != -1 || nt != -1)
	{
	  printf("%s: Incorrect lattice dimensions ",myname);
	  for(j=0;j<4;j++)
	    printf("%d ",gh->dims[j]); 
	  printf("\n");fflush(stdout);terminate(1);
	}
      else
	{
	  nx = gh->dims[0];
	  ny = gh->dims[1];
	  nz = gh->dims[2];
	  nt = gh->dims[3];
	  volume = nx*ny*nz*nt;
	}
    }


    /* Read date and time stamp */
    
    if(psread_data(parallel,fp,gh->time_stamp,sizeof(gh->time_stamp),
		   myname,"time stamp")!=0)terminate(1);
    
    /* Read header byte length */
    
    gh->header_bytes = sizeof(gh->magic_number) + sizeof(gh->dims) + 
      sizeof(gh->time_stamp) + sizeof(gh->order);
    
    /* Read data order */
    
    if(psread_byteorder(byterevflag,parallel,fp,&gh->order,sizeof(gh->order),
			myname,"order parameter")!=0)terminate(1);

  return byterevflag;
  
} /* read_u1_gauge_hdr */


