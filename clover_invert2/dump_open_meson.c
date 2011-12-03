/************* dump_open_meson.c ********************/
/* MIMD version 7 */
/* Print contents of an open meson correlator file */
/* C. DeTar 8 June 2008 */

#include <stdio.h>
#include <string.h>
#include <stdint.h>

#define PRECISION 2

typedef struct {   
  float real;	   
  float imag; 
} fcomplex;  

/* specific for double complex */
typedef struct {
   double real;
   double imag;
} dcomplex;

#if (PRECISION==1)
#define complex fcomplex
#else
#define complex dcomplex
#endif

typedef struct { fcomplex e[3][3]; } fsu3_matrix;
typedef struct { fcomplex c[3]; } fsu3_vector;
typedef struct { fsu3_vector d[4]; } fwilson_vector;
typedef struct { fsu3_vector h[2]; } fhalf_wilson_vector;
typedef struct { fwilson_vector c[3]; } fcolor_wilson_vector;
typedef struct { fwilson_vector d[4]; } fspin_wilson_vector;
typedef struct { fcolor_wilson_vector d[4]; } fwilson_matrix;
typedef struct { fspin_wilson_vector c[3]; } fwilson_propagator;

typedef struct { dcomplex e[3][3]; } dsu3_matrix;
typedef struct { dcomplex c[3]; } dsu3_vector;
typedef struct { dsu3_vector d[4]; } dwilson_vector;
typedef struct { dsu3_vector h[2]; } dhalf_wilson_vector;
typedef struct { dwilson_vector c[3]; } dcolor_wilson_vector;
typedef struct { dwilson_vector d[4]; } dspin_wilson_vector;
typedef struct { dcolor_wilson_vector d[4]; } dwilson_matrix;
typedef struct { dspin_wilson_vector c[3]; } dwilson_propagator;

#if (PRECISION==1)

#define wilson_vector       fwilson_vector
#define half_wilson_vector  fhalf_wilson_vector
#define color_wilson_vector fcolor_wilson_vector
#define spin_wilson_vector  fspin_wilson_vector
#define wilson_matrix       fwilson_matrix
#define wilson_propagator   fwilson_propagator

#else

#define wilson_vector       dwilson_vector
#define half_wilson_vector  dhalf_wilson_vector
#define color_wilson_vector dcolor_wilson_vector
#define spin_wilson_vector  dspin_wilson_vector
#define wilson_matrix       dwilson_matrix
#define wilson_propagator   dwilson_propagator

#endif

static int nt;

/*--------------------------------------------------------------------*/
typedef uint32_t u_int32type;
typedef uint64_t u_int64type;

int dobyterev;

#include <assert.h>

/* For doing byte reversal on n contiguous 32-bit words */

void byterevn(u_int32type w[], int n)
{
  u_int32type old,newv;
  int j;

  assert(sizeof(u_int32type) == 4);
  
  for(j=0; j<n; j++)
    {
      old = w[j];
      newv = old >> 24 & 0x000000ff;
      newv |= old >> 8 & 0x0000ff00;
      newv |= old << 8 & 0x00ff0000;
      newv |= old << 24 & 0xff000000;
      w[j] = newv;
    }
} /* byterevn */

/*--------------------------------------------------------------------*/
/* Do byte reversal on n contiguous 64-bit words */

void byterevn64(u_int32type w[], int n)
{
  u_int32type tmp;
  int j;

  assert(sizeof(u_int32type) == 4);
  
  /* First swap pairs of 32-bit words */
  for(j=0; j<n; j++){
    tmp = w[2*j];
    w[2*j] = w[2*j+1];
    w[2*j+1] = tmp;
  }

  /* Then swap bytes in 32-bit words */
  byterevn(w, 2*n);
}

/*--------------------------------------------------------------------*/
/* Open and read a FermiQCD header for the open meson correlator file */

typedef dcomplex mdp_complex;
static FILE* open_open_meson_file(char filename[]){

  /* Data members of the FermiQCD header class for the open meson correlator */
  struct {
    char  file_id[60];
    char  program_version[60];
    char  creation_date[60];
    u_int32type endianess;
    u_int32type ndim;
    u_int32type box[10];
    u_int32type bytes_per_site;
    u_int64type sites;
  } fermiQCD_header = 
      {
	"File Type: MDP FIELD",
	"milc_qcd version 7",
	"",
	0x87654321,
	1,
	{nt},
	144*sizeof(mdp_complex),
	nt
      };

  FILE *fp;
  int i;

  fp = fopen(filename,"rb");
  if(fp == NULL){
    printf("open_open_meson_file: Can't open %s for reading\n", filename);
    return NULL;
  }

  /* Read the header */
  if(fread(&fermiQCD_header, sizeof(fermiQCD_header), 1, fp) != 1){
    printf("open_open_meson_file: Can't read header from %s\n", filename);
    return NULL;
  }

  /* Do we need byte reversal? */
  if(fermiQCD_header.endianess != 0x87654321){
    u_int32type x = fermiQCD_header.endianess;
    byterevn(&x,1);
    if(x != 0x87654321){
      printf("Quitting: Can't figure out endianess: %x\n",fermiQCD_header.endianess);
      exit(1);
    } else {
      dobyterev = 1;
      printf("Reading with byte reversal\n");
    }
  } else {
    dobyterev = 0;
  }

  /* Fix the header data */
  if(dobyterev){
    byterevn(&fermiQCD_header.endianess,1);
    byterevn(&fermiQCD_header.ndim,1);
    byterevn(fermiQCD_header.box,fermiQCD_header.ndim);
    byterevn(&fermiQCD_header.bytes_per_site,1);
    byterevn64((u_int32type *)&fermiQCD_header.sites,1);
  }

  /* Set globals */
  nt = fermiQCD_header.box[0];

  /* Contents of header */
  fermiQCD_header.file_id[59] = '\0';
  printf("file_id:            %s\n", fermiQCD_header.file_id);
  fermiQCD_header.program_version[59] = '\0';
  printf("program_version:    %s\n", fermiQCD_header.program_version);
  fermiQCD_header.creation_date[59] = '\0';
  printf("creation_date:      %s\n", fermiQCD_header.creation_date);
  printf("endianess:          %x\n", fermiQCD_header.endianess);
  printf("ndim:               %d\n", fermiQCD_header.ndim);
  for(i = 0; i < fermiQCD_header.ndim; i++)
    printf("box[%d]:             %d\n", i, fermiQCD_header.box[i]);
  printf("bytes_per_site:     %d\n", fermiQCD_header.bytes_per_site);
  printf("sites:              %d\n", fermiQCD_header.sites);

  return fp;
}
/*--------------------------------------------------------------------*/
static void read_open_meson_prop(FILE *fp, wilson_propagator *wp)
{
  mdp_complex c[4][4][3][3];
  int c0,c1,s0,s1;
  
  /* Read the open propagator element c */
  if(fread(&c, sizeof(c), 1, fp) != 1){
    printf("read_open_meson_prop: Error writing open meson correlator\n");
    return;
  }

  if(dobyterev){
    byterevn64((u_int32type *)&c[0][0][0][0],4*4*3*3*2);
  }

  /* Remap Wilson propagator elements */
  for(c0=0;c0<3;c0++)
    for(s0=0;s0<4;s0++)
      for(s1=0;s1<4;s1++)
	for(c1=0;c1<3;c1++)
	  wp->c[c0].d[s0].d[s1].c[c1] = c[s0][s1][c0][c1];

}

/*--------------------------------------------------------------------*/
static void close_open_meson_file(FILE *fp){
  if(fp != NULL)fclose(fp);
}

/*--------------------------------------------------------------------*/
static dump_wprop(wilson_propagator *wp, int t){
  int ci,cf,si,sf;

  printf("==============================================================\n");
  printf("t = %d\n", t);
  printf("==============================================================\n");

  for(si = 0; si < 4; si++)
    for(sf = 0; sf < 4; sf++){
      printf("-----------------------------------------------------------\n");
      printf("spins ( %d , %d )\n",si,sf);
      printf("-----------------------------------------------------------\n");
      for(ci=0;ci<3;ci++){
	for(cf=0;cf<3;cf++)
	  printf("( %10.3e , %10.3e )  ",
		 wp->c[ci].d[si].d[sf].c[cf].real,
		 wp->c[ci].d[si].d[sf].c[cf].imag);
	printf("\n");
      }
      printf("\n");
    }
}

/*--------------------------------------------------------------------*/
int main(int argc, char *argv[]){

  FILE *corr_fp;
  int t, tselect = -1;
  char *filename;
  char scanfilename[512];
  wilson_propagator wprop;

  /* Process command line args */
  /* Or take values from stdin */
  if(argc < 2){
    if(scanf("%s %d",scanfilename,&tselect) != 2){
      fprintf(stderr, "Usage %s <filename> <tselect>\n", argv[0]);
      return 1;
    }
    filename = scanfilename;
  } else {
    filename = argv[1];
    if (argc > 2)
      tselect = atoi(argv[2]);
  }

  printf("Will read %s\n", filename);
  if(tselect >= 0)
    printf("and dump the correlator at time %d\n", tselect);

  corr_fp = open_open_meson_file(filename);

  for(t=0; t<nt; t++){
    read_open_meson_prop(corr_fp, &wprop);
    if(tselect < 0 || t == tselect)
      dump_wprop(&wprop, t);
  }
  
  close_open_meson_file(corr_fp);

  return 0;
}
