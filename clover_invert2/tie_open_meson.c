/************** tie_open_meson.c ***********************/
/* MIMD version 7 */
/* Print contents of an open meson correlator file */
/* C. DeTar 8 June 2008 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
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

#define CSUM(a,b) { (a).real += (b).real; (a).imag += (b).imag; }

#define GAMMAFIVE -1

/* Macros to multiply complex numbers by +-1 and +-i */
#define TIMESPLUSONE(a,b) { (b).real =  (a).real; (b).imag = (a).imag; }
#define TIMESMINUSONE(a,b) { (b).real =  -(a).real; (b).imag = -(a).imag; }
#define TIMESPLUSI(a,b) { (b).real = -(a).imag; (b).imag =  (a).real; }
#define TIMESMINUSI(a,b) { (b).real =  (a).imag; (b).imag = -(a).real; }

/*********************  gammatypes.h   **********************************
*									*
*  MIMD version 7 							*
*                                                                       *
*	A gamma matrix encoding scheme					*
*                                                                       *
*/

/*  DeGrand-Rossi convention

 gamma(XUP) 
 	    0  0  0  i
            0  0  i  0
            0 -i  0  0
           -i  0  0  0

 gamma(YUP)
 	    0  0  0 -1
            0  0  1  0
            0  1  0  0
           -1  0  0  0

 gamma(ZUP)
 	    0  0  i  0
            0  0  0 -i
           -i  0  0  0
            0  i  0  0

 gamma(TUP)
 	    0  0  1  0
            0  0  0  1
            1  0  0  0
            0  1  0  0

 gamma(FIVE) 
 	    1  0  0  0
            0  1  0  0
            0  0 -1  0
            0  0  0 -1
*/

/* Note: the order of the first five is fixed for compatibility with
   previous versions of mult_by_gamma */

enum gammatype {  GX, GY, GZ, GT, G5, 
		  GYZ, GZX, GXY, GXT, GYT, GZT, 
		  G5X, G5Y, G5Z, G5T, G1,
		  MAXGAMMA };

/* This encoding of gamma matrices works for representations in which
   matrix is just a permutation matrix with a phase change
   corresponding to the fourth roots of 1.  See the convention above.

   The matrix gamma_{r,c} (row index r, column index c) is encoded
   so that gamma.row[r].column = c is the only nonzero column in row r
   and gamma.row[r].phase specifies the phase of that matrix element

   */

typedef struct {
  int column;         /* encodes column with nonzero entry */
  int phase;          /* encodes phase of matrix element: 0,1,2,3 ->
                         1,i,-1,-i */
} gamma_row;

typedef struct {
  gamma_row row[4];
} gamma_matrix_t;

/* This list must agree exactly with the enum gammatype in gammatypes.h */

static char *gammalabel[MAXGAMMA] = {"GX", "GY", "GZ", "GT", "G5", 
				     "GYZ", "GZX", "GXY", "GXT", "GYT", "GZT", 
				     "G5X", "G5Y", "G5Z", "G5T", "G1" };

/* First four gamma matrices are initialized according to the
   conventions in gammatypes.h and the ones for gamma_x gamma_y
   gamma_z gamma_t are used to generate the rest */

static gamma_matrix_t gamma_matrix[MAXGAMMA] = {
   {{{3,1},  {2,1},  {1,3},  {0,3}}},
   {{{3,2},  {2,0},  {1,0},  {0,2}}},
   {{{2,1},  {3,3},  {0,3},  {1,1}}},
   {{{2,0},  {3,0},  {0,0},  {1,0}}}, 
};

static int gamma_initialized = 0;


/************* make_gammas.c **************************/

/* Starting from gamma_x, gamma_y, gamma_z, gamma_t, form the rest of the
   16 gamma matrices 
   Gamma matrices are represented via a "gamma_matrix" structure */

/* Does g3 = phase*g1 * g2  - used only to initialize the gamma matrices */
static void mult_gamma(int phase, gamma_matrix_t *g1, gamma_matrix_t *g2, gamma_matrix_t *g3)
{
  int r,s;
  for(r=0;r<4;r++)
    {
      s = g1->row[r].column;
      g3->row[r].column = g2->row[s].column;
      g3->row[r].phase  = (g1->row[r].phase + g2->row[s].phase + phase) % 4;
    }
}

/* For external use */
/* computes g3 = g1 * g2 */
void mult_gamma_by_gamma(gamma_matrix_t *g1, gamma_matrix_t *g2, 
			 gamma_matrix_t *g3)
{
  mult_gamma(0, g1, g2, g3);
}

static void make_gammas(void)
{
  int r;
  gamma_matrix_t *gamma = gamma_matrix;  

  /* gamma_yz = i * gamma_y * gamma_z = sigma_{yz}, etc */
  mult_gamma(1,&gamma[GY ],&gamma[GZ ],&gamma[GYZ]);
  mult_gamma(1,&gamma[GZ ],&gamma[GX ],&gamma[GZX]);
  mult_gamma(1,&gamma[GX ],&gamma[GY ],&gamma[GXY]);

  mult_gamma(1,&gamma[GX ],&gamma[GT ],&gamma[GXT]);
  mult_gamma(1,&gamma[GY ],&gamma[GT ],&gamma[GYT]);
  mult_gamma(1,&gamma[GZ ],&gamma[GT ],&gamma[GZT]);

  /* gamma_5t = gamma_5 * gamma_t = gamma_x * gamma_y * gamma_z */
  /* phase 3 -> -i compensates for the i in gamma_xy */
  mult_gamma(3,&gamma[GX ],&gamma[GYZ],&gamma[G5T]);
  /* gamma_5 = gamma_x * gamma_y * gamma_z * gamma_t */
  mult_gamma(0,&gamma[G5T],&gamma[GT ],&gamma[G5 ]);
  mult_gamma(0,&gamma[G5 ],&gamma[GX ],&gamma[G5X]);
  mult_gamma(0,&gamma[G5 ],&gamma[GY ],&gamma[G5Y]);
  mult_gamma(0,&gamma[G5 ],&gamma[GZ ],&gamma[G5Z]);
 
  /* G1 is the unit matrix */
  for(r=0;r<4;r++)
    {
      gamma[G1].row[r].column = r;
      gamma[G1].row[r].phase = 0;
    }
  gamma_initialized = 1;
}

/* Accessor */

gamma_matrix_t gamma_mat(enum gammatype i){
  if(gamma_initialized == 0)make_gammas();

  return gamma_matrix[i];
}


/************* mw_gamma.c (in su3.a) **************************/
/*
  Multiply a "Wilson vector" by a gamma matrix
  acting on the row index - equivalently, multiplying on the left
  usage:  mult_w_by_gamma( src, dest, dir)
	spin_wilson_vector *src,*dest;
	int dir;    dir = any of the gamma matrix types in gammatypes.h
*/


void mult_w_by_gamma_mat(wilson_vector * src, 
			 wilson_vector * dest, 
			 gamma_matrix_t *gm)
{
  register int c2,s2,s;	/* column indices, color and spin */

  if(gamma_initialized==0)make_gammas();

  for(s=0;s<4;s++){
    s2 = gm->row[s].column;
    switch (gm->row[s].phase){
    case 0:
      for(c2=0;c2<3;c2++){
	dest->d[s2].c[c2] = src->d[s].c[c2];}
      break;
    case 1:
      for(c2=0;c2<3;c2++){
	TIMESPLUSI( src->d[s].c[c2], dest->d[s2].c[c2] );}
      break;
    case 2:
      for(c2=0;c2<3;c2++){
	TIMESMINUSONE( src->d[s].c[c2], dest->d[s2].c[c2] );}
      break;
    case 3:
      for(c2=0;c2<3;c2++){
	TIMESMINUSI( src->d[s].c[c2], dest->d[s2].c[c2] );}
    }
  }
}

void mult_w_by_gamma(wilson_vector * src,
		     wilson_vector * dest, int dir)
{
  gamma_matrix_t gm;

  if(gamma_initialized==0)make_gammas();

  /* For compatibility */
  if(dir == GAMMAFIVE)dir = G5;

  if(dir >= MAXGAMMA)
    {
      printf("mult_w_by_gamma: Illegal gamma index %d\n",dir);
      exit(1);
    }

  gm = gamma_mat(dir);

  mult_w_by_gamma_mat(src, dest, &gm);
}

/************* msw_gamma_l.c (in su3.a) **************************/
/*
  Multiply a "spin Wilson vector" by a gamma matrix
  acting on the row index
  (This is the first index, or equivalently, multiplication on the left)
  usage:  mult_sw_by_gamma_l( src, dest, dir)
	spin_wilson_vector *src,*dest;
	int dir;    dir = any of the gamma matrix types in gammatypes.h


  Originally from MILC su3.a

  Modifications:
  UMH - modified for spin Wilson vector
  4/29/97 C. DeTar generalized to any gamma matrix.
*/


void mult_sw_by_gamma_mat_l(spin_wilson_vector * src, 
			    spin_wilson_vector * dest, 
			    gamma_matrix_t *gm)
{
  register int c2,s1,s2,s;	/* column indices, color and spin */

  if(gamma_initialized==0)make_gammas();

  for(s1=0;s1<4;s1++){
    s = gm->row[s1].column;
    switch (gm->row[s1].phase){
    case 0:
      for(s2=0;s2<4;s2++)for(c2=0;c2<3;c2++){
	dest->d[s1].d[s2].c[c2] = src->d[s].d[s2].c[c2];}
      break;
    case 1:
      for(s2=0;s2<4;s2++)for(c2=0;c2<3;c2++){
	TIMESPLUSI( src->d[s].d[s2].c[c2], dest->d[s1].d[s2].c[c2] );}
      break;
    case 2:
      for(s2=0;s2<4;s2++)for(c2=0;c2<3;c2++){
	TIMESMINUSONE( src->d[s].d[s2].c[c2], dest->d[s1].d[s2].c[c2] );}
      break;
    case 3:
      for(s2=0;s2<4;s2++)for(c2=0;c2<3;c2++){
	TIMESMINUSI( src->d[s].d[s2].c[c2], dest->d[s1].d[s2].c[c2] );}
    }
  }
}

void mult_sw_by_gamma_l(spin_wilson_vector * src,
			spin_wilson_vector * dest, int dir)
{
  gamma_matrix_t gm;

  if(gamma_initialized==0)make_gammas();

  /* For compatibility */
  if(dir == GAMMAFIVE)dir = G5;

  if(dir >= MAXGAMMA)
    {
      printf("mult_sw_by_gamma_l: Illegal gamma index %d\n",dir);
      exit(1);
    }

  gm = gamma_mat(dir);

  mult_sw_by_gamma_mat_l(src, dest, &gm);
}

/************* msw_gamma_r.c (in su3.a) **************************/


/* 
  Multiply a "Wilson matrix" (spin_wilson_vector) by a gamma matrix
  acting on the column index
  (This is the second index, or equivalently, multiplication on the right)
  usage:  mb_gamma_r( src, dest, dir)
	spin_wilson_vector *src,*dest;
	int dir;    dir = XUP, YUP, ZUP, TUP or GAMMAFIVE
*/

void mult_sw_by_gamma_mat_r(spin_wilson_vector * src,
			    spin_wilson_vector * dest, 
			    gamma_matrix_t *gm)
{
  register int c2,s1,s2,s;	/* column indices, color and spin */

  if(gamma_initialized==0)make_gammas();

  for(s=0;s<4;s++){
    s2 = gm->row[s].column;
    switch (gm->row[s].phase){
    case 0:
      for(s1=0;s1<4;s1++)for(c2=0;c2<3;c2++){
	dest->d[s1].d[s2].c[c2] = src->d[s1].d[s].c[c2];}
      break;
    case 1:
      for(s1=0;s1<4;s1++)for(c2=0;c2<3;c2++){
	TIMESPLUSI( src->d[s1].d[s].c[c2], dest->d[s1].d[s2].c[c2] );}
      break;
    case 2:
      for(s1=0;s1<4;s1++)for(c2=0;c2<3;c2++){
	TIMESMINUSONE( src->d[s1].d[s].c[c2], dest->d[s1].d[s2].c[c2] );}
      break;
    case 3:
      for(s1=0;s1<4;s1++)for(c2=0;c2<3;c2++){
	TIMESMINUSI( src->d[s1].d[s].c[c2], dest->d[s1].d[s2].c[c2] );}
    }
  }
}

void mult_sw_by_gamma_r(spin_wilson_vector * src,
			spin_wilson_vector * dest, int dir)
{
  gamma_matrix_t gm;

  if(gamma_initialized==0)make_gammas();

  /* For compatibility */
  if(dir == GAMMAFIVE)dir = G5;
  if(dir >= MAXGAMMA)
    {
      printf("mult_sw_by_gamma_r: Illegal gamma index %d\n",dir);
      exit(1);
    }

  gm = gamma_mat(dir);

  mult_sw_by_gamma_mat_r(src, dest, &gm);
}

/* Adjoint of gamma matrix */

void gamma_adj(gamma_matrix_t *dest, gamma_matrix_t *src){
  int r, c, p;
  int conj[4] = { 0, 3, 2, 1 };
  
  for(r = 0; r < 4; r++){
    c = src->row[r].column;
    p = src->row[r].phase;
    /* Transpose */
    dest->row[c].column = r; 
    /* Complex conjugate */
    dest->row[c].phase = conj[p];
  }
}

/* Transpose of gamma matrix */

void gamma_transp(gamma_matrix_t *dest, gamma_matrix_t *src){
  int r, c, p;
  
  for(r = 0; r < 4; r++){
    c = src->row[r].column;
    p = src->row[r].phase;
    /* Transpose */
    dest->row[c].column = r; 
    dest->row[c].phase = p;
  }
}

/* Complex conjugate */

void gamma_conj(gamma_matrix_t *dest, gamma_matrix_t *src){
  int r, c, p;
  int conj[4] = { 0, 3, 2, 1 };
  
  for(r = 0; r < 4; r++){
    c = src->row[r].column;
    p = src->row[r].phase;
    /* Complex conjugate */
    dest->row[r].column = c;
    dest->row[r].phase = conj[p];
  }
}

/* Map a label to the gamma index */

int gamma_index(char *label){
  int i;
  for(i = 0; i < MAXGAMMA; i++){
    if(strcmp(label,gammalabel[i]) == 0)return i;
  }
  return -1;  /* Error condition */
}

/* Map an index to the label */

char *gamma_label(int index){
  return gammalabel[index];
}

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
  printf("bytes_per_site:     %lu\n", fermiQCD_header.bytes_per_site);
  printf("sites:              %lu\n", fermiQCD_header.sites);

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
static void tie_wprop(wilson_propagator *wp, int t, int gam_snk){
  int ci,cf,si,sf;
  spin_wilson_vector sw1, sw2;
  complex c;

  c.real = 0.; c.imag = 0.;

  for(ci=0;ci<3;ci++){
    for(si = 0; si < 4; si++)
      for(sf = 0; sf < 4; sf++)
	for(cf=0; cf<3; cf++)
	  sw1.d[si].d[sf].c[cf] = wp->c[ci].d[si].d[sf].c[cf];

    /* Mult by sink gamma */
    mult_sw_by_gamma_l( &sw1, &sw2, gam_snk);

    /* Trace */
    for(sf = 0; sf < 4; sf++)
      CSUM(c, sw2.d[sf].d[sf].c[ci]);
  }

  printf("%d %14.6e %14.6e\n", t, c.real, c.imag);
}

/*--------------------------------------------------------------------*/
int main(int argc, char *argv[]){

  FILE *corr_fp;
  int t;
  char *gam_snk_lab;
  int gam_snk;
  char *filename;
  char scanfilename[256];
  char scangam_snk_lab[8];
  wilson_propagator wprop;

  /* Process command line args */
  if(argc < 2){
    fprintf(stderr, "Usage %s <filename> <gamma_label>\n", argv[0]);
    fprintf(stderr, "or just enter filename gamma_label now or Ctrl-d to abort\n");
    if(scanf("%s %s",scanfilename,scangam_snk_lab) != 2){
      return 1;
    }
    filename = scanfilename;
    gam_snk_lab = scangam_snk_lab;
  } else {
    filename = argv[1];
    gam_snk_lab = argv[2];
  }

  printf("Will read %s\n", filename);
  printf("and tie the correlator with gamma_label %s\n", gam_snk_lab);

  /* Decode gamma label */
  gam_snk = gamma_index(gam_snk_lab);
  if(gam_snk < 0){
    printf("\n%s is not a valid gamma matrix label\n",gam_snk_lab);
    return 1;
  }

  corr_fp = open_open_meson_file(filename);

  for(t=0; t<nt; t++){
    read_open_meson_prop(corr_fp, &wprop);
    tie_wprop(&wprop, t, gam_snk);
  }
  
  close_open_meson_file(corr_fp);

  return 0;
}
