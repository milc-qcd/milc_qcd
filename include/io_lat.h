#ifndef _IO_LAT_H
#define _IO_LAT_H
/************************* io_lat.h *************************************/

/* Macros, structures and prototypes for gauge configuration input and
   output */

/* Original by CD */
/* 
   2/17/98 added new gauge_info_keywords for hvy_qpot and ks_imp_dyn
*/

/* Definitions of restore and save lattice flags used in io_helpers.c */
/* Some are also used for Wilson propagators */
#define CONTINUE                         10
#define FRESH                            11
#define RELOAD_ASCII                     12
#define RELOAD_SERIAL                    13
#define RELOAD_MULTIDUMP                 18
#define RELOAD_PARALLEL                  19
#define FORGET                           40
#define SAVE_ASCII                       41
#define SAVE_SERIAL                      42
#define SAVE_CHECKPOINT                  43
#define SAVE_SERIAL_FM                   44
#define SAVE_SERIAL_FM_SC                45
#define SAVE_SERIAL_ILDG                 46
#define SAVE_PARALLEL_ILDG               47
#define SAVE_MULTIFILE_ILDG              48
#define SAVE_PARTFILE_ILDG               49
#define SAVE_SERIAL_SCIDAC               50
#define SAVE_PARALLEL_SCIDAC             51
#define SAVE_MULTIFILE_SCIDAC            52
#define SAVE_PARTFILE_SCIDAC             53
#define SAVE_PARALLEL                    54
#define SAVE_MULTIDUMP                   55
#define SAVE_SERIAL_ARCHIVE              56

/* Format for NERSC archive files */
#define ARCHIVE_3x2   0
#define ARCHIVE_3x3   1

/* For KS propagators */
#define SAVE_SERIAL_TSLICE 932
#define SAVE_ASCII_TSLICE 933

/* Helps in defining globals */
#ifdef CONTROL
#define EXTERN 
#else
#define EXTERN extern
#endif

#ifdef HAVE_UNISTD_H
#include <unistd.h>   /* For "write" and "close" "off_t" */
#endif
#include <sys/types.h> /* For "off_t" */
#include <stdio.h>
#include "../include/int32type.h"
#include "../include/macros.h"  /* For MAXFILENAME */
#include "../include/file_types.h"
#include "../include/su3.h"

/**********************************************************************/
/* Binary lattice formats                                             */
/**********************************************************************/
/* In version 5 we have two binary lattice file formats:
      serial files     Data written in coordinate natural order
      checkpoint files Data written in node dump order

      Further descriptive information is kept in a separate ASCII header
      file.  See below.

      */


/*--------------------------------------------------------------------*/
/* version 5 binary file format */

#define MAX_TIME_STAMP 64

/* 1. Header comes first    */

typedef struct {
  int32type magic_number;               /* Identifies file format */
  char   time_stamp[MAX_TIME_STAMP]; /* Date and time stamp - used to
					check consistency between the
					ASCII header file and the
					lattice file */
  int32type dims[4];                    /* Full lattice dimensions */
  int32type header_bytes;               /* NOT WRITTEN TO THE FILE but
					 helpful for finding the data */
  int32type order;                      /* 0 means no coordinate list is
				        attached and the values are in
				        coordinate serial order.
				        Nonzero means that a
				        coordinate list is attached,
				        specifying the order of values */
} gauge_header;


/* 2. Site list (ONLY for checkpoint files - i.e. node-dump order files)

      A listing of site coordinates for the data in this file
      in the order of appearance.  The number of coordinates must
      be exactly nx*ny*nz*nt.  The site coordinate is encoded
      as nx*(ny*(nz*t + z) + y) + x in a 32-bit integer.   

      */

/* 3. Next, the gauge field link matrices appear */

/* 4. Finally comes a checksum  */

typedef struct {
  u_int32type sum31;
  u_int32type sum29;
} gauge_check;

/**********************************************************************/
/* Info file format */

/* List of admissible keywords for version 5 ASCII lattice info file
   
   The information file describes the action that was used to create
   the lattice file in whatever detail you want.

   This list is intended to represent a union of the requirements
   for all applications.  A given application uses a subset.

   We keep this list to enforce some degree of consistency.
   Add more as needed.

   */

#ifdef CONTROL
char *gauge_info_keyword[] = {
      "magic_number",
      "time_stamp",
      "checksums",
      "nx",
      "ny",
      "nz",
      "nt",
      "action.description",

      "gauge.description",
      "gauge.beta11",
      "gauge.beta12",
      "gauge.tadpole.u0",
      "gauge.ssplaq",
      "gauge.stplaq",
      "gauge.nersc_linktr",
      "gauge.nersc_checksum",
      "gauge.nloops",
      "gauge.nreps",
      "gauge.previous.filename",
      "gauge.previous.time_stamp",
      "gauge.previous.checksums",
      "gauge.fix.description",
      "gauge.fix.tolerance",
      "gauge.translate",
      "gauge.smear.description",
      "gauge.smear.steps",
      "gauge.smear.factor",   /* (1-f)*link + f/6*staples */

      "quark.description",
      "quark.dyn_flavors",
      "quark.dyn_mass",
      "quark.flavors",
      "quark.flavors1",
      "quark.flavors2",
      "quark.mass",
      "quark.mass1",
      "quark.mass2",
      "quark.kappa",
      "quark.link.c1",
      "quark.link.c3",
      "quark.staple.w3",      /* link + w3*staples */
      "quark.clover.c0",
      "quark.clover.u0",
      ""       /* Last entry MUST be a zero-length keyword */
};
#else
extern char *gauge_info_keyword[];
#endif

/* Used to create info file name */

#define ASCII_INFO_EXT ".info"

/**********************************************************************/
/* 1996 Binary file format follows.  Note, this format was not in wide
  use.  We keep it for compatibility. */

#define MAX_GAUGE_FIELD_DESCRIPT 200
#define MAX_GAUGE_FIELD_PARAM 2
#define MAX_DIRAC_DESCRIPT 200
#define MAX_DIRAC_PARAM 3
#define MAX_SOURCE_DESCRIPT 200
#define MAX_SOURCE_PARAM 2
#define MAX_SOURCE_SPINS 4

/* Begin definition of header stuctures */

/* Note that an effort is made to make all Realing point and integer
   fields 32 bits long.  However, byte ordering may vary across
   platforms, and no effort is made in writing the file to produce
   a standard byte order.  The input routines attempt to compensate
   for byte reversal automatically, by examining the magic number at
   the beginning of the file */

/* 1. Header comes first */

typedef struct {
  int32type magic_number;          /* Identifies file format */
  int32type dims[4];               /* Full lattice dimensions */
  int32type header_bytes;          /* Number of bytes for data belonging to
				   this structure -- NOT necessarily 
				   the length of this structure! */
  int32type order;                 /* 0 means no coordinate list is attached
				   and the values are in coordinate serial order
				   Nonzero means that a coordinate list is attached,
				   specifying the order of values */
  struct {                      /* Gauge field parameters */
    int32type n_descript;          /* Number of bytes in character string */
    char   descript[MAX_GAUGE_FIELD_DESCRIPT];  /* Describes gauge field */
    int32type n_param;             /* Number of gauge field parameters */
    Real  param[MAX_GAUGE_FIELD_PARAM];        /* GF parameters */
  } gauge_field;
} gauge_header_1996  ;


/* Gauge Connection Header */

typedef struct {  /* Structure to hold header tokens */
  int ntoken;
  char **token;
  char **value;
} QCDheader ;

/* FNAL generic header */
typedef struct {
  int32type magic_number;          /* Identifies file format */
  int32type gmtime_stamp;             /* Used in FNAL header from call to time*/
  int32type size_of_element;	  /* bytes per data value, i.e, 4=single 						   precision */
  int32type elements_per_site;	   /* number of data elements stored at 
					each site */
  int32type dims[4];               /* Full lattice dimensions */
  int32type order;                 /* 0 means no coordinate list is attached
			   and the values are in coordinate serial order.
			  Nonzero is not a Fermilab option */
} FNALheader ;

/* 2. Parallel files only:

      Next comes a listing of site coordinates for the data in this file
      in the order of appearance.  The number of coordinates must
      be exactly nx*ny*nz*nt.  The site coordinate is encoded
      as nx*(ny*(nz*t + z) + y) + x in a 32-bit integer.   

      Serial files only:

      The site order of gauge field elements is required to be in subscript
      order (x,y,z,t) with x varying most rapidly, followed by y, etc. 
      so this list is omitted.

      */

/* 3. Next, a 32-bit checksum structure */

/* 4. Finally, the gauge field link matrices appear */


/**********************************************************************/
/* Versions 1-4 file format */

/* version number (32 bit integer) */

/* lattice dimensions nx,ny,nz,nt (4 32 bit integers) */

/* lattice constants c1, c2 (2 32 bit Realing number ) */

/* Link matrices in coordinate serial order */

/**********************************************************************/
typedef struct {
  FILE *         fp;            /* File pointer */
  gauge_header*  header;        /* Pointer to header for file */
  char *         filename;       /* Pointer to file name string */
  int            byterevflag;   /* Byte reverse flag - used only for reading */
  int32type *       rank2rcv;      /* File site list - used only for 
				   serial reading */ 
  int            parallel;      /* 1 if file was opened in parallel
				   0 if serial */
  int            dataformat;    /* 3x3 or 3x2 matrix (archive files only) */
  int            precision;     /* precision 1 or 2 (archive files for now) */
  gauge_check    check;         /* Checksum */
} gauge_file;

/**********************************************************************/
/* Globals specific to these file formats: */
/* Gauge Connection (Archive) format */

EXTERN  char ensemble_id[MAXFILENAME];
EXTERN  int sequence_number;

/**********************************************************************/
/* Prototypes for io_lat4.c */

int big_endian(void);

void read_lat_dim_gf(char *filename, int *ndim, int dims[]);
gauge_file *restore_ascii(char *filename);
gauge_file *save_ascii(char *filename);
gauge_file *restore_serial(char *filename);
gauge_file *save_serial(char *filename);
gauge_file *restore_parallel(char *filename);
gauge_file *save_parallel(char *filename);
gauge_file *save_checkpoint(char *filename);
gauge_file *save_serial_archive(char *filename);
gauge_file *save_parallel_archive(char *filename);
int write_gauge_info_item( FILE *fpout, /* ascii file pointer */
		       char *keyword,   /* keyword */
		       char *fmt,       /* output format -
					      must use s, d, f, or e */
		       char *src,       /* address of starting data */
		       int count,       /* number of data items if > 1 */
		       int stride);     /* byte stride of data if
					   count > 1 */
int sprint_gauge_info_item( 
  char *string,    /* character string */
  size_t nstring,     /* string length */			    
  char *keyword,   /* keyword */
  char *fmt,       /* output format -
		      must use s, d, e, f, or g */
  char *src,       /* address of starting data
		      floating point data must be
		      of type (Real) */
  int count,       /* number of data items if > 1 */
  int stride);     /* byte stride of data if
		      count > 1 */
gauge_file *setup_output_gauge_file(void);
gauge_file *setup_input_gauge_file(char *filename);

/* For compatibility */
gauge_file *save_old_binary(char *filename, Real c1, Real c2);

/**********************************************************************/
/* In gauge_info.c (application dependent) */
void write_appl_gauge_info(FILE *fp, gauge_file *gf);

/**********************************************************************/
/* Prototypes for io_scidac routines */
gauge_file *save_serial_scidac(char *filename);
gauge_file *save_parallel_scidac(char *filename);
gauge_file *save_multifile_scidac(char *filename);
gauge_file *save_partfile_scidac(char *filename);
gauge_file *save_serial_ildg(char *filename, char *stringLFN);
gauge_file *save_parallel_ildg(char *filename, char *stringLFN);
gauge_file *save_partfile_ildg(char *filename, char *stringLFN);
gauge_file *save_multifile_ildg(char *filename, char *stringLFN);
gauge_file *restore_serial_scidac(char *filename);
gauge_file *restore_parallel_scidac(char *filename);

/**********************************************************************/
/* Prototypes for parallel I/O routine interfaces:
   io_ansi.c, io_t3d.c, io_paragon2.c, etc  */

FILE *g_open(const char *filename, const char *mode);
int g_seek(FILE *stream, off_t offset, int whence);
size_t g_write(const void *ptr, size_t size, size_t nmemb,FILE *stream);
size_t g_read(void *ptr, size_t size, size_t nmemb, FILE *stream);
int g_close(FILE *stream);

/**********************************************************************/
/* Prototypes for io_lat_util.c routines */

int qcdhdr_get_str(char *s, QCDheader *hdr, char **q);
int qcdhdr_get_int(char *s,QCDheader *hdr,int *q);
int qcdhdr_get_int32x(char *s,QCDheader *hdr,u_int32type *q);
int qcdhdr_get_float(char *s, QCDheader *hdr, Real *q);
void error_exit(char *s);
void complete_U(float *u);
void complete_Ud(double *u);
QCDheader * qcdhdr_get_hdr(FILE *in);
void f2d_4mat(fsu3_matrix *a, su3_matrix *b);
void d2f_4mat(su3_matrix *a, fsu3_matrix *b);
void swrite_data(FILE* fp, void *src, size_t size, char *myname, char *descrip);
void pwrite_data(FILE* fp, void *src, size_t size, char *myname, char *descrip);
void pswrite_data(int parallel, FILE* fp, void *src, size_t size, 
		  char *myname, char *descrip);
int sread_data(FILE* fp, void *src, size_t size, char *myname, char *descrip);
int pread_data(FILE* fp, void *src, size_t size, char *myname, char *descrip);
int pread_byteorder(int byterevflag, FILE* fp, void *src, size_t size, char *myname, char *descrip);
int sread_byteorder(int byterevflag, FILE* fp, void *src, size_t size, char *myname, char *descrip);
int psread_data(int parallel, FILE* fp, void *src, size_t size, 
		char *myname, char *descrip);
int psread_byteorder(int byterevflag, int parallel, FILE* fp, 
		      void *src, size_t size, 
		     char *myname, char *descrip);
void pwrite_gauge_hdr(FILE *fp, gauge_header *gh);
void swrite_gauge_hdr(FILE *fp, gauge_header *gh);
int write_gauge_info_item( FILE *fpout,    /* ascii file pointer */
		       char *keyword,   /* keyword */
		       char *fmt,       /* output format -
					      must use s, d, e, f, or g */
		       char *src,       /* address of starting data
					   floating point data must be
					   of type (Real) */
		       int count,       /* number of data items if > 1 */
		       int stride);      /* byte stride of data if
                                           count > 1 */
int sprint_gauge_info_item( 
  char *string,    /* character string */
  size_t nstring,     /* string length */			    
  char *keyword,   /* keyword */
  char *fmt,       /* output format -
		      must use s, d, e, f, or g */
  char *src,       /* address of starting data
		      floating point data must be
		      of type (Real) */
  int count,       /* number of data items if > 1 */
  int stride);      /* byte stride of data if
		      count > 1 */
void write_generic_gauge_info(FILE *fp, gauge_file *gf);
void write_gauge_info_file(gauge_file *gf);
gauge_file *setup_input_gauge_file(char *filename);
gauge_file *setup_output_gauge_file(void);
void read_checksum(int parallel, gauge_file *gf, gauge_check *test_gc);
void write_checksum(int parallel, gauge_file *gf);
void read_site_list(int parallel,gauge_file *gf);
int read_v3_gauge_hdr(gauge_file *gf, int parallel, int *byterevflag);
int read_1996_gauge_hdr(gauge_file *gf, int parallel, int *byterevflag);
int read_fnal_gauge_hdr(gauge_file *gf, int parallel, int *byterevflag);
int read_gauge_hdr(gauge_file *gf, int parallel);
void write_site_list(FILE *fp, gauge_header *gh);
gauge_file *parallel_open(int order, char *filename);
fsu3_matrix *w_parallel_setup(gauge_file *gf, off_t *checksum_offset);
gauge_file *w_parallel_i(char *filename);
gauge_file *w_checkpoint_i(char *filename);
gauge_file *r_serial_i(char *filename);
void w_serial_f(gauge_file *gf);
void r_serial_f(gauge_file *gf);
void w_parallel_f(gauge_file *gf);
void r_parallel_f(gauge_file *gf);
char *read_info_file(char base_filename[]);

void byterevn(int32type w[], int n);


#endif /* _IO_LAT_H */
