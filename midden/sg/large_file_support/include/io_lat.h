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
#define CONTINUE 10
#define FRESH    11
#define RELOAD_ASCII  12
#define RELOAD_SERIAL  13
#define RELOAD_MULTIDUMP 18
#define RELOAD_PARALLEL  19
#define FORGET 20
#define SAVE_ASCII 21
#define SAVE_SERIAL 22
#define SAVE_CHECKPOINT 23
#define SAVE_MULTIDUMP 27
#define SAVE_PARALLEL 28
#define SAVE_OLD_BINARY 29

#ifndef NT
#include <unistd.h>   /* For "write" and "close" */
#endif

/* For portability we require a 32-bit integer type for writing
   integer data to binary files.  Machine-dependent. Add more as
   needed. Compile with -DSHORT32 for 64-bit integer machines.
   like the T3E */

#include "../include/type32.h"
#include <stdio.h>

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

#define GAUGE_VERSION_NUMBER 20103
#define MAX_TIME_STAMP 64

/* 1. Header comes first    */

typedef struct {
  type32 magic_number;               /* Identifies file format */
  char   time_stamp[MAX_TIME_STAMP]; /* Date and time stamp - used to
					check consistency between the
					ASCII header file and the
					lattice file */
  type32 dims[4];                    /* Full lattice dimensions */
  type32 header_bytes;               /* NOT WRITTEN TO THE FILE but
					 helpful for finding the data */
  type32 order;                      /* 0 means no coordinate list is
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
  u_type32 sum31;
  u_type32 sum29;
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
      "gauge.nloops",
      "gauge.nreps",
      "gauge.previous.filename",
      "gauge.previous.time_stamp",
      "gauge.previous.checksums",
      "gauge.fix.description",
      "gauge.fix.tolerance",
      "gauge.smear.description",
      "gauge.smear.steps",
      "gauge.smear.factor",   /* (1-f)*link + f/6*staples */

      "quark.description",
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

#define ASCII_GAUGE_INFO_EXT ".info"

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
#define GAUGE_VERSION_NUMBER_1996 53546

/* Begin definition of header stuctures */

/* Note that an effort is made to make all Realing point and integer
   fields 32 bits long.  However, byte ordering may vary across
   platforms, and no effort is made in writing the file to produce
   a standard byte order.  The input routines attempt to compensate
   for byte reversal automatically, by examining the magic number at
   the beginning of the file */

/* 1. Header comes first */

typedef struct {
  type32 magic_number;          /* Identifies file format */
  type32 dims[4];               /* Full lattice dimensions */
  type32 header_bytes;          /* Number of bytes for data belonging to
				   this structure -- NOT necessarily 
				   the length of this structure! */
  type32 order;                 /* 0 means no coordinate list is attached
				   and the values are in coordinate serial order
				   Nonzero means that a coordinate list is attached,
				   specifying the order of values */
  struct {                      /* Gauge field parameters */
    type32 n_descript;          /* Number of bytes in character string */
    char   descript[MAX_GAUGE_FIELD_DESCRIPT];  /* Describes gauge field */
    type32 n_param;             /* Number of gauge field parameters */
    Real  param[MAX_GAUGE_FIELD_PARAM];        /* GF parameters */
  } gauge_field;
} gauge_header_1996  ;

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

#define GAUGE_VERSION_NUMBER_V1 59354   /* Old number for compatibility */

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
  type32 *       rank2rcv;      /* File site list - used only for 
				   serial reading */ 
  int            parallel;      /* 1 if file was opened in parallel
				   0 if serial */
  gauge_check    check;         /* Checksum */
} gauge_file;

/**********************************************************************/
/* Prototypes for io_lat4.c */

gauge_file *restore_ascii(char *filename);
gauge_file *save_ascii(char *filename);
gauge_file *restore_serial(char *filename);
gauge_file *save_serial(char *filename);
gauge_file *restore_parallel(char *filename);
gauge_file *save_parallel(char *filename);
gauge_file *save_checkpoint(char *filename);
int write_gauge_info_item( FILE *fpout, /* ascii file pointer */
		       char *keyword,   /* keyword */
		       char *fmt,       /* output format -
					      must use s, d, f, or e */
		       char *src,       /* address of starting data */
		       int count,       /* number of data items if > 1 */
		       int stride);     /* byte stride of data if
					   count > 1 */

/* For compatibility */
gauge_file *save_old_binary(char *filename, Real c1, Real c2);

/**********************************************************************/
/* In gauge_info.c (application dependent) */
void write_appl_gauge_info(FILE *fp);

/**********************************************************************/
/* Prototypes for io_helpers.c */
gauge_file *save_lattice( int flag, char *filenamee, Real c1, Real c2 );
gauge_file *reload_lattice( int flag, char *filename);
int ask_starting_lattice( int prompt, int *flag, char *filename );
int ask_ending_lattice( int prompt, int *flag, char *filename );
void coldlat();
void funnylat();
int get_f( int prompt, char *variable_name_string, Real *value );
int get_i( int prompt, char *variable_name_string, int *value );
int get_s( int prompt, char *variable_name_string, char *value );
int get_prompt( int *value );


/**********************************************************************/
/* Prototypes for parallel I/O routine interfaces:
   io_ansi.c, io_t3d.c, io_paragon2.c, etc  */

FILE *g_open(const char *filename, const char *mode);
int g_seek(FILE *stream, off_t offset, int whence);
size_t g_write(const void *ptr, size_t size, size_t nmemb,FILE *stream);
size_t g_read(void *ptr, size_t size, size_t nmemb, FILE *stream);
int g_close(FILE *stream);

/**********************************************************************/
/* Prototypes for io_lat4.c routines also used in io_wb2.c */

void byterevn(type32 w[], int n);
void swrite_data(FILE* fp, void *src, size_t size, 
		 char *myname, char *descrip);
void pwrite_data(FILE* fp, void *src, size_t size, 
		 char *myname, char *descrip);
void pswrite_data(int parallel, FILE* fp, void *src, size_t size, 
		  char *myname, char *descrip);
int sread_data(FILE* fp, void *src, size_t size, 
	       char *myname, char *descrip);
int pread_data(FILE* fp, void *src, size_t size, 
	       char *myname, char *descrip);
int pread_byteorder(int byterevflag, FILE* fp, void *src, size_t size, 
		    char *myname, char *descrip);
int sread_byteorder(int byterevflag, FILE* fp, void *src, size_t size, 
		    char *myname, char *descrip);
int psread_data(int parallel, FILE* fp, void *src, size_t size, 
		char *myname, char *descrip);
int psread_byteorder(int byterevflag, int parallel, FILE* fp, 
		      void *src, size_t size, char *myname, char *descrip);


#endif /* _IO_LAT_H */
