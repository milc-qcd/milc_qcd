#ifndef _IO_WB_H
#define _IO_WB_H
/************************ io_wb.h *************************************/
/* This header file defines the binary file format for the propagator file and
   defines structures for file descriptors that include the file header information */

#include "../include/int32type.h"
#include "../include/macros.h"
#include "../include/file_types.h"
#include <stdio.h>

#ifdef CONTROL
#define EXTERN 
#else
#define EXTERN extern
#endif

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
#define MAX_SOURCE_SPINS 4

/* Begin definition of header stuctures */

/* Note that an effort is made to make all Realing point and integer
   fields 32 bits long.  However, byte ordering may vary across
   platforms, and no effort is made in writing the file to produce
   a standard byte order.  The input routines attempt to compensate
   for byte reversal automatically, by examining the magic number at
   the beginning of the file */

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

  int32type order;                      /* 0 file is in natural order
					  no coordinate list is attached.
				        1 file is in node-dump (checkpoint)
					  order.  Coordinate list is attached.
                                        2 file is in node-dump (checkpoint)
					  order but one file per node.
					  Coordinate list is attached
					  before each file. */

  int32type n_spins;                    /* Number of source spins in this file */
  int32type spins[MAX_SOURCE_SPINS];    /* List of source spin indices in file */

} w_prop_header;

/* 2. Site list (ONLY for checkpoint files - i.e. node-dump order files)

      A listing of site coordinates for the data in this file
      in the order of appearance.  The number of coordinates must
      be exactly nx*ny*nz*nt.  The site coordinate is encoded
      as nx*(ny*(nz*t + z) + y) + x in a 32-bit integer.   

      */

/* 3. Next comes a spin - color index and checksum  */

typedef struct {
  int32type spin;
  int32type color;
  u_int32type sum29;
  u_int32type sum31;
} w_prop_check;

/* 4. Finally the Wilson vectors appear */

/**********************************************************************/
/* Info file format */

/* List of admissible keywords for version 5 ASCII lattice info file */

#ifdef CONTROL
char *w_prop_info_keyword[] = {
      "magic_number",
      "time_stamp",
      "nx",
      "ny",
      "nz",
      "nt",
      "gauge.filename",
      "gauge.time_stamp",
      "gauge.checksums",
      "gauge.fix.description",
      "gauge.fix.tolerance",
      "gauge.fix.filename",
      "gauge.fix.time_stamp",
      "gauge.fix.checksums",
      "quark.description",
      "quark.kappa",
      "quark.clover.clov_c",
      "quark.clover.u0",
      "quark.boundary_condition",
      "source.description",
      "source.size",
      "source.x",
      "source.y",
      "source.z",
      "source.t",
      "source.n_spins",
      "source.spins",
      ""         /* Last entry MUST be a zero-length keyword */
};
#else
extern char *w_prop_info_keyword[];
#endif

/* Used to create info file name */

#define ASCII_W_PROP_INFO_EXT ".info"

/**********************************************************************/
/* 1996 Binary file format follows */
/* Kept for compatibility */

#define MAX_GAUGE_FIELD_DESCRIPT 200
#define MAX_GAUGE_FIELD_PARAM 2
#define MAX_DIRAC_DESCRIPT 200
#define MAX_DIRAC_PARAM 3
#define MAX_SOURCE_DESCRIPT 200
#define MAX_SOURCE_PARAM 2
#define IDENTITY_MAP -1
#define NO_MAP -2

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
  struct {                      /* Dirac operator parameters */
    int32type n_descript;          /* Number of bytes in character string */
    char   descript[MAX_DIRAC_DESCRIPT];        /* Describes Dirac operator */
    int32type n_param;             /* Number of Dirac operator parameters */
    Real  param[MAX_DIRAC_PARAM];              /* Dirac parameters */
  } dirac;
  struct {                      /* Source parameters */
    int32type n_descript;          /* Number of bytes in character string */
    char   descript[MAX_SOURCE_DESCRIPT];      /* Describes source */
    int32type n_param;             /* Number of source parameters */
    struct {                    /* Source parameters */
      int32type i1;
      Real  c1;
    } param;
    int32type n_spins;             /* Number of source spins in this file */
    int32type spins[MAX_SOURCE_SPINS];  /* List of source spin indices in file */
  } source;
} w_prop_header_1996 ;

/* 2. Parallel files only:

      Next comes a listing of site coordinates for the data in this file
      in the order of appearance.  The number of coordinates must
      be exactly nx*ny*nz*nt.  The site coordinate is encoded
      as nx*(ny*(nz*t + z) + y) + x in a 32-bit integer.   

      Serial files only:

      The site order of propagator elements is required to be in subscript
      order (x,y,z,t) with x varying most rapidly, followed by y, etc. 
      so this list is omitted.

      */

/*    Next, repeat Items 3 and 4 for each source spin and color */

/* 3. Next comes a check structure to introduce the propagator components
      for a given source spin and color and node number. */

EXTERN struct {
  int32type spin;
  int32type color;
  u_int32type checksum;
} w_prop_check_1996;

/* 4. Finally, the propagator Wilson vectors appear */

/*----------------------------------------------------------------------*/

/* File data structure */

typedef struct {
  FILE *         fp;            /* File pointer */
  w_prop_header* header;        /* Pointer to header for file */
  char *         filename;       /* Pointer to file name string */
  int            byterevflag;   /* Byte reverse flag - used only for reading */
  int32type *       rank2rcv;      /* File site list - used only for 
				   serial reading */ 
  int            parallel;      /* 0 if file was opened for serial reading
				   1 if opened for parallel reading */
  w_prop_check   check;         /* Current checksum, spin, color indices */
} w_prop_file;

/**********************************************************************/
/* Declarations for I/O routines in io_wb.c */

w_prop_file *r_ascii_w_i(char *filename);
int r_ascii_w(w_prop_file *wpf, int spin, int color, field_offset src);
void r_ascii_w_f(w_prop_file *wpf);

w_prop_file *r_serial_w_i(char *filename);
int r_serial_w(w_prop_file *wpf, int spin, int color, field_offset src);
void r_serial_w_f(w_prop_file *wpf);

w_prop_file *r_parallel_w_i(char *filename);
void r_parallel_w_o(w_prop_file *wpf);
int r_parallel_w(w_prop_file *wpf, int spin, int color, field_offset src);
void r_parallel_w_c(w_prop_file *wpf);
void r_parallel_w_f(w_prop_file *wpf);

w_prop_file *r_multidump_w_i(char *filename);
void r_multidump_w_o(w_prop_file *wpf);
int r_multidump_w(w_prop_file *wpf, int spin, int color, field_offset src);
void r_multidump_w_c(w_prop_file *wpf);
void r_multidump_w_f(w_prop_file *wpf);

w_prop_file *w_ascii_w_i(char *filename);
void w_ascii_w(w_prop_file *wpf, int spin, int color, field_offset src);
void w_ascii_w_f(w_prop_file *wpf);

w_prop_file *w_serial_w_i(char *filename);
void w_serial_w(w_prop_file *wpf, int spin, int color, field_offset src);
void w_serial_w_f(w_prop_file *wpf);

w_prop_file *w_parallel_w_i(char *filename);
void w_parallel_w_o(w_prop_file *wpf);
void w_parallel_w(w_prop_file *wpf, int spin, int color, field_offset src);
void w_parallel_w_c(w_prop_file *wpf);
void w_parallel_w_f(w_prop_file *wpf);

w_prop_file *w_checkpoint_w_i(char *filename);
void w_checkpoint_w_o(w_prop_file *wpf);
void w_checkpoint_w(w_prop_file *wpf, int spin, int color, field_offset src);
void w_checkpoint_w_c(w_prop_file *wpf);
void w_checkpoint_w_f(w_prop_file *wpf);

w_prop_file *w_multidump_w_i(char *filename);
void w_multidump_w_o(w_prop_file *wpf);
void w_multidump_w(w_prop_file *wpf, int spin, int color, field_offset src);
void w_multidump_w_c(w_prop_file *wpf);
void w_multidump_w_f(w_prop_file *wpf);

int write_w_prop_info_item( FILE *fpout,    /* ascii file pointer */
		       char *keyword,   /* keyword */
		       char *fmt,       /* output format -
					      must use s, d, f, or e */
		       char *src,       /* address of starting data */
		       int count,       /* number of data items if > 1 */
		       int stride);     /* byte stride of data if
                                           count > 1 */
/**********************************************************************/
/* In clover_info.c or wilson_info.c (application dependent) */
void write_appl_w_prop_info(FILE *fp);

/**********************************************************************/
/* Prototypes for io_helpers_w.c */
w_prop_file *r_open_prop(int flag, char *filename);
w_prop_file *w_open_prop(int flag, char *filename);
int reload_propagator( int flag, w_prop_file *wpf,
		       int spin, int color, field_offset dest, int timing);
void save_propagator( int flag, w_prop_file *wpf, 
		     int spin, int color, field_offset src, int timing);
int ask_starting_prop( int prompt, int *flag, char *filename );
int ask_ending_prop( int prompt, int *flag, char *filename );
void r_close_prop(int flag, w_prop_file *wpf);
void w_close_prop(int flag, w_prop_file *wpf);




#endif /* _IO_WB_H */
