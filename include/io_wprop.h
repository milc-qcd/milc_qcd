#ifndef _IO_WB_H
#define _IO_WB_H
/************************ io_wprop.h *************************************/
/* This header file defines the binary file format for the propagator file and
   defines structures for file descriptors that include the file header information */

#include "../include/int32type.h"
#include "../include/macros.h"
#include "../include/file_types.h"
#include "../include/generic_wilson.h"
#include <stdio.h>
#include <time.h>

#ifdef CONTROL
#define EXTERN 
#else
#define EXTERN extern
#endif

#ifdef HAVE_QIO
#include <qio.h>
#endif

/**********************************************************************/
/* Binary lattice formats  (deprecated)                               */
/**********************************************************************/
/* In version 5 we have two binary lattice file formats:
      serial files     Data written in coordinate natural order
      checkpoint files Data written in node dump order

      Further descriptive information is kept in a separate ASCII header
      file.  See below.

*/


/*--------------------------------------------------------------------*/
/* version 5 binary file format (deprecated) */

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
  time_t t_stamp;
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
  int32type elements_per_site;     /* For Fermilab propagator files */

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
      "header.time_stamp",
      "header.size_of_element",
      "header.elements_per_site",
      "header.dim[0]",
      "header.dim[1]",
      "header.dim[2]",
      "header.dim[3]",
      "header.site_order",
      "quark.description",
      "quark.kappa",
      "quark.clover.clov_c",
      "quark.clover.u0",
      "quark.boundary_condition",
      "record[0].checksum",
      "record[1].checksum",
      "record[2].checksum",
      "record[3].checksum",
      "record[4].checksum",
      "record[5].checksum",
      "record[6].checksum",
      "record[7].checksum",
      "record[8].checksum",
      "record[9].checksum",
      "record[10].checksum",
      "record[11].checksum",
      "source.description",
      "source.filename",
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

#define ASCII_INFO_EXT ".info"

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

/* 4. Finally, the propagator Wilson vectors appear */

/*----------------------------------------------------------------------*/

/* File data structure */

typedef struct {
  FILE *         fp;            /* File pointer */
  w_prop_header* header;        /* Pointer to header for file */
  char *         filename;      /* Pointer to file name string */
  int            byterevflag;   /* Byte reverse flag - used only for reading */
  int            parallel;      /* 0 if file was opened for serial reading
				   1 if opened for parallel reading */
  w_prop_check   check;         /* Current checksum, spin, color indices */
  wilson_propagator *prop;      /* Cache for the full propagator if needed */
  int            file_type;     /* File format */
  FILE *         info_fp;       /* File pointer for info file */
#ifdef HAVE_QIO
  QIO_Reader *   infile;        /* For QIO reading */
  QIO_Writer *   outfile;       /* For QIO writing */
#endif
} w_prop_file;

/**********************************************************************/
/* Declarations for I/O routines in io_source_w_fm.c */
void r_source_w_fm_to_site(char *filename, field_offset dest_site,
			   int spin, int color, int x0, int y0, int z0, int t0);
void r_source_w_fm_to_field(char *filename, wilson_vector *dest_field,
			    int spin, int color, int x0, int y0, int z0, int t0);

/* Declarations for I/O routines in io_prop_w.c */

w_prop_file *r_ascii_w_i(char *filename);
int r_ascii_w(w_prop_file *wpf, int spin, int color, field_offset src);
void r_ascii_w_f(w_prop_file *wpf);

w_prop_file *r_serial_w_i(char *filename);
int r_serial_w_to_site(w_prop_file *wpf, int spin, int color, 
		       field_offset dest_site);
int r_serial_w_to_field(w_prop_file *wpf, int spin, int color, 
			wilson_vector *dest_field);
void r_serial_w_f(w_prop_file *wpf);

w_prop_file *r_parallel_w_i(char *filename);
void r_parallel_w_o(w_prop_file *wpf);
int r_parallel_w_to_site(w_prop_file *wpf, int spin, int color, 
			 field_offset dest_site);
int r_parallel_w_to_field(w_prop_file *wpf, int spin, int color, 
			  wilson_vector *dest_field);
void r_parallel_w_c(w_prop_file *wpf);
void r_parallel_w_f(w_prop_file *wpf);

w_prop_file *r_multidump_w_i(char *filename);
void r_multidump_w_o(w_prop_file *wpf);
int r_multidump_w_to_site(w_prop_file *wpf, int spin, int color, 
			  field_offset dest_site);
int r_multidump_w_to_field(w_prop_file *wpf, int spin, int color, 
			   wilson_vector *dest_field);
void r_multidump_w_c(w_prop_file *wpf);
void r_multidump_w_f(w_prop_file *wpf);

w_prop_file *setup_input_w_prop_file(char *filename);
w_prop_file *setup_output_w_prop_file(void);
void clear_input_w_prop_file(w_prop_file *wpf);
void clear_output_w_prop_file(w_prop_file *wpf);

w_prop_file *w_ascii_w_i(char *filename);
void w_ascii_w(w_prop_file *wpf, int spin, int color, 
			 wilson_vector *src);
void w_ascii_w_f(w_prop_file *wpf);

w_prop_file *w_serial_w_i(char *filename);
void w_serial_w_from_site(w_prop_file *wpf, int spin, int color, 
			  field_offset src_site);
void w_serial_w_from_field(w_prop_file *wpf, int spin, int color, 
			   wilson_vector *src_field);
void w_serial_w_f(w_prop_file *wpf);

w_prop_file *w_parallel_w_i(char *filename);
void w_parallel_w_o(w_prop_file *wpf);
void w_parallel_w_from_site(w_prop_file *wpf, int spin, int color, 
			    field_offset src_site);
void w_parallel_w_from_field(w_prop_file *wpf, int spin, int color, 
			     wilson_vector *src_field);
void w_parallel_w_c(w_prop_file *wpf);
void w_parallel_w_f(w_prop_file *wpf);

w_prop_file *w_checkpoint_w_i(char *filename);
void w_checkpoint_w_o(w_prop_file *wpf);
void w_checkpoint_w_from_site(w_prop_file *wpf, int spin, int color, 
			      field_offset src_site);
void w_checkpoint_w_from_field(w_prop_file *wpf, int spin, int color, 
			       wilson_vector *src_field);
void w_checkpoint_w_c(w_prop_file *wpf);
void w_checkpoint_w_f(w_prop_file *wpf);

w_prop_file *w_multidump_w_i(char *filename);
void w_multidump_w_o(w_prop_file *wpf);
void w_multidump_w_from_site(w_prop_file *wpf, int spin, int color, 
			     field_offset src_site);
void w_multidump_w_from_field(w_prop_file *wpf, int spin, int color, 
			      wilson_vector *src_field);
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
int sprint_w_prop_info_item( 
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
/**********************************************************************/
/* In clover_info.c or wilson_info.c (application dependent) */
void write_appl_w_prop_info(FILE *fp);

/**********************************************************************/
/* Prototypes for io_helpers_w.c */

void interpret_usqcd_w_save_flag(int *volfmt, int *serpar, int flag);
int read_lat_dim_wprop(char *filename, int file_type, int *ndim, int dims[]);
w_prop_file *r_open_wprop(int flag, char *filename);
w_prop_file *w_open_wprop(int flag, char *filename, int source_type);
int reload_wprop_sc_to_field( int flag, w_prop_file *wpf, 
			      quark_source *wqs, int spin, int color, 
			      wilson_vector *dest, int timing);
int save_wprop_sc_from_field( int flag, w_prop_file *wpf, 
			      quark_source *wqs, int spin, int color, 
			      wilson_vector *src, char *recinfo, int timing);
int reload_wprop_to_field( int flag, char *filename, quark_source *wqs,
			   wilson_propagator *dest, int timing);
int reload_wprop_to_wp_field( int flag, char *filename, 
			      quark_source *wqs,
			      wilson_prop_field *dest, int timing);
int save_wprop_from_field( int flag, char *filename, quark_source *wqs,
			   wilson_propagator *src, char *recxml, int timing);
int save_wprop_from_wp_field( int flag, char *filename, 
			      quark_source *wqs,
			      wilson_prop_field *src, char *recxml, int timing);
int reload_wprop_sc_to_site( int flag, w_prop_file *wpf,
			     quark_source *wqs, int spin, int color, 
			     field_offset dest, int timing);
int save_wprop_sc_from_site( int flag, w_prop_file *wpf, 
			     quark_source *wqs, int spin, int color, 
			     field_offset src, int timing);
int reload_wprop_to_site( int flag, char *filename, quark_source *wqs,
			  field_offset dest, int timing);
int save_wprop_from_site( int flag, char *filename, quark_source *wqs,
			  field_offset src, char *recxml, int timing);
void r_close_wprop(int flag, w_prop_file *wpf);
void w_close_wprop(int flag, w_prop_file *wpf);
int ask_starting_wprop( FILE *fp, int prompt, int *flag, char *filename );
int ask_ending_wprop( FILE *fp, int prompt, int *flag, char *filename );
int convert_outflag_to_inflag_wprop(int out_flag);

/* Prototpyes for io_prop_w_fm.c */
w_prop_file *r_serial_w_fm_i(char *filename);
w_prop_file *r_serial_w_fm_sc_i(char *filename);
void r_serial_w_fm(w_prop_file *wpf, field_offset dest_site, 
		   wilson_propagator *dest_field);
void r_serial_w_fm_to_field(w_prop_file *wpf, wilson_propagator *dest_field);
void r_serial_w_fm_to_site(w_prop_file *wpf, field_offset dest_site);
void r_serial_w_fm_f(w_prop_file *wpf);
void r_prop_w_fm_to_site(char *filename, field_offset dest);
void r_prop_w_fm_to_field(char *filename, wilson_propagator *dest);

w_prop_file *w_serial_w_fm_i(char *filename);
w_prop_file *w_serial_w_fm_sc_i(char *filename);
void w_serial_w_fm_from_field(w_prop_file *wpf, wilson_propagator *src_field);
void w_serial_w_fm_from_site(w_prop_file *wpf, field_offset src_site);
void w_serial_w_fm_f(w_prop_file *wpf);

#endif /* _IO_WB_H */
