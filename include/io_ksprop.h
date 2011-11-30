#ifndef _IO_PROP_KS_H
#define _IO_PROP_KS_H

/************************* io_ksprop.h ******************************/

/* Macros, structures and prototypes for KS propagator I/O */

/* Original by MBW: Feb 2001, derived from io_lat.h 
    only implementing ASCII I/O for now 
    April 2002, MBW: adding binary I/O .. just serial
*/

/* For portability we require a 32-bit integer type for writing
   integer data to binary files.  Machine-dependent. Add more as
   needed. Compile with -DSHORT32 for 64-bit integer machines.
   like the T3E */

/* The read routines try to accommodate byte reversal */
/* No provision for 64-bit integers here - there should be no such files! */

#ifdef CONTROL
#define EXTERN 
#else
#define EXTERN extern
#endif

#include "../include/int32type.h"
#include "../include/macros.h"
#include "../include/io_lat.h"
#include "../include/generic.h"
//#include "../include/generic_ks.h"
#include <stdio.h>

#ifdef HAVE_QIO
#include <qio.h>
#endif

/* Used to create info file name */

#define ASCII_INFO_EXT ".info"

/* version numbers */
#include "../include/file_types.h"

/* Begin definition of stuctures */

typedef struct {
  int32type magic_number;               /* Identifies file format */
  int32type gmtime_stamp;             /* Used in FNAL header from call to time*/
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
} ks_prop_header;

typedef struct {
  int32type color;
  u_int32type sum31;
  u_int32type sum29;
} ks_prop_check;


/**********************************************************************/
typedef struct {
  FILE *           fp;            /* File pointer */
  ks_prop_header * header;        /* Pointer to header for file */
  char *           filename;      /* Pointer to file name string */
  int              byterevflag;   /* Byte reverse flag - used only for reading */
  int              parallel;      /* 1 if file was opened in parallel
				     0 if serial */
  ks_prop_check    check;         /* Checksum */
  su3_vector *     prop;          /* If we have to read the data in one lump*/
  int              file_type;     /* File format */
  FILE *           info_fp;       /* Pointer to info file */
  char *           info;          /* ASCII metadata */
#ifdef HAVE_QIO
  QIO_Reader *     infile;       /* For QIO reading */
  QIO_Writer *     outfile;      /* For QIO writing */
#endif
} ks_prop_file;

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
char *ks_prop_info_keyword[] = {
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
      "gauge.aniso.xi",
      "gauge.tadpole.u0",
      "gauge.tadpole.u0s",
      "gauge.tadpole.u0t",
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
      "quark.u0",
      "quark.lightspeed.c0",

      "ksprop.dim",
      "ksprop.magic_number",
      "ksprop.size-of-element",
      "ksprop.elements-per-site",
      "ksprop.site-order",
      "ksprop.time_stamp",
      ""       /* Last entry MUST be a zero-length keyword */
};
#else
extern char *ks_prop_info_keyword[];
#endif

/**********************************************************************/
/* Prototypes for io_prop_ks.c */

int write_ksprop_info_item( FILE *fpout, /* ascii file pointer */
			    char *keyword,   /* keyword */
			    char *fmt,       /* output format -
						must use s, d, f, or e */
			    char *src,       /* address of starting data */
			    int count,       /* number of data items if > 1 */
			    int stride);     /* byte stride of data if
						count > 1 */

ks_prop_file *create_input_ksprop_file_handle(char *filename);
ks_prop_file *create_output_ksprop_file_handle(void);
void destroy_ksprop_file_handle(ks_prop_file *kspf);

/**********************************************************************/
/* In ksprop_info.c (application dependent) */
void write_appl_ksprop_info(FILE *fp);

/**********************************************************************/
/* procedures for specific types of I/O                               */

/* 

   restore_ksprop_ascii: Node 0 opens and reads in prop in ASCII format
   save_ksprop_ascii:    Node 0 opens and writes prop in ASCII format

   w_serial_ks_i     Node 0 opens serial file for writing and writes header
   w_serial_ks       Node 0 writes propagator to the specified serial file
   w_serial_ks_f     Closes the file

   r_serial_ks_i     Node 0 Opens serial file for reading and reads header
   r_serial_ks       Node 0 reads propagator from specified serial file
   r_serial_ks_f     Closes the file

*/

ks_prop_file *r_serial_ks_i(char *filename);
int r_serial_ks_to_site(ks_prop_file *kspf, int color, field_offset src);
int r_serial_ks_to_field(ks_prop_file *kspf, int color, su3_vector *src);
void r_serial_ks_f(ks_prop_file *kspf);

ks_prop_file *w_serial_ks_i(char *filename);
void w_serial_ks_from_site(ks_prop_file *kspf, int color, field_offset src);
void w_serial_ks_from_field(ks_prop_file *kspf, int color, su3_vector *src);
void w_serial_ks_f(ks_prop_file *kspf);

void w_serial_ksprop_tt(char *filename, field_offset prop);
void w_ascii_ksprop_tt(char *filename, field_offset prop);

ks_prop_file *r_ascii_ks_i(char *filename);
int r_ascii_ks(ks_prop_file *kspf, int color, su3_vector *dest_field);
void r_ascii_ks_f(ks_prop_file *kspf);

ks_prop_file *w_ascii_ks_i(char *filename);
void w_ascii_ks(ks_prop_file *kspf, int color, su3_vector *src_field);
void w_ascii_ks_f(ks_prop_file *kspf);


/*** DEPRECATED: ***/
ks_prop_file *restore_ksprop_ascii(char *filename, field_offset prop);
ks_prop_file *save_ksprop_ascii(char *filename, field_offset prop);

/**********************************************************************/
/* Prototypes for io_prop_ks_fm.c */

void swrite_ks_fm_prop_hdr(FILE *fp, ks_prop_header *ksph);
void write_ks_fmprop_info_file(ks_prop_file *pf);
ks_prop_file *create_input_ks_fmprop_file_handle(char *filename);
ks_prop_file *create_output_ks_fmprop_file_handle(void);
ks_prop_file *w_serial_ks_fm_i(char *filename);
void w_serial_ks_fm_from_site(ks_prop_file *kspf, field_offset prop);
void w_serial_ks_fm_from_field(ks_prop_file *kspf, su3_vector *prop);
void w_serial_ks_fm_f(ks_prop_file *kspf);
ks_prop_file *r_serial_ks_fm_i(char *filename);
int r_serial_ks_fm_to_site(ks_prop_file *kspf, field_offset prop);
int r_serial_ks_fm_to_field(ks_prop_file *kspf, su3_vector *prop);
void r_serial_ks_fm_f(ks_prop_file *kspf);

/* Prototypes for io_helpers_ks.c */

void interpret_usqcd_ks_save_flag(int *volfmt, int *serpar, int flag);
int read_lat_dim_ksprop(char *filename, int file_type, int *ndim, int dims[]);

ks_prop_file *r_open_ksprop(int flag, char *filename);

void r_close_ksprop(int flag, ks_prop_file *kspf);
ks_prop_file *w_open_ksprop(int flag, char *filename, int source_type);
void w_close_ksprop(int flag, ks_prop_file *kspf);
int reload_ksprop_to_field3( int flag, char *filename, quark_source *ksqs,
			     su3_vector *dest, int timing);
int reload_ksprop_to_ksp_field( int flag, char *filename, 
				quark_source *ksqs,
				ks_prop_field *dest, int timing);
int reload_ksprop_to_site3( int flag, char *filename, quark_source *ksqs,
			   field_offset dest, int timing);
int reload_ksprop_c_to_field( int flag, ks_prop_file *kspf, 
			      quark_source *ksqs, int color, 
			      su3_vector *dest, int timing);
void save_ksprop_from_field3( int flag, char *filename, char *recxml, 
			      quark_source *ksqs,
			      su3_vector *src, int timing);
void save_ksprop_from_ksp_field( int flag, char *filename, char *recxml, 
				 quark_source *ksqs,
				 ks_prop_field *src, int timing);
void save_ksprop_from_site3( int flag, char *filename, char *recxml, 
			     quark_source *ksqs,
			     field_offset src, int timing);
int save_ksprop_c_from_field( int flag, ks_prop_file *kspf, 
			      quark_source *ksqs,
			      int color, su3_vector *src, 
			      char *recinfo, int timing);
int ask_starting_ksprop( FILE *fp, int prompt, int *flag, char *filename );
int ask_ending_ksprop( FILE *fp, int prompt, int *flag, char *filename );
int convert_outflag_to_inflag_ksprop(int out_flag);

/**********************************************************************/
/* Declarations for I/O routines in io_source_ks_fm.c */

ks_fm_source_file * r_source_ks_fm_i(char *filename);
void r_source_ks_fm(ks_fm_source_file *kssf, su3_vector *dest_field,
		    int x0, int y0, int z0, int t0);
void r_source_ks_fm_f(ks_fm_source_file *kssf);

#endif /* _IO_PROP_KS_H */
