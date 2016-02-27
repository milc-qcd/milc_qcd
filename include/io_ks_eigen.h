#ifndef _IO_KS_EIGEN_H
#define _IO_KS_EIGEN_H

/************************* io_ks_eigen.h ******************************/

/* Macros, structures and prototypes for KS eigenvector I/O */

/* H. Ohno: 11/20/2014, derived from io_ksprop.h */

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
#include "../include/imp_ferm_links.h"
#include <stdio.h>

/* Used to create info file name */
#define ASCII_INFO_EXT ".info"

/* version numbers */
//#include "../include/file_types.h"

/**********************************************************************/

/* Begin definition of stuctures */

typedef gauge_header ks_eigen_header;
typedef gauge_check ks_eigen_check;

typedef struct {
  FILE *fp;                    /* File pointer */
  ks_eigen_header *header;     /* Pointer to header for file */
  char            *filename;   /* Pointer to file name string */
  int             byterevflag; /* Byte reverse flag - used only for reading */
  int             parallel;    /* 1 if file was opened in parallel
				  0 if serial */
  ks_eigen_check  check;       /* Checksum */
  int             Nvecs;        /* number of eigenvectors */
  double          *eigVal;     /* eigenvalues */
  double          *resid;      /* residual norm */
  //int             file_type;   /* File format */
  FILE            *info_fp;    /* Pointer to info file */
  //char            *info;       /* ASCII metadata */
  int              parity; /* parity to save/load */
} ks_eigen_file;

/**********************************************************************/
/* Info file format */

/* List of admissible keywords for ASCII KS eigenvector info file

   This list is intended to represent a union of the requirements
   for all applications.  A given application uses a subset.

   We keep this list to enforce some degree of consistency.
   Add more as needed.
*/
#ifdef CONTROL
char *ks_eigen_info_keyword[] = {
  //      "magic_number",
      "time_stamp",
      "checksums",
      "nx",
      "ny",
      "nz",
      "nt",
      "parity",
      "Nvecs",
      "eigVal",
      ""       /* Last entry MUST be a zero-length keyword */
};
#else
extern char *ks_eigen_info_keyword[];
#endif

/**********************************************************************/

/* Prototypes for io_ks_eigen.c */
void swrite_ks_eigen_hdr(FILE *fp, ks_eigen_header *kseigh);
int read_ks_eigen_hdr(ks_eigen_file *kseigf, int parallel);

int write_ks_eigen_info_item(FILE *fpout, char *keyword, char *fmt, char *src, int count,
			     int stride);
void write_ks_eigen_info_file(ks_eigen_file *kseigf);

ks_eigen_file *create_input_ks_eigen_file_handle(char *filename);
ks_eigen_file *create_output_ks_eigen_file_handle(void);
void destroy_ks_eigen_file_handle(ks_eigen_file *kseigf);

/**********************************************************************/

/* procedures for specific types of I/O */
ks_eigen_file *w_serial_ks_eigen_i(char *filename, int parity);
void w_serial_ks_eigen(ks_eigen_file *kseigf, int Nvecs, double *eigVal,
		       su3_vector **eigVec, double *resid);
void w_serial_ks_eigen_f(ks_eigen_file *kseigf);

ks_eigen_file *r_serial_ks_eigen_i(char *filename);
int r_serial_ks_eigen(ks_eigen_file *kseigf, int Nvecs, double *eigVal,
		      su3_vector **eigVec);
void r_serial_ks_eigen_f(ks_eigen_file *kseigf);

ks_eigen_file *w_ascii_ks_eigen_i(char *filename, int parity);
void w_ascii_ks_eigen(ks_eigen_file *kseigf, int Nvecs, double *eigVal,
		      su3_vector **eigVec, double *resid);
void w_ascii_ks_eigen_f(ks_eigen_file *kseigf);

ks_eigen_file *r_ascii_ks_eigen_i(char *filename);
int r_ascii_ks_eigen(ks_eigen_file *kseigf, int Nvecs, double *eigVal,
		     su3_vector **eigVec);
void r_ascii_ks_eigen_f(ks_eigen_file *kseigf);

/**********************************************************************/

/* Prototypes for io_helpers_ks_eigen.c */
void restore_eigVec(int Nvecs, double *eigVal, su3_vector **eigVec, int parity,
		    imp_ferm_links_t *fm);
ks_eigen_file *r_open_ks_eigen(int flag, char *filename);
ks_eigen_file *w_open_ks_eigen(int flag, char *filename, int parity);

void r_close_ks_eigen(int flag, ks_eigen_file *kseigf);
void w_close_ks_eigen(int flag, ks_eigen_file *kseigf);
int reload_ks_eigen(int flag, char *eigfile, int *Nvecs, double *eigVal,
		    su3_vector **eigVec, int timing);
int save_ks_eigen(int flag, char *savefile, int Nvecs, double *eigVal, 
		  su3_vector **eigVec, double *resid, int timing);

int convert_outflag_to_inflag_ks_eigen(int out_flag);

int ask_starting_ks_eigen(FILE *fp, int prompt, int *flag, char *filename );
int ask_ending_ks_eigen(FILE *fp, int prompt, int *flag, char *filename );

#endif /* _IO_KS_EIGEN_H */
