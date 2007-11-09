#ifndef _IO_SCIDAC_W_H
#define _IO_SCIDAC_W_H

#include <qio.h>
#include "../include/macros.h"

/**********************************************************************/
/* In io_scidac_w.c */

void save_w_vector_scidac_from_site(char *filename, char *filexml, 
	   int volfmt, int serpar, field_offset src, int count);
void save_w_vector_scidac_from_field(char *filename, char *filexml, 
	    int volfmt, int serpar, wilson_vector *src, int count);

void restore_w_vector_scidac_to_site(char *filename, int serpar,
			      field_offset dest, int count);
void restore_w_vector_scidac_to_field(char *filename, int serpar,
			      wilson_vector *dest, int count);
int write_prop_usqcd_sc(QIO_Writer *outfile, wilson_vector *src, int spin, 
			int color, char *recinfo);
void close_usqcd_prop_write(QIO_Writer *outfile);
int io_detect_w_usqcd(char *filename);
QIO_Reader *open_usqcd_prop_read(char *filename, int serpar);
QIO_Writer *open_usqcd_prop_write(char *filename, int volfmt, 
    int serpar, int ildgstyle, char *stringLFN, int milc_type, char *fileinfo);
void read_propsource_C_usqcd(QIO_Reader *infile, char *srcinfo, int n,
			     complex *dest);
void read_propsource_D_usqcd(QIO_Reader *infile, char *srcinfo, int n,
			     wilson_vector *dest);
int write_propsource_C_usqcd(QIO_Writer *outfile, char *srcinfo, complex *src, 
			     int t0);
int write_propsource_D_usqcd(QIO_Writer *outfile, char *srcinfo, 
			     wilson_vector *src, int t0);
void read_proprecord_usqcd(QIO_Reader *infile, int *spin, int *color, 
			   wilson_vector *dest);
void close_usqcd_prop_read(QIO_Reader *infile);
/**********************************************************************/
/* In clover_info.c (application dependent) */

char *create_w_QCDML();
void free_w_QCDML(char *xml);

#endif /* _IO_SCIDAC_W_H */
