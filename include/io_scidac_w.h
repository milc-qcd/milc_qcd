#ifndef _IO_SCIDAC_W_H
#define _IO_SCIDAC_W_H

#include <qio.h>
#include "../include/macros.h"

/**********************************************************************/
/* In io_scidac_w.c */

QIO_Reader *r_open_w_vector_scidac_file(char *filename, int serpar);
void r_close_usqcd_wprop_file(QIO_Reader *infile);
void r_close_w_vector_scidac_file(QIO_Reader *infile);
QIO_Reader *r_open_w_vector_scidac_file_xml(char *filename, int serpar,
					    QIO_String *xml_file);
QIO_Reader *r_open_w_vector_scidac_file(char *filename, int serpar);
void restore_w_vector_scidac_to_field(char *filename, int serpar,
			      wilson_vector *dest, int count);
void restore_w_vector_scidac_to_site(char *filename, int serpar,
			      field_offset dest, int count);
QIO_Reader *r_open_usqcd_wprop_file(char *filename, int serpar);
int read_w_vector_scidac(QIO_Reader *infile, wilson_vector *dest, int count);
int read_w_vector_scidac_xml(QIO_Reader *infile, wilson_vector *dest, 
			     int count, QIO_String *recxml);
int read_wproprecord_usqcd(QIO_Reader *infile, int *spin, int *color, 
			   wilson_vector *dest);
int read_wpropsource_C_usqcd(QIO_Reader *infile, char *srcinfo, int n,
			     complex *dest);
int read_wpropsource_D_usqcd(QIO_Reader *infile, char *srcinfo, int n,
			     wilson_vector *dest);
int save_w_vector_scidac(QIO_Writer *outfile, char *filename, char *recinfo,
			 int volfmt, wilson_vector *src, int count);
void save_w_vector_scidac_from_field(char *filename, char *filexml, 
				     char *recxml, int volfmt, int serpar, 
				     wilson_vector *src, int count);
void save_w_vector_scidac_from_site(char *filename, char *fileinfo,
				    char *recinfo, int volfmt, int serpar,
				    field_offset src, int count);
void w_close_usqcd_wprop_file(QIO_Writer *outfile);
QIO_Writer *w_open_usqcd_wprop_file(char *filename, int volfmt, 
    int serpar, int ildgstyle, char *stringLFN, int milc_type, char *fileinfo);
void w_close_w_vector_scidac_file(QIO_Writer *outfile);
QIO_Writer *w_open_w_vector_scidac_file(char *filename, char *fileinfo, 
					int volfmt, int serpar);
int write_prop_usqcd_sc(QIO_Writer *outfile, wilson_vector *src, int spin, 
			int color, char *recinfo);
int write_wpropsource_C_usqcd(QIO_Writer *outfile, char *srcinfo, complex *src, 
			     int t0);
int write_wpropsource_C_usqcd_xml(QIO_Writer *outfile, QIO_String *recxml, 
				  complex *src, int t0);
int write_wpropsource_D_usqcd(QIO_Writer *outfile, char *srcinfo, 
			     wilson_vector *src, int t0);
int write_wpropsource_D_usqcd_xml(QIO_Writer *outfile, QIO_String *recxml, 
				  wilson_vector *src, int t0);
/**********************************************************************/
/* In clover_info.c (application dependent) */

char *create_w_QCDML(void);
void free_w_QCDML(char *xml);

#endif /* _IO_SCIDAC_W_H */
