#ifndef _IO_SCIDAC_H
#define _IO_SCIDAC_H

#include <qio.h>
#include "../include/macros.h"

/**********************************************************************/
/* In io_scidac.c */

void build_qio_layout(QIO_Layout *layout);

QIO_Writer *open_scidac_output(char *filename, int volfmt, 
			       int serpar, int ildgtype, 
			       char *stringLFN, QIO_Layout *layout,
			       char *xml_write_file);

QIO_Reader *open_scidac_input(char *filename, QIO_Layout *layout, 
			      int serpar);

QIO_Reader *open_input(char *filename, QIO_Layout *layout);

QIO_Writer *open_output(char *filename, int volfmt, QIO_Layout *layout,
			char *xml_write_file);
void close_input(QIO_Reader *infile);
void close_output(QIO_Writer *outfile);

int read_lat_dim_scidac(char *filename, int *ndim, int dims[]);

void save_color_matrix_scidac_from_site(char *filename, char *filexml,
                       char *recxml, int volfmt, field_offset src, int count);

void save_color_matrix_scidac_from_field(char *filename, char *filexml,
		       char *recxml, int volfmt, su3_matrix *src, int count);

void restore_color_matrix_scidac_to_site(char *filename, 
					field_offset dest, int count);
void restore_color_matrix_scidac_to_field(char *filename, 
					su3_matrix *dest, int count);
/**********************************************************************/
/* In gauge_info.c (application dependent) */

char *create_QCDML();
void free_QCDML(char *qcdml);

#endif /* _IO_SCIDAC_H */
