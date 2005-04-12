#ifndef _IO_SCIDAC_W_H
#define _IO_SCIDAC_W_H

#include <qio.h>
#include "../include/macros.h"

/**********************************************************************/
/* In io_scidac_w.c */

void save_w_vector_scidac_from_site(char *filename, char *filexml, 
	   int volfmt, field_offset src, int count);
void save_w_vector_scidac_from_field(char *filename, char *filexml, 
	    int volfmt, wilson_vector *src, int count);

void restore_w_vector_scidac_to_site(char *filename, 
			      field_offset dest, int count);
void restore_w_vector_scidac_to_field(char *filename, 
			      wilson_vector *dest, int count);
int read_F3_D_to_site(QIO_Reader *infile, field_offset dest, int count);
int read_F3_D_to_field(QIO_Reader *infile, wilson_vector *dest, int count);
int write_F3_D_from_site(QIO_Writer *outfile, char *xml_write_lattice, 
	       field_offset src, int count);
int write_F3_D_from_field(QIO_Writer *outfile, char *xml_write_lattice, 
	       wilson_vector *src, int count);

/**********************************************************************/
/* In clover_info.c (application dependent) */

char *create_w_XML();
void free_w_XML(char *xml);

#endif /* _IO_SCIDAC_W_H */
