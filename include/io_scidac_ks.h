#ifndef _IO_SCIDAC_KS_H
#define _IO_SCIDAC_KS_H

#include <qio.h>
#include "../include/macros.h"

/**********************************************************************/
/* In io_scidac_ks.c */

QIO_Writer *save_ks_vector_scidac(char *filename, char *filexml, int volfmt, 
				  field_offset src, int count);
QIO_Writer *save_ks_vector_scidac_from_temp(char *filename, char *recxml, 
			    int volfmt, su3_vector *src, int count);

QIO_Reader *restore_ks_vector_scidac(char *filename, field_offset dest,
				     int count);
QIO_Reader *restore_ks_vector_scidac_to_temp(char *filename, 
					     su3_vector *dest, int count);
/**********************************************************************/
/* In ks_info.c (application dependent) */

char *create_ks_XML();
void free_ks_XML(char *xml);

#endif /* _IO_SCIDAC_KS_H */
