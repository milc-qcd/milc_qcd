#ifndef _IO_SCIDAC_KS_H
#define _IO_SCIDAC_KS_H

#include <qio.h>
#include "../include/macros.h"

/**********************************************************************/
/* In io_scidac_ks.c */

void save_ks_vector_scidac_from_site(char *filename, 
				     char *filexml, int volfmt, int serpar, 
				     field_offset src, int count);
void save_ks_vector_scidac_from_field(char *filename, char *recxml, 
				      int volfmt, int serpar, 
				      su3_vector *src, int count);

void restore_ks_vector_scidac_to_site(char *filename, field_offset dest,
				      int serpar, int count);
void restore_ks_vector_scidac_to_field(char *filename, su3_vector *dest, 
				       int serpar, int count);
/**********************************************************************/
/* In ks_info.c (application dependent) */

char *create_ks_XML();
void free_ks_XML(char *xml);

#endif /* _IO_SCIDAC_KS_H */
