#ifndef _IO_SCIDAC_H
#define _IO_SCIDAC_H

#include <qio.h>
#include "../include/macros.h"

/**********************************************************************/
/* In io_scidac.c */

void build_layout(QIO_Layout *layout);

QIO_Reader *open_input(char *filename, QIO_Layout *layout);

QIO_Writer *open_output(char *filename, int volfmt, QIO_Layout *layout,
			char *xml_write_file);

QIO_Writer *save_color_matrix_scidac(char *filename, char *filexml,
                       char *recxml, int volfmt, field_offset src, int count);

QIO_Writer *save_color_matrix_scidac_from_temp(char *filename, char *filexml,
		       char *recxml, int volfmt, su3_matrix *src, int count);

QIO_Reader *restore_color_matrix_scidac(char *filename, 
					field_offset dest, int count);
/**********************************************************************/
/* In gauge_info.c (application dependent) */

char *create_QCDML();
void free_QCDML(char *qcdml);

#endif /* _IO_SCIDAC_H */
