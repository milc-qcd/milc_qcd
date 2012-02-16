#ifndef _IO_SCIDAC_H
#define _IO_SCIDAC_H

#include <qio.h>
#include "../include/macros.h"

/**********************************************************************/
/* In io_scidac.c */

void build_qio_layout(QIO_Layout *layout);
void build_qio_filesystem(QIO_Filesystem *fs);
void close_scidac_input(QIO_Reader *infile);
void close_scidac_output(QIO_Writer *outfile);
QIO_Writer *open_scidac_output(char *filename, int volfmt, 
			       int serpar, int ildgtype, 
			       char *stringLFN, QIO_Layout *layout,
			       QIO_Filesystem *fs, QIO_String *filexml);
QIO_Reader *open_scidac_input(char *filename, QIO_Layout *layout, 
			      QIO_Filesystem *fs, int serpar);
QIO_Reader *open_scidac_input_xml(char *filename, QIO_Layout *layout,
				  QIO_Filesystem *fs, int serpar,
				  QIO_String *xml_file_in);
int qio_status(int status);
int read_lat_dim_scidac(char *filename, int *ndim, int dims[]);
void r_close_complex_scidac_file(QIO_Reader *infile);
QIO_Reader *r_open_complex_scidac_file_xml(char *filename, int serpar,
					   QIO_String *xml_file);
QIO_Reader *r_open_complex_scidac_file(char *filename, int serpar);
void restore_color_matrix_scidac_to_site(char *filename, 
					field_offset dest, int count);
void restore_color_matrix_scidac_to_field(char *filename, 
					su3_matrix *dest, int count);
void restore_random_state_scidac_to_site(char *filename, field_offset dest);

int read_lat_dim_scidac(char *filename, int *ndim, int dims[]);
void restore_real_scidac_to_field(char *filename, Real *dest, int count);

void restore_real_scidac_to_site(char *filename, field_offset dest, int count);

int read_complex_scidac(QIO_Reader *infile, complex *dest, int count);
void r_close_complex_scidac_file(QIO_Reader *infile);
void restore_complex_scidac_to_field(char *filename, int serpar,
				     complex *dest, int count);
void save_color_matrix_scidac_from_site(char *filename, char *filexml,
                       char *recxml, int volfmt, field_offset src, int count);

void save_color_matrix_scidac_from_field(char *filename, char *filexml,
		       char *recxml, int volfmt, su3_matrix *src, int count);
void save_random_state_scidac_from_site(char *filename, 
	char *filexml, char *recxml, int volfmt, field_offset src);

void save_real_scidac_from_field(char *filename, 
    char *filexml, char *recxml, int volfmt, Real *src, int count);

void save_real_scidac_from_site(char *filename, 
    char *filexml, char *recxml, int volfmt, field_offset src, int count);
void save_complex_scidac_from_field(char *filename, char *fileinfo,
				    char *recinfo, int volfmt, int serpar, 
				    complex *src, int count);
void w_close_complex_scidac_file(QIO_Writer *outfile);
QIO_Writer *w_open_complex_scidac_file(char *filename, char *fileinfo, 
				       int volfmt, int serpar);


/**********************************************************************/
/* In io_scidac_types.c */

int read_F_R_to_site(QIO_Reader *infile, QIO_String *xml_record_in, 
      field_offset dest, int count);
int read_F_C_to_site(QIO_Reader *infile, QIO_String *xml_record_in, 
      field_offset dest, int count);
int read_F3_V_to_site(QIO_Reader *infile, QIO_String *xml_record_in, 
      field_offset dest, int count);
int read_F3_M_to_site(QIO_Reader *infile, QIO_String *xml_record_in, 
      field_offset dest, int count);
int read_F3_D_to_site(QIO_Reader *infile, QIO_String *xml_record_in, 
      field_offset dest, int count);
int read_S_to_site(QIO_Reader *infile, QIO_String *xml_record_in, 
      field_offset dest);

int read_F_R_to_field(QIO_Reader *infile, QIO_String *xml_record_in, 
      Real *dest, int count);
int read_F_C_to_field(QIO_Reader *infile, QIO_String *xml_record_in, 
      complex *dest, int count);
int read_F3_V_to_field(QIO_Reader *infile, QIO_String *xml_record_in, 
      su3_vector *dest, int count);
int read_F3_M_to_field(QIO_Reader *infile, QIO_String *xml_record_in, 
      su3_matrix *dest, int count);
int read_F3_D_to_field(QIO_Reader *infile, QIO_String *xml_record_in, 
      wilson_vector *dest, int count);

int read_D_C_to_site(QIO_Reader *infile, QIO_String *xml_record_in, 
      field_offset dest, int count);
int read_D3_V_to_site(QIO_Reader *infile, QIO_String *xml_record_in, 
      field_offset dest, int count);
int read_D3_M_to_site(QIO_Reader *infile, QIO_String *xml_record_in, 
      field_offset dest, int count);

int read_D_C_to_field(QIO_Reader *infile, QIO_String *xml_record_in, 
      complex *dest, int count);
int read_D3_V_to_field(QIO_Reader *infile, QIO_String *xml_record_in, 
      su3_vector *dest, int count);
int read_D3_M_to_field(QIO_Reader *infile, QIO_String *xml_record_in, 
      su3_matrix *dest, int count);
int read_D3_D_to_field(QIO_Reader *infile, QIO_String *xml_record_in, 
      wilson_vector *dest, int count);


int write_F_R_from_site(QIO_Writer *outfile, 
		 QIO_String *xml_record_out, field_offset src, int count);
int write_F_C_from_site(QIO_Writer *outfile, 
		 QIO_String *xml_record_out, field_offset src, int count);
int write_F3_V_from_site(QIO_Writer *outfile, 
		 QIO_String *xml_record_out, field_offset src, int count);
int write_F3_D_from_site(QIO_Writer *outfile, 
		 QIO_String *xml_record_out, field_offset src, int count);
int write_F3_M_from_site(QIO_Writer *outfile, 
		 QIO_String *xml_record_out, field_offset src, int count);
int write_S_from_site(QIO_Writer *outfile, QIO_String *xml_record_out,
		      field_offset src);

int write_D_C_from_site(QIO_Writer *outfile, 
		 QIO_String *xml_record_out, field_offset src, int count);
int write_D3_V_from_site(QIO_Writer *outfile, 
		 QIO_String *xml_record_out, field_offset src, int count);
int write_D3_D_from_site(QIO_Writer *outfile, 
		 QIO_String *xml_record_out, field_offset src, int count);

int write_F_R_timeslice_from_site(QIO_Writer *outfile, 
	 QIO_String *xml_record_out, field_offset src, int count, int t0);
int write_F_C_timeslice_from_site(QIO_Writer *outfile, 
	 QIO_String *xml_record_out, field_offset src, int count, int t0);
int write_F3_V_timeslice_from_site(QIO_Writer *outfile, 
	 QIO_String *xml_record_out, field_offset src, int count, int t0);
int write_F3_D_timeslice_from_site(QIO_Writer *outfile, 
	 QIO_String *xml_record_out, field_offset src, int count, int t0);

int write_D_C_timeslice_from_site(QIO_Writer *outfile, 
	 QIO_String *xml_record_out, field_offset src, int count, int t0);
int write_D3_V_timeslice_from_site(QIO_Writer *outfile, 
	 QIO_String *xml_record_out, field_offset src, int count, int t0);
int write_D3_D_timeslice_from_site(QIO_Writer *outfile, 
	 QIO_String *xml_record_out, field_offset src, int count, int t0);

int write_F_R_from_field(QIO_Writer *outfile, 
		 QIO_String *xml_record_out,  Real *src, int count);
int write_F_C_from_field(QIO_Writer *outfile, 
		 QIO_String *xml_record_out,  complex *src, int count);
int write_F3_V_from_field(QIO_Writer *outfile, 
		  QIO_String *xml_record_out, su3_vector *src, int count);
int write_F3_M_from_field(QIO_Writer *outfile, 
		  QIO_String *xml_record_out, su3_matrix *src, int count);
int write_F3_D_from_field(QIO_Writer *outfile, 
		  QIO_String *xml_record_out, wilson_vector *src, int count);

int write_D_C_from_field(QIO_Writer *outfile, 
		 QIO_String *xml_record_out,  complex *src, int count);
int write_D3_V_from_field(QIO_Writer *outfile, 
		  QIO_String *xml_record_out, su3_vector *src, int count);
int write_D3_D_from_field(QIO_Writer *outfile, 
		  QIO_String *xml_record_out, wilson_vector *src, int count);


int write_F_R_timeslice_from_field(QIO_Writer *outfile, 
		 QIO_String *xml_record_out,  Real *src, int count, int t0);
int write_F_C_timeslice_from_field(QIO_Writer *outfile, 
		 QIO_String *xml_record_out,  complex *src, int count, int t0);
int write_F3_V_timeslice_from_field(QIO_Writer *outfile, 
	 QIO_String *xml_record_out,  su3_vector *src, int count, int t0);
int write_F3_D_timeslice_from_field(QIO_Writer *outfile, 
	  QIO_String *xml_record_out, wilson_vector *src, int count, int t0);

int write_D_C_timeslice_from_field(QIO_Writer *outfile, 
		 QIO_String *xml_record_out,  complex *src, int count, int t0);
int write_D3_V_timeslice_from_field(QIO_Writer *outfile, 
	  QIO_String *xml_record_out, su3_vector *src, int count, int t0);
int write_D3_D_timeslice_from_field(QIO_Writer *outfile, 
	  QIO_String *xml_record_out, wilson_vector *src, int count, int t0);

/**********************************************************************/
/* In gauge_info.c (application dependent) */

char *create_QCDML(void);
void free_QCDML(char *qcdml);

#endif /* _IO_SCIDAC_H */
