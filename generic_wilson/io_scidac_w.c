/********************** io_scidac_w.c **********************************/
/* MILC/QIO interface for Wilson fields */
/* MIMD version 6 */
/* CD 2/2005
*/

#include "generic_wilson_includes.h"
#include <qio.h>
#include "../include/io_scidac.h"
#include "../include/io_scidac_w.h"
#include <string.h>

/* Factory function for moving Wilson vectors from the site structure
   to single precision output */
/* arg must point to the field_offset of the field */

void vget_F3_D_from_site(char *buf, size_t index, int count, void *arg)
{
  int c,d,j;
  /* Assume output vector is single precision */
  fwilson_vector *dest = (fwilson_vector *)buf;
  field_offset src = *((field_offset *)arg);
  site *s = &lattice[index];
  /* Source can be any precision */
  wilson_vector *src_vec = (wilson_vector *)F_PT(s,src);

  /* Copy, changing precision, if necessary */
  for(j = 0; j < count; j++)
    for(d = 0; d < 4; d++)for(c = 0; c < 3; c++){
      dest[j].d[d].c[c].real = src_vec[j].d[d].c[c].real;
      dest[j].d[d].c[c].imag = src_vec[j].d[d].c[c].imag;
    }
}

void vget_F3_D_from_field(char *buf, size_t index, int count, void *arg)
{
  int c,d,j;
  /* Assume output vector is single precision */
  fwilson_vector *dest = (fwilson_vector *)buf;
  wilson_vector *src = (wilson_vector *)arg;
  /* Source can be any precision */
  wilson_vector *src_vec = src + index * count;

  /* Copy, changing precision, if necessary */
  for(j = 0; j < count; j++)
    for(d = 0; d < 4; d++)for(c = 0; c < 3; c++){
      dest[j].d[d].c[c].real = src_vec[j].d[d].c[c].real;
      dest[j].d[d].c[c].imag = src_vec[j].d[d].c[c].imag;
    }
}

int write_F3_D_from_site(QIO_Writer *outfile, char *xml_write_lattice, 
			 field_offset src, int count)
{
  QIO_String *xml_record_out;
  int status;
  QIO_RecordInfo *rec_info;
  
  /* We assume output precision is single */
  char qdptype[] = "QDP_F3_DiracVector";
  char prec[] = "F";
  int datum_size = sizeof(fwilson_vector);
  int word_size = sizeof(float);

  /* Create the record info for the field */
  rec_info = QIO_create_record_info(QIO_FIELD, qdptype, prec, 3,
				    0, datum_size, count);
  /* Create the record XML for the field */
  xml_record_out = QIO_string_create();
  QIO_string_set(xml_record_out,xml_write_lattice);

  /* Write the record for the field */
  status = QIO_write(outfile, rec_info, xml_record_out, vget_F3_D_from_site, 
		     count*datum_size, word_size, (void *)&src);
  if(status != QIO_SUCCESS)return 1;

  QIO_destroy_record_info(rec_info);
  QIO_string_destroy(xml_record_out);

  return 0;
}

int write_F3_D_from_field(QIO_Writer *outfile, char *xml_write_lattice, 
			 wilson_vector *src, int count)
{
  QIO_String *xml_record_out;
  int status;
  QIO_RecordInfo *rec_info;
  
  /* We assume output precision is single */
  char qdptype[] = "QDP_F3_DiracVector";
  char prec[] = "F";
  int datum_size = sizeof(fwilson_vector);
  int word_size = sizeof(float);

  /* Create the record info for the field */
  rec_info = QIO_create_record_info(QIO_FIELD, qdptype, prec, 3,
				    0, datum_size, count);
  /* Create the record XML for the field */
  xml_record_out = QIO_string_create();
  QIO_string_set(xml_record_out,xml_write_lattice);

  /* Write the record for the field */
  status = QIO_write(outfile, rec_info, xml_record_out, vget_F3_D_from_field, 
		     count*datum_size, word_size, (void *)src);
  if(status != QIO_SUCCESS)return 1;

  QIO_destroy_record_info(rec_info);
  QIO_string_destroy(xml_record_out);

  return 0;
}

/* Factory function for moving single precision Wilson vector from input
   to site structure */
void vput_F3_D_to_site(char *buf, size_t index, int count, void *arg)
{
  int c,d,j;
  fwilson_vector *src = (fwilson_vector *)buf;
  field_offset dest = *((field_offset *)arg);
  site *s = &lattice[index];
  /* Destination can be any precision */
  wilson_vector *dest_vec = (wilson_vector *)F_PT(s,dest);
  
  /* Copy, changing precision, if necessary */
  for(j = 0; j < count; j++)
    for(d = 0; d<4; d++)for(c = 0; c < 3; c++){
      dest_vec[j].d[d].c[c].real = src[j].d[d].c[c].real;
      dest_vec[j].d[d].c[c].imag = src[j].d[d].c[c].imag;
    }
}

/* Factory function for moving single precision Wilson vector from input
   to site structure */
void vput_F3_D_to_field(char *buf, size_t index, int count, void *arg)
{
  int c,d,j;
  fwilson_vector *src = (fwilson_vector *)buf;
  wilson_vector *dest = (wilson_vector *)arg;
  /* Destination can be any precision */
  wilson_vector *dest_vec = dest + index * count;
  
  /* Copy, changing precision, if necessary */
  for(j = 0; j < count; j++)
    for(d = 0; d<4; d++)for(c = 0; c < 3; c++){
      dest_vec[j].d[d].c[c].real = src[j].d[d].c[c].real;
      dest_vec[j].d[d].c[c].imag = src[j].d[d].c[c].imag;
    }
}

/* Read a set of Wilson vectors */
int read_F3_D_to_site(QIO_Reader *infile, field_offset dest, int count)
{
  QIO_String *xml_record_in;
  QIO_RecordInfo rec_info;
  int status;
  /* We assume input precision is single */
  char qdptype[] = "QDP_F3_DiracVector";
  char prec[] = "F";
  int datum_size = sizeof(fwilson_vector);
  int word_size = sizeof(float);
  
  /* Read the field record */
  xml_record_in = QIO_string_create();
  status = QIO_read(infile, &rec_info, xml_record_in, 
		    vput_F3_D_to_site, count*datum_size, 
		    word_size, (void *)&dest);
  node0_printf("Record info \n\"%s\"\n",QIO_string_ptr(xml_record_in));
  if(status != QIO_SUCCESS)return 1;

  node0_printf("Checksums %x %x\n",
	       QIO_get_reader_last_checksuma(infile),
	       QIO_get_reader_last_checksumb(infile));

  QIO_string_destroy(xml_record_in);
  return 0;
}

/* Read a set of Wilson vectors */
int read_F3_D_to_field(QIO_Reader *infile, wilson_vector *dest, int count)
{
  QIO_String *xml_record_in;
  QIO_RecordInfo rec_info;
  int status;
  /* We assume input precision is single */
  char qdptype[] = "QDP_F3_DiracVector";
  char prec[] = "F";
  int datum_size = sizeof(fwilson_vector);
  int word_size = sizeof(float);
  
  /* Read the field record */
  xml_record_in = QIO_string_create();
  status = QIO_read(infile, &rec_info, xml_record_in, 
		    vput_F3_D_to_field, count*datum_size, 
		    word_size, (void *)dest);
  node0_printf("Record info \n\"%s\"\n",QIO_string_ptr(xml_record_in));
  if(status != QIO_SUCCESS)return 1;

  node0_printf("Checksums %x %x\n",
	       QIO_get_reader_last_checksuma(infile),
	       QIO_get_reader_last_checksumb(infile));

  QIO_string_destroy(xml_record_in);
  return 0;
}

/* Write Wilson vectors in SciDAC format, taking data from the site
   structure */
/* We don't have a MILC format for such a file */

void save_w_vector_scidac_from_site(char *filename, char *recxml, 
				   int volfmt, field_offset src, int count)
{
  QIO_Layout layout;
  QIO_Writer *outfile;
  int status;
  char *filexml;

  QIO_verbose(QIO_VERB_OFF);

  /* Build the layout structure */
  build_qio_layout(&layout);

  /* Create the file metadata */
  filexml = create_w_XML();

  /* Open file for writing */
  outfile = open_scidac_output(filename, volfmt, &layout, filexml);
  if(outfile == NULL)terminate(1);


  /* Write the lattice field */
  status = write_F3_D_from_site(outfile, recxml, src, count);
  if(status)terminate(1);
  
  /* Write information */
  if(volfmt == QIO_SINGLEFILE){
    node0_printf("Saved Wilson vector serially to binary file %s\n",
		 filename);
  }
  else if(volfmt == QIO_MULTIFILE){
    node0_printf("Saved Wilson vector as multifile to binary file %s\n",
	   filename);
  }
  else if(volfmt == QIO_PARTFILE){
    node0_printf("Saved Wilson vector in partition format to binary file %s\n",
	   filename);
  }

  node0_printf("Checksums %x %x\n",
	       QIO_get_writer_last_checksuma(outfile),
	       QIO_get_writer_last_checksumb(outfile));

  /* Close the file */
  QIO_close_write(outfile);

  free_w_XML(filexml);
}

/* Write Wilson vectors in SciDAC format, taking data from a field */
/* We don't have a MILC format for such a file */

void save_w_vector_scidac_from_field(char *filename, char *recxml, 
				   int volfmt, wilson_vector *src, int count)
{
  QIO_Layout layout;
  QIO_Writer *outfile;
  int status;
  char *filexml;

  QIO_verbose(QIO_VERB_OFF);

  /* Build the layout structure */
  build_qio_layout(&layout);

  /* Create the file metadata */
  filexml = create_w_XML();

  /* Open file for writing */
  outfile = open_scidac_output(filename, volfmt, &layout, filexml);
  if(outfile == NULL)terminate(1);


  /* Write the lattice field */
  status = write_F3_D_from_field(outfile, recxml, src, count);
  if(status)terminate(1);
  
  /* Write information */
  if(volfmt == QIO_SINGLEFILE){
    node0_printf("Saved Wilson vector serially to binary file %s\n",
		 filename);
  }
  else if(volfmt == QIO_MULTIFILE){
    node0_printf("Saved Wilson vector as multifile to binary file %s\n",
	   filename);
  }
  else if(volfmt == QIO_PARTFILE){
    node0_printf("Saved Wilson vector in partition format to binary file %s\n",
	   filename);
  }

  node0_printf("Checksums %x %x\n",
	       QIO_get_writer_last_checksuma(outfile),
	       QIO_get_writer_last_checksumb(outfile));

  /* Close the file */
  QIO_close_write(outfile);

  free_w_XML(filexml);
}

/* Read Wilson vectors in SciDAC format */

void restore_w_vector_scidac_to_site(char *filename, field_offset dest,
				     int count){
  QIO_Layout layout;
  QIO_Reader *infile;
  int status;

  QIO_verbose(QIO_VERB_OFF);

  /* Build the layout structure */
  build_qio_layout(&layout);

  /* Open file for reading */
  infile = open_scidac_input(filename, &layout);
  if(infile == NULL)terminate(1);

  /* Read the lattice field: one Wilson vector */
  status = read_F3_D_to_site(infile, dest, count);
  if(status)terminate(1);

  /* Close the file */
  QIO_close_read(infile);
}

/* Read Wilson vectors in SciDAC format */

void restore_w_vector_scidac_to_field(char *filename, 
		     wilson_vector *dest, int count){
  QIO_Layout layout;
  QIO_Reader *infile;
  int status;

  QIO_verbose(QIO_VERB_OFF);

  /* Build the layout structure */
  build_qio_layout(&layout);

  /* Open file for reading */
  infile = open_scidac_input(filename, &layout);
  if(infile == NULL)terminate(1);

  /* Read the lattice field: one Wilson vector */
  status = read_F3_D_to_field(infile, dest, count);
  if(status)terminate(1);

  /* Close the file */
  QIO_close_read(infile);
}

