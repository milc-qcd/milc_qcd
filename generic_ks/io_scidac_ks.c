/********************** io_scidac_ks.c **********************************/
/* MILC/QIO interface for KS fields */
/* MIMD version 6 */
/* CD 2/2005
*/

#include "generic_ks_includes.h"
#include <qio.h>
#include "../include/io_scidac.h"
#include "../include/io_scidac_ks.h"
#include <string.h>

/* Factory function for moving color vectors from the site structure
   to single precision output */
/* arg must point to the field_offset of the field */

void vget_F3_V(char *buf, size_t index, int count, void *arg)
{
  int i,j;
  /* Assume output vector is single precision */
  fsu3_vector *dest = (fsu3_vector *)buf;
  field_offset src = *((field_offset *)arg);
  site *s = &lattice[index];
  /* Source can be any precision */
  su3_vector *src_vec = (su3_vector *)F_PT(s,src);

  /* Copy, changing precision, if necessary */
  for(j = 0; j < count; j++)
    for(i = 0; i < 3; i++){
      dest[j].c[i].real = src_vec[j].c[i].real;
      dest[j].c[i].imag = src_vec[j].c[i].imag;
    }
}

/* Factory function for moving color vectors from the site structure
   to single precision output */
/* arg must point to the field_offset of the field */

void vget_F3_V_from_temp(char *buf, size_t index, int count, void *arg)
{
  int i,j;
  /* Assume output vector is single precision */
  fsu3_vector *dest = (fsu3_vector *)buf;
  /* Source can be any precision */
  su3_vector *src = (su3_vector *)arg;
  su3_vector *src_vec = src + index * count;

  /* Copy, changing precision, if necessary */
  for(j = 0; j < count; j++)
    for(i = 0; i < 3; i++){
      dest[j].c[i].real = src_vec[j].c[i].real;
      dest[j].c[i].imag = src_vec[j].c[i].imag;
    }
}

int write_F3_V(QIO_Writer *outfile, char *xml_write_lattice, 
		    field_offset src, int count)
{
  QIO_String *xml_record_out;
  int status;
  QIO_RecordInfo *rec_info;
  
  /* We assume output precision is single */
  char qdptype[] = "QDP_F3_ColorVector";
  char prec[] = "F";
  int datum_size = sizeof(fsu3_vector);
  int word_size = sizeof(float);

  /* Create the record info for the field */
  rec_info = QIO_create_record_info(QIO_FIELD, qdptype, prec, 3,
				    0, datum_size, count);
  /* Create the record XML for the field */
  xml_record_out = QIO_string_create();
  QIO_string_set(xml_record_out,xml_write_lattice);

  /* Write the record for the field */
  status = QIO_write(outfile, rec_info, xml_record_out, vget_F3_V, 
		     count*datum_size, word_size, (void *)&src);
  if(status != QIO_SUCCESS)return 1;

  QIO_destroy_record_info(rec_info);
  QIO_string_destroy(xml_record_out);

  return 0;
}

int write_F3_V_from_temp(QIO_Writer *outfile, char *xml_write_lattice, 
			 su3_vector *src, int count)
{
  QIO_String *xml_record_out;
  int status;
  QIO_RecordInfo *rec_info;
  
  /* We assume output precision is single */
  char qdptype[] = "QDP_F3_ColorVector";
  char prec[] = "F";
  int datum_size = sizeof(fsu3_vector);
  int word_size = sizeof(float);

  /* Create the record info for the field */
  rec_info = QIO_create_record_info(QIO_FIELD, qdptype, prec, 3,
				    0, datum_size, count);
  /* Create the record XML for the field */
  xml_record_out = QIO_string_create();
  QIO_string_set(xml_record_out,xml_write_lattice);

  /* Write the record for the field */
  status = QIO_write(outfile, rec_info, xml_record_out, vget_F3_V_from_temp, 
		     count*datum_size, word_size, (void *)src);
  if(status != QIO_SUCCESS)return 1;

  QIO_destroy_record_info(rec_info);
  QIO_string_destroy(xml_record_out);

  return 0;
}

/* Factory function for moving single precision KS vector from input
   to site structure */
void vput_F3_V(char *buf, size_t index, int count, void *arg)
{
  int i,j;
  fsu3_vector *src = (fsu3_vector *)buf;
  field_offset dest = *((field_offset *)arg);
  site *s = &lattice[index];
  /* Destination can be any precision */
  su3_vector *dest_vec = (su3_vector *)F_PT(s,dest);
  
  /* Copy, changing precision, if necessary */
  for(j = 0; j < count; j++)
    for(i = 0; i < 3; i++){
      dest_vec[j].c[i].real = src[j].c[i].real;
      dest_vec[j].c[i].imag = src[j].c[i].imag;
    }
}

/* Factory function for moving single precision KS vector from input
   to field-major vector */
void vput_F3_V_to_temp(char *buf, size_t index, int count, void *arg)
{
  int i,j;
  fsu3_vector *src = (fsu3_vector *)buf;
  su3_vector *dest = (su3_vector *)arg;
  /* Destination can be any precision */
  su3_vector *dest_vec = dest + index * count;
  
  /* Copy, changing precision, if necessary */
  for(j = 0; j < count; j++)
    for(i = 0; i < 3; i++){
      dest_vec[j].c[i].real = src[j].c[i].real;
      dest_vec[j].c[i].imag = src[j].c[i].imag;
    }
}

/* Read a set of color vectors */
int read_F3_V(QIO_Reader *infile, field_offset dest, int count)
{
  QIO_String *xml_record_in;
  QIO_RecordInfo rec_info;
  int status;
  /* We assume input precision is single */
  int datum_size = sizeof(fsu3_vector);
  int word_size = sizeof(float);
  
  /* Read the field record */
  xml_record_in = QIO_string_create();
  status = QIO_read(infile, &rec_info, xml_record_in, 
		    vput_F3_V, count*datum_size, word_size, (void *)&dest);
  node0_printf("Record info \n\"%s\"\n",QIO_string_ptr(xml_record_in));
  if(status != QIO_SUCCESS)return 1;

  node0_printf("Checksums %x %x\n",
	       QIO_get_reader_last_checksuma(infile),
	       QIO_get_reader_last_checksumb(infile));

  QIO_string_destroy(xml_record_in);
  return 0;
}

/* Read a set of color vectors */
int read_F3_V_to_temp(QIO_Reader *infile, su3_vector *dest, int count)
{
  QIO_String *xml_record_in;
  QIO_RecordInfo rec_info;
  int status;
  int word_size = sizeof(float);

  /* We assume input precision is single */
  int datum_size = sizeof(fsu3_vector);
  
  /* Read the field record */
  xml_record_in = QIO_string_create();
  status = QIO_read(infile, &rec_info, xml_record_in, 
		    vput_F3_V_to_temp, count*datum_size, 
		    word_size, (void *)dest);
  node0_printf("Record info \n\"%s\"\n",QIO_string_ptr(xml_record_in));
  if(status != QIO_SUCCESS)return 1;

  node0_printf("Checksums %x %x\n",
	       QIO_get_reader_last_checksuma(infile),
	       QIO_get_reader_last_checksumb(infile));

  QIO_string_destroy(xml_record_in);
  return 0;
}

/* Write color vectors in SciDAC format, taking data from the site
   structure */
/* We don't have a MILC format for such a file */

QIO_Writer *save_ks_vector_scidac(char *filename, char *recxml, int volfmt, 
				  field_offset src, int count)
{
  QIO_Layout layout;
  QIO_Writer *outfile;
  int status;
  char *filexml;

  QIO_verbose(QIO_VERB_OFF);

  /* Build the layout structure */
  build_layout(&layout);

  /* Create the file metadata */
  filexml = create_ks_XML();

  /* Open file for writing */
  outfile = open_output(filename, volfmt, &layout, filexml);
  if(outfile == NULL)terminate(1);


  /* Write the lattice field */
  status = write_F3_V(outfile, recxml, src, count);
  if(status)terminate(1);
  
  /* Close the file */
  QIO_close_write(outfile);

  /* Write information */
  if(volfmt == QIO_SINGLEFILE)
    node0_printf("Saved KS vector serially to binary file %s\n",
		 filename);
  else if(volfmt == QIO_MULTIFILE)
    node0_printf("Saved KS vector as multifile to binary file %s\n",
	   filename);
  else if(volfmt == QIO_PARTFILE)
    node0_printf("Saved KS vector in partition format to binary file %s\n",
	   filename);

  node0_printf("Checksums %x %x\n",
	       QIO_get_writer_last_checksuma(outfile),
	       QIO_get_writer_last_checksumb(outfile));

  free_ks_XML(filexml);
  return outfile;
}

QIO_Writer *save_ks_vector_scidac_from_temp(char *filename, char *recxml, 
    int volfmt, su3_vector *src, int count)
{
  QIO_Layout layout;
  QIO_Writer *outfile;
  int status;
  char *filexml;

  QIO_verbose(QIO_VERB_OFF);

  /* Build the layout structure */
  build_layout(&layout);

  /* Create the file metadata */
  filexml = create_ks_XML();

  /* Open file for writing */
  outfile = open_output(filename, volfmt, &layout, filexml);
  if(outfile == NULL)terminate(1);


  /* Write the lattice field */
  status = write_F3_V_from_temp(outfile, recxml, src, count);
  if(status)terminate(1);
  
  /* Close the file */
  QIO_close_write(outfile);

  /* Write information */
  if(volfmt == QIO_SINGLEFILE)
    node0_printf("Saved KS vector serially to binary file %s\n",
		 filename);
  else if(volfmt == QIO_MULTIFILE)
    node0_printf("Saved KS vector as multifile to binary file %s\n",
	   filename);
  else if(volfmt == QIO_PARTFILE)
    node0_printf("Saved KS vector in partition format to binary file %s\n",
	   filename);

  node0_printf("Checksums %0x %0x\n",
	       QIO_get_writer_last_checksuma(outfile),
	       QIO_get_writer_last_checksumb(outfile));

  free_ks_XML(filexml);
  return outfile;
}

/* Read color vectors in SciDAC format */

QIO_Reader *restore_ks_vector_scidac(char *filename, field_offset dest,
				     int count){
  QIO_Layout layout;
  QIO_Reader *infile;
  int status;

  QIO_verbose(QIO_VERB_OFF);

  /* Build the layout structure */
  build_layout(&layout);

  /* Open file for reading */
  infile = open_input(filename, &layout);
  if(infile == NULL)terminate(1);

  /* Read the lattice field: one color vector */
  status = read_F3_V(infile, dest, count);
  if(status)terminate(1);

  /* Close the file */
  QIO_close_read(infile);

  return infile;
}

/* Read color vectors in SciDAC format */

QIO_Reader *restore_ks_vector_scidac_to_temp(char *filename, 
			su3_vector *dest, int count){
  QIO_Layout layout;
  QIO_Reader *infile;
  int status;

  QIO_verbose(QIO_VERB_OFF);

  /* Build the layout structure */
  build_layout(&layout);

  /* Open file for reading */
  infile = open_input(filename, &layout);
  if(infile == NULL)terminate(1);

  /* Read the lattice field: one color vector */
  status = read_F3_V_to_temp(infile, dest, count);
  if(status)terminate(1);

  /* Close the file */
  QIO_close_read(infile);

  return infile;
}

