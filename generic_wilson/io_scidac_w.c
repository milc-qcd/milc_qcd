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

void vget_F3_D(char *buf, size_t index, int count, void *arg)
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

int write_F3_D(QIO_Writer *outfile, char *xml_write_lattice, 
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
  status = QIO_write(outfile, rec_info, xml_record_out, vget_F3_D, 
		     count*datum_size, word_size, (void *)&src);
  if(status != QIO_SUCCESS)return 1;

  QIO_destroy_record_info(rec_info);
  QIO_string_destroy(xml_record_out);

  return 0;
}

/* Factory function for moving single precision Wilson vector from input
   to site structure */
void vput_F3_D(char *buf, size_t index, int count, void *arg)
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

/* Read a set of Wilson vectors */
int read_F3_D(QIO_Reader *infile, field_offset dest, int count)
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
		    vput_F3_D, count*datum_size, word_size, (void *)&dest);
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

QIO_Writer *save_w_vector_scidac(char *filename, char *recxml, int volfmt, 
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
  filexml = create_w_XML();

  /* Open file for writing */
  outfile = open_output(filename, volfmt, &layout, filexml);
  if(outfile == NULL)terminate(1);


  /* Write the lattice field */
  status = write_F3_D(outfile, recxml, src, count);
  if(status)terminate(1);
  
  /* Close the file */
  QIO_close_write(outfile);

  /* Write information */
  if(volfmt == QIO_SINGLEFILE)
    node0_printf("Saved Wilson vector serially to binary file %s\n",
		 filename);
  else if(volfmt == QIO_MULTIFILE)
    node0_printf("Saved Wilson vector as multifile to binary file %s\n",
	   filename);
  else if(volfmt == QIO_PARTFILE)
    node0_printf("Saved Wilson vector in partition format to binary file %s\n",
	   filename);

  node0_printf("Checksums %x %x\n",
	       QIO_get_writer_last_checksuma(outfile),
	       QIO_get_writer_last_checksumb(outfile));

  free_w_XML(filexml);
  return outfile;
}

/* Read Wilson vectors in SciDAC format */

QIO_Reader *restore_w_vector_scidac(char *filename, field_offset dest,
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

  /* Read the lattice field: one Wilson vector */
  status = read_F3_D(infile, dest, count);
  if(status)terminate(1);

  /* Close the file */
  QIO_close_read(infile);

  return infile;
}

