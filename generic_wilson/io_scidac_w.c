/********************** io_scidac_w.c **********************************/
/* MILC/QIO interface for Wilson fields */
/* MIMD version 7 */
/* CD 2/2005
   CD 6/2007 Converted standard routines to macros and moved to
   generic/io_scidac.c
*/

#include "generic_wilson_includes.h"
#include <qio.h>
#include "../include/io_scidac.h"
#include "../include/io_scidac_w.h"
#include <string.h>

/********************************************************************/
/* Translate from MILC to USQCD file type */
int w_prop_usqcd_to_milc(int usqcd_type){

  int milc_type = -1;

  switch(usqcd_type){
  case QIO_USQCDPROPFILETYPE_C1D12:
    milc_type = FILE_TYPE_W_USQCD_C1D12;
    break;
  case QIO_USQCDPROPFILETYPE_DD_PAIRS:
    milc_type = FILE_TYPE_W_USQCD_DD_PAIRS;
    break;
  case QIO_USQCDPROPFILETYPE_CD_PAIRS:
    milc_type = FILE_TYPE_W_USQCD_CD_PAIRS;
    break;
  case QIO_USQCDPROPFILETYPE_LHPC:
    milc_type = FILE_TYPE_W_USQCD_LHPC;
    break;
  }
  return milc_type;
}

/********************************************************************/
/* Translate from MILC to USQCD file type */
int w_prop_milc_to_usqcd(int milc_type){

  int usqcd_type = -1;

  switch(milc_type){
  case FILE_TYPE_W_USQCD_C1D12:
    usqcd_type = QIO_USQCDPROPFILETYPE_C1D12;
    break;
  case FILE_TYPE_W_USQCD_DD_PAIRS:
    usqcd_type = QIO_USQCDPROPFILETYPE_DD_PAIRS;
    break;
  case FILE_TYPE_W_USQCD_CD_PAIRS:
    usqcd_type = QIO_USQCDPROPFILETYPE_CD_PAIRS;
    break;
  case FILE_TYPE_W_USQCD_LHPC:
    usqcd_type = QIO_USQCDPROPFILETYPE_LHPC;
    break;
  }
  return usqcd_type;
}

/********************************************************************/
/* Write Wilson vectors in SciDAC format, taking data from the site
   structure */
/* We don't have a MILC format for such a file */

void save_w_vector_scidac_from_site(char *filename, char *recinfo, 
				    int volfmt, int serpar,
				    field_offset src, int count)
{
  QIO_Layout layout;
  QIO_Filesystem fs;
  QIO_Writer *outfile;
  int status;
  char *fileinfo;
  QIO_String *filexml;
  QIO_String *recxml;

  QIO_verbose(QIO_VERB_OFF);

  /* Build the layout structure */
  build_qio_layout(&layout);

  /* Create the file metadata */
  fileinfo = create_w_QCDML();

  /* Define the I/O nodes */
  build_qio_filesystem(&fs);

  /* Open file for writing */
  filexml = QIO_string_create();
  QIO_string_set(filexml, fileinfo);
  outfile = open_scidac_output(filename, volfmt, serpar, QIO_ILDGNO,
			       NULL, &layout, &fs, filexml);
  if(outfile == NULL)terminate(1);
  QIO_string_destroy(filexml);

  /* Write the lattice field */
  recxml = QIO_string_create();
  QIO_string_set(recxml, recinfo);
  status = write_F3_D_from_site(outfile, recxml, src, count);
  if(status)terminate(1);
  QIO_string_destroy(recxml);
  
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

  free_w_QCDML(fileinfo);
}

/********************************************************************/
/* Write Wilson vectors in SciDAC format, taking data from a field */
/* We don't have a MILC format for such a file */

void save_w_vector_scidac_from_field(char *filename, char *recinfo, 
				     int volfmt, int serpar,
				     wilson_vector *src, int count)
{
  QIO_Layout layout;
  QIO_Filesystem fs;
  QIO_Writer *outfile;
  QIO_String *filexml;
  QIO_String *recxml;
  int status;
  char *fileinfo;

  QIO_verbose(QIO_VERB_OFF);

  /* Build the layout structure */
  build_qio_layout(&layout);

  /* Create the file metadata */
  fileinfo = create_w_QCDML();

  /* Define the I/O nodes */
  build_qio_filesystem(&fs);

  /* Open file for writing */
  filexml = QIO_string_create();
  QIO_string_set(filexml, fileinfo);
  outfile = open_scidac_output(filename, volfmt, serpar, QIO_ILDGNO,
			       NULL, &layout, &fs, filexml);
  if(outfile == NULL)terminate(1);
  QIO_string_destroy(filexml);


  /* Write the lattice field */
  recxml = QIO_string_create();
  QIO_string_set(recxml, recinfo);
  status = write_F3_D_from_field(outfile, recxml, src, count);
  if(status)terminate(1);
  QIO_string_destroy(recxml);
  
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

  free_w_QCDML(fileinfo);
}

/********************************************************************/
/* Read Wilson vectors in SciDAC format */

void restore_w_vector_scidac_to_site(char *filename, field_offset dest,
				     int serpar, int count){
  QIO_Layout layout;
  QIO_Filesystem fs;
  QIO_Reader *infile;
  QIO_String *recxml;
  int status;

  QIO_verbose(QIO_VERB_OFF);

  /* Build the layout structure */
  build_qio_layout(&layout);

  /* Define the I/O nodes */
  build_qio_filesystem(&fs);

  /* Open file for reading */
  infile = open_scidac_input(filename, &layout, &fs, serpar);
  if(infile == NULL)terminate(1);

  /* Read the lattice field: one Wilson vector */
  recxml = QIO_string_create();
  status = read_F3_D_to_site(infile, recxml, dest, count);
  if(status)terminate(1);
  
  /* Discard for now */
  QIO_string_destroy(recxml);

  /* Close the file */
  QIO_close_read(infile);
}

/********************************************************************/
/* Read Wilson vectors in SciDAC format */

void restore_w_vector_scidac_to_field(char *filename, int serpar,
		     wilson_vector *dest, int count){
  QIO_Layout layout;
  QIO_Filesystem fs;
  QIO_Reader *infile;
  QIO_String *recxml;
  int status;

  QIO_verbose(QIO_VERB_OFF);

  /* Build the layout structure */
  build_qio_layout(&layout);

  /* Define the I/O nodes */
  build_qio_filesystem(&fs);

  /* Open file for reading */
  infile = open_scidac_input(filename, &layout, &fs, serpar);
  if(infile == NULL)terminate(1);

  /* Read the lattice field: one Wilson vector */
  recxml = QIO_string_create();
  status = read_F3_D_to_field(infile, recxml, dest, count);
  if(status)terminate(1);

  /* Discard for now */
  QIO_string_destroy(recxml);

  /* Close the file */
  QIO_close_read(infile);
}

/********************************************************************/
/* Routines for SciDAC propagator formats */
/********************************************************************/

#define MAX_XML 512

/* Write the file header for the propagator */

QIO_Writer *open_usqcd_prop_write(char *filename, int volfmt, 
				  int serpar, int ildgstyle, 
				  char *stringLFN, int milc_type,
				  char *fileinfo){
  
  QIO_Layout layout;
  QIO_Filesystem fs;
  QIO_Writer *outfile;
  QIO_USQCDPropFileInfo *propfile_info;
  QIO_String *filexml;
  int usqcd_type = w_prop_milc_to_usqcd(milc_type);
  
  build_qio_layout(&layout);
  build_qio_filesystem(&fs);

  propfile_info = QIO_create_usqcd_propfile_info(usqcd_type, fileinfo);
  filexml = QIO_string_create();
  QIO_encode_usqcd_propfile_info(filexml, propfile_info);
  QIO_destroy_usqcd_propfile_info(propfile_info);

  outfile = open_scidac_output(filename, volfmt, serpar, ildgstyle, 
			       stringLFN, &layout, &fs, filexml);
  QIO_string_destroy(filexml);
  return outfile;
}

/********************************************************************/
/* Write a complex source field */

int write_propsource_C_usqcd(QIO_Writer *outfile, char *srcinfo, complex *src, 
			 int t0){
  QIO_USQCDPropSourceInfo *propsource_info;
  QIO_String *recxml;
  int status;

  propsource_info = QIO_create_usqcd_propsource_info(srcinfo);
  recxml = QIO_string_create();
  QIO_encode_usqcd_propsource_info(recxml, propsource_info);
  QIO_destroy_usqcd_propsource_info(propsource_info);

  status = write_F_C_timeslice_from_field(outfile, recxml, src, 1, t0);
  QIO_string_destroy(recxml);
  return status;
}

/********************************************************************/
/* Write a Wilson vector source field */

int write_propsource_D_usqcd(QIO_Writer *outfile, char *srcinfo, 
			 wilson_vector *src, int t0){
  QIO_USQCDPropSourceInfo *propsource_info;
  QIO_String *recxml;
  int status;

  propsource_info = QIO_create_usqcd_propsource_info(srcinfo);
  recxml = QIO_string_create();
  QIO_encode_usqcd_propsource_info(recxml, propsource_info);
  QIO_destroy_usqcd_propsource_info(propsource_info);

  status = write_F3_D_timeslice_from_field(outfile, recxml, src, 1, t0);
  QIO_string_destroy(recxml);
  return status;
}

/********************************************************************/
/* Write a Wilson vector solution field for a given source spin and color */

int write_prop_usqcd_sc(QIO_Writer *outfile, wilson_vector *src, int spin, 
			int color, char *recinfo)
{
  QIO_USQCDPropRecordInfo *proprecord_info;
  QIO_String *recxml;
  int status;

  proprecord_info = QIO_create_usqcd_proprecord_sc_info(spin, color, recinfo);
  recxml = QIO_string_create();
  QIO_encode_usqcd_proprecord_info(recxml, proprecord_info);
  QIO_destroy_usqcd_proprecord_info(proprecord_info);

  status = write_F3_D_from_field(outfile, recxml, src, 1);
  QIO_string_destroy(recxml);
  return status;
}

/********************************************************************/
/* Close the output Wilson propagator file */

void close_usqcd_prop_write(QIO_Writer *outfile){
  close_scidac_output(outfile);
}

/********************************************************************/
/* Open a Wilson propagator file and discover its format */

int io_detect_w_usqcd(char *filename){

  QIO_Layout layout;
  QIO_Filesystem fs;
  QIO_Reader *infile;
  QIO_String *xml_file_in;
  QIO_USQCDPropFileInfo prop_info;
  int status, usqcd_type, milc_type;

  /* Allocate for the file XML string */
  xml_file_in = QIO_string_create();

  /* Build the layout structure */
  build_qio_layout(&layout);

  /* Define the I/O nodes */
  build_qio_filesystem(&fs);

  /* Open file for reading */
  infile = open_scidac_input_xml(filename, &layout, &fs, QIO_SERIAL, 
				 xml_file_in);
  if(infile == NULL)terminate(1);

  /* Decode the file XML */

  status = QIO_decode_usqcd_propfile_info(&prop_info, xml_file_in);
  QIO_string_destroy(xml_file_in);
  if(status)terminate(1);

  milc_type = -1;
  if(QIO_defined_usqcd_propfile_type(&prop_info))
    {
      /* Translate the file type */
      usqcd_type = QIO_get_usqcd_propfile_type(&prop_info);
      milc_type = w_prop_usqcd_to_milc(usqcd_type);
    }

  close_scidac_input(infile);
  return milc_type;
}

/********************************************************************/
/* Read the file header for the propagator */

QIO_Reader *open_usqcd_prop_read(char *filename, int serpar){

  QIO_Layout layout;
  QIO_Filesystem fs;
  QIO_Reader *infile;
  QIO_String *filexml;

  filexml = QIO_string_create();
  build_qio_layout(&layout);
  build_qio_filesystem(&fs);

  infile = open_scidac_input(filename, &layout, &fs, serpar);
  if(infile == NULL)terminate(1);

  return infile;
}

/********************************************************************/
/* Read the source record and extract the info string */

void read_propsource_C_usqcd(QIO_Reader *infile, char *srcinfo, int n,
			     complex *dest)
{
  QIO_USQCDPropSourceInfo propsource_info;
  QIO_String *recxml;
  char *info;
  int status;

  recxml = QIO_string_create();
  status = read_F_C_to_field(infile, recxml, dest, 1);
  if(status)terminate(1);

  status = QIO_decode_usqcd_propsource_info(&propsource_info, recxml);
  if(status)terminate(1);
  QIO_string_destroy(recxml);

  info = QIO_get_usqcd_propsource_info(&propsource_info);

  strncpy(srcinfo, info, n);
}

/********************************************************************/
/* Read the source record and extract the info string */

void read_propsource_D_usqcd(QIO_Reader *infile, char *srcinfo, int n,
			     wilson_vector *dest)
{
  QIO_USQCDPropSourceInfo propsource_info;
  QIO_String *recxml;
  char *info;
  int status;

  recxml = QIO_string_create();
  status = read_F3_D_to_field(infile, recxml, dest, 1);
  if(status)terminate(1);

  status = QIO_decode_usqcd_propsource_info(&propsource_info, recxml);
  if(status)terminate(1);
  QIO_string_destroy(recxml);

  info = QIO_get_usqcd_propsource_info(&propsource_info);

  strncpy(srcinfo, info, n);
}

/********************************************************************/
/* Read the solution Dirac field and return its spin and color */

void read_proprecord_usqcd(QIO_Reader *infile, int *spin, int *color, 
			   wilson_vector *dest)
{
  QIO_USQCDPropRecordInfo proprecord_info;
  QIO_String *recxml;
  int status;

  recxml = QIO_string_create();
  status = read_F3_D_to_field(infile, recxml, dest, 1);
  if(status)terminate(1);

  status = QIO_decode_usqcd_proprecord_info(&proprecord_info, recxml);
  if(status)terminate(1);
  QIO_string_destroy(recxml);

  *spin  = QIO_get_usqcd_proprecord_spin(&proprecord_info);
  *color = QIO_get_usqcd_proprecord_color(&proprecord_info);
}

/********************************************************************/
/* Close the input file */

void close_usqcd_prop_read(QIO_Reader *infile){
  close_scidac_input(infile);
}
