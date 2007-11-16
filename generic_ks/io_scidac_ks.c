/********************** io_scidac_ks.c **********************************/
/* MILC/QIO interface for KS fields */
/* MIMD version 7 */
/* CD 2/2005
   CD 6/2007 Converted standard routines to macros and moved to
   generic/io_scidac.c
*/

#include "generic_ks_includes.h"
#ifndef HAVE_QIO
REQUIRES QIO
#else
#include <qio.h>
#endif
#include "../include/io_scidac.h"
#include "../include/io_scidac_ks.h"
#include <string.h>

/* Write color vectors in SciDAC format, taking data from the site
   structure */
/* We don't have a MILC format for such a file */

void save_ks_vector_scidac_from_site(char *filename, char *recinfo, 
			  int volfmt, int serpar, field_offset src, int count)
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
  fileinfo = create_ks_XML();

  /* Define the I/O nodes */
  build_qio_filesystem(&fs);

  /* Open file for writing */
  filexml = QIO_string_create();
  QIO_string_set(filexml,fileinfo);
  outfile = open_scidac_output(filename, volfmt, serpar, QIO_ILDGNO,
			       NULL, &layout, &fs, filexml);
  if(outfile == NULL)terminate(1);
  QIO_string_destroy(filexml);


  /* Write the lattice field */
  recxml = QIO_string_create();
  QIO_string_set(recxml, recinfo);
  status = write_F3_V_from_site(outfile, recxml, src, count);
  if(status)terminate(1);
  QIO_string_destroy(recxml);
  
  /* Write information */
  if(volfmt == QIO_SINGLEFILE){
    node0_printf("Saved KS vector serially to binary file %s\n",
		 filename);
  }
  else if(volfmt == QIO_MULTIFILE){
    node0_printf("Saved KS vector as multifile to binary file %s\n",
	   filename);
  }
  else if(volfmt == QIO_PARTFILE){
    node0_printf("Saved KS vector in partition format to binary file %s\n",
	   filename);
  }

  node0_printf("Checksums %x %x\n",
	       QIO_get_writer_last_checksuma(outfile),
	       QIO_get_writer_last_checksumb(outfile));

  /* Close the file */
  QIO_close_write(outfile);

  free_ks_XML(fileinfo);
}

void save_ks_vector_scidac_from_field(char *filename, char *recinfo, 
    int volfmt, int serpar, su3_vector *src, int count)
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
  fileinfo = create_ks_XML();

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
  status = write_F3_V_from_field(outfile, recxml, src, count);
  if(status)terminate(1);
  QIO_string_destroy(recxml);
  
  /* Write information */
  if(volfmt == QIO_SINGLEFILE){
    node0_printf("Saved KS vector serially to binary file %s\n",
		 filename);
  }
  else if(volfmt == QIO_MULTIFILE){
    node0_printf("Saved KS vector as multifile to binary file %s\n",
	   filename);
  }
  else if(volfmt == QIO_PARTFILE){
    node0_printf("Saved KS vector in partition format to binary file %s\n",
	   filename);
  }

  node0_printf("Checksums %x %x\n",
	       QIO_get_writer_last_checksuma(outfile),
	       QIO_get_writer_last_checksumb(outfile));

  /* Close the file */
  QIO_close_write(outfile);

  free_ks_XML(fileinfo);
}

/* Read color vectors in SciDAC format */

void restore_ks_vector_scidac_to_site(char *filename, field_offset dest,
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

  /* Read the lattice field: one color vector */
  recxml = QIO_string_create();
  status = read_F3_V_to_site(infile, recxml, dest, count);
  if(status)terminate(1);

  /* Discard for now */
  QIO_string_destroy(recxml);

  /* Close the file */
  QIO_close_read(infile);

}

/* Read color vectors in SciDAC format */

void restore_ks_vector_scidac_to_field(char *filename, 
			      su3_vector *dest, int serpar, int count){
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

  /* Read the lattice field: one color vector */
  recxml = QIO_string_create();
  status = read_F3_V_to_field(infile, recxml, dest, count);
  if(status)terminate(1);

  /* Discard for now */
  QIO_string_destroy(recxml);

  /* Close the file */
  QIO_close_read(infile);

}

