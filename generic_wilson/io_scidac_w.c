/********************** io_scidac_w.c **********************************/
/* MILC/QIO interface for Wilson fields */
/* MIMD version 7 */
/* CD 2/2005
   CD 6/2007 Converted standard routines to macros and moved to
   generic/io_scidac.c
*/

#include "generic_wilson_includes.h"
#ifndef HAVE_QIO
REQUIRES QIO
#else
#include <qio.h>
#endif
#include "../include/io_scidac.h"
#include "../include/io_scidac_w.h"
#include <string.h>

/********************************************************************/
/* Generic Wilson vector file (not USQCD)                           */
/* Write Wilson vectors in SciDAC format, taking data from the site
   structure */
/********************************************************************/
/* Generic Dirac vector file (not USQCD)                            */
/* Write Dirac vectors in SciDAC format, taking data from the site
   structure */

QIO_Writer *w_open_w_vector_scidac_file(char *filename, char *fileinfo, 
					int volfmt, int serpar)
{
  QIO_Layout layout;
  QIO_Filesystem fs;
  QIO_Writer *outfile;
  QIO_String *filexml;

  QIO_verbose(QIO_VERB_OFF);

  /* Build the layout structure */
  build_qio_layout(&layout);

  /* Define the I/O nodes */
  build_qio_filesystem(&fs);

  /* Open file for writing */
  filexml = QIO_string_create();
  QIO_string_set(filexml, fileinfo);
  outfile = open_scidac_output(filename, volfmt, serpar, QIO_ILDGNO,
			       NULL, &layout, &fs, filexml);
  QIO_string_destroy(filexml);
  return outfile;
}

int save_w_vector_scidac(QIO_Writer *outfile, char *filename, char *recinfo,
			 int volfmt, wilson_vector *src, int count)
{
  QIO_String *recxml;
  int status;

  recxml = QIO_string_create();
  QIO_string_set(recxml, recinfo);

  if(PRECISION == 1)
    status = write_F3_D_from_field(outfile, recxml, src, count);
  else
    status = write_D3_D_from_field(outfile, recxml, src, count);

  QIO_string_destroy(recxml);
  if(status)return status;
  
  /* Write information */
  if(volfmt == QIO_SINGLEFILE){
    node0_printf("Saved Dirac vector serially to binary file %s\n",
		 filename);
  }
  else if(volfmt == QIO_MULTIFILE){
    node0_printf("Saved Dirac vector as multifile to binary file %s\n",
	   filename);
  }
  else if(volfmt == QIO_PARTFILE){
    node0_printf("Saved Dirac vector in partition format to binary file %s\n",
	   filename);
  }

  node0_printf("Checksums %x %x\n",
	       QIO_get_writer_last_checksuma(outfile),
	       QIO_get_writer_last_checksumb(outfile));

  return status;
}

void w_close_w_vector_scidac_file(QIO_Writer *outfile)
{
  QIO_close_write(outfile);
}


/********************************************************************/
/* Write Wilson vectors in SciDAC format, taking data from a field */
/* We don't have a MILC format for such a file */

void save_w_vector_scidac_from_field(char *filename, char *fileinfo,
				     char *recinfo, 
				     int volfmt, int serpar,
				     wilson_vector *src, int count)
{
  QIO_Writer *outfile;
  int status;

  QIO_verbose(QIO_VERB_OFF);

  outfile = w_open_w_vector_scidac_file(filename, fileinfo, volfmt, serpar);
  if(outfile == NULL)terminate(1);

  status = save_w_vector_scidac(outfile, filename, recinfo, 
				volfmt, src, count);

  if(status)terminate(1);

  w_close_w_vector_scidac_file(outfile);

}

/********************************************************************/
/* Generic Dirac vector file (not USQCD)                            */
/* Write Dirac vectors in SciDAC format, taking data from the site
   structure */

void save_w_vector_scidac_from_site(char *filename, char *fileinfo,
				    char *recinfo, int volfmt, int serpar,
				    field_offset src, int count)
{
  wilson_vector *tmp;
  int i,j; site *s;

  tmp = (wilson_vector *)malloc(sites_on_node * count * sizeof(wilson_vector));
  if(tmp == NULL){
    printf("save_w_vector_scidac_from_site(%d): No room for tmp\n",
	   this_node);
    terminate(1);
  }

  FORALLSITES(i,s){
    for(j = 0; j < count; j++)
      copy_wvec((wilson_vector *)F_PT(s,src+j*sizeof(wilson_vector)),
		tmp+count*i+j);
  }
  save_w_vector_scidac_from_field(filename, fileinfo, recinfo, 
				  volfmt, serpar, tmp,count);
  free(tmp);
}

/********************************************************************/
/* Read Wilson vectors in SciDAC format (non-USQCD)                 */

QIO_Reader *r_open_w_vector_scidac_file_xml(char *filename, int serpar,
					    QIO_String *xml_file)
{
  QIO_Layout layout;
  QIO_Filesystem fs;
  QIO_Reader *infile;
  QIO_String *xml_file_in;

  QIO_verbose(QIO_VERB_OFF);

  /* Build the layout structure */
  build_qio_layout(&layout);

  /* Define the I/O nodes */
  build_qio_filesystem(&fs);

  /* Allocate for the file XML string */
  xml_file_in = QIO_string_create();

  /* Open file for reading */
  infile = open_scidac_input_xml(filename, &layout, &fs, serpar, xml_file);

  if(infile == NULL)return NULL;

  if(this_node==0){
    printf("Restoring binary SciDAC file %s\n",filename);
    printf("File info \n\"%s\"\n",QIO_string_ptr(xml_file_in));
  }

  QIO_string_destroy(xml_file_in);

  return infile;
}

QIO_Reader *r_open_w_vector_scidac_file(char *filename, int serpar)
{
  QIO_Reader *infile;
  QIO_String *xml_file;

  /* Open file for reading */
  xml_file = QIO_string_create();
  infile = r_open_w_vector_scidac_file_xml(filename, serpar, xml_file);
  QIO_string_destroy(xml_file);
  return infile;
}

int read_w_vector_scidac_xml(QIO_Reader *infile, wilson_vector *dest, 
			     int count, QIO_String *recxml)
{
  int status, typesize;
  QIO_RecordInfo recinfo;
  char myname[] = "read_w_vector_scidac_xml";

  /* Check the record type (double or single precision) */
  status = QIO_read_record_info(infile, &recinfo, recxml);
  if(status != QIO_SUCCESS){
    printf("read_w_vector_scidac_xml: Can't read file info\n");
    terminate(1);
  }
  typesize = QIO_get_typesize(&recinfo);

  /* Read "count" vectors per site.  Each has 12 complex values  */
  if(typesize == 24*4)
    /* Read them as a single precision record */
    return read_F3_D_to_field(infile, recxml, dest, count);
  else if(typesize == 24*8)
    /* Read them as a double precision record */
    return read_D3_D_to_field(infile, recxml, dest, count);
  else{
    node0_printf("%s: Incorrect data type size %d\n", myname, typesize);
    return 1;
  }

}

int read_w_vector_scidac(QIO_Reader *infile, wilson_vector *dest, int count)
{
  QIO_String *recxml;
  int status;

  /* Read the lattice field: "count" vectors per site */
  recxml = QIO_string_create();
  status = read_w_vector_scidac_xml(infile, dest, count, recxml);

  /* Discard for now */
  QIO_string_destroy(recxml);

  return status;
}

void r_close_w_vector_scidac_file(QIO_Reader *infile)
{
  QIO_close_read(infile);
}

/********************************************************************/
/* Read Wilson vectors in SciDAC format (non-USQCD)                 */

void restore_w_vector_scidac_to_field(char *filename, int serpar,
				      wilson_vector *dest, int count){
  QIO_Reader *infile;
  int status;

  infile = r_open_w_vector_scidac_file(filename, serpar);
  if(infile == NULL)terminate(1);

  /* Read the lattice field: "count" vectors per site */
  status = read_w_vector_scidac(infile, dest, count);
  if(status)terminate(1);

  /* Close the file */
  r_close_w_vector_scidac_file(infile);
}

/********************************************************************/
/* Read Wilson vectors in SciDAC format (non USQCD) */

void restore_w_vector_scidac_to_site(char *filename, int serpar,
				     field_offset dest, int count){
  wilson_vector *tmp;
  int i,j;
  site *s;

  tmp = (wilson_vector *)malloc(sites_on_node * count * sizeof(wilson_vector));
  if(tmp == NULL){
    printf("restore_w_vector_scidac_to_site(%d): No room for tmp\n",
	   this_node);
    terminate(1);
  }
			     
  restore_w_vector_scidac_to_field(filename, serpar, tmp, count);

  FORALLSITES(i,s){
    for(j = 0; j < count; j++)
      copy_wvec(tmp+count*i+j,
		(wilson_vector *)F_PT(s,dest+j*sizeof(wilson_vector)));
  }

  free(tmp);
}

/********************************************************************/
/* Routines for USQCD propagator formats */
/********************************************************************/

/* Write the file header for the propagator */

QIO_Writer *w_open_usqcd_wprop_file(char *filename, int volfmt, 
				    int serpar, int ildgstyle, 
				    char *stringLFN, int milc_type,
				    char *fileinfo){
  
  QIO_Layout layout;
  QIO_Filesystem fs;
  QIO_Writer *outfile;
  QIO_USQCDPropFileInfo *propfile_info;
  QIO_String *filexml;
  int usqcd_type = w_prop_milc_to_usqcd(milc_type);
  
  QIO_verbose(QIO_VERB_OFF);

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

int write_wpropsource_C_usqcd_xml(QIO_Writer *outfile, QIO_String *recxml, 
				  complex *src, int t0){
  int status;

  if(PRECISION == 1)
    status = write_F_C_timeslice_from_field(outfile, recxml, src, 1, t0);
  else
    status = write_D_C_timeslice_from_field(outfile, recxml, src, 1, t0);
  return status;
}

/********************************************************************/
/* Encode the record XML and write a complex source field */

int write_wpropsource_C_usqcd(QIO_Writer *outfile, char *srcinfo, 
			      complex *src, int t0){
  QIO_USQCDPropSourceInfo *propsource_info;
  QIO_String *recxml;
  int status;

  propsource_info = QIO_create_usqcd_propsource_info(srcinfo);
  recxml = QIO_string_create();
  QIO_encode_usqcd_propsource_info(recxml, propsource_info);
  QIO_destroy_usqcd_propsource_info(propsource_info);

  status = write_wpropsource_C_usqcd_xml(outfile, recxml, src, t0);
  QIO_string_destroy(recxml);
  return status;
}

/********************************************************************/
/* Write a Wilson vector source field */

int write_wpropsource_D_usqcd_xml(QIO_Writer *outfile, QIO_String *recxml, 
				  wilson_vector *src, int t0){
  int status;

  if(t0 == ALL_T_SLICES){
    if(PRECISION == 1)
      status = write_F3_D_from_field(outfile, recxml, src, 1);
    else
      status = write_D3_D_from_field(outfile, recxml, src, 1);
  }  else {
    if(PRECISION == 1)
      status = write_F3_D_timeslice_from_field(outfile, recxml, src, 1, t0);
    else
      status = write_D3_D_timeslice_from_field(outfile, recxml, src, 1, t0);
  }
  return status;
}

/********************************************************************/
/* Encode the record XML and write a Wilson vector source field */

int write_wpropsource_D_usqcd(QIO_Writer *outfile, char *srcinfo, 
			      wilson_vector *src, int t0){
  QIO_USQCDPropSourceInfo *propsource_info;
  QIO_String *recxml;
  int status;

  propsource_info = QIO_create_usqcd_propsource_info(srcinfo);
  recxml = QIO_string_create();
  QIO_encode_usqcd_propsource_info(recxml, propsource_info);
  QIO_destroy_usqcd_propsource_info(propsource_info);

  status = write_wpropsource_D_usqcd_xml(outfile, recxml, src, t0);
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

  QIO_verbose(QIO_VERB_OFF);

  proprecord_info = QIO_create_usqcd_proprecord_sc_info(spin, color, recinfo);
  recxml = QIO_string_create();
  QIO_encode_usqcd_proprecord_info(recxml, proprecord_info);
  QIO_destroy_usqcd_proprecord_info(proprecord_info);

  if(PRECISION == 1)
    status = write_F3_D_from_field(outfile, recxml, src, 1);
  else
    status = write_D3_D_from_field(outfile, recxml, src, 1);

  QIO_string_destroy(recxml);
  return status;
}

/********************************************************************/
/* Close the output Wilson propagator file */

void w_close_usqcd_wprop_file(QIO_Writer *outfile){
  close_scidac_output(outfile);
}

/********************************************************************/
/* Read the file header for the propagator */

QIO_Reader *r_open_usqcd_wprop_file(char *filename, int serpar){

  QIO_Layout layout;
  QIO_Filesystem fs;
  QIO_Reader *infile;

  QIO_verbose(QIO_VERB_OFF);

  build_qio_layout(&layout);
  build_qio_filesystem(&fs);

  infile = open_scidac_input(filename, &layout, &fs, serpar);
  if(infile == NULL)terminate(1);

  return infile;
}

/********************************************************************/
/* Read the source record and extract the info string */

int read_wpropsource_C_usqcd(QIO_Reader *infile, char *srcinfo, int n,
			     complex *dest)
{
  QIO_USQCDPropSourceInfo propsource_info;
  QIO_String *recxml;
  char *info;
  int status, typesize;
  QIO_RecordInfo recinfo;

  recxml = QIO_string_create();

  /* Check the record type (double or single precision) */
  status = QIO_read_record_info(infile, &recinfo, recxml);
  if(status)terminate(1);
  typesize = QIO_get_typesize(&recinfo);

  if(typesize == 8)
    status = read_F_C_to_field(infile, recxml, dest, 1);
  else
    status = read_D_C_to_field(infile, recxml, dest, 1);

  if(status != QIO_SUCCESS)return status;

  status = QIO_decode_usqcd_propsource_info(&propsource_info, recxml);
  if(status != QIO_SUCCESS)terminate(1);
  QIO_string_destroy(recxml);

  info = QIO_get_usqcd_propsource_info(&propsource_info);

  strncpy(srcinfo, info, n);
  return QIO_SUCCESS;
}

/********************************************************************/
/* Read the source record and extract the info string */

int read_wpropsource_D_usqcd(QIO_Reader *infile, char *srcinfo, int n,
			     wilson_vector *dest)
{
  QIO_USQCDPropSourceInfo propsource_info;
  QIO_String *recxml;
  char *info;
  int status, typesize;
  QIO_RecordInfo recinfo;

  recxml = QIO_string_create();

  /* Check the record type (double or single precision) */
  status = QIO_read_record_info(infile, &recinfo, recxml);
  if(status)terminate(1);
  typesize = QIO_get_typesize(&recinfo);
  
  if(typesize == 24*4)
    /* Read as a single precision record */
    status = read_F3_D_to_field(infile, recxml, dest, 1);
  else
    /* Read as a double precision record */
    status = read_D3_D_to_field(infile, recxml, dest, 1);

  if(status != QIO_SUCCESS)return QIO_SUCCESS;

  status = QIO_decode_usqcd_propsource_info(&propsource_info, recxml);
  if(status != QIO_SUCCESS)terminate(1);
  QIO_string_destroy(recxml);

  info = QIO_get_usqcd_propsource_info(&propsource_info);

  strncpy(srcinfo, info, n);
  return QIO_SUCCESS;
}

/********************************************************************/
/* Read the solution Dirac field and return its spin and color */

int read_wproprecord_usqcd(QIO_Reader *infile, int *spin, int *color, 
			   wilson_vector *dest)
{
  QIO_USQCDPropRecordInfo proprecord_info;
  QIO_String *recxml;
  int status, typesize;
  QIO_RecordInfo recinfo;

  recxml = QIO_string_create();

  /* Check the record type (double or single precision) */
  status = QIO_read_record_info(infile, &recinfo, recxml);
  if(status)terminate(1);
  typesize = QIO_get_typesize(&recinfo);

  if(typesize == 24*4)
    /* Read as a single precision record */
    status = read_F3_D_to_field(infile, recxml, dest, 1);
  else
    /* Read as a double precision record */
    status = read_D3_D_to_field(infile, recxml, dest, 1);

  if(status != QIO_SUCCESS)return status;

  status = QIO_decode_usqcd_proprecord_info(&proprecord_info, recxml);
  if(status)terminate(1);
  QIO_string_destroy(recxml);

  *spin  = QIO_get_usqcd_proprecord_spin(&proprecord_info);
  *color = QIO_get_usqcd_proprecord_color(&proprecord_info);
  return QIO_SUCCESS;
}

/********************************************************************/
/* Close the input file */

void r_close_usqcd_wprop_file(QIO_Reader *infile){
  close_scidac_input(infile);
}

