/********************** io_scidac_ks.c **********************************/
/* MILC/QIO interface for KS fields */
/* MIMD version 7 */
/* CD 2/2005
   CD 6/2007 Converted standard routines to macros and moved to
   generic/io_scidac.c
   CD 12/2007 Added proposed USQCD propagator formats
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

/********************************************************************/
/* Generic color vector file (not USQCD)                            */
/* Write color vectors in SciDAC format, taking data from the site
   structure */

QIO_Writer *w_open_ks_vector_scidac_file(char *filename, char *fileinfo, 
					 int volfmt, int serpar)
{
  QIO_Layout layout;
  QIO_Filesystem fs;
  QIO_Writer *outfile;
  QIO_String *filexml;

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

int save_ks_vector_scidac(QIO_Writer *outfile, char *filename, char *recinfo,
			  int volfmt, su3_vector *src, int count)
{
  QIO_String *recxml;
  int status;

  recxml = QIO_string_create();
  QIO_string_set(recxml, recinfo);
  if(PRECISION == 1)
    status = write_F3_V_from_field(outfile, recxml, src, count);
  else
    status = write_D3_V_from_field(outfile, recxml, src, count); 

  status = write_F3_V_from_field(outfile, recxml, src, count);

  QIO_string_destroy(recxml);
  if(status != QIO_SUCCESS)return status;
  
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

  return status;
}

void w_close_ks_vector_scidac_file(QIO_Writer *outfile)
{
  QIO_close_write(outfile);
}

/********************************************************************/
/* Generic color vector file (not USQCD)                            */
/* Write color vectors in SciDAC format, taking data from the site
   structure */

void save_ks_vector_scidac_from_field(char *filename, char *fileinfo,
				      char *recinfo, 
				      int volfmt, int serpar, 
				      su3_vector *src, int count)
{
  QIO_Writer *outfile;
  int status;

  QIO_verbose(QIO_VERB_OFF);

  outfile = w_open_ks_vector_scidac_file(filename, fileinfo, volfmt, serpar);
  if(outfile == NULL)terminate(1);

  status = save_ks_vector_scidac(outfile, filename, recinfo, 
				 volfmt, src, count);
  if(status)terminate(1);
  
  w_close_ks_vector_scidac_file(outfile);
}

/********************************************************************/
/* Generic color vector file (not USQCD)                            */
/* Write color vectors in SciDAC format, taking data from the site
   structure */

void save_ks_vector_scidac_from_site(char *filename, char *fileinfo,
				     char *recinfo, 
				     int volfmt, int serpar, 
				     field_offset src, int count)
{

  su3_vector *tmp;
  int i,j; site *s;

  tmp = (su3_vector *)malloc(sites_on_node * count * sizeof(su3_vector));
  if(tmp == NULL){
    printf("save_ks_vector_scidac_from_site(%d): No room for tmp\n",
	   this_node);
    terminate(1);
  }

  FORALLSITES(i,s){
    for(j = 0; j < count; j++)
      su3vec_copy((su3_vector *)F_PT(s,src+j*sizeof(su3_vector)),tmp+count*i+j);
  }

  save_ks_vector_scidac_from_field(filename, fileinfo, recinfo, 
				   volfmt, serpar, tmp,count);
  free(tmp);
}

/********************************************************************/
/* Read color vectors in SciDAC format (non-USQCD)                  */

QIO_Reader *r_open_ks_vector_scidac_file_xml(char *filename, int serpar,
					     QIO_String *xml_file)
{
  QIO_Layout layout;
  QIO_Filesystem fs;
  QIO_Reader *infile;

  /* Build the layout structure */
  build_qio_layout(&layout);

  /* Define the I/O nodes */
  build_qio_filesystem(&fs);

  /* Open file for reading */
  infile = open_scidac_input_xml(filename, &layout, &fs, serpar, xml_file);

  if(infile == NULL)return NULL;

  if(this_node==0){
    printf("Restoring binary SciDAC file %s\n",filename);
    printf("File info \n\"%s\"\n",QIO_string_ptr(xml_file));
  }

  return infile;
}

QIO_Reader *r_open_ks_vector_scidac_file(char *filename, int serpar)
{
  QIO_Reader *infile;
  QIO_String *xml_file;

  /* Open file for reading */
  xml_file = QIO_string_create();
  infile = r_open_ks_vector_scidac_file_xml(filename, serpar, xml_file);
  QIO_string_destroy(xml_file);
  return infile;
}

int read_ks_vector_scidac_xml(QIO_Reader *infile, su3_vector *dest, int count,
			      QIO_String *recxml)
{
  int status, typesize;
  QIO_RecordInfo recinfo;
  char myname[] = "read_ks_vector_scidac_xml";

  /* Check the record type (double or single precision) */
  status = QIO_read_record_info(infile, &recinfo, recxml);
  if(status != QIO_SUCCESS)return status;
  typesize = QIO_get_typesize(&recinfo);

  /* Read the lattice field: "count" color vectors per site */
  if(typesize == 6*4)
    /* Read them as a single precision record */
    return read_F3_V_to_field(infile, recxml, dest, count);
  else if(typesize == 6*8)
    /* Read them as a double precision record */
    return read_D3_V_to_field(infile, recxml, dest, count);
  else {
    node0_printf("%s: Incorrect data type size %d\n", myname, typesize);
    return 1;
  }

}

int read_ks_vector_scidac(QIO_Reader *infile, su3_vector *dest, int count)
{
  QIO_String *recxml;
  int status;

  /* Read the lattice field: "count" color vectors per site */
  recxml = QIO_string_create();
  status = read_ks_vector_scidac_xml(infile, dest, count, recxml);

  /* Discard for now */
  QIO_string_destroy(recxml);

  return status;
}

void r_close_ks_vector_scidac_file(QIO_Reader *infile)
{
  QIO_close_read(infile);
}

/********************************************************************/
/* Read color vectors in SciDAC format (non-USQCD)                  */
/* reads "count" vectors per site                                   */

void restore_ks_vector_scidac_to_field(char *filename, int serpar, 
				       su3_vector *dest, int count){
  QIO_Reader *infile;
  int status;

  QIO_verbose(QIO_VERB_OFF);

  infile = r_open_ks_vector_scidac_file(filename, serpar);
  if(infile == NULL)terminate(1);

  /* Read the lattice field: "count" color vectors per site */
  status = read_ks_vector_scidac(infile, dest, count);
  if(status != QIO_SUCCESS)terminate(1);

  r_close_ks_vector_scidac_file(infile);
}

/********************************************************************/
/* Read color vectors in SciDAC format (not USQCD) */

void restore_ks_vector_scidac_to_site(char *filename, int serpar,
				      field_offset dest, int count){
  su3_vector *tmp;
  int i,j;
  site *s;

  tmp = (su3_vector *)malloc(sites_on_node * count * sizeof(su3_vector));
  if(tmp == NULL){
    printf("restore_ks_vector_scidac_to_site(%d): No room for tmp\n",
	   this_node);
    terminate(1);
  }
			     
  restore_ks_vector_scidac_to_field(filename, serpar, tmp, count);

  FORALLSITES(i,s){
    for(j = 0; j < count; j++)
      su3vec_copy(tmp+count*i+j,(su3_vector *)F_PT(s,dest+j*sizeof(su3_vector)));
  }

  free(tmp);
}

/********************************************************************/
/* Routines for USQCD staggered propagator formats         */
/********************************************************************/

/* Write the file header for the propagator */

QIO_Writer *open_usqcd_ksprop_write(char *filename, int volfmt, 
				    int serpar, int ildgstyle, 
				    char *stringLFN, int milc_type,
				    char *fileinfo){
  
  QIO_Layout layout;
  QIO_Filesystem fs;
  QIO_Writer *outfile;
  QIO_USQCDKSPropFileInfo *propfile_info;
  QIO_String *filexml;
  int usqcd_type = ks_prop_milc_to_usqcd(milc_type);
  
  build_qio_layout(&layout);
  build_qio_filesystem(&fs);

  propfile_info = QIO_create_usqcd_kspropfile_info(usqcd_type, fileinfo);
  filexml = QIO_string_create();
  QIO_encode_usqcd_kspropfile_info(filexml, propfile_info);
  QIO_destroy_usqcd_kspropfile_info(propfile_info);

  outfile = open_scidac_output(filename, volfmt, serpar, ildgstyle, 
			       stringLFN, &layout, &fs, filexml);
  QIO_string_destroy(filexml);
  return outfile;
}

/********************************************************************/
/* Read the source record and extract the info string */

int read_kspropsource_C_usqcd(QIO_Reader *infile, char *srcinfo, int n,
			      complex *dest)
{
  QIO_USQCDKSPropSourceInfo propsource_info;
  QIO_String *recxml;
  char *info;
  int status, typesize;
  QIO_RecordInfo recinfo;

  recxml = QIO_string_create();

  /* Check the record type (double or single precision) */
  status = QIO_read_record_info(infile, &recinfo, recxml);
  if(status != QIO_SUCCESS)return status;
  typesize = QIO_get_typesize(&recinfo);

  if(typesize == 8)
    status = read_F_C_to_field(infile, recxml, dest, 1);
  else
    status = read_D_C_to_field(infile, recxml, dest, 1);

  if(status != QIO_SUCCESS)return status;

  status = QIO_decode_usqcd_kspropsource_info(&propsource_info, recxml);
  if(status != QIO_SUCCESS)terminate(1);
  QIO_string_destroy(recxml);

  info = QIO_get_usqcd_kspropsource_info(&propsource_info);

  strncpy(srcinfo, info, n);
  return QIO_SUCCESS;
}

/********************************************************************/
/* Read the next record as a color vector from an open file and
   extract the info string */

int read_kspropsource_V_usqcd(QIO_Reader *infile, char *srcinfo, int n,
			      su3_vector *dest)
{
  QIO_USQCDKSPropSourceInfo propsource_info;
  QIO_String *recxml;
  char *info;
  int status, typesize;
  QIO_RecordInfo recinfo;

  recxml = QIO_string_create();

  /* Check the record type (double or single precision) */
  status = QIO_read_record_info(infile, &recinfo, recxml);
  if(status != QIO_SUCCESS)return status;
  typesize = QIO_get_typesize(&recinfo);
  
  if(typesize == 6*4)
    /* Read as a single precision record */
    status = read_F3_V_to_field(infile, recxml, dest, 1);
  else
    /* Read as a double precision record */
    status = read_D3_V_to_field(infile, recxml, dest, 1);

  if(status != QIO_SUCCESS)return status;

  status = QIO_decode_usqcd_kspropsource_info(&propsource_info, recxml);
  if(status != QIO_SUCCESS)terminate(1);
  QIO_string_destroy(recxml);

  info = QIO_get_usqcd_kspropsource_info(&propsource_info);

  strncpy(srcinfo, info, n);
  return QIO_SUCCESS;
}

/********************************************************************/
/* Write a complex source field */

int write_kspropsource_C_usqcd_xml(QIO_Writer *outfile, QIO_String *recxml, 
			       complex *src, int t0){
  int status;
  if(PRECISION == 1)
    status = write_F_C_timeslice_from_field(outfile, recxml, src, 1, t0);
  else
    status = write_D_C_timeslice_from_field(outfile, recxml, src, 1, t0);
  return status;
}

/********************************************************************/
/* Encode the record info and write a complex source field */

int write_kspropsource_C_usqcd(QIO_Writer *outfile, char *srcinfo, 
			       complex *src, int t0){
  QIO_USQCDKSPropSourceInfo *propsource_info;
  QIO_String *recxml;
  int status;

  propsource_info = QIO_create_usqcd_kspropsource_info(srcinfo);
  recxml = QIO_string_create();
  QIO_encode_usqcd_kspropsource_info(recxml, propsource_info);
  QIO_destroy_usqcd_kspropsource_info(propsource_info);

  status = write_kspropsource_C_usqcd_xml(outfile, recxml, src, t0);
  QIO_string_destroy(recxml);
  return status;
}

/********************************************************************/
/* Write a KS vector source field */

int write_kspropsource_V_usqcd_xml(QIO_Writer *outfile, QIO_String *recxml,
				   su3_vector *src, int t0){
  int status;
  if(t0 == ALL_T_SLICES){
    if(PRECISION == 1)
      status = write_F3_V_from_field(outfile, recxml, src, 1);
    else
      status = write_D3_V_from_field(outfile, recxml, src, 1);
  }  else {
    if(PRECISION == 1)
      status = write_F3_V_timeslice_from_field(outfile, recxml, src, 1, t0);
    else
      status = write_D3_V_timeslice_from_field(outfile, recxml, src, 1, t0);
  }
  return status;
}

/********************************************************************/
/* Write a KS vector source field */

int write_kspropsource_V_usqcd(QIO_Writer *outfile, char *srcinfo, 
			       su3_vector *src, int t0){
  QIO_USQCDKSPropSourceInfo *propsource_info;
  QIO_String *recxml;
  int status;

  propsource_info = QIO_create_usqcd_kspropsource_info(srcinfo);
  recxml = QIO_string_create();
  QIO_encode_usqcd_kspropsource_info(recxml, propsource_info);
  QIO_destroy_usqcd_kspropsource_info(propsource_info);

  status = write_kspropsource_V_usqcd_xml(outfile, recxml, src, t0);
  QIO_string_destroy(recxml);
  return status;
}

/********************************************************************/
/* Write a KS vector solution field for a given source color */

int write_ksprop_usqcd_c(QIO_Writer *outfile, su3_vector *src, 
			 int color, char *recinfo)
{
  QIO_USQCDKSPropRecordInfo *proprecord_info;
  QIO_String *recxml;
  int status;

  proprecord_info = QIO_create_usqcd_ksproprecord_c_info(color, recinfo);
  recxml = QIO_string_create();
  QIO_encode_usqcd_ksproprecord_info(recxml, proprecord_info);
  QIO_destroy_usqcd_ksproprecord_info(proprecord_info);

  if(PRECISION == 1)
    status = write_F3_V_from_field(outfile, recxml, src, 1);
  else
    status = write_D3_V_from_field(outfile, recxml, src, 1);

  QIO_string_destroy(recxml);
  return status;
}

/********************************************************************/
/* Close the output KS propagator file */

void close_usqcd_ksprop_write(QIO_Writer *outfile){
  close_scidac_output(outfile);
}

/********************************************************************/
/* Read the file header for the propagator */

QIO_Reader *open_usqcd_ksprop_read(char *filename, int serpar, char **fileinfo){

  QIO_Layout layout;
  QIO_Filesystem fs;
  QIO_Reader *infile;
  QIO_String *xml_file;
  char *s;
  int len;

  xml_file = QIO_string_create();
  build_qio_layout(&layout);
  build_qio_filesystem(&fs);

  infile = open_scidac_input_xml(filename, &layout, &fs, serpar, xml_file);

  if(infile == NULL)terminate(1);

  /* Make a copy of the file metadata */
  s = QIO_string_ptr(xml_file);
  len = strlen(s);
  *fileinfo = (char *)malloc((len+1)*sizeof(char));
  if(*fileinfo == NULL){
    printf("open_usqcd_ksprop_read: No room for fileinfo\n");
    terminate(1);
  }
  strncpy(*fileinfo, s, len);
  (*fileinfo)[len] = '\0';
  QIO_string_destroy(xml_file);

  return infile;
}

/********************************************************************/
/* Read the next color vector record from an open file and return its
   color */

int read_ksproprecord_usqcd(QIO_Reader *infile, int *color, su3_vector *dest)
{
  QIO_USQCDKSPropRecordInfo proprecord_info;
  QIO_String *recxml;
  int status, typesize;
  QIO_RecordInfo recinfo;

  recxml = QIO_string_create();

  /* Check the record type (double or single precision) */
  status = QIO_read_record_info(infile, &recinfo, recxml);
  if(status != QIO_SUCCESS)return status;
  typesize = QIO_get_typesize(&recinfo);

  if(typesize == 6*4)
    /* Read as a single precision record */
    status = read_F3_V_to_field(infile, recxml, dest, 1);
  else
    /* Read as a double precision record */
    status = read_D3_V_to_field(infile, recxml, dest, 1);

  if(status != QIO_SUCCESS)return status;

  status = QIO_decode_usqcd_ksproprecord_info(&proprecord_info, recxml);
  if(status != QIO_SUCCESS)terminate(1);
  QIO_string_destroy(recxml);

  *color = QIO_get_usqcd_ksproprecord_color(&proprecord_info);
  return QIO_SUCCESS;
}

/********************************************************************/
/* Close the input file */

void close_usqcd_ksprop_read(QIO_Reader *infile){
  close_scidac_input(infile);
}
