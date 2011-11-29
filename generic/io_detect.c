/*********************** io_detect.c *************************/
/* MIMD version 7 */

/* Read the 32-bit magic number of a file and look it up in a table */
/* Returns the file type number from the table, which must be >= 0 */
/* Error return is -1 (not in table), -2 (open error) -3 (read error) */
/* Multidump files do not have a magic number, so can't be identified
   this way */
/* FNAL files have a common magic number, so it is necessary to
   read further to identify them. */

#include "generic_includes.h"
#include "../include/file_types.h"
#ifdef HAVE_QIO
#include "../include/io_scidac.h"
#include "../include/io_scidac_ks.h"
#include "../include/io_scidac_w.h"
#endif
#include <string.h>

int io_detect(char *filename, file_table ft[], int ntypes){
  FILE *fp;
  int i, status, words;
  int32type magic_no;
  int32type revmagic_no;
  char editfilename[513];

  /* Node 0 reads and checks */
  if(this_node == 0){
    fp = g_open(filename,"rb");
    if(fp == NULL){
      /* Special provision for partition or multifile format.  Try
	 adding the extension to the filename */
      strncpy(editfilename,filename,504);
      editfilename[504] = '\0';  /* Just in case of truncation */
      strcat(editfilename,".vol0000");
      fp = g_open(editfilename,"rb");
    }

    if(fp == NULL)status = -2;
    else
      {
	words = g_read(&magic_no, sizeof(int32type), 1, fp);
	g_close(fp);

	if(words != 1)status = -3;
	else
	  {
	    revmagic_no = magic_no;
	    byterevn(&revmagic_no, 1);

	    status = -1;
	    for(i = 0; i < ntypes; i++){
	      if(ft[i].magic_no == magic_no || 
		 ft[i].magic_no == revmagic_no)
		{
		  status = ft[i].type;
		  break;
		}
	    }
	  }
      }
  }

  /* Node 0 broadcasts the result */
  broadcast_bytes((char *)&status, sizeof(int));

  /* All nodes return the same value */
  return status;
}

#ifndef ONLY_GLUON_FILES
/********************************************************************/
/* Open a staggered propagator file and discover its format */

#ifdef HAVE_QIO
int io_detect_ks_usqcd(char *filename){

  QIO_Layout layout;
  QIO_Filesystem fs;
  QIO_Reader *infile;
  QIO_String *xml_file_in;
  QIO_USQCDKSPropFileInfo prop_info;
  int status, usqcd_type, milc_type;

  /* Allocate for the file XML string */
  xml_file_in = QIO_string_create();

  /* Build the layout structure */
  build_qio_layout(&layout);

  /* Define the I/O nodes */
  build_qio_filesystem(&fs);

  /* Open file for reading */
  layout.latdim = 0; /* Force discovery of dimensions */
  infile = open_scidac_input_xml(filename, &layout, &fs, QIO_SERIAL, 
				 xml_file_in);
  if(infile == NULL)terminate(1);

  /* Decode the file XML */

  status = QIO_decode_usqcd_kspropfile_info(&prop_info, xml_file_in);
  QIO_string_destroy(xml_file_in);
  if(status)return -1;

  milc_type = -1;
  if(QIO_defined_usqcd_kspropfile_type(&prop_info))
    {
      /* Translate the file type */
      usqcd_type = QIO_get_usqcd_kspropfile_type(&prop_info);
      milc_type = ks_prop_usqcd_to_milc(usqcd_type);
    }

  close_scidac_input(infile);
  return milc_type;
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
  layout.latdim = 0; /* Force discovery of dimensions */
  infile = open_scidac_input_xml(filename, &layout, &fs, QIO_SERIAL, 
				 xml_file_in);
  if(infile == NULL)terminate(1);

  /* Decode the file XML */

  status = QIO_decode_usqcd_propfile_info(&prop_info, xml_file_in);
  QIO_string_destroy(xml_file_in);
  if(status)return -1;

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

#endif

#endif
/********************************************************************/
/* For FNAL we base the detection on the number of elements per site */
int io_detect_fm(char *filename){
  FILE *fp;
  int status, words;
  int32type magic_no, revmagic_no, gmtime_stamp, size_of_element,
    elem_per_site;
  int byterevflag = 0;

  /* Node 0 reads and checks */
  if(this_node == 0){
    fp = g_open(filename,"rb");
    if(fp == NULL)status = -2;
    else
      {
	words = g_read(&magic_no, sizeof(int32type), 1, fp);

	if(words != 1)status = -3;
	else
	  {
	    revmagic_no = magic_no;
	    byterevn(&revmagic_no, 1);

	    status = -1;
	    if(revmagic_no == IO_UNI_MAGIC)byterevflag = 1;
	    if(magic_no == IO_UNI_MAGIC || revmagic_no == IO_UNI_MAGIC){
	      g_read(&gmtime_stamp, sizeof(gmtime_stamp), 1, fp);
	      g_read(&size_of_element, sizeof(int32type), 1, fp);
	      g_read(&elem_per_site, sizeof(int32type), 1, fp);

	      if(byterevflag){
		byterevn(&size_of_element, 1);
		byterevn(&elem_per_site, 1);
	      }

	      if(size_of_element != sizeof(float))status = -1;
	      else 
		if(elem_per_site == sizeof(fsu3_matrix)/sizeof(float))
		  status = FILE_TYPE_KS_FMPROP;
	      else 
		if(elem_per_site == sizeof(fwilson_propagator)/sizeof(float) ||
		   elem_per_site == sizeof(fwilson_vector)/sizeof(float))
		  status = FILE_TYPE_W_FMPROP;
	      else 
		if(elem_per_site == 4*sizeof(fsu3_matrix)/sizeof(float))
		  status = FILE_TYPE_GAUGE_FNAL;
	    }
	  }
      }
    g_close(fp);
  }
  /* Node 0 broadcasts the result */
  broadcast_bytes((char *)&status, sizeof(int));

  /* All nodes return the same value */
  return status;
}

/*---------------------------------------------------------------*/
/* Sniff out the file type */

#define N_BROAD_FILE_TYPES 6
static file_table broad_file_types[N_BROAD_FILE_TYPES] =
  { 
    {FILE_TYPE_LIME,          LIME_MAGIC_NO},
    {FILE_TYPE_FM,            IO_UNI_MAGIC},
    {FILE_TYPE_GAUGE_V5,      GAUGE_VERSION_NUMBER},
    {FILE_TYPE_GAUGE_ARCHIVE, GAUGE_VERSION_NUMBER_ARCHIVE},
    {FILE_TYPE_KS_PROP,       KSPROP_VERSION_NUMBER_V0},
    {FILE_TYPE_KS_PROP,       KSPROP_VERSION_NUMBER},
  };

int 
get_file_type(char *filename)
{
  int file_type = FILE_TYPE_UNKNOWN;
  char myname[] = "get_file_type";
  
  file_type = io_detect(filename, broad_file_types, N_BROAD_FILE_TYPES);
  if(file_type < 0){
    node0_printf("%s: Error opening file %s\n",myname,filename);
    terminate(1);
  }
  
  /* Look further for FNAL files */
  if(file_type == FILE_TYPE_FM){
    file_type = io_detect_fm(filename);
    if(file_type < 0){
      node0_printf("%s: Don't recognize FNAL subtype %s\n",myname,filename);
      terminate(1);
    }
  }
  
  /* For QIO(LIME) types, same thing */
  else if(file_type == FILE_TYPE_LIME){
#ifdef HAVE_QIO
#ifndef ONLY_GLUON_FILES
    file_type = io_detect_w_usqcd(filename);
    /* If this fails, try the KS propagator types */
    if(file_type < 0)
      file_type = io_detect_ks_usqcd(filename);
#endif
    if(file_type < 0){
      node0_printf("%s: Don't recognize QIO file type for %s\n",
		   myname,filename);
      terminate(1);
    }
#else
    node0_printf("This looks like a QIO file, but to read it requires QIO compilation\n");
#endif
  }
  return file_type;
} /* get_file_type */
