/************************ quark_source_io.c *****************************/
/* MIMD version 7 */

/* File-types and operations supported

   Complex source files
     Reading in USQCD (SciDAC) or FNAL format.
     Not writing in any format so far.

   Color vector source files
     Reading in USQCD or FNAL format
     Writing in USQCD format

   Dirac vector source files
     Reading in USQCD or FNAL format
     Writing in USQCD format

   In USQCD format the files may have any number of source fields.  In
   FNAL format, Dirac field files must have the traditional 12 color
   and spin source fields.

*/

/* Coding notes 

   Wrappers are not needed, so not provided here for
   reading complex, color vector, or Dirac FNAL files

*/

/* External entries

   choose_usqcd_ks_file_type
   choose_usqcd_w_file_type
   r_source_open
   r_source_close
   r_source_cmplx_scidac
   r_source_cmplx_scidac_open
   r_source_vector
   r_source_dirac
   w_source_open_ks
   w_source_open_dirac
   w_source_close
   w_source_ks
   w_source_dirac
   w_source_dirac_site
   ask_output_quark_source_file
 */

#include "generic_includes.h"
#include "../include/io_ksprop.h"
#include "../include/io_wprop.h"
#include <string.h>
#ifdef HAVE_QIO
#include <qio.h>
#include "../include/io_scidac.h"
#include "../include/io_scidac_ks.h"
#include "../include/io_scidac_w.h"
#endif

/*--------------------------------------------------------------------*/
/* USQCD propagator file types according to source type */
/*--------------------------------------------------------------------*/

/* Staggered propagator files */
int choose_usqcd_ks_file_type(int source_type){
  int file_type;

  switch(source_type){
    /* complex field sources */
  case COMPLEX_FIELD_FILE:
  case COMPLEX_FIELD_FM_FILE:
  case COMPLEX_FIELD_STORE:
  case CORNER_WALL:
  case EVEN_WALL:
  case EVENANDODD_WALL:
  case EVENMINUSODD_WALL:
  case GAUSSIAN:
  case POINT:
  case UNKNOWN:
    file_type = FILE_TYPE_KS_USQCD_C1V3;
    break;
    /* vector field sources */
  case FAT_COVARIANT_GAUSSIAN:
  case RANDOM_CORNER_COLOR_WALL:
  case RANDOM_COLOR_WALL:
  case VECTOR_FIELD_FILE:
  case VECTOR_FIELD_FM_FILE:
  case VECTOR_FIELD_STORE:
  case WAVEFUNCTION_FILE:
    file_type = FILE_TYPE_KS_USQCD_VV_PAIRS;
    break;
  default:
    file_type = -1;
  }
  return file_type;
}

/* Dirac propagator files */
int choose_usqcd_w_file_type(int source_type){
  int file_type;

  switch(source_type){
  case COMPLEX_FIELD_FILE:
  case COMPLEX_FIELD_FM_FILE:
  case COMPLEX_FIELD_STORE:
  case GAUSSIAN:
  case POINT:
  case WAVEFUNCTION_FILE:
  case UNKNOWN:
    file_type = FILE_TYPE_W_USQCD_C1D12;
    break;
  case COVARIANT_GAUSSIAN:
  case DIRAC_FIELD_FILE:
  case DIRAC_FIELD_FM_FILE:
  case DIRAC_FIELD_STORE:
  case DIRAC_PROPAGATOR_FILE:
  case FAT_COVARIANT_GAUSSIAN:
  case DERIV1:
  case DERIV2_D:
  case DERIV2_B:
  case DERIV3_A:
  case RANDOM_COLOR_WALL:
  case ROTATE_3D:
  case VECTOR_FIELD_FILE:
  case VECTOR_FIELD_FM_FILE:
  case VECTOR_PROPAGATOR_FILE:
    file_type = FILE_TYPE_W_USQCD_DD_PAIRS;
    break;
  default:
    file_type = -1;
  }
  return file_type;
}

#ifdef HAVE_QIO

#ifdef HAVE_DIRAC
/********************************************************************/
/* Parse the record XML to get the color and check it               */
/********************************************************************/
/* For the source file we borrow the XML encoding from the USQCD
   Wilson propagator file  */
static int check_color_spin(QIO_String *recxml, int color, int spin){
  int status;
  int input_color, input_spin;
  QIO_USQCDPropRecordInfo recinfo;
  char myname[] = "check_color_spin";

  status = QIO_decode_usqcd_proprecord_info(&recinfo, recxml);
  if(status != QIO_SUCCESS) 
    return qio_status(status);
  input_color = QIO_get_usqcd_proprecord_color(&recinfo);
  input_spin = QIO_get_usqcd_proprecord_spin(&recinfo);
  if(color != input_color || spin  != input_spin ){
    node0_printf("%s(%d): Error: expected color %d and spin %d got %d and %d\n",
		 myname, this_node, color, spin, 
		 input_color, input_spin);
    return 1;
  }
  return 0;
}

#endif /* HAVE_DIRAC */

#ifdef HAVE_KS

/********************************************************************/
/* Parse the record XML to get the color and check it               */
/********************************************************************/
/* For a KS source file we use the same XML encoding as with the KS
   propagator file  */
static int check_color(QIO_String *recxml, int color){
  int status;
  int input_color;
  QIO_USQCDKSPropRecordInfo recinfo;
  char myname[] = "check_color";

  status = QIO_decode_usqcd_ksproprecord_info(&recinfo, recxml);
  if(status != QIO_SUCCESS) 
    return qio_status(status);
  input_color = QIO_get_usqcd_ksproprecord_color(&recinfo);
  if(color != input_color){
    node0_printf("%s(%d): Error: expected color %d got %d\n",
		 myname, this_node, color, input_color);
    return 1;
  }
  return 0;
}

#endif /* HAVE_KS */

#endif /* HAVE_QIO */

/*--------------------------------------------------------------------*/
/* Routines for reading sources */
/*--------------------------------------------------------------------*/

/*--------------------------------------------------------------------*/
/* Open a complex, color vector or Dirac vector source file for reading */
/*--------------------------------------------------------------------*/

#ifdef HAVE_QIO

QIO_Reader *r_source_cmplx_scidac_open(char source_file[]){
  QIO_String *xml_file;
  QIO_Reader *infile;

  xml_file = QIO_string_create();
  infile = r_open_complex_scidac_file_xml(source_file, QIO_SERIAL,
					  xml_file);
  QIO_string_destroy(xml_file);
  return infile;
}
#endif

void r_source_open(quark_source *qs){

  char *source_file         = qs->source_file;
  int source_type           = qs->type;
#ifdef HAVE_QIO
  QIO_String *xml_file;
#endif
  char myname[] = "r_source_open";

  if(qs->source_file_initialized != 0){
    node0_printf("%s: file already opened\n",myname);
    return;
  }


  if(0);

#ifdef HAVE_KS

  else if(source_type == VECTOR_FIELD_FM_FILE){
    qs->kssf = r_source_ks_fm_i(source_file);
    if(qs->kssf == NULL){
      node0_printf("%s: Failed to open source %s\n", myname, source_file);
    }
  }
#endif

#ifdef HAVE_QIO

  else if(source_type == COMPLEX_FIELD_FILE){
    qs->infile = r_source_cmplx_scidac_open(source_file);
  }

#ifdef HAVE_KS

  else if(source_type == VECTOR_FIELD_FILE){
    int serpar;
    if(qs->sourceflag == RELOAD_PARALLEL)
      serpar = QIO_PARALLEL;
    else
      serpar = QIO_SERIAL;

    xml_file = QIO_string_create();
    qs->infile = r_open_ks_vector_scidac_file_xml(source_file, serpar,
						   xml_file);
    QIO_string_destroy(xml_file);
  }

#endif

#ifdef HAVE_DIRAC

  else if(source_type == DIRAC_FIELD_FILE){
    int serpar;
    if(qs->sourceflag == RELOAD_PARALLEL)
      serpar = QIO_PARALLEL;
    else
      serpar = QIO_SERIAL;

    xml_file = QIO_string_create();
    qs->infile = r_open_w_vector_scidac_file_xml(source_file, serpar,
						 xml_file);
    QIO_string_destroy(xml_file);
  }

#endif

  else {
    node0_printf("%s: bad source type %d\n",myname, source_type);
    return;
  }
#else
  else {
    node0_printf("%s: QIO compilation required for this source\n", myname);
    terminate(1);
  }
#endif

#ifdef HAVE_QIO
  if(qs->infile == NULL && qs->kssf == NULL)terminate(1);
#else
  if(qs->kssf == NULL)terminate(1);
#endif
  qs->source_file_initialized = 1;
} /* r_source_open */

/*--------------------------------------------------------------------*/
/* Close a complex, color vector or Dirac vector source after reading */
/*--------------------------------------------------------------------*/

void r_source_close(quark_source *qs){

  char myname[] = "r_source_close";

  if(qs->source_file_initialized == 0){
    node0_printf("%s: file has not been opened\n",myname);
    return;
  }

  if(0);

#ifdef HAVE_KS
  else if(qs->type == VECTOR_FIELD_FM_FILE){
    r_source_ks_fm_f(qs->kssf);
    qs->kssf = NULL;
  }
#endif
#ifdef HAVE_QIO
  else if(qs->type == COMPLEX_FIELD_FILE)
    r_close_complex_scidac_file(qs->infile);
#ifdef HAVE_KS
  else if(qs->type == VECTOR_FIELD_FILE)
    r_close_ks_vector_scidac_file(qs->infile);
#endif
#ifdef HAVE_DIRAC
  else if(qs->type == DIRAC_FIELD_FILE)
    r_close_w_vector_scidac_file(qs->infile);
#endif
  else
    node0_printf("%s: bad source type %d\n",myname, qs->type);

  qs->infile = NULL;
#else
  else {
    node0_printf("%s: QIO compilation required for this source\n", myname);
    terminate(1);
  }
#endif
  qs->source_file_initialized = 0;
} /* r_source_close */

#ifdef HAVE_QIO
/*--------------------------------------------------------------------*/
/* Read a USQCD complex source file */
/*--------------------------------------------------------------------*/

int r_source_cmplx_scidac(QIO_Reader *infile, complex *src,  
			  int x0, int y0, int z0, int t0){

  int status;

  status = qio_status(read_complex_scidac(infile, src, 1 ));

  return status;

} /* r_source_cmplx_scidac */

#endif

#if defined(HAVE_KS) & defined(HAVE_QIO)

/*--------------------------------------------------------------------*/
/* Read a USQCD color vector field source.  Result is cached. */
/*--------------------------------------------------------------------*/

int r_source_vector(quark_source *qs){

  //  int x0                      = qs->x0; 
  //  int y0                      = qs->y0; 
  //  int z0                      = qs->z0; 
  //  int t0                      = qs->t0;
  int file_initialized        = qs->source_file_initialized;

  QIO_String *recxml;
  int status = 0;
  //  int rshift[4] = {x0, y0, z0, t0};
  int color = qs->color;

  /* Source file contains a list of color vectors in USQCD format,
     i.e., a succession of records with one color vector field per
     record.  Each color vector field is considered to have support
     on the entire lattice, although the record might specify values
     for just one time slice. */
  /* We read and deliver the next color vector in the file. */
    
  if(file_initialized == 0)
    r_source_open(qs);
  
  recxml = QIO_string_create();

  /* Load the color vector from the file. */
  status = qio_status(
        read_ks_vector_scidac_xml(qs->infile, qs->v_src, 1, recxml) );

  /* Translate the origin as requested */
  // shift_su3_vector(qs->v_src, rshift); 

  /* Verify that the requested color matches the color encoded in
     the record XML */
  if(status == 0)
    status = check_color(recxml, color);

  QIO_string_destroy(recxml);

  return status;

} /* r_source_vector */
#endif

#if defined(HAVE_DIRAC) & defined(HAVE_QIO)

/*--------------------------------------------------------------------*/
/* Read a USQCD Dirac vector source.  Result is cached. */
/*--------------------------------------------------------------------*/

int r_source_dirac(quark_source *qs){

  //int x0                      = qs->x0; 
  //int y0                      = qs->y0; 
  //int z0                      = qs->z0; 
  //int t0                      = qs->t0;
  int file_initialized        = qs->source_file_initialized;

  QIO_String *recxml;
  int status = 0;
  //int rshift[4] = {x0, y0, z0, t0};

  /* Source file contains a list of Dirac vectors in USQCD format,
     i.e., a succession of records with one Dirac vector field per
     record.  Each Dirac vector field is considered to have support
     on the entire lattice, although the record might specify values
     for just one time slice. */
  /* We read and deliver the next Dirac vector in the file. */
  
  if(file_initialized == 0)
    r_source_open(qs);
  
  recxml = QIO_string_create();
  
  /* Load the Dirac vector from the file. */
  status = qio_status(
	read_w_vector_scidac_xml(qs->infile, qs->wv_src, 1, recxml) );
  
  /* Translate the origin as requested */
  //shift_wilson_vector(qs->wv_src, rshift);
  
  /* Verify that the requested color and spin matches the color encoded in
     the record XML */
  if(status == 0){
    int color = convert_ksource_to_color(qs->ksource);
    int spin  = convert_ksource_to_spin(qs->ksource);
    status = check_color_spin(recxml, color, spin);
  }
  
  QIO_string_destroy(recxml);

  return status;
} /* r_source_dirac */

#endif

/*--------------------------------------------------------------------*/
/* Routines for writing sources */
/*--------------------------------------------------------------------*/

#ifdef HAVE_KS

/*--------------------------------------------------------------------*/
/* Open a color vector source file for writing */
/*--------------------------------------------------------------------*/

int w_source_open_ks(quark_source *qs, char *fileinfo){

  char myname[] = "w_source_open_ks";

#ifdef HAVE_QIO
  int volfmt, serpar;
  char *source_file = qs->save_file;

  interpret_usqcd_ks_save_flag(&volfmt, &serpar, qs->saveflag);

  if(qs->save_file_initialized != 0){
    node0_printf("%s: file already opened\n",myname);
    return 0;
  }

  if(qs->savetype == VECTOR_FIELD_FILE || 
     qs->savetype == VECTOR_FIELD_STORE){
    qs->outfile = w_open_ks_vector_scidac_file(source_file, fileinfo,
					       volfmt, serpar);
    if(qs->outfile == NULL)return 1;
    qs->save_file_initialized = 1;
  }
  else {
    node0_printf("%s: bad source type %d\n",myname, qs->type);
    return 1;
  }
#else
  node0_printf("%s: QIO compilation required for this source\n", myname);
  terminate(1);
#endif
  return 0;
} /* w_source_open_ks */

#endif

#ifdef HAVE_DIRAC

/*--------------------------------------------------------------------*/
/* Open a Dirac vector source file for writing */
/*--------------------------------------------------------------------*/

int w_source_open_dirac(quark_source *qs, char *fileinfo){

  char myname[] = "w_source_open_dirac";

#ifdef HAVE_QIO
  int volfmt, serpar;
  char *source_file = qs->save_file;

  interpret_usqcd_w_save_flag(&volfmt, &serpar, qs->saveflag);

  if(qs->save_file_initialized != 0){
    node0_printf("%s: file already opened\n",myname);
    return 0;
  }

  if(qs->savetype == DIRAC_FIELD_FILE ||
     qs->savetype == DIRAC_FIELD_STORE){
    qs->outfile = w_open_w_vector_scidac_file(source_file, fileinfo,
					      volfmt, serpar);
    if(qs->outfile == NULL)return 1;
    qs->save_file_initialized = 1;
  } 

#ifdef HAVE_KS
  else if(qs->savetype == VECTOR_FIELD_FILE ||
	  qs->savetype == VECTOR_FIELD_STORE){
    /* In some cases the Dirac source is derived from a complex source
       or a color vector source, so we can write it as a KS vector
       source instead */
    qs->outfile = w_open_ks_vector_scidac_file(source_file, fileinfo,
					       volfmt, serpar);
    if(qs->outfile == NULL)return 1;
    qs->save_file_initialized = 1;
  }
#endif
    
  else {
    node0_printf("%s: bad source type %d\n",myname, qs->type);
    return 1;
  }
#else
  node0_printf("%s: QIO compilation required for this source\n", myname);
  terminate(1);
#endif
  return 0;
} /* w_source_open_dirac */

#endif

/*--------------------------------------------------------------------*/
/* C.ose either a color vector or Dirac vector source file after writing */
/*--------------------------------------------------------------------*/

void w_source_close(quark_source *qs){

  char myname[] = "w_source_close";

  if(qs->save_file_initialized == 0){
    node0_printf("%s: file has not been opened\n",myname);
    return;
  }

  if(0);
#ifdef HAVE_QIO
#ifdef HAVE_KS
  else if(qs->savetype == VECTOR_FIELD_FILE || 
     qs->savetype == VECTOR_FIELD_STORE){
    w_close_ks_vector_scidac_file(qs->outfile);
    qs->save_file_initialized = 0;
    qs->outfile = NULL;
  }
#endif
#ifdef HAVE_DIRAC
  else if(qs->savetype == DIRAC_FIELD_FILE || 
	  qs->savetype == DIRAC_FIELD_STORE){
    w_close_w_vector_scidac_file(qs->outfile);
    qs->save_file_initialized = 0;
    qs->outfile = NULL;
  }
#endif
  else
    node0_printf("%s: bad source type %d\n",myname, qs->savetype);
#else
  node0_printf("%s: QIO compilation required for this source\n", myname);
  terminate(1);
#endif
} /* w_source_close */

#ifdef HAVE_KS

/*--------------------------------------------------------------------*/
/* Write a color vector source record */
/*--------------------------------------------------------------------*/

int w_source_ks(su3_vector *src, quark_source *qs)
{
  char myname[] = "w_source_ks";
  
  int status = 0;
  double dtime = 0;
  
  /* Unpack structure */
  int source_type           = qs->savetype;
  int file_initialized      = qs->save_file_initialized;
#ifdef HAVE_QIO
  int t0                    = qs->t0;
  int color                 = qs->color;
  QIO_String *recxml;
  QIO_USQCDKSPropRecordInfo *recinfo;
#endif
  
  dtime = -dclock();

  if(source_type != VECTOR_FIELD_FILE && source_type != VECTOR_FIELD_STORE){
    node0_printf("%s: Don't know how to write source type %d\n", myname, source_type);
    return 1;
  }

  if(file_initialized == 0){
    node0_printf("%s: File has not been initialized\n", myname);
    return 1;
  }
#ifdef HAVE_QIO
  recxml = QIO_string_create();

  /* Construct the record XML */
  recinfo = QIO_create_usqcd_ksproprecord_c_info(color, "");
  QIO_encode_usqcd_ksproprecord_info(recxml, recinfo);
  QIO_destroy_usqcd_ksproprecord_info(recinfo);
  /* Write the SU(3) source vector to the file */
  status = qio_status(
      write_kspropsource_V_usqcd_xml(qs->outfile, recxml, src, t0) );
  node0_printf("Wrote source for color %d time slice %d\n", color, t0);
  QIO_string_destroy(recxml);
  dtime += dclock();
#ifdef IOTIME
  node0_printf("Time to save source color %d = %e\n", color,dtime);
#endif
#else
  node0_printf("%s: QIO compilation required for this operation\n", myname);
  terminate(1);
#endif
  return status;
} /* w_source_ks */

#endif

#ifdef HAVE_DIRAC

/*--------------------------------------------------------------------*/
/* Write a Dirac vector source record */
/*--------------------------------------------------------------------*/

int w_source_dirac(wilson_vector *src, quark_source *qs)
{
  char myname[] = "w_source_dirac";
  
  int status = 0;
  int spin = convert_ksource_to_spin(qs->ksource);
  int color = convert_ksource_to_color(qs->ksource);
  double dtime = 0;
  
  /* Unpack structure */
  int save_type             = qs->savetype;
  int file_initialized      = qs->save_file_initialized;
#ifdef HAVE_QIO
  int t0                    = qs->t0;
  QIO_String *recxml;
#endif
  
  if(save_type != DIRAC_FIELD_FILE && save_type != VECTOR_FIELD_FILE){
    node0_printf("%s: Unrecognized source type for writing\n", myname);
    return 1;
  }

  if(file_initialized == 0){
    node0_printf("%s: File has not been initialized\n", myname);
    return 1;
  }
#ifdef HAVE_QIO

  if(save_type == DIRAC_FIELD_FILE){
    QIO_USQCDPropRecordInfo *recinfo;

    dtime = -dclock();

    recxml = QIO_string_create();
    
    /* Construct the record XML - we borrow the USQCD prop record XML */
    recinfo = QIO_create_usqcd_proprecord_sc_info(spin, color, "");
    QIO_encode_usqcd_proprecord_info(recxml, recinfo);
    QIO_destroy_usqcd_proprecord_info(recinfo);
    
    /* Write the Dirac source vector to the file */
    status = qio_status(
			write_wpropsource_D_usqcd_xml(qs->outfile, recxml, src, t0) );
    node0_printf("Wrote source for color %d spin %d time slice %d\n", 
		 color, spin, t0);
    
    QIO_string_destroy(recxml);

    dtime += dclock();
#ifdef IOTIME
    node0_printf("Time to save source %d (spin %d color %d) = %e\n",
		 qs->ksource, spin, color, dtime);
#endif
  } else {

#ifdef HAVE_KS

    QIO_USQCDKSPropRecordInfo *recinfo;

    /* Save as a color vector source */

    /* We can do this only if we have cached an intermediate color vector source */

    if(qs->v_src == NULL){
      node0_printf("%s: Can't write this source as a KS vector source\n", myname);
      status = 1;
    } else {

      /* Do this only at spin 0 so we don't write four copies of the same color vector source */
      
      if(spin == 0){
	
	dtime = -dclock();
	
	recxml = QIO_string_create();
	
	/* Construct the record XML - we borrow the USQCD prop record XML */
	recinfo = QIO_create_usqcd_ksproprecord_c_info(color, "");
	QIO_encode_usqcd_ksproprecord_info(recxml, recinfo);
	QIO_destroy_usqcd_ksproprecord_info(recinfo);
	
	/* Write the color source vector to the file */
	status = qio_status(
			    write_kspropsource_V_usqcd_xml(qs->outfile, recxml, qs->v_src, t0) );
	node0_printf("Wrote source for color %d time slice %d\n",  color, t0);
	
	QIO_string_destroy(recxml);
	
	dtime += dclock();
#ifdef IOTIME
	node0_printf("Time to save source %d (color %d) = %e\n",
		     qs->ksource, color, dtime);
#endif
	
      }
    }
#endif /* HAVE_KS */
  }
    
#else
  node0_printf("%s: QIO compilation required for this operation\n", myname);
  terminate(1);
#endif
  return status;
} /* w_source_dirac */


/*--------------------------------------------------------------------*/
/* Deprecated */
/*--------------------------------------------------------------------*/

int w_source_dirac_site(field_offset src, quark_source *qs)
{
  wilson_vector *t_src;
  int status;

  t_src  = create_wv_field();
  status = w_source_dirac(t_src, qs);
  free(t_src);

  return status;

} /* w_source_dirac_site */

#endif

/*--------------------------------------------------------------------*/
/* Name and format for an independent output source file              */
/*--------------------------------------------------------------------*/

#define IF_OK if(status==0)

static void print_output_quark_source_choices(void){
  printf("'forget_source' or ");
  printf("'save_serial_scidac_ks_source' or ");
  printf("'save_parallel_scidac_ks_source' or ");
  printf("'save_multifile_scidac_ks_source' or ");
  printf("'save_partfile_scidac_ks_source' or");
  printf("'save_serial_scidac_w_source', or ");
  printf("'save_parallel_scidac_w_source', or ");
  printf("'save_multifile_scidac_w_source', or ");
  printf("'save_partfile_scidac_w_source'? ");
}

int ask_output_quark_source_file( FILE *fp, int prompt, 
				  int *flag, int *save_type,
				  int *t0, char *descrp, char *filename)
{
  char *savebuf;
  int status = 0;
  char myname[] = "ask_output_quark_source_file";

  filename[0] = '\0';  /* Set NULL default */

  if (prompt==1){
    print_output_quark_source_choices();
  }

  savebuf = get_next_tag(fp, "output quark source command", myname);
  if (savebuf == NULL)return 1;

  printf("%s ",savebuf);

  if(strcmp("forget_source",savebuf) == 0 ) {
    *flag=FORGET;
    *save_type = UNKNOWN;
    strcpy(descrp,"");
  }
  else if(strcmp("save_serial_scidac_ks_source",savebuf) == 0 ) {
#ifdef HAVE_QIO
    *flag=SAVE_SERIAL_SCIDAC;
    *save_type = VECTOR_FIELD_FILE;
    strcpy(descrp,"vector_field");
#else
    node0_printf("requires QIO compilation!\n");
    terminate(1);
#endif
  }
  else if(strcmp("save_parallel_scidac_ks_source",savebuf) == 0 ) {
#ifdef HAVE_QIO
    *flag=SAVE_PARALLEL_SCIDAC;
    *save_type = VECTOR_FIELD_FILE;
    strcpy(descrp,"vector_field");
#else
    node0_printf("requires QIO compilation!\n");
    terminate(1);
#endif
  }
  else if(strcmp("save_multifile_scidac_ks_source",savebuf) == 0 ) {
#ifdef HAVE_QIO
    *flag=SAVE_MULTIFILE_SCIDAC;
    *save_type = VECTOR_FIELD_FILE;
    strcpy(descrp,"vector_field");
#else
    node0_printf("requires QIO compilation!\n");
    terminate(1);
#endif
  }
  else if(strcmp("save_partfile_scidac_ks_source",savebuf) == 0 ) {
#ifdef HAVE_QIO
    *flag=SAVE_PARTFILE_SCIDAC;
    *save_type = VECTOR_FIELD_FILE;
    strcpy(descrp,"vector_field");
#else
    node0_printf("requires QIO compilation!\n");
    terminate(1);
#endif
  }
  else if(strcmp("save_serial_scidac_w_source",savebuf) == 0 ) {
#ifdef HAVE_QIO
    *flag=SAVE_SERIAL_SCIDAC;
    *save_type = DIRAC_FIELD_FILE;
    strcpy(descrp,"dirac_field");
#else
    node0_printf("requires QIO compilation!\n");
    terminate(1);
#endif
  }
  else if(strcmp("save_parallel_scidac_w_source",savebuf) == 0 ) {
#ifdef HAVE_QIO
    *flag=SAVE_PARALLEL_SCIDAC;
    *save_type = DIRAC_FIELD_FILE;
    strcpy(descrp,"dirac_field");
#else
    node0_printf("requires QIO compilation!\n");
    terminate(1);
#endif
  }
  else if(strcmp("save_multifile_scidac_w_source",savebuf) == 0 ) {
#ifdef HAVE_QIO
    *flag=SAVE_MULTIFILE_SCIDAC;
    *save_type = DIRAC_FIELD_FILE;
    strcpy(descrp,"dirac_field");
#else
    node0_printf("requires QIO compilation!\n");
    terminate(1);
#endif
  }
  else if(strcmp("save_partfile_scidac_w_source",savebuf) == 0 ) {
#ifdef HAVE_QIO
    *flag=SAVE_PARTFILE_SCIDAC;
    *save_type = DIRAC_FIELD_FILE;
    strcpy(descrp,"dirac_field");
#else
    node0_printf("requires QIO compilation!\n");
    terminate(1);
#endif
  }
  else{
    printf("\n%s: ERROR IN INPUT: source file command %s is invalid\n",
	   myname, savebuf); 
    printf("Choices are \n");
    print_output_quark_source_choices();
    return 1;
  }
  
  /* Get file name */
  if( *flag != FORGET ){
    if(prompt==1)printf("enter filename\n");
    if(scanf("%s",filename) != 1){
      printf("\n%s(%d): ERROR IN INPUT: Can't read filename\n",
	     myname, this_node); 
      status++;
    }
    else
      printf("%s\n",filename);
  } else {
    printf("\n");
  }

  /* Get time slice if requested */
  if(t0 != NULL && *flag != FORGET){
    char t0_string[5];
    IF_OK status += get_s(stdin, prompt, "t0", t0_string);
    /* Provide for writing all time slices */
    IF_OK {
      if(strcmp(t0_string, "all") == 0)
	*t0 = ALL_T_SLICES;
      else
	*t0 = atoi(t0_string);
    }
  }

  return status;
} /* ask_output_quark_source_file */



