/********************** io_helpers_ks.c *********************************/
/* MIMD version 7 */
/* Started 17 April 2002 by MBW, derived from io_helpers_w.c,
   pretty much cut-paste-query-replace.  Someday make one set of routines
   for lattices, Wilson, and KS props? */
/*
   General purpose high level propagator I/O routines, 
   to be used by any application that wants them.
*/

#include "generic_ks_includes.h"
#include "../include/io_lat.h"
#include "../include/io_ksprop.h"
#include "../include/file_types.h"
#include <string.h>
#ifdef HAVE_QIO
#include "../include/io_scidac.h"
#include "../include/io_scidac_ks.h"
#include <qio.h>
#endif


#ifdef HAVE_QIO
/*---------------------------------------------------------------*/
/* Set up a USQCD KS propagator file for reading */

static void
open_input_usqcd_ksprop_file(ks_prop_file *kspf)
{
  /* Open the file and read the header */
  kspf->infile = open_usqcd_ksprop_read(kspf->filename,QIO_SERIAL, &kspf->info);
} /* open_input_usqcd_ksprop_file */

/*---------------------------------------------------------------*/
/* Read an opened USQCD KS propagator file and read source and record
   according to file type */
static int 
read_usqcd_ksprop_record(ks_prop_file *kspf, 
			 int color, su3_vector *dest,
			 quark_source *ksqs )
{
  int status = 0;
  int input_color;
  char myname[] = "read_usqcd_ksprop_record";
  int file_type = kspf->file_type;

  if(file_type == FILE_TYPE_KS_USQCD_CV_PAIRS || 
     (file_type == FILE_TYPE_KS_USQCD_C1V3 && color == 0)){
    /* Read a complex source field into the source cache */
    alloc_cached_c_source(ksqs);
    status = qio_status(read_kspropsource_C_usqcd(kspf->infile, ksqs->descrp, 
				  MAXDESCRP, get_cached_c_source(ksqs)));
    if(status == 0)node0_printf("Read prop source %s\n",ksqs->descrp);
    ksqs->type = COMPLEX_FIELD_STORE;
  }
  else if(file_type == FILE_TYPE_KS_USQCD_VV_PAIRS){
    /* Read a color vector source field into the source cache */
    alloc_cached_v_source(ksqs);
    status = qio_status(read_kspropsource_V_usqcd(kspf->infile, ksqs->descrp, 
				  MAXDESCRP, get_cached_v_source(ksqs)));
    if(status == 0)node0_printf("Read prop source %s\n",ksqs->descrp);
    ksqs->type = VECTOR_FIELD_STORE;
  }

  /* Next, read the solution vector */
  status = qio_status(read_ksproprecord_usqcd(kspf->infile, &input_color, 
					      dest));

  /* Check color */
  if(status == 0 && input_color != color){
    node0_printf("%s: input color %d does not match requested color %d\n",
		 myname, input_color, color);
    status = 1;
  }

  return status;
} /* read_usqcd_ksprop_record */

/*---------------------------------------------------------------*/
/* Get the SciDAC format and mode parameters from the MILC I/O flag */
void 
interpret_usqcd_ks_save_flag(int *volfmt, int *serpar, int flag)
{
  switch(flag){
  case SAVE_PARALLEL_SCIDAC:   
    *serpar = QIO_PARALLEL;
    break;
  default:
    *serpar = QIO_SERIAL;
  }
  
  switch(flag){
  case SAVE_MULTIFILE_SCIDAC: 
    *volfmt = QIO_MULTIFILE;
    break;
  case SAVE_PARTFILE_SCIDAC:
    *volfmt = QIO_PARTFILE;
    break;
  default:
    *volfmt = QIO_SINGLEFILE;
  }
}
#endif  

/*---------------------------------------------------------------*/
/* read the lattice dimensions from a binary ks prop file */

int 
read_lat_dim_ksprop(char *filename, int file_type, int *ndim, int dims[])
{
  ks_prop_file *kspf;
  int i;

  switch(file_type){
  case FILE_TYPE_KS_FMPROP:
    *ndim = 4;
    nx = -1; ny = -1; nz = -1; nt = -1;
    kspf = r_serial_ks_fm_i(filename);
    for(i = 0; i < *ndim; i++)
      dims[i] = kspf->header->dims[i];
    r_serial_ks_fm_f(kspf);
    break;
    
  case FILE_TYPE_KS_PROP:
    *ndim = 4;
    nx = -1; ny = -1; nz = -1; nt = -1;
    kspf = r_serial_ks_i(filename);
    for(i = 0; i < *ndim; i++)
      dims[i] = kspf->header->dims[i];
    r_serial_ks_f(kspf);
    break;
    
  case FILE_TYPE_KS_USQCD_C1V3:
  case FILE_TYPE_KS_USQCD_VV_PAIRS:
  case FILE_TYPE_KS_USQCD_CV_PAIRS:
#ifdef HAVE_QIO
    read_lat_dim_scidac(filename, ndim, dims);
#else
    node0_printf("read_lat_dim_ksprop(%d): This looks like a QIO file, but to read it requires QIO compilation\n",this_node);
    return 1;
#endif
    break;

 default:
   node0_printf("read_lat_dim_ksprop(%d): Unsupported file type %d\n",
		this_node,file_type);
   return 1;
  }
  return 0;
} /* read_lat_dim_ksprop */

/*----------------------------------------------------------------*/
/* Open KS propagator file for reading one source color at a time */

ks_prop_file *
r_open_ksprop(int flag, char *filename)
{
  ks_prop_file *kspf = NULL;
  int file_type, color;
  su3_vector *ksp;
  char myname[] = "r_open_ksprop";

  /* ASCII file */
  if(flag == RELOAD_ASCII){
    kspf = r_ascii_ks_i(filename);
    return kspf;
  }

  /* No file */
  if(flag == FRESH)return NULL;

  /* Interpret non-ASCII file type */
  file_type = get_file_type(filename);
  if(file_type == FILE_TYPE_UNKNOWN){
    node0_printf("%s: unrecognized type file %s\n", myname, filename);
    return NULL;
  }

  /* If it is a standard propagator file, read the full prop and cache if
     so we can process one color at a time */

  if(file_type == FILE_TYPE_KS_PROP){
    kspf = r_serial_ks_i(filename);
    kspf->file_type = file_type;
    kspf->prop = (su3_vector *)malloc(sites_on_node*sizeof(su3_vector)*3);
    ksp = kspf->prop;
    if(ksp == NULL){
      printf("%s: Can't malloc for full input propagator\n", myname);
      terminate(1);
    }
    for(color = 0; color < 3; color++)
      r_serial_ks_to_field(kspf, color, ksp);
    
    /* Indicate propagator solution data is cached */
    kspf->file_type = file_type = FILE_TYPE_KS_STORE;
  }

  /* If it is an FNAL propagator file, read the full prop and cache it 
     because we process one color at a time */

  else if(file_type == FILE_TYPE_KS_FMPROP){
    kspf = r_serial_ks_fm_i(filename);
    kspf->file_type = file_type;
    kspf->prop = (su3_vector *)malloc(sites_on_node*sizeof(su3_vector)*3);
    ksp = kspf->prop;
    if(ksp == NULL){
      printf("%s: Can't malloc for full input propagator\n", myname);
      terminate(1);
    }
    r_serial_ks_fm_to_field(kspf, ksp);
    
    /* Indicate propagator data is cached */
    kspf->file_type = file_type = FILE_TYPE_KS_STORE;
  }
#ifdef HAVE_QIO
  /* SciDAC propagator format */
  else if(file_type == FILE_TYPE_KS_USQCD_C1V3 ||
	  file_type == FILE_TYPE_KS_USQCD_VV_PAIRS ||
	  file_type == FILE_TYPE_KS_USQCD_CV_PAIRS){
    /* Create a kspf structure. (No file movement here.) */
    kspf = create_input_ksprop_file_handle(filename);
    kspf->file_type = file_type;
    open_input_usqcd_ksprop_file(kspf);
  }
#endif
  else {
    node0_printf("%s: File %s is not a KS propagator file\n", myname, filename);
  }

  return kspf;
}

/*---------------------------------------------------------------*/
/* Open propagator file for writing an SU(3) vector for one
   source color at a time. */

ks_prop_file *
w_open_ksprop(int flag, char *filename, int source_type)
{
  ks_prop_file *kspf = NULL;
  su3_vector *ksp;
#ifdef HAVE_QIO
  char *fileinfo;
  int volfmt, serpar, file_type;
#endif
  
  switch(flag){

  case FORGET:
    kspf = NULL;
    break;

  case SAVE_ASCII:
    kspf = w_ascii_ks_i(filename);
    break;

  case SAVE_SERIAL:
    kspf = w_serial_ks_i(filename);
    /* Allocate space for the entire propagator */
    kspf->prop = (su3_vector *)
      malloc(sites_on_node*sizeof(su3_vector)*3);
    ksp = kspf->prop;
    if(ksp == NULL){
      printf("Can't malloc for full output propagator\n");
      terminate(1);
    }
    break;

  case SAVE_SERIAL_FM:
    kspf = w_serial_ks_fm_i(filename);
    /* Allocate space for the entire propagator */
    kspf->prop = (su3_vector *)
      malloc(sites_on_node*sizeof(su3_vector)*3);
    ksp = kspf->prop;
    if(ksp == NULL){
      printf("Can't malloc for full output propagator\n");
      terminate(1);
    }
    break;

  case SAVE_SERIAL_SCIDAC:
  case SAVE_PARALLEL_SCIDAC:   
  case SAVE_MULTIFILE_SCIDAC: 
  case SAVE_PARTFILE_SCIDAC:

#ifdef HAVE_QIO
    kspf = create_output_ksprop_file_handle();
    interpret_usqcd_ks_save_flag(&volfmt, &serpar, flag);
    file_type = choose_usqcd_ks_file_type(source_type);
    kspf->file_type = file_type;
    fileinfo = create_ks_XML();
    kspf->outfile = open_usqcd_ksprop_write(filename, volfmt, serpar, 
					    QIO_ILDGNO,  NULL, 
					    file_type, fileinfo);
    free_ks_XML(fileinfo);
    
#else
    node0_printf("w_open_ksprop: SciDAC formats require QIO compilation\n");
    terminate(1);
#endif
    break;

  default:
    node0_printf("w_open_ksprop: Unsupported save flag %d\n",flag);
    kspf = NULL;
  }

  return kspf;
} /* w_open_ksprop */

/*---------------------------------------------------------------*/
void 
r_close_ksprop(int flag, ks_prop_file *kspf)
{
  
  if(kspf == NULL)return;

  switch(flag){
  case RELOAD_ASCII:
    r_ascii_ks_f(kspf);
    break;
  case RELOAD_SERIAL:
    r_serial_ks_f(kspf);
    break;
  case RELOAD_PARALLEL:
    destroy_ksprop_file_handle(kspf);
    break;
  default:
    node0_printf("r_close_ksprop: Unrecognized read flag %d",flag);
  }
} /* r_close_ksprop */

/*---------------------------------------------------------------*/
void 
w_close_ksprop(int flag, ks_prop_file *kspf)
{
  su3_vector *ksp;
  int color;

  if(kspf == NULL)return;

  switch(flag){
  case FORGET:
    break;
  case SAVE_ASCII:
    w_ascii_ks_f(kspf);
    break;
  case SAVE_SERIAL:
    /* Dump accumulated propagator and free memory */
    ksp = kspf->prop;
    for(color = 0; color < 3; color++)
      w_serial_ks_from_field(kspf, color, ksp);
    w_serial_ks_f(kspf); 
    break;
  case SAVE_SERIAL_FM:
    /* Dump accumulated propagator and free memory */
    ksp = kspf->prop;
    w_serial_ks_fm_from_field(kspf, ksp);
    w_serial_ks_fm_f(kspf); 
    break;
  case SAVE_SERIAL_SCIDAC:
  case SAVE_PARALLEL_SCIDAC:   
  case SAVE_MULTIFILE_SCIDAC: 
  case SAVE_PARTFILE_SCIDAC:
#ifdef HAVE_QIO
    close_usqcd_ksprop_write(kspf->outfile);
    destroy_ksprop_file_handle(kspf);
#endif
    break;
  default:
    node0_printf("w_close_ksprop: Unrecognized save flag %d\n",flag);
  }
} /* w_close_ksprop */

/*---------------------------------------------------------------*/
/* reload a propagator for a single source color and spin in MILC or
   USQCD format, or cold propagator, or keep current propagator:
   FRESH, CONTINUE, RELOAD_ASCII, RELOAD_SERIAL, RELOAD_PARALLEL,
   RELOAD_MULTIDUMP
   */
int 
reload_ksprop_c_to_field( int flag, ks_prop_file *kspf, 
			  quark_source *ksqs, int color, 
			  su3_vector *dest, int timing)
{
  /* 0 normal exit value
     1 read error */

  double dtime = 0;
  int i,status;
  site *s;
  su3_vector *ksp;
  su3_vector *cv;
  int c0;
  int file_type = FILE_TYPE_UNKNOWN;  /* So the compiler doesn't say uninit */
  char myname[] = "reload_ksprop_c_to_field";

  if(timing)dtime = -dclock();
  status = 0;
  switch(flag){
  case CONTINUE:  /* do nothing */
    break;
  case FRESH:     /* zero initial guess */
    FORALLSITES(i,s)clearvec( &(dest[i]) );
    break;
  case RELOAD_ASCII:
    r_ascii_ks(kspf, color, dest);
    break;
  case RELOAD_SERIAL:
    ksp = kspf->prop;
    file_type = kspf->file_type;
    /* Special treatment for a cached propagator */
    if(file_type == FILE_TYPE_KS_STORE){
    /* Copy input KS vector for this color from buffer */
      FORALLSITES(i,s){
	cv = dest + i;
	for(c0=0;c0<3;c0++)
	  {
	    cv->c[c0].real = ksp[3*i+color].c[c0].real;
	    cv->c[c0].imag = ksp[3*i+color].c[c0].imag;
	  }
      }
      status = 0;
    }
#ifdef HAVE_QIO
    else if(file_type == FILE_TYPE_KS_USQCD_C1V3 ||
	    file_type == FILE_TYPE_KS_USQCD_VV_PAIRS ||
	    file_type == FILE_TYPE_KS_USQCD_CV_PAIRS){

      /* Read the propagator record */
      status = read_usqcd_ksprop_record(kspf, color, dest, ksqs);
    }
#endif
    else {
      status = r_serial_ks_to_field(kspf,color,dest); 
    }
    break;
  case RELOAD_PARALLEL:
#ifdef HAVE_QIO
    if(file_type == FILE_TYPE_KS_USQCD_C1V3 ||
       file_type == FILE_TYPE_KS_USQCD_VV_PAIRS ||
       file_type == FILE_TYPE_KS_USQCD_CV_PAIRS){
      status = read_usqcd_ksprop_record(kspf, color, dest, ksqs);
    }
#endif
    /* If not SciDAC */
    if(file_type != FILE_TYPE_KS_USQCD_C1V3 &&
       file_type != FILE_TYPE_KS_USQCD_VV_PAIRS &&
       file_type != FILE_TYPE_KS_USQCD_CV_PAIRS){
      node0_printf("%s: Parallel reading with this file type not supported\n",
		   myname);
    }
    break;
  default:
    node0_printf("%s: Unrecognized reload flag.\n", myname);
    terminate(1);
  }
  
  if(timing)
    {
      dtime += dclock();
      if(flag != FRESH && flag != CONTINUE)
	node0_printf("Time to reload ksprop color %d %e\n",
		     color,dtime);
    }

  return status;

} /* reload_ksprop_c_to_field */

/*---------------------------------------------------------------*/
/* save a propagator one source color at a time */
/* recinfo is for USQCD formats */
int 
save_ksprop_c_from_field( int flag, ks_prop_file *kspf, 
			  quark_source *ksqs,
			  int color, su3_vector *src, 
			  char *recinfo, int timing)
{
  double dtime = 0;
  int status;
  int i; site *s;
  su3_vector *ksp;
  su3_vector *cv;
  int c0;
#ifdef HAVE_QIO
  int file_type;
#endif
  char myname[] = "save_ksprop_c_from_field";
  
  if(timing)dtime = -dclock();
  status = 0;
  switch(flag){
  case FORGET:
    break;
  case SAVE_ASCII:
    w_ascii_ks(kspf, color, src);
    break;
  case SAVE_SERIAL:
  case SAVE_SERIAL_FM:
    ksp = kspf->prop;
    if(ksp == NULL){
      printf("%s{%d): Propagator field not allocated\n", myname, this_node);
      terminate(1);
    }
    /* Add output KS vector to propagator buffer */
    FORALLSITES(i,s){
      cv = src + i;
      for(c0=0;c0<3;c0++)
	{
	  ksp[3*i+color].c[c0].real = cv->c[c0].real;
	  ksp[3*i+color].c[c0].imag = cv->c[c0].imag;
	}
    }
    break;
  case SAVE_SERIAL_SCIDAC:
  case SAVE_PARALLEL_SCIDAC:   
  case SAVE_MULTIFILE_SCIDAC: 
  case SAVE_PARTFILE_SCIDAC:

#ifdef HAVE_QIO
    file_type = kspf->file_type;
    /* Save color source field */
    if(file_type == FILE_TYPE_KS_USQCD_CV_PAIRS ||
       (file_type == FILE_TYPE_KS_USQCD_C1V3 && color == 0))
      {
	complex *c_src = get_cached_c_source(ksqs);
	int null_src = (c_src == NULL);
	if(null_src){
	    node0_printf("%s complex source is missing\n",myname);
	    node0_printf("%s File will be written with a dummy zero source\n",
		   myname);
	    c_src = create_c_field();
	}
	status = write_kspropsource_C_usqcd(kspf->outfile, ksqs->descrp, 
				    c_src, ksqs->t0);
	if(null_src)free(c_src);
	if(status != 0)break;
      }
    /* Save color vector source field */
    else if(file_type == FILE_TYPE_KS_USQCD_VV_PAIRS){
      su3_vector *v_src = get_cached_v_source(ksqs);
      int null_src = (v_src == NULL);
      if(null_src){
	node0_printf("%s color vector source is missing\n",myname);
	node0_printf("%s File will be written with a dummy zero source\n",
		     myname);
	v_src = create_v_field();
      }
      status = write_kspropsource_V_usqcd(kspf->outfile, ksqs->descrp, 
				  v_src, ksqs->t0);
      if(null_src)free(v_src);
      if(status != 0)break;
    }
    /* Save solution field */
    if(status == 0)
      status = write_ksprop_usqcd_c(kspf->outfile, src, color, recinfo);
#else
    node0_printf("%s: SciDAC formats require QIO compilation\n",myname);
    terminate(1);
#endif
    
    break;
  default:
    node0_printf("%s: Unrecognized save flag.\n", myname);
    terminate(1);
  }
  
  if(timing)
    {
      dtime += dclock();
      if(flag != FORGET)
	node0_printf("Time to save prop color %d = %e\n", color, dtime);
    }
  return status;
} /* save_ksprop_c_from_field */

/*---------------------------------------------------------------*/
/* Reload a full three-color propagator: FRESH, CONTINUE,
   RELOAD_ASCII, RELOAD_SERIAL, RELOAD_PARALLEL

   dest takes a triplet of color vectors per site
*/
int 
reload_ksprop_to_field3( int flag, char *filename, quark_source *ksqs,
			 su3_vector *dest, int timing)
{
  /* 0 normal exit value
     1 read error */

  ks_prop_file *kspf;
  int i,color,status;
  site *s;
  su3_vector *ksp;

  ksp = create_v_field();

  kspf = r_open_ksprop(flag, filename);
  if(kspf == NULL)return 1;

  status = 0;
  for(color = 0; color < 3; color++){
    status = reload_ksprop_c_to_field(flag, kspf, ksqs, color, ksp, timing);
    if(status != 0)break;
    FORALLSITES(i,s){
      dest[3*i+color] = ksp[i];
    }
  }

  r_close_ksprop(flag, kspf);

  destroy_v_field(ksp);

  return status;

} /* reload_ksprop_to_field3 */

/*---------------------------------------------------------------*/
/* Reload a full nc-color propagator: FRESH, CONTINUE,
   RELOAD_ASCII, RELOAD_SERIAL, RELOAD_PARALLEL

   dest contains an array of nc field pointers, one for each color
*/
int 
reload_ksprop_to_ksp_field( int flag, char *filename, quark_source *ksqs,
			    ks_prop_field *dest, int timing)
{
  /* 0 normal exit value
     1 read error */

  ks_prop_file *kspf;
  int color,status;

  kspf = r_open_ksprop(flag, filename);
  if(kspf == NULL)return 1;

  status = 0;
  for(color = 0; color < dest->nc; color++){
    status = reload_ksprop_c_to_field(flag, kspf, ksqs, color, dest->v[color], 
				      timing);
    if(status != 0)break;
  }

  r_close_ksprop(flag, kspf);

  return status;

} /* reload_ksprop_to_ksprop_field */

/*---------------------------------------------------------------*/
/* Reload a full three-color propagator to the site structure: 
   FRESH, CONTINUE, RELOAD_ASCII, RELOAD_SERIAL, RELOAD_PARALLEL

   DEPRECATED.  KEPT FOR BACKWARD COMPATIBILITY.
*/
int 
reload_ksprop_to_site3( int flag, char *filename, quark_source *ksqs,
			field_offset dest, int timing)
{
  /* 0 normal exit value
     1 read error */

  int i,status;
  site *s;
  su3_vector *ksp, *vec;

  ksp = (su3_vector *)malloc(sites_on_node * 3 * sizeof(su3_vector));
  if(ksp == NULL){
    printf("reload_ksprop_to_site(%d): no room for ksp\n",this_node);
    return 1;
  }

  status = reload_ksprop_to_field3( flag, filename, ksqs, ksp, timing);

  FORALLSITES(i,s){
    vec = (su3_vector *)F_PT(s,dest);
    vec[0] = ksp[3*i];
    vec[1] = ksp[3*i+1];
    vec[2] = ksp[3*i+2];
  }

  free(ksp);
  return status;

} /* reload_ksprop_to_site */

/*---------------------------------------------------------------*/
/* save a three-source-color KS propagator to a file in
   various formats FORGET, SAVE_ASCII, SAVE_SERIAL_TSLICE,
   SAVE_SERIAL, SERIAL_SERIAL_FM, SAVE_XXX_SCIDAC

   src has a triplet of color vectors for each site
   */
void 
save_ksprop_from_field3( int flag, char *filename, char *recxml, 
			 quark_source *ksqs,
			 su3_vector *src, int timing)
{
  ks_prop_file *kspf;
  int i, color, status;
  site *s;
  su3_vector *ksp;

  ksp = create_v_field();

  kspf = w_open_ksprop(flag, filename, ksqs->type);
  if(kspf == NULL)return;

  status = 0;
  for(color = 0; color < 3; color++){
    FORALLSITES(i,s){
      ksp[i] = src[3*i+color];
    }
    status = save_ksprop_c_from_field(flag, kspf, ksqs, color, ksp, 
				      recxml, timing);
    if(status != 0)break;
  }

  w_close_ksprop(flag, kspf);

  destroy_v_field(ksp);

} /* save_ksprop_from_field3 */

/*---------------------------------------------------------------*/
/* save a three-source-color KS propagator to a file in
   various formats FORGET, SAVE_ASCII, SAVE_SERIAL_TSLICE,
   SAVE_SERIAL, SERIAL_SERIAL_FM, SAVE_XXX_SCIDAC

   src is an array of three color fields
   */
void 
save_ksprop_from_ksp_field( int flag, char *filename, char *recxml, 
			    quark_source *ksqs,
			    ks_prop_field *src, int timing)
{
  ks_prop_file *kspf;
  int  color, status;

  kspf = w_open_ksprop(flag, filename, ksqs->type);
  if(kspf == NULL)return;

  status = 0;
  for(color = 0; color < src->nc; color++){
    status = save_ksprop_c_from_field(flag, kspf, ksqs, color, src->v[color], 
				      recxml, timing);
    if(status != 0)break;
  }
  
  w_close_ksprop(flag, kspf);

} /* save_ksprop_from_field3 */

/*---------------------------------------------------------------*/
/* save a three-source-color KS propagator in the site structure to a
   file in various formats FORGET, SAVE_ASCII, SAVE_SERIAL_TSLICE,
   SAVE_SERIAL, SERIAL_SERIAL_FM, SAVE_XXX_SCIDAC

   DEPRECATED.  Kept for compatibility.

   src has type su3_vector[3]

   */
void 
save_ksprop_from_site3( int flag, char *filename, char *recxml, 
			quark_source *ksqs,
			field_offset src, int timing)
{
  int i, color;
  site *s;
  su3_vector *ksp, *vec;

  ksp = (su3_vector *)malloc(3 * sites_on_node * sizeof(su3_vector));
  if(ksp == NULL){
    printf("save_ksprop_from_site(%d): no room for ksp\n",this_node);
    return;
  }

  FORALLSITES(i,s){
    for(color = 0; color < 3; color++){
      vec = (su3_vector *)F_PT(s,src);
      ksp[3*i+color] = vec[color];
    }
  }

  save_ksprop_from_field3(flag, filename, recxml, ksqs, ksp, timing);

  free(ksp);

} /* save_ksprop_from_site3 */

/*---------------------------------------------------------------*/
/* Translate output flag to the appropriate input flag for restoring
   a propagator that was temporarily written to disk  */
int
convert_outflag_to_inflag_ksprop(int outflag){
  switch(outflag){
  case SAVE_ASCII:
    return RELOAD_ASCII;
  case SAVE_SERIAL_FM:
  case SAVE_SERIAL_SCIDAC:
  case SAVE_MULTIFILE_SCIDAC:            
  case SAVE_PARTFILE_SCIDAC:            
    return RELOAD_SERIAL;
  case SAVE_PARALLEL_SCIDAC:             
    return RELOAD_PARALLEL;
  default:
    return FRESH;  /* Error return */
  }
}

/*---------------------------------------------------------------*/

/* find out what if any KS propagator should be loaded.
   This routine is only called by node 0.
*/
int 
ask_starting_ksprop( FILE *fp, int prompt, int *flag, char *filename ){
  char *savebuf;
  int status;
  char myname[] = "ask_starting_ksprop";
  
  if (prompt==1) {
    printf("enter 'fresh_ksprop', ");
    printf("'reload_ascii_ksprop', 'reload_serial_ksprop', ");
    printf("or 'reload_parallel_ksprop' \n");
  }

  savebuf = get_next_tag(fp, "read ksprop command", myname);
  if (savebuf == NULL)return 1;
  
  printf("%s ",savebuf);
  if(strcmp("fresh_ksprop",savebuf) == 0 ){
    *flag = FRESH;
    printf("\n");
  }
  else if(strcmp("reload_ascii_ksprop",savebuf) == 0 ) {
    *flag = RELOAD_ASCII;
  }
  else if(strcmp("reload_serial_ksprop",savebuf) == 0 ) {
    *flag = RELOAD_SERIAL;
  }
  else if(strcmp("reload_parallel_ksprop",savebuf) == 0 ) {
    *flag = RELOAD_PARALLEL;
  }
  else{
    printf("ERROR IN INPUT: propagator_command is invalid\n"); return(1);
  }

  /*read name of file and load it */
  if( *flag != FRESH && *flag != CONTINUE ){
    if(prompt==1) printf("enter name of file containing ksprop\n");
    status = scanf("%s",filename);
    if(status != 1) {
      printf("\nask_starting_ksprop: ERROR IN INPUT: Can't read filename\n");
      return(1);
    }
    printf("%s\n",filename);
  }
  return(0);
  
} /* end ask_starting_ksprop() */

/*--------------------------------------------------------------------*/

/* find out what do to with lattice at end, and lattice name if
   necessary.  This routine is only called by node 0.
*/

static void 
print_options(void)
{
    node0_printf("'forget_ksprop', 'save_ascii_ksprop', ");
    node0_printf("'save_serial_ksprop', ");
    node0_printf("'save_serial_fm_ksprop', 'save_serial_scidac_ksprop', ");
    node0_printf("'save_parallel_scidac_ksprop', 'save_multifile_scidac_ksprop', ");
    node0_printf("'save_partfile_scidac_ksprop'");
}

int 
ask_ending_ksprop( FILE *fp, int prompt, int *flag, char *filename ){
  char *savebuf;
  int status;
  char myname[] = "ask_ending_ksprop";

  if (prompt==1) {
    print_options();
    printf(" ?\n");
  }

  savebuf = get_next_tag(fp, "write ksprop command", myname);
  if (savebuf == NULL)return 1;

  printf("%s ",savebuf);

  if(strcmp("save_ascii_ksprop",savebuf) == 0 )  {
    *flag=SAVE_ASCII;
  }
  else if(strcmp("save_serial_ksprop",savebuf) == 0 ) {
    *flag=SAVE_SERIAL;
  }
  else if(strcmp("save_serial_fm_ksprop",savebuf) == 0 ) {
    *flag=SAVE_SERIAL_FM;
  }
  else if(strcmp("save_serial_scidac_ksprop",savebuf) == 0 ) {
#ifdef HAVE_QIO
    *flag=SAVE_SERIAL_SCIDAC;
#else
    node0_printf("requires QIO compilation!\n");
    terminate(1);
#endif
  }
  else if(strcmp("save_parallel_scidac_ksprop",savebuf) == 0 ) {
#ifdef HAVE_QIO
    *flag=SAVE_PARALLEL_SCIDAC;
#else
    node0_printf("requires QIO compilation!\n");
    terminate(1);
#endif
  }
  else if(strcmp("save_multifile_scidac_ksprop",savebuf) == 0 ) {
#ifdef HAVE_QIO
    *flag=SAVE_MULTIFILE_SCIDAC;
#else
    node0_printf("requires QIO compilation!\n");
    terminate(1);
#endif
  }
  else if(strcmp("save_partfile_scidac_ksprop",savebuf) == 0 ) {
#ifdef HAVE_QIO
    *flag=SAVE_PARTFILE_SCIDAC;
#else
    node0_printf("requires QIO compilation!\n");
    terminate(1);
#endif
  }
  else if(strcmp("forget_ksprop",savebuf) == 0 ) {
    *flag=FORGET;
    printf("\n");
  }
  else {
    node0_printf("is not a valid save KS prop command. INPUT ERROR.\n");
    node0_printf("Choices are ");
    print_options();
    node0_printf("\n");
    return(1);
  }
  
  if( *flag != FORGET ){
    if(prompt==1)printf("enter filename\n");
    status = scanf("%s",filename);
    if(status != 1){
      printf("\n%s(%d): ERROR IN INPUT: Can't read filename\n",
	     myname, this_node); 
      return(1);
    }
    printf("%s\n",filename);
  }
  return(0);
  
} /* end ask_ending_ksprop() */
