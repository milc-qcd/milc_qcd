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
open_input_usqcd_ksprop_file(ks_prop_file *kspf, int serpar)
{
  /* Open the file and read the header */
  kspf->infile = open_usqcd_ksprop_read(kspf->filename, serpar, &kspf->info);
} /* open_input_usqcd_ksprop_file */

/*---------------------------------------------------------------*/
/* Read an opened USQCD KS propagator file.  Read source and record
   for the given source color */
static int 
read_usqcd_ksprop_record(ks_prop_file *kspf, int color, su3_vector *src,
			 su3_vector *prop, quark_source *ksqs )
{
  int status = 0;
  int input_color;
  char myname[] = "read_usqcd_ksprop_record";
  int file_type = kspf->file_type;

  if(file_type == FILE_TYPE_KS_USQCD_VV_PAIRS){
    /* Read a color vector source field into a temporary vector */
    /* (Really need to support skipping a record if we don't want to read it.) */
    su3_vector *tmp_src = create_v_field();
    status = qio_status(read_kspropsource_V_usqcd(kspf->infile, ksqs->descrp, 
				  MAXDESCRP, tmp_src));
    if(status == 0){node0_printf("Read prop source %s from %s\n",ksqs->descrp, kspf->filename);}
    else if(status == -1){
      node0_printf("%s: Unexpected EOF encountered on %s\n", myname, kspf->filename);
    }
    /* For source type VECTOR_PROPAGATOR_FILE, copy src from the file */
    if(ksqs->type == VECTOR_PROPAGATOR_FILE){
      copy_v_field(src, tmp_src);
      if(status == 0){node0_printf("%s: Read prop source %s from %s\n",
				   myname, ksqs->descrp, kspf->filename);}
    }
    destroy_v_field(tmp_src);
  }

  /* Next, read the propagator for this source */
  status = qio_status(read_ksproprecord_usqcd(kspf->infile, &input_color, 
					      prop));

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

/*---------------------------------------------------------------*/
/* Translate MILC flag to USQCD mode parameter */
static int
interpret_usqcd_ks_reload_flag(int flag)
{
  switch(flag){
  case RELOAD_PARALLEL:   
    return QIO_PARALLEL;
  case RELOAD_SERIAL:
    return QIO_SERIAL;
  default:
    printf("interpret_usqcd_ks_reload_flag: bad reload flag %d\n", flag);
    terminate(1);
    return 0;
  }

  /* (The volume format is determined by looking at the file) */
}
#endif  

/*---------------------------------------------------------------*/
/* read the lattice dimensions from a binary ks prop file */

int 
read_lat_dim_ksprop(const char *filename, int file_type, int *ndim, int dims[])
{
  ks_prop_file *kspf;
  int i;

  switch(file_type){
  case FILE_TYPE_KS_USQCD_VV_PAIRS:
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
r_open_ksprop(int flag, const char *filename)
{
  ks_prop_file *kspf = NULL;
  int file_type, color;
  su3_vector *ksp;
  char myname[] = "r_open_ksprop";

  /* No file */
  if(flag == FRESH)return NULL;

  /* ASCII file */
  if(flag == RELOAD_ASCII){
    kspf = r_ascii_ks_i(filename);
    return kspf;
  }

  /* Interpret non-ASCII file type */
  file_type = get_file_type(filename);
  if(file_type == FILE_TYPE_UNKNOWN){
    node0_printf("%s: unrecognized type file %s\n", myname, filename);
    return NULL;
  }

  if(file_type == FILE_TYPE_KS_USQCD_VV_PAIRS){ /* USQCD format */

#ifdef HAVE_QIO
    /* Create a kspf structure. (No file movement here.) */
    int serpar = interpret_usqcd_ks_reload_flag(flag);
    kspf = create_input_ksprop_file_handle(filename);
    open_input_usqcd_ksprop_file(kspf, serpar);
    kspf->file_type = file_type;
    if(kspf->infile == NULL){
      printf("r_open_ksprop: Failed to open %s for reading\n", filename);
      terminate(1);
    }
#else
    node0_printf("%s: %s looks like a QIO file, but to read it requires QIO compilation\n", myname, filename);
#endif

  } else {

    if(file_type == FILE_TYPE_FM){ /* Old MILC format */
      kspf = r_serial_ks_fm_i(filename);
      kspf->file_type = file_type;
    } else {
      node0_printf("%s: File %s with file_type %d is not a supported KS propagator file (see include/file_types.h)\n", myname, filename, file_type);
    }
    
  }

  return kspf;
} /* r_open_ksprop */

/*---------------------------------------------------------------*/
/* Open propagator file for writing an SU(3) vector for one
   source color at a time. */

ks_prop_file *
w_open_ksprop(int flag, const char *filename, int source_type)
{
  ks_prop_file *kspf = NULL;
  su3_vector *ksp;
  int file_type;

  switch(flag){

  case FORGET:
    kspf = NULL;
    break;

  case SAVE_ASCII:
    kspf = w_ascii_ks_i(filename);
    break;

  case SAVE_SERIAL_SCIDAC:
  case SAVE_PARALLEL_SCIDAC:   
  case SAVE_MULTIFILE_SCIDAC: 
  case SAVE_PARTFILE_SCIDAC:

#ifdef HAVE_QIO
    {    int volfmt, serpar;
      kspf = create_output_ksprop_file_handle();
      interpret_usqcd_ks_save_flag(&volfmt, &serpar, flag);
      file_type = choose_usqcd_ks_file_type(source_type);
      kspf->file_type = file_type;
      char *fileinfo = create_ks_XML();
      kspf->outfile = open_usqcd_ksprop_write(filename, volfmt, serpar, 
					      QIO_ILDGNO,  NULL, 
					      file_type, fileinfo);
      if(kspf->outfile == NULL){
	node0_printf("w_open_ksprop: Cannot open %s for writing\n",filename);
	terminate(1);
      }
      free_ks_XML(fileinfo);
    }
    
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
    if(kspf->file_type == FILE_TYPE_KS_USQCD_VV_PAIRS)
      r_serial_ks_f(kspf);
    else if(kspf->file_type == FILE_TYPE_FM)
      r_serial_ks_fm_f(kspf);
    else
      node0_printf("r_close_ksprop: ERROR: Can't close the file %s\n",
		   kspf->filename);
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

  if(kspf == NULL)return;

  switch(flag){
  case FORGET:
    break;
  case SAVE_ASCII:
    w_ascii_ks_f(kspf);
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
			  su3_vector *src, su3_vector *dest, int timing)
{
  /* 0 normal exit value
     1 read error */

  double dtime = 0;
  size_t i;
  int status;
  site *s;
  su3_vector *prop;
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
#ifdef HAVE_QIO
  case RELOAD_SERIAL:
    prop = kspf->prop;
    file_type = kspf->file_type;
    if(file_type == FILE_TYPE_KS_USQCD_VV_PAIRS){ /* USQCD format */
      /* Read the propagator record */
      status = read_usqcd_ksprop_record(kspf, color, src, dest, ksqs);
    } else {
      node0_printf("%s: Unsupported file type %d\n", myname, file_type);
      status = 1; /* Error status */
    }
    break;
  case RELOAD_PARALLEL:
    file_type = kspf->file_type;
    if(file_type == FILE_TYPE_KS_USQCD_VV_PAIRS){
      status = read_usqcd_ksprop_record(kspf, color, src, dest, ksqs);
    } else {
      node0_printf("%s: Unsupported file type %d.\n",
		   myname, file_type);
      status = 1;
    }
    break;

  default:
    node0_printf("%s: Unrecognized reload flag.\n", myname);
    terminate(1);
#endif
  }
  
#ifndef HAVE_QIO
  /* No QIO */
  file_type = kspf->file_type;
  if(file_type == FILE_TYPE_FM){ /* Old MILC format */
    r_serial_ks_to_field(kspf, color, dest);
  } else {
    node0_printf("%s: Unsupported file type %d. Recompile with QIO??\n", myname, file_type);
    status = 1; /* Error status */
  }
#endif

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
			  quark_source *ksqs, int color,
			  su3_vector *src, su3_vector *prop, 
			  char *recinfo, int timing)
{
  double dtime = 0;
  int status;
  int i; site *s;
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
    w_ascii_ks(kspf, color, prop);
    break;
  case SAVE_SERIAL_SCIDAC:
  case SAVE_PARALLEL_SCIDAC:   
  case SAVE_MULTIFILE_SCIDAC: 
  case SAVE_PARTFILE_SCIDAC:

#ifdef HAVE_QIO
    file_type = kspf->file_type;
    /* Save color vector source field */
    if(file_type == FILE_TYPE_KS_USQCD_VV_PAIRS){
      status = write_kspropsource_V_usqcd(kspf->outfile, ksqs->descrp, 
				  src, ksqs->t0);
    } else {
      node0_printf("%s: Unsupported file type %d.\n", myname, file_type);
      status = 1;
    }
    /* Save solution field */
    if(status == 0)
      status = write_ksprop_usqcd_c(kspf->outfile, prop, color, recinfo);
#else
    node0_printf("%s: SciDAC formats require QIO compilation\n",myname);
    terminate(1);
#endif
    break;
  default:
    node0_printf("%s: Unrecognized save flag.\n", myname);
    status = 1;
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
reload_ksprop_to_field3( int flag, const char *filename, quark_source *ksqs,
			 su3_vector *dest, int timing)
{
  /* 0 normal exit value
     1 read error */

  ks_prop_file *kspf;
  int i,color,status;
  site *s;

  if(flag == FRESH)return 0;
  
  kspf = r_open_ksprop(flag, filename);
  if(kspf == NULL)return 1;

  status = 0;
  su3_vector *prop = create_v_field();
  for(color = 0; color < 3; color++){
    /* Here src = NULL, because we don't read the source */
    status = reload_ksprop_c_to_field(flag, kspf, ksqs, color, NULL, prop, timing);
    if(status != 0)break;
    FORALLSITES(i,s){
      dest[3*i+color] = prop[i];
    }
  }

  r_close_ksprop(flag, kspf);

  destroy_v_field(prop);

  return status;

} /* reload_ksprop_to_field3 */

/*---------------------------------------------------------------*/
/* Reload a full nc-color propagator: FRESH, CONTINUE,
   RELOAD_ASCII, RELOAD_SERIAL, RELOAD_PARALLEL

   dest contains an array of nc field pointers, one for each color
*/
int 
reload_ksprop_to_ksp_field( int flag, const char *filename, quark_source *ksqs,
			    ks_prop_field *source, ks_prop_field *prop, int timing)
{
  /* 0 normal exit value
     1 read error */
  
  ks_prop_file *kspf;
  int color,status;
  
  if(flag == FRESH)return 0;
  
  kspf = r_open_ksprop(flag, filename);
  if(kspf == NULL)return 1;
  
  status = 0;
  for(color = 0; color < prop->nc; color++){
    status = reload_ksprop_c_to_field(flag, kspf, ksqs, color, source->v[color],
				      prop->v[color], timing);
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
reload_ksprop_to_site3( int flag, const char *filename, quark_source *ksqs,
			field_offset dest, int timing)
{
  /* 0 normal exit value
     1 read error */

  int i,status;
  site *s;
  su3_vector *prop, *vec;

  prop = (su3_vector *)malloc(sites_on_node * 3 * sizeof(su3_vector));
  if(prop == NULL){
    printf("reload_ksprop_to_site(%d): no room for prop\n",this_node);
    return 1;
  }

  status = reload_ksprop_to_field3( flag, filename, ksqs, prop, timing);

  FORALLSITES(i,s){
    vec = (su3_vector *)F_PT(s,dest);
    vec[0] = prop[3*i];
    vec[1] = prop[3*i+1];
    vec[2] = prop[3*i+2];
  }

  free(prop);
  return status;

} /* reload_ksprop_to_site */

/*---------------------------------------------------------------*/
/* save a three-source-color KS propagator to a file in
   various formats FORGET, SAVE_ASCII, SAVE_SERIAL_TSLICE,
   SAVE_XXX_SCIDAC

   src has a triplet of color vectors for each site
   */
void 
save_ksprop_from_field3( int flag, const char *filename, char *recxml, 
			 quark_source *ksqs,
			 su3_vector *src, int timing)
{
  ks_prop_file *kspf;
  int i, color, status;
  site *s;
  su3_vector *prop;

  prop = create_v_field();

  kspf = w_open_ksprop(flag, filename, ksqs->type);
  if(kspf == NULL)return;

  status = 0;
  for(color = 0; color < 3; color++){
    FORALLSITES(i,s){
      prop[i] = src[3*i+color];
    }
    status = save_ksprop_c_from_field(flag, kspf, ksqs, color, NULL,
				      prop, recxml, timing);
    if(status != 0)break;
  }

  w_close_ksprop(flag, kspf);

  destroy_v_field(prop);

} /* save_ksprop_from_field3 */

/*---------------------------------------------------------------*/
/* save an nc-source-color KS propagator to a file in
   various formats FORGET, SAVE_ASCII, SAVE_SERIAL_TSLICE,
   SAVE_XXX_SCIDAC

   src is an array of three color fields
   */
int
save_ksprop_from_ksp_field( int flag, const char *filename, char *recxml, 
			    quark_source *ksqs, ks_prop_field *source,
			    ks_prop_field *prop, int timing)
{
  ks_prop_file *kspf;
  int  color, status;

  if(flag == FORGET)return 0;

  kspf = w_open_ksprop(flag, filename, ksqs->type);
  if(kspf == NULL)return 1;

  status = 0;
  for(color = 0; color < prop->nc; color++){
    /* Dummy source (NULL) has no colors. Use prop instead */ 
    if (source == NULL) {
    status = save_ksprop_c_from_field(flag, kspf, ksqs, color, prop->v[color],
				      prop->v[color], recxml, timing);
    }
    else {
    status = save_ksprop_c_from_field(flag, kspf, ksqs, color, source->v[color],
                                      prop->v[color], recxml, timing);
    }
    if(status != 0)break;
  }
  
  w_close_ksprop(flag, kspf);

  return status;

} /* save_ksprop_from_field3 */

/*---------------------------------------------------------------*/
/* save a three-source-color KS propagator in the site structure to a
   file in various formats FORGET, SAVE_ASCII, SAVE_SERIAL_TSLICE,
   SAVE_XXX_SCIDAC

   DEPRECATED.  Kept for compatibility.

   src has type su3_vector[3]

   */
void 
save_ksprop_from_site3( int flag, const char *filename, char *recxml, 
			quark_source *ksqs,
			field_offset dest, int timing)
{
  int i, color;
  site *s;
  su3_vector *prop, *vec;

  prop = (su3_vector *)malloc(3 * sites_on_node * sizeof(su3_vector));
  if(prop == NULL){
    printf("save_ksprop_from_site(%d): no room for prop\n",this_node);
    return;
  }

  FORALLSITES(i,s){
    for(color = 0; color < 3; color++){
      vec = (su3_vector *)F_PT(s,dest);
      prop[3*i+color] = vec[color];
    }
  }

  save_ksprop_from_field3(flag, filename, recxml, ksqs, prop, timing);

  free(prop);

} /* save_ksprop_from_site3 */

/*---------------------------------------------------------------*/
/* Translate output flag to the appropriate input flag for restoring
   a propagator that was temporarily written to disk  */
int
convert_outflag_to_inflag_ksprop(int outflag){
  switch(outflag){
  case SAVE_ASCII:
    return RELOAD_ASCII;
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
  const char *savebuf;
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
    node0_printf("'save_serial_scidac_ksprop', ");
    node0_printf("'save_parallel_scidac_ksprop', 'save_multifile_scidac_ksprop', ");
    node0_printf("'save_partfile_scidac_ksprop'");
}

int 
ask_ending_ksprop( FILE *fp, int prompt, int *flag, char *filename ){
  const char *savebuf;
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
