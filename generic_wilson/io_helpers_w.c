/********************** io_helpers_w.c **********************************/
/* MIMD version 7 */
/* CD 10/97 from io_helpers.c DT 8/97 
   General purpose high level propagator I/O routines, 
   to be used by any application that wants them.
   */

#include "generic_wilson_includes.h"
#include "../include/io_lat.h"
#include "../include/io_wprop.h"
#include "../include/file_types.h"
#include <string.h>
#ifdef HAVE_QIO
#include <qio.h>
#include "../include/io_scidac.h"
#include "../include/io_scidac_w.h"
#endif

#ifdef HAVE_QIO
/*---------------------------------------------------------------*/
/* Set up a USQCD Wilson propagator file for reading */

static void
open_input_usqcd_prop_file(w_prop_file *wpf, int serpar)
{

  /* Open the file and read the header */
  wpf->infile = r_open_usqcd_wprop_file(wpf->filename, serpar);
} /* setup_input_usqcd_prop_file */

/*---------------------------------------------------------------*/
/* Read a USQCD Wilson propagator source and record according to file type */
static int 
read_usqcd_wprop_record(w_prop_file *wpf, int spin, int color,
			wilson_vector *src, wilson_vector *prop,
			quark_source *wqs )
{
  int status = 0;
  int input_spin, input_color;
  char myname[] = "read_usqcd_wprop_record";
  int file_type = wpf->file_type;

  if(file_type == FILE_TYPE_W_USQCD_DD_PAIRS){
    /* Read a Wilson vector source field */
    wilson_vector *tmp_src = create_wv_field();
    status = qio_status(read_wpropsource_D_usqcd(wpf->infile, wqs->descrp, 
						 MAXDESCRP, tmp_src));
    if(status == 0){node0_printf("Read prop source %s from %s\n",wqs->descrp, wpf->filename);}
    else if(status == -1){
      node0_printf("Unexpected EOF encountered on %s\n", wpf->filename);
    }
    /* For source type DIRAC_PROPAGATOR_FILE, copy src from the file */
    if(wqs->type == DIRAC_PROPAGATOR_FILE){
      copy_wv_field(src, tmp_src);
      if(status == 0){node0_printf("%s: Read prop source %s from %s\n",
				   myname, wqs->descrp, wpf->filename);}
    }
    destroy_wv_field(tmp_src);
  }

  /* Next, read the propagator for this source */
  if(status == 0)
    status = qio_status(read_wproprecord_usqcd(wpf->infile, 
					      &input_spin, &input_color, 
					      prop));

  /* Check spin and color */
  if(status == 0 && (input_spin != spin || input_color != color)){
    node0_printf("%s: input spin %d and color %d do not match requested spin %d and color %d\n",
		 myname, input_spin, input_color, spin, color);
    status = 1;
  }

  return status;
} /* read_usqcd_wprop_record */

/*---------------------------------------------------------------*/
/* Get the SciDAC format and mode parameters from the MILC I/O flag */
void 
interpret_usqcd_w_save_flag(int *volfmt, int *serpar, int flag)
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
} /* interpret_usqcd_w_save_flag */

/*---------------------------------------------------------------*/
/* Translate MILC flag to USQCD mode parameter */
static int
interpret_usqcd_w_reload_flag(int flag)
{
  switch(flag){
  case RELOAD_PARALLEL:   
    return QIO_PARALLEL;
  case RELOAD_SERIAL:
    return QIO_SERIAL;
  default:
    printf("interpret_usqcd_w_reload_flag: bad reload flag %d\n", flag);
    terminate(1);
  }

  /* (The volume format is determined by looking at the file) */
  return -1;
}
#endif  

/*---------------------------------------------------------------*/
/* read the lattice dimensions from a binary Wilson prop file */

int 
read_lat_dim_wprop(char *filename, int file_type, int *ndim, int dims[])
{
  w_prop_file *wpf;
  int i;

  switch(file_type){
  case FILE_TYPE_W_USQCD_DD_PAIRS:
#ifdef HAVE_QIO
    read_lat_dim_scidac(filename, ndim, dims);
#else
    node0_printf("read_lat_dim_wprop(%d): This looks like a QIO file, but to read it requires QIO compilation\n", this_node);
    return 1;
#endif
    break;
  default:
    node0_printf("read_lat_dim_wprop(%d): Unsupported file type %d\n",
		 this_node, file_type);
    return 1;
  }
  return 0;
} /* read_lat_dim_wprop */

/*---------------------------------------------------------------*/
/* Open Wilson propagator file for reading one source spin-color at a time */

w_prop_file *
r_open_wprop(int flag, char *filename)
{
  w_prop_file *wpf = NULL;
  int file_type;
  wilson_propagator *wp;
  char myname[] = "r_open_wprop";

  /* No file */
  if(flag == FRESH)return NULL;

  if(flag == RELOAD_ASCII){
    wpf = r_ascii_w_i(filename);
    return wpf;
  }

  /* Interpret non-ASCII file type */
  file_type = get_file_type(filename);
  if(file_type == FILE_TYPE_UNKNOWN){
    node0_printf("%s: unrecognized type file %s\n", myname, filename);
    return NULL;
  }

  /* If it is an FNAL propagator file, read the full prop and cache
     it.  (We have to do a spin rotation on source and sink spins to
     convert to MILC spin conventions.) */

  /* There are two types of Fermilab Wilson propagator files.  One is
     sorted by source color and spin, so the site datum is a Dirac
     vector.  The other has the full propagator on each site. */

  if(file_type == FILE_TYPE_W_USQCD_DD_PAIRS){

#ifdef HAVE_QIO
    /* Create a wpf structure. (No file movement here.) */
    int serpar = interpret_usqcd_w_reload_flag(flag);
    wpf = setup_input_w_prop_file(filename);
    wpf->file_type = file_type;
    open_input_usqcd_prop_file(wpf, serpar);
    if(wpf->infile == NULL){
      printf("r_open_wprop: Failed to open %s for reading\n", filename);
      terminate(1);
    }
#else
    node0_printf("%s: %s looks like a QIO file, but to read it requires QIO compilation\n",
		 myname, filename);
#endif
  } else {

    if(file_type == FILE_TYPE_FM){ /* Old MILC format */
      wpf = r_serial_w_fm_i(filename);
      wpf->file_type = file_type;
      wpf->rank2rcv = NULL;
    } else {
      node0_printf("%s: File %s with file_type %d is not a supported Dirac propagator file (see include/file_types.h)\n", myname, filename, file_type);
    }
  }

  return wpf;
} /* r_open_wprop */

/*---------------------------------------------------------------*/
/* Open propagator file for writing a Wilson vector for one
   source spin and color at a time. */

w_prop_file *
w_open_wprop(int flag, char *filename, int source_type)
{
  w_prop_file *wpf = NULL;
  wilson_propagator *wp;
#ifdef HAVE_QIO
  char *fileinfo;
  int volfmt, serpar, file_type;
#endif

  switch(flag){

  case FORGET:
    wpf = NULL;
    break;

  case SAVE_ASCII:
    wpf = w_ascii_w_i(filename);
    break;

  case SAVE_SERIAL_SCIDAC:
  case SAVE_PARALLEL_SCIDAC:   
  case SAVE_MULTIFILE_SCIDAC: 
  case SAVE_PARTFILE_SCIDAC:

#ifdef HAVE_QIO
    wpf = setup_output_w_prop_file();
    interpret_usqcd_w_save_flag(&volfmt, &serpar, flag);
    file_type = choose_usqcd_w_file_type(source_type);
    wpf->file_type = file_type;
    fileinfo = create_w_QCDML();
    wpf->outfile = w_open_usqcd_wprop_file(filename, volfmt, serpar, QIO_ILDGNO,
					   NULL, file_type, fileinfo);
    if(wpf->outfile == NULL){
      node0_printf("w_open_wprop: Cannot open %s for writing\n",filename);
      terminate(1);
    }
    free_w_QCDML(fileinfo);
    
#else
    node0_printf("w_open_wprop: SciDAC formats require QIO compilation\n");
    terminate(1);
#endif
    break;

  default:
    wpf = NULL;
  }

  return wpf;
} /* w_open_wprop */

/*---------------------------------------------------------------*/
void 
r_close_wprop(int flag, w_prop_file *wpf)
{

  if(wpf == NULL)return;

  switch(flag){
  case RELOAD_ASCII:
    r_ascii_w_f(wpf);
    break;
  case RELOAD_SERIAL:
    if(wpf->file_type == FILE_TYPE_W_USQCD_DD_PAIRS)
      clear_input_w_prop_file(wpf);
    else if(wpf->file_type == FILE_TYPE_FM)
      r_serial_w_fm_f(wpf);
    else
      node0_printf("r_close_wprop: ERROR: Can't close the file %s. Don't recognize the file_type %d\n",
		   wpf->filename,wpf->file_type);
    break;
  case RELOAD_PARALLEL:
    clear_input_w_prop_file(wpf);
    break;
  default:
    node0_printf("r_close_wprop: Unrecognized read flag %d",flag);
  }
}
/*---------------------------------------------------------------*/
void 
w_close_wprop(int flag, w_prop_file *wpf)
{
  wilson_propagator *wp;

  if(wpf == NULL)return;

  switch(flag){
  case FORGET:
    break;
  case SAVE_ASCII:
    w_ascii_w_f(wpf);
    break;
  case SAVE_SERIAL_SCIDAC:
  case SAVE_PARALLEL_SCIDAC:   
  case SAVE_MULTIFILE_SCIDAC: 
  case SAVE_PARTFILE_SCIDAC:
#ifdef HAVE_QIO
    clear_output_w_prop_file(wpf);
#endif
    break;
  default:
    node0_printf("w_close_wprop: Unrecognized save flag %d\n",flag);
  }
}

/*---------------------------------------------------------------*/
/* reload a Wilson vector field for a single source color and spin in
   FNAL or USQCD format, or cold propagator, or keep current
   propagator: FRESH, CONTINUE, RELOAD_ASCII, RELOAD_SERIAL,
   RELOAD_PARALLEL

   Return value is 

  -1 end of file
   0 normal exit value
   1 read error
*/

int 
reload_wprop_sc_to_field( int flag, w_prop_file *wpf, 
			  quark_source *wqs, int spin, int color, 
			  wilson_vector *src, wilson_vector *dest, int timing)
{

  double dtime = 0;
  int i,status;
  site *s;
  wilson_propagator *wp;
  wilson_vector *wv;
  int s0, c0;
  int file_type = FILE_TYPE_UNKNOWN;  /* So the compiler doesn't say uninit */
  char myname[] = "reload_wprop_sc_to_field";

  if(timing)dtime = -dclock();
  status = 0;
  switch(flag){
  case CONTINUE:  /* do nothing */
    break;
  case FRESH:     /* zero initial guess */
    FORALLSITES(i,s)clear_wvec( &(dest[i]) );
    break;
  case RELOAD_ASCII:
    node0_printf("Reloading ASCII to wprop field not supported\n");
    terminate(1);
    break;
  case RELOAD_SERIAL:

#ifdef HAVE_QIO
    file_type = wpf->file_type;
    if(file_type == FILE_TYPE_W_USQCD_DD_PAIRS){
      /* Read the propagator record */
      status = read_usqcd_wprop_record(wpf, spin, color, src, dest, wqs);
    }
    else {
      node0_printf("%s: File type %d not supported\n",myname, file_type);
      status = 1; /*Error status */
    }
#else

    file_type = wpf->file_type;
    if(file_type == FILE_TYPE_FM){
      status = r_serial_w_to_field(wpf, spin, color, dest);
    } else {
      node0_printf("%s: File type %d not supported\n",myname, file_type);
    }

#endif    

    break;

  case RELOAD_PARALLEL:

#ifdef HAVE_QIO
    file_type = wpf->file_type;
    if(file_type == FILE_TYPE_W_USQCD_DD_PAIRS){
      status = read_usqcd_wprop_record(wpf, spin, color, src, dest, wqs);
    } else {
      node0_printf("%s: Unsupported file type %d\n", myname, file_type);
      status = 1;
    }
#else

    /* No QIO */
    node0_printf("%s: Parallel reading with this file type not supported\n",
      myname);
    status = 1;

#endif

    break;

  default:
    node0_printf("%s: Unrecognized reload flag.\n", myname);
    terminate(1);
  }
  
  if(timing)
    {
      dtime += dclock();
      if(flag != FRESH && flag != CONTINUE)
	node0_printf("Time to reload wprop spin %d color %d %e\n",
		     spin,color,dtime);
    }

  return status;

} /* reload_wprop_sc_to_field */

/*---------------------------------------------------------------*/
/* save a propagator one source color and spin at a time */
/* recinfo is for USQCD formats */
int 
save_wprop_sc_from_field( int flag, w_prop_file *wpf, 
			  quark_source *wqs,
			  int spin, int color,
			  wilson_vector *src,
			  wilson_vector *prop, 
			  char *recinfo, int timing)
{
  double dtime = 0;
  int status;
  int i; site *s;
  wilson_vector *wv;
  int s0, c0;
#ifdef HAVE_QIO
  int file_type;
#endif
  char myname[] = "save_wprop_sc_from_field";
  
  if(timing)dtime = -dclock();
  status = 0;
  switch(flag){
  case FORGET:
    break;
  case SAVE_ASCII:
    w_ascii_w(wpf,spin,color,prop);
    break;
  case SAVE_SERIAL_SCIDAC:
  case SAVE_PARALLEL_SCIDAC:   
  case SAVE_MULTIFILE_SCIDAC: 
  case SAVE_PARTFILE_SCIDAC:

#ifdef HAVE_QIO
    file_type = wpf->file_type;
    /* Save Dirac source field */
    if(file_type == FILE_TYPE_W_USQCD_DD_PAIRS){

      /* If source is missing, write an empty source (zeros) */
      wilson_vector *src_out = NULL;
      if(src == NULL){
	src_out = create_wv_field();  /* created with zeros */
      } else {
	src_out = src;
      }
      status = write_wpropsource_D_usqcd(wpf->outfile, wqs->descrp, 
					 src_out, wqs->t0);
      if(src_out != NULL)
	destroy_wv_field(src_out);

    } else {

      node0_printf("%s: Unsupported file type %d.\n", myname, file_type);
      status = 1;

    }
    /* Save solution field */
    if(status == 0)
      status = write_prop_usqcd_sc(wpf->outfile, prop, spin, color, recinfo);
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
	node0_printf("Time to save prop spin %d color %d = %e\n",
		     spin,color,dtime);
    }
  return status;
} /* save_wprop_sc_from_field */

/*---------------------------------------------------------------*/
/* Reload a spin_wilson_vector field for a single source color and 4
   source spins in USQCD format, or keep current propagator: FRESH,
   CONTINUE, RELOAD_ASCII, RELOAD_SERIAL, RELOAD_PARALLEL

   Return value is 

  -1 end of file
   0 normal exit value
   1 read error
*/

int 
reload_wprop_c_to_field( int flag, w_prop_file *wpf, 
			 quark_source *wqs, int spin, int color,
			 spin_wilson_vector *source,
			 spin_wilson_vector *dest, int timing)
{

  double dtime = 0;
  int status;
  wilson_vector *psi, *src;

  if(flag == FRESH)return 0;

  status = 0;
  if(timing)dtime = -dclock();

  if(wqs->type == DIRAC_PROPAGATOR_FILE)
    src = create_wv_field();
  else
    src = NULL;
  psi = create_wv_field();
    
  for(int spin=0;spin<4;spin++){
    reload_wprop_sc_to_field(flag, wpf, wqs, spin, color, src, psi, 0);
    if(wqs->type == DIRAC_PROPAGATOR_FILE)
      insert_swv_from_wv(source, spin, src);
    insert_swv_from_wv(dest, spin, psi);
  }
    
  destroy_wv_field(psi);
  if(wqs->type == DIRAC_PROPAGATOR_FILE)
    destroy_wv_field(src);

  r_close_wprop(flag, wpf);

  if(timing)
    {
      dtime += dclock();
      if(flag != FRESH)
	node0_printf("Time to reload wprop %e\n",dtime);
    }
  
  return status;

} /* reload_wprop_c_to_field */

/*---------------------------------------------------------------*/
/* Reload a full propagator (3 source colors and 4 source spins) in
   most of the formats, or fresh propagator, or keep current propagator.
   Destination field is a wilson_prop_field.

   FRESH, CONTINUE,
   RELOAD_ASCII, RELOAD_SERIAL, RELOAD_PARALLEL, RELOAD_MULTIDUMP

   Return value is 

  -1 end of file
   0 normal exit value
   1 read error

*/
int 
reload_wprop_to_wp_field( int flag, char *filename, quark_source *wqs,
			  wilson_prop_field *source, wilson_prop_field *dest, int timing)
{

  double dtime = 0;
  int status;
  int spin, color;
  w_prop_file *wpf;
  wilson_vector *psi, *src;
  
  if(flag == FRESH)return 0;

  status = 0;
  if(timing)dtime = -dclock();

  wpf = r_open_wprop(flag, filename);
  if(wpf == NULL)return 1;
  
  if(wqs->type == DIRAC_PROPAGATOR_FILE)
    src = create_wv_field();
  else
    src = NULL;
  psi = create_wv_field();
    
  /* Loop over source colors */
  status = 0;
  for(color=0;color<dest->nc;color++){
    /* Loop over source spins */
    for(spin=0;spin<4;spin++){
      reload_wprop_sc_to_field(flag, wpf, wqs, spin, color, src, psi, 0);
      if(wqs->type == DIRAC_PROPAGATOR_FILE)
	copy_wp_from_wv(source, src, color, spin);
      copy_wp_from_wv(dest, psi, color, spin);
    }
  }
    
  r_close_wprop(flag, wpf);
    
  destroy_wv_field(psi);
  if(wqs->type == DIRAC_PROPAGATOR_FILE)
    destroy_wv_field(src);
    
  if(timing)
    {
      dtime += dclock();
      if(flag != FRESH)
	node0_printf("Time to reload wprop %e\n",dtime);
    }
  
  return status;
  
} /* reload_wprop_to_wp_field */

/*---------------------------------------------------------------*/
/* save the full propagator (src is wilson_prop_field type)
   FORGET,
   SAVE_ASCII, 
   SAVE_SERIAL_SCIDAC, SAVE_PARALLEL_SCIDAC, SAVE_PARTFILE_SCIDAC,
   SAVE_MULTFILE_SCIDAC
*/
int 
save_wprop_from_wp_field( int flag, char *filename, char *recxml,
			  quark_source *wqs, wilson_prop_field *source,
			  wilson_prop_field *prop,  int timing)
{
  w_prop_file *wpf;
  int spin, color, status = 0;
  double dtime = 0;
  wilson_vector *wvprop, *wvsrc;
  
  if(timing)dtime = -dclock();

  if(flag == FORGET)return status;
  wpf = w_open_wprop(flag, filename, DIRAC_FIELD_FILE);
  if(wpf == NULL)return 1;
  
  wvprop = create_wv_field();
  wvsrc = create_wv_field();

  for(color = 0; color < prop->nc; color++)
    for(spin = 0; spin < 4; spin++)
      {
	copy_wv_from_wp(wvsrc, prop, color, spin);
	copy_wv_from_wp(wvprop, prop, color, spin);
	status = save_wprop_sc_from_field (flag, wpf, wqs, spin, color, 
					   wvsrc, wvprop, recxml, 0);
	if(status != 0)break;
      }

  destroy_wv_field(wvsrc);
  destroy_wv_field(wvprop);

  if(timing)
    {
      dtime += dclock();
      if(flag != FORGET)
	node0_printf("Time to save wprop %e %s\n",dtime, filename);
    }
  w_close_wprop(flag, wpf);
  
  return status;

} /* save_wprop_from_field */

/*---------------------------------------------------------------*/
/* Temporary procedure to support legacy applications that read to the
   site structure */
/* reload a propagator for a single source color and spin

   Return value is 

  -1 end of file
   0 normal exit value
   1 read error
*/

int 
reload_wprop_sc_to_site( int flag, w_prop_file *wpf,
			 quark_source *wqs, int spin, int color, 
			 field_offset dest, int timing)
{
  int i,status = 0;
  site *s;
  wilson_vector *wv;

  if(flag == CONTINUE)return 0;
  if(flag == FRESH){
    FORALLSITES(i,s)clear_wvec( (wilson_vector *)F_PT(s,dest) );
    return 0;
  }

  wv = (wilson_vector *)malloc(sites_on_node*sizeof(wilson_vector));
  if(wv == NULL){
    printf("reload_wprop_sc_to_site(%d): Can't allocate wv\n", this_node);
    terminate(1);
  }

  status = reload_wprop_sc_to_field( flag, wpf, wqs, spin, color, NULL, wv, timing);
  if(status)return status;

  FORALLSITES(i,s){
    copy_wvec(wv+i, (wilson_vector *)F_PT(s,dest));
  }

  free(wv);
  return status;

} /* reload_wprop_sc_to_site */

/*---------------------------------------------------------------*/
/* Temporary procedure to support legacy applications that read to the
   site structure */
/* save a propagator one source color and spin at a time */
 
int 
save_wprop_sc_from_site( int flag, w_prop_file *wpf, quark_source *wqs, 
			 int spin, int color, field_offset src, int timing)
{
  int status;
  int i; site *s;
  wilson_vector *wv;
  char recinfo[] = "Dummy record info";

  wv = (wilson_vector *)malloc(sites_on_node*sizeof(wilson_vector));
  if(wv == NULL){
    printf("save_wprop_sc_from_site(%d): Can't allocate wv\n", this_node);
    terminate(1);
  }

  FORALLSITES(i,s){
    copy_wvec((wilson_vector *)F_PT(s,src),wv+i);
  }

  status = save_wprop_sc_from_field( flag, wpf, wqs, spin, color, NULL, wv, 
				     recinfo, timing);

  free(wv);
  return status;

} /* save_wprop_sc_from_site */

/*---------------------------------------------------------------*/
/* DEPRECATED */
/* Temporary procedure to support legacy applications that read to the
   site structure */
/* reload a full propagator (all source colors and spins) in any of
   the formats, or cold propagator, or keep current propagator:

   FRESH, CONTINUE,
   RELOAD_ASCII, RELOAD_SERIAL, RELOAD_PARALLEL, RELOAD_MULTIDUMP

   Return value is 

  -1 end of file
   0 normal exit value
   1 read error

*/
int 
reload_wprop_to_site( int flag, char *filename, quark_source *wqs,
		      field_offset dest, int timing )
{
  int i,status;
  site *s;
  wilson_propagator *wp;
  int spin, color;
  field_offset destcs;
  w_prop_file *wpf = NULL;
  
  status = 0;

  switch(flag){

  case CONTINUE:  /* do nothing */
    break;

  case FRESH:     /* zero initial guess */
    FORALLSITES(i,s){
      wp = (wilson_propagator *)F_PT(s,dest);
      for(color = 0; color < 3; color++)for(spin = 0; spin < 4; spin++)
	clear_wvec( &wp->c[color].d[spin] );
    }
    break;
    
  case RELOAD_ASCII:
    wpf = r_ascii_w_i(filename);
    for(color = 0; color < 3; color++)for(spin = 0; spin < 4; spin++)
      {
	destcs = dest + color*sizeof(spin_wilson_vector) + 
	  spin*sizeof(wilson_vector);
	if( r_ascii_w(wpf,spin,color,destcs) != 0)status = 1;
      }
    r_ascii_w_f(wpf);
    if(status == 0)
      node0_printf("Read wprop in ASCII format from file %s\n",wpf->filename);
    break;

  case RELOAD_SERIAL:
  case RELOAD_PARALLEL:
    wpf = r_open_wprop(flag, filename);
    for(color = 0; color < 3; color++)for(spin = 0; spin < 4; spin++)
      {
	destcs = dest + color*sizeof(spin_wilson_vector) + 
	  spin*sizeof(wilson_vector);
	if( reload_wprop_sc_to_site(flag, wpf, wqs, spin, color, destcs, timing)
	    != 0) status = 1;
      }
    r_close_wprop(flag, wpf);
    break;
    
  default:
    node0_printf("Bad reload flag %d\n",flag);
    terminate(1);
  }

  return status;

} /* reload_wprop_to_site */

/*---------------------------------------------------------------*/
/* DEPRECATED */
/* Temporary procedure to support legacy applications that read to the
   site structure */
/* save the full propagator
   FORGET,
   SAVE_ASCII, SAVE_SERIAL, SAVE_PARALLEL, SAVE_MULTIDUMP, SAVE_CHECKPOINT
*/
int 
save_wprop_from_site( int flag, char *filename, quark_source *wqs,
		      field_offset src, char *recxml, int timing)
{
  int spin, color;
  int source_type;
  w_prop_file *wpf;
  wilson_vector *srccs;
  site *s;
  int i;
  int status;

  status = 0;

  srccs = create_wv_field();

  switch(flag){

  case FORGET:
    break;

  case SAVE_ASCII:
    wpf = w_ascii_w_i(filename);
    for(color = 0; color < 3; color++)for(spin = 0; spin < 4; spin++)
      {
	FORALLSITES(i,s){
	  copy_wvec(&((wilson_propagator *)F_PT(s,src))->c[color].d[spin],
		    &srccs[i]);
	  w_ascii_w(wpf,spin,color,srccs);
	}
      }
    w_ascii_w_f(wpf);
    if(status == 0)
      node0_printf("Saved wprop in ASCII format to file %s\n",filename);
    break;
    
  case SAVE_SERIAL:
  case SAVE_PARALLEL:
    source_type = POINT;  /* For legacy code, this is the most likely choice */
    wpf = w_open_wprop(flag, filename, source_type);
    for(color = 0; color < 3; color++)for(spin = 0; spin < 4; spin++)
      {
	FORALLSITES(i,s){
	  copy_wvec(&((wilson_propagator *)F_PT(s,src))->c[color].d[spin],
		    &srccs[i]);
	  save_wprop_sc_from_field(flag, wpf, wqs, spin, color, NULL, srccs, 
				   recxml, timing);
	}
      }
    w_close_wprop(flag, wpf);
    if(status == 0)
      node0_printf("Saved wprop serially to file %s\n",filename);
    break;
  default:
    node0_printf("save_wprop_from_site: Unrecognized save flag.\n");
    terminate(1);
  }

  destroy_wv_field(srccs);
  
  return status;

} /* save_wprop_from_site */

/*---------------------------------------------------------------*/
/* Translate output flag to the appropriate input flag for restoring
   a propagator that was temporarily written to disk  */
int
convert_outflag_to_inflag_wprop(int outflag){
  switch(outflag){
  case SAVE_ASCII:
    return RELOAD_ASCII;
  case SAVE_SERIAL_FM:
  case SAVE_SERIAL_FM_SC:                
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
/* find out what kind of starting propagator to use, 
   and propagator name if necessary.  This routine is only 
   called by node 0.
   */
int 
ask_starting_wprop( FILE *fp, int prompt, int *flag, char *filename )
{
  const char *savebuf;
  int status;
  char myname[] = "ask_starting_wprop";
  
  if (prompt==1) {
    printf("loading wilson propagator:\n enter 'fresh_wprop', ");
    printf("'continue_wprop', 'reload_ascii_wprop', ");
    printf("'reload_serial_wprop', 'reload_parallel_wprop', ");
    printf("\n");
  }

  savebuf = get_next_tag(fp, "read wprop command", myname);
  if (savebuf == NULL)return 1;

  printf("%s ",savebuf);
  if(strcmp("fresh_wprop",savebuf) == 0 ){
    *flag = FRESH;
    printf("\n");
  }
  else if(strcmp("continue_wprop",savebuf) == 0 ) {
    *flag = CONTINUE;
    printf("(if possible)\n");
  }
  else if(strcmp("reload_ascii_wprop",savebuf) == 0 ) {
    *flag = RELOAD_ASCII;
  }
  else if(strcmp("reload_serial_wprop",savebuf) == 0 ) {
    *flag = RELOAD_SERIAL;
  }
  else if(strcmp("reload_parallel_wprop",savebuf) == 0 ) {
    *flag = RELOAD_PARALLEL;
  }
  else{
    printf("ERROR IN INPUT: starting propagator command %s is invalid\n",
	   savebuf); 
    return 1;
  }
  
  /*read name of file and load it */
  if( *flag != FRESH && *flag != CONTINUE ){
    if(prompt==1)printf("enter name of file containing props\n");
    status=scanf("%s",filename);
    if(status !=1) {
      printf("\n%s(%d): ERROR IN INPUT: Can't read file name.\n",
	     myname, this_node);
      return 1;
    }
    printf("%s\n",filename);
  }
  return 0;
} /* ask_starting_wprop */

/*---------------------------------------------------------------*/
/* find out what do to with propagator at end, and propagator name if
   necessary.  This routine is only called by node 0.
   */

static void print_output_wprop_choices(void){
    printf("'forget_wprop', 'save_ascii_wprop', ");
    printf("'save_serial_fm_wprop', 'save_serial_fm_sc_wprop', ");
    printf("'save_serial_scidac_wprop', 'save_parallel_scidac_wprop', ");
    printf("'save_partfile_scidac_wprop', 'save_multifile_scidac_wprop', ");
    printf("\n");
}

static int
parse_output_wprop_choices(int *flag, const char *savebuf){

  if(strcmp("save_ascii_wprop",savebuf) == 0 )  {
    *flag=SAVE_ASCII;
  }
  else if(strcmp("save_serial_fm_wprop",savebuf) == 0 ) {
    *flag=SAVE_SERIAL_FM;
  }
  else if(strcmp("save_serial_fm_sc_wprop",savebuf) == 0 ) {
    *flag=SAVE_SERIAL_FM_SC;
  }
  else if(strcmp("save_serial_scidac_wprop",savebuf) == 0 ) {
    *flag=SAVE_SERIAL_SCIDAC;
  }
  else if(strcmp("save_parallel_scidac_wprop",savebuf) == 0 ) {
    *flag=SAVE_PARALLEL_SCIDAC;
  }
  else if(strcmp("save_partfile_scidac_wprop",savebuf) == 0 ) {
    *flag=SAVE_PARTFILE_SCIDAC;
  }
  else if(strcmp("save_multifile_scidac_wprop",savebuf) == 0 ) {
    *flag=SAVE_MULTIFILE_SCIDAC;
  }
  else if(strcmp("forget_wprop",savebuf) == 0 ) {
    *flag=FORGET;
    printf("\n");
  }
  else
    return 1;

  return 0;
}

int 
ask_ending_wprop( FILE *fp, int prompt, int *flag, char *filename ){
  const char *savebuf;
  int status;
  char myname[] = "ask_ending_wprop";
  
  if (prompt==1) {
    printf("save wilson propagator:\n enter ");
    print_output_wprop_choices();
  }

  savebuf = get_next_tag(fp, "write wprop command", myname);
  if (savebuf == NULL)return 1;

  printf("%s ",savebuf);

  if(parse_output_wprop_choices(flag, savebuf) != 0){
    printf("ERROR IN INPUT: ending propagator command %s is invalid\n",
	   savebuf);
    printf("Choices are \n");
    print_output_wprop_choices();
    return 1;
  }
  
  if( *flag != FORGET ){
    if(prompt==1)printf("enter filename\n");
    status=scanf("%s",filename);
    if(status !=1){
      printf("\n%s(%d): ERROR IN INPUT. Can't read filename\n",
	     myname, this_node); 
      return 1;
    }
    printf("%s\n",filename);
  }
  return 0;
} /* ask_ending_wprop */

/* In this case we allow writing the Dirac propagator in propagator format or
   in source format suitable for use as an extended source on the whole lattice. */
int
ask_ending_wprop_or_wsource(FILE *fp, int prompt, int *flag, int *type, 
			    int *t0, char *descrp, char *filename){
  const char *savebuf;
  int status;
  char myname[] = "ask_ending_wprop_or_wsource";
  
  if (prompt==1) {
    printf("save wilson propagator:\n enter ");
    print_output_wprop_choices();
    print_output_quark_source_choices();
  }

  savebuf = get_next_tag(fp, "write wprop command", myname);
  if (savebuf == NULL)return 1;

  printf("%s ",savebuf);

  /* First see if a regular propagator file is requested */
  if(parse_output_wprop_choices(flag, savebuf) == 0){
    *type = DIRAC_PROPAGATOR_FILE;
  } else {
    /* If not, see if a source format file is requested */
    if(parse_output_quark_source_choices(flag, type, descrp, savebuf) != 0){
      printf("\n%s: ERROR IN INPUT: propagator or source file command %s is invalid\n",
	     myname, savebuf); 
      printf("Choices are \n");
      print_output_wprop_choices();
      print_output_quark_source_choices();
      return 1;
    }
  }

  if( *type == VECTOR_FIELD_FILE ){
    printf("\n%s: Can't convert a Dirac propagator to a color vector source\n", myname);
    return 1;
  }

  if( *flag != FORGET ){
    if(prompt==1)printf("enter filename\n");
    status=scanf("%s",filename);
    if(status !=1){
      printf("\n%s(%d): ERROR IN INPUT. Can't read filename\n",
	     myname, this_node); 
      return 1;
    }
    printf("%s\n",filename);
  }

  /* Read t0 line if we are looking for a t0 */
  if(t0 != NULL && *flag != FORGET){
    char t0_string[5];
    int status = 0;
    status += get_s(stdin, prompt, "t0", t0_string);
    /* Provide for writing all time slices */
    if(status == 0){
      if(strcmp(t0_string, "all") == 0)
	*t0 = ALL_T_SLICES;
      else
	*t0 = atoi(t0_string);
    }
  }

  return 0;
} /* ask_ending_wprop_or_wsource */
