/********************** io_helpers_w.c **********************************/
/* MIMD version 7 */
/* CD 10/97 from io_helpers.c DT 8/97 
   General purpose high level propagator I/O routines, 
   to be used by any application that wants them.
   */

#include "generic_wilson_includes.h"
#include "../include/io_lat.h"
#include "../include/io_wprop.h"
#include "../include/io_ksprop.h"
#include "../include/file_types.h"
#include <string.h>
#ifdef HAVE_QIO
#include "../include/io_scidac.h"
#include "../include/io_scidac_w.h"
#include <qio.h>
#endif

static file_type w_prop_list[N_WPROP_TYPES + N_KSPROP_TYPES] =
  { {FILE_TYPE_W_PROP,       W_PROP_VERSION_NUMBER},
    {FILE_TYPE_W_PROP_1996,  W_PROP_VERSION_NUMBER_1996},
    {FILE_TYPE_FM,           IO_UNI_MAGIC},
    {FILE_TYPE_LIME,         LIME_MAGIC_NO},
    {FILE_TYPE_KSPROP,       KSPROP_VERSION_NUMBER},
  };

#ifdef HAVE_QIO
/*---------------------------------------------------------------*/
/* Set up a USQCD SciDAC propagator file for reading */

static w_prop_file *setup_input_usqcd_prop_file(char *filename, int file_type){

  w_prop_file *wpf;

  /* Create a wpf structure. (No file movement here.) */
  wpf = setup_input_w_prop_file(filename);
  /* Open the file and read the header */
  wpf->infile = open_usqcd_prop_read(filename,QIO_SERIAL);
  return wpf;
}

/*---------------------------------------------------------------*/
/* Read a SciDAC propagator source and record according to file type */
static int read_usqcd_prop_record(w_prop_file *wpf, int file_type, 
				  int spin, int color, wilson_vector *dest,
				  wilson_quark_source *wqs ){
  int status;
  int input_spin, input_color;
  char descrp[MAXDESCRP];

  if(file_type == FILE_TYPE_W_USQCD_CD_PAIRS || 
     (file_type == FILE_TYPE_W_USQCD_C1D12 && spin == 0 && color == 0)){
    if(wqs->c_src == NULL){
      wqs->c_src = (complex *)malloc(sizeof(complex)*sites_on_node);
      if(wqs->c_src == NULL){
	printf("read_usqcd_prop_record(%d): No room for source field\n",
	       this_node);
	terminate(1);
      }
      memset(wqs->c_src, 0, sizeof(complex)*sites_on_node);
    }
    read_propsource_C_usqcd(wpf->infile, wqs->descrp, MAXDESCRP, wqs->c_src);
    printf("Read prop source %s\n",descrp);
  }
  else if(file_type == FILE_TYPE_W_USQCD_DD_PAIRS){
    if(wqs->wv_src == NULL){
      wqs->wv_src = 
	(wilson_vector *)malloc(sizeof(wilson_vector)*sites_on_node);
      if(wqs->wv_src == NULL){
	printf("read_usqcd_prop_record(%d): No room for source field\n",
	       this_node);
	terminate(1);
      }
      memset(wqs->wv_src, 0, sizeof(wilson_vector)*sites_on_node);
    }
    read_propsource_D_usqcd(wpf->infile, wqs->descrp, MAXDESCRP, 
			    wqs->wv_src);
    printf("Read prop source %s\n",descrp);
  }

  /* Next, read the solution vector */
  read_proprecord_usqcd(wpf->infile, &input_spin, &input_color, dest);

  /* Check spin and color */
  if(input_spin != spin || input_color != color){
    node0_printf("read_usqcd_prop_record: input spin %d and color %d do not match requested spin %d and color %d\n",
		 input_spin, input_color, spin, color);
    status = 1;
  }
  status = 0;

  return status;
}

/*---------------------------------------------------------------*/
static void interpret_usqcd_save_flag(int *volfmt, int *serpar, 
				 int flag){
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
  case SAVE_PARTITION_SCIDAC:
    *volfmt = QIO_PARTFILE;
    break;
  default:
    *volfmt = QIO_SINGLEFILE;
  }
}
#endif  

/*---------------------------------------------------------------*/
/* Open propagator file for reading one source spin-color at a time */
/* Supported only in MILC formats */

w_prop_file *r_open_wprop(int flag, char *filename)
{
  w_prop_file *wpf = NULL;
  int file_type;
  wilson_propagator *wp;

  switch(flag){
  case RELOAD_ASCII:
    wpf = r_ascii_w_i(filename);
    break;
  case RELOAD_SERIAL:

    /* Sniff out the file type */
    file_type = io_detect(filename, w_prop_list, N_WPROP_TYPES);
    if(file_type < 0){
      printf("Error opening %s\n",filename);
      terminate(1);
    }

    /* For FNAL types we need to look farther to distinguish Wilson
       prop files from KS prop files  */
    if(file_type == FILE_TYPE_FM)
      file_type = io_detect_fm(filename);

    /* For QIO(LIME) types, same thing */
    if(file_type == FILE_TYPE_LIME){
#ifdef HAVE_QIO
      file_type = io_detect_w_usqcd(filename);
#else
      node0_printf("This looks like a QIO file, but to read it requires QIO compilation\n");
      return NULL;
#endif
    }

    /* NOTE: We are not detecting any USQCD KS prop files until we have
       a standard for it */
    
    if(file_type < 0){
      node0_printf("r_open_wprop: Can't interpret file %s\n", filename);
      return NULL;
    }
    /* If it is an FNAL propagator file, read the full prop and cache it */
    else if(file_type == FILE_TYPE_W_FMPROP){
      /* There are two types of Fermilab propagator files.
	 One is sorted by source color and spin, so the site
	 datum is a Dirac vector.  The other has the full propagator
	 on each site. */
      wpf = r_serial_w_fm_i(filename);
      wpf->prop = (wilson_propagator *)
	malloc(sites_on_node*sizeof(wilson_propagator));
      wp = wpf->prop;
      if(wp == NULL){
	printf("Can't malloc for full input propagator\n");
	terminate(1);
      }
      /* For either format we must read the entire propagator so we
	 can convert the spin basis */
      r_serial_w_fm_to_field(wpf, wp);
      
      /* Convert from FNAL to MILC */
      convert_wprop_fnal_to_milc_field(wp);

      /* Indicate propagator data is cached */
      wpf->file_type = file_type = FILE_TYPE_W_STORE;
    }
    /* If it is a KS propagator file, we convert to naive Dirac and cache it */
    else if(file_type == FILE_TYPE_KSPROP || 
	    file_type == FILE_TYPE_KSFMPROP || 
	    file_type == FILE_TYPE_KSQIOPROP ){
      su3_vector *ksp;
      int status;
      int ks_source_r[4] = {0, 0, 0, 0};  /* Assumed for now */

      node0_printf("Converting KS prop %s to naive Dirac.\n", filename);
      node0_printf("Assuming KS source is %d %d %d %d\n",
		   ks_source_r[0],ks_source_r[1],
		   ks_source_r[2],ks_source_r[3]);

      /* Create a wpf structure. (No file movement here.) */
      wpf = setup_input_w_prop_file(filename);
      /* Make space for the naive propagator */
      wpf->prop = (wilson_propagator *)
	malloc(sites_on_node*sizeof(wilson_propagator));
      wp = wpf->prop;
      if(wp == NULL){
	printf("No room for naive propagator\n");
	terminate(1);
      }
      /* Make space for the KS propagator */
      ksp = (su3_vector *)malloc(3*sizeof(su3_vector)*sites_on_node);
      if(ksp == NULL){
	printf("No room for input KS propagator\n");
	terminate(1);
      }
      /* Read the KS propagator */
      status = reload_ksprop_to_field(RELOAD_SERIAL, filename, ksp, 1);
      
      /* Convert from KS to naive */
      convert_ksprop_to_wprop(wp, ksp, ks_source_r);
      free(ksp);

      /* Indicate propagator is now cached and not in a file */
      wpf->file_type = file_type = FILE_TYPE_W_STORE;
    }
#ifdef HAVE_QIO
    /* SciDAC propagator format */
    else if(file_type == FILE_TYPE_W_USQCD_C1D12 ||
	    file_type == FILE_TYPE_W_USQCD_DD_PAIRS ||
	    file_type == FILE_TYPE_W_USQCD_CD_PAIRS){
      wpf = setup_input_usqcd_prop_file(filename, file_type);
    }
#endif
    /* All other formats -- just open the file if we can */
    else {
      wpf = r_serial_w_i(filename);
    }
    wpf->file_type = file_type;
    break;
  case RELOAD_PARALLEL:
    /* Sniff out the file type */
    file_type = io_detect(filename, w_prop_list, N_WPROP_TYPES);

    /* For QIO(LIME) types, same thing */
    if(file_type == FILE_TYPE_LIME){
#ifdef HAVE_QIO
      file_type = io_detect_w_usqcd(filename);
#else
      node0_printf("This looks like a QIO file, but to read it requires QIO compilation\n");
      return NULL;
#endif
    }

#ifdef HAVE_QIO
    /* With QIO and MILC formats we can do parallel I/O */
    if(file_type == FILE_TYPE_W_USQCD_C1D12 ||
       file_type == FILE_TYPE_W_USQCD_DD_PAIRS ||
       file_type == FILE_TYPE_W_USQCD_CD_PAIRS){
      wpf = setup_input_usqcd_prop_file(filename, file_type);
    }
#endif
    if(file_type != FILE_TYPE_W_USQCD_C1D12 &&
       file_type != FILE_TYPE_W_USQCD_DD_PAIRS &&
       file_type != FILE_TYPE_W_USQCD_CD_PAIRS){

      wpf = r_parallel_w_i(filename);
      r_parallel_w_c(wpf);
    }
    break;
  case RELOAD_MULTIDUMP:
    wpf = r_multidump_w_i(filename);
    r_multidump_w_c(wpf);
    break;

  default:
    wpf = NULL;
  }

  return wpf;
}

/*---------------------------------------------------------------*/
/* Open propagator file for writing a wilson vector for one
   source spin and color at a time. */

w_prop_file *w_open_wprop(int flag, char *filename, int source_type)
{
  w_prop_file *wpf = NULL;
  wilson_propagator *wp;
#ifdef HAVE_QIO
  char *fileinfo;
  int volfmt, serpar, file_type;
#endif
  
  switch(flag){
  case SAVE_ASCII:
    wpf = w_ascii_w_i(filename);
    break;

  case SAVE_SERIAL:
    wpf = w_serial_w_i(filename);
    break;

  case SAVE_SERIAL_FM:
    wpf = w_serial_w_fm_i(filename);
    /* Allocate space for the entire propagator */
    wpf->prop = (wilson_propagator *)
      malloc(sites_on_node*sizeof(wilson_propagator));
    wp = wpf->prop;
    if(wp == NULL){
      printf("Can't malloc for full output propagator\n");
      terminate(1);
    }
    break;

  case SAVE_SERIAL_FM_SC:
    wpf = w_serial_w_fm_sc_i(filename);
    /* Allocate space for the entire propagator */
    wpf->prop = (wilson_propagator *)
      malloc(sites_on_node*sizeof(wilson_propagator));
    wp = wpf->prop;
    if(wp == NULL){
      printf("Can't malloc for full output propagator\n");
      terminate(1);
    }
    break;

  case SAVE_SERIAL_SCIDAC:
  case SAVE_PARALLEL_SCIDAC:   
  case SAVE_MULTIFILE_SCIDAC: 
  case SAVE_PARTITION_SCIDAC:

#ifdef HAVE_QIO
    wpf = setup_output_w_prop_file();
    interpret_usqcd_save_flag(&volfmt, &serpar, flag);
    file_type = choose_usqcd_file_type(source_type);
    wpf->file_type = file_type;
    fileinfo = create_w_QCDML();
    wpf->outfile = open_usqcd_prop_write(filename, volfmt, serpar, QIO_ILDGNO,
					 NULL, file_type, fileinfo);
    free(fileinfo);
    
#else
    node0_printf("w_open_wprop: SciDAC formats require QIO compilation\n");
    terminate(1);
#endif
    break;

  case SAVE_PARALLEL:
    wpf = w_parallel_w_i(filename);
    w_parallel_w_c(wpf);
    break;

  case SAVE_MULTIDUMP:
    wpf = w_multidump_w_i(filename);
    w_multidump_w_c(wpf);
    break;

  case SAVE_CHECKPOINT:
    wpf = w_checkpoint_w_i(filename);
    w_checkpoint_w_c(wpf);
    break;

  default:
    wpf = NULL;
  }

  return wpf;
}

/*---------------------------------------------------------------*/
/* reload a propagator for a single source color and spin in MILC
   format, or cold propagator, or keep current propagator: FRESH,
   CONTINUE, RELOAD_ASCII, RELOAD_SERIAL, RELOAD_PARALLEL,
   RELOAD_MULTIDUMP
   */
int reload_wprop_sc_to_site( int flag, w_prop_file *wpf,
		       int spin, int color, field_offset dest, int timing)
{
  /* 0 normal exit value
     1 read error */

  double dtime = 0;
  int i,status;
  int file_type;
  site *s;
  wilson_propagator *wp;
  wilson_vector *wv;
  int s0, c0;

  if(timing)dtime = -dclock();
  status = 0;
  switch(flag){
  case CONTINUE:  /* do nothing */
    break;
  case FRESH:     /* zero initial guess */
    FORALLSITES(i,s)clear_wvec( (wilson_vector *)F_PT(s,dest));
    break;
  case RELOAD_ASCII:
    node0_printf("Warning: no source spin and color checked when reading from an ASCII file\n");
    status = r_ascii_w(wpf,spin,color,dest);
    break;
  case RELOAD_SERIAL:
    file_type = wpf->file_type;
    wp = wpf->prop;
    /* Special treatment for a cached propagator */
    if(file_type == FILE_TYPE_W_STORE){
      /* Copy input Wilson vector for this color and spin from buffer */
      FORALLSITES(i,s){
	wv = (wilson_vector *)F_PT(s,dest);
	for(s0=0;s0<4;s0++)for(c0=0;c0<3;c0++)
	  {
	    wv->d[s0].c[c0].real = wp[i].c[color].d[spin].d[s0].c[c0].real;
	    wv->d[s0].c[c0].imag = wp[i].c[color].d[spin].d[s0].c[c0].imag;
	  }
      }
      status = 0;
    }
#ifdef HAVE_QIO
    else if(file_type == FILE_TYPE_W_USQCD_C1D12 ||
	    file_type == FILE_TYPE_W_USQCD_DD_PAIRS ||
	    file_type == FILE_TYPE_W_USQCD_CD_PAIRS){
      node0_printf("reload_wprop_sc_to_site: Reading SciDAC prop to site structure is not supported\n");
      terminate(1);
    }
#endif
    else {
      status = r_serial_w_to_site(wpf,spin,color,dest); 
    }
    break;
  case RELOAD_PARALLEL:
    /* Reopen, read, and close temporarily */
    r_parallel_w_o(wpf);
    status = r_parallel_w_to_site(wpf,spin,color,dest);
    r_parallel_w_c(wpf);
    break;
  case RELOAD_MULTIDUMP:
    /* Reopen, read, and close temporarily */
    r_multidump_w_o(wpf);
    status = r_multidump_w_to_site(wpf,spin,color,dest);
    r_multidump_w_c(wpf);
    break;
  default:
    node0_printf("reload_wprop_sc_to_site: Unrecognized reload flag.\n");
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

} /* reload_wprop_sc_to_site */

/*---------------------------------------------------------------*/
/* reload a propagator for a single source color and spin in MILC
   format, or cold propagator, or keep current propagator: FRESH,
   CONTINUE, RELOAD_ASCII, RELOAD_SERIAL, RELOAD_PARALLEL,
   RELOAD_MULTIDUMP
   */
int reload_wprop_sc_to_field( int flag, w_prop_file *wpf, 
			      wilson_quark_source *wqs, int spin, int color, 
			      wilson_vector *dest, int timing)
{
  /* 0 normal exit value
     1 read error */

  double dtime = 0;
  int i,status;
  site *s;
  wilson_propagator *wp;
  wilson_vector *wv;
  int s0, c0;
  int file_type = -999;

  if(timing)dtime = -dclock();
  status = 0;
  switch(flag){
  case CONTINUE:  /* do nothing */
    break;
  case FRESH:     /* zero initial guess */
    FORALLSITES(i,s)clear_wvec( &(dest[i]) );
    break;
  case RELOAD_ASCII:
    node0_printf("Reloading ASCII to temp wprop not supported\n");
    terminate(1);
    break;
  case RELOAD_SERIAL:
    wp = wpf->prop;
    file_type = wpf->file_type;
    /* Special treatment for a cached propagator */
    if(file_type == FILE_TYPE_W_STORE){
    /* Copy input Wilson vector for this color and spin from buffer */
      FORALLSITES(i,s){
	wv = dest + i;
	for(s0=0;s0<4;s0++)for(c0=0;c0<3;c0++)
	  {
	    wv->d[s0].c[c0].real = wp[i].c[color].d[spin].d[s0].c[c0].real;
	    wv->d[s0].c[c0].imag = wp[i].c[color].d[spin].d[s0].c[c0].imag;
	  }
      }
      status = 0;
    }
#ifdef HAVE_QIO
    else if(file_type == FILE_TYPE_W_USQCD_C1D12 ||
	    file_type == FILE_TYPE_W_USQCD_DD_PAIRS ||
	    file_type == FILE_TYPE_W_USQCD_CD_PAIRS){

      /* Read the propagator record */
      status = read_usqcd_prop_record(wpf, file_type, spin, color, dest, wqs);
    }
#endif
    else {
      status = r_serial_w_to_field(wpf,spin,color,dest); 
    }
    break;
  case RELOAD_PARALLEL:
#ifdef HAVE_QIO
    if(file_type == FILE_TYPE_W_USQCD_C1D12 ||
       file_type == FILE_TYPE_W_USQCD_DD_PAIRS ||
       file_type == FILE_TYPE_W_USQCD_CD_PAIRS){
      status = read_usqcd_prop_record(wpf, file_type, spin, color, dest, wqs);
    }
#endif
    /* If not SciDAC */
    if(file_type != FILE_TYPE_W_USQCD_C1D12 &&
       file_type != FILE_TYPE_W_USQCD_DD_PAIRS &&
       file_type != FILE_TYPE_W_USQCD_CD_PAIRS){
      /* Reopen, read, and close temporarily */
      r_parallel_w_o(wpf);
      status = r_parallel_w_to_field(wpf,spin,color,dest);
      r_parallel_w_c(wpf);
    }
    break;
  case RELOAD_MULTIDUMP:
    /* Reopen, read, and close temporarily */
    r_multidump_w_o(wpf);
    status = r_multidump_w_to_field(wpf,spin,color,dest);
    r_multidump_w_c(wpf);
    break;
  default:
    node0_printf("reload_wprop_sc_to_field: Unrecognized reload flag.\n");
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
/* reload a full propagator (all source colors and spins) in any of
   the formats, or cold propagator, or keep current propagator:

   FRESH, CONTINUE,
   RELOAD_ASCII, RELOAD_SERIAL, RELOAD_PARALLEL, RELOAD_MULTIDUMP
   */
int reload_wprop_to_site( int flag, char *filename,
			  field_offset dest, int timing)
{
  /* 0 normal exit value
     1 read error */
  
  double dtime = 0;
  int i,status;
  site *s;
  wilson_propagator *wp;
  int spin, color;
  field_offset destcs;
  int file_type;
  w_prop_file *wpf = NULL;
  
  if(timing)dtime = -dclock();
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
    /* Sniff out the file type */
    file_type = io_detect(filename, w_prop_list, N_WPROP_TYPES);

    /* For FNAL types we need to look farther to distinguish Wilson
       prop files from KS prop files  */
    if(file_type == FILE_TYPE_FM)
      file_type = io_detect_fm(filename);
    
    /* For QIO(LIME) types, same thing */
    if(file_type == FILE_TYPE_LIME){
      node0_printf("Reading a SciDAC propagator in one call not supported\n");
      return 1;
    }

    if(file_type < 0){
      node0_printf("reload_wprop_to_site: Can't read file %s\n", filename);
      return 1;
    }

    if(file_type == FILE_TYPE_W_PROP ||
       file_type == FILE_TYPE_W_PROP_1996)
      {
	if(flag == RELOAD_SERIAL){
	  node0_printf("Reading serially as a MILC Wilson prop file\n");
	  wpf = r_serial_w_i(filename);

	  /* Old MILC format has one record for each source spin, color */
	  for(color = 0; color < 3; color++)for(spin = 0; spin < 4; spin++)
	    {
	      destcs = dest + color*sizeof(spin_wilson_vector) + 
		spin*sizeof(wilson_vector);
	      if( r_serial_w_to_site(wpf,spin,color,destcs) != 0)status = 1; 
	    }
	  r_serial_w_f(wpf);
	}
	else{
	  node0_printf("Reading in parallel as a MILC Wilson prop file\n");
	  wpf = r_parallel_w_i(filename);

	  /* Old MILC format has one record for each source spin, color */
	  for(color = 0; color < 3; color++)for(spin = 0; spin < 4; spin++)
	    {
	      destcs = dest + color*sizeof(spin_wilson_vector) + 
		spin*sizeof(wilson_vector);
	      if( r_parallel_w_to_site(wpf,spin,color,destcs) != 0)status = 1; 
	    }
	  r_parallel_w_f(wpf);
	  node0_printf("Read wprop in parallel from file %s\n",filename);
	}
      }
    
    else if(file_type == FILE_TYPE_W_FMPROP){
      node0_printf("Reading as a Fermilab Wilson prop file\n");
      /*Read the entire file */
      r_prop_w_fm_to_site( filename, dest );
      convert_wprop_fnal_to_milc_site( dest );
    }
    else{
      node0_printf("Unsupported file type %d for reading to site\n",
		   file_type);
      return 1;
    }
    if(status == 0)
      node0_printf("Read wprop serially from file %s\n",filename);
    break;
    
  case RELOAD_MULTIDUMP:
    /* Reopen, read, and close temporarily */
    r_multidump_w_i(filename);
    /* Old MILC format has one record for each source spin, color */
    for(color = 0; color < 3; color++)for(spin = 0; spin < 4; spin++)
      {
	destcs = dest + color*sizeof(spin_wilson_vector) + 
	  spin*sizeof(wilson_vector);
	if(r_multidump_w_to_site(wpf,spin,color,destcs) != 0)status = 1;
      }
    r_multidump_w_f(wpf);
    if(status == 0)
      node0_printf("Read wprop in multidump format from file %s\n",filename);
    break;
    
  default:
    node0_printf("reload_wprop_to_site: Unrecognized reload flag.\n");
    terminate(1);
  }
  
  if(timing)
    {
      dtime += dclock();
      if(flag != FRESH && flag != CONTINUE)
	node0_printf("Time to reload wprop %e\n",dtime);
    }
  
  return status;

} /* reload_wprop_to_site */

/*---------------------------------------------------------------*/
/* reload a full propagator (all source colors and spins) in any of
   the formats, or cold propagator, or keep current propagator:

   FRESH, CONTINUE,
   RELOAD_ASCII, RELOAD_SERIAL, RELOAD_PARALLEL, RELOAD_MULTIDUMP
   */
int reload_wprop_to_field( int flag, char *filename,
			   wilson_propagator *dest, int timing)
{
  /* 0 normal exit value
     1 read error */
  
  double dtime = 0;
  int i,status;
  site *s;
  wilson_propagator *wp;
  wilson_vector *destcs;
  int spin, color;
  int file_type;
  w_prop_file *wpf;
  
  if(timing)dtime = -dclock();
  status = 0;
  switch(flag){

  case CONTINUE:  /* do nothing */
    break;

  case FRESH:     /* zero initial guess */
    FORALLSITES(i,s){
      wp = dest + i;
      for(color = 0; color < 3; color++)for(spin = 0; spin < 4; spin++)
	clear_wvec( &wp->c[color].d[spin] );
    }
    break;
    
  case RELOAD_ASCII:
    node0_printf("Reloading ASCII to temp wprop not supported\n");
    terminate(1);
    if(status == 0)
      node0_printf("Read wprop in ASCII format from file %s\n",filename);
    break;

  case RELOAD_SERIAL:
  case RELOAD_PARALLEL:
    /* Sniff out the file type */
    file_type = io_detect(filename, w_prop_list, N_WPROP_TYPES);

    /* For FNAL types we need to look farther to distinguish Wilson
       prop files from KS prop files  */
    if(file_type == FILE_TYPE_FM)
      file_type = io_detect_fm(filename);
    
    /* For QIO(LIME) types, same thing */
    if(file_type == FILE_TYPE_LIME){
#ifdef HAVE_QIO
      file_type = io_detect_w_usqcd(filename);
#else
      node0_printf("This looks like a QIO file, but to read it requires QIO compilation\n");
      return 0;
#endif
      node0_printf("Reading a SciDAC propagator in one call not supported\n");
      return 1;
    }
    
    if(file_type < 0){
      node0_printf("reload_wprop_to_field: Can't read file %s\n", filename);
      return 1;
    }

    if(file_type == FILE_TYPE_W_PROP ||
       file_type == FILE_TYPE_W_PROP_1996)
      {
	if(flag == RELOAD_SERIAL){
	  node0_printf("Reading serially as a MILC Wilson prop file\n");

	  /* Old MILC format has one record for each source spin, color */
	  /* So we allocate space for one wilson vector per site */
	  destcs = (wilson_vector *)malloc(sites_on_node*
					   sizeof(wilson_vector));
	  if(destcs == NULL){
	    node0_printf("Can't malloc space to read propagator\n");
	    status = 1;
	  }
	  
	  wpf = r_serial_w_i(filename);
	  for(color = 0; color < 3; color++)for(spin = 0; spin < 4; spin++)
	    {
	      if(status == 0){
		if( r_serial_w_to_field(wpf,spin,color,destcs) != 0)
		  status += 1; 
		FORALLSITES(i,s){
		  dest[i].c[color].d[spin] = destcs[i];
		}
	      }
	    }
	  
	  r_serial_w_f(wpf);
	  free(destcs); destcs = NULL;
	}
	else {
	  node0_printf("Reading in parallel as a MILC Wilson prop file\n");

	  /* Old MILC format has one record for each source spin, color */
	  /* So we allocate space for one wilson vector per site */
	  destcs = (wilson_vector *)malloc(sites_on_node*
					   sizeof(wilson_vector));
	  if(destcs == NULL){
	    node0_printf("Can't malloc space to read propagator\n");
	    status = 1;
	  }
	  
	  wpf = r_parallel_w_i(filename);
	  for(color = 0; color < 3; color++)for(spin = 0; spin < 4; spin++)
	    {
	      if(status == 0){
		if( r_parallel_w_to_field(wpf,spin,color,destcs) != 0)
		  status += 1; 
		FORALLSITES(i,s){
		  dest[i].c[color].d[spin] = destcs[i];
		}
	      }
	    }
	  
	  r_parallel_w_f(wpf);
	  free(destcs); destcs = NULL;
	}
      }
    
    else if(file_type == FILE_TYPE_W_FMPROP){
      node0_printf("Reading as a Fermilab Wilson prop file\n");
      r_prop_w_fm_to_field( filename, dest );
      convert_wprop_fnal_to_milc_field( dest );
    }

    /* If it is a KS propagator file, we convert it to naive Dirac */
    else if(file_type == FILE_TYPE_KSPROP || 
	    file_type == FILE_TYPE_KSFMPROP || 
	    file_type == FILE_TYPE_KSQIOPROP ){
      su3_vector *ksp;
      int status;
      int ks_source_r[4] = {0, 0, 0, 0};  /* Assumed for now */

      node0_printf("Converting KS prop %s to naive Dirac\n",filename);
      /* Make space for the KS propagator */
      ksp = (su3_vector *)malloc(3*sizeof(su3_vector)*sites_on_node);
      if(ksp == NULL){
	printf("No room for input KS propagator\n");
	terminate(1);
      }
      /* Read the KS propagator */
      status = reload_ksprop_to_field(RELOAD_SERIAL, filename, ksp, 1);
      
      /* Convert from KS to naive */
      convert_ksprop_to_wprop(dest, ksp, ks_source_r);
      free(ksp);
    }
    else{
      node0_printf("Unsupported file type %d\n",file_type);
      return 1;
    }
    if(status == 0)
      node0_printf("Read wprop serially from file %s\n",filename);
    break;
    
  case RELOAD_MULTIDUMP:
    /* Old MILC format has one record for each source spin, color */
    /* So we allocate space for one wilson vector per site */
    destcs = (wilson_vector *)malloc(sites_on_node*
				     sizeof(wilson_vector));
    if(destcs == NULL){
      node0_printf("Can't malloc space to read propagator\n");
      status = 1;
    }

    wpf = r_multidump_w_i(filename);
    for(color = 0; color < 3; color++)for(spin = 0; spin < 4; spin++)
      {
	if(status == 0){
	  if(r_multidump_w_to_field(wpf,spin,color,destcs) != 0)
	    status += 1;
	  FORALLSITES(i,s){
	    dest[i].c[color].d[spin] = destcs[i];
	  }
	}
      }
    free(destcs); destcs = NULL;
    r_multidump_w_f(wpf);

    if(status == 0)
      node0_printf("Read wprop in multidump format from file %s\n",filename);
    break;
    
  default:
    node0_printf("reload_wprop_to_field: Unrecognized reload flag.\n");
    terminate(1);
  }
  
  if(timing)
    {
      dtime += dclock();
      if(flag != FRESH && flag != CONTINUE)
	node0_printf("Time to reload wprop %e\n",dtime);
    }
  
  return status;

} /* reload_wprop_to_field */

/*---------------------------------------------------------------*/
/* save a propagator one source color and spin at a time MILC formats only:
   FORGET,
   SAVE_ASCII, SAVE_SERIAL, SAVE_PARALLEL, SAVE_MULTIDUMP, SAVE_CHECKPOINT
*/
int save_wprop_sc_from_site( int flag, w_prop_file *wpf, 
		      int spin, int color, field_offset src, int timing)
{
  double dtime = 0;
  int status;
  int i; site *s;
  wilson_propagator *wp;
  wilson_vector *wv;
  int s0, c0;
  
  if(timing)dtime = -dclock();
  status = 0;
  switch(flag){
  case FORGET:
    break;
  case SAVE_ASCII:
    node0_printf("Warning: no source spin and color recorded in an ASCII file\n");
    w_ascii_w(wpf,spin,color,src);
    break;
  case SAVE_SERIAL:
    w_serial_w_from_site(wpf,spin,color,src);
    break;
  case SAVE_SERIAL_FM:
  case SAVE_SERIAL_FM_SC:
    wp = wpf->prop;
    if(wp == NULL){
      printf("save_wprop_sc_from_site: Propagator field not allocated\n");
      terminate(1);
    }
    /* Add output Wilson vector to propagator buffer */
    FORALLSITES(i,s){
      wv = (wilson_vector *)F_PT(s,src);
      for(s0=0;s0<4;s0++)for(c0=0;c0<3;c0++)
	{
	  wp[i].c[color].d[spin].d[s0].c[c0].real = wv->d[s0].c[c0].real;
	  wp[i].c[color].d[spin].d[s0].c[c0].imag = wv->d[s0].c[c0].imag;
	}
    }
    break;
  case SAVE_SERIAL_SCIDAC:
  case SAVE_PARALLEL_SCIDAC:   
  case SAVE_MULTIFILE_SCIDAC: 
  case SAVE_PARTITION_SCIDAC:

    node0_printf("w_open_wprop: Saving SciDAC formats from site structure data is not supported\n");
    terminate(1);
    break;

  case SAVE_PARALLEL:
    w_parallel_w_o(wpf);
    w_parallel_w_from_site(wpf,spin,color,src);
    w_parallel_w_c(wpf);
    break;
  case SAVE_CHECKPOINT:
    w_checkpoint_w_o(wpf);
    w_checkpoint_w_from_site(wpf,spin,color,src);
    w_checkpoint_w_c(wpf);
    break;
  case SAVE_MULTIDUMP:
    w_multidump_w_o(wpf);
    w_multidump_w_from_site(wpf,spin,color,src);
    w_multidump_w_c(wpf);
    break;
  default:
    node0_printf("save_wprop_sc_from_site: Unrecognized save flag.\n");
    terminate(1);
  }
  
  if(timing)
    {
      dtime += dclock();
      if(flag != FORGET)
	node0_printf("Time to save prop spin %d color %d = %e\n",
		     spin,color,dtime);
    }

  return status;
} /* save_wprop_sc_from_site */

/*---------------------------------------------------------------*/
/* save a propagator one source color and spin at a time */
/* recinfo is for USQCD formats */
int save_wprop_sc_from_field( int flag, w_prop_file *wpf, 
			      wilson_quark_source *wqs,
			      int spin, int color, wilson_vector *src, 
			      char *recinfo, int timing)
{
  double dtime = 0;
  int status;
  int i; site *s;
  wilson_propagator *wp;
  wilson_vector *wv;
  int s0, c0;
#ifdef HAVE_QIO
  int volfmt, serpar;
  int file_type;
#endif
  
  if(timing)dtime = -dclock();
  status = 0;
  switch(flag){
  case FORGET:
    break;
  case SAVE_ASCII:
    node0_printf("Reading to field from ASCII is not supported\n");
    terminate(1);
    break;
  case SAVE_SERIAL:
    w_serial_w_from_field(wpf,spin,color,src);
    break;
  case SAVE_SERIAL_FM:
  case SAVE_SERIAL_FM_SC:
    wp = wpf->prop;
    if(wp == NULL){
      printf("save_wprop_sc_from_field: Propagator field not allocated\n");
      terminate(1);
    }
    /* Add output Wilson vector to propagator buffer */
    FORALLSITES(i,s){
      wv = src + i;
      for(s0=0;s0<4;s0++)for(c0=0;c0<3;c0++)
	{
	  wp[i].c[color].d[spin].d[s0].c[c0].real = wv->d[s0].c[c0].real;
	  wp[i].c[color].d[spin].d[s0].c[c0].imag = wv->d[s0].c[c0].imag;
	}
    }
    break;
  case SAVE_SERIAL_SCIDAC:
  case SAVE_PARALLEL_SCIDAC:   
  case SAVE_MULTIFILE_SCIDAC: 
  case SAVE_PARTITION_SCIDAC:

#ifdef HAVE_QIO
    file_type = wpf->file_type;
    interpret_usqcd_save_flag(&volfmt, &serpar, flag);
    /* Save color source field */
    if(file_type == FILE_TYPE_W_USQCD_CD_PAIRS ||
       (file_type == FILE_TYPE_W_USQCD_C1D12 && spin == 0 && color == 0))
      {
	if(wqs->c_src == NULL){
	  printf("save_wprop_sc_from_field(%d): source field must be defined before calling.\n",this_node);
	  terminate(1);
	}
	status = write_propsource_C_usqcd(wpf->outfile, wqs->descrp, wqs->c_src, wqs->t0);
      }
    /* Save Dirac source field */
    else if(file_type == FILE_TYPE_W_USQCD_DD_PAIRS){
	if(wqs->wv_src == NULL){
	  printf("save_wprop_sc_from_field(%d): source field must be defined before calling.\n",this_node);
	  terminate(1);
	}
      status = write_propsource_D_usqcd(wpf->outfile, wqs->descrp, 
				    wqs->wv_src, wqs->t0);
      if(status != 0)break;
    }
    /* Save solution field */
    if(status == 0)
      status = write_prop_usqcd_sc(wpf->outfile, src, spin, color, recinfo);
#else
    node0_printf("w_open_wprop: SciDAC formats require QIO compilation\n");
    terminate(1);
#endif
    
    break;
  case SAVE_PARALLEL:
    w_parallel_w_o(wpf);
    w_parallel_w_from_field(wpf,spin,color,src);
    w_parallel_w_c(wpf);
    break;
  case SAVE_CHECKPOINT:
    w_checkpoint_w_o(wpf);
    w_checkpoint_w_from_field(wpf,spin,color,src);
    w_checkpoint_w_c(wpf);
    break;
  case SAVE_MULTIDUMP:
    w_multidump_w_o(wpf);
    w_multidump_w_from_field(wpf,spin,color,src);
    w_multidump_w_c(wpf);
    break;
  default:
    node0_printf("save_wprop_sc_from_field: Unrecognized save flag.\n");
    terminate(1);
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
/* save the full propagator
   FORGET,
   SAVE_ASCII, SAVE_SERIAL, SAVE_PARALLEL, SAVE_MULTIDUMP, SAVE_CHECKPOINT
*/
int save_wprop_from_site( int flag, char *filename, char *recxml,
			    field_offset src, int timing)
{
  double dtime = 0;
  int spin, color;
  w_prop_file *wpf;
  wilson_propagator *wp;
  field_offset srccs;
  int i;
  site *s;
  int status;

  if(timing)dtime = -dclock();
  status = 0;
  switch(flag){
  case FORGET:
    break;

  case SAVE_ASCII:
    wpf = w_ascii_w_i(filename);
    for(color = 0; color < 3; color++)for(spin = 0; spin < 4; spin++)
      {
	srccs = src + color*sizeof(spin_wilson_vector) + 
	  spin*sizeof(wilson_vector);
	w_ascii_w(wpf,spin,color,srccs);
      }
    w_ascii_w_f(wpf);
    if(status == 0)
      node0_printf("Saved wprop in ASCII format to file %s\n",filename);
    break;

  case SAVE_SERIAL:
    wpf = w_serial_w_i(filename);
    for(color = 0; color < 3; color++)for(spin = 0; spin < 4; spin++)
      {
	srccs = src + color*sizeof(spin_wilson_vector) + 
	  spin*sizeof(wilson_vector);
	w_serial_w_from_site(wpf,spin,color,srccs);
      }
    w_serial_w_f(wpf);
    if(status == 0)
      node0_printf("Saved wprop serially to file %s\n",filename);
    break;

  case SAVE_SERIAL_FM:
  case SAVE_SERIAL_FM_SC:
    if(flag == SAVE_SERIAL_FM)
      wpf = w_serial_w_fm_i(filename);
    else
      wpf = w_serial_w_fm_sc_i(filename);

    wpf->prop = (wilson_propagator *)malloc(sites_on_node*
					    sizeof(wilson_propagator));
    wp = wpf->prop;
    if(wp == NULL){
      node0_printf("Can't malloc space to save propagator\n");
      status = 1;
      break;
    }
    
    /* Copy full propagator to buffer */
    FORALLSITES(i,s){
      wp[i] = *((wilson_propagator *)F_PT(s,src));
    }

    /* Convert to FNAL spin conventions */
    convert_wprop_milc_to_fnal_field(wp);

    /* Write the converted propagator */
    
    w_serial_w_fm_from_field(wpf, wp);
    w_serial_w_fm_f(wpf);
    break;
  case SAVE_SERIAL_SCIDAC:
  case SAVE_PARALLEL_SCIDAC:   
  case SAVE_MULTIFILE_SCIDAC: 
  case SAVE_PARTITION_SCIDAC:

    node0_printf("Currently saving the entire propagator in SciDAC format is not supported\n");
    status = 1;

    break;
  case SAVE_PARALLEL:
    wpf = w_parallel_w_i(filename);
    for(color = 0; color < 3; color++)for(spin = 0; spin < 4; spin++)
      {
	srccs = src + color*sizeof(spin_wilson_vector) + 
	  spin*sizeof(wilson_vector);
	w_parallel_w_from_site(wpf,spin,color,srccs);
      }
    w_parallel_w_f(wpf);
    if(status == 0)
      node0_printf("Saved wprop in parallel to file %s\n",filename);
    break;
  case SAVE_CHECKPOINT:
    wpf = w_checkpoint_w_i(filename);
    for(color = 0; color < 3; color++)for(spin = 0; spin < 4; spin++)
      {
	srccs = src + color*sizeof(spin_wilson_vector) + 
	  spin*sizeof(wilson_vector);
	w_checkpoint_w_from_site(wpf,spin,color,srccs);
      }
    w_checkpoint_w_f(wpf);
    if(status == 0)
      node0_printf("Saved wprop in checkpoint format to file %s\n",filename);
    break;
  case SAVE_MULTIDUMP:
    wpf = w_multidump_w_i(filename);
    for(color = 0; color < 3; color++)for(spin = 0; spin < 4; spin++)
      {
	srccs = src + color*sizeof(spin_wilson_vector) + 
	  spin*sizeof(wilson_vector);
	w_multidump_w_from_site(wpf,spin,color,srccs);
      }
    w_multidump_w_f(wpf);
    if(status == 0)
      node0_printf("Saved wprop in multidump format to file %s\n",filename);
    break;
  default:
    node0_printf("save_wprop_from_site: Unrecognized save flag.\n");
    terminate(1);
  }
  
  if(timing)
    {
      dtime += dclock();
      if(flag != FORGET)
	node0_printf("Time to save wprop = %e\n",dtime);
    }
  return status;
} /* save_wsprop_from_site */
/*---------------------------------------------------------------*/
/* read the lattice dimensions from a binary Wilson prop file */

int read_lat_dim_wprop(char *filename, int file_type, int *ndim, int dims[])
{
  w_prop_file *wpf;
  int i;

  switch(file_type){
  case FILE_TYPE_W_PROP:
  case FILE_TYPE_W_PROP_1996:
    *ndim = 4;
    nx = -1; ny = -1; nz = -1; nt = -1;
    wpf = r_serial_w_i(filename);
    for(i = 0; i < *ndim; i++)
      dims[i] = wpf->header->dims[i];
    r_serial_w_f(wpf);

    break;
  case FILE_TYPE_W_FMPROP:
    *ndim = 4;
    nx = -1; ny = -1; nz = -1; nt = -1;
    wpf = r_serial_w_fm_i(filename);
    for(i = 0; i < *ndim; i++)
      dims[i] = wpf->header->dims[i];
    r_serial_w_fm_f(wpf);
    break;

  case FILE_TYPE_W_USQCD_C1D12:
  case FILE_TYPE_W_USQCD_DD_PAIRS:
  case FILE_TYPE_W_USQCD_CD_PAIRS:
  case FILE_TYPE_W_USQCD_LHPC:
#ifdef HAVE_QIO
    read_lat_dim_scidac(filename, ndim, dims);
#else
    node0_printf("This looks like a QIO file, but to read it requires QIO compilation\n");
    return 1;
#endif
    break;
  default:
    node0_printf("Unsupported file type %d\n",file_type);
    return 1;
  }
  return 0;
}

/*---------------------------------------------------------------*/
/* save the full propagator
   FORGET,
   SAVE_ASCII, SAVE_SERIAL, SAVE_PARALLEL, SAVE_MULTIDUMP, SAVE_CHECKPOINT
   SAVE_SERIAL_SCIDAC, SAVE_PARALLEL_SCIDAC, SAVE_PARTITION_SCIDAC,
   SAVE_MULTFILE_SCIDAC
*/
int save_wprop_from_field( int flag, char *filename, char *recxml,
			     wilson_propagator *src, int timing)
{
  double dtime = 0;
  int spin, color;
  w_prop_file *wpf = NULL;
  wilson_vector *srccs;
  wilson_propagator *wp;
  int status;
  site *s;
  int i;
  
  if(timing)dtime = -dclock();
  status = 0;
  switch(flag){
  case FORGET:
    break;

  case SAVE_ASCII:
    node0_printf("Writing from field to ASCII not supported\n");
    terminate(1);
    break;

  case SAVE_SERIAL:
    srccs = (wilson_vector *)malloc(sites_on_node*
				    sizeof(wilson_vector));
    if(srccs == NULL){
      node0_printf("Can't malloc space to read propagator\n");
      status = 1;
    }

    wpf = w_serial_w_i(filename);
    for(color = 0; color < 3; color++)for(spin = 0; spin < 4; spin++)
      {
	FORALLSITES(i,s){
	  srccs[i] = src[i].c[color].d[spin];
	}
	w_serial_w_from_field(wpf,spin,color,srccs);
      }
    free(srccs); srccs = NULL;
    w_serial_w_f(wpf);
    if(status == 0)
      node0_printf("Saved wprop serially to file %s\n",filename);
    break;

  case SAVE_SERIAL_FM:
  case SAVE_SERIAL_FM_SC:
    if(flag == SAVE_SERIAL_FM)
      wpf = w_serial_w_fm_i(filename);
    else
      wpf = w_serial_w_fm_sc_i(filename);
    wpf->prop = (wilson_propagator *)malloc(sites_on_node*
					    sizeof(wilson_propagator));
    wp = wpf->prop;
    if(wp == NULL){
      node0_printf("Can't malloc space to save propagator\n");
      status = 1;
      break;
    }
    
    /* Copy full propagator to buffer */
    FORALLSITES(i,s){
      wp[i] = src[i];
    }

    /* Convert to FNAL spin conventions */
    convert_wprop_milc_to_fnal_field(wp);

    /* Write the converted propagator */
    w_serial_w_fm_from_field(wpf, wp);
    w_serial_w_fm_f(wpf);
    break;

  case SAVE_SERIAL_SCIDAC:
  case SAVE_PARALLEL_SCIDAC:   
  case SAVE_MULTIFILE_SCIDAC: 
  case SAVE_PARTITION_SCIDAC:

    node0_printf("Currently saving the entire propagator in SciDAC format is not supported\n");
    status = 1;

    break;
  case SAVE_PARALLEL:
    wpf = w_parallel_w_i(filename);

    srccs = (wilson_vector *)malloc(sites_on_node*
				    4*sizeof(wilson_vector));
    if(srccs == NULL){
      node0_printf("Can't malloc space to read propagator\n");
      status = 1;
    }

    for(color = 0; color < 3; color++)for(spin = 0; spin < 4; spin++)
      {
	FORALLSITES(i,s){
	  for(spin = 0; spin < 4; spin++)
	    srccs[4*i + spin] = src[i].c[color].d[spin];
	}
	w_parallel_w_from_field(wpf,spin,color,srccs);
      }
    w_parallel_w_f(wpf);
    free(srccs); srccs = NULL;
    if(status == 0)
      node0_printf("Saved wprop in parallel to file %s\n",filename);
    break;

  case SAVE_CHECKPOINT:
    srccs = (wilson_vector *)malloc(sites_on_node*
				    sizeof(wilson_vector));
    if(srccs == NULL){
      node0_printf("Can't malloc space to read propagator\n");
      status = 1;
    }

    wpf = w_checkpoint_w_i(filename);
    for(color = 0; color < 3; color++)for(spin = 0; spin < 4; spin++)
      {
	FORALLSITES(i,s){
	  srccs[i] = src[i].c[color].d[spin];
	}
	w_checkpoint_w_from_field(wpf,spin,color,srccs);
      }
    free(srccs); srccs = NULL;
    w_checkpoint_w_f(wpf);
    if(status == 0)
      node0_printf("Saved wprop in checkpoint format to file %s\n",filename);
    break;

  case SAVE_MULTIDUMP:
    srccs = (wilson_vector *)malloc(sites_on_node*
				    sizeof(wilson_vector));
    if(srccs == NULL){
      node0_printf("Can't malloc space to read propagator\n");
      status = 1;
    }

    wpf = w_multidump_w_i(filename);
    for(color = 0; color < 3; color++)for(spin = 0; spin < 4; spin++)
      {
	FORALLSITES(i,s){
	  srccs[i] = src[i].c[color].d[spin];
	}
	w_multidump_w_from_field(wpf,spin,color,srccs);
      }
    free(srccs); srccs = NULL;
    w_multidump_w_f(wpf);
    if(status == 0)
      node0_printf("Saved wprop in multidump format to file %s\n",filename);
    break;
  default:
    node0_printf("save_wprop_from_field: Unrecognized save flag.\n");
    terminate(1);
  }
  
  if(timing)
    {
      dtime += dclock();
      if(flag != FORGET)
	node0_printf("Time to save wprop = %e\n",dtime);
    }
  return status;

} /* save_wprop_from_field */
/*---------------------------------------------------------------*/
void r_close_wprop(int flag, w_prop_file *wpf)
{
  
  switch(flag){
  case RELOAD_ASCII:
    r_ascii_w_f(wpf);
    break;
  case RELOAD_SERIAL:
    if(wpf->prop != NULL)
      free(wpf->prop); wpf->prop = NULL;
    r_serial_w_f(wpf);
    break;
  case RELOAD_PARALLEL:
    if(wpf->file_type == FILE_TYPE_W_PROP ||
       wpf->file_type == FILE_TYPE_W_PROP_1996)
      r_parallel_w_f(wpf); 
    else
      free(wpf);
    break;
  case RELOAD_MULTIDUMP:
    r_multidump_w_f(wpf); 
  }
}
/*---------------------------------------------------------------*/
void w_close_wprop(int flag, w_prop_file *wpf)
{
  wilson_propagator *wp;

  switch(flag){
  case SAVE_ASCII:
    w_ascii_w_f(wpf);
    break;
  case SAVE_SERIAL:
    w_serial_w_f(wpf); 
    break;
  case SAVE_SERIAL_FM:
  case SAVE_SERIAL_FM_SC:
    /* Dump accumulated propagator and free memory */
    wp = wpf->prop;
    if(wp != NULL){
      /* Convert from MILC to FNAL */
      convert_wprop_milc_to_fnal_field(wp);
      w_serial_w_fm_from_field(wpf, wp);
      free(wp);  wp = NULL;
    }
    w_serial_w_fm_f(wpf); 
    break;
  case SAVE_SERIAL_SCIDAC:
  case SAVE_PARALLEL_SCIDAC:   
  case SAVE_MULTIFILE_SCIDAC: 
  case SAVE_PARTITION_SCIDAC:
#ifdef HAVE_QIO
    close_usqcd_prop_write(wpf->outfile);
    if(wpf->prop != NULL)
      free(wpf->prop); 
    free(wpf);
#endif
    break;
  case SAVE_PARALLEL:
    w_parallel_w_f(wpf); 
    break;
  case SAVE_CHECKPOINT:
    w_checkpoint_w_f(wpf); 
    break;
  case SAVE_MULTIDUMP:
    w_multidump_w_f(wpf); 
  }
}

/*---------------------------------------------------------------*/
/* find out what kind of starting propagator to use, 
   and propagator name if necessary.  This routine is only 
   called by node 0.
   */
int ask_starting_wprop( FILE *fp, int prompt, int *flag, char *filename ){
  char *savebuf;
  int status;
  char myname[] = "ask_starting_wprop";
  
  if (prompt!=0) {
    printf("loading wilson propagator:\n enter 'fresh_wprop', ");
    printf("'continue_wprop', 'reload_ascii_wprop', ");
    printf("'reload_serial_wprop', 'reload_parallel_wprop', ");
    printf("or 'reload_multidump_wprop'\n");
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
  else if(strcmp("reload_multidump_wprop",savebuf) == 0 ) {
    *flag = RELOAD_MULTIDUMP;
  }
  else if(strcmp("reload_parallel_wprop",savebuf) == 0 ) {
    *flag = RELOAD_PARALLEL;
  }
  else{
    printf("ERROR IN INPUT: propagator_command is invalid\n"); return(1);
  }
  
  /*read name of file and load it */
  if( *flag != FRESH && *flag != CONTINUE ){
    if(prompt!=0)printf("enter name of file containing props\n");
    status=scanf("%s",filename);
    if(status !=1) {
      printf("\n%s(%d): ERROR IN INPUT: Can't read file name.\n",
	     myname, this_node);
      return(1);
    }
    printf("%s\n",filename);
  }
  return(0);
}

/*---------------------------------------------------------------*/
/* find out what do to with propagator at end, and propagator name if
   necessary.  This routine is only called by node 0.
   */
int ask_ending_wprop( FILE *fp, int prompt, int *flag, char *filename ){
  char *savebuf;
  int status;
  char myname[] = "ask_ending_wprop";
  
  if (prompt!=0) {
    printf("save wilson propagator:\n enter ");
    printf("'forget_wprop', 'save_ascii_wprop', 'save_serial_wprop', ");
    printf("'save_serial_fm_wprop', 'save_serial_fm_sc_wprop', ");
    printf("'save_serial_scidac_wprop', 'save_parallel_scidac_wprop', ");
    printf("'save_partfile_scidac_wprop', 'save_multifile_scidac_wprop', ");
    printf("'save_parallel_wprop', 'save_multidump_wprop', ");
    printf("or 'save_checkpoint_wprop'\n");
  }

  savebuf = get_next_tag(fp, "write wprop command", myname);
  if (savebuf == NULL)return 1;

  printf("%s ",savebuf);

  if(strcmp("save_ascii_wprop",savebuf) == 0 )  {
    *flag=SAVE_ASCII;
  }
  else if(strcmp("save_serial_wprop",savebuf) == 0 ) {
    *flag=SAVE_SERIAL;
  }
  else if(strcmp("save_parallel_wprop",savebuf) == 0 ) {
    *flag=SAVE_PARALLEL;
  }
  else if(strcmp("save_multidump_wprop",savebuf) == 0 ) {
    *flag=SAVE_MULTIDUMP;
  }
  else if(strcmp("save_checkpoint_wprop",savebuf) == 0 ) {
    *flag=SAVE_CHECKPOINT;
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
    *flag=SAVE_PARTITION_SCIDAC;
  }
  else if(strcmp("save_multifile_scidac_wprop",savebuf) == 0 ) {
    *flag=SAVE_MULTIFILE_SCIDAC;
  }
  else if(strcmp("forget_wprop",savebuf) == 0 ) {
    *flag=FORGET;
    printf("\n");
  }
  else {
    printf(" is not a valid save wprop command. INPUT ERROR.\n");
    return(1);
  }
  
  if( *flag != FORGET ){
    if(prompt!=0)printf("enter filename\n");
    status=scanf("%s",filename);
    if(status !=1){
      printf("\n%s(%d): ERROR IN INPUT. Can't read filename\n",
	     myname, this_node); 
      return(1);
    }
    printf("%s\n",filename);
  }
  return(0);
}
