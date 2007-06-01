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
#include "../include/io_scidac.h"
#include "../include/io_scidac_w.h"
#include <qio.h>
#endif

static file_type w_prop_list[N_WPROP_TYPES] =
  { {FILE_TYPE_W_PROP,       W_PROP_VERSION_NUMBER},
    {FILE_TYPE_W_PROP_1996,  W_PROP_VERSION_NUMBER_1996},
    {FILE_TYPE_W_FMPROP,     W_FMPROP_VERSION_NUMBER},
    {FILE_TYPE_W_QIOPROP,    LIME_MAGIC_NO}
  };

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
      node0_printf("r_open_wprop: Can't read file %s\n", filename);
      return NULL;
    }
    /* If it is an FNAL propagator file, special treatment */
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
      /* For either format we read the entire file */
      r_serial_w_fm_to_field(wpf, wp);

      /* Convert from FNAL to MILC */
      convert_wprop_fnal_to_milc_field(wp);
    }
    /* Otherwise, just open the file */
    else {
      wpf = r_serial_w_i(filename);
    }
    wpf->file_type = file_type;
    break;
  case RELOAD_PARALLEL:
    wpf = r_parallel_w_i(filename);
    r_parallel_w_c(wpf);
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
   source spin and color at a time.
   Supported only in MILC formats */

w_prop_file *w_open_wprop(int flag, char *filename)
{
  w_prop_file *wpf;
  wilson_propagator *wp;
  
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
    wp = wpf->prop;
    /* Special treatment for FNAL formatted file */
    if(wpf->file_type == FILE_TYPE_W_FMPROP){
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
		       int spin, int color, wilson_vector *dest, int timing)
{
  /* 0 normal exit value
     1 read error */

  double dtime = 0;
  int i,status;
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
    FORALLSITES(i,s)clear_wvec( &(dest[i]) );
    break;
  case RELOAD_ASCII:
    node0_printf("Reloading ASCII to temp wprop not supported\n");
    terminate(1);
    break;
  case RELOAD_SERIAL:
    wp = wpf->prop;
    /* Special treatment for FNAL formatted file */
    if(wpf->file_type == FILE_TYPE_W_FMPROP){
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
    else {
      status = r_serial_w_to_field(wpf,spin,color,dest); 
    }
    break;
  case RELOAD_PARALLEL:
    /* Reopen, read, and close temporarily */
    r_parallel_w_o(wpf);
    status = r_parallel_w_to_field(wpf,spin,color,dest);
    r_parallel_w_c(wpf);
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
#ifdef HAVE_QIO
  field_offset destc;
#endif
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

    else if(file_type == FILE_TYPE_W_QIOPROP)
      {
#ifdef HAVE_QIO
	/* In this format we have three records, one for each
	   source color.  So one spin_wilson_vector field per record */
	QIO_Layout layout;
	QIO_Filesystem fs;
	QIO_Reader *infile;
	
	build_qio_layout(&layout);
	build_qio_filesystem(&fs);

	if(flag == RELOAD_SERIAL){
	  node0_printf("Reading serially as a SciDAC Wilson prop file\n");
	  infile = open_scidac_input(filename, &layout, &fs, QIO_SERIAL);
	}
	else{
	  node0_printf("Reading in parallel as a SciDAC Wilson prop file\n");
	  infile = open_scidac_input(filename, &layout, &fs, QIO_PARALLEL);
	}

	if(infile == NULL)return 1;
	
	for(color = 0; color < 3; color++)
	  {
	    destc = dest + color*sizeof(spin_wilson_vector);
	    if(read_F3_D_to_site(infile, destc, 4) != QIO_SUCCESS)status = 1;
	  }
	QIO_close_read(infile);
#else
	node0_printf("This looks like a QIO file, but to read it requires QIO compilation\n");
	return 1;
#endif
      }

    else{
      node0_printf("Unsupported file type %d\n",file_type);
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

    else if(file_type == FILE_TYPE_W_QIOPROP)
      {
#ifdef HAVE_QIO
	/* In this format we have three records, one for each
	   source color.  So one spin_wilson_vector field per record */
	QIO_Layout layout;
	QIO_Filesystem fs;
	QIO_Reader *infile;

	build_qio_layout(&layout);
	build_qio_filesystem(&fs);
	if(flag == RELOAD_SERIAL){
	  node0_printf("Reading serially as a SciDAC Wilson prop file\n");
	  infile = open_scidac_input(filename, &layout, &fs, QIO_SERIAL);
	}
	else{
	  node0_printf("Reading in parallel as a SciDAC Wilson prop file\n");
	  infile = open_scidac_input(filename, &layout, &fs, QIO_PARALLEL);
	}

	if(infile == NULL)return 1;

	destcs = (wilson_vector *)malloc(sites_on_node*
					 4*sizeof(wilson_vector));
	if(destcs == NULL){
	  node0_printf("Can't malloc space to read propagator\n");
	  status = 1;
	}

	for(color = 0; color < 3; color++)
	  {
	    if(status == 0){
	      if(read_F3_D_to_field(infile, destcs, 4) != QIO_SUCCESS)
		status += 1;
	      FORALLSITES(i,s){
		for(spin = 0; spin < 4; spin++)
		  dest[i].c[color].d[spin] = destcs[4*i + spin];
	      }
	    }
	  }
	free(destcs); destcs = NULL;
	QIO_close_read(infile);
#else
	node0_printf("This looks like a QIO file, but to read it requires QIO compilation\n");
	return 1;
#endif
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
/* save a propagator one source color and spin at a time MILC formats only:
   FORGET,
   SAVE_ASCII, SAVE_SERIAL, SAVE_PARALLEL, SAVE_MULTIDUMP, SAVE_CHECKPOINT
*/
int save_wprop_sc_from_field( int flag, w_prop_file *wpf, 
		      int spin, int color, wilson_vector *src, int timing)
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
#ifdef HAVE_QIO
  field_offset srcc;
  int volfmt = QIO_SINGLEFILE;
#endif
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
  case SAVE_PARTITION_SCIDAC:
  case SAVE_MULTIFILE_SCIDAC:
#ifdef HAVE_QIO
    {
      QIO_Layout layout;
      QIO_Filesystem fs;
      QIO_Writer *outfile;

      build_qio_layout(&layout);
      build_qio_filesystem(&fs);
      /* In this format we have three records, one for each
	 source color.  So one spin_wilson_vector field per record */
      if(flag == SAVE_SERIAL_SCIDAC)volfmt = QIO_SINGLEFILE;
      else if(flag == SAVE_PARALLEL_SCIDAC)volfmt = QIO_SINGLEFILE;
      else if(flag == SAVE_PARTITION_SCIDAC)volfmt = QIO_PARTFILE;
      else if(flag == SAVE_MULTIFILE_SCIDAC)volfmt = QIO_MULTIFILE;

      build_qio_layout(&layout);
      build_qio_filesystem(&fs);
      if(flag == SAVE_PARALLEL_SCIDAC)
	outfile = open_scidac_output(filename, volfmt, QIO_PARALLEL, 
				     QIO_ILDGNO, NULL, &layout, &fs,
				     "MILC Wilson propagator");
      else
	outfile = open_scidac_output(filename, volfmt, QIO_SERIAL, QIO_ILDGNO,
				     NULL, &layout, &fs, 
				     "MILC Wilson propagator");

      if(outfile == NULL)break;

      for(color = 0; color < 3; color++)
	{
	  srcc = src + color*sizeof(spin_wilson_vector);
	  if(write_F3_D_from_site(outfile, recxml, srcc, 4) != QIO_SUCCESS)
	    status += 1;
	}
      close_output(outfile);
    }
#else
    node0_printf("To write a SciDAC file requires QIO compilation\n");
#endif
    if(status == 0)
      node0_printf("Saved wprop serially to file %s\n",filename);
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

  if(file_type == FILE_TYPE_W_PROP ||
     file_type == FILE_TYPE_W_PROP_1996){
    *ndim = 4;
    nx = -1; ny = -1; nz = -1; nt = -1;
    wpf = r_serial_w_i(filename);
    for(i = 0; i < *ndim; i++)
      dims[i] = wpf->header->dims[i];
    r_serial_w_f(wpf);
  }
  else if(file_type == FILE_TYPE_W_FMPROP){
    *ndim = 4;
    nx = -1; ny = -1; nz = -1; nt = -1;
    wpf = r_serial_w_fm_i(filename);
    for(i = 0; i < *ndim; i++)
      dims[i] = wpf->header->dims[i];
    r_serial_w_fm_f(wpf);
  }
  else if(file_type == FILE_TYPE_W_QIOPROP){
#ifdef HAVE_QIO
    read_lat_dim_scidac(filename, ndim, dims);
#else
    node0_printf("This looks like a QIO file, but to read it requires QIO compilation\n");
    return 1;
#endif
  }
  else{
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
#ifdef HAVE_QIO
  int volfmt = QIO_SINGLEFILE;
#endif
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
  case SAVE_PARTITION_SCIDAC:
  case SAVE_MULTIFILE_SCIDAC:
#ifdef HAVE_QIO
    {
      QIO_Layout layout;
      QIO_Filesystem fs;
      QIO_Writer *outfile;

      build_qio_layout(&layout);
      build_qio_filesystem(&fs);
      /* In this format we have three records, one for each
	 source color.  So one spin_wilson_vector field per record */
      if(flag == SAVE_SERIAL_SCIDAC)volfmt = QIO_SINGLEFILE;
      else if(flag == SAVE_PARALLEL_SCIDAC)volfmt = QIO_SINGLEFILE;
      else if(flag == SAVE_PARTITION_SCIDAC)volfmt = QIO_PARTFILE;
      else if(flag == SAVE_MULTIFILE_SCIDAC)volfmt = QIO_MULTIFILE;

      build_qio_layout(&layout);
      build_qio_filesystem(&fs);
      if(flag == SAVE_PARALLEL_SCIDAC)
	outfile = open_scidac_output(filename, volfmt, QIO_PARALLEL, 
				     QIO_ILDGNO, NULL, &layout, &fs,
				     "MILC Wilson propagator");
      else
	outfile = open_scidac_output(filename, volfmt, QIO_SERIAL, 
				     QIO_ILDGNO, NULL, &layout, &fs,
				     "MILC Wilson propagator");
      if(outfile == NULL)break;

      srccs = (wilson_vector *)malloc(sites_on_node*
				       4*sizeof(wilson_vector));
      if(srccs == NULL){
	node0_printf("Can't malloc space to read propagator\n");
	status = 1;
      }

      for(color = 0; color < 3; color++)
	{
	  FORALLSITES(i,s){
	    for(spin = 0; spin < 4; spin++)
	      srccs[4*i + spin] = src[i].c[color].d[spin];
	  }
	  if(write_F3_D_from_field(outfile, recxml, srccs, 4) 
	     != QIO_SUCCESS)break;
	}
      close_output(outfile);
    }
    free(srccs); srccs = NULL;
#else
    node0_printf("To write a SciDAC file requires QIO compilation\n");
#endif
    if(status == 0)
      node0_printf("Saved wprop serially to file %s\n",filename);
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
    r_parallel_w_f(wpf); 
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
    printf("'save_partfile_scidac_wprop', 'save_multfile_scidac_wprop', ");
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
      printf("\n%s(%d). ERROR IN INPUT. Can't read filename\n",
	     myname, this_node); 
      return(1);
    }
    printf("%s\n",filename);
  }
  return(0);
}
