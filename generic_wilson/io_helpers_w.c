/********************** io_helpers_w.c **********************************/
/* MIMD version 6 */
/* CD 10/97 from io_helpers.c DT 8/97 
   General purpose high level propagator I/O routines, 
   to be used by any application that wants them.
   */

#include "generic_wilson_includes.h"
#include "../include/io_lat.h"
#include "../include/io_wb.h"
#include "../include/file_types.h"
#include "../include/io_scidac.h"
#include "../include/io_scidac_w.h"
#include <string.h>
#ifdef HAVE_QIO
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

w_prop_file *r_open_prop(int flag, char *filename)
{
  w_prop_file *wpf;

  switch(flag){
  case RELOAD_ASCII:
    wpf = r_ascii_w_i(filename);
    break;
  case RELOAD_SERIAL:
    wpf = r_serial_w_i(filename);
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

w_prop_file *w_open_prop(int flag, char *filename)
{
  w_prop_file *wpf;
  
  switch(flag){
  case SAVE_ASCII:
    wpf = w_ascii_w_i(filename);
    break;

  case SAVE_SERIAL:
    wpf = w_serial_w_i(filename);
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
int reload_propagator( int flag, w_prop_file *wpf,
		       int spin, int color, field_offset dest, int timing)
{
  /* 0 normal exit value
     1 read error */

  double dtime;
  int i,status;
  site *s;

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
    status = r_serial_w(wpf,spin,color,dest); 
    break;
  case RELOAD_PARALLEL:
    /* Reopen, read, and close temporarily */
    r_parallel_w_o(wpf);
    status = r_parallel_w(wpf,spin,color,dest);
    r_parallel_w_c(wpf);
    break;
  case RELOAD_MULTIDUMP:
    /* Reopen, read, and close temporarily */
    r_multidump_w_o(wpf);
    status = r_multidump_w(wpf,spin,color,dest);
    r_multidump_w_c(wpf);
    break;
  default:
    node0_printf("reload_propagator: Unrecognized reload flag.\n");
    terminate(1);
  }
  
  if(timing)
    {
      dtime += dclock();
      if(flag != FRESH && flag != CONTINUE)
	node0_printf("Time to reload prop spin %d color %d %e\n",
		     spin,color,dtime);
    }

  return status;

} /* reload_propagator */

/*---------------------------------------------------------------*/
/* reload a full propagator (all source colors and spins) in any of
   the formats, or cold propagator, or keep current propagator:

   FRESH, CONTINUE,
   RELOAD_ASCII, RELOAD_SERIAL, RELOAD_PARALLEL, RELOAD_MULTIDUMP
   */
int reload_full_propagator( int flag, char *filename,
			    field_offset dest, int timing)
{
  /* 0 normal exit value
     1 read error */
  
  double dtime;
  int i,status;
  site *s;
  wilson_propagator *wp;
  int spin, color;
  field_offset destcs;
#ifdef HAVE_QIO
  field_offset destc;
#endif
  int file_type;
  w_prop_file *wpf;
  
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
    break;

  case RELOAD_SERIAL:
    /* Sniff out the file type */
    file_type = io_detect(filename, w_prop_list, N_WPROP_TYPES);

    if(file_type < 0){
      node0_printf("reload_ksprop: Can't read file %s\n", filename);
      return 1;
    }

    if(file_type == FILE_TYPE_W_PROP ||
       file_type == FILE_TYPE_W_PROP_1996)
      {
	node0_printf("Reading as a MILC Wilson prop file\n");
	/* Old MILC format has one record for each source spin, color */
	wpf = r_serial_w_i(filename);
	for(color = 0; color < 3; color++)for(spin = 0; spin < 4; spin++)
	  {
	    destcs = dest + color*sizeof(spin_wilson_vector) + 
	      spin*sizeof(wilson_vector);
	    if( r_serial_w(wpf,spin,color,dest) != 0)status = 1; 
	  }
	r_serial_w_f(wpf);
      }
    
    else if(file_type == FILE_TYPE_W_FMPROP){
      node0_printf("Reading as a Fermilab Wilson prop file\n");
      /* FNAL format has a full propagator in one record */
      r_prop_w_fm( filename, dest );
    }

    else if(file_type == FILE_TYPE_W_QIOPROP)
      {
#ifdef HAVE_QIO
	/* In this format we have three records, one for each
	   source color.  So one spin_wilson_vector field per record */
	QIO_Layout layout;
	QIO_Reader *infile;
	
	node0_printf("Reading as a SciDAC Wilson prop file\n");
	build_layout(&layout);
	infile = open_input(filename, &layout);
	if(infile == NULL)return 1;
	
	for(color = 0; color < 3; color++)
	  {
	    destc = dest + color*sizeof(spin_wilson_vector);
	    if(read_F3_D(infile, destc, 4) != QIO_SUCCESS)status = 1;
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
    break;
    
  case RELOAD_PARALLEL:
    /* Reopen, read, and close temporarily */
    r_parallel_w_i(filename);
    /* Old MILC format has one record for each source spin, color */
    for(color = 0; color < 3; color++)for(spin = 0; spin < 4; spin++)
      {
	destcs = dest + color*sizeof(spin_wilson_vector) + 
	  spin*sizeof(wilson_vector);
	if( r_parallel_w(wpf,spin,color,destcs) != 0)status = 1;
      }
    r_parallel_w_f(wpf);
    break;
    
  case RELOAD_MULTIDUMP:
    /* Reopen, read, and close temporarily */
    r_multidump_w_i(filename);
    /* Old MILC format has one record for each source spin, color */
    for(color = 0; color < 3; color++)for(spin = 0; spin < 4; spin++)
      {
	destcs = dest + color*sizeof(spin_wilson_vector) + 
	  spin*sizeof(wilson_vector);
	if(r_multidump_w(wpf,spin,color,destcs) != 0)status = 1;
      }
    r_multidump_w_f(wpf);
    break;
    
  default:
    node0_printf("reload_propagator: Unrecognized reload flag.\n");
    terminate(1);
  }
  
  if(timing)
    {
      dtime += dclock();
      if(flag != FRESH && flag != CONTINUE)
	node0_printf("Time to reload prop spin %d color %d %e\n",
		     spin,color,dtime);
    }
  
  return status;

} /* reload_full_propagator */

/*---------------------------------------------------------------*/
/* save a propagator one source color and spin at a time MILC formats only:
   FORGET,
   SAVE_ASCII, SAVE_SERIAL, SAVE_PARALLEL, SAVE_MULTIDUMP, SAVE_CHECKPOINT
*/
void save_propagator( int flag, w_prop_file *wpf, 
		      int spin, int color, field_offset src, int timing)
{
  double dtime;
  
  if(timing)dtime = -dclock();
  switch(flag){
  case FORGET:
    break;
  case SAVE_ASCII:
    node0_printf("Warning: no source spin and color recorded in an ASCII file\n");
    w_ascii_w(wpf,spin,color,src);
    break;
  case SAVE_SERIAL:
    w_serial_w(wpf,spin,color,src);
    break;
  case SAVE_PARALLEL:
    w_parallel_w_o(wpf);
    w_parallel_w(wpf,spin,color,src);
    w_parallel_w_c(wpf);
    break;
  case SAVE_CHECKPOINT:
    w_checkpoint_w_o(wpf);
    w_checkpoint_w(wpf,spin,color,src);
    w_checkpoint_w_c(wpf);
    break;
  case SAVE_MULTIDUMP:
    w_multidump_w_o(wpf);
    w_multidump_w(wpf,spin,color,src);
    w_multidump_w_c(wpf);
    break;
  default:
    node0_printf("save_propagator: Unrecognized save flag.\n");
    terminate(1);
  }
  
  if(timing)
    {
      dtime += dclock();
      if(flag != FORGET)
	node0_printf("Time to save prop spin %d color %d = %e\n",
		     spin,color,dtime);
    }
} /* save_propagator */

void save_full_propagator( int flag, char *filename,
			   field_offset src, int timing)
{
  double dtime;
  int spin, color;
  w_prop_file *wpf;
  field_offset srccs;
#ifdef HAVE_QIO
  field_offset srcc;
  int volfmt;
#endif
  
  if(timing)dtime = -dclock();
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
    break;

  case SAVE_SERIAL:
    wpf = w_serial_w_i(filename);
    for(color = 0; color < 3; color++)for(spin = 0; spin < 4; spin++)
      {
	srccs = src + color*sizeof(spin_wilson_vector) + 
	  spin*sizeof(wilson_vector);
	w_serial_w(wpf,spin,color,srccs);
      }
    w_serial_w_f(wpf);
    break;

  case SAVE_SERIAL_SCIDAC:
  case SAVE_PARTITION_SCIDAC:
  case SAVE_MULTIFILE_SCIDAC:
#ifdef HAVE_QIO
    {
      QIO_Layout layout;
      QIO_Writer *outfile;
      char recxml[64];

      build_layout(&layout);
      /* In this format we have three records, one for each
	 source color.  So one spin_wilson_vector field per record */
      if(flag == SAVE_SERIAL_SCIDAC)volfmt = QIO_SINGLEFILE;
      else if(flag == SAVE_PARTITION_SCIDAC)volfmt = QIO_PARTFILE;
      else if(flag == SAVE_MULTIFILE_SCIDAC)volfmt = QIO_MULTIFILE;

      build_layout(&layout);
      outfile = open_output(filename, volfmt, &layout,
			    "MILC Wilson propagator");
      if(outfile == NULL)break;

      for(color = 0; color < 3; color++)
	{
	  srcc = src + color*sizeof(spin_wilson_vector);
	  sprintf(recxml,"Source color %d\n",color);
	  if(write_F3_D(outfile, recxml, srcc, 4) != QIO_SUCCESS)break;
	}
      QIO_close_write(outfile);
    }
#else
    node0_printf("To write a SciDAC file requires QIO compilation\n");
#endif
    break;
  case SAVE_PARALLEL:
    wpf = w_parallel_w_i(filename);
    for(color = 0; color < 3; color++)for(spin = 0; spin < 4; spin++)
      {
	srccs = src + color*sizeof(spin_wilson_vector) + 
	  spin*sizeof(wilson_vector);
	w_parallel_w(wpf,spin,color,srccs);
      }
    w_parallel_w_f(wpf);
    break;
  case SAVE_CHECKPOINT:
    wpf = w_checkpoint_w_i(filename);
    for(color = 0; color < 3; color++)for(spin = 0; spin < 4; spin++)
      {
	srccs = src + color*sizeof(spin_wilson_vector) + 
	  spin*sizeof(wilson_vector);
	w_checkpoint_w(wpf,spin,color,srccs);
      }
    w_checkpoint_w_f(wpf);
    break;
  case SAVE_MULTIDUMP:
    wpf = w_multidump_w_i(filename);
    for(color = 0; color < 3; color++)for(spin = 0; spin < 4; spin++)
      {
	srccs = src + color*sizeof(spin_wilson_vector) + 
	  spin*sizeof(wilson_vector);
	w_multidump_w(wpf,spin,color,srccs);
      }
    w_multidump_w_f(wpf);
    break;
  default:
    node0_printf("save_propagator: Unrecognized save flag.\n");
    terminate(1);
  }
  
  if(timing)
    {
      dtime += dclock();
      if(flag != FORGET)
	node0_printf("Time to save prop spin %d color %d = %e\n",
		     spin,color,dtime);
    }
} /* save_full_propagator */
/*---------------------------------------------------------------*/
void r_close_prop(int flag, w_prop_file *wpf)
{
  
  switch(flag){
  case RELOAD_ASCII:
    r_ascii_w_f(wpf);
    break;
  case RELOAD_SERIAL:
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
void w_close_prop(int flag, w_prop_file *wpf)
{
  switch(flag){
  case SAVE_ASCII:
    w_ascii_w_f(wpf);
    break;
  case SAVE_SERIAL:
    w_serial_w_f(wpf); 
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
int ask_starting_prop( int prompt, int *flag, char *filename ){
  char savebuf[256];
  int status;
  
  if (prompt!=0) 
    printf("loading wilson propagator:\n enter 'fresh_prop','continue_prop', 'reload_ascii_prop', 'reload_serial_prop', 'reload_parallel_prop', or 'reload_multidump_prop'\n");
  status=scanf("%s",savebuf);
  if(status !=1) {
    printf("ERROR IN INPUT: 'fresh_prop' or 'reload_prop'\n");
    return(1);
  }

  printf("%s ",savebuf);
  if(strcmp("fresh_prop",savebuf) == 0 ){
    *flag = FRESH;
    printf("\n");
  }
  else if(strcmp("continue_prop",savebuf) == 0 ) {
    *flag = CONTINUE;
    printf("(if possible)\n");
  }
  else if(strcmp("reload_ascii_prop",savebuf) == 0 ) {
    *flag = RELOAD_ASCII;
  }
  else if(strcmp("reload_serial_prop",savebuf) == 0 ) {
    *flag = RELOAD_SERIAL;
  }
  else if(strcmp("reload_multidump_prop",savebuf) == 0 ) {
    *flag = RELOAD_MULTIDUMP;
  }
  else if(strcmp("reload_parallel_prop",savebuf) == 0 ) {
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
      printf("ERROR IN INPUT: file name read\n"); return(1);
    }
    printf("%s\n",filename);
  }
  return(0);
}

/*---------------------------------------------------------------*/
/* find out what do to with propagator at end, and propagator name if
   necessary.  This routine is only called by node 0.
   */
int ask_ending_prop( int prompt, int *flag, char *filename ){
  char savebuf[256];
  int status;
  
  if (prompt!=0) 
    printf("save wilson propagator:\n enter 'forget_prop', 'save_ascii_prop', 'save_serial_prop' , 'save_serial_scidac_prop', 'save_partfile_scidac_prop', 'save_multfile_scidac_prop', 'save_parallel_prop', 'save_multidump_prop', or 'save_checkpoint_prop'\n");
  status=scanf("%s",savebuf);
  if(status !=1) {
    printf("ERROR IN INPUT: 'save_prop' or 'forget_prop'\n");
    return(1);
  }
  printf("%s ",savebuf);
  if(strcmp("save_ascii_prop",savebuf) == 0 )  {
    *flag=SAVE_ASCII;
  }
  else if(strcmp("save_serial_prop",savebuf) == 0 ) {
    *flag=SAVE_SERIAL;
  }
  else if(strcmp("save_parallel_prop",savebuf) == 0 ) {
    *flag=SAVE_PARALLEL;
  }
  else if(strcmp("save_multidump_prop",savebuf) == 0 ) {
    *flag=SAVE_MULTIDUMP;
  }
  else if(strcmp("save_checkpoint_prop",savebuf) == 0 ) {
    *flag=SAVE_CHECKPOINT;
  }
  else if(strcmp("save_serial_scidac_prop",savebuf) == 0 ) {
    *flag=SAVE_SERIAL_SCIDAC;
  }
  else if(strcmp("save_partfile_scidac_prop",savebuf) == 0 ) {
    *flag=SAVE_PARTITION_SCIDAC;
  }
  else if(strcmp("save_multifile_scidac_prop",savebuf) == 0 ) {
    *flag=SAVE_MULTIFILE_SCIDAC;
  }
  else if(strcmp("forget_prop",savebuf) == 0 ) {
    *flag=FORGET;
    printf("\n");
  }
  else {
    printf("ERROR IN INPUT: %s is not a save propagator command\n",savebuf);
    return(1);
  }
  
  if( *flag != FORGET ){
    if(prompt!=0)printf("enter filename\n");
    status=scanf("%s",filename);
    if(status !=1){
      printf("ERROR IN INPUT: save filename\n"); return(1);
    }
    printf("%s\n",filename);
  }
  return(0);
}
