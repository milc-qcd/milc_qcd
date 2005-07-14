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

static file_type ksprop_list[N_KSPROP_TYPES] =
  { {FILE_TYPE_KSPROP,       KSPROP_VERSION_NUMBER},
    {FILE_TYPE_KSFMPROP,     KSFMPROP_VERSION_NUMBER},
    {FILE_TYPE_KSQIOPROP,    LIME_MAGIC_NO}
  };

/*---------------------------------------------------------------*/
/* reload a binary propagator in any of the formats */

int reload_serial_ksprop_to_site( int flag, int file_type, char *filename, 
			  field_offset dest, int timing)
{
  /* 0 normal exit value
     1 read error */

  ks_prop_file *kspf;
  int status,color;
  field_offset destc;

  if(file_type == FILE_TYPE_KSPROP){
    node0_printf("Reading as a standard KS prop file\n");
    kspf = r_serial_ks_i(filename);
    for(color = 0; color < 3; color++){
      destc = dest + color*sizeof(su3_vector);
      status = r_serial_ks_to_site(kspf,color,destc); 
    }
    r_serial_ks_f(kspf);
  }
  else if(file_type == FILE_TYPE_KSFMPROP){
    node0_printf("Reading as a Fermilab KS prop file\n");
    kspf = r_serial_ks_fm_i(filename);
    r_serial_ks_fm_to_site(kspf,dest);
    r_serial_ks_fm_f(kspf);
  }
  else if(file_type == FILE_TYPE_KSQIOPROP){
#ifdef HAVE_QIO
    node0_printf("Reading as a QIO KS prop file\n");
    restore_ks_vector_scidac_to_site(filename, dest, 3);
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

} /* reload_serial_ksprop_to_site */

/*---------------------------------------------------------------*/
/* reload a binary propagator in any of the formats */

int reload_serial_ksprop_to_field( int flag, int file_type, char *filename, 
				   su3_vector *dest, int timing)
{
  /* 0 normal exit value
     1 read error */

  ks_prop_file *kspf;
  int status,color;

  if(file_type == FILE_TYPE_KSPROP){
    node0_printf("Reading as a standard KS prop file\n");
    kspf = r_serial_ks_i(filename);
    for(color = 0; color < 3; color++){
      status = r_serial_ks_to_field(kspf,color,dest); 
    }
    r_serial_ks_f(kspf);
  }
  else if(file_type == FILE_TYPE_KSFMPROP){
    node0_printf("Reading as a Fermilab KS prop file\n");
    kspf = r_serial_ks_fm_i(filename);
    r_serial_ks_fm_to_field(kspf,dest);
    r_serial_ks_fm_f(kspf);
  }
  else if(file_type == FILE_TYPE_KSQIOPROP){
#ifdef HAVE_QIO
    node0_printf("Reading as a QIO KS prop file\n");
    restore_ks_vector_scidac_to_field(filename, dest, 3);
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

} /* reload_serial_ksprop_to_field */

/*---------------------------------------------------------------*/
/* reload a propagator in any of the formats, or cold propagator, or keep
   current propagator:
   FRESH, CONTINUE, RELOAD_ASCII, RELOAD_SERIAL
   */
int reload_ksprop_to_site( int flag, char *filename, field_offset dest, 
			   int timing)
{
  /* 0 normal exit value
     1 read error */

  double dtime;
  int i,status,color;
  field_offset destc;
  int file_type;
  ks_prop_file *kspf;
  site *s;

  if(timing)dtime = -dclock();
  status = 0;
  switch(flag){
  case CONTINUE:  /* do nothing */
    break;
  case FRESH:     /* zero initial guess */
    FORALLSITES(i,s) {
      clearvec( (su3_vector *)F_PT(s,dest));
      clearvec( (su3_vector *)F_PT(s,dest)+1);
      clearvec( (su3_vector *)F_PT(s,dest)+2);
    }
    break;
  case RELOAD_ASCII:
    kspf = r_ascii_ks_i(filename);
    for(color = 0; color < 3; color++){
      destc = dest + color*sizeof(su3_vector);
      status = r_ascii_ks(kspf,color,destc);
    }
    r_ascii_ks_f(kspf);
    break;
  case RELOAD_SERIAL:
    file_type = io_detect(filename, ksprop_list, N_KSPROP_TYPES);
    if(file_type < 0){
      node0_printf("reload_ksprop: Can't read file %s\n", filename);
      return 1;
    }

    status = reload_serial_ksprop_to_site(flag, file_type, filename, dest, timing);
    if(status != 0)return status;

    break;
  default:
    node0_printf("reload_ksprop: Unrecognized reload flag.\n");
    terminate(1);
  }
  
  if(timing)
    {
      dtime += dclock();
      if(flag != FRESH && flag != CONTINUE)
	node0_printf("Time to reload ksprop = %e\n", dtime);
    }

  return status;

} /* reload_ksprop */

/*---------------------------------------------------------------*/
/* reload a propagator in any of the formats, or cold propagator, or keep
   current propagator:
   FRESH, CONTINUE, RELOAD_ASCII, RELOAD_SERIAL
   */
int reload_ksprop_to_field( int flag, char *filename, 
			    su3_vector *dest, int timing)
{
  /* 0 normal exit value
     1 read error */

  double dtime;
  int i,status,color;
  field_offset destc;
  int file_type;
  ks_prop_file *kspf;
  site *s;

  if(timing)dtime = -dclock();
  status = 0;
  switch(flag){
  case CONTINUE:  /* do nothing */
    break;
  case FRESH:     /* zero initial guess */
    FORALLSITES(i,s) {
      for(color = 0; color < 3; color++)
	clearvec( dest + 3*i + color );
    }
    break;
  case RELOAD_ASCII: node0_printf("reload_ksprop_to_field: ASCII to field not supported\n");
    terminate(1);
    break;
  case RELOAD_SERIAL:
    file_type = io_detect(filename, ksprop_list, N_KSPROP_TYPES);
    if(file_type < 0){
      node0_printf("reload_ksprop_to_field: Can't read file %s\n", filename);
      return 1;
    }
    status = reload_serial_ksprop_to_field(flag, file_type, filename, 
					   dest, timing);
    if(status != 0)return status;

    break;
  default:
    node0_printf("reload_ksprop_to_field: Unrecognized reload flag.\n");
    terminate(1);
  }
  
  if(timing)
    {
      dtime += dclock();
      if(flag != FRESH && flag != CONTINUE)
	node0_printf("Time to reload ksprop = %e\n", dtime);
    }

  return status;

} /* reload_ksprop_to_field */

/*---------------------------------------------------------------*/
/* read the lattice dimensions from a binary ks prop file */

int read_lat_dim_ksprop(char *filename, int file_type, int *ndim, int dims[])
{
  ks_prop_file *kspf;
  int i;

  if(file_type == FILE_TYPE_KSPROP){
    *ndim = 4;
    nx = -1; ny = -1; nz = -1; nt = -1;
    kspf = r_serial_ks_i(filename);
    for(i = 0; i < *ndim; i++)
      dims[i] = kspf->header->dims[i];
    r_serial_ks_f(kspf);
  }
  else if(file_type == FILE_TYPE_KSFMPROP){
    *ndim = 4;
    nx = -1; ny = -1; nz = -1; nt = -1;
    kspf = r_serial_ks_fm_i(filename);
    for(i = 0; i < *ndim; i++)
      dims[i] = kspf->header->dims[i];
    r_serial_ks_fm_f(kspf);
  }
  else if(file_type == FILE_TYPE_KSQIOPROP){
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
/* write a single su3_vector field to an open file in various formats
   FORGET, SAVE_ASCII, SAVE_SERIAL_TSLICE, SAVE_SERIAL,
   SERIAL_SERIAL_FM, SAVE_XXX_SCIDAC
   */
void save_ksprop_from_site( int flag, char *filename, char *recxml, 
		  field_offset src, int timing)
{
  double dtime;
  int color;
  field_offset srcc;
  ks_prop_file *kspf;
  
  if(timing)dtime = -dclock();
  switch(flag){
  case FORGET:
    break;
  case SAVE_ASCII:
    kspf = w_ascii_ks_i(filename);
    for(color = 0; color < 3; color++){
      srcc = src + color*sizeof(su3_vector);
      w_ascii_ks(kspf,color,srcc);
    }
    w_ascii_ks_f(kspf);
    break;
  case SAVE_SERIAL_TSLICE:
    w_ascii_ksprop_tt(filename, src);
    break;
  case SAVE_SERIAL:
    kspf = w_serial_ks_i(filename);
    for(color = 0; color < 3; color++){
      srcc = src + color*sizeof(su3_vector);
      w_serial_ks_from_site(kspf,color,srcc);
    }
    w_serial_ks_f(kspf);
    break;
  case SAVE_SERIAL_FM:
    kspf = w_serial_ks_fm_i(filename);
    w_serial_ks_fm_from_site(kspf,src);
    w_serial_ks_fm_f(kspf);
    break;
  case SAVE_SERIAL_SCIDAC:
#ifdef HAVE_QIO
    save_ks_vector_scidac_from_site(filename, recxml, QIO_SINGLEFILE, src, 3);
#else
    node0_printf("Need QIO compilation to save in SciDAC format\n");
#endif
    break;
  case SAVE_PARTITION_SCIDAC:
#ifdef HAVE_QIO
    save_ks_vector_scidac_from_site(filename, recxml, QIO_PARTFILE, src, 3);
#else
    node0_printf("Need QIO compilation to save in SciDAC format\n");
#endif
    break;
  case SAVE_MULTIFILE_SCIDAC:
#ifdef HAVE_QIO
    save_ks_vector_scidac_from_site(filename, recxml, QIO_MULTIFILE, src, 3);
#else
    node0_printf("Need QIO compilation to save in SciDAC format\n");
#endif
    break;
  default:
    node0_printf("save_ksprop: Unrecognized save flag.\n");
    terminate(1);
  }
  
  if(timing)
    {
      dtime += dclock();
      if(flag != FORGET)
	node0_printf("Time to save ksprop = %e\n", dtime);
    }

} /* save_ksprop_from_site */

/*---------------------------------------------------------------*/
/* write a single su3_vector field to an open file in various formats
   FORGET, SAVE_ASCII, SAVE_SERIAL_TSLICE, SAVE_SERIAL,
   SERIAL_SERIAL_FM, SAVE_XXX_SCIDAC
   */
void save_ksprop_from_field( int flag, char *filename, char *recxml, 
		  su3_vector *src, int timing)
{
  double dtime;
  int color;
  field_offset srcc;
  ks_prop_file *kspf;
  
  if(timing)dtime = -dclock();
  switch(flag){
  case FORGET:
    break;
  case SAVE_ASCII:
    node0_printf("Saving ASCII from fields not supported\n");
    break;
  case SAVE_SERIAL_TSLICE:
    node0_printf("Saving ASCII TSLICE from fields not supported\n");
    break;
  case SAVE_SERIAL:
    kspf = w_serial_ks_i(filename);
    for(color = 0; color < 3; color++){
      w_serial_ks_from_field(kspf,color,src);
    }
    w_serial_ks_f(kspf);
    break;
  case SAVE_SERIAL_FM:
    kspf = w_serial_ks_fm_i(filename);
    w_serial_ks_fm_from_field(kspf,src);
    w_serial_ks_fm_f(kspf);
    break;
  case SAVE_SERIAL_SCIDAC:
#ifdef HAVE_QIO
    save_ks_vector_scidac_from_field(filename, recxml, QIO_SINGLEFILE, src, 3);
#else
    node0_printf("Need QIO compilation to save in SciDAC format\n");
#endif
    break;
  case SAVE_PARTITION_SCIDAC:
#ifdef HAVE_QIO
    save_ks_vector_scidac_from_field(filename, recxml, QIO_PARTFILE, src, 3);
#else
    node0_printf("Need QIO compilation to save in SciDAC format\n");
#endif
    break;
  case SAVE_MULTIFILE_SCIDAC:
#ifdef HAVE_QIO
    save_ks_vector_scidac_from_field(filename, recxml, QIO_MULTIFILE, src, 3);
#else
    node0_printf("Need QIO compilation to save in SciDAC format\n");
#endif
    break;
  default:
    node0_printf("save_ksprop: Unrecognized save flag.\n");
    terminate(1);
  }
  
  if(timing)
    {
      dtime += dclock();
      if(flag != FORGET)
	node0_printf("Time to save ksprop = %e\n", dtime);
    }

} /* save_ksprop_from_site */


/* find out what if any KS propagator should be loaded.
   This routine is only called by node 0.
*/
int ask_starting_ksprop( int prompt, int *flag, char *filename ){
    char savebuf[256];
    int status;

    if (prompt!=0) 
      printf( "enter 'fresh_ks', 'reload_ascii_ksprop', 'reload_serial_ksprop' \n");
    status=scanf("%s",savebuf);
    if (status == EOF){
      printf("ask_starting_ksprop: EOF on STDIN.\n");
      return(1);
    }
    if(status !=1) {
        printf("ask_starting_ksprop: ERROR IN INPUT: error reading starting ksprop command\n");
        return(1);
    }

    printf("%s ",savebuf);
    if(strcmp("fresh_ks",savebuf) == 0 ){
       *flag = FRESH;
    printf("\n");
    }
    else if(strcmp("reload_ascii_ksprop",savebuf) == 0 ) {
       *flag = RELOAD_ASCII;
    }
    else if(strcmp("reload_serial_ksprop",savebuf) == 0 ) {
       *flag = RELOAD_SERIAL;
    }
    else{
      printf("is not a valid starting ksprop command: INPUT ERROR\n");
	return(1);
    }

    /*read name of file and load it */
    if( *flag != FRESH ){
        if(prompt!=0) printf("enter name of file containing ksprop\n");
        status = scanf("%s",filename);
        if(status != 1) {
	  printf("\nask_starting_ksprop: ERROR IN INPUT: Can't read filename\n");
	    return(1);
        }
	printf("%s\n",filename);
    }
    return(0);

} /* end ask_starting_ksprop() */


/* find out what do to with lattice at end, and lattice name if
   necessary.  This routine is only called by node 0.
*/
int ask_ending_ksprop( int prompt, int *flag, char *filename ){
    char savebuf[256];
    int status;

    if (prompt!=0) printf(
        "'forget_ksprop', 'save_ascii_ksprop', 'save_serial_ksprop_tslice', 'save_serial_ksprop', 'save_serial_fm_ksprop', 'save_serial_scidac_ksprop', 'save_multifile_scidac_ksprop', 'save_partition_scidac_ksprop' ?\n");

    status=scanf("%s",savebuf);
    if (status == EOF){
      printf("ask_ending_ksprop: EOF on STDIN.\n");
      return(1);
    }
    if(status !=1) {
        printf("ask_ending_ksprop: ERROR IN INPUT: ending ksprop command\n");
        return(1);
    }
    printf("%s ",savebuf);

    if(strcmp("save_ascii_ksprop",savebuf) == 0 )  {
        *flag=SAVE_ASCII;
    }
    else if(strcmp("save_serial_ksprop_tslice",savebuf) == 0 ) {
        *flag=SAVE_SERIAL_TSLICE;
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
    else if(strcmp("save_multifile_scidac_ksprop",savebuf) == 0 ) {
#ifdef HAVE_QIO
        *flag=SAVE_MULTIFILE_SCIDAC;
#else
	node0_printf("requires QIO compilation!\n");
	terminate(1);
#endif
    }
    else if(strcmp("save_partition_scidac_ksprop",savebuf) == 0 ) {
#ifdef HAVE_QIO
        *flag=SAVE_PARTITION_SCIDAC;
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
      return(1);
    }

    if( *flag != FORGET ){
        if(prompt!=0)printf("enter filename\n");
        status = scanf("%s",filename);
        if(status != 1){
    	    printf("ask_ending_ksprop: ERROR IN INPUT: error reading filename\n"); 
	    return(1);
        }
	printf("%s\n",filename);
    }
    return(0);

} /* end ask_ending_ksprop() */
