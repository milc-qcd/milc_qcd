/********************** io_helpers_ks.c *********************************/
/* MIMD version 6 */
/* Started 17 April 2002 by MBW, derived from io_helpers_w.c,
   pretty much cut-paste-query-replace.  Someday make one set of routines
   for lattices, Wilson, and KS props? */
/*
   General purpose high level propagator I/O routines, 
   to be used by any application that wants them.
*/

#include "generic_ks_includes.h"
#include "../include/io_lat.h"
#include "../include/io_prop_ks.h"
#include "../include/file_types.h"
#include <string.h>
#ifdef HAVE_QIO
#include <qio.h>
#endif

static file_type ksprop_list[N_KSPROP_TYPES] =
  { {FILE_TYPE_KSPROP,       KSPROP_VERSION_NUMBER},
    {FILE_TYPE_KSFMPROP,     KSFMPROP_VERSION_NUMBER},
    {FILE_TYPE_KSQIOPROP,    LIME_MAGIC_NO}
  };

/*---------------------------------------------------------------*/
/* reload a binary propagator in any of the formats */

int reload_serial_ksprop( int flag, int file_type, char *filename, 
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
      status = r_serial_ks(kspf,color,destc); 
    }
    r_serial_ks_f(kspf);
  }
  else if(file_type == FILE_TYPE_KSFMPROP){
    node0_printf("Reading as a Fermilab KS prop file\n");
    kspf = r_serial_ks_fm_i(filename);
    r_serial_ks_fm(kspf,dest);
    r_serial_ks_fm_f(kspf);
  }
  else if(file_type == FILE_TYPE_KSQIOPROP){
#ifdef HAVE_QIO
    node0_printf("Reading as a QIO KS prop file\n");
    restore_ks_vector_scidac(filename, dest, 3);
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

} /* reload_serial_ksprop */

/*---------------------------------------------------------------*/
/* reload a propagator in any of the formats, or cold propagator, or keep
   current propagator:
   FRESH, CONTINUE, RELOAD_ASCII, RELOAD_SERIAL
   */
int reload_ksprop( int flag, char *filename, field_offset dest, int timing)
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
    FORALLSITES(i,s) clearvec( (su3_vector *)F_PT(s,dest));
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

    status = reload_serial_ksprop(flag, file_type, filename, dest, timing);
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
    node0_printf("Reading as a QIO KS prop file\n");
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
void save_ksprop( int flag, char *filename, char *recxml, 
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
      w_serial_ks(kspf,color,srcc);
    }
    w_serial_ks_f(kspf);
    break;
  case SAVE_SERIAL_FM:
    kspf = w_serial_ks_fm_i(filename);
    w_serial_ks_fm(kspf,src);
    w_serial_ks_fm_f(kspf);
    break;
  case SAVE_SERIAL_SCIDAC:
#ifdef HAVE_QIO
    save_ks_vector_scidac(filename, recxml, QIO_SINGLEFILE, src, 3);
#else
    node0_printf("Need QIO compilation to save in SciDAC format\n");
#endif
    break;
  case SAVE_PARTITION_SCIDAC:
#ifdef HAVE_QIO
    save_ks_vector_scidac(filename, recxml, QIO_PARTFILE, src, 3);
#else
    node0_printf("Need QIO compilation to save in SciDAC format\n");
#endif
    break;
  case SAVE_MULTIFILE_SCIDAC:
#ifdef HAVE_QIO
    save_ks_vector_scidac(filename, recxml, QIO_MULTIFILE, src, 3);
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

} /* save_ksprop */


/* find out what if any KS propagator should be loaded.
   This routine is only called by node 0.
*/
int ask_starting_ksprop( int prompt, int *flag, char *filename ){
    char savebuf[256];
    int status;

    if (prompt!=0) 
      printf( "enter 'fresh_ks', 'reload_ks_ascii', 'reload_ks_serial', \n");
    status=scanf("%s",savebuf);
    if (status == EOF){
      printf("ask_starting_ksprop: EOF on STDIN.\n");
      return(1);
    }
    if(status !=1) {
        printf("ask_starting_ksprop: ERROR IN INPUT: error reading starting prop command\n");
        return(1);
    }

    printf("%s ",savebuf);
    if(strcmp("fresh_ks",savebuf) == 0 ){
       *flag = FRESH;
    printf("\n");
    }
    else if(strcmp("reload_ks_ascii",savebuf) == 0 ) {
       *flag = RELOAD_ASCII;
    }
    else if(strcmp("reload_ks_serial",savebuf) == 0 ) {
       *flag = RELOAD_SERIAL;
    }
    else{
    	printf("ask_starting_ksprop: ERROR IN INPUT: command %s is invalid\n",
	       savebuf); 
	return(1);
    }

    /*read name of file and load it */
    if( *flag != FRESH ){
        if(prompt!=0) printf("enter name of file containing ksprop\n");
        status = scanf("%s",filename);
        if(status != 1) {
	    printf("ask_starting_ksprop: ERROR IN INPUT: error reading filename\n"); 
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
        "'forget_ks', 'save_ks_ascii', 'save_ks_serial_tslice', 'save_ks_serial', 'save_ks_serial_fm', 'save_ks_serial_scidac', 'save_ks_multifile_scidac', 'save_ks_partition_scidac' ?\n");
    status=scanf("%s",savebuf);
    if(status !=1) {
        printf("ask_ending_ksprop: ERROR IN INPUT: ending ksprop command\n");
        return(1);
    }
    printf("%s ",savebuf);
    if(strcmp("save_ks_ascii",savebuf) == 0 )  {
        *flag=SAVE_ASCII;
    }
    else if(strcmp("save_ks_serial_tslice",savebuf) == 0 ) {
        *flag=SAVE_SERIAL_TSLICE;
    }
    else if(strcmp("save_ks_serial",savebuf) == 0 ) {
        *flag=SAVE_SERIAL;
    }
    else if(strcmp("save_ks_serial_fm",savebuf) == 0 ) {
        *flag=SAVE_SERIAL_FM;
    }
    else if(strcmp("save_ks_serial_scidac",savebuf) == 0 ) {
#ifdef HAVE_QIO
        *flag=SAVE_SERIAL_SCIDAC;
#else
	node0_printf("requires QIO compilation!\n");
	terminate(1);
#endif
    }
    else if(strcmp("save_ks_multifile_scidac",savebuf) == 0 ) {
#ifdef HAVE_QIO
        *flag=SAVE_MULTIFILE_SCIDAC;
#else
	node0_printf("requires QIO compilation!\n");
	terminate(1);
#endif
    }
    else if(strcmp("save_ks_partition_scidac",savebuf) == 0 ) {
#ifdef HAVE_QIO
        *flag=SAVE_PARTITION_SCIDAC;
#else
	node0_printf("requires QIO compilation!\n");
	terminate(1);
#endif
    }
    else if(strcmp("forget_ks",savebuf) == 0 ) {
        *flag=FORGET;
	printf("\n");
    }
    else {
      printf("ask_ending_ksprop: ERROR IN INPUT: %s is not a valid save KS prop command\n",
	     savebuf);
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
