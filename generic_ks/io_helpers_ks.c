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

/*---------------------------------------------------------------*/
/* Open propagator file for reading */
ks_prop_file *r_open_ksprop(int flag, char *filename)
{
  ks_prop_file *kspf;
  
  switch(flag){
  case RELOAD_ASCII:
    kspf = r_ascii_ks_i(filename);
    break;
  case RELOAD_SERIAL:
    kspf = r_serial_ks_i(filename);
    break;
  default:
    kspf = NULL;
  }
  return kspf;
}

/*---------------------------------------------------------------*/
/* Open propagator file for writing */
ks_prop_file *w_open_ksprop(int flag, char *filename)
{
  ks_prop_file *kspf;
  
  switch(flag){
  case SAVE_ASCII:
    kspf = w_ascii_ks_i(filename);
    break;
  case SAVE_SERIAL:
    kspf = w_serial_ks_i(filename);
    break;
  case SAVE_SERIAL_FM:
    kspf = w_serial_ks_fm_i(filename);
    break;
  case SAVE_SERIAL_TSLICE:
    kspf = w_serial_ks_i(filename);
    break;
  default:
    kspf = NULL;
  }

  return kspf;

}

/*---------------------------------------------------------------*/
/* reload a propagator in any of the formats, or cold propagator, or keep
   current propagator:
   FRESH, CONTINUE,
   RELOAD_ASCII, RELOAD_SERIAL, RELOAD_PARALLEL, RELOAD_MULTIDUMP
   */
int reload_ksprop( int flag, ks_prop_file *kspf, int color,
		   field_offset dest, int timing)
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
    FORALLSITES(i,s) clearvec( (su3_vector *)F_PT(s,dest));
    break;
  case RELOAD_ASCII:
    status = r_ascii_ks(kspf,color,dest);
    break;
  case RELOAD_SERIAL:
    status = r_serial_ks(kspf,color,dest); 
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
/* save a ksprop in any of the formats:
   FORGET,
   SAVE_ASCII, SAVE_SERIAL, SAVE_PARALLEL, SAVE_MULTIDUMP, SAVE_CHECKPOINT
   */
void save_ksprop( int flag, ks_prop_file *kspf, int color, 
		  field_offset src, int timing)
{
  double dtime;
  
  if(timing)dtime = -dclock();
  switch(flag){
  case FORGET:
    break;
  case SAVE_ASCII:
    w_ascii_ks(kspf,color,src);
    break;
  case SAVE_SERIAL:
    w_serial_ks(kspf,color,src);
    break;
  case SAVE_SERIAL_TSLICE:
    w_serial_ks(kspf,color,src);
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

/*---------------------------------------------------------------*/
void r_close_ksprop(int flag, ks_prop_file *kspf)
{
  
  switch(flag){
  case RELOAD_ASCII:
    r_ascii_ks_f(kspf);
    break;
  case RELOAD_SERIAL:
    r_serial_ks_f(kspf);
    break;
  default:
    node0_printf("r_close_ksprop: Unrecognized flag.\n");
    terminate(1);    
  }
}

/*---------------------------------------------------------------*/
void w_close_ksprop(int flag, ks_prop_file *kspf)
{
  switch(flag){
  case FORGET:
    break;
  case SAVE_ASCII:
    w_ascii_ks_f(kspf);
    break;
  case SAVE_SERIAL:
  case SAVE_SERIAL_FM:
  case SAVE_SERIAL_TSLICE:
    w_serial_ks_f(kspf); 
    break;
  default:
    node0_printf("w_close_ksprop: Unrecognized flag.\n");
    terminate(1);    
  }
}

/* find out what if any KS propagator should be loaded.
   This routine is only called by node 0.
*/
int ask_starting_ksprop( int prompt, int *flag, char *filename ){
    char savebuf[256];
    int status;

    if (prompt!=0) 
      printf( "enter 'fresh_ks', 'reload_ks_ascii', 'reload_ks_serial'\n");
    status=scanf("%s",savebuf);
    if(status !=1) {
        printf("ask_starting_ksprop: ERROR IN INPUT: starting prop command\n");
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
	    printf("ask_starting_ksprop: ERROR IN INPUT: file name read\n"); 
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
        "'forget_ks', 'save_ks_ascii', 'save_ks_serial', 'save_ks_serial_fm', save_ks_serial_tslice' ?\n");
    status=scanf("%s",savebuf);
    if(status !=1) {
        printf("ask_ending_ksprop: ERROR IN INPUT: ending ksprop command\n");
        return(1);
    }
    printf("%s ",savebuf);
    if(strcmp("save_ks_ascii",savebuf) == 0 )  {
        *flag=SAVE_ASCII;
    }
    else if(strcmp("save_ks_serial",savebuf) == 0 ) {
        *flag=SAVE_SERIAL;
    }
    else if(strcmp("save_ks_serial_fm",savebuf) == 0 ) {
        *flag=SAVE_SERIAL_FM;
    }
    else if(strcmp("save_ks_serial_tslice",savebuf) == 0 ) {
        *flag=SAVE_SERIAL_TSLICE;
    }
    else if(strcmp("forget_ks",savebuf) == 0 ) {
        *flag=FORGET;
	printf("\n");
    }
    else {
      printf("ask_ending_ksprop: ERROR IN INPUT: %s is not a valid command\n",
	     savebuf);
      return(1);
    }

    if( *flag != FORGET ){
        if(prompt!=0)printf("enter filename\n");
        status = scanf("%s",filename);
        if(status != 1){
    	    printf("ask_ending_ksprop: ERROR IN INPUT: save filename\n"); 
	    return(1);
        }
	printf("%s\n",filename);
    }
    return(0);

} /* end ask_ending_ksprop() */
