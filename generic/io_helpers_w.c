/********************** io_helpers_w.c **********************************/
/* MIMD version 6 */
/* CD 10/97 from io_helpers.c DT 8/97 
   General purpose high level propagator I/O routines, 
   to be used by any application that wants them.
   */

#include "generic_includes.h"
#include "../include/io_lat.h"
#include "../include/io_wb.h"


/*---------------------------------------------------------------*/
/* Open propagator file for reading */
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
    /* Read header and close temporarily */
    wpf = r_parallel_w_i(filename);
    r_parallel_w_c(wpf);
    break;
  case RELOAD_MULTIDUMP:
    /* Read header and close temporarily */
    wpf = r_multidump_w_i(filename);
    r_multidump_w_c(wpf);
    break;
  default:
    wpf = NULL;
  }
  return wpf;
}

/*---------------------------------------------------------------*/
/* Open propagator file for writing */
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
    /* Write header and close temporarily */
    wpf = w_parallel_w_i(filename);
    w_parallel_w_c(wpf);
    break;
  case SAVE_MULTIDUMP:
    /* Write header and close temporarily */
    wpf = w_multidump_w_i(filename);
    w_multidump_w_c(wpf);
    break;
  case SAVE_CHECKPOINT:
    /* Write header and close temporarily */
    wpf = w_checkpoint_w_i(filename);
    w_checkpoint_w_c(wpf);
    break;
  default:
    wpf = NULL;
  }
  return wpf;
}

/*---------------------------------------------------------------*/
/* reload a propagator in any of the formats, or cold propagator, or keep
   current propagator:
   FRESH, CONTINUE,
   RELOAD_ASCII, RELOAD_SERIAL, RELOAD_PARALLEL, RELOAD_MULTIDUMP
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
/* save a propagator in any of the formats:
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
    printf("save wilson propagator:\n enter 'forget_prop', 'save_ascii_prop', 'save_serial_prop' , 'save_parallel_prop', 'save_multidump_prop', or 'save_checkpoint_prop'\n");
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
