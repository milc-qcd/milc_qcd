/********************** io_helpers_ks_utils.c ***************************/
/* MIMD version 7 */

#include "ks_imp_includes.h"
#include "../include/io_lat.h"

/* find out the name of the input or output color vector source
   file. This routine is called only by node 0.
*/

int ask_color_vector( int prompt, int *flag, char *filename ){
  char savebuf[256];
  int status;
  
  if (prompt==1) 
    printf("enter 'fresh', 'reload_serial_ks_vector' 'save_serial_ks_vector' 'save_partfile_scidac_ks_vector'\n");
  status=scanf("%s",savebuf);
  if(status !=1) {
    printf("\nask_color_vector: ERROR IN INPUT: command \"%s\" is invalid\n",savebuf);
    return(1);
  }

  printf("%s ",savebuf);
  if(strcmp("fresh",savebuf) == 0 ){
       *flag = FRESH;
       printf("\n");
  }
  else if(strcmp("forget",savebuf) == 0 ){
       *flag = FORGET;
       printf("\n");
  }
  else if(strcmp("reload_serial_ks_vector",savebuf) == 0 ) {
    *flag = RELOAD_SERIAL;
  }
  else if(strcmp("save_serial_ks_vector",savebuf) == 0 ) {
    *flag = SAVE_SERIAL;
  }
  else if(strcmp("save_partfile_scidac_ks_vector",savebuf) == 0 ) {
    *flag = SAVE_PARTFILE_SCIDAC;
  }
  else{
    printf("\nask_color_vector: ERROR IN INPUT: color vector command \"%s\" is invalid\n",savebuf); return(1);
  }

  /*read name of file and load it */
  if( *flag != FRESH && *flag != FORGET ){
    if(prompt==1)printf("enter name of file containing the color vector\n");
    status=scanf("%s",filename);
    if(status !=1) {
      printf("\nask_color_vector: ERROR IN INPUT: error reading file name\n"); return(1);
    }
    printf("%s\n",filename);
  }
  return(0);
}

int ask_color_matrix( int prompt, int *flag, char *filename ){
  char savebuf[256];
  int status;
  
  if (prompt==1) 
    printf("enter 'fresh', 'reload_serial_color_matrix' 'save_serial_color_matrix' 'save_partfile_scidac_color_matrix'\n");
  status=scanf("%s",savebuf);
  if(status !=1) {
    printf("\nask_color_matrix: ERROR IN INPUT: command \"%s\" is invalid\n",savebuf);
    return(1);
  }

  printf("%s ",savebuf);
  if(strcmp("fresh",savebuf) == 0 ){
       *flag = FRESH;
       printf("\n");
  }
  else if(strcmp("forget",savebuf) == 0 ){
       *flag = FORGET;
       printf("\n");
  }
  else if(strcmp("reload_serial_color_matrix",savebuf) == 0 ) {
    *flag = RELOAD_SERIAL;
  }
  else if(strcmp("save_serial_color_matrix",savebuf) == 0 ) {
    *flag = SAVE_SERIAL;
  }
  else if(strcmp("save_partfile_scidac_color_matrix",savebuf) == 0 ) {
    *flag = SAVE_PARTFILE_SCIDAC;
  }
  else{
    printf("\nask_color_matrix: ERROR IN INPUT: color matrix command \"%s\" is invalid\n",savebuf); return(1);
  }

  /*read name of file and load it */
  if( *flag != FRESH && *flag != FORGET ){
    if(prompt==1)printf("enter name of file containing the color matrix\n");
    status=scanf("%s",filename);
    if(status !=1) {
      printf("\nask_color_matrix: ERROR IN INPUT: error reading file name\n"); return(1);
    }
    printf("%s\n",filename);
  }
  return(0);
}

