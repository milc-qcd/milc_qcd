/********************** remap_stdio_from_args.c ***************************/
/* MIMD version 7 */
/* CD 11/04
   AB 2/11 Add append option

   In case input and output redirection is not an option

   Remap stdin, stdout, and stderr from command line arguments.  This
   routine must be called only after all preceding command line
   arguments have been consumed and removed.

   CD 10/09 This version supports multiple jobs with separate stdin
   stdout and stderr for each.

   Usage:

   cmd [stdin [stdout [stderr]]]

*/

//AB Temporarily put the define here, later move to Makefile
//   Allow to append to stdout, stderr files, if defined
//#define REMAP_STDIO_APPEND

#include "generic_includes.h"
#include <string.h>

int remap_stdio_from_args(int argc, char *argv[]){
  FILE *fp;
  int num_jobs = numjobs();
  int jobid = myjobid();

  /* stdin is remapped only on node 0 on any machine */
  if(argc > 1 && mynode() == 0){
    char *newfile = (char *)malloc(strlen(argv[1])+16);
    if(num_jobs == 1){
      sprintf(newfile,"%s",argv[1]);
    } else {
      sprintf(newfile,"%s.j%02d",argv[1],jobid);
#ifdef JOB_DEBUG
      printf("(%d:%d) input file is %s\n",jobid,mynode(),newfile);
#endif
    }
    fp = freopen(newfile,"r",stdin);
    if(fp == NULL){
      node0_printf("Can't open stdin file %s for reading.\n",newfile);
      return 1;
    }
    free(newfile);
  }

  if(argc > 2){
    char *newfile = (char *)malloc(strlen(argv[2])+16);
    if(num_jobs == 1){
      sprintf(newfile,"%s",argv[2]);
    } else {
      sprintf(newfile,"%s.j%02d",argv[2],jobid);
#ifdef JOB_DEBUG
      printf("(%d:%d) output file is %s\n",jobid,mynode(),newfile);
#endif
    }

#ifdef REMAP_STDIO_APPEND
    fp = freopen(newfile,"a",stdout);
#else
    fp = freopen(newfile,"w",stdout);
#endif
    if(fp == NULL){
      node0_printf("Can't open stdout file %s for writing\n",argv[2]);
      return 1;
    }
    free(newfile);
  }
  
  if(argc > 3){
    char *newfile = (char *)malloc(strlen(argv[2])+16);
    if(num_jobs == 1){
      sprintf(newfile,"%s",argv[3]);
    } else {
      sprintf(newfile,"%s.j%02d",argv[3],jobid);
#ifdef JOB_DEBUG
      printf("(%d:%d) error file is %s\n",jobid,mynode(),newfile);
#endif
    }
#ifdef REMAP_STDIO_APPEND
    fp = freopen(newfile,"a",stderr);
#else
    fp = freopen(newfile,"w",stderr);
#endif
    if(fp == NULL){
      node0_printf("Can't open stderr file %s for writing\n",newfile);
      return 1;
    }
    free(newfile);
  }
  return 0;
}
