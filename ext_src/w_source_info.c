/*********************** w_source_info.c *************************/
/* MIMD version 7 */

/* For ext_src */

/* Create the file XML for an extended clover source file */

#include "ext_src_includes.h"
#include <string.h>

#define INFOSTRING_MAX 2048


char *create_ws_XML(char *filename, quark_source *wqs){

  size_t bytes = 0;
  char *info = (char *)malloc(INFOSTRING_MAX);
  size_t max = INFOSTRING_MAX;
  //  char begin[] = "<?xml version=\"1.0\" encoding=\"UTF-8\"?><info>";
  //  char end[] = "</info>";

  //  snprintf(info+bytes, max-bytes,"%s",begin);
  //  bytes = strlen(info);
  
  sprint_w_prop_info_item(info+bytes, max-bytes,"quark.description","%s",
			  "Extended clover source",0,0);

  bytes = strlen(info);
  sprint_w_prop_info_item(info+bytes, max-bytes,"source.filename","%s",
			  filename,0,0);

  bytes = strlen(info);
  sprint_w_prop_info_item(info+bytes, max-bytes,"source.t",
			  "%d",(char *)&wqs->t0,0,0);
  
  // bytes = strlen(info);
  //  snprintf(info+bytes, max-bytes,"%s",end);

  return info;
}

void free_ws_XML(char *info){
  if(info != NULL)free(info);
}

