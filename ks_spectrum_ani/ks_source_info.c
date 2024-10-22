/*********************** ks_source_info.c *************************/
/* MIMD version 7 */

/* For ks_spectrum */

/* Create the file XML for an extended clover source file */

/* NEEDS WORK: PLEASE REWRITE THIS TO CONFORM TO ks_info.c */

#include "ks_spectrum_includes.h"

#define MAX_XML 2049

char *create_kss_XML(char *filename, quark_source *ksqs)
{
  char *xml;

  xml = (char *)malloc(MAX_XML);
  
  snprintf(xml,MAX_XML,"\nsource.description KS source\nquark.filename %s\nsource.t %d\n\n",
	   filename, ksqs->t0);

  xml[MAX_XML-1] = '\0';
  return xml;
}


void free_kss_XML(char *xml){
  if(xml != NULL)free(xml);
}
