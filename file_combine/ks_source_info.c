/*********************** ks_source_info.c *************************/
/* MIMD version 7 */

/* For file_combine */

/* Read and write file and record info for KS propagator source file formats */

#include <qio_config.h>
#include <stdio.h>
#include <string.h>
#include <qio.h>
#include <qioxml.h>
#include <qio_string.h>
#include <qio_stdint.h>
#include <sys/types.h>
#include <time.h>
#include <stdlib.h>
#include "qioxml_usqcd_ksprop_source.h"


int QIO_insert_usqcdkspropsourcefile_info_tag_string(QIO_USQCDKSPropSourceFileInfoWrapper *wrapper,
						     char *fileinfo_tags){
  wrapper->usqcdkspropsourcefileinfo_tags.occur = 0;
  if(!fileinfo_tags)return QIO_BAD_ARG;
  strncpy(wrapper->usqcdkspropsourcefileinfo_tags.value, fileinfo_tags,
          QIO_MAXVALUESTRING-1);
  wrapper->usqcdkspropsourcefileinfo_tags.value[QIO_MAXVALUESTRING-1] = '\0';
  wrapper->usqcdkspropsourcefileinfo_tags.occur = 1;
  if(strlen(fileinfo_tags) >= QIO_MAXVALUESTRING)return QIO_ERR_ALLOC;
  else return QIO_SUCCESS;
}

char *QIO_get_usqcd_kspropsourcefile_info_tag_string(QIO_USQCDKSPropSourceFileInfoWrapper *wrapper){
  return wrapper->usqcdkspropsourcefileinfo_tags.value;
}

void QIO_encode_usqcd_kspropsourcefile_info(QIO_String *file_string, 
				      QIO_USQCDKSPropSourceFileInfo *file_info)
{
  char *buf;
  int remainder,n;
  char fileinfo_tags[QIO_MAXVALUESTRING];
  QIO_USQCDKSPropSourceFileInfoWrapper wrapper = QIO_USQCD_KSPROPSOURCEFILE_INFO_WRAPPER;

  /* Start by creating string of inner tags */
  buf = fileinfo_tags;
  remainder = QIO_MAXVALUESTRING;

  /* Build inner tag string by appending tags */
  *buf = '\0';
  buf = QIO_encode_as_string(buf,&file_info->version, &remainder);
  buf = QIO_encode_as_string(buf,&file_info->type, &remainder);
  buf = QIO_encode_as_string(buf,&file_info->info, &remainder);

  /* Insert inner tag string into file wrapper structure */
  QIO_insert_usqcdkspropsourcefile_info_tag_string(&wrapper, fileinfo_tags);

  /* Now build final XML string */
  QIO_string_realloc(file_string, QIO_STRINGALLOC);
  buf  = QIO_string_ptr(file_string);
  remainder = QIO_string_length(file_string);

  /* Begin with xml info stuff */
  strncpy(buf,QIO_XMLINFO,remainder);
  buf[remainder-1] = '\0';
  n = strlen(buf);
  remainder -= n;
  buf += n;
  if(remainder < 0){
    printf("QIO_encode_usqcd_kspropsourcefile_info: file_string overflow\n");
  }
  else{
    /* Conclude by appending the wrapped tag string */
    buf = QIO_encode_as_string (buf,&wrapper.usqcdkspropsourcefileinfo_tags, &remainder);
  }
}

int QIO_decode_usqcd_kspropsourcefile_info(QIO_USQCDKSPropSourceFileInfo *file_info,
				   QIO_String *file_string)
{
  char *parse_pt = QIO_string_ptr(file_string);
  char *tmp_pt;
  char tag[QIO_MAXTAG];
  char tags_string[QIO_MAXVALUESTRING];
  char value_string[QIO_MAXVALUESTRING];
  int errors = 0;
  QIO_USQCDKSPropSourceFileInfoWrapper wrapper = QIO_USQCD_KSPROPSOURCEFILE_INFO_WRAPPER;
  QIO_USQCDKSPropSourceFileInfo templ = QIO_USQCD_KSPROPSOURCEFILE_INFO_TEMPLATE;
  char *left_angle;

  /* Initialize file info structure from a template */
  memcpy(file_info, &templ, sizeof(QIO_USQCDKSPropSourceFileInfo));
  
  /* Start parsing file_string */
  /* Check leading tag, which is probably the info phrase "<?xml ...?>" */
  /* We ignore it if it is there */
  tmp_pt = QIO_next_tag(parse_pt, tag, &left_angle);
  if(strcmp(tag,QIO_QUESTXML)==0){
    /* Found ?xml, so resume parsing after the closing ">", ignoring
       the field. Otherwise, leave the parse_pt at its initial value */
    parse_pt = tmp_pt;
  }

  /* Open top-level tag (wrapper) and extract string containing tags */
  parse_pt = QIO_get_tag_value(parse_pt, tag, tags_string);
  QIO_decode_as_string (tag, tags_string, &wrapper.usqcdkspropsourcefileinfo_tags);

  /* If outer wrapper has bad tag, exit with error status */
  if(QIO_check_string_occur(&wrapper.usqcdkspropsourcefileinfo_tags))
    return QIO_BAD_XML;
  /* Otherwise start parsing the string of tags */
  parse_pt = QIO_get_usqcd_kspropsourcefile_info_tag_string(&wrapper);
  /* Scan string until null character is reached */
  while(*parse_pt){
    parse_pt = QIO_get_tag_value(parse_pt, tag, value_string);

    QIO_decode_as_string(tag,value_string,&file_info->version);
    QIO_decode_as_string(tag,value_string,&file_info->type);
    QIO_decode_as_string(tag,value_string,&file_info->info);
  }

  /* Check for completeness */

  errors += QIO_check_string_occur(&file_info->version);
  errors += QIO_check_string_occur(&file_info->type);
  errors += QIO_check_string_occur(&file_info->info);

  return errors;
}

/********************************************************************/
/* Support for USQCD KS propagator source file info */

/* Accessors */

/* Return integer code or negative value for failure */
int QIO_get_usqcd_kspropsourcefile_type(QIO_USQCDKSPropSourceFileInfo *file_info)
{
  char *string = file_info->type.value;

  if(strcmp(string,QIO_USQCDKSPROPSOURCEFILETYPESTRING) == 0)
    return QIO_USQCDKSPROPSOURCEFILETYPE;
  else
    return QIO_ERR_FILE_INFO;
}

char *QIO_get_usqcd_kspropsourcefile_info(QIO_USQCDKSPropSourceFileInfo *file_info)
{
  return file_info->info.value;
}

int QIO_defined_usqcd_kspropsourcefile_type(QIO_USQCDKSPropSourceFileInfo *file_info)
{
  return file_info->type.occur;
}

int QIO_defined_usqcd_kspropsourcefile_info(QIO_USQCDKSPropSourceFileInfo *file_info)
{
  return file_info->info.occur;
}

int QIO_insert_usqcdkspropsourcefile_version(QIO_USQCDKSPropSourceFileInfo *file_info, char *version)
{
  file_info->version.occur = 0;
  if(!version)return QIO_BAD_ARG;
  strncpy(file_info->version.value, version, QIO_MAXVALUESTRING-1);
  file_info->version.value[QIO_MAXVALUESTRING-1] = '\0';
  file_info->version.occur = 1;
  if(strlen(version) >= QIO_MAXVALUESTRING)return QIO_ERR_ALLOC;
  else return QIO_SUCCESS;
}

/* Takes one of the integer type codes and translates it to the string */
int QIO_insert_usqcdkspropsourcefile_type( QIO_USQCDKSPropSourceFileInfo *file_info)
{
  file_info->type.occur = 0;
  if(!file_info)return QIO_BAD_ARG;
  strncpy(file_info->type.value, 
	  QIO_USQCDKSPROPSOURCEFILETYPESTRING, QIO_MAXVALUESTRING-1);
  file_info->type.value[QIO_MAXVALUESTRING-1] = '\0';
  file_info->type.occur = 1;
  return QIO_SUCCESS;
}

int QIO_insert_usqcdkspropsourcefile_info( QIO_USQCDKSPropSourceFileInfo *file_info, char *info)
{
  file_info->info.occur = 0;
  if(!file_info)return QIO_BAD_ARG;
  strncpy(file_info->info.value, info, QIO_MAXVALUESTRING-1);
  file_info->info.value[QIO_MAXVALUESTRING-1] = '\0';
  file_info->info.occur = 1;
  if(strlen(info) >= QIO_MAXVALUESTRING)return QIO_ERR_ALLOC;
  else return QIO_SUCCESS;
}

QIO_USQCDKSPropSourceFileInfo *QIO_create_usqcd_kspropsourcefile_info(char *info)
{
  
  QIO_USQCDKSPropSourceFileInfo templ = QIO_USQCD_KSPROPSOURCEFILE_INFO_TEMPLATE;
  QIO_USQCDKSPropSourceFileInfo *file_info;

  file_info = (QIO_USQCDKSPropSourceFileInfo *)malloc(sizeof(QIO_USQCDKSPropSourceFileInfo));
  if(!file_info)return NULL;

  memcpy(file_info, &templ, sizeof(QIO_USQCDKSPropSourceFileInfo));
  QIO_insert_usqcdkspropsourcefile_version(file_info,QIO_USQCDKSPROPSOURCEFILEFORMATVERSION);  
  
  QIO_insert_usqcdkspropsourcefile_type( file_info );
  QIO_insert_usqcdkspropsourcefile_info( file_info, info );

  return file_info;
}

#include "ext_src_includes.h"

#define MAX_XML 2049

char *create_kss_XML(char *filename, quark_source *ksqs)
{
  char *xml;

  xml = (char *)malloc(MAX_XML);
  
  snprintf(xml,MAX_XML,"\nsource.description Extended KS source\nquark.filename %s\nsource.t %d\n\n",
	   filename, ksqs->t0);

  xml[MAX_XML-1] = '\0';
  return xml;
}


void free_kss_XML(char *xml){
  if(xml != NULL)free(xml);
}
