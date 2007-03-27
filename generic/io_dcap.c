/******************  io_dcap.c *****************************************/
/* I/O routines for the SU3 program. DCAP version */
/*   MIMD version 7. */

/* Wrappers for unformatted file I/O. */
/* This version uses DCAP to support dcache I/O */

#include "generic_includes.h"
#include <sys/types.h>
#include <errno.h>
#include <fcntl.h>
#include <unistd.h>
#include <dcap.h>
#include <string.h>

FILE *g_open(const char *filename, const char *mode)
{
  if(strstr(mode,"+") != NULL && strstr(mode,"w") != NULL){
    fprintf(stderr,"dc_fopen does not support appending to a file\n");
    return NULL;
  }
  return dc_fopen64(filename,mode);
}

int g_seek(FILE *stream, off_t offset, int whence)
{
  return dc_fseeko64(stream,offset,whence);
}

size_t g_write(const void *ptr, size_t size, size_t nmemb,FILE *stream)
{
  return dc_fwrite(ptr,size,nmemb,stream);
}

size_t g_read(void *ptr, size_t size, size_t nmemb, FILE *stream)
{
  return dc_fread(ptr,size,nmemb,stream);
}

int g_close(FILE *stream)
{
  return dc_fclose(stream);
}

