/******************  io_ansi.c *****************************************/
/* Parallel I/O routines for the SU3 program */
/*   MIMD version 7. */

/* Wrappers for parallel file access. */
/* This version uses ANSI standard I/O which works for the T3D 
   and, of course, the scalar vanilla code */
/* FOR THE T3E USE io_nonansi.c! */

/* These are patterned after stdio fopen, fseek, fwrite, fread, fclose */
/* We need them because some systems (e.g. Intel Paragon) 
   use home-made routines, rather than using ANSI standard calls. */

#include "generic_includes.h"
#include <sys/types.h>
#include <errno.h>
#include <fcntl.h>

FILE *g_open(const char *filename, const char *mode)
{
  return fopen(filename,mode);
}

int g_seek(FILE *stream, off_t offset, int whence)
{
#if HAVE_FSEEKO
  return fseeko(stream,offset,whence);
#else
  return fseek(stream,(long)offset,whence);
#endif
}

size_t g_write(const void *ptr, size_t size, size_t nmemb,FILE *stream)
{
  return fwrite(ptr,size,nmemb,stream);
}

size_t g_read(void *ptr, size_t size, size_t nmemb, FILE *stream)
{
  return fread(ptr,size,nmemb,stream);
}

int g_close(FILE *stream)
{
  return fclose(stream);
}

