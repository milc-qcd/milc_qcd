/******************  io_paragon3.c *****************************************/
/* Parallel I/O routines for the SU3 program */
/*   MIMD version 6. */
/* NOTE: This file is NOT an upgraded version of io_paragon.c! */

/* Wrappers for parallel file access. */
/* This version calls Paragon I/O routines */

/* These are patterned after stdio fopen, fseek, fwrite, fread, fclose */
/* We need them because some systems (e.g. Intel Paragon) 
   use home-made routines, rather than using ANSI standard calls. */

/* Modifications

   9/11/97 Removed obsolete save_checkpoint reload_checkpoint C.D.
   8/8/97 Created by splitting from com_paragon.c C.D.
   */

#include "generic_includes.h"
#include <sys/types.h>
#include <errno.h>
#include <fcntl.h>

/**********************************************************************/
/* Wrappers for parallel file access. */
/* These are patterned after stdio fopen, fseek, fwrite, fread, fclose */
/* We need them because in systems (e.g. Intel Paragon)
   use home-made routines, rather than using ANSI standard calls.

   Note: Once a file is opened with g_open, it should not be accessed
   with ANSI routines, even though it pretends to return a FILE * pointer */

FILE *g_open(const char *filename, const char *mode)
{
  int fd, oflg;
  int *fp;

  /* Here we support only three ways to open the file:
     wb : create a new binary file for writing
     rb : open an existing binary file for reading
     rb+: open an existing binary file for reading and/or writing  */

  /* Decide on the mode for opening */

  if(mode[0] == 'a')
    {
      printf("g_open: Node %d. Append not supported\n",this_node);
      return NULL;
    }
  
  else if(mode[0] == 'w')
    oflg = O_WRONLY | O_CREAT;
  
  else if(mode[0] == 'r')
    {
      oflg = O_RDONLY;
      if(strchr(mode,'+') != NULL)
	oflg = O_RDWR;
    }

  else
    {
      printf("g_open: Node %d. mode %s not recognized\n",this_node,mode);
      return NULL;
    }

  /* Now open the file */

  if( (fd = gopen(filename, oflg, 0, 0644)) < 0)
    {
      printf("g_open: Node %d error %d opening %s\n",
	     this_node,errno,filename);
      return NULL;
    }
  
  /** Note, we may have to use this approach and distinguish between open for read and write: **/
  /**  fp = creat(filename, 0666); all nodes open file **/

  /* Set up a dummy file pointer structure */
  /* Since the structure FILE is system dependent, we use a pointer
     to an integer instead, and pretend that it is a FILE 
     No routines other than those in this package are supposed
     to look for "members" of "FILE" */

  fp = (int *)malloc(sizeof(int));
  if(fp == NULL)
    {
      printf("g_open: Node %d can't malloc fp\n",this_node);
      fflush(stdout); terminate(1);
    }

  *fp = fd;

  return (FILE *)fp;
}

int g_seek(FILE *stream, long int offset, int whence)
{
  int fd;
  fd = *((int *)stream);

  return lseek(fd, (off_t)offset, whence );
}

size_t g_write(const void *ptr, size_t size,size_t  nmemb, FILE *stream)
{
  int fd;
  fd = *((int *)stream);

  return write( fd, ptr, (int)(size*nmemb) )/size;
}

size_t g_read(const void *ptr, size_t size, size_t nmemb, FILE *stream)
{
  int fd;
  fd = *((int *)stream);

  return read( fd, ptr, (int)(size*nmemb) )/size;
}

int g_close(FILE *stream)
{
  int fd, status;
  fd = *((int *)stream);

  status = close(fd);
  free(stream);
  return status;
}
