/****************** io_piofs.c *****************************************/
/* MIMD version 7 */
/* Parallel I/O routines for the SU3 program */
/* NOT MAINTAINED.  TEST BEFORE USE. */
/*   MIMD version 7. */

/* Wrappers for parallel file access. */
/* This version uses IBM/PIOFS for the SP2 */

/* Wrappers for parallel file access. */
/* These are patterned after stdio fopen, fseek, fwrite, fread, fclose */
/* We need them because some systems (e.g. Intel Paragon and SP2) 
   use home-made routines, rather than using ANSI standard calls.

   NOTE: Once a file is opened with g_open, it should not be accessed
   with ANSI routines, even though it pretends to return a FILE * pointer */

#include "generic_includes.h"
#include <sys/types.h>
#include <errno.h>
#include <piofs/piofs_ioctl.h>   /* For IBM SP2 PIOFS */
#include <fcntl.h>

FILE *g_open(const char *filename, const char *mode)
{
  int fd,oflg;
  int *fp;
  piofs_create_t mycreate;

  /* Here we support only three ways to open the file:
     wb : create a new binary file for writing
     rb : open an existing binary file for reading
     rb+: open an existing binary file for reading and/or writing  */

  /* Decide on the mode for opening */

  if(mode[0] == 'a')
    {
      printf("g_open: Node %d. Append not supported in PIOFS\n",this_node);
      return NULL;
    }
  
  else if(mode[0] == 'w')
    {
     if(this_node==0)
       {
	 /* Node 0 creates a new parallel file */
	 strcpy(mycreate.name, filename);
	 mycreate.bsu          = 4096;
	 mycreate.cells        = 8;
	 mycreate.permissions  = 0644;
	 mycreate.flags        = REPLACE;
	 mycreate.base_node    = -1;
	 if(piofsioctl(0,PIOFS_CREATE, &mycreate) != 0)
	   {
	     printf("g_open: Node %d Error %d in PIOFS_CREATE\n",
		    this_node,errno);
	     return NULL;
	   }
       }

     /* Nodes must wait for file creation before proceeding to open it */
     g_sync();
     oflg = O_WRONLY | O_CREAT;
   }
  
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

  if( (fd = open(filename, oflg)) == -1)
    {
      printf("g_open: Node %d error %d opening %s\n",
	     this_node,errno,filename);
      return NULL;
    }
  
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

size_t g_write(const void *ptr, size_t size, size_t nmemb,FILE *stream)
{
  int fd;
  fd = *((int *)stream);

  return write( fd, ptr, (int)(size*nmemb) )/size;
}

size_t g_read(void *ptr, size_t size, size_t nmemb, FILE *stream)
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
