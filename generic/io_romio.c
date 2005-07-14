/******************  io_romio.c *****************************************/
/* Parallel I/O routines for the SU3 program */
/*   MIMD version 7. */

/* Wrappers for parallel file access. */
/* This version uses MPI-2 standard routines */

/* These interface routines are patterned after stdio fopen, fseek,
   fwrite, fread, fclose */

/* NOTE: Once a file is opened with g_open, it should not be accessed
   with ANSI routines, even though it pretends to return a FILE * pointer */

#include "generic_includes.h"
#include <sys/types.h>
#include <errno.h>
#include <fcntl.h>
#include <string.h>
#include <mpi.h>
/*#define IO_ROMIO_DEBUG */
#ifdef IO_ROMIO_DEBUG
#define DEBUG 1
#else
#define DEBUG 0
#endif

FILE *g_open(const char *filename, const char *mode)
{
  int amode;
  int status;
  MPI_File *mpifile;
  char datarep[] = "native";
  int debug = DEBUG;

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
    {
      amode = MPI_MODE_WRONLY | MPI_MODE_CREATE;
    }
  
  else if(mode[0] == 'r')
    {
      amode = MPI_MODE_RDONLY;
      if(strchr(mode,'+') != NULL)
	amode = MPI_MODE_RDWR;
    }

  else
    {
      printf("g_open: Node %d. mode %s not recognized\n",this_node,mode);
      return NULL;
    }

  mpifile = (MPI_File *)malloc(sizeof(MPI_File));

  status = MPI_File_open(MPI_COMM_WORLD, (char *)filename, amode, MPI_INFO_NULL, mpifile);

  if(debug)printf("g_open: Node %d. Return from MPI_File_open %d\n",this_node,status);
  status = MPI_File_set_view(*mpifile, 0, MPI_BYTE,
		    MPI_BYTE, datarep, MPI_INFO_NULL);
  if(debug)printf("g_open: Node %d. Return from MPI_File_set_view %d\n",this_node,status);
  return (FILE *)mpifile;
}

int g_seek(FILE *stream, off_t offset, int whence)
{
  MPI_Offset seek = offset, current;
  int status;
  int mpi_whence;
  int debug = DEBUG;

  /* Interpret "whence" code */

  if(whence == SEEK_CUR){
    mpi_whence = MPI_SEEK_CUR;
  }
  else if(whence == SEEK_SET){
    mpi_whence = MPI_SEEK_SET;
  }
  else if(whence == SEEK_END){
    mpi_whence = MPI_SEEK_END;
  }
  else{
      printf("g_seek: Node %d. whence code %d not supported\n",this_node,whence);
      return 1;
  }

  if(debug)printf("g_seek: Node %d. Calling MPI_File_seek offset %d\n",this_node,offset);
  status = MPI_File_seek(*((MPI_File *)stream), 
				(MPI_Offset) offset, mpi_whence);
  if(debug)printf("g_seek: Node %d. Return from MPI_File_seek %d\n",this_node,status);
  return status;
}

size_t g_write(const void *ptr, size_t size, size_t nmemb,FILE *stream)
{
  MPI_Status status;
  int debug = DEBUG;

  if(debug)printf("g_write: Node %d calling MPI_File_write for %d bytes\n",
	 this_node,(int)(size*nmemb));
  if(MPI_File_write(*((MPI_File *)stream), (void *)ptr, size*nmemb, 
		    MPI_BYTE, &status) == 0)
    return (size_t)nmemb;
  else
    return 0;
}

size_t g_read(void *ptr, size_t size, size_t nmemb, FILE *stream)
{
  MPI_Status status;
  int count;
  int debug = DEBUG;

  if(debug)printf("g_read: Node %d calling MPI_File_read\n",this_node);

  MPI_File_read(*((MPI_File *)stream), ptr, size*nmemb, MPI_BYTE, &status);

  MPI_Get_count(&status, MPI_BYTE, &count);

  if(debug)printf("g_read: Node %d size %d nmemb %d count %d datatype %d\n",this_node,
	 (int)size, (int)nmemb, count, MPI_BYTE);

  /*  return (size_t)(count/size);*/
  return (size_t)nmemb;
}

int g_close(FILE *stream)
{
  int debug = DEBUG;

  if(debug)printf("g_close: Node %d\n",this_node);
  return MPI_File_close((MPI_File *)stream);
}

