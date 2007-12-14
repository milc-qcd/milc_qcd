/*********************** io_detect.c *************************/
/* MIMD version 7 */

/* Read the 32-bit magic number of a file and look it up in a table */
/* Returns the file type number from the table, which must be >= 0 */
/* Error return is -1 (not in table), -2 (open error) -3 (read error) */
/* Multidump files do not have a magic number, so can't be identified
   this way */
/* FNAL files have a common magic number, so it is necessary to
   read further to identify them. */

#include "generic_includes.h"
#include "../include/file_types.h"
#include <string.h>

int io_detect(char *filename, file_type ft[], int ntypes){
  FILE *fp;
  int i, status, words;
  int32type magic_no;
  int32type revmagic_no;
  char editfilename[513];

  /* Node 0 reads and checks */
  if(this_node == 0){
    fp = g_open(filename,"rb");
    if(fp == NULL){
      /* Special provision for partition or multifile format.  Try
	 adding the extension to the filename */
      strncpy(editfilename,filename,504);
      editfilename[504] = '\0';  /* Just in case of truncation */
      strcat(editfilename,".vol0000");
      fp = g_open(editfilename,"rb");
    }

    if(fp == NULL)status = -2;
    else
      {
	words = g_read(&magic_no, sizeof(int32type), 1, fp);
	g_close(fp);

	if(words != 1)status = -3;
	else
	  {
	    revmagic_no = magic_no;
	    byterevn(&revmagic_no, 1);

	    status = -1;
	    for(i = 0; i < ntypes; i++){
	      if(ft[i].magic_no == magic_no || 
		 ft[i].magic_no == revmagic_no)
		{
		  status = ft[i].type;
		  break;
		}
	    }
	  }
      }
  }

  /* Node 0 broadcasts the result */
  broadcast_bytes((char *)&status, sizeof(int));

  /* All nodes return the same value */
  return status;
}

/* For FNAL we base the detection on the number of elements per site */
int io_detect_fm(char *filename){
  FILE *fp;
  int status, words;
  int32type magic_no, revmagic_no, gmtime_stamp, size_of_element,
    elem_per_site;
  int byterevflag = 0;

  /* Node 0 reads and checks */
  if(this_node == 0){
    fp = g_open(filename,"rb");
    if(fp == NULL)status = -2;
    else
      {
	words = g_read(&magic_no, sizeof(int32type), 1, fp);

	if(words != 1)status = -3;
	else
	  {
	    revmagic_no = magic_no;
	    byterevn(&revmagic_no, 1);

	    status = -1;
	    if(revmagic_no == IO_UNI_MAGIC)byterevflag = 1;
	    if(magic_no == IO_UNI_MAGIC || revmagic_no == IO_UNI_MAGIC){
	      g_read(&gmtime_stamp, sizeof(gmtime_stamp), 1, fp);
	      g_read(&size_of_element, sizeof(int32type), 1, fp);
	      g_read(&elem_per_site, sizeof(int32type), 1, fp);

	      if(byterevflag){
		byterevn(&size_of_element, 1);
		byterevn(&elem_per_site, 1);
	      }

	      if(size_of_element != sizeof(float))status = -1;
	      else 
		if(elem_per_site == sizeof(fsu3_matrix)/sizeof(float))
		  status = FILE_TYPE_KSFMPROP;
	      else 
		if(elem_per_site == sizeof(fwilson_propagator)/sizeof(float) ||
		   elem_per_site == sizeof(fwilson_vector)/sizeof(float))
		  status = FILE_TYPE_W_FMPROP;
	      else 
		if(elem_per_site == 4*sizeof(fsu3_matrix)/sizeof(float))
		  status = FILE_TYPE_GAUGE_FNAL;
	    }
	  }
      }
    g_close(fp);
  }
  /* Node 0 broadcasts the result */
  broadcast_bytes((char *)&status, sizeof(int));

  /* All nodes return the same value */
  return status;
}
