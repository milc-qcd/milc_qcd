/*********************** io_detect.c *************************/
/* MIMD version 7 */

/* Read the 32-bit magic number of a file and look it up in a table */
/* Returns the file type number from the table, which must be >= 0 */
/* Error return is -1 (not in table), -2 (open error) -3 (read error) */
/* Multidump files do not have a magic number, so can't be identified
   this way */

#include "generic_includes.h"
#include <string.h>

int io_detect(char *filename, file_type ft[], int ntypes){
  FILE *fp;
  int i, status, words;
  int32type magic_no;
  int32type revmagic_no;
  char editfilename[513];

  /* Node 0 reads and checks */
  if(this_node == 0){
    fp = fopen(filename,"rb");
    if(fp == NULL){
      /* Special provision for partition or multifile format.  Try
	 adding the extension to the filename */
      strncpy(editfilename,filename,504);
      editfilename[504] = '\0';  /* Just in case of truncation */
      strcat(editfilename,".vol0000");
      fp = fopen(filename,"rb");
    }

    if(fp == NULL)status = -2;
    else
      {
	words = fread(&magic_no, sizeof(int32type), 1, fp);
	fclose(fp);

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
