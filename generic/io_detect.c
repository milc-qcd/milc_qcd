/*********************** io_detect.c *************************/
/* MIMD version 7 */

/* Read the 32-bit magic number of a file and look it up in a table */
/* Returns the file type number from the table, which must be >= 0 */
/* Error return is -1 (not in table), -2 (open error) -3 (read error) */

#include "generic_includes.h"

int io_detect(char *filename, file_type ft[], int ntypes){
  FILE *fp;
  int i, status;
  int32type magic_no;
  int32type revmagic_no;
  
  fp = fopen(filename,"r");
  if(fp == NULL)return -2;

  status = fread(&magic_no, sizeof(int32type), 1, fp);
  fclose(fp);

  if(status != 1)return -3;

  revmagic_no = magic_no;
  byterevn(&revmagic_no, 1);

  for(i = 0; i < ntypes; i++)
    if(ft[i].magic_no == magic_no || ft[i].magic_no == revmagic_no)
      return ft[i].type;

  return -1;
}
