/*
   Read a set of three point functions from a binary correlator file
   and write it to specified files

*/

#include "read_hl_form.h"
#include "prop_form_utilities.h"


int main(int argc, char *argv[])
{
  char *form_file ;
  int nselect;
  int nfile;
  threept_oneselect fileparam[MAX_NO_FILE];

  if( argc != 2 )
  {
    printf("usage::  %s  [propagating form factor 3pt file] \n",argv[0]);
    exit(1);
  }
  form_file = argv[1];

  read_multiselect_form_param(&nfile,fileparam);

  /*** read and dump the three point functions to the specified files  ******/
  read_multiselect_form_corr(form_file, nfile, fileparam);

  return 0 ;
}
