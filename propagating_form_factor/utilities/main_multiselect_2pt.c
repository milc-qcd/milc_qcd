/*
   Read a set of two point functions from a binary correlator file
   and write it to specified files

*/

#include "read_hl_form.h"
#include "prop_form_utilities.h"


int main(int argc, char *argv[])
{
  char *twopt_file ;
  int nselect;
  int nfile;
  twopt_oneselect fileparam[MAX_NO_FILE];

  /***--------------------------------------------------*****/

  if( argc != 2 )
  {
    printf("usage::  %s  [2pt file] \n",argv[0]);
    exit(1);
  }

  twopt_file = argv[1];

  read_multiselect_twopt_param(&nfile,fileparam);

  read_multiselect_twopt(twopt_file, nfile, fileparam);
  
  return 0 ;
}
