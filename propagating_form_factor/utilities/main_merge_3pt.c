/*
   Read three point functions from two files and merge
   sequential quark lists to one (with selections hard coded below)
   The two files must be identical except for sequential quark
   lists.

*/

#include "read_hl_form.h"
#include "prop_form_utilities.h"


int main(int argc, char *argv[])
{
  char *input1_file, *input2_file, *output_file ;
    
  int32type  *p_momentum  ;
  int32type  *q_momentum  ;
  complex  *corr ;
  int *corr_oper_list;
  int *corr_copy_list;
  int nt ;
  int hl_flag ;
  int no_q_values ;
  int no_p_values ;
  int no_oper ;
  int no_zonked ;
  int no_spectator ;
  int no_sequential ; 
  int no_copies ;

  /***--------------------------------------------------*****/

  printf("=============================================================\n");
  printf(" Merging three point functions (propagating form factor code)  \n");
  printf("=============================================================\n");

  if( argc != 4 )
  {
    printf("usage::  %s  [1st input file] [2nd input file] [output file]\n",argv[0]);
    exit(1);
  }

  input1_file = argv[1];
  input2_file = argv[2];
  output_file = argv[3];

  /*** read and reserve memory for the first three point function file ******/
  read_merge1_corr(&corr, &corr_oper_list, &corr_copy_list,
    &p_momentum, &q_momentum, 
    &nt, &no_p_values, &no_q_values, &no_oper,
    &no_spectator, &no_sequential, &no_zonked,    
    &hl_flag, input1_file, &no_copies );

  /*** write out some titles *****/
  titles_hl(nt, no_p_values, no_q_values, no_oper ,p_momentum, q_momentum , hl_flag) ;



  /*** read the second three point function file ******/
    read_merge2_corr(&corr, &corr_oper_list, &corr_copy_list,
	&p_momentum, &q_momentum, 
	&nt, &no_p_values, &no_q_values, &no_oper,
	&no_spectator, &no_sequential, &no_zonked,    
	&hl_flag, input2_file, &no_copies );
	
  /*** write out some titles *****/
  titles_hl(nt, no_p_values, no_q_values, no_oper ,p_momentum, q_momentum , hl_flag) ;



  /*** write the merged three point functions to disk ****/

    write_3pt(corr,output_file,q_momentum,p_momentum,hl_flag,no_zonked,
	no_spectator,no_sequential,no_p_values,no_q_values,
	no_oper,nt,no_copies);
	
  return 0 ;
}
