/*
   Read the three point function from a file with two or three source momenta
   and split it into separate files according to source momentum

*/

#include "read_hl_form.h"
#include "prop_form_utilities.h"


int main(int argc, char *argv[])
{
  char *input_3pt_file, *output1_file, *output2_file, *output3_file ;

  int32type  *p_momentum  ;
  int32type  *q_momentum  ;
  complex  *corr ;
  int *corr_oper_list;
  int *corr_copy_list;
  int nt ;
  int hl_flag ;
  int no_q_mom ;
  int no_p_mom ;
  int no_oper ;

  int no_output_files ;

  int no_zonked ;
  int no_spectator ;

  int no_sequential ; 
  int nocopies ;
  int p_mom_select;

  /***--------------------------------------------------*****/

  printf("=============================================================\n");
  printf(" Splitting the three point function (propagating form factor code)  \n");
  printf("=============================================================\n");

  if( argc < 4 )
  {
    printf("usage::  %s  <input 3pt file> <output file for 1st p> <output file for 2nd p> [< output file for 3rd p>]\n",argv[0]);
    exit(1);
  }

  no_output_files = argc - 2;
  input_3pt_file = argv[1];
  output1_file = argv[2];
  output2_file = argv[3];
  output3_file = argv[4];

  /*** read and reserve memory for the three point function ******/
  read_propagating_form_corr(&corr,&corr_oper_list,&corr_copy_list,
			     &p_momentum , &q_momentum , 
			     &nt, &no_p_mom , &no_q_mom , &no_oper,
			     &no_spectator , &no_sequential, &no_zonked, 
			     &hl_flag , input_3pt_file, &nocopies );



  /*** write out some titles *****/
  titles_hl(nt, no_p_mom, no_q_mom, no_oper ,p_momentum, q_momentum , hl_flag) ;

  /*** write the three point functions to disk ****/

  p_mom_select = 0;
  write_3pt_onemom(corr, 
		   output1_file, p_mom_select, q_momentum, p_momentum, hl_flag,
		   no_zonked , no_spectator,  no_sequential,
		   no_p_mom, no_q_mom, no_oper, nt, nocopies)  ;

  if(no_p_mom < 2)
    {
      printf("ERROR: Can't write 2nd file. Only one p momentum available\n");
      return 1;
    }

  p_mom_select = 1;
  write_3pt_onemom(corr, 
		   output2_file, p_mom_select, q_momentum, p_momentum, hl_flag,
		   no_zonked , no_spectator,  no_sequential,
		   no_p_mom, no_q_mom, no_oper, nt, nocopies)  ;

  if(no_output_files == 2)return 0;

  if(no_p_mom < 3)
    {
      printf("ERROR: Can't write 3rd file. Only two p momenta available\n");
      return 1;
    }
    
  p_mom_select = 2;
  write_3pt_onemom(corr, 
		   output3_file, p_mom_select, q_momentum, p_momentum, hl_flag,
		   no_zonked , no_spectator,  no_sequential,
		   no_p_mom, no_q_mom, no_oper, nt, nocopies)  ;

  return 0 ;
}

