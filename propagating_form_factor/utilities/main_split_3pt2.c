/*
   Read the three point function from a file with two or three source momenta
   and split it into separate files according to source momentum

   2/14/99 This version requires less memory.  It reads only one
   momentum set at a time from the input file. C.D.

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
  p_mom_select = 0;
  read_3pt_onemom(&corr,&corr_oper_list,&corr_copy_list,
		  &p_momentum , &q_momentum , 
		  &nt, &no_p_mom , &no_q_mom , &no_oper,
		  &no_spectator , &no_sequential, &no_zonked, 
		  &hl_flag , input_3pt_file, &nocopies,p_mom_select );


  /*** write out some titles *****/
  titles_hl(nt, no_p_mom, no_q_mom, no_oper ,p_momentum, q_momentum , hl_flag) ;

  /*** write the three point functions to disk ****/

  p_mom_select = 0;
  write_3pt_onemom(corr, 
		   output1_file, 0, q_momentum, p_momentum, hl_flag,
		   no_zonked , no_spectator,  no_sequential,
		   no_p_mom, no_q_mom, no_oper, nt, nocopies)  ;


  free(corr); free(corr_oper_list); free(corr_copy_list);
  free(q_momentum); free(p_momentum);

  p_mom_select = 1;
  read_3pt_onemom(&corr,&corr_oper_list,&corr_copy_list,
		  &p_momentum , &q_momentum , 
		  &nt, &no_p_mom , &no_q_mom , &no_oper,
		  &no_spectator , &no_sequential, &no_zonked, 
		  &hl_flag , input_3pt_file, &nocopies,p_mom_select );


  /*** write out some titles *****/
  titles_hl(nt, no_p_mom, no_q_mom, no_oper ,p_momentum, q_momentum , hl_flag) ;

  write_3pt_onemom(corr, 
		   output2_file, 0, q_momentum, p_momentum, hl_flag,
		   no_zonked , no_spectator,  no_sequential,
		   no_p_mom, no_q_mom, no_oper, nt, nocopies)  ;

  free(corr); free(corr_oper_list); free(corr_copy_list);
  free(q_momentum); free(p_momentum);

  if(no_output_files == 2)return 0;

  p_mom_select = 2;
  read_3pt_onemom(&corr,&corr_oper_list,&corr_copy_list,
		  &p_momentum , &q_momentum , 
		  &nt, &no_p_mom , &no_q_mom , &no_oper,
		  &no_spectator , &no_sequential, &no_zonked, 
		  &hl_flag , input_3pt_file, &nocopies,p_mom_select );


  /*** write out some titles *****/
  titles_hl(nt, no_p_mom, no_q_mom, no_oper ,p_momentum, q_momentum , hl_flag) ;

  write_3pt_onemom(corr, 
		   output3_file, 0, q_momentum, p_momentum, hl_flag,
		   no_zonked , no_spectator,  no_sequential,
		   no_p_mom, no_q_mom, no_oper, nt, nocopies)  ;

  free(corr); free(corr_oper_list); free(corr_copy_list);
  free(q_momentum); free(p_momentum);

  return 0 ;
}

