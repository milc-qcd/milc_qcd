/*
   Read the three point function from a file with two or three
   spectator kappas and split it into separate files according to
   spectator quark kappa values

   4/23/99  C.D.

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
  int spect_select;

  /***--------------------------------------------------*****/

  printf("=============================================================\n");
  printf(" Splitting the three point function (propagating form factor code)  \n");
  printf("=============================================================\n");

  if( argc < 4 )
  {
    printf("usage::  %s  <input 3pt file> <output file for 1st sp> <output file for 2nd sp> [< output file for 3rd sp>]\n",argv[0]);
    exit(1);
  }

  no_output_files = argc - 2;
  input_3pt_file = argv[1];
  output1_file = argv[2];
  output2_file = argv[3];
  output3_file = argv[4];

  /*** read and reserve memory for the three point function ******/
  spect_select = 0;
  read_3pt_onespect(&corr,&corr_oper_list,&corr_copy_list,
		  &p_momentum , &q_momentum , 
		  &nt, &no_p_mom , &no_q_mom , &no_oper,
		  &no_spectator , &no_sequential, &no_zonked, 
		  &hl_flag , input_3pt_file, &nocopies,spect_select );


  /*** write out some titles *****/
  titles_hl(nt, no_p_mom, no_q_mom, no_oper ,p_momentum, q_momentum , hl_flag) ;

  /*** write the three point functions to disk ****/

  spect_select = 0;
  write_3pt_onespect(corr, 
		   output1_file, 0, q_momentum, p_momentum, hl_flag,
		   no_zonked , no_spectator,  no_sequential,
		   no_p_mom, no_q_mom, no_oper, nt, nocopies)  ;


  free(corr); free(corr_oper_list); free(corr_copy_list);
  free(q_momentum); free(p_momentum);

  spect_select = 1;
  read_3pt_onespect(&corr,&corr_oper_list,&corr_copy_list,
		  &p_momentum , &q_momentum , 
		  &nt, &no_p_mom , &no_q_mom , &no_oper,
		  &no_spectator , &no_sequential, &no_zonked, 
		  &hl_flag , input_3pt_file, &nocopies,spect_select );


  /*** write out some titles *****/
  titles_hl(nt, no_p_mom, no_q_mom, no_oper ,p_momentum, q_momentum , hl_flag) ;

  write_3pt_onespect(corr, 
		   output2_file, 0, q_momentum, p_momentum, hl_flag,
		   no_zonked , no_spectator,  no_sequential,
		   no_p_mom, no_q_mom, no_oper, nt, nocopies)  ;

  free(corr); free(corr_oper_list); free(corr_copy_list);
  free(q_momentum); free(p_momentum);

  if(no_output_files == 2)return 0;

  spect_select = 2;
  read_3pt_onespect(&corr,&corr_oper_list,&corr_copy_list,
		  &p_momentum , &q_momentum , 
		  &nt, &no_p_mom , &no_q_mom , &no_oper,
		  &no_spectator , &no_sequential, &no_zonked, 
		  &hl_flag , input_3pt_file, &nocopies,spect_select );


  /*** write out some titles *****/
  titles_hl(nt, no_p_mom, no_q_mom, no_oper ,p_momentum, q_momentum , hl_flag) ;

  write_3pt_onespect(corr, 
		   output3_file, 0, q_momentum, p_momentum, hl_flag,
		   no_zonked , no_spectator,  no_sequential,
		   no_p_mom, no_q_mom, no_oper, nt, nocopies)  ;

  free(corr); free(corr_oper_list); free(corr_copy_list);
  free(q_momentum); free(p_momentum);

  return 0 ;
}

