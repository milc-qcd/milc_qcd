/*
   Read a two point function from a file with two or three spectator
   kappas and split it into separate files according to spectator kappa

*/

#include "read_hl_form.h"
#include "prop_form_utilities.h"

int main(int argc, char *argv[])
{
  char *input_2pt_file, *output1_file, *output2_file, *output3_file ;

  int32type  *p_momentum  ;
  complex  *corr ;
  int *corr_oper_list;
  int *corr_copy_list;
  int nt ;
  int hl_flag ;
  int no_p_mom ;
  int no_oper ;

  int no_output_files ;

  int no_zonked = 1 ;
  int no_spectator = 1 ;
  int nocopies ; 

  int spect_select;

  /***--------------------------------------------------*****/

  printf("==================================================================\n");
  printf(" Splitting the two point function (propagating form factor code)  \n");
  printf("==================================================================\n");

  if( argc < 4 )
  {
    printf("usage::  %s  <input 2pt file> <output1> <output2>  [<output3>] \n",argv[0]);
    exit(1);
  }

  no_output_files = argc - 2;
  input_2pt_file = argv[1];
  output1_file = argv[2];
  output2_file = argv[3];
  output3_file = argv[4];

  /*** read and reserve memory for the three point function ******/
  read_prop_twopt(&corr, &corr_oper_list, &corr_copy_list,
		  input_2pt_file ,   &p_momentum , 
		  &hl_flag , &no_spectator , &no_zonked,  
		  &no_oper, &nt, &no_p_mom , &nocopies );


  /*** write out some titles *****/
  titles_2pt(nt, no_p_mom, no_oper,p_momentum, hl_flag, nocopies) ;

  /*** write the two point functions to disk ****/

  if( hl_flag ==  HL_2PT_BAG  )
  {
    printf("---------------------------------------\n"); 
    printf("Heavy-light SHELL smearing function\n"); 
    printf("---------------------------------------\n"); 
  }
  else if( hl_flag == HL_REL_2PT )
  {
    printf("---------------------------------------\n"); 
    printf("Heavy-light RELATIVE smearing function\n"); 
    printf("---------------------------------------\n"); 

  }
  else if( hl_flag == HL_LOCAL_SINK )
  {
    printf("-------------------------------------------------------\n"); 
    printf("Heavy-light LOCAL-SOURCE SMEARED-SINK smearing function\n"); 
    printf("-------------------------------------------------------\n"); 
  }
  else if( hl_flag == LL_2PT )
  {
    printf("---------------------------------------\n"); 
    printf("Light-light RELATIVE smearing function\n"); 
    printf("---------------------------------------\n"); 

    if( nocopies != 1 )
    {
      printf("ERROR:main: nocopies = %d <> 1 \n",nocopies); 
      exit(2); 
    }
  }
  else
    {
      printf("ERROR hl_flag =%d is out of range\n",hl_flag); 
      exit(1) ; 
    }
  
  spect_select = 0;
  write_prop_twopt_onespect(corr, 
	      output1_file, spect_select, p_momentum, hl_flag,
              no_zonked , no_spectator,  no_p_mom, no_oper, nt, nocopies)  ;

  spect_select = 1;
  write_prop_twopt_onespect(corr, 
	      output2_file, spect_select, p_momentum, hl_flag,
	      no_zonked , no_spectator,  no_p_mom, no_oper, nt, nocopies)  ;

  if(no_output_files == 2)return 0;
  
  spect_select = 2;
  write_prop_twopt_onespect(corr, 
	      output3_file, spect_select, p_momentum, hl_flag,
	      no_zonked , no_spectator,  no_p_mom, no_oper, nt, nocopies)  ;
  
  return 0 ;
}
