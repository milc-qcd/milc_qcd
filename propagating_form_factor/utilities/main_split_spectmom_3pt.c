/*
   Read the three point function from a file with two or three source momenta
   and split it into separate files according to source momentum

   2/14/99 This version requires less memory.  It reads only one
   momentum set at a time from the input file. C.D.

*/

#include "read_hl_form.h"
#include "prop_form_utilities.h"
#define MAX 16


int main(int argc, char *argv[])
{
  char *input_3pt_file, *output_file[MAX];

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
  int spect_select;
  int k_output,j;

  /***--------------------------------------------------*****/

  printf("=============================================================\n");
  printf(" Splitting the three point function (propagating form factor code)  \n");
  printf("=============================================================\n");

  if( argc < 4 )
  {
    printf("usage::  %s  <input 3pt file> <output p 0 sp 0> <output p 0 sp 1> < output p 0 sp 2> <output p 1 sp 0> etc\n",argv[0]);
    return 1;
  }

  no_output_files = argc - 2;
  if(no_output_files >= MAX){
    fprintf(stderr,"Number of output files exceeds %d\nFix the program\n",
	    MAX);
    return 1;
  }
  input_3pt_file = argv[1];
  for(k_output = 0,j = 2; k_output < no_output_files; k_output++,j++)
    output_file[k_output] = argv[j];

  p_mom_select = 0;
  spect_select = 0;
  k_output = 0;
  no_p_mom = no_spectator = 1;

  while(p_mom_select < no_p_mom){
    
    /*** read the three point function ******/
    read_3pt_onemom_onespect(&corr,&corr_oper_list,&corr_copy_list,
			     &p_momentum , &q_momentum , 
			     &nt, &no_p_mom , &no_q_mom , &no_oper,
			     &no_spectator , &no_sequential, &no_zonked, 
			     &hl_flag , input_3pt_file, &nocopies,
			     p_mom_select, spect_select);

    if(no_output_files != no_p_mom*no_spectator){
      fprintf(stderr,"File %s has %d p momenta and %d spectators but I have only %d output files\n",input_3pt_file,no_p_mom,no_spectator,no_output_files);
      return 1;
    }
    
    /*** write out some titles *****/
    titles_hl(nt, 1, no_q_mom, no_oper ,p_momentum, q_momentum , hl_flag) ;
    
    /*** write the three point functions to disk ****/
    
    write_3pt_onemom(corr, 
		     output_file[k_output], 0, q_momentum, p_momentum, hl_flag,
		     no_zonked , 1,  no_sequential,
		     1, no_q_mom, no_oper, nt, nocopies)  ;
    
    
    free(corr); free(corr_oper_list); free(corr_copy_list);
    free(q_momentum); free(p_momentum);
    
    k_output++;
    spect_select++;
    if(spect_select >= no_spectator){
      spect_select = 0;
      p_mom_select++;
    }
  }
    
    return 0 ;
}

