/*
   Read the two point function from a file, fix its normalization
   and write it to a new file.

*/

#include "read_hl_form.h"
#include "prop_form_utilities.h"


int main(int argc, char *argv[])
{
  char *input_2pt_file, *output_2pt_file;
  int32type  *k_momentum  ;
  complex  *corr ;
  int nt ;
  int hl_flag ;
  int no_k_mom = 1 ;
  int no_q_values = 1 ;
  int number_of_operators = 1 ;

  int no_zonked = 1 ;
  int no_spectator = 1 ;
  int no_copies = 1 ; 

  int zonked_pt, spect_pt, q_pt, oper_pt, copy_pt;

  int *corr_oper_list;
  int *corr_copy_list;

  int corr_stride_output;
  int i,dim;

  /***--------------------------------------------------*****/

  printf("==================================================================\n");
  printf(" Fixing the normalization of the two point function (propagating form factor code)  \n");
  printf("==================================================================\n");

  if( argc != 3 )
  {
    printf("usage::  %s  [input 2pt file] [output 2pt file]\n",argv[0]);
    exit(1);
  }

  input_2pt_file = argv[1];
  output_2pt_file = argv[2];

  /*** read and reserve memory for the three point function ******/
  read_prop_twopt(&corr, &corr_oper_list, &corr_copy_list,
		  input_2pt_file ,   &k_momentum , 
		  &hl_flag , &no_spectator , &no_zonked,  
		  &number_of_operators, &nt, &no_k_mom , &no_copies );

  no_q_values = no_k_mom;

  /*** write out some titles *****/
  titles_2pt(nt, no_k_mom, number_of_operators, k_momentum, hl_flag, 
	     no_copies) ;

  /* Renormalize the 2 pt correlators */
  
  /* compute dim = number of complex values in corr */
  corr_stride_output = 
    TWOPT_FORM_WHERE(nt,no_zonked-1,no_spectator-1,no_q_values-1,
			number_of_operators-1);

  dim = corr_stride_output * no_copies;

  /* We divide every number by 2 here */
  for(i = 0; i < dim; i++)
    {
      corr[i].real /= 2;
      corr[i].imag /= 2;
    }

  /*** write the two point functions to disk ****/

  write_prop_twopt(corr, 
		   output_2pt_file, k_momentum, 
		   hl_flag, no_spectator,  no_zonked , 
		   number_of_operators, nt, no_k_mom, no_copies)  ;
  
  return 0 ;
}
