/*
   Read the three point function from a file and write
   it to the screen


*/

#include "read_hl_form.h"
#include "prop_form_utilities.h"


int main(int argc, char *argv[])
{
  char hqet_file[80] ;
  int32type  *p_momentum  ;
  int32type  *q_momentum  ;
  complex  *corr ;
  int corr_oper_list[1];
  int corr_copy_list[1];
  int nt ;
  int hl_flag ;
  int no_q_mom = 1 ;
  int no_p_mom = 1 ;
  int no_oper = 1 ;

  int no_zonked = 1 ;
  int no_spectator = 1 ;

  int no_sequential = 1 ; 
  int nocopies = 1 ;

  int nselect;
  threept_select threeselect[MAX_SELECT_PER_FILE];

  /***--------------------------------------------------*****/

  printf("=============================================================\n");
  printf(" Raw output of the propagating heavy --> light form factors  \n");
  printf("=============================================================\n");

  if( argc != 2 )
  {
    printf("usage::  %s  [propagating form factor data file] \n",argv[0]);
    exit(1);
  }
  strcpy(hqet_file,argv[1]);

  read_select_form_param(&nselect,threeselect);

  /*** read and reserve memory for the three point function ******/
  read_select_form_corr(&corr, hqet_file, nselect, threeselect,
			&p_momentum, &q_momentum,
			&hl_flag , &nt );

  /* For the dump we report only the first selection in the list */
  corr_oper_list[0] = threeselect[0].oper;
  corr_copy_list[0] = threeselect[0].copy;


  /*** write out some titles *****/
  titles_hl(nt, no_p_mom, no_q_mom, no_oper ,p_momentum, q_momentum , hl_flag) ;



  /*** write the two point function to the terminal ****/

  dump_prop_form(corr, corr_oper_list, corr_copy_list,
		 no_zonked , no_spectator, no_sequential, 
		 no_q_mom, no_p_mom, no_oper, nt, nocopies)  ;

  return 0 ;
}
