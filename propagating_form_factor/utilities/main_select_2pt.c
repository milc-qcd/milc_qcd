/*
   Read the two point function from a file and write
   it to the screen


*/

#include "read_hl_form.h"
#include "prop_form_utilities.h"


int main(int argc, char *argv[])
{
  char hqet_file[80] ;
  int32type  *q_momentum  ;
  complex  *corr ;
  int nt ;
  int hl_flag ;
  int no_q_mom = 1 ;
  int no_oper = 1 ;

  int no_zonked = 1 ;
  int no_spectator = 1 ;
  int nocopies = 1 ; 

  int corr_oper_list[1];
  int corr_copy_list[1];

  int nselect;
  twopt_select twoselect[MAX_SELECT_PER_FILE];

  /***--------------------------------------------------*****/

  printf("==================================================================\n");
  printf(" Raw output of the two point function (propagating form factor code)  \n");
  printf("==================================================================\n");

  if( argc != 2 )
  {
    printf("usage::  %s  [propagating form factor 2pt file] \n",argv[0]);
    exit(1);
  }
  strcpy(hqet_file,argv[1]);

  read_select_twopt_param(&nselect,twoselect);

  /*** read and reserve memory for the two point function ******/
  read_select_twopt(&corr, hqet_file , nselect, twoselect, &q_momentum,
		    &hl_flag , &nt );
  /* For the dump we report only the first selection in the list */
  corr_oper_list[0] = twoselect[0].oper;
  corr_copy_list[0] = twoselect[0].copy;

  /*** write out some titles *****/
  titles_2pt(nt, no_q_mom, no_oper, q_momentum, hl_flag, nocopies) ;

  /*** write the two point function to the terminal ****/

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
  
  dump_twopt(corr, corr_oper_list, corr_copy_list,
	     no_zonked , no_spectator, no_q_mom, no_oper, nt, nocopies)  ;
  
  return 0 ;
}
