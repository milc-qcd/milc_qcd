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
  int *corr_oper_list;
  int *corr_copy_list;
  int nt ;
  int hl_flag ;
  int no_q_mom ;
  int no_oper ;

  int no_zonked = 1 ;
  int no_spectator = 1 ;
  int nocopies = 1; 

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


  /*** read and reserve memory for the three point function ******/
  read_prop_twopt(&corr, &corr_oper_list, &corr_copy_list,
		  hqet_file ,   &q_momentum , 
		  &hl_flag , &no_spectator , &no_zonked,  
		  &no_oper, &nt, &no_q_mom , &nocopies );



  /*** write out some titles *****/
  titles_2pt(nt, no_q_mom, no_oper,q_momentum, hl_flag, nocopies) ;

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
