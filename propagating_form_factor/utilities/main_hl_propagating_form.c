/*
   Read the three point function from a file and write
   it to the screen


*/

#include"read_hl_form.h"



/**** required functions prototypes ******/


void read_propagating_form_corr(complex **corr, int **corr_oper_list,
  int **corr_copy_list,
  int32type **p_momstore, int32type **q_momstore,
  int *nt, int *no_p_values, int *no_q_values, int *no_oper,
  int *no_spectator,   int *no_sequential, int *no_zonked , int *hl_flag,
  char filename[80], int *nocopies ) ;

void dump_hl_prop_form(complex *corr,  int corr_oper_list[],
      int corr_copy_list[],
      int no_zonked_light, int no_spectator, int no_sequential,
      int no_q_values, int no_p_values, int no_oper, int nt, int nocopies);


/***** end of function prototypes *****/


int main(int argc, char *argv[])
{
  char hqet_file[80] ;
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

  int no_zonked = 1 ;
  int no_spectator = 1 ;

  int no_sequential = 1 ; 
  int nocopies ;

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


  /*** read and reserve memory for the three point function ******/
  read_propagating_form_corr(&corr,&corr_oper_list,&corr_copy_list,
			     &p_momentum , &q_momentum , 
			     &nt, &no_p_mom , &no_q_mom , &no_oper,
			     &no_spectator , &no_sequential, &no_zonked, 
			     &hl_flag , hqet_file, &nocopies );



  /*** write out some titles *****/
  titles_hl(nt, no_p_mom, no_q_mom, no_oper ,p_momentum, q_momentum , hl_flag) ;



  /*** write the two point function to the terminal ****/

  if( hl_flag == HEAVY_TO_LIGHT  )
  {
    dump_hl_prop_form(corr, corr_oper_list, corr_copy_list,
		      no_zonked , no_spectator, no_sequential, 
		      no_q_mom, no_p_mom, no_oper, nt,  nocopies)  ;
  }
  else if( hl_flag == HEAVY_TO_HEAVY )
  {
    dump_hh_prop_form(corr, corr_oper_list, corr_copy_list,
		no_zonked , no_spectator, no_sequential, no_q_mom, no_p_mom, 
		      no_oper, nt, nocopies)  ;
  }


  return 0 ;
}
