/*
   Read the three point functions from two file and write
   them and their difference to the screen.


*/

#include"read_hl_form.h"



/**** required functions prototypes ******/

void titles(int nt, int novel, int no_mom, int no_oper,
Real  *velocity,int  *momentum)  ;

void read_hqet_form_corr(complex **corr, int **corr_oper_list, 
  int **corr_copy_list,
  int **q_momstore, Real **velocity,
  int *nt, int *novel, int *no_q_values, int *no_oper,
  int *no_spectator , int *no_zonked , 
  char filename[80] )  ;


void diff_hl_form(complex *corr_one,  complex *corr_two, 
  int corr_oper_list[], int corr_copy_list[],
  int no_zonked, int no_spectator, int no_q_values, int novel,
  int no_oper, int nt)  ;


/***** end of function prototypes *****/


int main(int argc, char *argv[])
{
  char hqet_file[80] ;
  char hqet_file_one[80] ;

  Real  *velocity ,  *velocity_one  ;
  int  *momentum ,  *momentum_one ;
  complex  *corr ,   *corr_one  ;
  int *corr_oper_list;
  int *corr_copy_list;
  int nt , nt_one ;
  int novel , novel_one  ;
  int no_mom , no_mom_one ;
  int no_oper , no_oper_one;

  int no_zonked = 1 ;
  int no_spectator = 1 ;


  int no_zonked_one = 1 ;
  int no_spectator_one = 1 ;



  /***--------------------------------------------------*****/

  printf("==================================================\n");
  printf(" Comparison  of the hqet heavy --> light three pt  \n");
  printf("==================================================\n");

  if( argc != 3 )
  {
    printf("usage::  %s  [hqet 3pt file ] [hqet 3pt file ] \n",argv[0]);
    exit(1);
  }
  strcpy(hqet_file,argv[1]);
  strcpy(hqet_file_one,argv[2]);


  /*** read and reserve memory for the three point function ******/
  read_hqet_form_corr(&corr,&corr_oper_list,&corr_copy_list,
		      &momentum , &velocity,
		      &nt, &novel, &no_mom , &no_oper,
		      &no_spectator , &no_zonked , hqet_file );


  read_hqet_form_corr(&corr_one,&corr_oper_list, &corr_copy_list,
		      &momentum_one , &velocity_one,
		      &nt_one, &novel_one, &no_mom_one , &no_oper_one,
		      &no_spectator_one , &no_zonked_one , hqet_file_one );


  /*** check the parameters of the two files ****/
  if( nt != nt_one  )
  {
    printf("ERROR: mismatch between the timeslices of the two files %d != %d \n",nt,nt_one);
    exit(1); 
  }

  if( novel != novel_one  )
  {
    printf("ERROR: mismatch between the velocites of the two files %d != %d \n",novel,novel_one);
    exit(1); 
  }

  if( no_oper != no_oper_one  )
  {
    printf("ERROR: mismatch between the number of operators of the two files %d != %d \n",no_oper,no_oper_one);
    exit(1); 
  }

  if( no_mom != no_mom_one  )
  {
    printf("ERROR: mismatch between the momentum of the two files %d != %d \n",no_mom,no_mom_one);
    exit(1); 
  }



  /*** write out some titles *****/
  titles(nt, novel,no_mom, no_oper,velocity,momentum) ;

  /*** write the two point function to the terminal ****/
  diff_hl_form(corr,corr_one,corr_oper_list,corr_copy_list,
	       no_zonked,no_spectator,no_mom,novel,no_oper,nt) ;



  /*** free up the reserved memory  **/
  free(corr);
  free(corr_one); 
  free(corr_oper_list);
  free(corr_copy_list);
  free(velocity) ;
  free(velocity_one) ;
  free(momentum) ;
  free(momentum_one) ;

  return 0 ;
}
