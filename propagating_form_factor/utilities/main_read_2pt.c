/*
   Read the two point function from a file and write
   it to the screen


*/

#include "read_hl_form.h"
#include "prop_form_utilities.h"



/**** required functions prototypes ******/


void read_prop_twopt(complex **corr, int **corr_oper_list,
		     int **corr_copy_list,
		     char filename[80], int32type **q_momstore,
int *hl_flag, int *no_k_one , int *no_k_two, int *numer_of_operators, 
int *nt , int *no_q_values, int *nocopies) ;


void dump_ll_twopt(complex *corr,  int corr_oper_list[],
  int corr_copy_list[],
  int no_zonked_light, int no_spectator, 
  int no_q_values, int no_oper, int nt) ;


void dump_hl_twopt(complex *corr,  int corr_oper_list[],
  int corr_copy_list[],
  int no_zonked_heavy, int no_spectator, 
  int no_q_values, int no_oper, int nt, int nocopies) ;


void titles(int nt, int no_q_mom, int no_oper,int *q_momentum, int hl_flag) ;

/***** end of function prototypes *****/


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
  int nocopies ; 

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
  titles(nt, no_q_mom, no_oper,q_momentum, hl_flag) ;

  /*** write the two point function to the terminal ****/

  if( hl_flag ==  HL_2PT_BAG  )
  {
    printf("---------------------------------------\n"); 
    printf("Heavy-light SHELL smearing function\n"); 
    printf("---------------------------------------\n"); 
    dump_hl_twopt(corr, corr_oper_list, corr_copy_list,
		  no_zonked , no_spectator, no_q_mom, no_oper, nt, nocopies)  ;
  }
  else if( hl_flag == HL_REL_2PT )
  {
    printf("---------------------------------------\n"); 
    printf("Heavy-light RELATIVE smearing function\n"); 
    printf("---------------------------------------\n"); 

    dump_hl_twopt(corr, corr_oper_list, corr_copy_list,
		  no_zonked , no_spectator, no_q_mom, no_oper, nt, nocopies)  ;
  }
  else if( hl_flag == HL_LOCAL_SINK )
  {
    printf("-------------------------------------------------------\n"); 
    printf("Heavy-light LOCAL-SOURCE SMEARED-SINK smearing function\n"); 
    printf("-------------------------------------------------------\n"); 

    dump_hl_twopt(corr, corr_oper_list, corr_copy_list,
		  no_zonked , no_spectator, no_q_mom, no_oper, nt, nocopies)  ;
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

    dump_ll_twopt(corr, corr_oper_list, corr_copy_list,
		  no_zonked , no_spectator, no_q_mom,  no_oper, nt)  ;
  }
  else
  {
    printf("ERROR hl_flag =%d is out of range\n",hl_flag); 
    exit(1) ; 
  }



  return 0 ;
}
