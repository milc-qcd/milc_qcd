/*
 *
 *  Dump the heavy --> light two  point functions to 
 *  the screen.
 *
 *
 *
 */


#include "read_hl_form.h"

char *two_oper_name(int n ) ;
char *oper_name_and_corrections(int n , int copy_pt ) ;


void dump_twopt(complex *corr, int corr_oper_list[], int corr_copy_list[],
		int no_zonked, int no_spectator, 
		int no_q_values, int no_oper, int nt, int nocopies)
{
  int where ;
  int t ;
  int zonk_pt = no_zonked -1 ;
  int spect_pt = no_spectator -1 ; 
  int q_pt = no_q_values -1 ; 
  int oper_pt = no_oper -1 ;

  int copy_pt ; 
  int corr_stride ;
  /**************************************************/

  corr_stride = TWOPT_FORM_WHERE(nt,zonk_pt ,spect_pt,q_pt, oper_pt) ; ; 

  printf("Here is the heavy-light two point function\n");


  for(zonk_pt = 0 ; zonk_pt < no_zonked ; ++zonk_pt )
  {
    printf("zonked quark = %d\n",zonk_pt);
    for(spect_pt = 0 ; spect_pt < no_spectator ; ++spect_pt )
    {
      printf("  spectator quark = %d\n",spect_pt);
      for(q_pt = 0 ; q_pt < no_q_values ; ++q_pt ) 
      {
	printf("   Q--Momentum pointer = %d \n",q_pt  ); 
	for(copy_pt = 0 ; copy_pt < nocopies ; ++copy_pt)
	  for(oper_pt = 0 ; oper_pt < no_oper ; ++oper_pt)
	  {
	    printf("       Operator pointer = %d  name = %s\n",
		   corr_oper_list[oper_pt], 
		   oper_name_and_corrections(corr_oper_list[oper_pt],
					     corr_copy_list[copy_pt]) );
	    for(t = 0 ; t < nt ; ++t)
	    {
	      where = corr_stride*copy_pt + TWOPT_FORM_WHERE(t,zonk_pt ,spect_pt,q_pt, oper_pt) ;
	      printf("         corr[ %d ] = ( %g , %g )\n",t , (corr+where)->real, (corr+where)->imag);
	    }

		

	  }
      }
    }
  }
  


}  /*** end of the dump of the heavy-light 2pt functions ****/


