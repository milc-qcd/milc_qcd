/*
 *   Write out the titles for the code that reads 
 *   the heavy-light hqet form factors
 * 
 *
 */


#include<stdio.h>
#include<math.h>

#include "read_hl_form.h"
#include "prop_form_utilities.h"

void titles_2pt(int nt, int no_q_mom, int no_oper,int32type *q_momentum, int hl_flag, int no_copies)
{
  int imom ;
  int i ;


  printf("---------------------------------------\n"); 
  if( hl_flag ==  HL_2PT_BAG  )
  {
    printf("Heavy-light SHELL smearing function\n"); 
  }
  else if( hl_flag == HL_REL_2PT )
  {
    printf("Heavy-light RELATIVE smearing function\n"); 
  }
  else if( hl_flag == LL_2PT )
  {
    printf("Light-light  smearing function\n"); 
  }
  else if( hl_flag == HL_LOCAL_SINK )
  {
    printf("Heavy-light local-source smeared-sink smearing function\n"); 
  }
  else
  {
    printf("ERROR hl_flag =%d is out of range\n",hl_flag); 
    exit(1) ; 
  }
  printf("---------------------------------------\n"); 





  printf("\nNumber of timeslices = %d\n",nt);
  printf("The number of operators = %d\n",no_oper);
  printf("Number of q-momenta     = %d\n",no_q_mom   );
  printf("Number of copies        = %d\n",no_copies   );




  printf("\nHere is a list of the q-momentum  \n");
  for(imom=0 ; imom < no_q_mom ; ++imom)
  {
    printf("%d   %d   %d\n",*(q_momentum + 3*imom),
	   *(q_momentum + 1 + 3*imom),
	   *(q_momentum + 2 + 3*imom)  ) ;
  }





}
