/*
 *   Write out the titles for the code that reads 
 *   the heavy-light hqet form factors
 * 
 *
 */


#include<stdio.h>
#include<math.h>

#include"read_hl_form.h"
#include "prop_form_utilities.h"

void titles_hl(int nt, int no_p_mom, int no_q_mom, int no_oper,
int32type  *p_momentum, int32type  *q_momentum, int hl_flag)
{
  int imom ;
  int i ;

  if( hl_flag == HEAVY_TO_HEAVY ) 
  {
    printf("HEAVY_TO_HEAVY  form factor\n"); 
  }
  else if (hl_flag == HEAVY_TO_LIGHT )
  {
    printf("HEAVY_TO_LIGHT  form factor\n"); 
  }
  else
  {
    printf("ERROR:: hl_flag = %d is out of range\n",hl_flag) ; 
    exit(1) ; 
  }


  printf("\nNumber of timeslices = %d\n",nt);
  printf("The number of operators = %d\n",no_oper);
  printf("Number of p-momentum    = %d\n",no_p_mom   );
  printf("Number of q-momentum    = %d\n",no_q_mom   );


  printf("\nHere is a list of the p-momentum  \n");
  for(imom=0 ; imom < no_p_mom ; ++imom)
  {
    printf("%d   %d   %d\n",*(p_momentum + 3*imom),
	   *(p_momentum + 1 + 3*imom),
	   *(p_momentum + 2 + 3*imom)  ) ;
  }




  printf("\nHere is a list of the q-momentum  \n");
  for(imom=0 ; imom < no_q_mom ; ++imom)
  {
    printf("%d   %d   %d\n",*(q_momentum + 3*imom),
	   *(q_momentum + 1 + 3*imom),
	   *(q_momentum + 2 + 3*imom)  ) ;
  }





}
