/* 
 *  Return the name of the three point function
 *
 *  See K&R section 5.8, 113
 *
 */

#include "read_hl_form.h"
#include <assert.h>
#include "../opertypes.h"
#include "../../generic_form/gammatypes.h"
#include "../corrlist.h"

char *two_oper_name(int n )
{

  assert( n >= 0 && n <  MAX_TWOPT ); 

  return twopt_name[two_pt[n].oper] ; 
}


/* 
 *  Return the name of the two point function
 *  Subroutine arguments:
 *    On input
 *      n       :: number of dirac matrix in the operator
 *      copy_pt :: the type of operator insertion
 */

char *oper_name_and_corrections(int n , int copy_pt )
{
  static char name[80] ; 

  assert( n >= 0 && n < NO_TWOPT_OPERS ); 


  switch( copy_pt ) 
  {
   case 0 :  
     sprintf(name,"%s",two_oper_name(n))  ;
     break ;
   case 1 :  
     sprintf(name,"%s_3D_on_spect", two_oper_name(n) )  ;
     break ;
   case 2 :  
     sprintf(name,"%s_3D_on_zonked",two_oper_name(n)) ; 
     break ;
   case 3 :  
     sprintf(name,"%s_4D_on_spect",two_oper_name(n)) ;
     break ;
   case 4 :  
     sprintf(name,"%s_4D_on_zonked",two_oper_name(n) )  ;
     break ;
   default :
     printf("ERROR:oper_name_and_corrections: copy_pt = %d is out of range\n",copy_pt); 
     exit(1) ;
     break ;
  }



  return &name[0] ;


}

