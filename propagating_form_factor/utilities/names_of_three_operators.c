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

char *three_oper_name(int n , int copy_pt)
{
  static char name[80] ; 

  assert( n >= 0 && n <   MAX_THREEPT ); 



  switch( copy_pt ) 
  {
   case 0 :  
     sprintf(name,"%s",name_3pt[three_pt[n].oper])  ;
     break ;
   case 1 :  
     sprintf(name,"%s_3D_on_seq",name_3pt[three_pt[n].oper])  ;
     break ;
   case 2 :  
     sprintf(name,"%s_3D_on_zonked",name_3pt[three_pt[n].oper]) ; 
     break ;
   case 3 :  
     sprintf(name,"%s_4D_on_seq",name_3pt[three_pt[n].oper] ) ;
     break ;
   case 4 :  
     sprintf(name,"%s_4D_on_zonked",name_3pt[three_pt[n].oper] )  ;
     break ;
   default :
     printf("ERROR:three_oper_name: copy_pt = %d is out of range\n",copy_pt); 
     exit(1) ;
     break ;
  }




  return &name[0] ;
}

int oper_chirality(int oper)
{
  int gamma_chirality[MAXGAMMA];    /* Sign change under gamma5 conjugation */

  gamma_chirality[GX ] = -1;
  gamma_chirality[GY ] = -1;
  gamma_chirality[GZ ] = -1;
  gamma_chirality[GT ] = -1;
  gamma_chirality[G5 ] =  1;
  gamma_chirality[GYZ] =  1;
  gamma_chirality[GZX] =  1;
  gamma_chirality[GXY] =  1;
  gamma_chirality[GXT] =  1;
  gamma_chirality[GYT] =  1;
  gamma_chirality[GZT] =  1;
  gamma_chirality[G5X] = -1;
  gamma_chirality[G5Y] = -1;
  gamma_chirality[G5T] = -1;
  gamma_chirality[G5Z] = -1;
  gamma_chirality[G1 ] =  1;

  return gamma_chirality[three_pt[oper].gin];
}

int oper_phase(int oper)
{
  int gamma_phase[MAXGAMMA];    /* Phase of operator */

  gamma_phase[GX ] =  0;
  gamma_phase[GY ] =  0;
  gamma_phase[GZ ] =  0;
  gamma_phase[GT ] =  0;
  gamma_phase[G5 ] =  0;
  gamma_phase[GYZ] =  0;
  gamma_phase[GZX] =  0;
  gamma_phase[GXY] =  0;
  gamma_phase[GXT] =  0;
  gamma_phase[GYT] =  0;
  gamma_phase[GZT] =  0;
  gamma_phase[G5X] =  0;
  gamma_phase[G5Y] =  0;
  gamma_phase[G5T] =  0;
  gamma_phase[G5Z] =  0;
  gamma_phase[G1 ] =  0;

  return gamma_phase[three_pt[oper].gin];
}



