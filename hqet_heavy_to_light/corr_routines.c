/************* corr_routines.c **************************/
/* MIMD version 4  ***/

/*  
 *  Routines to set up and save the various types of   
 *  two and three point functions
*
 *
 *
*/

/* Modifications
   C. McNeile 1997 Original version
   C. DeTar 5/24/97  Counting and naming of 3 pt form factors 
                     now provides for SW rotations
*/


#include "hqet_light_includes.h"
#include <assert.h>
#include <string.h>
#include "opertypes.h"
#include "corrlist.h"

/* function prototypes */
void dump_3pt_text(complex *corr);
void dump_light_light_2pt_text(complex *corr) ;
void dump_hqet_light_2pt_text(complex *corr) ;
/***********************/

#ifdef DEBUGDEF
#include DEBUGDEF
#endif


/*
 *  Set up the memory for the heavy-light correlator
 *
 *
 */

void setup_hlcorr(complex **hl_corr, int *hl_corr_dim )
{
  int i ; 
  int v_pt = novel - 1 ; 
  int no_oper = MAX_THREEPT*2 - 1 ; 
  int q_pt = no_q_values - 1 ;
  int k_spectator = no_spectator - 1 , k_zonked_light = no_zonked_light - 1 ;


  *hl_corr_dim = HQET_FORM_WHERE(nt,k_zonked_light,k_spectator,q_pt , v_pt ,no_oper ) ;

  if( ( *hl_corr = (complex *) calloc( (size_t) *hl_corr_dim, sizeof(complex) )  ) == NULL )
  {
    printf("node[%d] ERROR: could not reserve buffer space for the H-L correlators\n",this_node);
    terminate(1);
  }


  for(i = 0 ; i  < *hl_corr_dim ; ++i)
  {
    (*hl_corr + i)->real = 0.0 ; 
    (*hl_corr + i)->imag = 0.0 ; 
  }




}




/*
 *  A routine to finish up the calculation of the 
 *  heavy-light form factors
 *
 */

void finish_hlcorr(complex *hl_corr, int hl_corr_dim )
{
  int i ; 

  /***  Sum the heavy-light correlators over the nodes ****/
  for(i = 0 ; i  < hl_corr_dim ; ++i)
    g_complexsum(hl_corr + i ) ;


  /*** write the correlators to disk ****/
  IF_MASTER
  {
    write_hqet_form_corr(hl_corr, heavy_light_out, hl_corr_dim) ;
    dump_3pt_text(hl_corr  ) ; 
  }

}




/*
 *  Set up the memory for the required two point
 *  functions.
 *
 */

void setup_twocorr(complex **two_corr, int *two_corr_dim )
{
  int i ; 
  int no_oper = MAX_TWOPT*2 - 1 ; 
  int q_pt = no_q_values - 1 ;
  int k_spectator = no_spectator - 1 , k_zonked_light = no_zonked_light - 1 ;


  *two_corr_dim = TWOPT_FORM_WHERE(nt,k_zonked_light,k_spectator,q_pt  ,no_oper ) ;

  if( ( *two_corr = (complex *) calloc( (size_t) *two_corr_dim, sizeof(complex) )  ) == NULL )
  {
    printf("node[%d] ERROR: could not reserve buffer space for the two point correlators\n",this_node);
    terminate(1);
  }


  for(i = 0 ; i  < *two_corr_dim ; ++i)
  {
    (*two_corr + i)->real = 0.0 ; 
    (*two_corr + i)->imag = 0.0 ; 
  }

}






/*
 *  A routine to finish up the calculation of the 
 *  two point functions
 *
 */

void finish_twocorr(complex *two_corr, int two_corr_dim )
{
  int i ; 

  /***  Sum the heavy-light correlators over the nodes ****/
  for(i = 0 ; i  < two_corr_dim ; ++i)
    g_complexsum(two_corr + i ) ;


  /*** write the correlators to disk ****/
  IF_MASTER
  {
    write_twopt(two_corr, twopt_out, two_corr_dim) ;  
    dump_light_light_2pt_text(two_corr ) ;

  }

}




/*
 *  Set up the memory for the correlators from the 
 *  sequential source corerelators.
 *
 */

void setup_seq_corr(complex **seq_corr, int *seq_corr_dim )
{
  int i ; 
  int v_pt = novel - 1 ; 
  int k_spectator = no_spectator - 1 ; 


  *seq_corr_dim =   TWOPT_SEQ_WHERE(nt, k_spectator, v_pt ) ;

  if( ( *seq_corr = (complex *) calloc( (size_t) *seq_corr_dim, sizeof(complex) )  ) == NULL )
  {
    printf("node[%d] ERROR: could not reserve buffer space for the sequential correlators\n",this_node);
    terminate(1);
  }


  for(i = 0 ; i  < *seq_corr_dim ; ++i)
  {
    (*seq_corr + i)->real = 0.0 ; 
    (*seq_corr + i)->imag = 0.0 ; 
  }




}




/*
 *  A routine to finish up the calculation of the 
 *  sequential two point functions.
 *
 */

void finish_seqcorr(complex *seq_corr, int seq_corr_dim )
{
  int i ; 

  /***  Sum the heavy-light correlators over the nodes ****/
  for(i = 0 ; i  < seq_corr_dim ; ++i)
    g_complexsum(seq_corr + i ) ;


  /*** write the correlators to disk ****/
  IF_MASTER
  {
    write_seq_twopt(seq_corr, seq_out , seq_corr_dim) ;  
    dump_hqet_light_2pt_text( seq_corr); 
  }



}





/* 
 *  Return the name of the three point function
 *
 */



char *three_oper_name(int n )
{

  static char opername[MAXNAME];
  register int i;


  if(n < MAX_THREEPT)
    return name_3pt[hqet_to_light[n].oper] ;

  else
    {
      /* Modified name for 3 pt function with rotated propagator */
      /* Appends suffix to specify rotation */

      for(i=0; i<MAXNAME; i++)opername[i] = '\0';
      strncpy(opername,name_3pt[hqet_to_light[n-MAX_THREEPT].oper],MAXNAME-1);
      strncat(opername,SW_suffix,MAXNAME-strlen(opername)-1);
      return opername;
    }
}



/*
 *
 *  Dump the heavy --> light hqet three point function to the
 *  terminal screen.
 *
 *  This is scalar code. It must be run on node 0.
 */



void dump_3pt_text(complex *corr)
{
  int where ;
  int t,zonk_pt,spect_pt,q_pt,vel_pt,oper_pt ;
  const int no_oper = MAX_THREEPT*2 ; 

  printf("\nHQET-->light three point functions\n\n") ;

  for(zonk_pt = 0 ; zonk_pt < no_zonked_light ; ++zonk_pt )
  {
    printf("zonked quark kappa = %f\n",  kappa_zonked_light[zonk_pt]);
    for(spect_pt = 0 ; spect_pt < no_spectator ; ++spect_pt )
    {
      printf("spectator quark kappa = %f\n", kappa_spectator[spect_pt]);
      for(q_pt = 0 ; q_pt < no_q_values ; ++q_pt ) 
      {
	printf("Momentum %d %d %d\n", q_momstore[q_pt][0],q_momstore[q_pt][1], q_momstore[q_pt][2]  ); 
	for(vel_pt = 0 ; vel_pt < novel ; ++vel_pt)
	{
	  printf("velocity xyzt = %f %f %f    %f\n", 
		 velocity[vel_pt][XUP],velocity[vel_pt][YUP],velocity[vel_pt][ZUP],
		 velocity[vel_pt][TUP]) ; 
	  for(oper_pt = 0 ; oper_pt < no_oper ; ++oper_pt)
	  {
	    printf("Operator pointer = %d  name = %s\n",oper_pt, 
		   three_oper_name(oper_pt) );
	    for(t = 0 ; t < nt ; ++t)
	    {
	      where = HQET_FORM_WHERE(t,zonk_pt,spect_pt,q_pt,vel_pt,oper_pt) ;
	      printf("%d    %g  %g\n",t , (corr+where)->real, (corr+where)->imag);
	    }


	  }
	}
      }
    }
  }




}  /*** end of the text dump of the hqet--light three point function ***/






/* 
 *  Return the name of the two point function
 *
 */


char *oper_name(int n )
{
  /*  When the global array "twopt_opertype" is updated, 
   *  this routine must also up updated.
   */

  assert( n >= 0 && n < NO_TWOPT_OPERS ); 

  return twopt_name[n] ; 

}



/*
 *
 *  Dump the light--light two point function to the 
 *  terminal screen.
 *
 *
 *
 */



void dump_light_light_2pt_text(complex *corr)
{
  int where ;
  int t,zonk_pt,spect_pt,q_pt, oper_pt ;
  const int no_oper = MAX_TWOPT ; 
  

  printf("\nHere are the light--light two point functions\n\n") ; 


  for(zonk_pt = 0 ; zonk_pt < no_zonked_light ; ++zonk_pt )
  {
    printf("zonked quark kappa = %f\n", kappa_zonked_light[zonk_pt] );
    for(spect_pt = 0 ; spect_pt < no_spectator ; ++spect_pt )
    {
      printf("spectator quark kappa = %f\n",kappa_spectator[spect_pt]);
      for(q_pt = 0 ; q_pt < no_q_values ; ++q_pt ) 
      {
	printf("Momentum %d %d %d\n", q_momstore[q_pt][0],q_momstore[q_pt][1], q_momstore[q_pt][2]  ); 
	for(oper_pt = 0 ; oper_pt < no_oper ; ++oper_pt)
	{
	  printf("Operator pointer = %d oper = %s\n",oper_pt, oper_name(oper_pt) );
	  for(t = 0 ; t < nt ; ++t)
	  {
	    where = TWOPT_FORM_WHERE(t,zonk_pt,spect_pt,q_pt,oper_pt) ;
	    printf("%d   %g  %g \n",t , (corr+where)->real, (corr+where)->imag);
	  }
	  
	  
	}
      }
    }
  }
  
  


} /*** end of the TEXT dump of the light-light 2pt functions ****/


/*
 *
 *  Dump the heavy --> light hqet sequential two point functions 
 *  to the terminal screen.
 *
 *
 *
 */



void dump_hqet_light_2pt_text(complex *corr)
{
  int where ;
  int t,  spect_pt, vel_pt ; 
  

  printf("Here are the HQET-spectatator_quark two point function \n\n") ; 


  for(spect_pt = 0 ; spect_pt < no_spectator ; ++spect_pt )
  {
    printf("spectator quark kappa = %f\n", kappa_spectator[spect_pt]);
    for(vel_pt = 0 ; vel_pt < novel ; ++vel_pt ) 
    {
      printf("velocity xyzt = %f %f %f    %f\n", 
	     velocity[vel_pt][XUP],velocity[vel_pt][YUP],velocity[vel_pt][ZUP],velocity[vel_pt][TUP]) ; 
      for(t = 0 ; t < nt ; ++t)
      {
	where = TWOPT_SEQ_WHERE(t, spect_pt, vel_pt) ; 
	printf("%d    %g  %g\n",t , (corr+where)->real, (corr+where)->imag);
      }

    }
  }
  
  


}

