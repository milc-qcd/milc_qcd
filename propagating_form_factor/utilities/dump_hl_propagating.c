/*
 *
 *  Dump the heavy --> light hqet from factor to the
 *  terminal screen.
 *
 *
 *
 */


#include"read_hl_form.h"

char *three_oper_name(int n , int copy_pt) ;


void dump_hl_prop_form(complex *corr,  int corr_oper_list[], 
      int corr_copy_list[],
      int no_zonked_light, int no_spectator, int no_sequential,
      int no_q_values, int no_p_values, int no_oper, int nt, int nocopies)
{
  int where ;
  int t = nt ;
  int zonk_pt  = no_zonked_light - 1 ; 
  int spect_pt = no_spectator -1 ;
  int q_pt = no_q_values -1 ;
  int p_pt = no_p_values -1 ;
  int oper_pt = no_oper - 1 ;
  int seq_pt = no_sequential - 1 ; 

  int copy_pt ;
  int corr_stride ; 
  /**************************************************/
  corr_stride = LIGHT_FORM_WHERE(t,zonk_pt,seq_pt,spect_pt,q_pt,p_pt,oper_pt) ; 


  for(zonk_pt = 0 ; zonk_pt < no_zonked_light ; ++zonk_pt )
  {
    printf("zonked quark = %d\n",zonk_pt);
    for(spect_pt = 0 ; spect_pt < no_spectator ; ++spect_pt )
    {
      printf("  spectator quark = %d\n",spect_pt);
	for(seq_pt = 0 ; seq_pt <no_sequential  ; ++seq_pt)
	{
	  printf("  sequential quark = %d\n",seq_pt);
	  for(q_pt = 0 ; q_pt < no_q_values ; ++q_pt ) 
	  {
	    printf("   Q--Momentum pointer = %d \n",q_pt  ); 
	    
	    for(p_pt = 0 ; p_pt < no_p_values ; ++p_pt)
	    {
	      printf("   P--Momentum pointer = %d \n",q_pt  ); 
	      for( copy_pt = 0 ; copy_pt < nocopies ; ++copy_pt ) 
	      for(oper_pt = 0 ; oper_pt < no_oper ; ++oper_pt)
	      {
		printf("       Operator pointer = %d  name = %s\n",
		       corr_oper_list[oper_pt], 
		       three_oper_name(corr_oper_list[oper_pt],
				       corr_copy_list[copy_pt]) );
		for(t = 0 ; t < nt ; ++t)
		{
		  where = corr_stride*copy_pt + LIGHT_FORM_WHERE(t,zonk_pt,seq_pt,spect_pt,q_pt,p_pt,oper_pt) ;
		  printf("         corr[ %d ] = ( %g , %g )\n",t , (corr+where)->real, (corr+where)->imag);
		}

		
	      }
	    }
	  }
	}
    }
  }



}


/*
 *   Dump the heavy to heavy form factors to 
 *   the screen.
 *
 *
 *
 */


void dump_hh_prop_form(complex *corr,  int corr_oper_list[],
       int corr_copy_list[],
       int no_zonked_heavy, int no_spectator, int no_sequential,
       int no_q_values, int no_p_values, int no_oper, int nt, int nocopies)
{
  int where ;

  int t = nt ;
  int zonk_pt  = no_zonked_heavy - 1 ; 
  int spect_pt = no_spectator -1 ;
  int q_pt = no_q_values -1 ;
  int p_pt = no_p_values -1 ;
  int oper_pt = no_oper - 1 ;
  int seq_pt = no_sequential - 1 ; 

  int copy_pt ; 
  int corr_stride ; 
  /************************************************************/
  corr_stride = HEAVY_FORM_WHERE(t,zonk_pt,seq_pt,spect_pt,q_pt,p_pt,oper_pt)  ; 

  for(zonk_pt = 0 ; zonk_pt < no_zonked_heavy ; ++zonk_pt )
  {
    printf("zonked quark = %d\n",zonk_pt);
    for(spect_pt = 0 ; spect_pt < no_spectator ; ++spect_pt )
    {
      printf("  spectator quark = %d\n",spect_pt);
	for(seq_pt = 0 ; seq_pt <no_sequential  ; ++seq_pt)
	{
	  printf("  sequential quark = %d\n",seq_pt);
	  for(q_pt = 0 ; q_pt < no_q_values ; ++q_pt ) 
	  {
	    printf("   Q--Momentum pointer = %d \n",q_pt  ); 
	    
	    for(p_pt = 0 ; p_pt < no_p_values ; ++p_pt)
	    {
	      printf("   P--Momentum pointer = %d \n",q_pt  ); 
	      for( copy_pt = 0 ; copy_pt < nocopies ; ++copy_pt ) 
		for(oper_pt = 0 ; oper_pt < no_oper ; ++oper_pt)
		{
		  printf("       Operator pointer = %d  name = %s\n",
			 corr_oper_list[oper_pt], 
			 three_oper_name(corr_oper_list[oper_pt],
					 corr_copy_list[copy_pt]) );
		  for(t = 0 ; t < nt ; ++t)
		  {
		    where = copy_pt*corr_stride  + HEAVY_FORM_WHERE(t,zonk_pt,seq_pt,spect_pt,q_pt,p_pt,oper_pt)  ;
		    printf("         corr[ %d ] = ( %g , %g )\n",t , (corr+where)->real, (corr+where)->imag);
		  }
		  
		  
		}
	    }
	  }
	}
    }
  }



}



#undef NONONON_DO_NOT_INCLUDE_THIS_CODE
#ifdef NONONON_DO_NOT_INCLUDE_THIS_CODE


/*
 *  Dump the difference of two heavy-light correlators
 *
 *
 */


void diff_hl_form(complex *corr_one,  complex *corr_two, int corr_oper_list[], 
 int corr_copy_list[],
int no_zonked_light, int no_spectator, int no_q_values, int novel,
int no_oper, int nt)
{
  int where ;
  int t,zonk_pt,spect_pt,q_pt,vel_pt,oper_pt ;
  complex diff ; 
  Real maxdiff = 0.0 ;
  Real diff_mag ; 

  printf("Here is the DIFFERENCE hqet, heavy-->light,  form factor \n") ; 


  for(zonk_pt = 0 ; zonk_pt < no_zonked_light ; ++zonk_pt )
  {
    printf("zonked quark = %d\n",zonk_pt);
    for(spect_pt = 0 ; spect_pt < no_spectator ; ++spect_pt )
    {
      printf("  spectator quark = %d\n",spect_pt);
      for(q_pt = 0 ; q_pt < no_q_values ; ++q_pt ) 
      {
	printf("    Momentum pointer = %d \n",q_pt  ); 
	for(vel_pt = 0 ; vel_pt < novel ; ++vel_pt)
	{
	  printf("      velocity pointer = %d\n",+vel_pt); 
	  for(oper_pt = 0 ; oper_pt < no_oper ; ++oper_pt)
	  {
	    printf("       Operator pointer = %d  name = %s\n",
		   corr_oper_list[oper_pt], 
		   three_oper_name(corr_oper_list[oper_pt]) );
	    for(t = 0 ; t < nt ; ++t)
	    {
	      where = HQET_FORM_WHERE(t,zonk_pt,spect_pt,q_pt,vel_pt,oper_pt) ;

	      diff.real = (corr_one+where)->real - (corr_two+where)->real  ;
	      diff.imag = (corr_one+where)->imag - (corr_two+where)->imag  ;

	      diff_mag = sqrt( diff.real*diff.real +  diff.imag*diff.imag)  ; 
	      if( diff_mag > maxdiff)   maxdiff = diff_mag ;
	      

	      printf("         corr[ %d ] :: ( %g , %g ) -  ( %g , %g ) =  ( %g , %g ) \n",t , 
		     (corr_one + where)->real, (corr_one + where)->imag , 
		     (corr_two + where)->real, (corr_two + where)->imag , 
		     diff.real, diff.imag );
	    }


	  }
	}
      }
    }
  }

  printf("\n>>> Maximum magnitude of the difference = %e\n",maxdiff); 


}


#endif
