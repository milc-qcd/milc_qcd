/*********************** setup_smearing.c *************************/
/* MIMD version 7 */
/*
 **   Load in the parameters for the smaearing functions from
 **   standard input, and create the smearing function
 **/

#include "w_static_includes.h"
#include <string.h>

/****
    This functions reads the options in the standrad input for the type
    of smearing functions to be used. It stores this information in the array
    "smear_code", which is eventually sent to all the nodes. ((This routine only
     runs on the master node))

*****/
int get_smearing_funcs_code(char savebuf[], params *par_buf )
{
  int ismear ; 
  int status = 0;
  Real decay, A, B ,D ; 

  /** 
    load in the names of the file to read the smearing functions from 
    ***/

  for(ismear = 0 ; ismear < par_buf->nosmear ; ++ismear)
  {
    if( scanf("%s",savebuf) != 1 )
    {
      printf("error in input: smear_func required\n"); 
      return status++;
    }
    if(strcmp("smear_func",savebuf) != 0 ) 
    {
      printf("error in input: smear_func required\n"); 
      return status++;
    }


    /*** now choose the type of smearing function to be  created ***/
    if( scanf("%s",savebuf) != 1 )
    {
      printf("error in input: smearing type required (eg local, wall, ....)\n"); 
      return status++;
    }
    

    if(strcmp("local",savebuf) == 0 ) 
    {
      par_buf->smear_code[ismear][0] = LOCALsmear ;
    }
    else if(strcmp("wall",savebuf) == 0 ) 
    {
      par_buf->smear_code[ismear][0] = WALLsmear ;
    }
    else if(strcmp("expon",savebuf) == 0 ) 
    {
      read_expon_params( &decay, savebuf); 
      par_buf->smear_code[ismear][0] = EXPON ;
      par_buf->smear_code[ismear][1] = decay ; 
    }
    else if(strcmp("2S_smear",savebuf) == 0 ) 
    {
      read_twoS_params(&decay, &A, savebuf );
      par_buf->smear_code[ismear][0] = twoSpoly ; 
      par_buf->smear_code[ismear][1] = decay ; 
      par_buf->smear_code[ismear][2] = A ;

    }
    else if(strcmp("3S_smear",savebuf) == 0 ) 
    {
      read_threeS_params(&decay, &B, &D,  savebuf );
      par_buf->smear_code[ismear][0] = threeSpoly ; 
      par_buf->smear_code[ismear][1] = decay ; 
      par_buf->smear_code[ismear][2] = B ;
      par_buf->smear_code[ismear][3] = D ;
    }
    else
    {
      printf("ERROR: I do not recognize the smearing type %s \n",savebuf) ; 
      return status++;
    }

    printf("The smearing code is %s\n",savebuf);

  } /*** end of the loop over the smearing functions ***/


/*******  DEBUG ****DEBUG ****DEBUG ****DEBUG ****DEBUG ****
  for(ismear = 0 ; ismear < par_buf->nosmear ; ++ismear )
  {
    printf("DEBUG_smear %g %g %g %g\n",par_buf->smear_code[ismear][0],par_buf->smear_code[ismear][1],
	   par_buf->smear_code[ismear][2],par_buf->smear_code[ismear][3] ) ;

  }
  exit(1); 

*******  DEBUG ****DEBUG ****DEBUG ****DEBUG ****DEBUG ****/


  return status;

} /*** end of the get_smearing_funcs_code function   ****/







/***
     Read the parameters for the exponential smearing function
****/


void read_expon_params(Real *decay, char savebuf[]  )
{

#if PRECISION == 1
  if( scanf("%f",decay) != 1 )
#else
  if( scanf("%lf",decay) != 1 )
#endif
  {
    printf("ERROR:: The decay parameter for the smearing function is required \n");
    terminate(1); 
  }


  if( scanf("%s",savebuf) != 1 )
  {
    printf("error in input: \"end_param\" required\n"); 
    terminate(1);
  }

  if(strcmp("end_param",savebuf) != 0 ) 
  {
    printf("error in input: \"end_param\" required\n"); 
    terminate(1);
  }


}






/***
     Read the parameters for the 2S smearing function
****/


void read_twoS_params(Real *decay, Real *A, char savebuf[]  )
{

#if PRECISION == 1
  if( scanf("%f",decay) != 1 )
#else
  if( scanf("%lf",decay) != 1 )
#endif
  {
    printf("ERROR:: The decay parameter for the 2S smearing function is required \n");
    terminate(1); 
  }


#if PRECISION == 1
  if( scanf("%f",A) != 1 )
#else
  if( scanf("%lf",A) != 1 )
#endif
  {
    printf("ERROR:: The linear parameter for the 2S smearing function is required \n");
    terminate(1); 
  }


  if( scanf("%s",savebuf) != 1 )
  {
    printf("error in input: \"end_param\" required\n"); 
    terminate(1);
  }

  if(strcmp("end_param",savebuf) != 0 ) 
  {
    printf("error in input: \"end_param\" required\n"); 
    terminate(1);
  }


}






/***
     Read the parameters for the 3S smearing function
****/


void read_threeS_params(Real *decay, Real *B, Real *D, char savebuf[]  )
{

#if PRECISION == 1
  if( scanf("%f",decay) != 1 )
#else
  if( scanf("%lf",decay) != 1 )
#endif
  {
    printf("ERROR:: The decay parameter for the 2S smearing function is required \n");
    terminate(1); 
  }


#if PRECISION == 1
  if( scanf("%f",B) != 1 )
#else
  if( scanf("%lf",B) != 1 )
#endif
  {
    printf("ERROR:: The linear parameter for the 3S smearing function is required \n");
    terminate(1); 
  }

#if PRECISION == 1
  if( scanf("%f",D) != 1 )
#else
  if( scanf("%lf",D) != 1 )
#endif
  {
    printf("ERROR:: The quadratic parameter for the 3S smearing function is required \n");
    terminate(1); 
  }


  if( scanf("%s",savebuf) != 1 )
  {
    printf("error in input: \"end_param\" required\n"); 
    terminate(1);
  }

  if(strcmp("end_param",savebuf) != 0 ) 
  {
    printf("error in input: \"end_param\" required\n"); 
    terminate(1);
  }


}




/*
 *  Construct the smearing functions using the 
 *  information stored in the "smear_code" array.
 *
 */ 

void create_smearing_funcs(void)
{
  int ismear ; 
  Real decay, A, B ,D ; 
  int which_smear ; 

  /** 
    load in the names of the file to read the smearing functions from 
    ***/

  for(ismear = 0 ; ismear < nosmear ; ++ismear)
  {
    which_smear = smear_code[ismear][0]   ;

    

    if( which_smear ==  LOCALsmear   )
    {
      create_local_oper(ismear, smearfile_in[ismear] );
    }
    else if(which_smear ==  WALLsmear   )
    {
      create_wall_oper(ismear, smearfile_in[ismear] );
    }
    else if( which_smear == EXPON )
    {
      decay = smear_code[ismear][1]  ; 
      create_expon_oper(ismear, decay,smearfile_in[ismear] );
    }
    else if ( which_smear == twoSpoly  )
    {
      decay  = smear_code[ismear][1]  ; 
      A      = smear_code[ismear][2]  ;

      create_2S_polyexpon_oper(ismear, decay,A,smearfile_in[ismear] );
    }
    else if(which_smear == threeSpoly  )
    {
      decay =  smear_code[ismear][1]  ; 
      B     =  smear_code[ismear][2] ;
      D     =  smear_code[ismear][3] ;

      create_3S_polyexpon_oper(ismear, decay,B,D,smearfile_in[ismear] );
    }
    else
    {
      printf("ERROR: I do not recognize the smearing type \n") ; 
      terminate(1);
    }

    IF_MASTER
      printf("The smearing function name is %s\n",smearfile_in[ismear]);

  } /*** end of the loop over the smearing functions ***/



} /*** end of the create_smearing_funcs function   ****/


