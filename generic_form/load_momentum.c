/**************************** load_momentum.c ****************************/
/* MIMD version 7 */
/*
 *  Scalar code to load in the spatial momentum from a disk file.
 *
 *  The total number of momentum components is returned
 */

#include "generic_form_includes.h"
#include <string.h>

#ifndef IF_OK
#define IF_OK if(status==0)
#endif

int load_momentum_from_disk(int mom_in[][3], char filename[], int max_mom)
{
  int countmom  = -1 ;
  int i , j ;
  FILE *fp;
  int v1,v2,v3 ;
  /***--------------------------------------------------***/
  for(i=0 ; i < max_mom ; ++i)
    for(j=0; j < 4 ; ++j)
      mom_in[i][j] = 0  ;

  /*** open the file *****/
  if( (fp = fopen(filename ,"r")) == NULL )
  {
    printf("ERROR::load_momentum_from_disk::Could not open the file %s\n",filename);
    terminate(1);
  }


  /*** read in the momentum and compute the fourth component ****/
  while( fscanf(fp,"%d %d %d",&v1,&v2,&v3) == 3 && countmom < max_mom )
  {
    ++countmom ;
    mom_in[countmom][XUP ] = v1 ; 
    mom_in[countmom][YUP ] = v2 ; 
    mom_in[countmom][ZUP ] = v3 ; 

  }


  if( countmom < 0  )
  {
    printf("ERROR: No momenta have been read from %s\n",filename);
    terminate(1); 
  }
  else if( countmom >= max_mom )
  {
    printf("ERROR: Not enough space for the momentum in %s \n",filename);
    terminate(1); 
  }

  ++countmom ;

  /*** close the file ****/
  if( fclose(fp) != 0 )
  {
    printf("There was an error during the closing of %s \n",filename);
    terminate(1);
  }

  printf("The number of momentum read in = %d\n",countmom);
  printf("p_x      p_y     p_z\n");
  for(i=0 ; i < countmom ; ++i)
  {
    printf("%d    %d    %d\n", mom_in[i][0],mom_in[i][1],mom_in[i][2] );
  }


  return countmom ;
}




/*
 *  Routine to load in the momentum values
 *  from the input file.
 *
 *  Subroutine arguments
 *    prompt  :: whether to print stuff before the input is required
 *    no_mom  :: return the number of momentum values to run
 *    mom_in  :: array containing the momentum (as integers)
 *    max_mom :: the maxium number of momentum values allowed by the
 *               code
 */


int load_momentum(int prompt, char *label, int *no_mom, int mom_in[][3], int max_mom)
{
  int status = 0 ; 
  char savebuf[128] ;
  char checklabel[4] ;
  int what_to_do ; 
  enum what_to_do_choices { read_momentum = 50 ,  stop_reading_momentum  } ; 
  int n_mom ; 

  /***------------------------------------------------------------****/

  /*** read in all the momentum that the user wants to calculate ***/

  IF_OK if (prompt==1) 
    printf("enter 'start_of_momentum %s' \n",label);
  IF_OK scanf("%s",savebuf);
  IF_OK scanf("%s",checklabel);
  IF_OK printf("%s ",savebuf);
  IF_OK printf("%s\n",checklabel);
  IF_OK {
    if(strcmp("start_of_momentum",savebuf) == 0 
       && strcmp(label,checklabel) == 0)
    {
      what_to_do = read_momentum ;
    }
    else
    {
      printf("error in input: start_of_momentum %s required, but got %s %s\n",
	     label,savebuf,checklabel);
      status++ ;
      what_to_do = stop_reading_momentum ; 
    }
  }



  n_mom = 0 ; 

  /*** now read the momentum to compute ***/
  while( what_to_do == read_momentum   ) 
  {

    IF_OK if (prompt==1) 
      printf("enter three integers for the momentum\n");
    IF_OK scanf("%s",savebuf);

    if(strcmp("end_momentum",savebuf) == 0 )
    {
      what_to_do = stop_reading_momentum ; 
    }
    else if( n_mom >= max_mom )
    {
      ++status ; 
      what_to_do = stop_reading_momentum ; 
      printf("Too many momentum values (maximum = %d)\n",max_mom) ; 
    }
    else
    {
      mom_in[n_mom][0]  = atoi(savebuf) ; 

      if( scanf("%d %d",&mom_in[n_mom][1] ,&mom_in[n_mom][2] ) == 2 )
      {
	printf("%s= %d   %d   %d\n",label,mom_in[n_mom][0], mom_in[n_mom][1], mom_in[n_mom][2] ) ; 
	++n_mom ; 
      }
      else
      {
	++status ; 
	what_to_do = stop_reading_momentum ; 
	printf("Error reading momentum values\n"); 

      }


    }  


  }  /*** end the loop over the momentum  (while loop) ***/





  *no_mom = n_mom ; 


  return status ; 
}





