/**************************** load_velocity.c ****************************/
/* MIMD version 6 */
/*
 *  Scalar code to load in the spatial velocities from a disk file.
 *  The fourth component of the velocity is calculated.
 *
 *  The total number of velocity components is returned
 *
 *  The smearing functions for the load are also read in.
 *  Format of the input file:
 *      0.0 0.0  0.0    zero_smear
 *      0.1 0.0  0.0    one_smear
 *      0.3 0.0  0.0    two_smear
 *
 * 
 */

#include "hqet_light_includes.h"
#include <string.h>

int load_velocity_from_disk(params *par_buf, char filename[], int maxvel)
{
  int countvel  = -1 ;
  int i , j ;
  FILE *fp;
  Real v1,v2,v3 ;
  char buff[MAXFILENAME] ;

  /***--------------------------------------------------***/
  for(i=0 ; i < maxvel ; ++i)
    for(j=0; j < 4 ; ++j)
      par_buf->velocity[i][j] = 0.0 ;

  /*** open the file *****/
  if( (fp = fopen(filename ,"r")) == NULL )
  {
    printf("Could not open the file %s\n",filename);
    terminate(1);
  }


  /*** read in the velocities and compute the fourth component ****/
  while( fscanf(fp,"%lfHELP %lfHELP %lfHELP",&v1,&v2,&v3) == 3 && countvel < maxvel )
  {
    ++countvel ;
    par_buf->velocity[countvel][XUP ] = v1 ; 
    par_buf->velocity[countvel][YUP ] = v2 ; 
    par_buf->velocity[countvel][ZUP ] = v3 ; 

    par_buf->velocity[countvel][TUP] = 
      (Real) sqrt( (double)  1.0 + v1*v1 + v2*v2 + v3*v3  );


    if(  fscanf(fp,"%s %s",buff,par_buf->hqet_smear_file[countvel] ) != 2  )
    {
      printf("Error reading the velocity filename\n");
      terminate(1);
    }

    if(strcmp("smear_file",buff) != 0 ) 
    {
      printf("Error in the file %s containg the velocites, \"smear_file\" required \n",filename);
      terminate(1);
    }





  } 


  if( countvel < 0  )
  {
    printf("ERROR: No velocities have been read from %s\n",filename);
    terminate(1); 
  }
  else if( countvel >= maxvel )
  {
    printf("ERROR: Not enough space for the velocities in %s \n",filename);
    terminate(1); 
  }

  ++countvel ;

  /*** close the file ****/
  if( fclose(fp) != 0 )
  {
    printf("There was an error during the closing of %s \n",filename);
    terminate(1);
  }

  printf("The number of velocities read in = %d\n",countvel);
  printf("v0      v1     v2    v3   v4\n");
  for(i=0 ; i < countvel ; ++i)
  {
    printf("%f    %f  %f  %f\n", 
	   par_buf->velocity[i][0],
	   par_buf->velocity[i][1],
	   par_buf->velocity[i][2],
	   par_buf->velocity[i][3] );
  }


  printf("Here are the names of the files to read the smearing functions from \n");
  for(i=0 ; i < countvel ; ++i)
  {
    printf("vel_smear[%d] = %s\n",i,par_buf->hqet_smear_file[ i ] );
  }

  return countvel ;
}


