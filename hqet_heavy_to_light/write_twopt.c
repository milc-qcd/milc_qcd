/************* write_twopt.c **************************/
/* MIMD version 4 ***/

/* 
 *  Write out the hqet heavy---> light form factors to 
 *  a disk file.
 *
 *  NOTE this is scalar code -- it should only be run on 
 *  one node.
 *
 *  Subroutine argument
 *   filename  :: the name of the disk file
 *   corr      :: the hqet form factors correlators
 *   dim       :: the amount of data
 */

/* Modifications
   C. McNeile  1997 Original version magic number was 23631
   C. DeTar 5/24/97 Now dumps character names of two-pt propagators -
                    changed magic number to 23731
*/

/** $Header: /lqcdproj/detar/cvsroot/milc_qcd/hqet_heavy_to_light/write_twopt.c,v 1.1 2005/02/23 00:05:07 detar Exp $ ***/

#include "hqet_light_includes.h"
#include "opertypes.h"
#include "corrlist.h"
#include <string.h>

#ifdef DEBUGDEF
#include DEBUGDEF
#endif


void write_twopt(complex *corr, char filename[], int dim)
{
  FILE *fp ;
  size_t nobj = (size_t) dim ;
  /* Hack to distinguish single and double precision files */
  const int magic_number = 23631 + 8*(sizeof(Real) - 4) ; 
  const int version_number  = 1  ; /** update this flag whwn the data format changes ***/
  int i ;
  size_t  name_len = 80 ;
  size_t how_many ;
  int check_sum ; 
  char out[ MAXNAME ] ; 


  /* Memory for some header information ***/
#define HEADER_DIM_WRITE_CORR 11
  size_t  header_size = (size_t) HEADER_DIM_WRITE_CORR ;
  int header_data[HEADER_DIM_WRITE_CORR] ;
  char *oper_name(int n );


  /*** calculate the checksum of the form factors ****/
  check_sum =  bsd_sum((char*)corr, sizeof(complex)*dim)  ;

  /***** load up the header *****/
  header_data[0] = magic_number ;
  header_data[1] = version_number ;
  header_data[2] = nt ;
  header_data[3] = check_sum ;
  header_data[4] = no_q_values ;
  header_data[5] = MAX_TWOPT  ;
  header_data[6] = no_spectator  ;
  header_data[7] = no_zonked_light  ;
  header_data[8] = dim ;
  header_data[9] = MAXNAME;
  header_data[10] = tf;

  /**** open the file ******/
  if( (fp = fopen(filename ,"wb")) == NULL )
  {
    printf("Could not open the file %s\n",filename);
    terminate(1);
  }

  /** write the header information ***/
  if( fwrite(header_data,sizeof(int),header_size,fp) != header_size )
  {
    printf("There was an error during the writing of the form factor HEADER \n");
    terminate(1);
  }



  /*** write the external momentum to file ***/
  how_many = (size_t) 3*no_q_values ;
  if( fwrite(q_momstore  ,sizeof(int),how_many,fp) != how_many   )
  {
    printf("There was an error during the writing of the momentum table n");
    terminate(1);
  }


  /*** write the names of the correlators to file ***/

  for(i=0; i<MAX_TWOPT; i++)
    {
      strncpy(out, oper_name(i) , MAXNAME ) ;
      if( fwrite(out, MAXNAME,1,fp) != 1)
	{
	  printf("There was an error during the writing of the 2 pt op names \n");
	  terminate(1);
	}
    }

  /** write the two point function to disk *****/
  if( fwrite(corr,sizeof(complex),nobj,fp) != nobj   )
  {
    printf("There was an error during the writing of the two point names\n");
    terminate(1);
  }


  /*** close the file ****/
  if( fclose(fp) != 0 )
  {
    printf("There was an error during the closing of %s \n",filename);
    terminate(1);
  }


  printf("I have written some form factors to the file %s\n",filename);



}


