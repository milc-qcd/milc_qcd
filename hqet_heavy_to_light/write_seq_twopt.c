/********************** write_seq_twopt.c ****************************/ 
/* MIMD version 6 */
/*  Write out the sequential two point functions to 
 *  a disk file.
 *
 *  NOTE this is scalar code -- it should only be run on 
 *  one node.
 *
 *  Subroutine argument
 *   filename  :: the name of the disk file
 *   corr      :: the two point  correlators
 *   dim       :: the amount of data
 */

/** $Header: /lqcdproj/detar/cvsroot/milc_qcd/hqet_heavy_to_light/write_seq_twopt.c,v 1.1 2005/02/23 00:05:07 detar Exp $ ***/

#include "hqet_light_includes.h"

#ifdef DEBUGDEF
#include DEBUGDEF
#endif



void write_seq_twopt(complex *corr, char filename[], int dim)
{
  FILE *fp ;
  size_t nobj = (size_t) dim ;
  /* Hack to distinguish single and double precision files */
  const int magic_number = 23991 + 8*(sizeof(Real) - 4) ; 
  const int version_number  = 1  ; /** update this flag when the data format changes ***/
  int i ;
  size_t  name_len = 80 ;
  size_t how_many ;
  int check_sum ; 

  /* Memory for some header information ***/
#define HEADER_DIM_WRITE_CORR 7
  size_t  header_size = (size_t) HEADER_DIM_WRITE_CORR ;
  int header_data[HEADER_DIM_WRITE_CORR] ;


  /*** calculate the checksum of the form factors ****/
  check_sum =  bsd_sum((char*)corr, sizeof(complex)*dim)  ;

  /***** load up the header *****/
  header_data[0] = magic_number ;
  header_data[1] = version_number ;
  header_data[2] = nt ;
  header_data[3] = check_sum ;
  header_data[4] = novel ; 
  header_data[5] = no_spectator  ;
  header_data[6] = dim ;

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



  /*** write the velocities out to the file ***/
  how_many = (size_t) 4*novel ;
  if( fwrite(velocity,sizeof(Real),how_many,fp) != how_many   )
  {
    printf("There was an error during the writing of the velocity table \n");
    terminate(1);
  }

  /** write the three point function to disk *****/
  if( fwrite(corr,sizeof(complex),nobj,fp) != nobj   )
  {
    printf("There was an error during the writing of the form factors\n");
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


