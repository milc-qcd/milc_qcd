/****************************** write_propagating_two_pt.c ****************/
/* MIMD version 6 */
/* NOT MAINTAINED.  TEST BEFORE USE. Consider using write_propagating_two_point.c */
/* 
 *  Write out the propagating two point functions
 *  a disk file.
 *
 *
 *
 *  NOTE this is scalar code -- it should only be run on 
 *  one node.
 *
 *  Subroutine argument
 *   filename  :: the name of the disk file
 *   corr      :: the two point correlators
 *   no_k_one :: the number of  kappa values
 *   no_k_two :: the number of  kappa values
 *   hl_flag   :: flag saying where the heavy-heavy or  heavy-light correlators
 *                are to be written.
 *   numer_of_operators  :: 
 *   dim       :: the amount of data

 */

/** $Header: /lqcdproj/detar/cvsroot/milc_qcd/propagating_form_factor/write_propagating_two_pt.c,v 1.1 2005/02/23 00:05:53 detar Exp $ ***/

#include "prop_form_includes.h"

#ifdef DEBUGDEF
#include "prop_form.h"
#endif



void write_prop_twopt(complex *corr, char filename[], 
int hl_flag, int no_k_one , int no_k_two, int numer_of_operators, int no_copies,int dim)
{
  FILE *fp ;
  size_t nobj = (size_t) dim ;
  /* Hack to distinguish single and double precision files */
  const int magic_number = 44451189 + 8*(sizeof(Real) - 4)); 
  const int version_number  = 1  ; /** update this flag when the data format changes ***/
  size_t writecount ;
  int check_sum ; 
  int32type q_momstore_local[MAXMOM][3] ;
  int i , p ; 


  /* Memory for some header information ***/
#define HEADER_DIM_WRITE_CORR 11
  size_t  header_size = (size_t) HEADER_DIM_WRITE_CORR ;
  int32type header_data[HEADER_DIM_WRITE_CORR] ;
#undef HEADER_DIM_WRITE_CORR

  /*** calculate the checksum of the two point function ****/
  check_sum =  bsd_sum((char*)corr, sizeof(complex)*dim)  ;


  /***** load up the header *****/
  header_data[0] = (int32type) magic_number ;
  header_data[1] = (int32type) version_number ;
  header_data[2] = (int32type) nt ;
  header_data[3] = (int32type) check_sum ;
  header_data[4] = (int32type) no_k_values ;
  header_data[5] = (int32type) no_k_one ;
  header_data[6] = (int32type) no_k_two ;
  header_data[7] = (int32type) numer_of_operators ;
  header_data[8] = (int32type) dim ;
  header_data[9] = (int32type) hl_flag ;
  header_data[10] = (int32type) no_copies ;


  /**** open the file ******/
  if( (fp = fopen(filename ,"wb")) == NULL )
  {
    printf("ERROR::write_prop_twopt::Could not open the file %s\n",filename);
    terminate(1);
  }

  /** write the header information ***/
  if( fwrite(header_data,sizeof(int32type),header_size,fp) != header_size )
  {
    printf("There was an error during the writing of the two point HEADER to %s\n",
	   filename);
    terminate(1);
  }


  /*** write the external k--momentum to file ***/
  for(i = 0 ; i < 3 ; ++i)
    for(p = 0 ; p < no_k_values ; ++p )
      q_momstore_local[p][i] =  k_momstore[p][i] ; 


  writecount = (size_t) 3*no_k_values ;
  if( fwrite(q_momstore_local  ,sizeof(int32type),writecount,fp) != writecount   )
  {
    printf("There was an error during the writing of the Q momentum table to %s\n",
	   filename);
    terminate(1);
  }



  /** write the two point function to disk *****/
  if( fwrite(corr,sizeof(complex),nobj,fp) != nobj   )
  {
    printf("There was an error during the writing of the two point functions to %s\n",
	   filename);
    terminate(1);
  }


  /*** close the file ****/
  if( fclose(fp) != 0 )
  {
    printf("There was an error during the closing of %s \n",filename);
    terminate(1);
  }


  printf("I have written the two point functions to the file %s\n",filename);

}


