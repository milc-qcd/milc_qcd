/* 
 *  Write the propagating two point functions
 *  to a disk file.
 *
 *  NOTE this is scalar code -- it should only be run on 
 *  one node.
 *
 *  Subroutine argument
 *   filename  :: the name of the disk file
 *   corr      :: the two point correlators
 *   no_zonked :: the number of  kappa values
 *   no_spectator :: the number of  kappa values
 *   hl_flag   :: flag saying where the heavy-heavy or  heavy-light correlators
 *                are to be written.
 *   number_of_operators  :: 
 *   dim       :: the amount of data

 */


#include "read_hl_form.h"
#include "prop_form_utilities.h"



/* Functions declarations used in the routines in this file ***/

void write_prop_twopt(
  complex *corr,            /* the two point correlators */
  char filename[80],        /* the name of the disk file */
  int32type *q_momstore,       /* list of momenta */
  int hl_flag,              /* flag identifying type of correlator */
  int no_spectator,         /* the number of  kappa values */
  int no_zonked ,           /* the number of  kappa values */
  int number_of_operators,  /* number of operators */
  int nt ,                  /* time values */
  int no_q_values,          /* momentum values */
  int no_copies             /* rotations */
  )
{
  size_t nobj ; 
  FILE *fp ;
  const int32type magic_number = 44451189 ; 
  const int32type version_number  = 1  ; /** update this flag when the data format changes ***/
  int i ;
  size_t  name_len = 80 ;
  size_t howmany ;
  int check_sum ; 
  int dim ; 
  size_t corr_stride_output;

  /* Memory for some header information ***/
#define HEADER_DIM_WRITE_CORR 11
  size_t  header_size = (size_t) HEADER_DIM_WRITE_CORR ;
  int32type header_data[HEADER_DIM_WRITE_CORR] ;
#undef HEADER_DIM_WRITE_CORR

  /**** open the file ******/
  if( (fp = fopen(filename ,"wb")) == NULL )
  {
    printf("ERROR::write_2pt::Could not open the file %s\n",filename);
    exit(1);
  }

  corr_stride_output = 
    TWOPT_FORM_WHERE(nt,no_zonked-1,no_spectator-1,no_q_values-1,
			number_of_operators-1);

  dim = corr_stride_output * no_copies;
  
  /*** calculate the checksum of the form factors ****/
  check_sum =  bsd_sum((char*)(corr), sizeof(complex)*dim)  ;

  /*** pack and write the header ******/

  header_data[0]  = magic_number;
  header_data[1]  = version_number;
  header_data[2]  = nt ;
  header_data[3]  = check_sum ;
  header_data[4]  = no_q_values ;   /* no_q_values */
  header_data[5]  = no_spectator;
  header_data[6]  = no_zonked;
  header_data[7]  = number_of_operators;
  header_data[8]  = dim; 
  header_data[9]  = hl_flag;
  header_data[10] = no_copies;

  /** write the header information ***/
  if( fwrite(header_data,sizeof(int32type),header_size,fp) != header_size )
  {
    printf("There was an error during the writing of the two point (propagating) HEADER to %s\n",
	   filename);
    exit(1);
  }

  /*** write the external q--momentum to the file ***/
  howmany = (size_t) 3*(no_q_values) ;
  if( fwrite((q_momstore),sizeof(int32type),howmany,fp) != howmany   )
    {
      printf("There was an error during the writing of the Q momentum table to %s\n",
	     filename);
      exit(1);
    }

  /** write the two point functions to disk *****/
  nobj = (size_t) dim ; 

  if( fwrite((const void *)corr,sizeof(complex),nobj,fp) != nobj   )
  {
    printf("There was an error during the writing of he two point finctions to %s\n",
	   filename);
    exit(1);
  }

  /*** close the file ****/
  if( fclose(fp) != 0 )
  {
    printf("There was an error during the closing of %s \n",filename);
    exit(1);
  }

  printf("I have written the two point functions to the file %s\n",filename);

}


