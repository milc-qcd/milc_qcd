/* 
 *  Read the propagating two point functions
 *  from a disk file.
 *
 *  NOTE this is scalar code -- it should only be run on 
 *  one node.
 *
 *  Subroutine argument
 *   filename  :: 
 *   corr      :: 
 *   no_k_one :: 
 *   no_k_two :: 
 *   hl_flag   :: flag saying where the heavy-heavy or  heavy-light correlators
 *                are to be written.
 *   number_of_operators  :: 
 *   dim       :: the amount of data

 */


#include "read_hl_form.h"
#include "prop_form_utilities.h"



/* Functions declarations used in the routines in this file ***/

void read_prop_twopt(
  complex **corr,           /* the two point correlators */
  int **corr_oper_list,     /* List of operator numbers read */
  int **corr_copy_list,     /* List of copy numbers read */
  char filename[80],        /* the name of the disk file */
  int32type **q_momstore,      /* list of momenta */
  int *hl_flag,             /* flag identifying type of correlator */
  int *no_k_one ,           /* the number of  kappa values */
  int *no_k_two,            /* the number of  kappa values */
  int *number_of_operators,  /* number of operators */
  int *nt ,                 /* time values */
  int *no_q_values,         /* momentum values */
  int *nocopies             /* rotations */
  )
{
  size_t nobj ; 
  FILE *fp ;
  int32type magic_number = 44451189 ; 
  int32type version_number  = 1  ; /** update this flag when the data format changes ***/
  int i ;
  size_t  name_len = 80 ;
  size_t howmany ;
  int check_sum ; 
  int check_sum_in ; 
  int dim ; 
  int byte_rev_flag  =  do_nothing ;

  /* Memory for some header information ***/
#define HEADER_DIM_WRITE_CORR 11
  size_t  header_size = (size_t) HEADER_DIM_WRITE_CORR ;
  int32type header_data[HEADER_DIM_WRITE_CORR] ;
#undef HEADER_DIM_WRITE_CORR


  /**** open the file ******/
  if( (fp = fopen(filename ,"rb")) == NULL )
  {
    printf("ERROR::read_prop_form::Could not open the file %s\n",filename);
    exit(1);
  }

  /** read the header information ***/
  if( fread(header_data,sizeof(int32type),header_size,fp) != header_size )
  {
    printf("There was an error during the reading of the two point (propagating) HEADER to %s\n",
	   filename);
    exit(1);
  }

  /*** unpack and check the header ******/

  /***** check the header ****/
  if( header_data[0] != magic_number )
  {
    printf("ERROR: magic number mismatch between code %d and file %d\n",magic_number, header_data[0]  );

    printf("I will try byte reversing the data\n"); 
    byte_rev_flag = do_byte_rev ;

  }


  if(  byte_rev_flag == do_byte_rev    )
  {
    byte_rev_array( header_data, header_size); 
  }


  /***** check the header ****/
  if( header_data[0] != magic_number )
  {
    printf("ERROR: magic number mismatch between code %d and file %d\n",magic_number, header_data[0]  );
    exit(1); 
  }


  if( header_data[1] != version_number )
  {
    printf("ERROR: version number mismatch between code %d and file %d\n",version_number, header_data[1]);
    exit(1); 
  }





  *nt  = header_data[2] ;
  check_sum_in  = header_data[3]  ;
  *no_q_values = header_data[4]  ;
  *no_k_one = header_data[5]  ;
  *no_k_two = header_data[6]  ;
  *number_of_operators = header_data[7]  ;
  dim = header_data[8]  ; 
  *hl_flag = header_data[9]  ;
  *nocopies = header_data[10]  ;

  /*** read the external q--momentum from the file ***/
  howmany = (size_t) 3*(*no_q_values) ;
  if( ( (*q_momstore) = (int32type *)calloc( howmany , sizeof(int32type))   )  == NULL) 
  {
    printf("There was an error in allocating q-momentum \n");
    exit(1);
  }


  if( fread((*q_momstore),sizeof(int32type),howmany,fp) != howmany   )
  {
    printf("There was an error during the reading of the Q momentum table to %s\n",
	   filename);
    exit(1);
  }

  if(  byte_rev_flag == do_byte_rev    )
  {
    byte_rev_array((*q_momstore) , (int) howmany );
  }


  /** read the two point functions from disk *****/
  nobj = (size_t) dim ; 
  if( ((*corr) = (complex *)calloc( nobj , sizeof(complex) ) )  == NULL) 
  {
    printf("There was an error in allocating the two point function, nobj = %d\n",nobj);
    exit(1);
  }


  if( fread((*corr),sizeof(complex),nobj,fp) != nobj   )
  {
    printf("There was an error during the reading of he two point finctions from %s\n",
	   filename);
    exit(1);
  }


  if(  byte_rev_flag == do_byte_rev    )
  {
    byte_rev_array((int32type*) (*corr), dim*2 );
  }



  /*** close the file ****/
  if( fclose(fp) != 0 )
  {
    printf("There was an error during the closing of %s \n",filename);
    exit(1);
  }

 


  /*** calculate the checksum of the form factors ****/
  check_sum =  bsd_sum((char*)(*corr), sizeof(complex)*dim)  ;

  if( check_sum !=  check_sum_in )
  {
    printf("ERROR: checksum mismatch bewteen the file %d and the code %d \n",check_sum,check_sum_in) ; 
    printf("ERROR: checksum mismatch bewteen the file %d and the code %d \n",check_sum,check_sum_in) ; 
    printf("ERROR: checksum mismatch bewteen the file %d and the code %d \n",check_sum,check_sum_in) ; 

/**    exit(1) ;  ****/
  }

  /* Create operator list */
  if( ((*corr_oper_list) = (int *)calloc( *number_of_operators, sizeof(int) ) )  == NULL) 
    {
      printf("There was an error in allocating the operator list\n");
      exit(1);
    }

  for(i = 0; i < *number_of_operators; i++)
    (*corr_oper_list)[i] = i;

  /* Create copies list */
  if( ((*corr_copy_list) = (int *)calloc( *nocopies, sizeof(int) ) )  == NULL) 
    {
      printf("There was an error in allocating the copies list\n");
      exit(1);
    }

  for(i = 0; i < *nocopies; i++)
    (*corr_copy_list)[i] = i;

  printf("I have read the two point functions from the file %s\n",filename);

}


