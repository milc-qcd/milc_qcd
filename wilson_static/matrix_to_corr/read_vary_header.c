
/*
 *  Read the header from a static variational matrix.
 *  This version does not do any clever byte reversing stuff.
 *
 */

#include LATDEF

FILE *(read_vary_header(char filename[], int local_src, int which_source, char src_name[],
int *nosmear, int *nt,  int *byte_rev_flag  ))
{
  FILE *fp ;
  const int magic_number = 10213 ;
  int magic_in ;

  /* Memory for some header information ***/
#define HEADER_DIM 3 
  size_t h_dim = (size_t) HEADER_DIM ;
  int header_data[ HEADER_DIM ] ;
#undef HEADER_DIM

  int i ;
  char smear_name[NAME_LEN] ;
  size_t  name_len = NAME_LEN ;

  /**------------------------------------------------------------------**/

  *byte_rev_flag  =  do_nothing ;


  /**** open the file ******/
  if( (fp = fopen(filename ,"rb")) == NULL )
  {
    printf("Could not open the file %s\n",filename);
    exit(1);
  }

  /** read the header information ***/
  if( fread(header_data,sizeof(int),h_dim,fp) !=h_dim )
  {
    printf("There was an error during the reading of the variational matrix HEADER \n");
    exit(1);
  }


  /***** check the header ****/
  if( header_data[0] != magic_number )
  {
    printf("Magic number mismatch between code %d and file %d\n",magic_number, header_data[0]  );
    printf("I will try byte reversing the data\n"); 

    *byte_rev_flag = do_byte_rev ;

  }

  if(  *byte_rev_flag == do_byte_rev    )
  {
    byte_rev_array( header_data , h_dim ); 
  }



  magic_in = header_data[0] ;
  *nt      = header_data[1] ;
  *nosmear = header_data[2]  ;

  /** check the magic number *****/
  if( magic_in == magic_number )
    printf("Magic number agreement between file and code \n");
  else
  {
    printf("ERROR:: Magic number mismatch between file and header \n");
    printf("magic code = %d magic file = %d\n",magic_number,magic_in);
    exit(2);
  }



  /*** read the names of the smearing functions from the file  ***/
  for(i= 0 ; i < *nosmear ; ++i)
  {
    if( fread(smear_name,sizeof(char),name_len,fp) != name_len )
    {
      printf("There was an error during the writing of the file name number %d \n",i);
      exit(1);  
    }


    if( i == which_source)
      strcpy(src_name ,smear_name);

    if( i == local_src)
      printf("local_source  = %s\n",smear_name);

  }



  return fp ;
}



