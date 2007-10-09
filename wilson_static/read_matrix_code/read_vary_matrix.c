/*   $Header: /lqcdproj/detar/cvsroot/milc_qcd/wilson_static/read_matrix_code/read_vary_matrix.c,v 1.3 2007/10/09 21:02:55 detar Exp $
 *  This is a acalar utility program to read the 
 *  static variational matrix from disk and
 *  print it to the screen --- this is useful for
 *  debugging the code.
 *  
 *  This program takes the name of the matrix
 *  file as a command line argument.
 *
 */



#include <stdio.h>
#include <math.h>
#include <string.h>
#include "../../include/int32type.h"
#include "../../include/precision.h"

void dump_vary(Real *vary_matrix,int nt, int nosmear);
void byte_rev_array(int32type buf[], int words) ;

#define NAME_LEN 80 

int main(int argc, char *argv[])
{
  FILE *fp ;
  int nt,nosmear ;
  size_t nobj   ;
  Real *vary_matrix ;
  char vary_out[80] ;
  int32type magic_number = 10213 + 8*(sizeof(Real) - 4) ;
  int32type magic_in ;
  int i ;
  char smear_name[NAME_LEN] ;
  size_t  name_len = NAME_LEN ;


  /* Memory for some header information ***/
  size_t  h_dim = (size_t) 3 ;
  int32type header_data[3] ;

  enum rev_choice { do_byte_rev = 10 , do_nothing } ; 
  int byte_rev_flag =  do_nothing ; 
  /*********..........**********........*****************/

  printf("==================================================\n");
  printf("Raw output of the static variational matrix \n");
  printf("==================================================\n");

  if( argc != 2 )
  {
    printf("usage::  %s  [matrix data file] \n",argv[0]);
    exit(1);
  }
  strcpy(vary_out,argv[1]);


  /**** open the file ******/
  if( (fp = fopen(vary_out ,"rb")) == NULL )
  {
    printf("Could not open the file %s\n",vary_out);
    exit(1);
  }

  /** read the header information ***/
  if( fread(header_data,sizeof(int32type),h_dim,fp) != h_dim )
  {
    printf("There was an error during the reading of the variational matrix HEADER \n");
    exit(1);
  }

  magic_in = header_data[0]  ;
  nt      = header_data[1]  ;
  nosmear = header_data[2]   ;

  /** check the magic number *****/
  if( magic_in == magic_number )
    printf("Magic number agreement between file and code \n");
  else
  {
    printf("ERROR:: Magic number mismatch between file and header \n");
    printf("magic code = %d magic file = %d\n",magic_number,magic_in);
    printf("I will try byte reversing"); 

    byte_rev_flag = do_byte_rev ;

    byte_rev_array( header_data , h_dim ); 

    magic_in = header_data[0]  ;
    nt      = header_data[1]  ;
    nosmear = header_data[2]   ;
    

    if( magic_in == magic_number )
      printf("Magic number agreement between file and code \n");
    else
      {
	printf("ERROR:: Magic number mismatch between file and header \n");
	printf("magic code = %d magic file = %d\n",magic_number,magic_in);
	exit(0) ;
      }


  }

  /*** read the names of the smearing functions from the file  ***/
  for(i= 0 ; i < nosmear ; ++i)
  {
    if( fread(smear_name,sizeof(char),name_len,fp) != name_len )
    {
      printf("There was an error during the reading of the file name number %d \n",i);
      exit(1);  
    }
    printf("smearing_function_file[%d] = %s\n",i,smear_name);
  }
    
  nobj    = (size_t) nt*nosmear*nosmear  ;

  printf("The number of timeslices = %d\n",nt);
  printf("The number of smearing functions = %d\n",nosmear);

  /* reserve memory for the variational matrix *****/
  vary_matrix = (Real *) calloc( (size_t) nt*nosmear*nosmear, sizeof(Real) );

  /** read the variational matrix from disk *****/
  if( fread(vary_matrix,sizeof(Real),nobj,fp) != nobj   )
  {
    printf("There was an error during the reading of the variational matrix \n");
    exit(1);
  }

  if( byte_rev_flag == do_byte_rev )
    {
      byte_rev_array((int32type*) vary_matrix, (int) nobj); 
    }

  /*** close the file ****/
  if( fclose(fp) != 0 )
  {
    printf("There was an error during the closing of %s \n",vary_out);
    exit(1);
  }


  printf("I have read the variational matrix from the file %s\n",vary_out);

  /* write the variational matrix to the screen  ***/  
  dump_vary(vary_matrix,nt, nosmear);

  /* Free up the memory used in the calculation ***/
  free(vary_matrix);

  return 0 ;

}

/*
 *     Print out the variational matrix to the 
 *     screen.
 *
 */


void dump_vary(Real *vary_matrix,int nt, int nosmear)
{
  int t ;
  int aloop,bloop ;
  int pt ;



  printf("Here is the variational matrix\n\n");
  for( t= 0 ; t < nt ;++t)
  {
    printf("Time slice = %d \n",t);
    for(aloop= 0 ; aloop < nosmear ; ++aloop)
      for(bloop= 0 ; bloop < nosmear ; ++bloop)
      {
	pt = t + nt*(aloop + nosmear*(bloop)); 

	printf("%d   %d  VM=  %e\n",aloop,bloop,vary_matrix[pt]);
      }

  }



}



/*
   Byte reverse an array of integers.
   This is a minor modification of the 
   code byterev.c, written by  Doug Toussaint.

   Usage:
       int buff[dim-1] ; byte_rev_array(buf, dim)) ;

       Real buff[dim-1] ; byte_rev_array((int*)buf, dim)) ;

     I have not tried using "double"
*/


void byte_rev_array(int32type buf[], int words)
{
  register int i ;
  register int32type old,new;

  
  for(i=0;i<words;i++)
  {
    old = buf[i];
    new = old >> 24 & 0x000000ff;
    new |= old >> 8 & 0x0000ff00;
    new |= old << 8 & 0x00ff0000;
    new |= old << 24 & 0xff000000;
    buf[i] = new;
  }


}
