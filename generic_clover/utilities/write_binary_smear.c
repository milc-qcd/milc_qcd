/*
 * Write s spatial smearing function in a format that
 * can be loaded into the static variational code.
 *
 */



#include <stdio.h>
#include <math.h>

/*
 *  Ths routine writes a  single smearing function,
 *  over the spatial lattice, to a disk file.
 *
 *  This is scalar code
 *
 *  The smearing function is written to a binary file.
 *  Some header information is written to the start of the file.
 *  
 *
 *  Function arguments
 *  ------------------
 *   On input 
 *     dim          :: the total amount of data
 *     nx           :: linear dimension of the lattice
 *     data[1..dim] :: the data to be written to disk
 *     filename     :: the name of the file to write the data too
 *
 */

void write_smear_wave(Real *data, int nx, int dim, char filename[])
{
  FILE *fp ;
  size_t nobj = (size_t) dim  ;

  const int magic_number = 43241 ;
  char filenamemod[80] ;

  /* Memory for some header information ***/
  size_t  three_object = (size_t) 3 ;
  int header_data[3] ;

  header_data[0] = magic_number ;
  header_data[1] = dim ;
  header_data[2] = nx  ;

  /** append lattice dimension to filename */
  sprintf(filenamemod,"%s_ns%d",filename,nx);

  /** open the file ****/
  if( (fp = fopen(filenamemod,"wb")) == NULL )
  {
    printf("Could not open the file %s\n",filenamemod);
    exit(1);
  }
  
  /** write the header information to disk *******/
  if( fwrite(header_data,sizeof(int),three_object,fp) != three_object )
  {
    printf("There was an error during the writing of the smearing function HEADER \n");
    exit(1);
  }


  /*** write the smearing function  ***/
  if( fwrite(data,sizeof(Real),nobj,fp) != nobj   )
  {
    printf("There was an error during the writing of the smearing function \n");
    exit(1);
  }

  /*** close the file ***/
  if( fclose(fp) != 0 )
  {
    printf("There was an error during the closing of %s \n",filenamemod);
    exit(1);
  }

  printf("I have written a smearing function to the file %s\n",filenamemod);



}


