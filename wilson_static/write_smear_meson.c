/**************************** write_smear_meson.c ****************************/
/* MIMD version 7 */
/*
 *  This routine writes the variational matrix 
 *  out to a binary disk file.
 *
 *  Some header information is dumped out at the start of
 *  the code.
 *
 *  NOTE this is scalar code -- it should only be run on 
 *  one node.
 *
 */

/** $Header: /lqcdproj/detar/cvsroot/milc_qcd/wilson_static/write_smear_meson.c,v 1.3 2007/10/09 21:02:32 detar Exp $  **/

#include "w_static_includes.h"

void write_smear_mesonx(complex *meson)
{
  FILE *fp ;
  size_t nobj = (size_t) nt*nosmear*144  ;
  /* Hack to distinguish single and double precision files */
  const int magic_number = 564132 + 8*(sizeof(Real) - 4) ;
  int i ;
  size_t  name_len = 80 ;

  /* Memory for some header information ***/
  size_t  three_oject = (size_t) 3 ;

  int32type header_data[3] ;


  header_data[0] = (int32type) magic_number ;
  header_data[1] = (int32type) nt ;
  header_data[2] = (int32type) nosmear ;

  /**** open the file ******/
  if( (fp = fopen(smear_meson_out ,"wb")) == NULL )
  {
    printf("Could not open the file %s\n",smear_meson_out );
    terminate(1);
  }

  /** write the header information ***/
  if( fwrite(header_data,sizeof(header_data[0]),three_oject,fp) != three_oject )
  {
    printf("There was an error during the writing of the smeared meson HEADER \n");
    terminate(1);
  }


  /*** write the names of the smearing functions to the file ****/
  for(i=0 ; i< nosmear ;++i)
  {
    if( fwrite(smearfile_in[i],sizeof(char),name_len,fp) != name_len )
    {
      printf("There was an error during the writing of the file name number %d \n",i);
      terminate(1); 
    }
  }

  /** write the smeraed correlators out to disk *****/
  if( fwrite(meson,sizeof(complex),nobj,fp) != nobj   )
  {
    printf("There was an error during the writing of the smeared meson correlators \n");
    terminate(1);
  }


  /*** close the file ****/
  if( fclose(fp) != 0 )
  {
    printf("There was an error during the closing of %s \n",vary_out);
    terminate(1);
  }


  printf("I have written the smeared mesons to the file %s\n",smear_meson_out);


}



