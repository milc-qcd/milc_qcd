/*
 *   Write the computed matrix elements to disk.
 *
 *  MORE WORK, need to add errorbars
 *
 */


#include"read_hl_form.h"


void write_data(Real *data,Real *data_err, int dim ,char filename[80] )
{
  FILE *fp;
  int i ;

  /*** open the file ******/
 if( (fp = fopen(filename ,"w")) == NULL )
 {
   printf("ERROR:write_data: Could not open the file %s\n",filename);
   exit(1);
 }


  for(i = 0 ; i <= dim ; ++i)
    fprintf(fp,"%d   %e   %e\n",i,data[i],data_err[i]);


  /*** close the file ****/
  if( fclose(fp) != 0 )
  {
    printf("ERROR:write_data:There was an error during the closing of %s \n",filename);
    exit(1);
  }

  printf("Some data has been written to the file %s\n",filename) ; 

}

