/* 
 *  Write the propagating two point functions
 *  to a disk file FOR ONE CHOICE OF MOMENTUM.
 *
 *  NOTE this is scalar code -- it should only be run on 
 *  one node.
 *
 */


#include "read_hl_form.h"
#include "prop_form_utilities.h"

void write_prop_twopt_onemom(
  complex *corr,            /* the two point correlator */
  char filename[80],        /* the name of the output file */
  int q_mom_select,         /* Index of selected momentum */
  int32type *q_momstore,       /* Array of momentum values */
  int hl_flag,              /* File identifier */
  int no_zonked ,     /* Zonked quarks */
  int no_spectator,         /* Spectator quarks */
  int no_q_corr_values,     /* Number of momentum values in corr */
  int no_oper,              /* Number of operators */
  int nt ,                  /* Number of t values */
  int no_copies             /* Number of rotations */
  )
{
  size_t nobj ; 
  FILE *fp ;
  const int32type magic_number = 44451189 ; 
  const int32type version_number  = 1  ; /** update this flag when the data format changes ***/
  int i ;
  size_t howmany ;
  int dim ; 
  int corr_stride;
  int corr_stride_mom_op;
  int no_q_values;

  int oper_pt;         /* operator index */
  int copy_pt;         /* copy index */
  int where;


  /* Memory for some header information ***/
#define HEADER_DIM_WRITE_CORR 11
  size_t  header_size = (size_t) HEADER_DIM_WRITE_CORR ;
  int32type header_data[HEADER_DIM_WRITE_CORR] ;
#undef HEADER_DIM_WRITE_CORR


  /**** open the file ******/
  if( (fp = fopen(filename ,"wb")) == NULL )
  {
    printf("ERROR::write_prop_form::Could not open the file %s\n",filename);
    exit(1);
  }

  /* Calculate the size of the correlator to be written */

  /** Output dimensions - only one momentum value in corr array **/
  no_q_values = 1;
  corr_stride = TWOPT_FORM_WHERE(nt,no_zonked-1 ,no_spectator-1,
				 0, no_oper-1) ;
  dim = corr_stride * no_copies;

  /*** pack and write the header ******/

  header_data[0]  = magic_number;
  header_data[1]  = version_number;
  header_data[2]  = nt ;
  header_data[3]  = 0 ;
  header_data[4]  = 1 ;   /* no_q_values */
  header_data[5]  = no_spectator;
  header_data[6]  = no_zonked;
  header_data[7]  = no_oper;
  header_data[8]  = dim; 
  header_data[9]  = hl_flag;
  header_data[10]  = no_copies;

  /** write the header information ***/
  if( fwrite(header_data,sizeof(int32type),header_size,fp) != header_size )
  {
    printf("There was an error during the writing of the two point (propagating) HEADER to %s\n",
	   filename);
    exit(1);
  }

  /*** write the external q--momentum to the file ***/
  howmany = (size_t) 3 ;

  if( fwrite(&q_momstore[3*q_mom_select],sizeof(int32type),howmany,fp) != howmany   )
  {
    printf("There was an error during the writing of the Q momentum table to %s\n",
	   filename);
    exit(1);
  }

  /* Stride per momentum, operator, and copy in corr array */
  no_q_values = no_q_corr_values;
  corr_stride_mom_op = 
     TWOPT_FORM_WHERE(nt,no_zonked-1 ,no_spectator-1,0,0);

  /* The momentum index is 4th in the corr array, so we
     block on the first 3 indices, select the 4th, and 
     iterate explicitly over the remaining 2 indices */

  for(copy_pt = 0; copy_pt < no_copies; copy_pt++)
    for(oper_pt = 0; oper_pt < no_oper; oper_pt++){

      /** locate the required correlator in the corr array **/
      where = corr_stride*copy_pt + 
	TWOPT_FORM_WHERE(0,0,0, q_mom_select, oper_pt) ;
      
      /** write the two point functions to disk *****/
      nobj = (size_t) corr_stride_mom_op ; 
      
      if( fwrite(&corr[where],sizeof(complex),nobj,fp) != nobj   )
	{
	  printf("There was an error during the writing of he two point finctions to %s\n",
		 filename);
	  exit(1);
	}
    }

  /*** close the file ****/
  if( fclose(fp) != 0 )
  {
    printf("There was an error during the closing of %s \n",filename);
    exit(1);
  }

  printf("I have written the two point functions to the file %s\n",filename);

}


