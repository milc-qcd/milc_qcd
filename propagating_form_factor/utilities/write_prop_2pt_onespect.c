/* 
 *  Write the propagating two point functions
 *  to a disk file FOR ONE CHOICE OF SPECTATOR QUARK
 *
 *  NOTE this is scalar code -- it should only be run on 
 *  one node.
 *
 */


#include "read_hl_form.h"
#include "prop_form_utilities.h"

void write_prop_twopt_onespect(
  complex *corr,            /* the two point correlator */
  char filename[80],        /* the name of the output file */
  int spect_select,         /* Index of selected spectator kappa */
  int32type *q_momstore,       /* Array of momentum values */
  int hl_flag,              /* File identifier */
  int no_zonked ,           /* Zonked quarks */
  int no_spectator,         /* Spectator quarks */
  int no_q_corr_values,     /* Number of momentum values in corr */
  int no_oper,              /* Number of operators */
  int nt ,                  /* Number of t values */
  int no_copies             /* Number of rotations */
  )
{
  size_t nobj ; 
  FILE *fp ;
  int32type magic_number = 44451189 ; 
  int32type version_number  = 1  ; /** update this flag when the data format changes ***/
  int i ;
  size_t howmany ;
  int dim ; 
  int corr_stride;
  int corr_stride_spect_mom_op;
  int no_q_values;

  int mom_pt;          /* momentum index */
  int oper_pt;         /* operator index */
  int copy_pt;         /* copy index */
  int where;
  char myname[] = "write_prop_twopt_onespect";

  /* Memory for some header information ***/
#define HEADER_DIM_WRITE_CORR 11
  size_t  header_size = (size_t) HEADER_DIM_WRITE_CORR ;
  int32type header_data[HEADER_DIM_WRITE_CORR] ;
#undef HEADER_DIM_WRITE_CORR


  /**** open the file ******/
  if( (fp = fopen(filename ,"wb")) == NULL )
  {
    printf("ERROR:: %s: Could not open the file %s\n",myname,filename);
    exit(1);
  }

  /* Calculate the size of the correlator to be written */

  /** Output dimensions - only one momentum value in corr array **/
  no_q_values = no_q_corr_values;
  corr_stride = TWOPT_FORM_WHERE(nt,no_zonked-1 , 0,
				 no_q_values-1, no_oper-1) ;
  dim = corr_stride * no_copies;

  /*** pack and write the header ******/

  header_data[0]  = magic_number;
  header_data[1]  = version_number;
  header_data[2]  = nt ;
  header_data[3]  = 0 ;
  header_data[4]  = no_q_corr_values ;   /* no_q_values */
  header_data[5]  = 1;                   /* no of spectators (output) */
  header_data[6]  = no_zonked;
  header_data[7]  = no_oper;
  header_data[8]  = dim; 
  header_data[9]  = hl_flag;
  header_data[10]  = no_copies;

  /** write the header information ***/
  if( fwrite(header_data,sizeof(int32type),header_size,fp) != header_size )
  {
    printf("%s: There was an error during the writing of the two point (propagating) HEADER to %s\n",
	   myname,filename);
    exit(1);
  }

  /*** write the external q--momentum to the file ***/
  howmany = (size_t) no_q_corr_values*3 ;

  if( fwrite(q_momstore,sizeof(int32type),howmany,fp) != howmany   )
  {
    printf("%s: There was an error during the writing of the Q momentum table to %s\n",
	   myname, filename);
    exit(1);
  }

  /* Stride per spectator, momentum, operator, and copy in corr array */
  no_q_values = no_q_corr_values;
  corr_stride_spect_mom_op = 
     TWOPT_FORM_WHERE(nt,no_zonked-1 ,0,0,0);

  /* The spectator index is 3rd in the corr array, so we
     block on the first 2 indices, select the 3rd, and 
     iterate explicitly over the remaining 3 indices */

  for(copy_pt = 0; copy_pt < no_copies; copy_pt++)
    for(oper_pt = 0; oper_pt < no_oper; oper_pt++)
      for(mom_pt = 0; mom_pt < no_q_values; mom_pt++){
	
	/** locate the required correlator in the corr array **/
	where = corr_stride*copy_pt + 
	  TWOPT_FORM_WHERE(0,0,spect_select, mom_pt, oper_pt) ;
	
	/** write the two point functions to disk *****/
	nobj = (size_t) corr_stride_spect_mom_op ; 
	
	if( fwrite(&corr[where],sizeof(complex),nobj,fp) != nobj   )
	  {
	    printf("%s: There was an error during the writing of he two point finctions to %s\n",
		   myname,filename);
	    exit(1);
	  }
      }
  
  /*** close the file ****/
  if( fclose(fp) != 0 )
  {
    printf("%s: There was an error during the closing of %s \n",
	   myname,filename);
    exit(1);
  }

  printf("I have written the two point functions to the file %s\n",filename);

}


