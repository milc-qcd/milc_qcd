/* 
 *  Write the propagating three point functions
 *  to a disk file FOR ONE CHOICE OF B MESON MOMENTUM.
 *
 *  NOTE this is scalar code -- it should only be run on 
 *  one node.
 *
 */


#include "read_hl_form.h"
#include "prop_form_utilities.h"

void write_3pt_onespect(
  complex *corr,            /* the three point correlator */
  char filename[80],        /* the name of the output file */
  int spect_select,         /* Index of selected spectator */
  int32type *q_momstore,       /* Array of q momentum values */
  int32type *p_momstore,       /* Array of p momentum values */
  int hl_flag,              /* File identifier */
  int no_zonked ,           /* Zonked quarks */
  int no_spectator_corr,    /* Spectator quarks */
  int no_sequential,        /* Sequential quarks */
  int no_p_values,          /* Number of p momentum values in corr */
  int no_q_values,          /* Number of q momentum values */
  int no_oper,              /* Number of operators */
  int nt ,                  /* Number of t values */
  int no_copies             /* Number of rotations */
  )
{
  size_t nobj ; 
  FILE *fp ;
  const int magic_number = 14567332 ; 
  const int32type version_number  = 1  ; /** update this flag when the data format changes ***/
  int i ;
  size_t howmany ;
  int dim ; 
  int corr_stride_input, corr_stride_output, corr_stride_mom_op;
  int corr_stride_spect;
  int no_spectator;

  int oper_pt;         /* operator index */
  int copy_pt;         /* copy index */
  int p_pt, q_pt;      /* momentum indices */
  int where;


  /* Memory for some header information ***/
#define HEADER_DIM_WRITE_CORR 13
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

  /** Output dimensions - only one spectator value in corr array **/
  no_spectator = 1;
  corr_stride_output = 
    FORM_WHERE(nt,no_zonked-1,no_sequential-1,
				 no_spectator-1,no_q_values-1,no_p_values-1,
				 no_oper-1) ;
  dim = corr_stride_output * no_copies;

  /*** pack and write the header ******/

  header_data[0]  = magic_number;
  header_data[1]  = version_number;
  header_data[2]  = nt ;
  header_data[3]  = 0 ;
  header_data[4]  = no_p_values ;   /* no_p_values */
  header_data[5]  = no_q_values ;   /* no_q_values */
  header_data[6]  = no_oper;
  header_data[7]  = 1;              /* no_spectator */
  header_data[8]  = no_zonked;
  header_data[9]  = dim; 
  header_data[10] = hl_flag;
  header_data[11] = no_sequential;
  header_data[12] = no_copies;

  /** write the header information ***/
  if( fwrite(header_data,sizeof(int32type),header_size,fp) != header_size )
  {
    printf("There was an error during the writing of the three point (propagating) HEADER to %s\n",
	   filename);
    exit(1);
  }

  /*** write the external q momenta to the  file ***/
  howmany = (size_t) 3*(no_q_values) ;
  if( fwrite((q_momstore)  ,sizeof(int32type),howmany,fp) != howmany   )
  {
    printf("There was an error during the writing of the q-momentum table to %s\n",filename);
    exit(1);
  }

  /*** write the external p--momentum to the file ***/
  howmany = (size_t) 3 * no_p_values ;

  if( fwrite(p_momstore,sizeof(int32type),howmany,fp) != howmany   )
  {
    printf("There was an error during the writing of the P momentum table to %s\n",
	   filename);
    exit(1);
  }

  /*#define FORM_WHERE(t,zonk_pt,seq_pt,spect_pt,q_pt,p_pt,oper_pt) \
         (t) + nt*((zonk_pt) + no_zonked*((seq_pt) + no_sequential*((spect_pt) + no_spectator*((q_pt) + no_q_values*((p_pt) + no_p_values*(oper_pt))))))
  */

  /* Stride per spectator, q_value, momentum, operator, and copy in corr array */
  no_spectator = no_spectator_corr;
  corr_stride_spect = 
    FORM_WHERE(nt,no_zonked-1,no_sequential-1,0,0,0,0);
  /* Stride for all values in one copy in corr array */
  corr_stride_input = 
    FORM_WHERE(nt,no_zonked-1,no_sequential-1,
				 no_spectator-1,no_q_values-1,no_p_values-1,
				 no_oper-1) ;

  /* The spectator index is 4th in the corr array, so we
     block on the first 3 indices, select the 4th, and 
     iterate explicitly over the remaining 4 indices */

  for(copy_pt = 0; copy_pt < no_copies; copy_pt++)
    for(oper_pt = 0; oper_pt < no_oper; oper_pt++)
      for(p_pt = 0; p_pt < no_p_values; p_pt++)
	for(q_pt = 0; q_pt < no_q_values; q_pt++){

	  /** locate the required correlator in the corr array **/
	  no_spectator = no_spectator_corr;
	  where = corr_stride_input*copy_pt + 
	    FORM_WHERE(0, 0, 0, spect_select, q_pt, p_pt, oper_pt) ;
	  
	  /** write the three point functions to disk *****/
	  nobj = (size_t) corr_stride_spect ; 
	  
	  if( fwrite(&corr[where],sizeof(complex),nobj,fp) != nobj   )
	    {
	      printf("There was an error during the writing of the three point finctions to %s\n",
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

  printf("I have written the three point functions to the file %s\n",filename);

}


