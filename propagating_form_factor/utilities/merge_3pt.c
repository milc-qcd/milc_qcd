/* 
 *  Read and MERGE the propagating heavy---> light three point functions
 *  from disk files.
 *
 *  NOTE this is scalar code -- it should only be run on 
 *  one node.
 *
 */


#include "read_hl_form.h"
#include "prop_form_utilities.h"

void read_merge1_corr(
  complex **corr,           /* the two point correlators */
  int **corr_oper_list,     /* List of operator numbers read */
  int **corr_copy_list,     /* List of copy numbers read */
  int32type **p_momstore,      /* list of B meson momenta */
  int32type **q_momstore,      /* list of momentum transfers */
  int *pt_nt,               /* time values */
  int *pt_no_p_values,      /* momentum values */
  int *pt_no_q_values,      /* momentum values */
  int *pt_no_oper,          /* number of operators */
  int *pt_no_spectator,     /* the number of  kappa values */
  int *pt_no_sequential,    /* the number of  kappa values */
  int *pt_no_zonked ,       /* the number of  kappa values */
  int *hl_flag,             /* flag identifying type of correlator */
  char filename[80],        /* the name of the disk file */
  int *pt_no_copies         /* rotations */
)
{
  FILE *fp ;
  size_t nobj, skip ;
  int nobj_skip ; 
  const int32type magic_number = 14567332 ; 
  const int32type version_number  = 1  ; /** update this flag when the data format changes ***/
  int i ;
  size_t  name_len = 80 ;
  size_t howmany ;
  int check_sum ; 
  int check_sum_in ; 
  int dim ; 
  int byte_rev_flag  =  do_nothing ;
  int nt, no_p_values, no_q_values, no_oper, no_spectator, no_sequential,
    no_zonked, no_copies;
  int spect_pt, q_pt, p_pt, oper_pt, copy_pt;
  int corr_stride;
  int where;

  /* Memory for some header information ***/
#define HEADER_DIM_WRITE_CORR 13
  size_t  header_size = (size_t) HEADER_DIM_WRITE_CORR ;
  int32type header_data[HEADER_DIM_WRITE_CORR] ;


  /**** open the file ******/
  if( (fp = fopen(filename ,"rb")) == NULL )
  {
    printf("Could not open the file %s\n",filename);
    exit(1);
  }

  /** read the header information ***/
  if( fread(header_data,sizeof(int32type),header_size,fp) != header_size )
  {
    printf("There was an error during the raeding of the form factor HEADER from %s\n",filename);
    exit(1);
  }




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


  *pt_nt = header_data[2]  ;
  check_sum_in = header_data[3] ;
  *pt_no_p_values = header_data[4] ;
  *pt_no_q_values = header_data[5] ;
  *pt_no_oper = header_data[6] ;
  *pt_no_spectator = header_data[7]   ;
  *pt_no_zonked = header_data[8]  ;
  dim = header_data[9]  ;
  *hl_flag  = header_data[10]  ;
  *pt_no_sequential = header_data[11]  ;
  *pt_no_copies = header_data[12]  ;

  nt = *pt_nt;
  no_p_values = *pt_no_p_values;
  no_q_values = *pt_no_q_values;
  no_oper = *pt_no_oper;
  no_spectator = *pt_no_spectator;
  no_zonked = *pt_no_zonked;
  no_sequential = *pt_no_sequential;
  no_copies = *pt_no_copies;

  /*******  read some more data from the file **********/

  /*** read the external momentum from the  file ***/
  howmany = (size_t) 3*(no_q_values) ;
  if( ( (*q_momstore) = (int32type *)calloc( howmany , sizeof(int32type))   )  == NULL) 
  {
    printf("There was an error in allocating q-momentum \n");
    exit(1);
  }

  if( fread((*q_momstore)  ,sizeof(int32type),howmany,fp) != howmany   )
  {
    printf("There was an error during the reading of the q-momentum table from %s\n",filename);
    exit(1);
  }

  if(  byte_rev_flag == do_byte_rev    )
  {
    byte_rev_array((*q_momstore) , (int) howmany );
  }




  /*** read the external momentum from the  file ***/
  howmany = (size_t) 3*(no_p_values) ;
  if( ( (*p_momstore) = (int32type *)calloc( howmany , sizeof(int32type))   )  == NULL) 
  {
    printf("There was an error in allocating p--momentum \n");
    exit(1);
  }

  if( fread((*p_momstore)  ,sizeof(int32type),howmany,fp) != howmany   )
  {
    printf("There was an error during the reading of the p-momentum table from %s\n",filename);
    exit(1);
  }

  if(  byte_rev_flag == do_byte_rev    )
  {
    byte_rev_array((*p_momstore) , (int) howmany );
  }




  /*** reserve memory for the three point function ****/
  nobj = (size_t) dim ; 
  if( ((*corr) = (complex *)calloc( nobj , sizeof(complex) ) )  == NULL) 
  {
    printf("There was an error in allocating the three point function, nobj = %d\n",nobj);
    exit(1);
  }

  /** read the three point function to disk *****/
  /* Here we read data for only the FIRST SEQUENTIAL QUARK */

  /* stride between copies */
  corr_stride = FORM_WHERE(nt,no_zonked-1,no_sequential-1,
			   no_spectator-1,no_q_values-1,no_p_values-1,
			   no_oper-1) ; 

  /* No initial skip required for the first sequential quark */

  /* Number of objects to read each time */
  nobj = nt*no_zonked;

  /* Number of objects to skip each time */
  
  nobj_skip = nt*no_zonked*(no_sequential-1);
  skip = nobj_skip*sizeof(complex);
  
  for(copy_pt = 0 ; copy_pt  < no_copies   ; copy_pt++)
    for(oper_pt = 0 ; oper_pt  < no_oper     ; oper_pt++)
      for(p_pt = 0    ; p_pt     < no_p_values ; p_pt++)
	for(q_pt = 0    ; q_pt     < no_q_values ; q_pt++)
	  for(spect_pt = 0; spect_pt < no_spectator; spect_pt++){

	    where = copy_pt*corr_stride + 
	      FORM_WHERE(0,0,0,spect_pt,q_pt,p_pt,oper_pt) ;
	    
	    if( fread((*corr+where),sizeof(complex),nobj,fp) != nobj   )
	      {
		printf("There was an error during the reading of the form factor data from %s\n",filename);
		exit(1);
	      }
	    
	    if(  byte_rev_flag == do_byte_rev    )
	      {
		byte_rev_array((int32type*)(*corr+where), nobj*2 );
	      }
	    /* Position the file for next read */
	    if(fseek(fp,skip,SEEK_CUR) != 0)
	      {
		printf("There was an error seeking the two point function on %s\n",
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

  /* Create operator list */
  if( ((*corr_oper_list) = (int *)calloc( no_oper, sizeof(int) ) )  == NULL) 
    {
      printf("There was an error in allocating the operator list\n");
      exit(1);
    }

  for(i = 0; i < no_oper; i++)
    (*corr_oper_list)[i] = i;

  /* Create copy list */
  if( ((*corr_copy_list) = (int *)calloc( no_copies, sizeof(int) ) )  == NULL) 
    {
      printf("There was an error in allocating the copy list\n");
      exit(1);
    }

  for(i = 0; i < no_copies; i++)
    (*corr_copy_list)[i] = i;

  printf("I have read propagating three point functions from the file %s\n",filename);

}


void read_merge2_corr(
  complex **corr,           /* the two point correlators */
  int **corr_oper_list,     /* List of operator numbers read */
  int **corr_copy_list,     /* List of copy numbers read */
  int32type **p_momstore,      /* list of B meson momenta */
  int32type **q_momstore,      /* list of momentum transfers */
  int *pt_nt,               /* time values */
  int *pt_no_p_values,      /* momentum values */
  int *pt_no_q_values,      /* momentum values */
  int *pt_no_oper,          /* number of operators */
  int *pt_no_spectator,     /* the number of  kappa values */
  int *pt_no_sequential,    /* the number of  kappa values */
  int *pt_no_zonked ,       /* the number of  kappa values */
  int *hl_flag,             /* flag identifying type of correlator */
  char filename[80],        /* the name of the disk file */
  int *pt_no_copies         /* rotations */
)
{
  FILE *fp ;
  size_t nobj, skip ;
  int nobj_skip ; 
  const int32type magic_number = 14567332 ; 
  const int32type version_number  = 1  ; /** update this flag when the data format changes ***/
  int i ;
  size_t  name_len = 80 ;
  size_t howmany ;
  int check_sum ; 
  int check_sum_in ; 
  int dim ; 
  int nt = *pt_nt;
  int no_p_values = *pt_no_p_values;
  int no_q_values = *pt_no_q_values;
  int no_oper = *pt_no_oper;
  int no_spectator = *pt_no_spectator;
  int no_sequential = *pt_no_sequential;
  int no_zonked = *pt_no_zonked;
  int no_copies = *pt_no_copies;
  int byte_rev_flag  =  do_nothing ;
  int32type *q_momstore2, *p_momstore2;
  int spect_pt, q_pt, p_pt, oper_pt, copy_pt;
  int corr_stride;
  int where;
  
  /* Memory for some header information ***/
#define HEADER_DIM_WRITE_CORR 13
  size_t  header_size = (size_t) HEADER_DIM_WRITE_CORR ;
  int32type header_data[HEADER_DIM_WRITE_CORR] ;


  /**** open the file ******/
  if( (fp = fopen(filename ,"rb")) == NULL )
  {
    printf("Could not open the file %s\n",filename);
    exit(1);
  }

  /** read the header information ***/
  if( fread(header_data,sizeof(int32type),header_size,fp) != header_size )
  {
    printf("There was an error during the raeding of the form factor HEADER from %s\n",filename);
    exit(1);
  }




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


  if( header_data[2] != nt )
  {
    printf("ERROR: nt mismatch between code %d and file %d\n",nt, header_data[2]);
    exit(1); 
  }

  if( header_data[4] != no_p_values )
  {
    printf("ERROR: no_p_values mismatch between code %d and file %d\n",no_p_values, header_data[4]);
    exit(1); 
  }

  if( header_data[5] != no_q_values )
  {
    printf("ERROR: no_q_values mismatch between code %d and file %d\n",no_q_values, header_data[5]);
    exit(1); 
  }

  if( header_data[6] != no_oper )
  {
    printf("ERROR: no_oper mismatch between code %d and file %d\n",no_oper, header_data[6]);
    exit(1); 
  }

  if( header_data[7] != no_spectator )
  {
    printf("ERROR: no_spectator mismatch between code %d and file %d\n",no_spectator, header_data[7]);
    exit(1); 
  }

  if( header_data[8] != no_zonked )
  {
    printf("ERROR: no_zonked mismatch between code %d and file %d\n",no_zonked, header_data[8]);
    exit(1); 
  }

  if( header_data[10] != *hl_flag )
  {
    printf("ERROR: hl_flag mismatch between code %d and file %d\n",*hl_flag, header_data[10]);
    exit(1); 
  }

  /* Note: we retain no_sequential = larger merged value from first file */

  if( header_data[11] != no_sequential-1 )
  {
    printf("ERROR: no_sequential mismatch between code %d and file %d\n",no_sequential-1, header_data[11]);
    exit(1); 
  }

  if( header_data[12] != no_copies )
  {
    printf("ERROR: no_copies mismatch between code %d and file %d\n",no_copies, header_data[12]);
    exit(1); 
  }

  /*******  read some more data from the file **********/

  /*** read the external momentum from the  file ***/
  howmany = (size_t) 3*(no_q_values) ;
  if( ( (q_momstore2) = (int32type *)calloc( howmany , sizeof(int32type))   )  == NULL) 
  {
    printf("There was an error in allocating q-momentum \n");
    exit(1);
  }

  if( fread((q_momstore2)  ,sizeof(int32type),howmany,fp) != howmany   )
  {
    printf("There was an error during the reading of the q-momentum table from %s\n",filename);
    exit(1);
  }

  if(  byte_rev_flag == do_byte_rev    )
  {
    byte_rev_array((q_momstore2) , (int) howmany );
  }


  for(i = 0; i < no_q_values*3; i++)
    {
      if(q_momstore2[i] != (*q_momstore)[i])
	{
	  printf("ERROR: q_mom mismatch between code %d and file %d\n",(*q_momstore)[i], q_momstore2[i]);
	  exit(1); 
	}
    }

  /*** read the external momentum from the  file ***/
  howmany = (size_t) 3*(no_p_values) ;
  if( ( (p_momstore2) = (int32type *)calloc( howmany , sizeof(int32type))   )  == NULL) 
  {
    printf("There was an error in allocating p--momentum \n");
    exit(1);
  }

  if( fread((p_momstore2)  ,sizeof(int32type),howmany,fp) != howmany   )
  {
    printf("There was an error during the reading of the p-momentum table from %s\n",filename);
    exit(1);
  }

  if(  byte_rev_flag == do_byte_rev    )
  {
    byte_rev_array((p_momstore2) , (int) howmany );
  }

  for(i = 0; i < no_p_values*3; i++)
    {
      if(p_momstore2[i] != (*p_momstore)[i])
	{
	  printf("ERROR: p_mom mismatch between code %d and file %d\n",(*p_momstore)[i], p_momstore2[i]);
	  exit(1); 
	}
    }

  /** read the three point function to disk *****/
  /* Here we read data for all sequential quarks in this file, which
     is the number for the first file minus one.  But we begin storing
     in the place for the SECOND sequential quark */

  /* Stride for corr array between "copies" */
  corr_stride = FORM_WHERE(nt,no_zonked-1,no_sequential-1,
			   no_spectator-1,no_q_values-1,no_p_values-1,
			   no_oper-1) ; 

  /* Number of objects to read each time */
  nobj = nt*no_zonked*(no_sequential-1);

  for(copy_pt = 0 ; copy_pt  < no_copies   ; copy_pt++)
    for(oper_pt = 0 ; oper_pt  < no_oper     ; oper_pt++)
      for(p_pt = 0    ; p_pt     < no_p_values ; p_pt++)
	for(q_pt = 0    ; q_pt     < no_q_values ; q_pt++)
	  for(spect_pt = 0; spect_pt < no_spectator; spect_pt++){

	    /* We start storing at seq_pt = 1 */
	    where = copy_pt*corr_stride + 
	      FORM_WHERE(0,0,1,spect_pt,q_pt,p_pt,oper_pt) ; 
	    
	    if( fread((*corr+where),sizeof(complex),nobj,fp) != nobj   )
	      {
		printf("There was an error during the reading of the form factor data from %s\n",filename);
		exit(1);
	      }
	    
	    if(  byte_rev_flag == do_byte_rev    )
	      {
		byte_rev_array((int32type*)(*corr+where), nobj*2 );
	      }
	  }

  /*** close the file ****/
  if( fclose(fp) != 0 )
  {
    printf("There was an error during the closing of %s \n",filename);
    exit(1);
  }

  printf("I have read propagating three point functions from the file %s\n",filename);

}


