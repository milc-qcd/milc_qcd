/* 
 *  Read in the propagating heavy---> light three point functions
 *  from a disk file.
 *
 *  NOTE this is scalar code -- it should only be run on 
 *  one node.
 *
 */


#include "read_hl_form.h"
#include "prop_form_utilities.h"
#include <errno.h>

/* Functions declarations used in the routines in this file ***/

void read_propagating_form_corr(
  complex **corr,           /* the two point correlators */
  int **corr_oper_list,     /* List of operator numbers read */
  int **corr_copy_list,     /* List of copy numbers read */
  int32type **p_momstore,      /* list of B meson momenta */
  int32type **q_momstore,      /* list of momentum transfers */
  int *nt,                  /* time values */
  int *no_p_values,         /* momentum values */
  int *no_q_values,         /* momentum values */
  int *no_oper,             /* number of operators */
  int *no_spectator,        /* the number of  kappa values */
  int *no_sequential,       /* the number of  kappa values */
  int *no_zonked ,          /* the number of  kappa values */
  int *hl_flag,             /* flag identifying type of correlator */
  char filename[80],        /* the name of the disk file */
  int *nocopies             /* rotations */
)
{
  char myname[] = "read_propagating_form_corr";
  FILE *fp ;
  size_t nobj ; 
  const int32type magic_number = 14567332 ; 
  const int32type version_number  = 1  ; /** update this flag when the data format changes ***/
  int i ;
  size_t  name_len = 80 ;
  size_t howmany ;
  int check_sum ; 
  int check_sum_in ; 
  int dim ; 
  int byte_rev_flag  =  do_nothing ;
  
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
    printf("%s: error %d reading the form factor HEADER from %s\n",
	   myname,errno,filename);
    exit(1);
  }




  /***** check the header ****/
  if( header_data[0] != magic_number )
  {
    printf("%s: ERROR: magic number mismatch between code %d and file %d\n",
	   myname,magic_number, header_data[0]  );
    printf("I will try byte reversing the data\n"); 
    byte_rev_flag = do_byte_rev ;
  }

  if(  byte_rev_flag == do_byte_rev    )
  {
    byte_rev_array( header_data, header_size); 
  }



  if( header_data[0] != magic_number )
  {
    printf("%s: ERROR: magic number mismatch between code %d and file %d\n",
	   myname,magic_number, header_data[0]  );
    exit(1); 
  }

  if( header_data[1] != version_number )
  {
    printf("%s: ERROR: version number mismatch between code %d and file %d\n",
	   myname,version_number, header_data[1]);
    exit(1); 
  }


  *nt = header_data[2]  ;
  check_sum_in = header_data[3] ;
  *no_p_values = header_data[4] ;
  *no_q_values = header_data[5] ;
  *no_oper = header_data[6] ;
  *no_spectator = header_data[7]   ;
  *no_zonked = header_data[8]  ;
  dim = header_data[9]  ;
  *hl_flag  = header_data[10]  ;
  *no_sequential = header_data[11]  ;
  *nocopies = header_data[12]  ;

  /*******  read some more data from the file **********/

  /*** read the external momentum from the  file ***/
  howmany = (size_t) 3*(*no_q_values) ;
  if( ( (*q_momstore) = (int32type *)calloc( howmany , sizeof(int32type))   )  == NULL) 
  {
    printf("%s: There was an error in allocating q-momentum \n",myname);
    exit(1);
  }

  if( fread((*q_momstore)  ,sizeof(int32type),howmany,fp) != howmany   )
  {
    printf("%s: error %d reading the q-momentum table from %s\n",
	   myname,errno,filename);
    exit(1);
  }

  if(  byte_rev_flag == do_byte_rev    )
  {
    byte_rev_array((*q_momstore) , (int) howmany );
  }




  /*** read the external momentum from the  file ***/
  howmany = (size_t) 3*(*no_p_values) ;
  if( ( (*p_momstore) = (int32type *)calloc( howmany , sizeof(int32type))   )  == NULL) 
  {
    printf("%s: There was an error in allocating p--momentum \n",myname);
    exit(1);
  }

  if( fread((*p_momstore)  ,sizeof(int32type),howmany,fp) != howmany   )
  {
    printf("%s: error %d reading the p-momentum table from %s\n",
	   myname,errno,filename);
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
    printf("%s: There was an error in allocating the three point function, nobj = %d\n",
	   myname,nobj);
    exit(1);
  }

  /** read the three point function to disk *****/
  if( fread((*corr),sizeof(complex),nobj,fp) != nobj   )
  {
    printf("%s: error %d reading the form factor data from %s\n",
	   myname,errno,filename);
    exit(1);
  }

  if(  byte_rev_flag == do_byte_rev    )
  {
    byte_rev_array((int32type*) (*corr), dim*2 );
  }


  /*** calculate the checksum of the form factors ****/
  check_sum =  bsd_sum((char*)(*corr), sizeof(complex)*dim)  ;

  if( check_sum !=  check_sum_in )
  {
    printf("%s: ERROR: checksum mismatch bewteen the file %d and the code %d \n",
	   myname,check_sum,check_sum_in) ; 
    printf("%s: ERROR: checksum mismatch bewteen the file %d and the code %d \n",
	   myname,check_sum,check_sum_in) ; 
    printf("%s: ERROR: checksum mismatch bewteen the file %d and the code %d \n",
	   myname,check_sum,check_sum_in) ; 

/**    exit(1) ;  ****/
  }


  /*** close the file ****/
  if( fclose(fp) != 0 )
  {
    printf("%s: error %d closing %s \n",myname,errno,filename);
    exit(1);
  }

  /* Create operator list */
  if( ((*corr_oper_list) = (int *)calloc( *no_oper, sizeof(int) ) )  == NULL) 
    {
      printf("%s: There was an error in allocating the operator list\n",
	     myname);
      exit(1);
    }

  for(i = 0; i < *no_oper; i++)
    (*corr_oper_list)[i] = i;

  /* Create copy list */
  if( ((*corr_copy_list) = (int *)calloc( *nocopies, sizeof(int) ) )  == NULL) 
    {
      printf("%s: There was an error in allocating the copy list\n",
	     myname);
      exit(1);
    }

  for(i = 0; i < *nocopies; i++)
    (*corr_copy_list)[i] = i;

  printf("I have read propagating three point functions from the file %s\n",filename);

}

/* Same as above, but reads data for only one B meson momentum */
void read_3pt_onemom(
  complex **corr,           /* the two point correlators */
  int **corr_oper_list,     /* List of operator numbers read */
  int **corr_copy_list,     /* List of copy numbers read */
  int32type **p_momstore,      /* list of B meson momenta */
  int32type **q_momstore,      /* list of momentum transfers */
  int *nt,                  /* time values */
  int *no_p_values,         /* momentum values */
  int *no_q_values,         /* momentum values */
  int *no_oper,             /* number of operators */
  int *no_spectator,        /* the number of  kappa values */
  int *no_sequential,       /* the number of  kappa values */
  int *no_zonked ,          /* the number of  kappa values */
  int *hl_flag,             /* flag identifying type of correlator */
  char filename[80],        /* the name of the disk file */
  int *nocopies,            /* rotations */
  int p_mom_select          /* Momentum index selected */
)
{
  FILE *fp ;
  size_t nobj ; 
  const int32type magic_number = 14567332 ; 
  const int32type version_number  = 1  ; /** update this flag when the data format changes ***/
  int i ;
  size_t  name_len = 80 ;
  size_t howmany ;
  int check_sum ; 
  int check_sum_in ; 
  int dim ; 
  int byte_rev_flag  =  do_nothing ;
  size_t base,skip;
  int corr_stride_copy, corr_stride_mom_op;
  int copy_pt,oper_pt;
  int where,where_onemom;
  char myname[] = "read_3pt_onemom";
  
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
    printf("%s: error %d reading the form factor HEADER from %s\n",
	   myname,errno,filename);
    exit(1);
  }




  /***** check the header ****/
  if( header_data[0] != magic_number )
  {
    printf("%s: ERROR: magic number mismatch between code %d and file %d\n",
	   myname,magic_number, header_data[0]  );
    printf("I will try byte reversing the data\n"); 
    byte_rev_flag = do_byte_rev ;
  }

  if(  byte_rev_flag == do_byte_rev    )
  {
    byte_rev_array( header_data, header_size); 
  }



  if( header_data[0] != magic_number )
  {
    printf("%s: ERROR: magic number mismatch between code %d and file %d\n",
	   myname,magic_number, header_data[0]  );
    exit(1); 
  }

  if( header_data[1] != version_number )
  {
    printf("%s: ERROR: version number mismatch between code %d and file %d\n",
	   myname,version_number, header_data[1]);
    exit(1); 
  }


  *nt = header_data[2]  ;
  check_sum_in = header_data[3] ;
  *no_p_values = header_data[4] ;
  *no_q_values = header_data[5] ;
  *no_oper = header_data[6] ;
  *no_spectator = header_data[7]   ;
  *no_zonked = header_data[8]  ;
  dim = header_data[9]  ;
  *hl_flag  = header_data[10]  ;
  *no_sequential = header_data[11]  ;
  *nocopies = header_data[12]  ;

  /*******  read some more data from the file **********/

  /*** read the external momentum from the  file ***/
  howmany = (size_t) 3*(*no_q_values) ;
  if( ( (*q_momstore) = (int32type *)calloc( howmany , sizeof(int32type))   )  == NULL) 
  {
    printf("%s: There was an error in allocating q-momentum \n",myname);
    exit(1);
  }

  if( fread((*q_momstore)  ,sizeof(int32type),howmany,fp) != howmany   )
  {
    printf("%s: error %d during reading the q-momentum table from %s\n",
	   myname,errno,filename);
    exit(1);
  }

  if(  byte_rev_flag == do_byte_rev    )
  {
    byte_rev_array((*q_momstore) , (int) howmany );
  }




  /*** read the external momentum from the  file ***/
  howmany = (size_t) 3*(*no_p_values) ;
  if( ( (*p_momstore) = (int32type *)calloc( howmany , sizeof(int32type))   )  == NULL) 
  {
    printf("%s: There was an error in allocating p--momentum \n",myname);
    exit(1);
  }

  if( fread((*p_momstore)  ,sizeof(int32type),howmany,fp) != howmany   )
  {
    printf("%s: error %d reading the p-momentum table from %s\n",
	   myname,errno,filename);
    exit(1);
  }

  if(  byte_rev_flag == do_byte_rev    )
  {
    byte_rev_array((*p_momstore) , (int) howmany );
  }

  /*** reserve memory for the three point function ****/
  /* REDUCE THE STATED DIMENSION BY THE NUMBER OF P-VALUES BECAUSE
     WE ARE SELECTING ONLY ONE */
  dim /= *no_p_values;
  nobj = (size_t) dim ; 
  if( ((*corr) = (complex *)calloc( nobj , sizeof(complex) ) )  == NULL) 
  {
    printf("%s: There was an error in allocating the three point function, nobj = %d\n",
	   myname,nobj);
    exit(1);
  }

  /* Stride per momentum, operator, and copy in corr array */
  corr_stride_mom_op = 
    (*nt)*(*no_zonked)*(*no_sequential)*(*no_spectator)*(*no_q_values);

  /* Stride for all values in one copy in corr array */
  corr_stride_copy = 
    corr_stride_mom_op*(*no_p_values)*(*no_oper);

  /* Get current file pointer location (at beginning of data) */
  base = ftell(fp);
  
  /* The momentum index is 6th in the corr array, so we
     block on the first 5 indices, select the 6th, and 
     iterate explicitly over the remaining 2 indices */

  for(copy_pt = 0; copy_pt < *nocopies; copy_pt++)
    for(oper_pt = 0; oper_pt < *no_oper; oper_pt++){

      /** locate the required correlator in the file corr array **/

      where = corr_stride_copy*copy_pt + 
	corr_stride_mom_op*(p_mom_select + (*no_p_values)*oper_pt);

      /** locate the required correlator in the RAM corr array **/

      where_onemom = corr_stride_copy*copy_pt/(*no_p_values) +
	corr_stride_mom_op*oper_pt;

      /** read a portion of the the three point functions from disk *****/
      nobj = (size_t) corr_stride_mom_op ; 
      
      /* Position the file for reading */
      skip = where*sizeof(complex);
      fseek(fp,base+skip,SEEK_SET);
      
      /** read the three point function to disk *****/
      if( fread(&((*corr)[where_onemom]),sizeof(complex),nobj,fp) != nobj   )
	{
	  printf("%s: error %d reading the form factor data from %s\n",
		 myname,errno,filename);
	  exit(1);
	}

      if(  byte_rev_flag == do_byte_rev    )
	{
	  byte_rev_array((int32type*)(&((*corr)[where_onemom])), nobj*2 );
	}
    }

  /*** calculate the checksum of the form factors ****/
  /**  check_sum =  bsd_sum((char*)(*corr), sizeof(complex)*dim)  ;

  if( check_sum !=  check_sum_in )
  {
    printf("%s: ERROR: checksum mismatch bewteen the file %d and the code %d \n",
    myname,check_sum,check_sum_in) ; 
    printf("%s: ERROR: checksum mismatch bewteen the file %d and the code %d \n",
    myname,check_sum,check_sum_in) ; 
    printf("%s: ERROR: checksum mismatch bewteen the file %d and the code %d \n",
    myname,check_sum,check_sum_in) ; 

  } **/


  /*** close the file ****/
  if( fclose(fp) != 0 )
  {
    printf("%s: error %d closing %s \n",myname,errno,filename);
    exit(1);
  }

  /* Create operator list */
  if( ((*corr_oper_list) = (int *)calloc( *no_oper, sizeof(int) ) )  == NULL) 
    {
      printf("%s: There was an error in allocating the operator list\n",
	     myname);
      exit(1);
    }

  for(i = 0; i < *no_oper; i++)
    (*corr_oper_list)[i] = i;

  /* Create copy list */
  if( ((*corr_copy_list) = (int *)calloc( *nocopies, sizeof(int) ) )  == NULL) 
    {
      printf("%s: There was an error in allocating the copy list\n",myname);
      exit(1);
    }

  for(i = 0; i < *nocopies; i++)
    (*corr_copy_list)[i] = i;

  /* Now pretend that only one momentum was read from this file */
  (*p_momstore)[0] = (*p_momstore)[3*p_mom_select];
  (*p_momstore)[1] = (*p_momstore)[3*p_mom_select+1];
  (*p_momstore)[2] = (*p_momstore)[3*p_mom_select+2];
  *no_p_values = 1;

  printf("I have read three point functions for B momentum %d %d %d from %s\n",
	 (*p_momstore)[0],(*p_momstore)[1],
	 (*p_momstore)[2],filename);
}

/* Same as above, but reads data for only one B meson momentum AND
   only one spectator kappa */
void read_3pt_onemom_onespect(
  complex **corr,           /* the two point correlators */
  int **corr_oper_list,     /* List of operator numbers read */
  int **corr_copy_list,     /* List of copy numbers read */
  int32type **p_momstore,      /* list of B meson momenta */
  int32type **q_momstore,      /* list of momentum transfers */
  int *nt,                  /* time values */
  int *no_p_values,         /* momentum values */
  int *no_q_values,         /* momentum values */
  int *no_oper,             /* number of operators */
  int *no_spectator,        /* the number of  kappa values */
  int *no_sequential,       /* the number of  kappa values */
  int *no_zonked ,          /* the number of  kappa values */
  int *hl_flag,             /* flag identifying type of correlator */
  char filename[80],        /* the name of the disk file */
  int *nocopies,            /* rotations */
  int p_mom_select,         /* Momentum index selected */
  int spect_select          /* Spectator kappa index selected */
)
{
  FILE *fp ;
  size_t nobj ; 
  const int32type magic_number = 14567332 ; 
  const int32type version_number  = 1  ; /** update this flag when the data format changes ***/
  int i ;
  size_t  name_len = 80 ;
  size_t howmany ;
  int check_sum ; 
  int check_sum_in ; 
  int dim ; 
  int byte_rev_flag  =  do_nothing ;
  size_t base,skip;
  int corr_stride_copy, corr_stride_mom_op, corr_stride_spect_op;
  int copy_pt,oper_pt,q_pt;
  int where,where_ram;
  char myname[] = "read_3pt_onemom";
  
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
    printf("%s: error %s reading the form factor HEADER from %s\n",
	   myname,errno,filename);
    exit(1);
  }




  /***** check the header ****/
  if( header_data[0] != magic_number )
  {
    printf("%s: ERROR: magic number mismatch between code %d and file %d\n",
	   myname,magic_number, header_data[0]  );
    printf("I will try byte reversing the data\n"); 
    byte_rev_flag = do_byte_rev ;
  }

  if(  byte_rev_flag == do_byte_rev    )
  {
    byte_rev_array( header_data, header_size); 
  }



  if( header_data[0] != magic_number )
  {
    printf("%s: ERROR: magic number mismatch between code %d and file %d\n",
	   myname,magic_number, header_data[0]  );
    exit(1); 
  }

  if( header_data[1] != version_number )
  {
    printf("%s: ERROR: version number mismatch between code %d and file %d\n",
	   myname,version_number, header_data[1]);
    exit(1); 
  }


  *nt = header_data[2]  ;
  check_sum_in = header_data[3] ;
  *no_p_values = header_data[4] ;
  *no_q_values = header_data[5] ;
  *no_oper = header_data[6] ;
  *no_spectator = header_data[7]   ;
  *no_zonked = header_data[8]  ;
  dim = header_data[9]  ;
  *hl_flag  = header_data[10]  ;
  *no_sequential = header_data[11]  ;
  *nocopies = header_data[12]  ;

  /*******  read some more data from the file **********/

  /*** read the external momentum from the  file ***/
  howmany = (size_t) 3*(*no_q_values) ;
  if( ( (*q_momstore) = (int32type *)calloc( howmany , sizeof(int32type))   )  == NULL) 
  {
    printf("%s: There was an error in allocating q-momentum \n",myname);
    exit(1);
  }

  if( fread((*q_momstore)  ,sizeof(int32type),howmany,fp) != howmany   )
  {
    printf("%s: error %d during reading the q-momentum table from %s\n",
	   myname,errno,filename);
    exit(1);
  }

  if(  byte_rev_flag == do_byte_rev    )
  {
    byte_rev_array((*q_momstore) , (int) howmany );
  }




  /*** read the external momentum from the  file ***/
  howmany = (size_t) 3*(*no_p_values) ;
  if( ( (*p_momstore) = (int32type *)calloc( howmany , sizeof(int32type))   )  == NULL) 
  {
    printf("%s: There was an error in allocating p--momentum \n",myname);
    exit(1);
  }

  if( fread((*p_momstore)  ,sizeof(int32type),howmany,fp) != howmany   )
  {
    printf("%s: error %d reading the p-momentum table from %s\n",
	   myname,errno,filename);
    exit(1);
  }

  if(  byte_rev_flag == do_byte_rev    )
  {
    byte_rev_array((*p_momstore) , (int) howmany );
  }

  /*** reserve memory for the three point function ****/
  /* REDUCE THE STATED DIMENSION BY THE NUMBER OF P-VALUES AND
     THE NUMBER OF SPECTATOR VALUES BECAUSE
     WE ARE SELECTING ONLY ONE */
  dim /= *no_p_values;
  dim /= *no_spectator;
  nobj = (size_t) dim ; 
  if( ((*corr) = (complex *)calloc( nobj , sizeof(complex) ) )  == NULL) 
  {
    printf("%s: There was an error in allocating the three point function, nobj = %d\n",
	   myname,nobj);
    exit(1);
  }

  /* The indexing is as follows:
     linear(t,zonk_pt,seq_pt,spect_pt,q_pt,p_pt,oper_pt,copy) =
     t + nt*(zonk_pt + no_zonked*(seq_pt + no_sequential*(spect_pt +
         no_spectator*(q_pt + no_q_values*(p_pt + 
         no_p_values*(oper_pt + no_oper*copy_pt))))))
  */

  /* Stride between spectator quark values */
  /* linear(0,0,0,spect_pt+1,0,0,0,0) - linear(0,0,0,spect_pt,0,0,0,0) */
  corr_stride_spect_op = 
    (*nt)*(*no_zonked)*(*no_sequential);

  /* Get current file pointer location (at beginning of data) */
  base = ftell(fp);

  /* The spect_pt index is 4th in the corr array, and the momentum
     index is 6th so we block on the 1st 3 indices, select the 4th,
     iterate explicitly over the 5th, select the 6th and iterate over
     the 7th and 8th. */

  for(copy_pt = 0; copy_pt < *nocopies; copy_pt++)
    for(oper_pt = 0; oper_pt < *no_oper; oper_pt++)
      for(q_pt = 0; q_pt < *no_q_values; q_pt++){

	/** locate the required correlator in the file corr array **/

	where = corr_stride_spect_op*(spect_select +
          (*no_spectator)*(q_pt + (*no_q_values)*(p_mom_select + 
	   (*no_p_values)*(oper_pt + (*no_oper)*copy_pt))));

	/** locate the required correlator in the RAM corr array **/
	/* Because no_spectator = no_p_values = 1 the indexing
	   is as follows:

	   linear(t,zonk_pt,seq_pt,spect_pt,q_pt,p_pt,oper_pt,copy) =
	   t + nt*(zonk_pt + no_zonked*(seq_pt + 
           no_sequential*(q_pt + no_q_values*(oper_pt + no_oper*copy_pt))))
	*/

	where_ram = corr_stride_spect_op*(q_pt + 
   	     (*no_q_values)*(oper_pt + (*no_oper)*copy_pt));

	/** read a portion of the the three point functions from disk *****/
	nobj = (size_t) corr_stride_spect_op ; 
	
	/* Position the file for reading */
	skip = where*sizeof(complex);
	fseek(fp,base+skip,SEEK_SET);
	
	/** read the three point function to disk *****/
	if( fread(&((*corr)[where_ram]),sizeof(complex),nobj,fp) != nobj   )
	  {
	    printf("%s: error %d reading the form factor data from %s\n",
		   myname,errno,filename);
	    exit(1);
	  }
	
	if(  byte_rev_flag == do_byte_rev    )
	  {
	    byte_rev_array((int32type*)(&((*corr)[where_ram])), nobj*2 );
	  }
      }
  
  /*** calculate the checksum of the form factors ****/
  /**  check_sum =  bsd_sum((char*)(*corr), sizeof(complex)*dim)  ;
       
  if( check_sum !=  check_sum_in )
  {
    printf("%s: ERROR: checksum mismatch bewteen the file %d and the code %d \n",
    myname,check_sum,check_sum_in) ; 
    printf("%s: ERROR: checksum mismatch bewteen the file %d and the code %d \n",
    myname,check_sum,check_sum_in) ; 
    printf("%s: ERROR: checksum mismatch bewteen the file %d and the code %d \n",
    myname,check_sum,check_sum_in) ; 

  } **/


  /*** close the file ****/
  if( fclose(fp) != 0 )
    {
      printf("%s: error %d closing %s \n",myname,errno,filename);
      exit(1);
    }
  
  /* Create operator list */
  if( ((*corr_oper_list) = (int *)calloc( *no_oper, sizeof(int) ) )  == NULL) 
    {
      printf("%s: There was an error in allocating the operator list\n",
	     myname);
      exit(1);
    }
  
  for(i = 0; i < *no_oper; i++)
    (*corr_oper_list)[i] = i;
  
  /* Create copy list */
  if( ((*corr_copy_list) = (int *)calloc( *nocopies, sizeof(int) ) )  == NULL) 
    {
      printf("%s: There was an error in allocating the copy list\n",myname);
      exit(1);
    }
  
  for(i = 0; i < *nocopies; i++)
    (*corr_copy_list)[i] = i;

  /* Now pretend that only one momentum was read from this file */
  (*p_momstore)[0] = (*p_momstore)[3*p_mom_select];
  (*p_momstore)[1] = (*p_momstore)[3*p_mom_select+1];
  (*p_momstore)[2] = (*p_momstore)[3*p_mom_select+2];

  printf("I have read three point functions for B momentum %d %d %d and spectator index %d from %s\n",
	 (*p_momstore)[0],(*p_momstore)[1],
	 (*p_momstore)[2],spect_select,filename);
}



