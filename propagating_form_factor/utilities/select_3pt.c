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

void read_select_form_corr(
  complex **corr,      /* pointer to complex array: the two point correlator */
  char filename[80],   /* the name of the disk file */
  int nselect,         /* Number of sets of selection parameters */
  threept_select *threeselect, /* Selection parameter sets */
  int32type **p_momselect,   /* momentum list */	
  int32type **q_momselect,   /* momentum list */	
  int *hl_flag,        /* identifies meson source/sink smearing */
  int *ntime           /* number of time steps */
  )
{
  FILE *fp ;
  size_t nobj, ngot ; 
  int32type magic_number = 14567332 ; 
  int32type version_number  = 1  ; /** update this flag when the data format changes ***/
  int i ;
  int t ;
  int jselect;
  size_t  name_len = 80 ;
  size_t howmany ;
  size_t where,skip,corr_stride,base;
  int check_sum ; 
  int check_sum_in ; 
  int dim ; 
  int byte_rev_flag  =  do_nothing ;
  int nt,no_q_values,no_p_values,no_spectator,no_sequential;
  int no_zonked,no_oper,nocopies;
  int32type *q_momstore,*p_momstore;
  int zonked_pt;       /* zonked quark index */
  int seq_pt;          /* sequential quark index */
  int spect_pt;        /* spectator quark index */
  int q_pt;            /* momentum index */
  int p_pt;            /* momentum index */
  int oper_pt;         /* operator index */
  int copy_pt;         /* copy index */
  Real wt;            /* weight */
  complex *corr_tmp;   /* temporary storage for correlator */
  char myname[] = "read_select_form_corr";
  
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
    printf("%s: Error %d reading header of size %d from %s\n",
	   myname,errno,header_size,filename);
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


  nt = header_data[2]  ;
  *ntime = nt;
  check_sum_in = header_data[3] ;
  no_p_values = header_data[4] ;
  no_q_values = header_data[5] ;
  no_oper = header_data[6] ;
  no_spectator = header_data[7]   ;
  no_zonked = header_data[8]  ;
  dim = header_data[9]  ;
  *hl_flag  = header_data[10]  ;
  no_sequential = header_data[11]  ;
  nocopies = header_data[12]  ;

  /*******  read some more data from the file **********/

  /*** read the external momentum from the  file ***/
  howmany = (size_t) 3*(no_q_values) ;
  if( ( q_momstore = (int32type *)calloc( howmany , sizeof(int32type))   )  == NULL) 
  {
    printf("There was an error in allocating q-momentum \n");
    exit(1);
  }

  if( fread((q_momstore)  ,sizeof(int32type),howmany,fp) != howmany   )
  {
    printf("There was an error during the reading of the q-momentum table from %s\n",filename);
    exit(1);
  }

  if(  byte_rev_flag == do_byte_rev    )
  {
    byte_rev_array((q_momstore) , (int) howmany );
  }




  /*** read the external momentum from the  file ***/
  howmany = (size_t) 3*(no_p_values) ;
  if( ( (p_momstore) = (int32type *)calloc( howmany , sizeof(int32type))   )  == NULL) 
  {
    printf("There was an error in allocating p--momentum \n");
    exit(1);
  }

  if( fread((p_momstore)  ,sizeof(int32type),howmany,fp) != howmany   )
  {
    printf("There was an error during the reading of the p-momentum table from %s\n",filename);
    exit(1);
  }

  if(  byte_rev_flag == do_byte_rev    )
  {
    byte_rev_array((p_momstore) , (int) howmany );
  }

  /* Pretend we read only one three-vector - and only the first in the list */
  *q_momselect = &q_momstore[3*threeselect[0].q];

  /* Pretend we read only one three-vector - and only the first in the list */
  *p_momselect = &p_momstore[3*threeselect[0].p];

  /*** reserve memory for the three point function ****/
  nobj = (size_t) nt ; 
  if( ((*corr) = (complex *)calloc( nobj , sizeof(complex) ) )  == NULL) 
  {
    printf("There was an error in allocating the three point function, nobj = %d\n",nobj);
    exit(1);
  }

  if( (corr_tmp = (complex *)calloc( nobj , sizeof(complex) ) )  == NULL) 
  {
    printf("There was an error in allocating corr_tmp, nobj = %d\n",nobj);
    exit(1);
  }

  for(t = 0; t < nt; t++)
    {
      (*corr)[t].real = (*corr)[t].imag = 0.;
    }

  /* Iterate over sets of selection parameters, reading data for each
     and accumulating results in the correlation array */

  printf("Selection from %s\n",filename);

  base = ftell(fp);
  for(jselect = 0; jselect < nselect; jselect++)
    {
      /* Unload selection parameters */
      
      spect_pt  = threeselect[jselect].spect;
      zonked_pt = threeselect[jselect].zonked;
      seq_pt    = threeselect[jselect].seq;
      q_pt      = threeselect[jselect].q;
      p_pt      = threeselect[jselect].p;
      oper_pt   = threeselect[jselect].oper;
      copy_pt   = threeselect[jselect].copy;
      wt        = threeselect[jselect].wt;
      
      printf(" WT %f SP %d ZK %d SQ %d Q %d %d %d P %d %d %d OP %s\n",
	     wt,spect_pt,zonked_pt,seq_pt,
	     q_momstore[3*q_pt],q_momstore[3*q_pt+1],q_momstore[3*q_pt+2],
	     p_momstore[3*p_pt],p_momstore[3*p_pt+1],p_momstore[3*p_pt+2],
	     three_oper_name(oper_pt,copy_pt));

      /** locate the required correlator on disk **/
      
      corr_stride = FORM_WHERE(nt,no_zonked-1,no_sequential-1,
	     no_spectator-1,no_q_values-1,no_p_values-1,
	     no_oper-1) ; 
      where = corr_stride*copy_pt + 
	FORM_WHERE(0,zonked_pt,seq_pt,spect_pt,q_pt,p_pt,oper_pt) ;
      
      skip = where*sizeof(complex);
      
      if(fseek(fp,base+skip,SEEK_SET) != 0)
	{
	  printf("There was an error seeking the two point function on %s\n",
		 filename);
	  exit(1);
	}
      
      
      /** read the three point function from disk *****/
      if( (ngot = fread(corr_tmp,sizeof(complex),nobj,fp)) != nobj   )
	{
	  printf("%s: error %d reading %d form factor data items from %s. Got %d.\n",
		 myname,errno,nobj,filename,ngot);
	  printf("Error after seeking %d bytes\n",base+skip);
	  exit(1);
	}
      
      if(  byte_rev_flag == do_byte_rev    )
	{
	  byte_rev_array((int32type*) corr_tmp, nobj*2 );
	}

      /* Accumulate result in correlation array */
      for(t = 0; t < nt; t ++)
	{
	  (*corr)[t].real += corr_tmp[t].real*wt;
	  (*corr)[t].imag += corr_tmp[t].imag*wt;
	}

    } /* jselect */


  /*** close the file ****/
  if( fclose(fp) != 0 )
  {
    printf("There was an error during the closing of %s \n",filename);
    exit(1);
  }

  free(corr_tmp);

  printf("3pt file: %s\n",filename);

}


