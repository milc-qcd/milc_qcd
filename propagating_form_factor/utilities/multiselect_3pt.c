/* 
 *  Read in the propagating heavy---> light three point functions
 *  from a disk file and append to the set of specified files
 *  following specified selection criteria
 *
 *  NOTE this is scalar code -- it should only be run on 
 *  one node.
 *
 */


#include "read_hl_form.h"
#include "prop_form_utilities.h"
#include <errno.h>

/* Functions declarations used in the routines in this file ***/
void read_multiselect_form_corr(
  char filename[80],   /* the name of the input correlator file */
  int nfile,           /* Number of output files */
  threept_oneselect *fp /* Selection parameter sets */
  )
{
  FILE *corrfp ;
  FILE *dumpfp[MAX_NO_FILE];
  complex *corr;
  size_t nobj, ngot ; 
  const int32type magic_number = 14567332 ; 
  const int32type version_number  = 1  ; /** update this flag when the data format changes ***/
  int i ;
  int ifile ;
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
  char *nextfield,*cfgid;   /* Identifying characters for binary
                              correlator file */
  char myname[] = "read_multiselect_form_corr";
  
  /* Memory for some header information ***/
#define HEADER_DIM_WRITE_CORR 13
  size_t  header_size = (size_t) HEADER_DIM_WRITE_CORR ;
  int32type header_data[HEADER_DIM_WRITE_CORR] ;

  /**** open the ASCII output files for APPENDING *****/
  for(ifile = 0; ifile < nfile; ifile++)
    {
      if( (dumpfp[ifile] = fopen(fp[ifile].filename,"a")) == NULL )
	{
	  printf("Could not open the file %s\n",fp[ifile].filename);
	  exit(1);
	}
    }

  /**** open the binary correlator file ******/
  if( (corrfp = fopen(filename ,"rb")) == NULL )
  {
    printf("Could not open the file %s\n",filename);
    exit(1);
  }

  /** read the header information ***/
  if( fread(header_data,sizeof(int32type),header_size,corrfp) != header_size )
  {
    printf("%s: Error %d reading header of size %d from %s\n",
	   myname,errno,header_size,filename);
    exit(1);
  }




  /***** check the header ****/
  if( header_data[0] != magic_number )
  {
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
  check_sum_in = header_data[3] ;
  no_p_values = header_data[4] ;
  no_q_values = header_data[5] ;
  no_oper = header_data[6] ;
  no_spectator = header_data[7]   ;
  no_zonked = header_data[8]  ;
  dim = header_data[9]  ;
/**  hl_flag  = header_data[10]  ; **/
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

  if( fread((q_momstore)  ,sizeof(int32type),howmany,corrfp) != howmany   )
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

  if( fread((p_momstore)  ,sizeof(int32type),howmany,corrfp) != howmany   )
  {
    printf("There was an error during the reading of the p-momentum table from %s\n",filename);
    exit(1);
  }

  if(  byte_rev_flag == do_byte_rev    )
  {
    byte_rev_array((p_momstore) , (int) howmany );
  }

  /*** reserve memory for the three point function ****/
  nobj = (size_t) nt ; 
  if( (corr = (complex *)calloc( nobj , sizeof(complex) ) )  == NULL) 
  {
    printf("There was an error in allocating the three point function, nobj = %d\n",nobj);
    exit(1);
  }

  if( (corr_tmp = (complex *)calloc( nobj , sizeof(complex) ) )  == NULL) 
  {
    printf("There was an error in allocating corr_tmp, nobj = %d\n",nobj);
    exit(1);
  }

  base = ftell(corrfp);

  /* Iterate over sets of selection parameters, reading data for each
     and accumulating results in the correlation array */

  for(ifile = 0; ifile < nfile; ifile++)
    {
      for(t = 0; t < nt; t++)
	{
	  corr[t].real = corr[t].imag = 0.;
	}
      
      for(jselect = 0; jselect < fp[ifile].nselect; jselect++)
	{
	  /* Unload selection parameters */
	  
	  spect_pt  = fp[ifile].select[jselect].spect;
	  zonked_pt = fp[ifile].select[jselect].zonked;
	  seq_pt    = fp[ifile].select[jselect].seq;
	  q_pt      = fp[ifile].select[jselect].q;
	  p_pt      = fp[ifile].select[jselect].p;
	  oper_pt   = fp[ifile].select[jselect].oper;
	  copy_pt   = fp[ifile].select[jselect].copy;
	  wt        = fp[ifile].select[jselect].wt;
	  
	  /** locate the required correlator on disk **/
	  
	  corr_stride = FORM_WHERE(nt,no_zonked-1,no_sequential-1,
				   no_spectator-1,no_q_values-1,no_p_values-1,
				   no_oper-1) ; 
	  where = corr_stride*copy_pt + 
	    FORM_WHERE(0,zonked_pt,seq_pt,spect_pt,q_pt,p_pt,oper_pt) ;
	  
	  skip = where*sizeof(complex);
	  
	  if(fseek(corrfp,base+skip,SEEK_SET) != 0)
	    {
	      printf("There was an error seeking the two point function on %s\n",
		     filename);
	      exit(1);
	    }
	  
	  
	  /** read the three point function from disk *****/
	  if( (ngot = fread(corr_tmp,sizeof(complex),nobj,corrfp)) != nobj   )
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
	      corr[t].real += corr_tmp[t].real*wt;
	      corr[t].imag += corr_tmp[t].imag*wt;
	    }
	  
	} /* jselect */

      /* Write the result to the ASCII file */
      /* Use the last field after the last period in the filename to
         identify the data set */
      for(nextfield=filename; nextfield!=(char *)NULL+1; 
	  nextfield=strstr(nextfield,".")+1)cfgid = nextfield;

      for(t = 0; t < nt; t++)
	fprintf(dumpfp[ifile],"%2d %12.6g %12.6g %% %s\n",
		t,corr[t].real,corr[t].imag,cfgid);

    } /* ifile */


  /*** close the file ****/
  if( fclose(corrfp) != 0 )
  {
    printf("There was an error during the closing of %s \n",filename);
    exit(1);
  }

  /**** close the ASCII output files *****/
  for(ifile = 0; ifile < nfile; ifile++)
    {
      if( fclose(dumpfp[ifile]) != 0 )
	{
	  printf("Error closing the file %s\n",fp[ifile].filename);
	  exit(1);
	}
    }

  free(corr_tmp);

}


