/* 
 *  Read a SPECIFIED two point function from a file
 *
 *  NOTE this is scalar code -- it should only be run on 
 *  one node.
 *
 *   
 */


#include "read_hl_form.h"
#include "prop_form_utilities.h"
#include <errno.h>

/* Functions declarations used in the routines in this file ***/

void read_select_twopt(
  complex **corr,      /* pointer to complex array: the two point correlator */
  char filename[80],   /* the name of the disk file */
  int nselect,         /* Number of sets of selection parameters */
  twopt_select *twoselect, /* Selection parameter sets */
  int32type **q_momselect,   /* momentum list */		       
  int *hl_flag,        /* identifies meson source/sink smearing */
  int *ntime           /* number of time steps */
  )
{
  size_t nobj ; 
  FILE *fp ;
  int32type magic_number = 44451189 ; 
  int32type version_number  = 1  ; /** update this flag when the data format changes ***/
  int i ;
  int t ;
  int ngot ;
  int jselect ;
  size_t  name_len = 80 ;
  size_t howmany ;
  size_t where,skip,corr_stride,base;
  int check_sum ; 
  int check_sum_in ; 
  int dim ; 
  int byte_rev_flag  =  do_nothing ;
  int nt,no_q_values,no_spectator,no_zonked,no_oper,nocopies;
  int32type *q_momstore;
  int zonked_pt;       /* zonked quark index */
  int spect_pt;        /* spectator quark index */
  int q_pt;            /* momentum index */
  int oper_pt;         /* operator index */
  int copy_pt;         /* copy index */
  Real wt;            /* weight */
  complex *corr_tmp;   /* temporary storage for correlator */
  char myname[] = "read_select_twopt";

  /* Memory for some header information ***/
#define HEADER_DIM_WRITE_CORR 11
  size_t  header_size = (size_t) HEADER_DIM_WRITE_CORR ;
  int32type header_data[HEADER_DIM_WRITE_CORR] ;
#undef HEADER_DIM_WRITE_CORR

  /**** open the file ******/
  if( (fp = fopen(filename ,"rb")) == NULL )
  {
    fprintf(stderr,"%s ERROR: Could not open the file %s\n",
	   myname,filename);
    exit(1);
  }

  /** read the header information ***/
  if( fread(header_data,sizeof(int32type),header_size,fp) != header_size )
  {
    fprintf(stderr,"%s Error %d reading the HEADER from %s\n",
	   myname,errno,filename);
    exit(1);
  }

  /*** unpack and check the header ******/

  /***** check the header ****/
  if( header_data[0] != magic_number )
  {
    byte_rev_flag = do_byte_rev ;
  }


  if(  byte_rev_flag == do_byte_rev    )
  {
    byte_rev_array( header_data, header_size); 
  }


  /***** check the header ****/
  if( header_data[0] != magic_number )
  {
    fprintf(stderr,"%s ERROR: magic number mismatch between code %x and file %x\n",
	   myname,magic_number, header_data[0]  );
    exit(1); 
  }


  if( header_data[1] != version_number )
  {
    fprintf(stderr,"%s ERROR: version number mismatch between code %d and file %d\n",
	    myname,version_number, header_data[1]);
    exit(1); 
  }


  nt  = header_data[2] ;
  *ntime = nt;
  check_sum_in  = header_data[3]  ;
  no_q_values = header_data[4]  ;
  no_spectator = header_data[5]  ;
  no_zonked = header_data[6]  ;
  no_oper = header_data[7]  ;
  dim = header_data[8]  ; 
  *hl_flag = header_data[9]  ;
  nocopies = header_data[10]  ;

  /*** read the external q--momentum from the file ***/
  howmany = (size_t) 3*(no_q_values) ;
  if( ( q_momstore = (int32type *)calloc( howmany , sizeof(int32type))   )  == NULL) 
  {
    fprintf(stderr,"%s Can't malloc q-momentum \n",myname);
    exit(1);
  }


  if( fread(q_momstore,sizeof(int32type),howmany,fp) != howmany   )
  {
    fprintf(stderr,"%s: error %d reading the Q momentum table from %s\n",
	   myname,errno,filename);
    exit(1);
  }

  if(  byte_rev_flag == do_byte_rev    )
  {
    byte_rev_array(q_momstore , (int) howmany );
  }

  /* Pretend we read only one three-vector - and only the first in the list */
  *q_momselect = &q_momstore[3*twoselect[0].mom];

  /** read the two point functions from disk *****/
  nobj = (size_t) nt ; 
  if( ((*corr) = (complex *)calloc( nobj , sizeof(complex) ) )  == NULL) 
  {
    fprintf(stderr,"%s: Can't calloc corr, nobj = %d\n",myname,nobj);
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
      zonked_pt       = twoselect[jselect].other;
      spect_pt        = twoselect[jselect].spect;
      q_pt            = twoselect[jselect].mom;
      oper_pt         = twoselect[jselect].oper;
      copy_pt         = twoselect[jselect].copy;
      wt              = twoselect[jselect].wt;

      printf(" WT %f ZK %d SP %d K %d %d %d OP %s\n",
	     wt,zonked_pt,spect_pt,
	     q_momstore[3*q_pt],q_momstore[3*q_pt+1],q_momstore[3*q_pt+2],
	     oper_name_and_corrections(oper_pt,copy_pt));

      /** locate the required correlator on disk **/

      corr_stride = TWOPT_FORM_WHERE(nt,no_zonked-1 ,no_spectator-1,
				     no_q_values-1, no_oper-1) ; ; 
      where = corr_stride*copy_pt + 
	TWOPT_FORM_WHERE(0,zonked_pt ,spect_pt,q_pt, oper_pt) ;


      skip = where*sizeof(complex);

      if(fseek(fp,base+skip,SEEK_SET) != 0)
	{
	  fprintf(stderr,"%s: Error %d seeking the two point function on %s\n",
		  myname,errno,filename);
	  exit(1);
	}
      
      /** read the two point function from disk *****/
      if( (ngot = fread(corr_tmp,sizeof(complex),nobj,fp)) != nobj   )
	{
	  printf("%s: error %d reading %d form factor data items from %s. Got %d.\n",
		 myname,errno,nobj,filename,ngot);
	  printf("Error after seeking %d bytes\n",base+skip);
	  exit(1);
	}
      
      
      if(  byte_rev_flag == do_byte_rev    )
	{
	  byte_rev_array((int32type*) (corr_tmp), nobj*2 );
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
    fprintf(stderr,"%s error %d closing %s \n",myname,errno,filename);
    exit(1);
  }

  free(corr_tmp);

  printf("2pt file: %s\n",filename);

}


