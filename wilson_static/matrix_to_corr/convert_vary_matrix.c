/*    $Header: /lqcdproj/detar/cvsroot/milc_qcd/wilson_static/matrix_to_corr/convert_vary_matrix.c,v 1.1 2005/02/23 00:06:11 detar Exp $
 *  This is a scalar utility program to read the 
 *  static variational matrix from disk and
 *  and create a file of LS and SS correlators, suitable for 
 *  use in the effective mass code.
 *  
 *  This program takes the name of the matrix
 *  file as a command line argument; and the source number.
 *
 *  WARNING :: there is an assumption that the local operator
 *  is the first in the set.
 *
 */


#include LATDEF

int main(int argc, char *argv[])
{
  FILE *fp ;
  int nt,nosmear ;
  size_t nobj   ;
  Real *vary_matrix ;
  Real *coeff ;
  char vary_out[80] ;
  int i ;
  const int local_src = 0 ;
  char src_name[NAME_LEN] ;
  size_t  name_len = NAME_LEN ;
  char coeff_file[NAME_LEN] ;
  int which_source ; 
  int arg_pt ;
  enum which_corr { NOT_SET = -10 , SINGLE_SRC , VARY_SRC} ;
  int what_src = NOT_SET ;

  int what_fold = NO_TIME_AVERAGE ; /** default is not to average over time ***/
   int byte_rev_flag ; 

  /*********..........**********........*****************/

  printf("=============================================================\n");
  printf(" Extraction of LS and SS correlators from variational matrix \n");
  printf("=============================================================\n");

  if( argc < 3  )
  {
    dump_usage(argv[0]) ;
  }
  strcpy(vary_out,argv[1]);
  which_source = atoi(argv[2]);

  /** parse the input arguments *******/
  arg_pt = 2 ;
  while( arg_pt < argc )
  {

    if( strcmp(argv[arg_pt],"-s" ) == 0 )
    {
      if( what_src != NOT_SET ) dump_usage(argv[0]) ;
      which_source = atoi(argv[arg_pt + 1]);
      arg_pt += 2 ;
      what_src = SINGLE_SRC  ;
    }
    else if( strcmp(argv[arg_pt],"-c" ) == 0 )
    {
      if( what_src != NOT_SET ) dump_usage(argv[0]) ;
      strcpy(coeff_file ,argv[arg_pt+1]);
      arg_pt += 2 ;
      what_src = VARY_SRC  ;
    }
    else if( strcmp(argv[arg_pt],"-fold" ) == 0 )
    {
      what_fold = AVERAGE_TIME  ;
      arg_pt += 1 ;
    }
    else
    {
      dump_usage(argv[0]);
    }


  }


  /**** open the file and read the header  ******/
  fp = read_vary_header(vary_out, local_src,which_source,
			src_name,&nosmear, &nt,&byte_rev_flag ) ;


  printf("The number of timeslices = %d\n",nt);
  printf("The number of smearing functions = %d\n",nosmear);

  /* reserve memory for the variational matrix *****/
  nobj        =  (size_t) nt*nosmear*nosmear  ;
  vary_matrix =  (Real *)  calloc( nobj , sizeof(Real) );

  /** read the variational matrix from disk *****/
  if( fread(vary_matrix,sizeof(Real),nobj,fp) != nobj   )
  {
    printf("There was an error during the reading of the variational matrix \n");
    exit(1);
  }

  /*** close the file ****/
  if( fclose(fp) != 0 )
  {
    printf("There was an error during the closing of %s \n",vary_out);
    exit(1);
  }

  /***  byte_reverse if required ***/

  if(  byte_rev_flag == do_byte_rev    )
  {
    byte_rev_array((int*) vary_matrix, (int) nobj); 
  }




  printf("I have read the variational matrix from the file %s\n",vary_out);

  /* write the correlatorsto file *****/ 
  if( what_src == SINGLE_SRC  )
  {
    if( which_source < 0 || which_source >= nosmear)
    {
      printf("ERROR which_source = %d is out of range 0 to %d\n",which_source,nosmear);
      exit(2);
    }
    printf("smearing_function_file  = %s\n",src_name);

    dump_LSandSS_corr(vary_matrix,which_source, vary_out, local_src,what_fold,nt, nosmear);
  }
  else if( what_src == VARY_SRC  )
  {
    coeff = (Real *) calloc( (size_t) nosmear, sizeof(Real) );
    read_coeff_file(coeff_file ,coeff, nosmear);
    dump_vary_LSandSS_corr(vary_matrix,coeff, vary_out, local_src,what_fold,nt, nosmear);
    free(coeff);
  }
  else
  {
    printf("Bad option chosen \n");
    dump_usage(argv[0]);
  }
  
  if(what_fold == AVERAGE_TIME   )
  {
    printf("\nThe correlators will be averaged over positive and negative times\n");
  }


  /* Free up the memory used in the calculation ***/
  free(vary_matrix);

  return 0 ;

}

/* Usage output for bad input ****/
void dump_usage(char exec_file[])
{
    /** this should be a function call *****/
    printf("usage::  %s  [matrix data file] \n",exec_file);
    printf("                               -s [smearing function number]\n");
    printf("                               -c [coefficient file]\n");
    exit(1);


}


/*
 *  Write to file the LS and SS orrelators of a
 *  specific source.
 *
 *  Subroutine arguments::
 *
 *   vary_matrix[] :: the smearing matrix
 *   which_source  :: number of source to for the smearing function
 *   filebase[]    :: name of output file
 *   local_src     :: location of the local operator
 *   what_fold     :: a flag which controls whether positive and negative
 *                    correlators are averaged together
 *   nt            :: tyhe number of time slices
 *   nosmear       :: the number of smearing functions
 */


void dump_LSandSS_corr(Real *vary_matrix,int which_source, char file_base[], int local_src, int what_fold, int nt, int nosmear)
{
  int t ;
  int pt ;
  char filename[NAME_LEN] ;
  Real *corr ;
  /*****------------------------------**********/

  corr = (Real *) calloc( (size_t) nt, sizeof(Real) );

  /*** calculate the LS correlators  ********/

  for( t= 0 ; t < nt ;++t)
  {
    pt = t + nt*(local_src + nosmear*(which_source)); 
    corr[t] = vary_matrix[pt] ;
  }

  sprintf(filename,"%s_LS",file_base);
  write_corr_out(filename, corr, what_fold, nt);



  /*** extract the SS correlators  ********/

  for( t= 0 ; t < nt ;++t)
  {
    pt = t + nt*(which_source+ nosmear*(which_source)); 
    corr[t] = vary_matrix[pt] ;
  }

  sprintf(filename,"%s_SS",file_base);
  write_corr_out(filename, corr, what_fold, nt);


  free(corr);


}






/*
 *  Write to file the LS and SS orrelators of a
 *  source which is a linear combination of the basis wave functions.
 *
 *  Subroutine arguments::
 *
 *   vary_matrix[] :: the smearing matrix
 *   coeff         :: vector coefficients of the required source
 *   filebase[]    :: name of output file
 *   local_src     :: location of the local operator
 *   what_fold     :: a flag which controls whether positive and negative
 *                    correlators are averaged together
 *   nt            :: tyhe number of time slices
 *   nosmear       :: the number of smearing functions
 */


void dump_vary_LSandSS_corr(Real *vary_matrix,Real *coeff , char file_base[], int local_src, int what_fold,int nt, int nosmear)
{
  int t ;
  int pt ;
  char filename[NAME_LEN] ;
  Real *corr ;
  int ismear,iss,ils ;
  /*****------------------------------**********/

  corr = (Real *) calloc( (size_t) nt, sizeof(Real) );

  /*** calculate the LS correlators  ********/


  for( t= 0 ; t < nt ;++t)
  {
    corr[t] = 0.0 ;
    for(ismear=0 ; ismear < nosmear ; ++ismear)
    {
      pt = t + nt*(local_src + nosmear*(ismear)); 
      corr[t] += coeff[ismear]*vary_matrix[pt] ;
    }
  }

  sprintf(filename,"%s_LS",file_base);
  write_corr_out(filename, corr, what_fold,nt);



  /*** extract the SS correlators  ********/

  for( t= 0 ; t < nt ;++t)
  {
    corr[t] = 0.0 ;
    for(iss=0 ; iss < nosmear ; ++iss)
      for(ils=0 ; ils < nosmear ; ++ils)  
      {
	pt = t + nt*(ils+ nosmear*(iss)); 
	corr[t] += vary_matrix[pt] * coeff[ ils ] * coeff[ iss ]  ;

/**	printf("DEBUG ils=%d ss=%d m= %g corr = %g coeff_ls = %g coeff_ss = %g \n",  
        ils,iss,vary_matrix[pt],corr[t],coeff[ ils ],coeff[ iss ]    );  ***/


      }
  }

  sprintf(filename,"%s_SS",file_base);
  write_corr_out(filename, corr,what_fold, nt);


  free(corr);


}


/*
 *  Write an individual source correlators to a 
 *  file.
 *
 *  dim            :: the amount of data
 *  what_fold      :: a flag which controls whether positive and negative
 *                    correlators are averaged together
 *  filename       :: the name of the output file 
 *  corr[0..dim-1] :: the correlators
 */


void write_corr_out(char filename[], Real *corr, int what_fold, int dim)
{
  FILE *fp ;
  int t ;
  Real corr_av ;

  /**** open the file ******/
  if( (fp = fopen(filename ,"w")) == NULL )
  {
    printf("Could not open the file %s\n",filename);
    exit(1);
  }


  if( what_fold == NO_TIME_AVERAGE )
  {
    for(t=0; t < dim  ; ++t)
      fprintf(fp,"%e\n",corr[t] );
  }
  else if ( what_fold == AVERAGE_TIME  )
  {
    fprintf(fp,"%e\n",corr[0] );

    for(t=1; t < dim/2  ; ++t)
    {
      corr_av = 0.5*(corr[t] + corr[dim - t] ) ;
      fprintf(fp,"%e\n",corr_av );
    }

  }



  /*** close the file ****/
  if( fclose(fp) != 0 )
  {
    printf("There was an error during the closing of %s \n",filename);
    exit(1);
  }




  printf("I have written a correlator to the file %s\n",filename);


}


/*
 *  Read the coefficients for the variational source from a 
 *  simple file
 *
 *  Subroutine arguments
 *   filename[]      :: the name of the input file
 *   dim             :: the number of coffecients
 *   coeff[0..dim-1] :: the coffecients
 */


void read_coeff_file(char filename[], Real *coeff, int dim)
{
  FILE *fp ;
  int i ;
  Real coeff_tmp  ;
  int pt ;

  /** zsero the coefficeint store ****/
  for(i=0 ; i < dim ;++i)
    *(coeff + i) = 0.0 ;


  /**** open the file ******/
  if( (fp = fopen(filename ,"r")) == NULL )
  {
    printf("Could not open the file %s\n",filename);
    exit(1);
  }

  /** read the coefficients ***/
  while( fscanf(fp,"%d %lfHELP",&pt,&coeff_tmp) == 2  )
  {
    if( pt < 0 || pt >= dim)
    {
      printf("read_coeff_file:: index %d out of range in %s\n",pt,filename);
      exit(1);
    }
    else
    {
      *(coeff + pt ) = coeff_tmp ;
    }


  }


  /*** close the file ****/
  if( fclose(fp) != 0 )
  {
    printf("There was an error during the closing of %s \n",filename);
    exit(1);
  }

  printf("I have read the VARIATIONAL coefficients from the file %s\n",filename);

  printf("Here are the coefficients of the variational source\n");
  for(i=0 ; i < dim ; ++i)
  {
    printf("coeff[ %d ] = %g \n",i,coeff[i]);
  }


  


}

