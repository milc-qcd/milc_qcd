/*  Sun Jul  6 08:17:16 MDT 1997
    UPDATE the comments for they are out of date

   For the heavy-light and heavy-heavy project, calculate the ratio of 
   two and three point functions, that is proportional to the 
   required matrix element.

   See for example UKQCD,Phys. REv D52, 5067, 1995.
   equation (22).


  matrix_element(t)  = c3(t,tf)
                       c2_L(t) * c2_hqet(tf - t)

  This code calculates matrix_element(t) with the jackknife errors.
  The two point function depends on velocity and momentum


  MORE WORK, Fri Mar 28 14:15:43 MST 1997 :: 
    1) Think about the operator pointer
              for the two and three point functions.

    2) Think about time reversal symmetrry
    3) Think about effective mass plos for the two pint functions

**/


#include "read_hl_form.h"

/**** function protypes ***/

int load_filenames(char filename[80], char file_list[MAX_NO_FILE][80], int maxfile) ;


void read_input_param(int *w_p , int *w_zonked , int *w_spectate , int *w_q , 
int *w_oper_2pt, int *w_oper_3pt , int *w_seq ,
char two_names[MAX_NO_FILE][80] , char three_names[MAX_NO_FILE][80] ,
char seq_names[MAX_NO_FILE][80] , int *nofile , 
char filename[80] ) ;



void read_propagating_form_corr(complex **corr, int **corr_oper_list,
  int **corr_copy_list,
  int **p_momstore, int **q_momstore,
  int *nt, int *no_p_values, int *no_q_values, int *no_oper,
  int *no_spectator,   int *no_sequential, int *no_zonked , int *hl_flag,
  char filename[80] )  ;


void write_data(Real *data,Real *data_err, int dim ,char filename[80] ) ;

Real jack_sample_twopt_light(int jomit, int nosample, complex *two_corr[MAX_NO_FILE ] ,
			int t, int nt , int q_pt, int no_q_values, int w_spectate, int no_spectator,
			int zonked_pt, int no_zonked_light, int oper_pt) ;


Real jack_sample_twopt_hl(int jomit, int nosample, complex *two_corr[MAX_NO_FILE ] ,
			int t, int nt , int q_pt, int no_q_values, int w_spectate, int no_spectator,
			int zonked_pt, int no_zonked_light, int oper_pt)  ;




Real jack_sample_threept_hh(int jomit, int nosample, complex *three_corr[MAX_NO_FILE ] ,
			int t, int nt , int q_pt, int no_q_values, 
			int p_pt, int no_p_values, 
			int w_spectate, int no_spectator,
                        int w_seq, int no_sequential , 
			int zonked_pt, int no_zonked_heavy, int oper_pt);


Real jack_sample_threept_hl(int jomit, int nosample, complex *three_corr[MAX_NO_FILE ] ,
			int t, int nt , int q_pt, int no_q_values, 
			int p_pt, int no_p_values, 
			int w_spectate, int no_spectator,
                        int w_seq, int no_sequential , 
			int zonked_pt, int no_zonked_light, int oper_pt ) ;




Real jackerror(Real mean,Real *data,int nodata) ;
Real jackmean(Real *data,int nodata,int jomit) ;

char *three_oper_name(int n ) ;
char *two_oper_name(int n ) ;


/**** end of function protypes ***/


int main(int argc, char *argv[])
{
  const int max_no_file  = MAX_NO_FILE ; 
  int nofile ; 
  int nofile_tmp  ; 

  char two_names[MAX_NO_FILE][80] ;   
  char three_names[MAX_NO_FILE][80] ; 
  char seq_names[MAX_NO_FILE][80] ; 

  char input_two[]   = "input_2pt" ; 
  char input_three[] = "input_3pt" ; 
  char input_seq[]   = "input_seq_pt" ; 
  char fileout[80]   = "matrix_out" ; 

  complex *three_corr[MAX_NO_FILE ] ; /*** The store of the 3pt correlators ****/
  complex *seq_corr[MAX_NO_FILE ] ; /*** The store of the seq correlators ****/
  complex *two_corr[MAX_NO_FILE ] ; /*** The store of the 2pt correlators ****/
  int *three_corr_oper_list;
  int *seq_corr_oper_list;
  int *two_corr_oper_list;
  int *three_corr_copy_list;
  int *seq_corr_copy_list;
  int *two_corr_copy_list;

  Real j_sample[MAX_NO_FILE ]  ; 
  int  j_denom ;
  enum j_denom_options { all_non_zero , some_zero }  ;
  int sample ; 
  Real j_mean ;

  int *p_momstore;
  int *q_momstore;

  Real top, bot ;
  int nt ;
  int no_p_values ;
  int no_q_values ;
  int no_oper_3pt ; 
  int no_oper_2pt ; 
  int no_oper_seq ; 
  int no_spectator ;
  int no_zonked  ;
  int no_sequential  ;
  int hl_flag ,  hl_flag_seq , hl_flag_two ;
  int t , trev ;
  int no_k_one , no_k_two ; 
  int seq_no_k_one , seq_no_k_two ; 
  int w_seq; 

  int w_p , w_zonked , w_spectate , w_q , tf , tf_two ; 
  int w_oper_2pt , w_oper_3pt ;

  char input_param[80] = "input_param" ;  
  int i ;

  int dim_elem ; 
  Real *matrix_elem ; 
  Real *matrix_elem_err ; 

  int maxname  ;
/**  char *three_oper_name[ MAX_NAME ] ;
  char *two_oper_name[ MAX_NAME ] ;  **/

  /*** print some generic titles ******/

  printf("============================================================\n"); 
  printf("Plot the propagating matrix element with jackknife errors \n") ; 
  printf("============================================================\n"); 

  if( argc == 1 ) 
  {
    printf("By default I will read the INPUT parameters from %s\n",input_param); 
  }
  else if( argc == 2 ) 
  {
    strcpy(input_param, argv[1]);
    printf("I will read the INPUT parameters from %s\n",input_param); 
  }
  else
  {
    printf("ERROR:: too many arguments argc = %d \n",argc) ; 
    exit(1); 
  }



  /*** load the names of the files containing the data **/
  read_input_param( &w_p , &w_zonked , &w_spectate , &w_q ,
		   &w_oper_2pt, &w_oper_3pt, &w_seq , 
		   two_names,three_names,seq_names,&nofile , input_param  ) ; 
  

  /*** more work::: move somewhere ******/
  printf("The number of input files = %d\n",nofile); 


  /*** read all the three and two point data from disk *****/
  for(i = 0 ; i < nofile ; ++i)
  {
    read_propagating_form_corr(&three_corr[i] , &three_corr_oper_list,
			&three_corr_copy_list,
			&p_momstore, &q_momstore,
			&nt, &no_p_values, &no_q_values, &no_oper_3pt,
			&no_spectator, &no_sequential, &no_zonked, &hl_flag,
		        three_names[i]) ;

    if( i == 0 )
    {
      tf = nt /2 -1 ; /*** more work ***/


      printf("\tParameters read from the the file %s\n",three_names[i] ) ; 
      printf("Number of timeslices = %d\n",nt); 
      printf("The EXTERNAL mesons are fixed at 0 and %d\n",tf); 
      printf("Total number of operator momentum P values = %d\n",no_p_values );
      printf("Total number of operator momentum Q values = %d\n",no_q_values );
      printf("The number of spectator kappa values = %d\n", no_spectator); 
      printf("The number of zonked kappa values = %d\n", no_zonked); 
      printf("The number of sequential kappa values = %d\n", no_sequential); 
      printf("P momentum = %d %d %d\n",
	     *(p_momstore + 3*w_p), *(p_momstore + 1 + 3*w_p), *(p_momstore + 2 + 3*w_p) ); 
      printf("Q momentum = %d %d %d\n",
	     *(q_momstore + 3*w_q), *(q_momstore + 1 + 3*w_q), *(q_momstore + 2 + 3*w_q) ); 
	    
	    


      if(hl_flag == HEAVY_TO_LIGHT) 
      {
	printf("*****> Heavy-light three point function <*****\n"); 

      }
      else if (hl_flag == HEAVY_TO_HEAVY) 
      {
	printf("*****> Heavy-heavy three point function <*****\n"); 

      }
      else
      {
	printf("ERROR: hl_flag = %d is out of range\n",hl_flag); 
	exit(10);
      }


      /***  more work, print the choosen velocity and momentum ***/


    }

    free(q_momstore  ) ; /*** this is dumb *****/
    free(p_momstore) ;      /*** this is dumb *****/

    /***  load the two point functions ***/

    read_prop_twopt(&two_corr[i], &two_corr_oper_list,
		    &two_corr_copy_list,
		    two_names[i], &p_momstore,
		    &hl_flag_two, &no_k_one , &no_k_two,  &no_oper_2pt,
		    &nt , &no_q_values) ;

		       
    free(q_momstore  ) ; /*** this is dumb *****/

    /*** more work, need to check the parameters ****/

    if( i== 0 )
    {
      printf("Three point operator = %s\n", three_oper_name(w_oper_3pt) ) ; 
      printf("Two point operator = %s\n", two_oper_name(w_oper_2pt) ) ; 
    }


    /***  load the sequential two point functions *******/

    read_prop_twopt(&seq_corr[i], &seq_corr_oper_list,
		    &seq_corr_copy_list,
		    seq_names[i], &q_momstore,
		    &hl_flag_seq, &no_k_one , &no_k_two,  &no_oper_seq,
		    &nt , &no_q_values) ;





    free(q_momstore) ;      /*** this is dumb *****/

    /*** more work, need to check the parameters ****/


  }

  printf("..... All the data has been read in\n"); 

  /****************** 
    calculate the full and jackknife samples of the ratio of 2 and 3 functions 
    *****************/

  /*** just look at the full sample for now ****/
  dim_elem = nt/2 ;  /*** more work ***/

  if( ( (matrix_elem) = calloc( dim_elem , sizeof(int))   )  == NULL) 
  {
    printf("There was an error in allocating \"matrix_elem\" \n");
    exit(1);
  }



  if( ( (matrix_elem_err) = calloc( dim_elem , sizeof(int))   )  == NULL) 
  {
    printf("There was an error in allocating \"matrix_elem_err\" \n");
    exit(1);
  }



  /***
    calculate the full sample correlators 
   ***/

  for( t = 0 ; t < dim_elem ; ++t)
  {
    trev = tf - t ;

    if(hl_flag == HEAVY_TO_LIGHT) 
    {
      top  = jack_sample_threept_hl(nofile , nofile , three_corr,
				    t, nt , w_q, no_q_values, w_p, no_p_values, 
				    w_spectate, no_spectator, w_seq ,  no_sequential ,
				    w_zonked, no_zonked,w_oper_3pt);


      bot = jack_sample_twopt_hl(nofile, nofile, seq_corr,
			      trev, nt , w_p, no_p_values, w_spectate, 
				 no_spectator,w_zonked, no_zonked, w_oper_2pt) * 
				jack_sample_twopt_light(nofile,nofile, two_corr, t, nt , w_q, no_q_values, 
						  w_spectate, no_spectator,
							w_zonked, no_zonked, w_oper_2pt) ;

      }
      else if (hl_flag == HEAVY_TO_HEAVY) 
      {
	top  = jack_sample_threept_hh(nofile , nofile , three_corr,
				   t, nt , w_q, no_q_values, w_p, no_p_values, 
				      w_spectate, no_spectator, w_seq , no_sequential ,
				   w_zonked, no_zonked,w_oper_3pt);


	/***** this section requires more thought ******/

      bot = jack_sample_twopt_hl(nofile, nofile, seq_corr,
			      trev, nt , w_p, no_p_values, w_spectate, 
				 no_spectator,w_zonked, no_zonked, w_oper_2pt) * 
				jack_sample_twopt_light(nofile,nofile, two_corr, t, nt , w_q, no_q_values, 
						  w_spectate, no_spectator,
							w_zonked, no_zonked, w_oper_2pt) ;


      }



    if( bot != 0 )
      matrix_elem[t] = top / bot ;
    else 
      matrix_elem[t] = -1 ;



    /***** jackknife analysis ******/
    j_denom = all_non_zero ; 

    for( sample = 0 ; sample < nofile ; ++sample )
    {

      if(hl_flag == HEAVY_TO_LIGHT) 
	{
	  top  = jack_sample_threept_hl(sample , nofile , three_corr,
					t, nt , w_q, no_q_values, w_p, no_p_values, 
					w_spectate, no_spectator, w_seq , no_sequential ,
					w_zonked, no_zonked,w_oper_3pt);


	  bot = jack_sample_twopt_hl(sample, nofile, seq_corr,
				     trev, nt , w_p, no_p_values, w_spectate, 
				     no_spectator,w_zonked, no_zonked, w_oper_2pt) * 
				       jack_sample_twopt_light(nofile,nofile, two_corr, t, nt , w_q, no_q_values, 
							       w_spectate, no_spectator,
							       w_zonked, no_zonked, w_oper_2pt) ;

	}
      else if (hl_flag == HEAVY_TO_HEAVY) 
	{
	  top  = jack_sample_threept_hh(sample , nofile , three_corr,
					t, nt , w_q, no_q_values, w_p, no_p_values, 
					w_spectate, no_spectator, w_seq , no_sequential ,
					w_zonked, no_zonked,w_oper_3pt);


	  /***** this section requires more thought ******/

	  bot = jack_sample_twopt_hl(sample, nofile, seq_corr,
				     trev, nt , w_p, no_p_values, w_spectate, 
				     no_spectator,w_zonked, no_zonked, w_oper_2pt) * 
				       jack_sample_twopt_light(nofile,nofile, two_corr, t, nt , w_q, no_q_values, 
							       w_spectate, no_spectator,
							       w_zonked, no_zonked, w_oper_2pt) ;
      }




    if( bot != 0 )
     j_sample[ sample ] = top / bot ;
    else 
      j_denom =   some_zero ;


    }  /*** end the loop over samples ****/

    if( j_denom ==  all_non_zero  )
    {
      j_mean = 0.0 ; 
      for(sample = 0 ; sample < nofile ;++sample)
	j_mean +=  j_sample[ sample ] ; 

      j_mean /= nofile ;

      matrix_elem_err[t ] = jackerror( j_mean ,  j_sample,nofile ) ;

    }
    else
    {
      matrix_elem_err[t] = 0.0 ; 
    }





  }  /*** end of the loop over time ****/






  /**** write the data to disk file ****/

  /*** create the filename ****/
  /*** more work, need to target the filename with useful information ***/
  write_data(matrix_elem ,matrix_elem_err  , dim_elem , fileout); 


  /*** free up the reserved memory in the code ****/
  free(matrix_elem); 
  free(matrix_elem_err); 


  return 0 ;
}
