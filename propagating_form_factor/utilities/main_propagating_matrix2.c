/*  5/29/98 C.D. Reads only specified data from specified
    correlator files, rather than the whole file */

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
#include "prop_form_utilities.h"

/**** function protypes ***/

void write_data(Real *data,Real *data_err, int dim ,char filename[80] ) ;


/**** end of function protypes ***/

Real matel_ratio(int t, int tf, int nt, int jack_pt, int nofile, 
		  int *zeros, int top_sign, int top_phase,
		  int top_fb, int bot1_fb, int bot2_fb,
		  complex *three_corr[MAX_NO_FILE], 
		  complex *seq_corr[MAX_NO_FILE], 
		  complex *rec_corr[MAX_NO_FILE])
{
  complex ff3_forw,ff3_back;
  Real top_forw,top_back;
  Real bot1_forw, bot2_forw, ratio_forw ;
  Real bot1_back, bot2_back, ratio_back ;
  Real temp;

  int trev ;

  *zeros = 0;
  
  trev = tf - t ;
  ff3_back.real = ff3_back.imag = top_back = bot1_back = bot2_back = 0.;
  
  ff3_forw  = jack_sample_threept(jack_pt, nofile, three_corr, t);
  bot1_forw = jack_sample_twopt(jack_pt, nofile, seq_corr, trev).real;
  bot2_forw = jack_sample_twopt(jack_pt, nofile, rec_corr, t).real;

  if(t > 0 && t < nt/2)
    {
      ff3_back  = jack_sample_threept(jack_pt, nofile, three_corr, nt - t);
      bot1_back = jack_sample_twopt(jack_pt, nofile, seq_corr, nt - trev).real;
      bot2_back = jack_sample_twopt(jack_pt, nofile, rec_corr, nt - t).real;
    }
  else
    {
      ff3_back = ff3_forw;
      bot1_back = bot1_forw;
      bot2_back = bot2_forw;
    }

  printf("%2d top %.3g %.3g %.3g %.3g bot1 %.3g %.3g bot2 %.3g %.3g\n",
	 t,ff3_forw.real,ff3_back.real,ff3_forw.imag,ff3_back.imag,bot1_forw,bot1_back,bot2_forw,bot2_back);


  if(bot1_fb == FOLD)
    {
      bot1_back = (bot1_back + bot1_forw)/2.;
      bot1_forw = bot1_back;
    }
  else if (bot1_fb == BACKWARD)
    {
      temp = bot1_forw;
      bot1_forw = bot1_back;
      bot1_back = temp;
    }
  if(bot2_fb == FOLD)
    {
      bot2_back = (bot2_back + bot2_forw)/2.;
      bot2_forw = bot2_back;
    }
  else if (bot2_fb == BACKWARD)
    {
      temp = bot2_forw;
      bot2_forw = bot2_back;
      bot2_back = temp;
    }

  if(top_phase==1 || top_phase==3)top_sign = -top_sign;
  if(top_phase==0 || top_phase==2)top_forw = ff3_forw.real; 
  else top_forw = ff3_forw.imag;
  if(top_phase==0 || top_phase==2)top_back = ff3_back.real; 
  else top_back = ff3_back.imag;
  if(top_phase==2 || top_phase==3)top_forw = -top_forw;
  if(top_phase==2 || top_phase==3)top_back = -top_back;
  

  if( bot1_forw != 0 && bot2_forw != 0 )
    ratio_forw = top_forw / (bot1_forw * bot2_forw) ;
  else 
    {
      ratio_forw = 0;
      *zeros = 1;
    }
  
  if( bot1_back != 0 && bot2_back != 0 )
    ratio_back = top_back/(bot1_back * bot2_back);
  else 
    {
      ratio_back = 0;
      *zeros = 1;
    }
  
  if(top_fb == FORWARD)  return ratio_forw;
  if(top_fb == BACKWARD) return ratio_back;
  if(top_fb == FOLD)     return (ratio_forw + top_sign*ratio_back)/2.;

}    

int main(int argc, char *argv[])
{
  const int max_no_file  = MAX_NO_FILE ; 
  int nofile ; 
  int nofile_tmp  ; 

  three_list threept;
  two_list twopt_recoil;
  two_list twopt_sequential;

  char fileout[80]   = "matrix_out" ; 

  complex *three_corr[MAX_NO_FILE ] ; /*** The store of the 3pt correlators ****/
  complex *seq_corr[MAX_NO_FILE ] ; /*** The store of the seq correlators ****/
  complex *rec_corr[MAX_NO_FILE ] ; /*** The store of the 2pt correlators ****/

  Real j_sample[MAX_NO_FILE ]  ; 
  int  j_denom ;
  enum j_denom_options { all_non_zero , some_zero }  ;
  int sample ; 
  Real j_mean ;

  int32type *p_momstore;
  int32type *q_momstore;

  int nt ;
  int hl_flag ,  hl_flag_seq , hl_flag_two ;
  int t ;
  int tf ; 
  int q_pt, p_pt, zeros;
  int top_sign, top_phase;

  char input_param[80] = "input_param" ;  
  int i,j ;

  int dim_elem ; 
  Real *matrix_elem ; 
  Real *matrix_elem_err ; 

  int maxname  ;


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


  /*** load file names and selection parameters **/
  read_input_param(input_param, &threept, &twopt_recoil, &twopt_sequential);


  /*** more work::: move somewhere ******/
  printf("The number of three point files = %d\n",
	 threept.nfile); 
  nofile = threept.nfile;
  if(nofile != twopt_recoil.nfile)
    printf("WARNING: The number of two point recoil files = %d\n",
	 twopt_recoil.nfile); 

  if(nofile != twopt_sequential.nfile)
    printf("WARNING: The number of two point sequential files = %d\n",
	   twopt_sequential.nfile); 


  /*** read all the three and two point data from disk *****/
  for(i = 0 ; i < threept.nfile ; ++i)
  {
    /*** read and reserve memory for the three point function ******/
    read_select_form_corr(&three_corr[i], threept.filename[i], 
			  threept.nselect, threept.select, 
			  &p_momstore, &q_momstore,
			  &hl_flag , &nt );
    
    if( i == 0 )
    {
      tf = nt/2 ; /*** more work ***/

      printf("\tParameters read from the the file %s\n",threept.filename[i]); 
      printf("Number of timeslices = %d\n",nt); 
      printf("The EXTERNAL mesons are fixed at 0 and %d\n",tf); 
      /* Needs more work: only first selection reported here */
      p_pt = threept.select[0].p;
      q_pt = threept.select[0].q;
      printf("P momentum = %d %d %d\n",
	     *(p_momstore + 3*p_pt), *(p_momstore + 1 + 3*p_pt), *(p_momstore + 2 + 3*p_pt) ); 
      printf("Q momentum = %d %d %d\n",
	     *(q_momstore + 3*q_pt), *(q_momstore + 1 + 3*q_pt), *(q_momstore + 2 + 3*q_pt) ); 
	    
	    
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


      /*** more work, need to check the parameters ****/


      printf("Three op is %d\n",threept.select[0].oper);
      printf("Three point operator = %s\n", 
	     three_oper_name(threept.select[i].oper,
			     threept.select[0].copy) ) ; 
      printf("Sequential two point operator = %s\n", 
	     two_oper_name(twopt_sequential.select[i].oper,
			   twopt_sequential.select[0].copy) ) ; 
      printf("Recoil two point operator = %s\n", 
	     two_oper_name(twopt_recoil.select[i].oper,
			   twopt_recoil.select[0].copy) ) ; 
    }


    /***  load the recoil two point function ***/

    read_select_twopt(&rec_corr[i], twopt_recoil.filename[i],  
		      twopt_recoil.nselect, twopt_recoil.select,
		      &q_momstore, &hl_flag_two , &nt );

    /***  load the sequential two point function *******/

    read_select_twopt(&seq_corr[i], twopt_sequential.filename[i],  
		      twopt_sequential.nselect, twopt_sequential.select,
		      &q_momstore, &hl_flag_seq , &nt );
		       

    /*** more work, need to check the parameters ****/


  }

  printf("..... All the data has been read in\n"); 

  /****************** 
    calculate the full and jackknife samples of the ratio of 2 and 3 functions 
    *****************/

  /*** just look at the full sample for now ****/
  dim_elem = nt/2 ;  /*** more work ***/

  if( ( (matrix_elem) = (Real *)calloc( dim_elem , sizeof(Real))   )  == NULL) 
  {
    printf("There was an error in allocating \"matrix_elem\" \n");
    exit(1);
  }



  if( ( (matrix_elem_err) = (Real *)calloc( dim_elem+1, sizeof(Real))   )  == NULL) 
  {
    printf("There was an error in allocating \"matrix_elem_err\" \n");
    exit(1);
  }



  /***
    calculate the full sample correlators 
   ***/

  top_sign = oper_chirality(threept.select[0].oper);
  top_phase = oper_phase(threept.select[0].oper);

  for( t = 0 ; t <= dim_elem ; ++t)
    {
      matrix_elem[t] = matel_ratio(t,tf,nt,nofile,nofile,&zeros,
				   top_sign,top_phase,
				   threept.forwback, 
				   twopt_sequential.forwback,
				   twopt_recoil.forwback,
				   three_corr,seq_corr,rec_corr);
      
      /***** jackknife analysis ******/
      j_denom = all_non_zero ; 
      
      for( sample = 0 ; sample < nofile ; ++sample )
	{
	  j_sample [ sample ] = matel_ratio(t,tf,nt,sample,nofile, &zeros,
					    top_sign,top_phase,
					    threept.forwback, 
					    twopt_sequential.forwback,
					    twopt_recoil.forwback,
					    three_corr,seq_corr,rec_corr);
	  if(zeros) j_denom = some_zero;
	  
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
