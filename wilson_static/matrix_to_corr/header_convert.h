#ifndef HEADER_CONVERT_INC
#define HEADER_CONVERT_INC

/*
 *  Header file for the static variational conversion routine
 */

/**
 **  >>>>>>>>>>  include files <<<<<<<<<<
 **/


#include <stdio.h>
#include <math.h>
#include<string.h>
#include<stdlib.h>


/*
 *   >>>>>> function prototypes <<<<<<
 */


void dump_LSandSS_corr(Real *vary_matrix,int which_source, char file_base[], int local_src, int what_fold, int nt, int nosmear);

void dump_vary_LSandSS_corr(Real *vary_matrix,Real *coeff , char file_base[], int local_src, int what_fold,int nt, int nosmear);

void write_corr_out(char filename[], Real *corr, int what_fold, int dim);

void dump_usage(char exec_file[]);

void read_coeff_file(char filename[], Real *coeff, int dim);

FILE *(read_vary_header(char filename[], int local_src, int which_source, char src_name[],
int *nosmear, int *nt,  int *byte_rev_flag  )) ;


void byte_rev_array(int buf[], int words) ;

/**
 **   Macro definitions
 **/


enum which_fold { NO_TIME_AVERAGE , AVERAGE_TIME } ;
#define NAME_LEN 80 

enum byte_rev_option { do_byte_rev = 10  , do_nothing  } ; 



/** end of the "header_convert" file ***/

#endif
