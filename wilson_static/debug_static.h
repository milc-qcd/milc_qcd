#ifndef DEBUG_STATIC_INCLUDE
#define DEBUG_STATIC_INCLUDE

/*
 *  Include files for the prototypes of the debug functions
 *
 */



void random_smear_func();

void dump_gauge();

void test_gauge_config();

void dump_strip_quark();

Real ran1(long *idum)  ;

void dump_vary_matrix(Real *vary_matrix);

void dump_static_prop();


void dump_smear_line();


void dump_smear_func();



void write_smear_func() ;

void dump_convolve();


void test_gamma_left();


void copy_w_right_vec_mat(wilson_vector *v, wilson_matrix *m, int spin, int color);
void copy_w_right_mat_vec(wilson_matrix *m, wilson_vector *v, int spin, int color) ;

void dump_complex_array(complex *data, int nodata);
void dump_psi();
void dump_psi_smear();

void dump_smearedmeson(complex *data);

void time_gauge_write() ;

void light_quark_pion(int flag);
void random_wilson_vec(int colour, int spin); 

void check_calc_matrix() ;
void zero_w_line() ;
void zero_strip_quark() ;

/******* end of the include file for the debug functions ******/

#endif
