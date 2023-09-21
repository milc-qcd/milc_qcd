/*************** gauge_utilities_includes.h ******************************/
/* Include files and prototypes common to most files in the application  */

/* Include files */
#include "../include/config.h"  /* Keep this first */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "defines.h"
#include "../include/complex.h"
#include "../include/su3.h"
#include "lattice.h"
#include "../include/macros.h"
#include "../include/comdefs.h"
#include "../include/io_lat.h"
#include "../include/generic.h"
#include "../include/dirs.h"

/* Prototypes for functions in high level code */
int setup();
int readin(int prompt);
void run_gradient_flow();
void flow_step();

//void stout_step_rk();
void staple();
void fmunu_fmunu(double *time, double *space, double *charge);
void initialize_integrator();
void gauge_action_w_s( double *wl1x1s, double *wl1x1t,
                       double *wl1x2s, double *wl1x2t );
// various integrators, compile-time choice
void integrate_RK_2N();
void integrate_RKMK3();
void integrate_RKMK_generic();
void integrate_adapt_RK_2N();
void integrate_adapt_bs();

// some flow helper functions
void clear_anti_hermitian( anti_hermitmat *dest );
void set_identity( su3_matrix *dest );
void mult_ah_su3_nn( anti_hermitmat *a, su3_matrix *b, su3_matrix *dest );
void scalar_mult_ah( anti_hermitmat *a, Real c, anti_hermitmat *dest );
void scalar_mult_add_ah( anti_hermitmat *a, anti_hermitmat *b, Real c,
                         anti_hermitmat *dest );
void anti_hermitian_traceless_proj( su3_matrix *a, anti_hermitmat *dest );
void exp_anti_hermitian( anti_hermitmat *a, su3_matrix *dest, int n );
void ahmat_copy( anti_hermitmat *a, anti_hermitmat *b );
void commutator_ah( anti_hermitmat *a, anti_hermitmat *b, anti_hermitmat *c );
void dexpinv( anti_hermitmat *u, anti_hermitmat *v, int q, anti_hermitmat *d );
Real su3mat_distance( su3_matrix *a, su3_matrix *b );

#ifdef DEBUG_FIELDS
void repack_site( su3_matrix *a, MATRIX_TYPE *b );
void dump_double_lattice();
#endif


#ifdef SPHALERON
void run_gradient_flow_region( int region_flag );
void flow_step_region();
void staple_region( int region_flag );

// various integrators, compile-time choice
void integrate_RK_2N_region();
void integrate_RKMK3_region();
void integrate_RKMK_generic_region();
void integrate_adapt_RK_2N_region();
void integrate_adapt_bs_region();

void clear_links ( su3_matrix **link, size_t nlink );
void destroy_links ( su3_matrix **link );
su3_matrix ** new_links( int region_flag, size_t nlink );

su3_matrix ** new_half_links( void );
void renew_half_links( su3_matrix ** link_half   );

void clear_last_flow_links ( su3_matrix *** link_last_flow  );
void destroy_last_flow_links ( su3_matrix *** link_last_flow  );
su3_matrix *** new_last_flow_links( void );
void update_last_flow_links( su3_matrix *** link_last_flow );


void clear_fieldstrength ( su3_matrix **fieldstrength, size_t ncomp );
void destroy_fieldstrength ( su3_matrix **fieldstrength );
su3_matrix ** new_fieldstrength ( size_t ncomp );


void 
make_fieldstrength_region( int region_flag,
                          field_offset link_src, 
                          su3_matrix **fieldstrength );

void 
make_fieldstrength_bdry ( field_offset link_src, 
                          su3_matrix ***link_last_flow, 
                          su3_matrix **fieldstrength );

void 
make_fieldstrength_bulk ( field_offset link_src, 
                          su3_matrix **fieldstrength );

void 
make_fieldstrength_half ( field_offset link_src, 
                          su3_matrix **fieldstrength );

void 
make_fieldstrength_full ( field_offset link_src, 
                          su3_matrix **fieldstrength );


void 
fmunu_fmunu_bdry( su3_matrix ***link_last_flow, 
                  double *time, double *space, double *charge );

void
fmunu_fmunu_bulk( double *time, double *space, double *charge );

void
fmunu_fmunu_half( double *time, double *space, double *charge );

void
fmunu_fmunu_full( double *time, double *space, double *charge );


void 
gauge_action_w_s_bdry( su3_matrix ***link_last_flow,  
                       double *wl1x1s, double *wl1x1t,
                       double *wl1x2s, double *wl1x2t );

void 
gauge_action_w_s_bulk( double *wl1x1s, double *wl1x1t,
                       double *wl1x2s, double *wl1x2t );

void 
gauge_action_w_s_half( double *wl1x1s, double *wl1x1t,
                       double *wl1x2s, double *wl1x2t );

void
gauge_action_w_s_full( double *wl1x1t, double *wl1x1s,
                       double *wl1x2t, double *wl1x2s );


void prepare_bulk_links( void );
void report_bulk( Real time );
void spatial_blocking( void );

#ifdef DEBUG_BLOCKING
void test_blocking ( void );
void make_test_lat (void);
#endif

#endif