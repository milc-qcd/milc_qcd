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
/* setup.c */
int setup();
int readin(int prompt);
void initialize_integrator();



#ifndef REGIONS
//void stout_step_rk();
/* integrate.c */
void run_gradient_flow();
void flow_step();
/* staple.c */
void staple();
#else
/* integrate_region.c */
void run_gradient_flow( int region_flag );
/* staple_region.c */
#ifndef USE_FIELD
void staple( int region_flag );
#else
void
staple( int region_flag, su3_matrix **stap );
#endif

/* integrate.c or integrate_region.c */
void flow_step();
// various integrators, compile-time choice
void integrate_RK_2N();
void integrate_RKMK3();
void integrate_RKMK_generic();
void integrate_adapt_RK_2N();
void integrate_adapt_bs();

#endif

// printing observables
void print_observables( char *TAG, double flowtime, 
                        double *Et_WS, double *Es_WS, 
                        double *Et_C, double *Es_C, 
                        double *charge );
// generic helper functions
void clear_field ( su3_matrix **this, size_t ncomp );
void destroy_field ( su3_matrix ***this );
su3_matrix ** new_field( size_t ncomp );
su3_matrix ** new_links_from_site( int region_flag, size_t nlink );

// some flow helper functions
void clear_anti_hermitian_field ( anti_hermitmat **this, size_t ncomp );
void destroy_anti_hermitian_field ( anti_hermitmat ***this );
anti_hermitmat ** new_anti_hermitian_field( size_t ncomp );
void destroy_anti_hermitian_twodim_field ( anti_hermitmat ****this );
anti_hermitmat *** new_anti_hermitian_twodim_field( size_t dim1, size_t dim2 );

// more flow helper functions
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

// blocking
#ifdef BLOCKING
void spatial_blocking( void );
#ifdef DEBUG_BLOCKING
void test_blocking ( void );
void make_test_lat (void);
#endif
#endif
// automatize setup of blocked directions
void setup_blocked_dirs( int *dir, int *dirb, int *diro );

#ifndef REGIONS
// making generic/regionalized fieldstrengths
/* fmunu.c */
void fmunu_fmunu(double *time, double *space, double *charge);
void gauge_action_w_s( double *wl1x1s, double *wl1x1t,
                       double *wl1x2s, double *wl1x2t );
#endif
#ifdef REGIONS
// regionalized observables
/* fieldstrength.c and fmunu_region.c */
// making generic/regionalized fieldstrengths
void 
make_fieldstrength_region( int region_flag,
                          su3_matrix **fieldstrength );
void 
make_fieldstrength_bulk ( su3_matrix **fieldstrength );
void 
make_fieldstrength_full ( su3_matrix **fieldstrength );

// making generic/regionalized fmunu squared and topological charges
void
fmunu_fmunu_bulk( double *time, double *space, double *charge );
void
fmunu_fmunu_full( double *time, double *space, double *charge );

// making generic/regionalized gauge loops
void 
gauge_action_w_s_bulk( double wlt[2], double wls[2] );
void
gauge_action_w_s_full( double wlt[2], double wls[2] );

#endif


#ifdef SPHALERON
// general helpers for SPHALERON application
void prepare_bulk_links( void );
void report_bulk( Real time, Real *q_bulk );

void bulk_flow( Real *q_bulk );
void bdry_flow( Real *q_bulk );

// making half-integer time links
su3_matrix ** new_half_links( void );
void renew_half_links( su3_matrix ** link_half   );

// making last flow time links
void update_last_flow_links( su3_matrix ** link_last_flow );

// making special fieldstrengths
void 
make_fieldstrength_bdry ( field_offset link_src, 
                          su3_matrix **link_last_flow, 
                          su3_matrix **fieldstrength );
void 
make_fieldstrength_lwr_bdry ( field_offset link_src, 
                          su3_matrix **link_last_flow, 
                          su3_matrix **fieldstrength );
void 
make_fieldstrength_upr_bdry ( field_offset link_src, 
                          su3_matrix **link_last_flow, 
                          su3_matrix **fieldstrength );
void 
make_fieldstrength_half ( field_offset link_src, 
                          su3_matrix **fieldstrength );

// making special fmunu squared and topological charges
void 
fmunu_fmunu_bdry( su3_matrix **link_last_flow, 
                  double *time, double *space, double *charge );
void
fmunu_fmunu_half( double *time, double *space, double *charge );

// making special gauge loops
void 
gauge_action_w_s_lwr_bdry( su3_matrix **link_last_flow,  
                       double wlt[2], double wls[2] );
void 
gauge_action_w_s_upr_bdry( su3_matrix **link_last_flow,  
                       double wlt[2], double wls[2] );
void 
gauge_action_w_s_bdry( su3_matrix **link_last_flow,  
                       double wlt[2], double wls[2] );
void 
gauge_action_w_s_half( double wlt[2], double wls[2] );
#endif