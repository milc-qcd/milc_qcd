/****************** hvy_qpot_includes.h ******************************/
/*
*  Include files for heavy quark potential application
*/

#ifndef HVY_QPOT_INCLUDES_H
#define HVY_QPOT_INCLUDES_H

/* Include files */
#include "../include/config.h"  /* Keep this first */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include "../include/complex.h"
#include "../include/su3.h"
#include "../include/macros.h"
#include "lattice.h"
#include "../include/comdefs.h"
#include "../include/io_lat.h"
#include "../include/generic.h"
#include "../include/dirs.h"

#define site_coord(s,mu) \
        (((short*)&(s->x))[mu])

#define NLLXUP (X3UP)
#define NLLYUP (Y3UP)
#define NLLZUP (Z3UP)
#define NLLTUP (T3UP)
#define NLLTDOWN (T3DOWN)
#define NLLZDOWN (Z3DOWN)
#define NLLYDOWN (Y3DOWN)
#define NLLXDOWN (X3DOWN)
#define OPP_NLL_DIR(dir) (23-(dir))
#define DIR_NLL(dir) ((dir)+1*NDIRS)

/* buffer handling for all hqp routines */
double * hqp_alloc_dble_buffer( int mysize );
su3_matrix * hqp_alloc_su3mat_buffer( int mysize );

void hqp_free_dble_buffer( double *buf );
void hqp_free_su3mat_buffer( su3_matrix *buf );

void hqp_output_corr(char cname[], char smtag[], int disp[], double corr);

/* new data type that contains information about the geometry */
typedef struct {
  int llat[4];
  int maxlen;
#ifndef ANISOTROPY
  int max_r2;
#else
  Real max_r2;
#endif
} hqp_geom;

int hqp_disp_rsq_ok ( int disp[], hqp_geom *hqp );
void hqp_free_geometry( hqp_geom *hqp ); 
hqp_geom *hqp_geometry( char myname[] );
void local_lattice_size (int *llat);

int loop_rt ( hqp_geom *hqp, 
              int mysize, 
              int mi, 
              int disp[], 
              double *wils_loop);

void make_loops ( int t, int mysize,
                  double *wils_loop, 
                  su3_matrix *s_link, 
                  su3_matrix *s_link_f, 
                  su3_matrix *t_link_f ) ;

void ax_gauge();
int setup();
int readin(int prompt);
void gball_simp(int tot_smear);
void smearing(void);
void w_loop1(int tot_smear);
void w_loop2(int tot_smear);
void hybrid_loop1(int tot_smear);

#define MAXSMTAG 4

#if (defined HYP_3D_SMEARING || defined HYP_4D_SMEARING) 
#define HYP_SMEARING
#define SMEARING
#endif

#if ( defined APE_1D2_SMEARING /* Smear one spatial and the temporal direction with a single direction of staples */\
   || defined APE_1D_SMEARING  /* Smear one spatial direction with a single direction of staples */\
   || defined APE_2D_SMEARING  /* Smear two spatial directions with the opposed direction of staples */\
   || defined APE_3D_SMEARING  /* Smear three spatial directions with the two opposed directions of staples */\
   || defined APE_4D_SMEARING) /* Smear all four directions with each three directions of staples */
#define APE_SMEARING
#define SMEARING
#if ( (defined APE_1D2_SMEARING && defined APE_1D_SMEARING ) \
   || (defined APE_1D2_SMEARING && defined APE_2D_SMEARING ) \
   || (defined APE_1D2_SMEARING && defined APE_3D_SMEARING ) \
   || (defined APE_1D2_SMEARING && defined APE_4D_SMEARING ) \
   || (defined APE_1D_SMEARING && defined APE_2D_SMEARING ) \
   || (defined APE_1D_SMEARING && defined APE_3D_SMEARING ) \
   || (defined APE_1D_SMEARING && defined APE_4D_SMEARING ) \
   || (defined APE_2D_SMEARING && defined APE_3D_SMEARING ) \
   || (defined APE_2D_SMEARING && defined APE_4D_SMEARING ) \
   || (defined APE_3D_SMEARING && defined APE_4D_SMEARING ) )
BOMB THE COMPILE
#endif
#endif
#if ( (defined APE_4D_SMEARING || defined HYP_4D_SMEARING || defined APE_1D2_SMEARING) && (defined AX_GAUGE) )
BOMB THE COMPILE
#endif
#if ( (defined APE_2D_SMEARING || defined APE_1D_SMEARING || defined APE_1D2_SMEARING ) )
#define STAP_APE
#endif

#if (defined SMEARING && defined GFIXONLY) 
BOMB THE COMPILE
#endif

#if (defined DIQUARK && ! defined COULOMB )
BOMB THE COMPILE
#endif

#if (!(defined COULOMB)&&!(defined PLOOPCOR_MEAS)&&!(defined PLANE_AVERAGED_PLC))
#define WLOOP_MEAS
#endif

#if (!(defined COULOMB)&&!(defined WLOOP_MEAS) &&(defined PLOOPCOR_MEAS||defined PLANE_AVERAGED_PLC))
#define PLOOPCOR_EXCLUSIVE
#endif

#ifdef NEW_HVY_POT
enum { HQPALG_AUTO=0, HQPALG_ONE=1, HQPALG_TWO=2, 
       HQPALG_THREE=3, HQPALG_FOUR=4, HQPALG_FIVE=5, 
       HQPALG_TIME=6 };
#define LLS_THRESH 512

void hqp_cycle_spatial( su3_matrix *links, int hqp_alg );
void hqp_switch( su3_matrix *links, int hqp_alg );
void hvy_pot_alg_old( su3_matrix *links );
void hvy_pot_alg_new( su3_matrix *links );
#endif

#ifdef PLANE_AVERAGED_PLC
void plane_averaged_plc( su3_matrix *links );
#endif

#ifdef OCTET_WLOOP
void free_gen(su3_matrix *genT);
su3_matrix *set_gen(void);
#endif

#endif
