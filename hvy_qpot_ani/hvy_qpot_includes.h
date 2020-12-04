/****************** hvy_qpot_includes.h ******************************/
/*
*  Include files for heavy quark potential application
*/

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

void ax_gauge();
int setup();
int readin(int prompt);
void gball_simp(int tot_smear);
void smearing(void);
void w_loop1(int tot_smear);
void w_loop2(int tot_smear);
void hybrid_loop1(int tot_smear);

#if (defined HYP_3D_SMEARING || defined HYP_4D_SMEARING) 
#define HYP_SMEARING
#define SMEARING
#endif

#if ( defined APE_1D2_SMEARING /* Smear one spatial and the temporal direction with a single direction of staples */\
   || defined APE_1D_SMEARING  /* Smear one spatial direction with a single direction of staples */\
   || defined APE_2D_SMEARING  /* Smear two spatial directions with the opposed direction of staples */\
   || defined APE_3D_SMEARING  /* Smear three spatial directions with the two opposed directions of staples */\
   || defined APE_4D_SMEARING) /* Smear all four directions with each three directions of staples */
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


#define site_coord(s,mu) \
        (((short*)&(s->x))[mu])

#ifdef NEW_HVY_POT
void local_lattice_size (int *llat);
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
#endif

#ifdef NEW_HVY_POT
void new_hvy_pot( su3_matrix *links );
#endif

#ifdef PLANE_AVERAGED_PLC
void plane_averaged_plc( su3_matrix *links );
#endif

#ifdef OCTET_WLOOP
void free_gen(su3_matrix *genT);
su3_matrix *set_gen(void);
#endif

