/*************************************/
/* Function declarations of D-3 type */
/*************************************/

#ifndef _MGRID_D3_H
#define _MGRID_D3_H

#ifdef __cplusplus
extern "C" {
#endif

// The following strategem is intended to provide a strong chain of typedefs for both C and C++ compilation
// Above the API the data types are opaque (incomplete).
// If mGrid_internal.h is used, it must precede this header file

#if defined( _MGRID_INTERNAL_H) && defined( __cplusplus)
#define STRUCT
#else
#define STRUCT struct
#endif

#include "../include/milc_datatypes.h"


typedef STRUCT GRID_D3_ColorVector_struct		GRID_D3_ColorVector;
typedef STRUCT GRID_D3_ColorVectorBlock_struct		GRID_D3_ColorVectorBlock;
typedef STRUCT GRID_D3_ColorMatrix_struct	        GRID_D3_ColorMatrix;
typedef STRUCT GRID_D3_FermionLinksAsqtad_struct	GRID_D3_FermionLinksAsqtad;

// create color vectors
GRID_D3_ColorVector *GRID_D3_create_V( int milc_parity,
				       GRID_4Dgrid *grid_full,GRID_4DRBgrid *grid_rb );

GRID_D3_ColorVectorBlock *GRID_D3_create_nV( int n, int milc_parity,
					     GRID_5Dgrid *grid_5D, GRID_5DRBgrid *grid_5Drb, 
					     GRID_4Dgrid *grid_full,GRID_4DRBgrid *grid_rb );
					       
// free color vectors
void GRID_D3_destroy_V(GRID_D3_ColorVector *V);
void GRID_D3_destroy_nV(GRID_D3_ColorVectorBlock *V);

// create color vectors from MILC type
GRID_D3_ColorVector  *GRID_D3_create_V_from_vec( su3_vector *src, int milc_parity,
						 GRID_4Dgrid *grid_full,GRID_4DRBgrid *grid_rb );

GRID_D3_ColorVectorBlock  *GRID_D3_create_nV_from_vecs( su3_vector *src[], int n, int milc_parity,
							GRID_5Dgrid *grid_5D,GRID_5DRBgrid *grid_5Drb, 
							GRID_4Dgrid *grid_full, GRID_4DRBgrid *grid_rb );

// copy color vectors from Grid structure to MILC type
void GRID_D3_extract_V_to_vec( su3_vector *dest, GRID_D3_ColorVector *src, int milc_parity);
void GRID_D3_extract_nV_to_vecs( su3_vector *dest[], int n, GRID_D3_ColorVectorBlock *src, int milc_parity);

  /*********************/
  /*  FN routines  */
  /*********************/

  /* fermion matrix link routines */

// link fattening
void GRID_D3_fn_links(GRID_info_t *info,
		      GRID_D3_FermionLinksAsqtad *out,
		      su3_matrix *in,
		      GRID_4Dgrid *grid_full);
  
// create asqtad fermion links from MILC
GRID_D3_FermionLinksAsqtad  *GRID_D3_asqtad_create_L_from_MILC( su3_matrix *thn, su3_matrix *fat, 
								su3_matrix *lng, GRID_4Dgrid *grid_full);

// extract MILC links from the asqtad link structure
void GRID_D3_extract_MILC_from_L( su3_matrix *fat, su3_matrix *lng, GRID_D3_FermionLinksAsqtad  *fn,
				  GRID_4Dgrid *grid_full );

// free asqtad fermion links
void GRID_D3_asqtad_destroy_L(GRID_D3_FermionLinksAsqtad *L);

// dslash
void GRID_D3_asqtad_dslash (GRID_D3_FermionLinksAsqtad *asqtad,
			    GRID_D3_ColorVector *out,
			    GRID_D3_ColorVector *in,
			    double mass,
			    int milc_parity);

// inverter
void GRID_D3_asqtad_invert (GRID_info_t *info,
			    GRID_D3_FermionLinksAsqtad *asqtad,
			    GRID_invert_arg_t *inv_arg,
			    GRID_resid_arg_t *res_arg,
			    double mass,
			    GRID_D3_ColorVector *out,
			    GRID_D3_ColorVector *in,
			    GRID_4Dgrid *grid_full, GRID_4DRBgrid *grid_rb);

// multi-mass inverter
void GRID_D3_asqtad_invert_multi (GRID_info_t *info,
				  GRID_D3_FermionLinksAsqtad *asqtad,
				  GRID_invert_arg_t *inv_arg,
				  GRID_resid_arg_t *res_arg[],
				  double *mass, int nmass,
				  GRID_D3_ColorVector *out[],
				  GRID_D3_ColorVector *in,
				  GRID_4Dgrid *grid_full, GRID_4DRBgrid *grid_rb);

  // mixed-precision CG inverter
void GRID_D3_asqtad_invert_mixed(GRID_info_t *info,
				 GRID_F3_FermionLinksAsqtad *asqtad_f,
				 GRID_D3_FermionLinksAsqtad *asqtad_d,
				 GRID_invert_arg_t *inv_arg,
				 GRID_resid_arg_t *res_arg,
				 double mass,
				 GRID_D3_ColorVector *out,
				 GRID_D3_ColorVector *in,
				 GRID_4Dgrid *grid_full, GRID_4DRBgrid *grid_rb);

  // block CG inverter
void GRID_D3_asqtad_invert_block (GRID_info_t *info,
				  GRID_D3_FermionLinksAsqtad *asqtad,
				  GRID_invert_arg_t *inv_arg,
				  GRID_resid_arg_t *res_arg,
				  double mass, int nrhs,
				  GRID_D3_ColorVectorBlock *out,
				  GRID_D3_ColorVectorBlock *in,
				  GRID_5Dgrid *grid_5D, GRID_5DRBgrid *grid_5Drb, 
				  GRID_4Dgrid *grid_full, GRID_4DRBgrid *grid_rb);

  // block mixed-precision CG inverter
void GRID_D3_asqtad_invert_mixed_block(GRID_info_t *info,
				       GRID_F3_FermionLinksAsqtad *asqtad_f,
				       GRID_D3_FermionLinksAsqtad *asqtad_d,
				       GRID_invert_arg_t *inv_arg,
				       GRID_resid_arg_t *res_arg,
				       double mass, int nrhs,
				       GRID_D3_ColorVectorBlock *out,
				       GRID_D3_ColorVectorBlock *in,
				       GRID_5Dgrid *grid_5D, GRID_5DRBgrid *grid_5Drb, 
				       GRID_4Dgrid *grid_full, GRID_4DRBgrid *grid_rb);

  /*********************/
  /*  HISQ routines  */
  /*********************/
void GRID_D3_hisq_links(GRID_info_t *info,
			double path_coeff[],
			su3_matrix *fat,
			su3_matrix *lng,
			su3_matrix *in,
			GRID_4Dgrid *grid_full);

void GRID_D3_hisq_aux_links(GRID_info_t *info,
			    double path_coeff[],
			    su3_matrix *U, su3_matrix *V, su3_matrix *W,
			    GRID_4Dgrid *grid_full);

void GRID_D3_reunit_deriv(GRID_info_t *info, su3_matrix *V, su3_matrix *dW,
			  su3_matrix *Q, GRID_4Dgrid * grid_full);

/* implicitly restarted Lanczos */

typedef STRUCT GRID_D3_ColorVectorArray_struct GRID_D3_ColorVectorArray;

void GRID_D3_implicitly_restarted_lanczos(
  GRID_D3_ColorVectorArray * eigVecs,
  double * eigVals,
  GRID_D3_FermionLinksAsqtad * asqtad,
  GRID_eig_arg_t * eig_arg,
  double mass,
  GRID_4Dgrid * grid_full,
  GRID_4DRBgrid * grid_rb );
  
/* create color vector array */
GRID_D3_ColorVectorArray * GRID_D3_create_V_array(
  int n,
  int milc_parity,
  GRID_4Dgrid * grid_full,
  GRID_4DRBgrid * grid_rb );

/* free color vector array */
void GRID_D3_destroy_V_array( GRID_D3_ColorVectorArray * V );
  
/* create color vector array from MILC type */
GRID_D3_ColorVectorArray * GRID_D3_create_V_array_from_vec_array(
  su3_vector ** src,
  int n,
  int milc_parity,
  GRID_4Dgrid * grid_full,
  GRID_4DRBgrid * grid_rb );

/* copy color vector array from Grid structure to MILC type */
void GRID_D3_extract_V_array_to_vec_array(
  su3_vector ** dest,
  int n,
  GRID_D3_ColorVectorArray * src,
  int milc_parity );
  
  /**************************************************/
  /* Mapping of generic names to specific precision */
  /**************************************************/

#if GRID_Precision == 'D'
#include "mGrid_d3_generic.h"
#endif


#ifdef __cplusplus
}
#endif

#endif /* _MGRID_D3_H */
