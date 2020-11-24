/*************************************/
/* Function declarations of D-3 type */
/*************************************/

#ifndef QPHIXJ_D3
#define QPHIXJ_D3

#ifdef __cplusplus
extern "C" {
#endif

// This strategem is intended to provide a strong chain of typedefs for both C and C++ compilation 
// Above the API the data types are opaque (incomplete).
// If qphixj_internal.h is used, it must precede this header file
#if defined( QPHIXJ_INTERNAL_H) && defined( __cplusplus)
#define STRUCT
#else
#define STRUCT struct
#endif

#include "../include/generic_clover.h"

typedef STRUCT QPHIXJ_D3_DiracFermion_struct   QPHIXJ_D3_DiracFermion;
typedef STRUCT QPHIXJ_D3_FermionLinksWilson_struct  QPHIXJ_D3_FermionLinksWilson;

// create Dirac vectors from MILC to qphixj structure
QPHIXJ_D3_DiracFermion  *QPHIXJ_D3_create_D_from_wvec( wilson_vector *src, QPHIXJ_evenodd_t evenodd );

// copy Dirac vectors from qphixj structure to MILC
void QPHIXJ_D3_extract_D_to_wvec( wilson_vector *dest, QPHIXJ_D3_DiracFermion  *src, QPHIXJ_evenodd_t evenodd );

// create wilson fermion links from MILC
QPHIXJ_D3_FermionLinksWilson  *QPHIXJ_D3_wilson_create_L_from_MILC( su3_matrix *backlinks,
								    clover *clov, Real kappa,
								    QPHIXJ_evenodd_t evenodd );

// free Dirac vectors
void  QPHIXJ_D3_destroy_D( QPHIXJ_D3_DiracFermion *D );

// free wilson fermion links
void  QPHIXJ_D3_wilson_destroy_L( QPHIXJ_D3_FermionLinksWilson *L );

// dslash
void QPHIXJ_D3_wilson_dslash( QPHIXJ_D3_FermionLinksWilson *wilson,
			      QPHIXJ_D3_DiracFermion *out,
			      QPHIXJ_D3_DiracFermion *in,
			      QPHIXJ_evenodd_t parity );

// inverter
void QPHIXJ_D3_wilson_invert(QPHIXJ_info_t *info,
			     QPHIXJ_D3_FermionLinksWilson *links,
			     QPHIXJ_invert_arg_t *inv_arg,
			     QPHIXJ_resid_arg_t *res_arg,
			     QPHIXJ_D_Real kappa,
			     QPHIXJ_D3_DiracFermion *out_pt,
			     QPHIXJ_D3_DiracFermion *in_pt);
  
// multi-mass inverter
void QPHIXJ_D3_wilson_invert_multi(QPHIXJ_info_t *info,
				   QPHIXJ_D3_FermionLinksWilson *links,
				   QPHIXJ_invert_arg_t *inv_arg,
				   QPHIXJ_resid_arg_t **res_arg[],
				   QPHIXJ_D_Real *kappas[],
				   int nkappa[],
				   QPHIXJ_D3_DiracFermion **out_pt[],
				   QPHIXJ_D3_DiracFermion *in_pt[],
				   int nsrc);

  /**************************************************/
  /* Mapping of generic names to specific precision */
  /**************************************************/

#if QPHIXJ_Precision == 'D'
#include "../include/qphixj/qphixj_d3_generic.h"
#endif

#ifdef __cplusplus
};
#endif

#endif // QPHIXJ_D3
