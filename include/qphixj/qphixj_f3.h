/*************************************/
/* Function declarations of F-3 type */
/*************************************/

#ifndef QPHIXJ_F3
#define QPHIXJ_F3

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

typedef STRUCT QPHIXJ_F3_DiracFermion_struct   QPHIXJ_F3_DiracFermion;
typedef STRUCT QPHIXJ_F3_FermionLinksWilson_struct  QPHIXJ_F3_FermionLinksWilson;

// create Dirac vectors from MILC to qphixj structure
QPHIXJ_F3_DiracFermion  *QPHIXJ_F3_create_D_from_wvec( wilson_vector *src, QPHIXJ_evenodd_t evenodd );

// copy Dirac vectors from qphixj structure to MILC
void QPHIXJ_F3_extract_D_to_wvec( wilson_vector *dest, QPHIXJ_F3_DiracFermion  *src, QPHIXJ_evenodd_t evenodd );

// create wilson fermion links from MILC
QPHIXJ_F3_FermionLinksWilson  *QPHIXJ_F3_wilson_create_L_from_MILC( su3_matrix *backlinks,
								    clover *clov, Real kappa,
								    QPHIXJ_evenodd_t evenodd );
// free Dirac vectors
void  QPHIXJ_F3_destroy_D( QPHIXJ_F3_DiracFermion *D );

// free wilson fermion links
void  QPHIXJ_F3_wilson_destroy_L( QPHIXJ_F3_FermionLinksWilson *L );

// dslash
void QPHIXJ_F3_wilson_dslash( QPHIXJ_F3_FermionLinksWilson *wilson,
			      QPHIXJ_F3_DiracFermion *out,
			      QPHIXJ_F3_DiracFermion *in,
			      QPHIXJ_evenodd_t parity );

// inverter
void QPHIXJ_F3_wilson_invert(QPHIXJ_info_t *info,
			     QPHIXJ_F3_FermionLinksWilson *links,
			     QPHIXJ_invert_arg_t *inv_arg,
			     QPHIXJ_resid_arg_t *res_arg,
			     QPHIXJ_F_Real kappa,
			     QPHIXJ_F3_DiracFermion *out_pt,
			     QPHIXJ_F3_DiracFermion *in_pt);
  
// multi-mass inverter
void QPHIXJ_F3_wilson_invert_multi(QPHIXJ_info_t *info,
				   QPHIXJ_F3_FermionLinksWilson *links,
				   QPHIXJ_invert_arg_t *inv_arg,
				   QPHIXJ_resid_arg_t **res_arg[],
				   QPHIXJ_F_Real *kappas[],
				   int nkappa[],
				   QPHIXJ_F3_DiracFermion **out_pt[],
				   QPHIXJ_F3_DiracFermion *in_pt[],
				   int nsrc);

  /**************************************************/
  /* Mapping of generic names to specific precision */
  /**************************************************/

#if QPHIXJ_Precision == 'F'
#include "../include/qphixj/qphixj_f3_generic.h"
#endif

#ifdef __cplusplus
};
#endif

#endif // QPHIXJ_F3
