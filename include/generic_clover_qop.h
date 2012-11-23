#ifndef _GENERIC_CLOVER_QOP_H
#define _GENERIC_CLOVER_QOP_H
/******************** generic_clover_qop.h *********************************
*  MIMD version 7 		  				       *
*/

#include <qop.h>

/* d_bicgilu_cl_qop_F.c */

int 
bicgilu_cl_milc2qop_F( int prop_type, wilson_vector *milc_src, wilson_vector *milc_dest, 
		       quark_invert_control *qic, void *dmp);

void
bicgilu_cl_qop_generic_F( int prop_type, QOP_info_t *info,
			  QOP_F3_FermionLinksWilson *qop_links, 
			  QOP_invert_arg_t *qop_invert_arg,
			  QOP_resid_arg_t  ***qop_resid_arg,
			  void *dmps[], int nkappa[],
			  QOP_F3_DiracFermion **qop_sol[], 
			  QOP_F3_DiracFermion *qop_src[], 
			  int nsrc );
  
/* d_bicgilu_cl_qop_D.c */

int 
bicgilu_cl_milc2qop_D( int prop_type, wilson_vector *milc_src, wilson_vector *milc_dest, 
		       quark_invert_control *qic, void *dmp);

void
bicgilu_cl_qop_generic_D( int prop_type, QOP_info_t *info,
			  QOP_D3_FermionLinksWilson *qop_links, 
			  QOP_invert_arg_t *qop_invert_arg,
			  QOP_resid_arg_t  ***qop_resid_arg,
			  void *dmps[], int nkappa[],
			  QOP_D3_DiracFermion **qop_sol[], 
			  QOP_D3_DiracFermion *qop_src[], 
			  int nsrc );
  

#endif /* _GENERIC_CLOVER_QOP_H */
