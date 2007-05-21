#ifndef _GENERIC_CLOVER_QOP_H
#define _GENERIC_CLOVER_QOP_H
/******************** generic_clover_qop.h *********************************
*  MIMD version 7 		  				       *
*/

#include <qop.h>

/* d_bicgilu_cl_qop_F.c */

int 
bicgilu_cl_milc2qop_F( wilson_vector *milc_src, wilson_vector *milc_dest, 
		       quark_invert_control *qic, void *dmp);

/* d_bicgilu_cl_qop_D.c */

int 
bicgilu_cl_milc2qop_D( wilson_vector *milc_src, wilson_vector *milc_dest, 
		       quark_invert_control *qic, void *dmp);

#endif /* _GENERIC_CLOVER_QOP_H */
