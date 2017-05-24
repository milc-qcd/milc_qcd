/* QPhiX-MILC Interface */
#ifndef _GENERIC_QPHIXJ_H
#define _GENERIC_QPHIXJ_H

#include <../include/qphixj/qphixj.h>
#include "../include/config.h"
#include "../include/generic_ks.h"
#include "../include/generic.h"
#include "../include/su3.h"
#include "../include/dirs.h"
#include "../include/fermion_links.h"
#include <stdbool.h>

/* milc_to_qphixj_utilities.c */
QPHIXJ_evenodd_t milc2qphixj_parity(int milc_parity);
int qphixj2milc_parity(QPHIXJ_evenodd_t qphixj_parity);
QPHIXJ_status_t initialize_qphixj(void);

/* d_bicgilu_qphixj_*.c */

#include "../include/generic_clover.h"

void bicgilu_cl_qphixj_inner_F(clover *milc_clov, Real kappa, wilson_vector r[], 
			       wilson_vector dest[], quark_invert_control *qic);

void bicgilu_cl_qphixj_inner_D(clover *milc_clov, Real kappa, wilson_vector r[], 
			       wilson_vector dest[], quark_invert_control *qic);

/* Generic mapping */

#if ( QPHIXJ_PrecisionInt == 1 )

#define create_qphixj_raw4_G  create_qphixj_raw4_F_G
#define create_qphixj_raw4_F  create_qphixj_raw4_F_F
#define create_qphixj_raw_V   create_qphixj_raw_F_V
#define create_qphixj_raw_D   create_qphixj_raw_F_D

#define destroy_qphixj_raw4_G  destroy_qphixj_raw4_F_G
#define destroy_qphixj_raw4_F  destroy_qphixj_raw4_F_F
#define destroy_qphixj_raw_V   destroy_qphixj_raw_F_V
#define destroy_qphixj_raw_D   destroy_qphixj_raw_F_D

#define create_qphixj_L_from_sites  create_qphixj_F_L_from_sites
#define create_qphixj_L_from_fieldback_and_sites  create_qphixj_F_L_from_fieldback_and_sites
#define create_qphixj_L_from_fields    create_qphixj_F_L_from_fields

#define create_qphixj_raw4_G_from_site create_qphixj_raw4_F_G_from_site

#define create_qphixj_raw4_G_from_field create_qphixj_raw4_F_G_from_field
#define create_qphixj_raw4_F_from_field create_qphixj_raw4_F_F_from_field
#define create_qphixj_raw_V_from_field  create_qphixj_raw_F_V_from_field
#define create_qphixj_raw_D_from_field  create_qphixj_raw_F_D_from_field

#define create_qphixj_F_from_field4  create_qphixj_F_F_from_field4
#define create_qphixj_G_from_field4  create_qphixj_F_G_from_field4
#define create_qphixj_V_from_field   create_qphixj_F_V_from_field
#define create_qphixj_D_from_field   create_qphixj_F_D_from_field

#define unload_qphixj_raw4_G_to_field unload_qphixj_raw4_F_G_to_field
#define unload_qphixj_raw4_F_to_field unload_qphixj_raw4_F_F_to_field
#define unload_qphixj_raw_V_to_field  unload_qphixj_raw_F_V_to_field
#define unload_qphixj_raw_D_to_field  unload_qphixj_raw_F_D_to_field

#define unload_qphixj_F_to_field4  unload_qphixj_F_F_to_field4
#define unload_qphixj_G_to_field4  unload_qphixj_F_G_to_field4

#define unload_qphixj_V_to_field unload_qphixj_F_V_to_field
#define unload_qphixj_D_to_field unload_qphixj_F_D_to_field

#define map_milc_clov_to_qphixj_raw map_milc_clov_to_qphixj_raw_F

#else

#define create_qphixj_raw4_G  create_qphixj_raw4_D_G
#define create_qphixj_raw4_F  create_qphixj_raw4_D_F
#define create_qphixj_raw_V   create_qphixj_raw_D_V
#define create_qphixj_raw_D   create_qphixj_raw_D_D

#define destroy_qphixj_raw4_G  destroy_qphixj_raw4_D_G
#define destroy_qphixj_raw4_F  destroy_qphixj_raw4_D_F
#define destroy_qphixj_raw_V   destroy_qphixj_raw_D_V
#define destroy_qphixj_raw_D   destroy_qphixj_raw_D_D

#define create_qphixj_L_from_sites  create_qphixj_D_L_from_sites
#define create_qphixj_L_from_fieldback_and_sites  create_qphixj_D_L_from_fieldback_and_sites
#define create_qphixj_L_from_fields    create_qphixj_D_L_from_fields

#define create_qphixj_raw4_G_from_site create_qphixj_raw4_D_G_from_site

#define create_qphixj_raw4_G_from_field create_qphixj_raw4_D_G_from_field
#define create_qphixj_raw4_F_from_field create_qphixj_raw4_D_F_from_field
#define create_qphixj_raw_V_from_field  create_qphixj_raw_D_V_from_field
#define create_qphixj_raw_D_from_field  create_qphixj_raw_D_D_from_field

#define create_qphixj_F_from_field4   create_qphixj_D_F_from_field4
#define create_qphixj_G_from_field4   create_qphixj_D_G_from_field4
#define create_qphixj_V_from_field   create_qphixj_D_V_from_field
#define create_qphixj_D_from_field   create_qphixj_D_D_from_field

#define unload_qphixj_raw4_G_to_field unload_qphixj_raw4_D_G_to_field
#define unload_qphixj_raw4_F_to_field unload_qphixj_raw4_D_F_to_field
#define unload_qphixj_raw_V_to_field  unload_qphixj_raw_D_V_to_field
#define unload_qphixj_raw_D_to_field  unload_qphixj_raw_D_D_to_field


#define unload_qphixj_F_to_field4  unload_qphixj_D_F_to_field4
#define unload_qphixj_G_to_field4  unload_qphixj_D_G_to_field4

#define unload_qphixj_V_to_field unload_qphixj_D_V_to_field
#define unload_qphixj_D_to_field unload_qphixj_D_D_to_field

#define map_milc_clov_to_qphixj_raw map_milc_clov_to_qphixj_raw_D

#endif

fsu3_matrix *create_qphixj_raw4_F_G_from_site(field_offset src, int parity);

QPHIXJ_F3_DiracFermion *create_qphixj_F_D_from_field(wilson_vector *src, int parity);

fwilson_vector *create_qphixj_raw_F_D_from_field(wilson_vector *src, int milc_parity);

void unload_qphixj_F_D_to_field(wilson_vector *dest, QPHIXJ_F3_DiracFermion* src, 
				int parity);

void destroy_qphixj_raw4_F_G(fsu3_matrix *raw);

void destroy_qphixj_raw_F_D(fwilson_vector *raw);

QPHIXJ_F3_FermionLinksWilson *
create_qphixj_F_L_from_sites (float *clov, int parity);

QPHIXJ_F3_FermionLinksWilson *
create_qphixj_F_L_from_fieldback_and_sites(float *clov, 
					   fsu3_matrix *fieldback, int parity);

dsu3_matrix *create_qphixj_raw4_D_G_from_site(field_offset src, int parity);

QPHIXJ_D3_DiracFermion *create_qphixj_D_D_from_field(wilson_vector *src, int parity);

dwilson_vector *create_qphixj_raw_D_D_from_field(wilson_vector *src, int milc_parity);

void unload_qphixj_D_D_to_field(wilson_vector *dest, QPHIXJ_D3_DiracFermion* src, 
				int parity);

void destroy_qphixj_raw4_D_G(dsu3_matrix *raw);

void destroy_qphixj_raw_D_D(dwilson_vector *raw);

QPHIXJ_D3_FermionLinksWilson *
create_qphixj_D_L_from_sites (float *clov, int parity);

QPHIXJ_D3_FermionLinksWilson *
create_qphixj_D_L_from_fieldback_and_sites(double *clov, 
					   dsu3_matrix *fieldback, int parity);

#endif // _GENERIC_QPHIXJ_H
