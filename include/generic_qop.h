#ifndef _GENERIC_QOP_H
#define _GENERIC_QOP_H
/******************** generic_qop.h *********************************
*  MIMD version 7 	 				            *
*/

#include <qop.h>
#include "../include/su3.h"

/* milc_to_qop_utilities.c */
QOP_evenodd_t milc2qop_parity(int milc_parity);
int qop2milc_parity(QOP_evenodd_t qop_parity);
QOP_status_t initialize_qop(void);

/* map_milc_to_qop*.c */

/* Generic mapping */

#if ( QOP_Precision == 1 )

#define create_raw4_G  create_raw4_F_G
#define create_raw4_F  create_raw4_F_F
#define create_raw_V   create_raw_F_V
#define create_raw_D   create_raw_F_D

#define destroy_raw4_G  destroy_raw4_F_G
#define destroy_raw4_F  destroy_raw4_F_F
#define destroy_raw_V   destroy_raw_F_V
#define destroy_raw_D   destroy_raw_F_D

#define create_raw4_G_from_site create_raw4_F_G_from_site
#define create_raw4_F_from_site create_raw4_F_F_from_site
#define create_raw_V_from_site  create_raw_F_V_from_site
#define create_raw_D_from_site  create_raw_F_D_from_site

#define create_raw4_G_from_field create_raw4_F_G_from_field
#define create_raw4_F_from_field create_raw4_F_F_from_field
#define create_raw_V_from_field  create_raw_F_V_from_field
#define create_raw_D_from_field  create_raw_F_D_from_field

#define create_F_from_site4  create_F_F_from_site4
#define create_G_from_site4  create_F_G_from_site4
#define create_V_from_site   create_F_V_from_site
#define create_D_from_site   create_F_D_from_site

#define create_F_from_field4  create_F_F_from_field4
#define create_G_from_field4  create_F_G_from_field4
#define create_V_from_field   create_F_V_from_field
#define create_D_from_field   create_F_D_from_field

#define unload_raw4_G_to_site unload_raw4_F_G_to_site
#define unload_raw4_F_to_site unload_raw4_F_F_to_site
#define unload_raw_V_to_site  unload_raw_F_V_to_site
#define unload_raw_D_to_site  unload_raw_F_D_to_site

#define unload_raw4_G_to_field unload_raw4_F_G_to_field
#define unload_raw4_F_to_field unload_raw4_F_F_to_field
#define unload_raw_V_to_field  unload_raw_F_V_to_field
#define unload_raw_D_to_field  unload_raw_F_D_to_field

#define unload_F_to_site4  unload_F_F_to_site4
#define unload_G_to_site4  unload_F_G_to_site4

#define unload_F_to_field4  unload_F_F_to_field4
#define unload_G_to_field4  unload_F_G_to_field4

#define unload_V_to_site  unload_F_V_to_site
#define unload_V_to_field unload_F_V_to_field
#define unload_D_to_site  unload_F_D_to_site
#define unload_D_to_field unload_F_D_to_field

#define map_milc_clov_to_qop_raw map_milc_clov_to_qop_raw_F

#else

#define create_raw4_G  create_raw4_D_G
#define create_raw4_F  create_raw4_D_F
#define create_raw_V   create_raw_D_V
#define create_raw_D   create_raw_D_D

#define destroy_raw4_G  destroy_raw4_D_G
#define destroy_raw4_F  destroy_raw4_D_F
#define destroy_raw_V   destroy_raw_D_V
#define destroy_raw_D   destroy_raw_D_D

#define create_raw4_G_from_site create_raw4_D_G_from_site
#define create_raw4_F_from_site create_raw4_D_F_from_site
#define create_raw_V_from_site  create_raw_D_V_from_site
#define create_raw_D_from_site  create_raw_D_D_from_site

#define create_raw4_G_from_field create_raw4_D_G_from_field
#define create_raw4_F_from_field create_raw4_D_F_from_field
#define create_raw_V_from_field  create_raw_D_V_from_field
#define create_raw_D_from_field  create_raw_D_D_from_field

#define create_F_from_site4  create_D_F_from_site4
#define create_G_from_site4  create_D_G_from_site4
#define create_V_from_site   create_D_V_from_site
#define create_D_from_site   create_D_D_from_site

#define create_F_from_field4  create_D_F_from_field4
#define create_G_from_field4  create_D_G_from_field4
#define create_V_from_field   create_D_V_from_field
#define create_D_from_field   create_D_D_from_field

#define unload_raw4_G_to_site unload_raw4_D_G_to_site
#define unload_raw4_F_to_site unload_raw4_D_F_to_site
#define unload_raw_V_to_site  unload_raw_D_V_to_site
#define unload_raw_D_to_site  unload_raw_D_D_to_site

#define unload_raw4_G_to_field unload_raw4_D_G_to_field
#define unload_raw4_F_to_field unload_raw4_D_F_to_field
#define unload_raw_V_to_field  unload_raw_D_V_to_field
#define unload_raw_D_to_field  unload_raw_D_D_to_field

#define unload_F_to_site4  unload_D_F_to_site4
#define unload_G_to_site4  unload_D_G_to_site4

#define unload_F_to_field4  unload_D_F_to_field4
#define unload_G_to_field4  unload_D_G_to_field4

#define unload_V_to_site  unload_D_V_to_site
#define unload_V_to_field unload_D_V_to_field
#define unload_D_to_site  unload_D_D_to_site
#define unload_D_to_field unload_D_D_to_field

#define map_milc_clov_to_qop_raw map_milc_clov_to_qop_raw_D

#endif

QOP_F3_ColorVector* create_F_V_from_site(field_offset src, int parity);
QOP_F3_ColorVector* create_F_V_from_field(su3_vector *src, int parity);
QOP_F3_DiracFermion* create_F_D_from_site(field_offset src,	int parity);
QOP_F3_DiracFermion* create_F_D_from_field(wilson_vector *src, 
					       int parity);
QOP_F3_GaugeField *create_F_G_from_site4(field_offset link, int parity);
QOP_F3_Force *create_F_F_from_site4(field_offset mom, int parity);

void unload_F_V_to_site(field_offset dest, QOP_F3_ColorVector* src, 
			int parity);
void unload_F_V_to_field(su3_vector *dest, QOP_F3_ColorVector* src, 
			int parity);
void unload_F_D_to_site(field_offset dest, QOP_F3_DiracFermion* src, 
			int parity);
void unload_F_D_to_field(wilson_vector *dest, QOP_F3_DiracFermion* src, 
			int parity);
void unload_F_F_to_site4(field_offset link, QOP_F3_Force *qop, 
			 int parity);
void unload_F_G_to_site4(field_offset link, QOP_F3_GaugeField *qop, 
			 int parity);

QOP_F3_FermionLinksAsqtad* create_F_L_from_sites( field_offset fat, 
			 field_offset lng, int parity);

QOP_F3_FermionLinksAsqtad* create_F_L_from_fields( su3_matrix *fat, 
			  su3_matrix *lng, int parity);

void unload_F_L_to_fields( su3_matrix *fat, su3_matrix *lng, 
			  QOP_F3_FermionLinksAsqtad* qop, int parity);

void unload_F_hisq_L_to_fields( su3_matrix *fat, su3_matrix *lng, 
			  QOP_F3_FermionLinksHisq* qop, int parity);


QOP_D3_ColorVector* create_D_V_from_site(field_offset src, int parity);
QOP_D3_ColorVector* create_D_V_from_field(su3_vector *src, int parity);
QOP_D3_DiracFermion* create_D_D_from_site(field_offset src,	int parity);
QOP_D3_DiracFermion* create_D_D_from_field(wilson_vector *src, 
					       int parity);
QOP_D3_GaugeField *create_D_G_from_site4(field_offset link, int parity);
QOP_D3_Force *create_D_F_from_site4(field_offset mom, int parity);

void unload_D_V_to_site(field_offset dest, QOP_D3_ColorVector* src, 
			int parity);
void unload_D_V_to_field(su3_vector *dest, QOP_D3_ColorVector* src, 
			int parity);

void unload_D_D_to_site(field_offset dest, QOP_D3_DiracFermion* src, 
			int parity);
void unload_D_D_to_field(wilson_vector *dest, QOP_D3_DiracFermion* src, 
			int parity);
void unload_D_F_to_site4(field_offset link, QOP_D3_Force *qop, 
			 int parity);
void unload_D_G_to_site4(field_offset link, QOP_D3_GaugeField *qop,
			 int parity);

QOP_D3_FermionLinksAsqtad* create_D_L_from_sites( field_offset fat, 
			 field_offset lng, int parity);

QOP_D3_FermionLinksAsqtad* create_D_L_from_fields( su3_matrix *fat, 
			  su3_matrix *lng, int parity);

void unload_D_L_to_fields( su3_matrix *fat, su3_matrix *lng, 
			  QOP_D3_FermionLinksAsqtad* qop, int parity);

void unload_D_hisq_L_to_fields( su3_matrix *fat, su3_matrix *lng, 
			  QOP_D3_FermionLinksHisq* qop, int parity);

#include "../include/generic_clover.h"

void map_milc_clov_to_qop_raw_F(float *raw_clov, clover *milc_clov );
void map_milc_clov_to_qop_raw_D(double *raw_clov, clover *milc_clov );

#endif /* GENERIC_QOP_H */
