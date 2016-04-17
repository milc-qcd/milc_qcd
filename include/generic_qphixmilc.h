#ifndef _GENERIC_QPHIXMILC_H
#define _GENERIC_QPHIXMILC_H
/******************** generic_qphixmilc.h *****************************
*  MIMD version 7 	 				            *
*/

#include "../include/su3.h"

/* map_milc_to_qphixmilc.c */

fsu3_matrix * create_qphix_raw4_F_G (void);
fsu3_matrix ** create_qphix_raw4_F_F (void);
fsu3_vector * create_qphix_raw_F_V(void);
fwilson_vector * create_qphix_raw_F_D(void);

void destroy_qphix_raw4_F_G (fsu3_matrix *raw);
void destroy_qphix_raw4_F_F (fsu3_matrix *raw[]);
void destroy_qphix_raw_F_V (fsu3_vector *raw);
void destroy_qphix_raw_F_D (fwilson_vector *raw);

QPHIX_F3_FermionLinksAsqtad *create_qphix_F_L_from_fields(su3_matrix *fat, 
   su3_matrix *lng, su3_matrix *fatback, su3_matrix *lngback, int milc_parity);

fsu3_matrix * create_qphix_raw4_F_G_from_site(field_offset src, int milc_parity);
fsu3_matrix ** create_qphix_raw4_F_F_from_site(field_offset src, int milc_parity);
fsu3_vector * create_qphix_raw_F_V_from_site(field_offset src, int milc_parity);
fwilson_vector * create_qphix_raw_F_D_from_site(field_offset src, int milc_parity);

fsu3_matrix * create_qphix_raw4_F_G_from_field(su3_matrix *src, int milc_parity);
fsu3_matrix ** create_qphix_raw4_F_F_from_field(anti_hermitmat *src, 
					  int milc_parity);
fsu3_vector * create_qphix_raw_F_V_from_field(su3_vector *src, int milc_parity);
fwilson_vector * create_qphix_raw_F_D_from_field(wilson_vector *src, int milc_parity);

void unload_qphix_raw4_F_G_to_site(field_offset dest, fsu3_matrix *raw, 
			     int milc_parity);
void unload_qphix_raw4_F_F_to_site(field_offset dest, fsu3_matrix *raw[], 
			     int milc_parity);
void unload_qphix_raw_F_V_to_site(field_offset dest, fsu3_vector *raw, 
			    int milc_parity);
void unload_qphix_raw_F_D_to_site(field_offset dest, fwilson_vector *raw, 
			    int milc_parity);

void unload_qphix_raw4_F_G_to_field(su3_matrix *dest, fsu3_matrix *raw, 
			      int milc_parity);
void unload_qphix_raw4_F_F_to_field(anti_hermitmat *dest, fsu3_matrix *raw, 
			      int milc_parity);
void unload_qphix_raw_F_V_to_field(su3_vector *dest, fsu3_vector *raw, 
			     int milc_parity);
void unload_qphix_raw_F_D_to_field(wilson_vector *dest, fwilson_vector *raw, 
			     int milc_parity);

dsu3_matrix * create_qphix_raw4_D_G (void);
dsu3_matrix ** create_qphix_raw4_D_F (void);
dsu3_vector * create_qphix_raw_D_V(void);
dwilson_vector * create_qphix_raw_D_D(void);

void destroy_qphix_raw4_D_G (dsu3_matrix *raw);
void destroy_qphix_raw4_D_F (dsu3_matrix *raw);
void destroy_qphix_raw_D_V (dsu3_vector *raw);
void destroy_qphix_raw_D_D (dwilson_vector *raw);

QPHIX_D3_FermionLinksAsqtad *create_qphix_D_L_from_fields(su3_matrix *fat, 
   su3_matrix *lng, su3_matrix *fatback, su3_matrix *lngback, int milc_parity);

dsu3_matrix * create_qphix_raw4_D_G_from_site(field_offset src, int milc_parity);
dsu3_vector * create_qphix_raw_D_V_from_site(field_offset src, int milc_parity);
dwilson_vector * create_qphix_raw_D_D_from_site(field_offset src, int milc_parity);

dsu3_matrix * create_qphix_raw4_D_G_from_field(su3_matrix *src, int milc_parity);
dsu3_matrix ** create_qphix_raw4_D_F_from_site(field_offset src, int milc_parity);
dsu3_vector * create_qphix_raw_D_V_from_site(field_offset src, int milc_parity);
dwilson_vector * create_qphix_raw_D_D_from_site(field_offset src, int milc_parity);

dsu3_matrix * create_qphix_raw4_D_G_from_field(su3_matrix *src, int milc_parity);
dsu3_matrix ** create_qphix_raw4_D_F_from_field(anti_hermitmat *src, 
					  int milc_parity);
dsu3_vector * create_qphix_raw_D_V_from_field(su3_vector *src, int milc_parity);
dwilson_vector * create_qphix_raw_D_D_from_field(wilson_vector *src, int milc_parity);

void unload_qphix_raw4_D_G_to_site(field_offset dest, dsu3_matrix *raw, 
			     int milc_parity);
void unload_qphix_raw4_D_G_to_field(su3_matrix *dest, dsu3_matrix *raw, 
			      int milc_parity);
void unload_qphix_raw4_D_F_to_site(field_offset dest, dsu3_matrix *raw[], 
			     int milc_parity);
void unload_qphix_raw4_D_F_to_field(anti_hermitmat *dest, dsu3_matrix *raw[], 
			      int milc_parity);
void unload_qphix_raw_D_V_to_site(field_offset dest, dsu3_vector *raw, 
			    int milc_parity);
void unload_qphix_raw_D_V_to_field(su3_vector *dest, dsu3_vector *raw, 
			     int milc_parity);
void unload_qphix_raw_D_D_to_site(field_offset dest, dwilson_vector *raw, 
			    int milc_parity);
void unload_qphix_raw_D_D_to_field(wilson_vector *dest, dwilson_vector *raw, 
			     int milc_parity);

#endif /* GENERIC_QPHIXMILC_H */
