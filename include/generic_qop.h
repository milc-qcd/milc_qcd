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
QOP_status_t initialize_qop();

/* map_milc_to_qop.c */

/* Generic mapping */

#if ( QOP_Precision == 1 )

#define create_raw4_G  create_raw4_F_G
#define create_raw4_F  create_raw4_F_F
#define create_raw_V   create_raw_F_V

#define destroy_raw4_G  destroy_raw4_F_G
#define destroy_raw4_F  destroy_raw4_F_F
#define destroy_raw_V   destroy_raw_F_V

#define create_raw4_G_from_site create_raw4_F_G_from_site
#define create_raw4_F_from_site create_raw4_F_F_from_site
#define create_raw_V_from_site  create_raw_F_V_from_site

#define create_raw4_G_from_field create_raw4_F_G_from_field
#define create_raw4_F_from_field create_raw4_F_F_from_field
#define create_raw_V_from_field  create_raw_F_V_from_field

#define unload_raw4_G_to_site unload_raw4_F_G_to_site
#define unload_raw4_F_to_site unload_raw4_F_F_to_site
#define unload_raw_V_to_site  unload_raw_F_V_to_site

#define unload_raw4_G_to_field unload_raw4_F_G_to_field
#define unload_raw4_F_to_field unload_raw4_F_F_to_field
#define unload_raw_V_to_field  unload_raw_F_V_to_field

#define load_links_and_mom_site load_F_links_and_mom_site
#define load_links_and_mom_field load_F_links_and_mom_field

#define unload_links_and_mom_site unload_F_links_and_mom_site
#define unload_links_and_mom_field unload_F_links_and_mom_field

#define load_V_from_site  load_F_V_from_site
#define load_V_from_field load_F_V_from_field
#define unload_V_to_site  unload_F_V_to_site
#define unload_V_to_field unload_F_V_to_field

#else

#define create_raw4_G  create_raw4_D_G
#define create_raw4_F  create_raw4_D_F
#define create_raw_V   create_raw_D_V

#define destroy_raw4_G  destroy_raw4_D_G
#define destroy_raw4_F  destroy_raw4_D_F
#define destroy_raw_V   destroy_raw_D_V

#define create_raw4_G_from_site create_raw4_D_G_from_site
#define create_raw4_F_from_site create_raw4_D_F_from_site
#define create_raw_V_from_site  create_raw_D_V_from_site

#define create_raw4_G_from_field create_raw4_D_G_from_field
#define create_raw4_F_from_field create_raw4_D_F_from_field
#define create_raw_V_from_field  create_raw_D_V_from_field

#define unload_raw4_G_to_site unload_raw4_D_G_to_site
#define unload_raw4_F_to_site unload_raw4_D_F_to_site
#define unload_raw_V_to_site  unload_raw_D_V_to_site

#define unload_raw4_G_to_field unload_raw4_D_G_to_field
#define unload_raw4_F_to_field unload_raw4_D_F_to_field
#define unload_raw_V_to_field  unload_raw_D_V_to_field

#define load_links_and_mom_site load_D_links_and_mom_site
#define load_links_and_mom_field load_D_links_and_mom_field

#define unload_links_and_mom_site unload_D_links_and_mom_site
#define unload_links_and_mom_field unload_D_links_and_mom_field

#define load_V_from_site  load_D_V_from_site
#define load_V_from_field load_D_V_from_field
#define unload_V_to_site  unload_D_V_to_site
#define unload_V_to_field unload_D_V_to_field

#endif

fsu3_matrix ** create_raw4_F_G (void);
fsu3_matrix ** create_raw4_F_F (void);
fsu3_vector * create_raw_F_V(void);

void destroy_raw4_F_G (fsu3_matrix *raw[]);
void destroy_raw4_F_F (fsu3_matrix *raw[]);
void destroy_raw_F_V (fsu3_vector *raw);

fsu3_matrix ** create_raw4_F_G_from_site(field_offset src, int milc_parity);
fsu3_matrix ** create_raw4_F_F_from_site(field_offset src, int milc_parity);
fsu3_vector * create_raw_F_V_from_site(field_offset src, int milc_parity);

fsu3_matrix ** create_raw4_F_G_from_field(su3_matrix *src, int milc_parity);
fsu3_matrix ** create_raw4_F_F_from_field(anti_hermitmat *src, 
					  int milc_parity);
fsu3_vector * create_raw_F_V_from_field(su3_vector *src, int milc_parity);

void unload_raw4_F_G_to_site(field_offset dest, fsu3_matrix *raw[], 
			     int milc_parity);
void unload_raw4_F_F_to_site(field_offset dest, fsu3_matrix *raw[], 
			     int milc_parity);
void unload_raw_F_V_to_site(field_offset dest, fsu3_vector *raw, 
			    int milc_parity);

void unload_raw4_F_G_to_field(su3_matrix *dest, fsu3_matrix *raw[], 
			      int milc_parity);
void unload_raw4_F_F_to_field(anti_hermitmat *dest, fsu3_matrix *raw[], 
			      int milc_parity);
void unload_raw_F_V_to_field(su3_vector *dest, fsu3_vector *raw, 
			     int milc_parity);

void load_F_links_and_mom_site(QOP_F3_GaugeField **links, QOP_F3_Force **mom, 
			       fsu3_matrix ***rawlinks, fsu3_matrix ***rawmom);
void load_F_links_and_mom_field(QOP_F3_GaugeField **links, QOP_F3_Force **mom,
				fsu3_matrix ***rawlinks, fsu3_matrix ***rawmom,
				su3_matrix *srclink, 
				anti_hermitmat *srcmom);

void unload_F_links_and_mom_site(QOP_F3_GaugeField **links, QOP_F3_Force **mom,
				 fsu3_matrix ***rawlinks, 
				 fsu3_matrix ***rawmom);
void unload_F_links_and_mom_field(su3_matrix *dstlink, anti_hermitmat *dstmom,
	  QOP_F3_GaugeField **links, QOP_F3_Force **mom, 
          fsu3_matrix ***rawlinks, fsu3_matrix ***rawmom);

void load_F_V_from_site(QOP_F3_ColorVector** dest, field_offset src, 
			int parity);
void load_F_V_from_field(QOP_F3_ColorVector** dest, su3_vector *src, 
			int parity);
void unload_F_V_to_site(field_offset dest, QOP_F3_ColorVector* src, 
			int parity);
void unload_F_V_to_field(su3_vector *dest, QOP_F3_ColorVector* src, 
			int parity);


dsu3_matrix ** create_raw4_D_G (void);
void destroy_raw4_D_G (dsu3_matrix *raw[]);
dsu3_matrix ** create_raw4_D_F (void);
void destroy_raw4_D_F (dsu3_matrix *raw[]);
dsu3_vector * create_raw_D_V(void);
void destroy_raw_D_V (dsu3_vector *raw);
dsu3_matrix ** create_raw4_D_G_from_site(field_offset src, int milc_parity);
dsu3_matrix ** create_raw4_D_G_from_field(su3_matrix *src, int milc_parity);
dsu3_matrix ** create_raw4_D_F_from_site(field_offset src, int milc_parity);
dsu3_matrix ** create_raw4_D_F_from_field(anti_hermitmat *src, 
					  int milc_parity);
dsu3_vector * create_raw_D_V_from_site(field_offset src, int milc_parity);
dsu3_vector * create_raw_D_V_from_field(su3_vector *src, int milc_parity);
void unload_raw4_D_G_to_site(field_offset dest, dsu3_matrix *raw[], 
			     int milc_parity);
void unload_raw4_D_G_to_field(su3_matrix *dest, dsu3_matrix *raw[], 
			      int milc_parity);
void unload_raw4_D_F_to_site(field_offset dest, dsu3_matrix *raw[], 
			     int milc_parity);
void unload_raw4_D_F_to_field(anti_hermitmat *dest, dsu3_matrix *raw[], 
			      int milc_parity);
void unload_raw_D_V_to_site(field_offset dest, dsu3_vector *raw, 
			    int milc_parity);
void unload_raw_D_V_to_field(su3_vector *dest, dsu3_vector *raw, 
			     int milc_parity);

void load_D_links_and_mom_site(QOP_D3_GaugeField **links, QOP_D3_Force **mom,
			       dsu3_matrix ***rawlinks, dsu3_matrix ***rawmom);
void load_D_links_and_mom_field(QOP_D3_GaugeField **links, QOP_D3_Force **mom,
				dsu3_matrix ***rawlinks, dsu3_matrix ***rawmom,
				su3_matrix *srclink, anti_hermitmat *srcmom);
void unload_D_links_and_mom_site(QOP_D3_GaugeField **links, QOP_D3_Force **mom,
			 dsu3_matrix ***rawlinks, dsu3_matrix ***rawmom);
void unload_D_links_and_mom_field(su3_matrix *dstlink, anti_hermitmat *dstmom,
			QOP_D3_GaugeField **links, QOP_D3_Force **mom, 
			dsu3_matrix ***rawlinks, dsu3_matrix ***rawmom);

void load_D_V_from_site(QOP_D3_ColorVector** dest, field_offset src, 
			int parity);
void load_D_V_from_field(QOP_D3_ColorVector** dest, su3_vector *src, 
			int parity);
void unload_D_V_to_site(field_offset dest, QOP_D3_ColorVector* src, 
			int parity);
void unload_D_V_to_field(su3_vector *dest, QOP_D3_ColorVector* src, 
			int parity);

#endif /* GENERIC_QOP_H */
