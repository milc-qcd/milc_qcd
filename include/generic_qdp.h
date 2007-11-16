#ifndef _GENERIC_QDP_H
#define _GENERIC_QDP_H
/******************** generic_qdp.h *********************************
*  MIMD version 7 							*
*/

#include <qdp.h>
#include "../include/macros.h"
#include "../include/su3.h"

/* Generic mapping */

#if ( QDP_Precision == 'F' )

#define set_V_from_site    set_F_V_from_site  
#define set_H_from_site	   set_F_H_from_site 
#define set_D_from_site	   set_F_D_from_site 
#define set_M_from_site	   set_F_M_from_site 
#define set4_M_from_site   set4_F_M_from_site

#define set_site_from_V    set_site_from_F_V
#define set_site_from_H	   set_site_from_F_H
#define set_site_from_D	   set_site_from_F_D
#define set_site_from_M	   set_site_from_F_M

#define set_V_from_field   set_F_V_from_field	
#define set_H_from_field   set_F_H_from_field	
#define set_M_from_field   set_F_M_from_field 

#define set_field_from_V   set_field_from_F_V
#define set_field_from_M   set_field_from_F_M

#define set4_field_from_M  set4_field_from_F_M

#define set4_V_from_field  set4_F_V_from_field
#define set4_M_from_field  set4_F_M_from_field

#else

#define set_V_from_site    set_D_V_from_site  
#define set_H_from_site	   set_D_H_from_site 
#define set_D_from_site	   set_D_D_from_site 
#define set_M_from_site	   set_D_M_from_site 
#define set4_M_from_site   set4_D_M_from_site

#define set_site_from_V    set_site_from_D_V
#define set_site_from_H	   set_site_from_D_H
#define set_site_from_D	   set_site_from_D_D
#define set_site_from_M	   set_site_from_D_M

#define set_V_from_field   set_D_V_from_field	
#define set_H_from_field   set_D_H_from_field	
#define set_M_from_field   set_D_M_from_field 

#define set_field_from_V   set_field_from_D_V
#define set_field_from_M   set_field_from_D_M

#define set4_field_from_M  set4_field_from_D_M

#define set4_V_from_field  set4_D_V_from_field
#define set4_M_from_field  set4_D_M_from_field

#endif

/* Mappings for single precision QDP types */

void set_F_V_from_site(QDP_F_ColorVector *dest, field_offset src, int parity);
void set_F_H_from_site(QDP_F_HalfFermion *dest, field_offset src, int parity);
void set_F_D_from_site(QDP_F_DiracFermion *dest, field_offset src, int parity);
void set_F_M_from_site(QDP_F_ColorMatrix *dest, field_offset src, int parity);
void set4_F_M_from_site(QDP_F_ColorMatrix *dest[], field_offset src, int parity);

void set_site_from_F_V(field_offset dest, QDP_F_ColorVector *src, int parity);
void set_site_from_F_H(field_offset dest, QDP_F_HalfFermion *src, int parity);
void set_site_from_F_D(field_offset dest, QDP_F_DiracFermion *src, int parity);
void set_site_from_F_M(field_offset dest, QDP_F_ColorMatrix *src, int parity);

void set_F_V_from_field(QDP_F_ColorVector *dest, su3_vector *src, int parity);
void set_F_H_from_field(QDP_F_HalfFermion *dest, half_wilson_vector *src, int parity);
void set_F_M_from_field(QDP_F_ColorMatrix *dest, su3_matrix *src, int parity);

void set_field_from_F_V(su3_vector *dest, QDP_F_ColorVector *src, int parity);
void set_field_from_F_M(su3_matrix *dest, QDP_F_ColorMatrix *src, int parity);

void set4_field_from_F_M(su3_matrix *dest, QDP_F_ColorMatrix *src[], int parity);

void set4_F_V_from_field(QDP_F_ColorVector *dest[], su3_vector *src, int parity);
void set4_F_M_from_field(QDP_F_ColorMatrix *dest[], su3_matrix *src, int parity);

/* Mappings for double precision QDP types */

void set_D_V_from_site(QDP_D_ColorVector *dest, field_offset src, int parity);
void set_D_H_from_site(QDP_D_HalfFermion *dest, field_offset src, int parity);
void set_D_D_from_site(QDP_D_DiracFermion *dest, field_offset src, int parity);
void set_D_M_from_site(QDP_D_ColorMatrix *dest, field_offset src, int parity);
void set4_D_M_from_site(QDP_D_ColorMatrix *dest[], field_offset src, int parity);

void set_site_from_D_V(field_offset dest, QDP_D_ColorVector *src, int parity);
void set_site_from_D_H(field_offset dest, QDP_D_HalfFermion *src, int parity);
void set_site_from_D_D(field_offset dest, QDP_D_DiracFermion *src, int parity);
void set_site_from_D_M(field_offset dest, QDP_D_ColorMatrix *src, int parity);

void set_D_V_from_field(QDP_D_ColorVector *dest, su3_vector *src, int parity);
void set_D_H_from_field(QDP_D_HalfFermion *dest, half_wilson_vector *src, int parity);
void set_D_M_from_field(QDP_D_ColorMatrix *dest, su3_matrix *src, int parity);

void set_field_from_D_V(su3_vector *dest, QDP_D_ColorVector *src, int parity);
void set_field_from_D_M(su3_matrix *dest, QDP_D_ColorMatrix *src, int parity);

void set4_field_from_D_M(su3_matrix *dest, QDP_D_ColorMatrix *src[], int parity);

void set4_D_V_from_field(QDP_D_ColorVector *dest[], su3_vector *src, int parity);
void set4_D_M_from_field(QDP_D_ColorMatrix *dest[], su3_matrix *src, int parity);

#endif /* _GENERIC_QDP_H */

