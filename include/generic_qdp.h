#include "../include/generic.h"
#include <qdp.h>

void set_V_from_field(QDP_ColorVector *dest, field_offset src);
void set_H_from_field(QDP_HalfFermion *dest, field_offset src);
void set_D_from_field(QDP_DiracFermion *dest, field_offset src);
void set_M_from_field(QDP_ColorMatrix *dest, field_offset src);

void set_field_from_V(field_offset dest, QDP_ColorVector *src);
void set_field_from_H(field_offset dest, QDP_HalfFermion *src);
void set_field_from_D(field_offset dest, QDP_DiracFermion *src);
void set_field_from_M(field_offset dest, QDP_ColorMatrix *src);

void set_V_from_temp(QDP_ColorVector *dest, su3_vector *src);
void set_M_from_temp(QDP_ColorMatrix *dest, su3_matrix *src);

void set_temp_from_V(su3_vector *dest, QDP_ColorVector *src);
void set_temp_from_M(su3_matrix *dest, QDP_ColorMatrix *src);

void set4_V_from_temp(QDP_ColorVector *dest[], su3_vector *src);
void set4_M_from_temp(QDP_ColorMatrix *dest[], su3_matrix *src);
