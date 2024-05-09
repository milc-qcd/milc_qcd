#ifndef __blind_data_h__
#define __blind_data_h__

#include "../include/complex.h"

extern const char* blind_short_tag;

void blind_vfloat(float v[],int veclen);

void blind_vdouble(double v[],int veclen);

//void blind_vfcomplex(float complex v[],int veclen); // original
void blind_vfcomplex(complex v[],int veclen); // compiles

//void blind_vdcomplex(double complex v[],int veclen); // original
void blind_vdcomplex(dcomplex v[],int veclen); // compiles

#endif
