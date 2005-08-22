//#include <stdio.h>
#include <qdp.h>
//#include <math.h>
//#include "RG_Shamir_includes.h"
//#include "RG_include.h"

typedef struct {
  QDP_Subset sub;
  int fact;
} QDP_Sub_Block;

void SQDP_M_eq_M(QDP_ColorMatrix *a ,QDP_ColorMatrix *b, QDP_Sub_Block s);

void SQDP_M_eq_c(QDP_ColorMatrix *a ,QDP_Complex *c, QDP_Sub_Block s);

void SQDP_M_eq_sM(QDP_ColorMatrix *a ,QDP_ColorMatrix *b, QDP_Shift shift, QDP_ShiftDir dir, QDP_Sub_Block s);

void SQDP_M_eq_fun(QDP_ColorMatrix *a, void (*func)(QLA_ColorMatrix *gl, int coord[]), QDP_Sub_Block s);

void RG_create_block(QDP_Sub_Block *block, int n);

#ifdef RG_BLOCK_SUB
#define EXTERN_OP
#else
#define EXTERN_OP extern
#endif

EXTERN_OP int fact;

