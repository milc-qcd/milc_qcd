/************************* RG_operations.c *******************************/
/* MIMD version 7 */
/* F. Maresca Jul 2005 */
#include <stdio.h>
#include <qdp.h>
#include <math.h>
#include "RG_Shamir_includes.h"
#define RG_BLOCK_SUB
#include "RG_include.h"

int intpow(int base,int exp)
{
int y;

  y = 1;

  if (exp > 0)
   while ( exp > 0)
   {
    y = y * base;
    exp = exp - 1;
   }
  if (exp < 0)
   {
   node0_printf("Error. Positive exponent expected\n");
   y = -1;
   }


 return y;

}


int func_block (int x[], void *arg)
{
int i,c=0;


  for(i=0;i<QDP_ndim();i++)
   c += x[i] % fact;

  if ( c == 0) return 0;
  else return 1;

}

void RG_create_block(QDP_Sub_Block *block, int n)
{
QDP_Subset *test;

  block->fact = fact = intpow(2,nrg-n);

  printf("Create a subset for lattice %d x a, this node %d\n",fact,this_node); fflush(stdout);

  if ( n == nrg) block->sub = QDP_all;
  else
  {
   test = QDP_create_subset(func_block,NULL,2);
   block->sub = test[0];
  }

return;
}


void RG_space(QDP_ColorMatrix *a)
{

  if (a == NULL) 
   {
    printf("No space: this node %d\n",this_node);
    terminate(-1);
   }
return;
}

void SQDP_M_eq_M(QDP_ColorMatrix *a ,QDP_ColorMatrix *b, QDP_Sub_Block s)
{

fact = s.fact;
QDP_M_eq_M(a,b,s.sub);

return;
}

void SQDP_M_peq_M(QDP_ColorMatrix *a ,QDP_ColorMatrix *b, QDP_Sub_Block s)
{

fact = s.fact;
QDP_M_peq_M(a,b,s.sub);

return;
}

void SQDP_M_eq_c(QDP_ColorMatrix *a ,QLA_Complex *c, QDP_Sub_Block s)
{

fact = s.fact;
QDP_M_eq_c(a,c,s.sub);

return;
}


void SQDP_M_eq_sM(QDP_ColorMatrix *a ,QDP_ColorMatrix *b, QDP_Shift shift, QDP_ShiftDir dir, QDP_Sub_Block s)
{

fact = s.fact;
QDP_M_eq_sM(a, b, shift, dir, s.sub);

return;
}

void SQDP_M_eq_func(QDP_ColorMatrix *a, void (*func)(QLA_ColorMatrix *gl, int coord[]), QDP_Sub_Block s)
{

fact = s.fact;
QDP_M_eq_func(a,func,s.sub);

return;
}

void SQDP_M_eq_zero(QDP_ColorMatrix *a, QDP_Sub_Block s)
{

fact = s.fact;
QDP_M_eq_zero(a,s.sub);

return;
}

void SQDP_M_eq_M_times_M(QDP_ColorMatrix *a,QDP_ColorMatrix *b,QDP_ColorMatrix *c, QDP_Sub_Block s)
{

fact = s.fact;
QDP_M_eq_M_times_M(a,b,c,s.sub);
 
return;
}

void SQDP_M_eq_M_times_Ma(QDP_ColorMatrix *a,QDP_ColorMatrix *b,QDP_ColorMatrix *c, QDP_Sub_Block s)
{

fact = s.fact;
QDP_M_eq_M_times_Ma(a,b,c,s.sub);
 
return;
}

void SQDP_M_peq_M_times_M(QDP_ColorMatrix *a,QDP_ColorMatrix *b,QDP_ColorMatrix *c, QDP_Sub_Block s)
{

fact = s.fact;
QDP_M_peq_M_times_M(a,b,c,s.sub);
 
return;
}


void SQDP_M_eq_Ma_times_M(QDP_ColorMatrix *a,QDP_ColorMatrix *b,QDP_ColorMatrix *c, QDP_Sub_Block s)
{

fact = s.fact;
QDP_M_eq_Ma_times_M(a,b,c,s.sub);
 
return;
}

void SQDP_M_eq_M_minus_M(QDP_ColorMatrix *a,QDP_ColorMatrix *b,QDP_ColorMatrix *c, QDP_Sub_Block s)
{

fact = s.fact;
QDP_M_eq_M_minus_M(a,b,c,s.sub);
 
return;
}


void SQDP_M_eq_M_times_sM(QDP_ColorMatrix *a,QDP_ColorMatrix *b,QDP_ColorMatrix *c, QDP_Shift shift, QDP_ShiftDir dir,QDP_Sub_Block s)
{

fact = s.fact;
QDP_M_eq_M_times_sM(a,b,c,shift,dir,s.sub);
 
return;
}

void SQDP_M_peq_M_times_sM(QDP_ColorMatrix *a,QDP_ColorMatrix *b,QDP_ColorMatrix *c, QDP_Shift shift, QDP_ShiftDir dir,QDP_Sub_Block s)
{

fact = s.fact;
QDP_M_peq_M_times_sM(a,b,c,shift,dir,s.sub);
 
return;
}

void SQDP_M_eq_r_times_M(QDP_ColorMatrix *a,QLA_Real *b, QDP_ColorMatrix *c, QDP_Sub_Block s)
{

fact = s.fact;
QDP_M_eq_r_times_M(a,b,c,s.sub);

return;
}

void SQDP_M_peq_r_times_M(QDP_ColorMatrix *a,QLA_Real *b, QDP_ColorMatrix *c, QDP_Sub_Block s)
{

fact = s.fact;
QDP_M_peq_r_times_M(a,b,c,s.sub);

return;
}

void SQDP_M_eq_r_times_M_plus_M(QDP_ColorMatrix *a,QLA_Real *b, QDP_ColorMatrix *c,QDP_ColorMatrix *d,QDP_Sub_Block s)
{

fact = s.fact;
QDP_M_eq_r_times_M_plus_M(a,b,c,d,s.sub);

return;
}

void SQDP_R_eq_re_M_dot_M(QLA_Real *r, QDP_ColorMatrix *a, QDP_ColorMatrix *b, QDP_Sub_Block s)
{

fact = s.fact;
QDP_r_eq_re_M_dot_M(r, a, b, s.sub);

return;
}

void SQDP_V_eq_V(QDP_ColorVector *a ,QDP_ColorVector *b, QDP_Sub_Block s)
{

  fact = s.fact;
  QDP_V_eq_V(a,b,s.sub);
  
  return;
}

void SQDP_V_peq_V(QDP_ColorVector *a ,QDP_ColorVector *b, QDP_Sub_Block s)
{

fact = s.fact;
QDP_V_peq_V(a,b,s.sub);

return;
}


void SQDP_V_eq_sV(QDP_ColorVector *a ,QDP_ColorVector *b, QDP_Shift shift, QDP_ShiftDir dir, QDP_Sub_Block s)
{

fact = s.fact;
QDP_V_eq_sV(a, b, shift, dir, s.sub);

return;
}

void SQDP_V_eq_func(QDP_ColorVector *a, void (*func)(QLA_ColorVector *gl, int coord[]), QDP_Sub_Block s)
{

fact = s.fact;
QDP_V_eq_func(a,func,s.sub);

return;
}

void SQDP_V_eq_zero(QDP_ColorVector *a, QDP_Sub_Block s)
{

fact = s.fact;
QDP_V_eq_zero(a,s.sub);

return;
}

void SQDP_V_eq_M_times_sV(QDP_ColorVector *a, QDP_ColorMatrix *b, QDP_ColorVector *c, QDP_Shift shift, QDP_ShiftDir dir, QDP_Sub_Block s)
{

fact = s.fact;
QDP_V_eq_M_times_sV(a,b,c,shift,dir,s.sub);

return;
}

void SQDP_V_peq_r_times_V(QDP_ColorVector *a, QLA_Real *r, QDP_ColorVector *b, QDP_Sub_Block s)
{

fact = s.fact;
QDP_V_peq_r_times_V(a,r,b,s.sub);

return;
}

void SQDP_V_eq_r_times_V_minus_V(QDP_ColorVector *a, QLA_Real *r, QDP_ColorVector *b, QDP_ColorVector *c, QDP_Sub_Block s)
{

fact = s.fact;
QDP_V_eq_r_times_V_minus_V(a,r,b,c,s.sub);

return;
}

void SQDP_V_eq_Ma_times_V(QDP_ColorVector *a, QDP_ColorMatrix *b, QDP_ColorVector *c, QDP_Sub_Block s)
{

fact = s.fact;
QDP_V_eq_Ma_times_V(a,b,c,s.sub);

return;
}

void SQDP_V_eq_r_times_V(QDP_ColorVector *a, QLA_Real *r, QDP_ColorVector *b, QDP_Sub_Block s)
{

fact = s.fact;
QDP_V_eq_r_times_V(a,r,b,s.sub);

return;
}

void SQDP_r_eq_re_M_dot_M(QLA_Real *r, QDP_ColorMatrix *a, QDP_ColorMatrix *b, QDP_Sub_Block s)
{

fact = s.fact;
QDP_r_eq_re_M_dot_M(r,a,b,s.sub);

return;
}
