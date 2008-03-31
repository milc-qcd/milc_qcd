#define _inline_sse_mult_su3_mat_vec_sum_4dir(aa,bb0,bb1,bb2,bb3,cc) \
{ \
__asm__ __volatile__ ("movlps %0, %%xmm14 \n\t" \
                      "movaps %%xmm14, %%xmm0 \n\t" \
                      "unpcklps %%xmm0, %%xmm0 \n\t" \
                      "movhlps %%xmm0, %%xmm11 \n\t" \
                      "movlhps %%xmm0, %%xmm0 \n\t" \
                      "movhps %1, %%xmm14 \n\t" \
                      "movaps %%xmm14, %%xmm6 \n\t" \
                      "unpckhps %%xmm6, %%xmm6 \n\t" \
                      "movhlps %%xmm6, %%xmm12 \n\t" \
                      "movlhps %%xmm6, %%xmm6 \n\t" \
                      "movhps %2, %%xmm15 \n\t" \
                      "movaps %%xmm15, %%xmm7 \n\t" \
                      "unpckhps %%xmm7, %%xmm7 \n\t" \
                      "movhlps %%xmm7, %%xmm13 \n\t" \
                      "movlhps %%xmm7, %%xmm7 \n\t" \
                      "movlps %3, %%xmm1 \n\t" \
                      "movhps %4, %%xmm1 \n\t" \
                      "mulps %%xmm1, %%xmm0 \n\t" \
                      "movlps %5, %%xmm9" \
                      : \
                      : \
                      "m" ((bb0)->c[0]), \
                      "m" ((bb0)->c[1]), \
                      "m" ((bb0)->c[2]), \
                      "m" ((aa)[0].e[0][0]), \
                      "m" ((aa)[0].e[1][0]), \
                      "m" ((aa)[0].e[0][1])\
                      : \
                      "%xmm0", \
                      "%xmm1", \
                      "%xmm11", \
                      "%xmm12", \
                      "%xmm13", \
                      "%xmm14", \
                      "%xmm15", \
                      "%xmm6", \
                      "%xmm7", \
                      "%xmm9", \
                      "memory"); \
__asm__ __volatile__ ("movhps %0, %%xmm9 \n\t" \
                      "mulps %%xmm9, %%xmm6 \n\t" \
                      "movlps %1, %%xmm10 \n\t" \
                      "movhps %2, %%xmm10 \n\t" \
                      "mulps %%xmm10, %%xmm7 \n\t" \
                      "addps %%xmm6, %%xmm0 \n\t" \
                      "addps %%xmm7, %%xmm0 \n\t" \
                      "movlhps %%xmm11, %%xmm11 \n\t" \
                      "movlhps %%xmm12, %%xmm12 \n\t" \
                      "movlhps %%xmm13, %%xmm13 \n\t" \
                      "mulps %%xmm11, %%xmm1 \n\t" \
                      "mulps %%xmm12, %%xmm9 \n\t" \
                      "mulps %%xmm13, %%xmm10 \n\t" \
                      "addps %%xmm9, %%xmm1 \n\t" \
                      "addps %%xmm10, %%xmm1 \n\t" \
                      "movlps %3, %%xmm6 \n\t" \
                      "unpcklps %%xmm6, %%xmm6 \n\t" \
                      "movaps %%xmm6, %%xmm2 \n\t" \
                      "movlps %4, %%xmm3 \n\t" \
                      "unpcklps %%xmm3, %%xmm3 \n\t" \
                      "movlhps %%xmm3, %%xmm2 \n\t" \
                      "movhlps %%xmm6, %%xmm3 \n\t" \
                      "mulps %%xmm14, %%xmm2 \n\t" \
                      "mulps %%xmm14, %%xmm3 \n\t" \
                      "movlps %5, %%xmm4" \
                      : \
                      : \
                      "m" ((aa)[0].e[1][1]), \
                      "m" ((aa)[0].e[0][2]), \
                      "m" ((aa)[0].e[1][2]), \
                      "m" ((aa)[0].e[2][0]), \
                      "m" ((aa)[0].e[2][1]), \
                      "m" ((aa)[0].e[2][2])\
                      : \
                      "%xmm0", \
                      "%xmm1", \
                      "%xmm10", \
                      "%xmm11", \
                      "%xmm12", \
                      "%xmm13", \
                      "%xmm14", \
                      "%xmm2", \
                      "%xmm3", \
                      "%xmm4", \
                      "%xmm6", \
                      "%xmm7", \
                      "%xmm9", \
                      "memory"); \
__asm__ __volatile__ ("unpcklps %%xmm4, %%xmm4 \n\t" \
                      "movhlps %%xmm15, %%xmm15 \n\t" \
                      "mulps %%xmm15, %%xmm4 \n\t" \
                      "movlps %0, %%xmm14 \n\t" \
                      "movaps %%xmm14, %%xmm5 \n\t" \
                      "unpcklps %%xmm5, %%xmm5 \n\t" \
                      "movhlps %%xmm5, %%xmm11 \n\t" \
                      "movlhps %%xmm5, %%xmm5 \n\t" \
                      "movhps %1, %%xmm14 \n\t" \
                      "movaps %%xmm14, %%xmm6 \n\t" \
                      "unpckhps %%xmm6, %%xmm6 \n\t" \
                      "movhlps %%xmm6, %%xmm12 \n\t" \
                      "movlhps %%xmm6, %%xmm6 \n\t" \
                      "movhps %2, %%xmm15 \n\t" \
                      "movaps %%xmm15, %%xmm7 \n\t" \
                      "unpckhps %%xmm7, %%xmm7 \n\t" \
                      "movhlps %%xmm7, %%xmm13 \n\t" \
                      "movlhps %%xmm7, %%xmm7 \n\t" \
                      "movlps %3, %%xmm8 \n\t" \
                      "movhps %4, %%xmm8 \n\t" \
                      "mulps %%xmm8, %%xmm5 \n\t" \
                      "movlps %5, %%xmm9" \
                      : \
                      : \
                      "m" ((bb1)->c[0]), \
                      "m" ((bb1)->c[1]), \
                      "m" ((bb1)->c[2]), \
                      "m" ((aa)[1].e[0][0]), \
                      "m" ((aa)[1].e[1][0]), \
                      "m" ((aa)[1].e[0][1])\
                      : \
                      "%xmm11", \
                      "%xmm12", \
                      "%xmm13", \
                      "%xmm14", \
                      "%xmm15", \
                      "%xmm4", \
                      "%xmm5", \
                      "%xmm6", \
                      "%xmm7", \
                      "%xmm8", \
                      "%xmm9", \
                      "memory"); \
__asm__ __volatile__ ("movhps %0, %%xmm9 \n\t" \
                      "mulps %%xmm9, %%xmm6 \n\t" \
                      "movlps %1, %%xmm10 \n\t" \
                      "movhps %2, %%xmm10 \n\t" \
                      "mulps %%xmm10, %%xmm7 \n\t" \
                      "addps %%xmm6, %%xmm5 \n\t" \
                      "addps %%xmm7, %%xmm5 \n\t" \
                      "addps %%xmm5, %%xmm0 \n\t" \
                      "movlhps %%xmm11, %%xmm11 \n\t" \
                      "movlhps %%xmm12, %%xmm12 \n\t" \
                      "movlhps %%xmm13, %%xmm13 \n\t" \
                      "mulps %%xmm11, %%xmm8 \n\t" \
                      "mulps %%xmm12, %%xmm9 \n\t" \
                      "mulps %%xmm13, %%xmm10 \n\t" \
                      "addps %%xmm9, %%xmm8 \n\t" \
                      "addps %%xmm10, %%xmm8 \n\t" \
                      "addps %%xmm8, %%xmm1 \n\t" \
                      "movlps %3, %%xmm6 \n\t" \
                      "unpcklps %%xmm6, %%xmm6 \n\t" \
                      "movaps %%xmm6, %%xmm7 \n\t" \
                      "movlps %4, %%xmm5 \n\t" \
                      "unpcklps %%xmm5, %%xmm5 \n\t" \
                      "movlhps %%xmm5, %%xmm7 \n\t" \
                      "movhlps %%xmm6, %%xmm5 \n\t" \
                      "mulps %%xmm14, %%xmm7 \n\t" \
                      "mulps %%xmm14, %%xmm5 \n\t" \
                      "addps %%xmm7, %%xmm2 \n\t" \
                      "addps %%xmm5, %%xmm3 \n\t" \
                      "movlps %5, %%xmm11" \
                      : \
                      : \
                      "m" ((aa)[1].e[1][1]), \
                      "m" ((aa)[1].e[0][2]), \
                      "m" ((aa)[1].e[1][2]), \
                      "m" ((aa)[1].e[2][0]), \
                      "m" ((aa)[1].e[2][1]), \
                      "m" ((aa)[1].e[2][2])\
                      : \
                      "%xmm0", \
                      "%xmm1", \
                      "%xmm10", \
                      "%xmm11", \
                      "%xmm12", \
                      "%xmm13", \
                      "%xmm14", \
                      "%xmm2", \
                      "%xmm3", \
                      "%xmm5", \
                      "%xmm6", \
                      "%xmm7", \
                      "%xmm8", \
                      "%xmm9", \
                      "memory"); \
__asm__ __volatile__ ("unpcklps %%xmm11, %%xmm11 \n\t" \
                      "movhlps %%xmm15, %%xmm15 \n\t" \
                      "mulps %%xmm15, %%xmm11 \n\t" \
                      "addps %%xmm11, %%xmm4 \n\t" \
                      "movlps %0, %%xmm14 \n\t" \
                      "movaps %%xmm14, %%xmm5 \n\t" \
                      "unpcklps %%xmm5, %%xmm5 \n\t" \
                      "movhlps %%xmm5, %%xmm11 \n\t" \
                      "movlhps %%xmm5, %%xmm5 \n\t" \
                      "movhps %1, %%xmm14 \n\t" \
                      "movaps %%xmm14, %%xmm6 \n\t" \
                      "unpckhps %%xmm6, %%xmm6 \n\t" \
                      "movhlps %%xmm6, %%xmm12 \n\t" \
                      "movlhps %%xmm6, %%xmm6 \n\t" \
                      "movhps %2, %%xmm15 \n\t" \
                      "movaps %%xmm15, %%xmm7 \n\t" \
                      "unpckhps %%xmm7, %%xmm7 \n\t" \
                      "movhlps %%xmm7, %%xmm13 \n\t" \
                      "movlhps %%xmm7, %%xmm7 \n\t" \
                      "movlps %3, %%xmm8 \n\t" \
                      "movhps %4, %%xmm8 \n\t" \
                      "mulps %%xmm8, %%xmm5 \n\t" \
                      "movlps %5, %%xmm9" \
                      : \
                      : \
                      "m" ((bb2)->c[0]), \
                      "m" ((bb2)->c[1]), \
                      "m" ((bb2)->c[2]), \
                      "m" ((aa)[2].e[0][0]), \
                      "m" ((aa)[2].e[1][0]), \
                      "m" ((aa)[2].e[0][1])\
                      : \
                      "%xmm11", \
                      "%xmm12", \
                      "%xmm13", \
                      "%xmm14", \
                      "%xmm15", \
                      "%xmm4", \
                      "%xmm5", \
                      "%xmm6", \
                      "%xmm7", \
                      "%xmm8", \
                      "%xmm9", \
                      "memory"); \
__asm__ __volatile__ ("movhps %0, %%xmm9 \n\t" \
                      "mulps %%xmm9, %%xmm6 \n\t" \
                      "movlps %1, %%xmm10 \n\t" \
                      "movhps %2, %%xmm10 \n\t" \
                      "mulps %%xmm10, %%xmm7 \n\t" \
                      "addps %%xmm6, %%xmm5 \n\t" \
                      "addps %%xmm7, %%xmm5 \n\t" \
                      "addps %%xmm5, %%xmm0 \n\t" \
                      "movlhps %%xmm11, %%xmm11 \n\t" \
                      "movlhps %%xmm12, %%xmm12 \n\t" \
                      "movlhps %%xmm13, %%xmm13 \n\t" \
                      "mulps %%xmm11, %%xmm8 \n\t" \
                      "mulps %%xmm12, %%xmm9 \n\t" \
                      "mulps %%xmm13, %%xmm10 \n\t" \
                      "addps %%xmm9, %%xmm8 \n\t" \
                      "addps %%xmm10, %%xmm8 \n\t" \
                      "addps %%xmm8, %%xmm1 \n\t" \
                      "movlps %3, %%xmm6 \n\t" \
                      "unpcklps %%xmm6, %%xmm6 \n\t" \
                      "movaps %%xmm6, %%xmm7 \n\t" \
                      "movlps %4, %%xmm5 \n\t" \
                      "unpcklps %%xmm5, %%xmm5 \n\t" \
                      "movlhps %%xmm5, %%xmm7 \n\t" \
                      "movhlps %%xmm6, %%xmm5 \n\t" \
                      "mulps %%xmm14, %%xmm7 \n\t" \
                      "mulps %%xmm14, %%xmm5 \n\t" \
                      "addps %%xmm7, %%xmm2 \n\t" \
                      "addps %%xmm5, %%xmm3 \n\t" \
                      "movlps %5, %%xmm11" \
                      : \
                      : \
                      "m" ((aa)[2].e[1][1]), \
                      "m" ((aa)[2].e[0][2]), \
                      "m" ((aa)[2].e[1][2]), \
                      "m" ((aa)[2].e[2][0]), \
                      "m" ((aa)[2].e[2][1]), \
                      "m" ((aa)[2].e[2][2])\
                      : \
                      "%xmm0", \
                      "%xmm1", \
                      "%xmm10", \
                      "%xmm11", \
                      "%xmm12", \
                      "%xmm13", \
                      "%xmm14", \
                      "%xmm2", \
                      "%xmm3", \
                      "%xmm5", \
                      "%xmm6", \
                      "%xmm7", \
                      "%xmm8", \
                      "%xmm9", \
                      "memory"); \
__asm__ __volatile__ ("unpcklps %%xmm11, %%xmm11 \n\t" \
                      "movhlps %%xmm15, %%xmm15 \n\t" \
                      "mulps %%xmm15, %%xmm11 \n\t" \
                      "addps %%xmm11, %%xmm4 \n\t" \
                      "movlps %0, %%xmm14 \n\t" \
                      "movaps %%xmm14, %%xmm5 \n\t" \
                      "unpcklps %%xmm5, %%xmm5 \n\t" \
                      "movhlps %%xmm5, %%xmm11 \n\t" \
                      "movlhps %%xmm5, %%xmm5 \n\t" \
                      "movhps %1, %%xmm14 \n\t" \
                      "movaps %%xmm14, %%xmm6 \n\t" \
                      "unpckhps %%xmm6, %%xmm6 \n\t" \
                      "movhlps %%xmm6, %%xmm12 \n\t" \
                      "movlhps %%xmm6, %%xmm6 \n\t" \
                      "movhps %2, %%xmm15 \n\t" \
                      "movaps %%xmm15, %%xmm7 \n\t" \
                      "unpckhps %%xmm7, %%xmm7 \n\t" \
                      "movhlps %%xmm7, %%xmm13 \n\t" \
                      "movlhps %%xmm7, %%xmm7 \n\t" \
                      "movlps %3, %%xmm8 \n\t" \
                      "movhps %4, %%xmm8 \n\t" \
                      "mulps %%xmm8, %%xmm5 \n\t" \
                      "movlps %5, %%xmm9" \
                      : \
                      : \
                      "m" ((bb3)->c[0]), \
                      "m" ((bb3)->c[1]), \
                      "m" ((bb3)->c[2]), \
                      "m" ((aa)[3].e[0][0]), \
                      "m" ((aa)[3].e[1][0]), \
                      "m" ((aa)[3].e[0][1])\
                      : \
                      "%xmm11", \
                      "%xmm12", \
                      "%xmm13", \
                      "%xmm14", \
                      "%xmm15", \
                      "%xmm4", \
                      "%xmm5", \
                      "%xmm6", \
                      "%xmm7", \
                      "%xmm8", \
                      "%xmm9", \
                      "memory"); \
__asm__ __volatile__ ("movhps %0, %%xmm9 \n\t" \
                      "mulps %%xmm9, %%xmm6 \n\t" \
                      "movlps %1, %%xmm10 \n\t" \
                      "movhps %2, %%xmm10 \n\t" \
                      "mulps %%xmm10, %%xmm7 \n\t" \
                      "addps %%xmm6, %%xmm5 \n\t" \
                      "addps %%xmm7, %%xmm5 \n\t" \
                      "addps %%xmm5, %%xmm0 \n\t" \
                      "movlhps %%xmm11, %%xmm11 \n\t" \
                      "movlhps %%xmm12, %%xmm12 \n\t" \
                      "movlhps %%xmm13, %%xmm13 \n\t" \
                      "mulps %%xmm11, %%xmm8 \n\t" \
                      "mulps %%xmm12, %%xmm9 \n\t" \
                      "mulps %%xmm13, %%xmm10 \n\t" \
                      "addps %%xmm9, %%xmm8 \n\t" \
                      "addps %%xmm10, %%xmm8 \n\t" \
                      "addps %%xmm8, %%xmm1 \n\t" \
                      "xorps %3, %%xmm1 \n\t" \
                      "shufps $0xb1, %%xmm1, %%xmm1 \n\t" \
                      "addps %%xmm1, %%xmm0 \n\t" \
                      : \
                      : \
                      "m" ((aa)[3].e[1][1]), \
                      "m" ((aa)[3].e[0][2]), \
                      "m" ((aa)[3].e[1][2]), \
                      "m" (_sse_sgn24)\
                      : \
                      "%xmm0", \
                      "%xmm1", \
                      "%xmm10", \
                      "%xmm11", \
                      "%xmm12", \
                      "%xmm13", \
                      "%xmm5", \
                      "%xmm6", \
                      "%xmm7", \
                      "%xmm8", \
                      "%xmm9", \
                      "memory"); \
__asm__ __volatile__ ("movlps %%xmm0, %0 \n\t" \
                      "movhps %%xmm0, %1 \n\t" \
                      : \
                      "=m" ((cc)->c[0]), \
                      "=m" ((cc)->c[1])\
                      : \
                      : \
                      "%xmm0", \
                      "memory"); \
__asm__ __volatile__ ("movlps %0, %%xmm6 \n\t" \
                      "unpcklps %%xmm6, %%xmm6 \n\t" \
                      "movaps %%xmm6, %%xmm7 \n\t" \
                      "movlps %1, %%xmm5 \n\t" \
                      "unpcklps %%xmm5, %%xmm5 \n\t" \
                      "movlhps %%xmm5, %%xmm7 \n\t" \
                      "movhlps %%xmm6, %%xmm5 \n\t" \
                      "mulps %%xmm14, %%xmm7 \n\t" \
                      "mulps %%xmm14, %%xmm5 \n\t" \
                      "addps %%xmm7, %%xmm2 \n\t" \
                      "addps %%xmm5, %%xmm3 \n\t" \
                      "xorps %2, %%xmm3 \n\t" \
                      "shufps $0xb1, %%xmm3, %%xmm3 \n\t" \
                      "addps %%xmm3, %%xmm2 \n\t" \
                      "movlps %3, %%xmm11 \n\t" \
                      "unpcklps %%xmm11, %%xmm11 \n\t" \
                      "movhlps %%xmm15, %%xmm15 \n\t" \
                      "mulps %%xmm15, %%xmm11 \n\t" \
                      "addps %%xmm11, %%xmm4 \n\t" \
                      "xorps %4, %%xmm4 \n\t" \
                      "shufps $0xb4, %%xmm4, %%xmm4 \n\t" \
                      "addps %%xmm4, %%xmm2 \n\t" \
                      "movhlps %%xmm2, %%xmm10 \n\t" \
                      "addps %%xmm10, %%xmm2 \n\t" \
                      : \
                      : \
                      "m" ((aa)[3].e[2][0]), \
                      "m" ((aa)[3].e[2][1]), \
                      "m" (_sse_sgn24), \
                      "m" ((aa)[3].e[2][2]), \
                      "m" (_sse_sgn4)\
                      : \
                      "%xmm10", \
                      "%xmm11", \
                      "%xmm14", \
                      "%xmm15", \
                      "%xmm2", \
                      "%xmm3", \
                      "%xmm4", \
                      "%xmm5", \
                      "%xmm6", \
                      "%xmm7", \
                      "memory"); \
__asm__ __volatile__ ("movlps %%xmm2, %0 \n\t" \
                      : \
                      "=m" ((cc)->c[2])\
                      : \
                      : \
                      "%xmm2", \
                      "memory"); \
}
