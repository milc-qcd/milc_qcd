#define _inline_sse_mult_adj_su3_mat_vec_4dir(aa,bb,cc) \
{ \
__asm__ __volatile__ ("movlps %0, %%xmm14 \n\t" \
                      "movaps %%xmm14, %%xmm8 \n\t" \
                      "unpcklps %%xmm8, %%xmm8 \n\t" \
                      "movhlps %%xmm8, %%xmm11 \n\t" \
                      "movlhps %%xmm8, %%xmm8 \n\t" \
                      "movhps %1, %%xmm14 \n\t" \
                      "movaps %%xmm14, %%xmm9 \n\t" \
                      "unpckhps %%xmm9, %%xmm9 \n\t" \
                      "movhlps %%xmm9, %%xmm12 \n\t" \
                      "movlhps %%xmm9, %%xmm9 \n\t" \
                      "movhps %2, %%xmm15 \n\t" \
                      "movaps %%xmm15, %%xmm10 \n\t" \
                      "unpckhps %%xmm10, %%xmm10 \n\t" \
                      "movhlps %%xmm10, %%xmm13 \n\t" \
                      "movlhps %%xmm10, %%xmm10 \n\t" \
                      "movaps %%xmm8, %%xmm0 \n\t" \
                      "movlps %3, %%xmm3 \n\t" \
                      "movhps %4, %%xmm3 \n\t" \
                      "mulps %%xmm3, %%xmm0 \n\t" \
                      "movaps %%xmm9, %%xmm1 \n\t" \
                      "movlps %5, %%xmm4" \
                      : \
                      : \
                      "m" ((bb)->c[0]), \
                      "m" ((bb)->c[1]), \
                      "m" ((bb)->c[2]), \
                      "m" ((aa)[0].e[0][0]), \
                      "m" ((aa)[0].e[0][1]), \
                      "m" ((aa)[0].e[1][0])\
                      : \
                      "%xmm0", \
                      "%xmm1", \
                      "%xmm10", \
                      "%xmm11", \
                      "%xmm12", \
                      "%xmm13", \
                      "%xmm14", \
                      "%xmm15", \
                      "%xmm3", \
                      "%xmm4", \
                      "%xmm8", \
                      "%xmm9", \
                      "memory"); \
__asm__ __volatile__ ("movhps %0, %%xmm4 \n\t" \
                      "mulps %%xmm4, %%xmm1 \n\t" \
                      "movaps %%xmm10, %%xmm2 \n\t" \
                      "movlps %1, %%xmm5 \n\t" \
                      "movhps %2, %%xmm5 \n\t" \
                      "mulps %%xmm5, %%xmm2 \n\t" \
                      "addps %%xmm1, %%xmm0 \n\t" \
                      "addps %%xmm2, %%xmm0 \n\t" \
                      "movlhps %%xmm11, %%xmm11 \n\t" \
                      "movlhps %%xmm12, %%xmm12 \n\t" \
                      "movlhps %%xmm13, %%xmm13 \n\t" \
                      "mulps %%xmm11, %%xmm3 \n\t" \
                      "mulps %%xmm12, %%xmm4 \n\t" \
                      "mulps %%xmm13, %%xmm5 \n\t" \
                      "addps %%xmm4, %%xmm3 \n\t" \
                      "addps %%xmm5, %%xmm3 \n\t" \
                      "xorps %3, %%xmm0 \n\t" \
                      "shufps $0xb1, %%xmm3, %%xmm3 \n\t" \
                      "addps %%xmm3, %%xmm0 \n\t" \
                      : \
                      : \
                      "m" ((aa)[0].e[1][1]), \
                      "m" ((aa)[0].e[2][0]), \
                      "m" ((aa)[0].e[2][1]), \
                      "m" (_sse_sgn24)\
                      : \
                      "%xmm0", \
                      "%xmm1", \
                      "%xmm10", \
                      "%xmm11", \
                      "%xmm12", \
                      "%xmm13", \
                      "%xmm2", \
                      "%xmm3", \
                      "%xmm4", \
                      "%xmm5", \
                      "memory"); \
__asm__ __volatile__ ("movlps %%xmm0, %0 \n\t" \
                      "movhps %%xmm0, %1 \n\t" \
                      : \
                      "=m" ((cc)[0].c[0]), \
                      "=m" ((cc)[0].c[1])\
                      : \
                      : \
                      "%xmm0", \
                      "memory"); \
__asm__ __volatile__ ("movlps %0, %%xmm1 \n\t" \
                      "unpcklps %%xmm1, %%xmm1 \n\t" \
                      "movaps %%xmm1, %%xmm6 \n\t" \
                      "movlps %1, %%xmm2 \n\t" \
                      "unpcklps %%xmm2, %%xmm2 \n\t" \
                      "movlhps %%xmm2, %%xmm6 \n\t" \
                      "movhlps %%xmm1, %%xmm2 \n\t" \
                      "mulps %%xmm14, %%xmm6 \n\t" \
                      "mulps %%xmm14, %%xmm2 \n\t" \
                      "shufps $0xb1, %%xmm2, %%xmm2 \n\t" \
                      "xorps %2, %%xmm2 \n\t" \
                      "addps %%xmm2, %%xmm6 \n\t" \
                      "movlps %3, %%xmm7 \n\t" \
                      "unpcklps %%xmm7, %%xmm7 \n\t" \
                      "movhlps %%xmm15, %%xmm15 \n\t" \
                      "mulps %%xmm15, %%xmm7 \n\t" \
                      "shufps $0xb4, %%xmm7, %%xmm7 \n\t" \
                      "xorps %4, %%xmm7 \n\t" \
                      "addps %%xmm7, %%xmm6 \n\t" \
                      "movhlps %%xmm6, %%xmm5 \n\t" \
                      "addps %%xmm5, %%xmm6 \n\t" \
                      : \
                      : \
                      "m" ((aa)[0].e[0][2]), \
                      "m" ((aa)[0].e[1][2]), \
                      "m" (_sse_sgn24), \
                      "m" ((aa)[0].e[2][2]), \
                      "m" (_sse_sgn4)\
                      : \
                      "%xmm1", \
                      "%xmm14", \
                      "%xmm15", \
                      "%xmm2", \
                      "%xmm5", \
                      "%xmm6", \
                      "%xmm7", \
                      "memory"); \
__asm__ __volatile__ ("movlps %%xmm6, %0 \n\t" \
                      "movaps %%xmm8, %%xmm0 \n\t" \
                      : \
                      "=m" ((cc)[0].c[2])\
                      : \
                      : \
                      "%xmm0", \
                      "%xmm6", \
                      "%xmm8", \
                      "memory"); \
__asm__ __volatile__ ("movlps %0, %%xmm3 \n\t" \
                      "movhps %1, %%xmm3 \n\t" \
                      "mulps %%xmm3, %%xmm0 \n\t" \
                      "movaps %%xmm9, %%xmm1 \n\t" \
                      "movlps %2, %%xmm4 \n\t" \
                      "movhps %3, %%xmm4 \n\t" \
                      "mulps %%xmm4, %%xmm1 \n\t" \
                      "movaps %%xmm10, %%xmm2 \n\t" \
                      "movlps %4, %%xmm5 \n\t" \
                      "movhps %5, %%xmm5" \
                      : \
                      : \
                      "m" ((aa)[1].e[0][0]), \
                      "m" ((aa)[1].e[0][1]), \
                      "m" ((aa)[1].e[1][0]), \
                      "m" ((aa)[1].e[1][1]), \
                      "m" ((aa)[1].e[2][0]), \
                      "m" ((aa)[1].e[2][1])\
                      : \
                      "%xmm0", \
                      "%xmm1", \
                      "%xmm10", \
                      "%xmm2", \
                      "%xmm3", \
                      "%xmm4", \
                      "%xmm5", \
                      "%xmm9", \
                      "memory"); \
__asm__ __volatile__ ("mulps %%xmm5, %%xmm2 \n\t" \
                      "addps %%xmm1, %%xmm0 \n\t" \
                      "addps %%xmm2, %%xmm0 \n\t" \
                      "mulps %%xmm11, %%xmm3 \n\t" \
                      "mulps %%xmm12, %%xmm4 \n\t" \
                      "mulps %%xmm13, %%xmm5 \n\t" \
                      "addps %%xmm4, %%xmm3 \n\t" \
                      "addps %%xmm5, %%xmm3 \n\t" \
                      "xorps %0, %%xmm0 \n\t" \
                      "shufps $0xb1, %%xmm3, %%xmm3 \n\t" \
                      "addps %%xmm3, %%xmm0 \n\t" \
                      : \
                      : \
                      "m" (_sse_sgn24)\
                      : \
                      "%xmm0", \
                      "%xmm1", \
                      "%xmm11", \
                      "%xmm12", \
                      "%xmm13", \
                      "%xmm2", \
                      "%xmm3", \
                      "%xmm4", \
                      "%xmm5", \
                      "memory"); \
__asm__ __volatile__ ("movlps %%xmm0, %0 \n\t" \
                      "movhps %%xmm0, %1 \n\t" \
                      : \
                      "=m" ((cc)[1].c[0]), \
                      "=m" ((cc)[1].c[1])\
                      : \
                      : \
                      "%xmm0", \
                      "memory"); \
__asm__ __volatile__ ("movlps %0, %%xmm1 \n\t" \
                      "unpcklps %%xmm1, %%xmm1 \n\t" \
                      "movaps %%xmm1, %%xmm6 \n\t" \
                      "movlps %1, %%xmm2 \n\t" \
                      "unpcklps %%xmm2, %%xmm2 \n\t" \
                      "movlhps %%xmm2, %%xmm6 \n\t" \
                      "movhlps %%xmm1, %%xmm2 \n\t" \
                      "mulps %%xmm14, %%xmm6 \n\t" \
                      "mulps %%xmm14, %%xmm2 \n\t" \
                      "shufps $0xb1, %%xmm2, %%xmm2 \n\t" \
                      "xorps %2, %%xmm2 \n\t" \
                      "addps %%xmm2, %%xmm6 \n\t" \
                      "movlps %3, %%xmm7 \n\t" \
                      "unpcklps %%xmm7, %%xmm7 \n\t" \
                      "mulps %%xmm15, %%xmm7 \n\t" \
                      "shufps $0xb4, %%xmm7, %%xmm7 \n\t" \
                      "xorps %4, %%xmm7 \n\t" \
                      "addps %%xmm7, %%xmm6 \n\t" \
                      "movhlps %%xmm6, %%xmm5 \n\t" \
                      "addps %%xmm5, %%xmm6 \n\t" \
                      : \
                      : \
                      "m" ((aa)[1].e[0][2]), \
                      "m" ((aa)[1].e[1][2]), \
                      "m" (_sse_sgn24), \
                      "m" ((aa)[1].e[2][2]), \
                      "m" (_sse_sgn4)\
                      : \
                      "%xmm1", \
                      "%xmm14", \
                      "%xmm15", \
                      "%xmm2", \
                      "%xmm5", \
                      "%xmm6", \
                      "%xmm7", \
                      "memory"); \
__asm__ __volatile__ ("movlps %%xmm6, %0 \n\t" \
                      "movaps %%xmm8, %%xmm0 \n\t" \
                      : \
                      "=m" ((cc)[1].c[2])\
                      : \
                      : \
                      "%xmm0", \
                      "%xmm6", \
                      "%xmm8", \
                      "memory"); \
__asm__ __volatile__ ("movlps %0, %%xmm3 \n\t" \
                      "movhps %1, %%xmm3 \n\t" \
                      "mulps %%xmm3, %%xmm0 \n\t" \
                      "movaps %%xmm9, %%xmm1 \n\t" \
                      "movlps %2, %%xmm4 \n\t" \
                      "movhps %3, %%xmm4 \n\t" \
                      "mulps %%xmm4, %%xmm1 \n\t" \
                      "movaps %%xmm10, %%xmm2 \n\t" \
                      "movlps %4, %%xmm5 \n\t" \
                      "movhps %5, %%xmm5" \
                      : \
                      : \
                      "m" ((aa)[2].e[0][0]), \
                      "m" ((aa)[2].e[0][1]), \
                      "m" ((aa)[2].e[1][0]), \
                      "m" ((aa)[2].e[1][1]), \
                      "m" ((aa)[2].e[2][0]), \
                      "m" ((aa)[2].e[2][1])\
                      : \
                      "%xmm0", \
                      "%xmm1", \
                      "%xmm10", \
                      "%xmm2", \
                      "%xmm3", \
                      "%xmm4", \
                      "%xmm5", \
                      "%xmm9", \
                      "memory"); \
__asm__ __volatile__ ("mulps %%xmm5, %%xmm2 \n\t" \
                      "addps %%xmm1, %%xmm0 \n\t" \
                      "addps %%xmm2, %%xmm0 \n\t" \
                      "mulps %%xmm11, %%xmm3 \n\t" \
                      "mulps %%xmm12, %%xmm4 \n\t" \
                      "mulps %%xmm13, %%xmm5 \n\t" \
                      "addps %%xmm4, %%xmm3 \n\t" \
                      "addps %%xmm5, %%xmm3 \n\t" \
                      "xorps %0, %%xmm0 \n\t" \
                      "shufps $0xb1, %%xmm3, %%xmm3 \n\t" \
                      "addps %%xmm3, %%xmm0 \n\t" \
                      : \
                      : \
                      "m" (_sse_sgn24)\
                      : \
                      "%xmm0", \
                      "%xmm1", \
                      "%xmm11", \
                      "%xmm12", \
                      "%xmm13", \
                      "%xmm2", \
                      "%xmm3", \
                      "%xmm4", \
                      "%xmm5", \
                      "memory"); \
__asm__ __volatile__ ("movlps %%xmm0, %0 \n\t" \
                      "movhps %%xmm0, %1 \n\t" \
                      : \
                      "=m" ((cc)[2].c[0]), \
                      "=m" ((cc)[2].c[1])\
                      : \
                      : \
                      "%xmm0", \
                      "memory"); \
__asm__ __volatile__ ("movlps %0, %%xmm1 \n\t" \
                      "unpcklps %%xmm1, %%xmm1 \n\t" \
                      "movaps %%xmm1, %%xmm6 \n\t" \
                      "movlps %1, %%xmm2 \n\t" \
                      "unpcklps %%xmm2, %%xmm2 \n\t" \
                      "movlhps %%xmm2, %%xmm6 \n\t" \
                      "movhlps %%xmm1, %%xmm2 \n\t" \
                      "mulps %%xmm14, %%xmm6 \n\t" \
                      "mulps %%xmm14, %%xmm2 \n\t" \
                      "shufps $0xb1, %%xmm2, %%xmm2 \n\t" \
                      "xorps %2, %%xmm2 \n\t" \
                      "addps %%xmm2, %%xmm6 \n\t" \
                      "movlps %3, %%xmm7 \n\t" \
                      "unpcklps %%xmm7, %%xmm7 \n\t" \
                      "mulps %%xmm15, %%xmm7 \n\t" \
                      "shufps $0xb4, %%xmm7, %%xmm7 \n\t" \
                      "xorps %4, %%xmm7 \n\t" \
                      "addps %%xmm7, %%xmm6 \n\t" \
                      "movhlps %%xmm6, %%xmm5 \n\t" \
                      "addps %%xmm5, %%xmm6 \n\t" \
                      : \
                      : \
                      "m" ((aa)[2].e[0][2]), \
                      "m" ((aa)[2].e[1][2]), \
                      "m" (_sse_sgn24), \
                      "m" ((aa)[2].e[2][2]), \
                      "m" (_sse_sgn4)\
                      : \
                      "%xmm1", \
                      "%xmm14", \
                      "%xmm15", \
                      "%xmm2", \
                      "%xmm5", \
                      "%xmm6", \
                      "%xmm7", \
                      "memory"); \
__asm__ __volatile__ ("movlps %%xmm6, %0 \n\t" \
                      "movaps %%xmm8, %%xmm0 \n\t" \
                      : \
                      "=m" ((cc)[2].c[2])\
                      : \
                      : \
                      "%xmm0", \
                      "%xmm6", \
                      "%xmm8", \
                      "memory"); \
__asm__ __volatile__ ("movlps %0, %%xmm3 \n\t" \
                      "movhps %1, %%xmm3 \n\t" \
                      "mulps %%xmm3, %%xmm0 \n\t" \
                      "movaps %%xmm9, %%xmm1 \n\t" \
                      "movlps %2, %%xmm4 \n\t" \
                      "movhps %3, %%xmm4 \n\t" \
                      "mulps %%xmm4, %%xmm1 \n\t" \
                      "movaps %%xmm10, %%xmm2 \n\t" \
                      "movlps %4, %%xmm5 \n\t" \
                      "movhps %5, %%xmm5" \
                      : \
                      : \
                      "m" ((aa)[3].e[0][0]), \
                      "m" ((aa)[3].e[0][1]), \
                      "m" ((aa)[3].e[1][0]), \
                      "m" ((aa)[3].e[1][1]), \
                      "m" ((aa)[3].e[2][0]), \
                      "m" ((aa)[3].e[2][1])\
                      : \
                      "%xmm0", \
                      "%xmm1", \
                      "%xmm10", \
                      "%xmm2", \
                      "%xmm3", \
                      "%xmm4", \
                      "%xmm5", \
                      "%xmm9", \
                      "memory"); \
__asm__ __volatile__ ("mulps %%xmm5, %%xmm2 \n\t" \
                      "addps %%xmm1, %%xmm0 \n\t" \
                      "addps %%xmm2, %%xmm0 \n\t" \
                      "mulps %%xmm11, %%xmm3 \n\t" \
                      "mulps %%xmm12, %%xmm4 \n\t" \
                      "mulps %%xmm13, %%xmm5 \n\t" \
                      "addps %%xmm4, %%xmm3 \n\t" \
                      "addps %%xmm5, %%xmm3 \n\t" \
                      "xorps %0, %%xmm0 \n\t" \
                      "shufps $0xb1, %%xmm3, %%xmm3 \n\t" \
                      "addps %%xmm3, %%xmm0 \n\t" \
                      : \
                      : \
                      "m" (_sse_sgn24)\
                      : \
                      "%xmm0", \
                      "%xmm1", \
                      "%xmm11", \
                      "%xmm12", \
                      "%xmm13", \
                      "%xmm2", \
                      "%xmm3", \
                      "%xmm4", \
                      "%xmm5", \
                      "memory"); \
__asm__ __volatile__ ("movlps %%xmm0, %0 \n\t" \
                      "movhps %%xmm0, %1 \n\t" \
                      : \
                      "=m" ((cc)[3].c[0]), \
                      "=m" ((cc)[3].c[1])\
                      : \
                      : \
                      "%xmm0", \
                      "memory"); \
__asm__ __volatile__ ("movlps %0, %%xmm1 \n\t" \
                      "unpcklps %%xmm1, %%xmm1 \n\t" \
                      "movaps %%xmm1, %%xmm6 \n\t" \
                      "movlps %1, %%xmm2 \n\t" \
                      "unpcklps %%xmm2, %%xmm2 \n\t" \
                      "movlhps %%xmm2, %%xmm6 \n\t" \
                      "movhlps %%xmm1, %%xmm2 \n\t" \
                      "mulps %%xmm14, %%xmm6 \n\t" \
                      "mulps %%xmm14, %%xmm2 \n\t" \
                      "shufps $0xb1, %%xmm2, %%xmm2 \n\t" \
                      "xorps %2, %%xmm2 \n\t" \
                      "addps %%xmm2, %%xmm6 \n\t" \
                      "movlps %3, %%xmm7 \n\t" \
                      "unpcklps %%xmm7, %%xmm7 \n\t" \
                      "mulps %%xmm15, %%xmm7 \n\t" \
                      "shufps $0xb4, %%xmm7, %%xmm7 \n\t" \
                      "xorps %4, %%xmm7 \n\t" \
                      "addps %%xmm7, %%xmm6 \n\t" \
                      "movhlps %%xmm6, %%xmm5 \n\t" \
                      "addps %%xmm5, %%xmm6 \n\t" \
                      : \
                      : \
                      "m" ((aa)[3].e[0][2]), \
                      "m" ((aa)[3].e[1][2]), \
                      "m" (_sse_sgn24), \
                      "m" ((aa)[3].e[2][2]), \
                      "m" (_sse_sgn4)\
                      : \
                      "%xmm1", \
                      "%xmm14", \
                      "%xmm15", \
                      "%xmm2", \
                      "%xmm5", \
                      "%xmm6", \
                      "%xmm7", \
                      "memory"); \
__asm__ __volatile__ ("movlps %%xmm6, %0 \n\t" \
                      : \
                      "=m" ((cc)[3].c[2])\
                      : \
                      : \
                      "%xmm6", \
                      "memory"); \
}
