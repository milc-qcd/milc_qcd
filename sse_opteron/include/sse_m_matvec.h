#define _inline_sse_mult_su3_mat_vec(aa,bb,cc) \
{ \
__asm__ __volatile__ ("movlps %0, %%xmm15 \n\t" \
                      "movaps %%xmm15, %%xmm2 \n\t" \
                      "unpcklps %%xmm2, %%xmm2 \n\t" \
                      "movhlps %%xmm2, %%xmm6 \n\t" \
                      "movlhps %%xmm2, %%xmm2 \n\t" \
                      "movhps %1, %%xmm15 \n\t" \
                      "movaps %%xmm15, %%xmm3 \n\t" \
                      "unpckhps %%xmm3, %%xmm3 \n\t" \
                      "movhlps %%xmm3, %%xmm7 \n\t" \
                      "movlhps %%xmm3, %%xmm3 \n\t" \
                      "movhps %2, %%xmm11 \n\t" \
                      "movaps %%xmm11, %%xmm4 \n\t" \
                      "unpckhps %%xmm4, %%xmm4 \n\t" \
                      "movhlps %%xmm4, %%xmm9 \n\t" \
                      "movlhps %%xmm4, %%xmm4 \n\t" \
                      "movlps %3, %%xmm5 \n\t" \
                      "movhps %4, %%xmm5 \n\t" \
                      "mulps %%xmm5, %%xmm2 \n\t" \
                      "movlps %5, %%xmm8" \
                      : \
                      : \
                      "m" ((bb)->c[0]), \
                      "m" ((bb)->c[1]), \
                      "m" ((bb)->c[2]), \
                      "m" ((aa)->e[0][0]), \
                      "m" ((aa)->e[1][0]), \
                      "m" ((aa)->e[0][1])\
                      : \
                      "%xmm11", \
                      "%xmm15", \
                      "%xmm2", \
                      "%xmm3", \
                      "%xmm4", \
                      "%xmm5", \
                      "%xmm6", \
                      "%xmm7", \
                      "%xmm8", \
                      "%xmm9", \
                      "memory"); \
__asm__ __volatile__ ("movhps %0, %%xmm8 \n\t" \
                      "mulps %%xmm8, %%xmm3 \n\t" \
                      "movlps %1, %%xmm10 \n\t" \
                      "movhps %2, %%xmm10 \n\t" \
                      "mulps %%xmm10, %%xmm4 \n\t" \
                      "addps %%xmm3, %%xmm2 \n\t" \
                      "addps %%xmm4, %%xmm2 \n\t" \
                      "movlhps %%xmm6, %%xmm6 \n\t" \
                      "movlhps %%xmm7, %%xmm7 \n\t" \
                      "movlhps %%xmm9, %%xmm9 \n\t" \
                      "mulps %%xmm5, %%xmm6 \n\t" \
                      "mulps %%xmm8, %%xmm7 \n\t" \
                      "mulps %%xmm10, %%xmm9 \n\t" \
                      "addps %%xmm7, %%xmm6 \n\t" \
                      "addps %%xmm9, %%xmm6 \n\t" \
                      "xorps %3, %%xmm6 \n\t" \
                      "shufps $0xb1, %%xmm6, %%xmm6 \n\t" \
                      "addps %%xmm6, %%xmm2 \n\t" \
                      : \
                      : \
                      "m" ((aa)->e[1][1]), \
                      "m" ((aa)->e[0][2]), \
                      "m" ((aa)->e[1][2]), \
                      "m" (_sse_sgn24)\
                      : \
                      "%xmm10", \
                      "%xmm2", \
                      "%xmm3", \
                      "%xmm4", \
                      "%xmm5", \
                      "%xmm6", \
                      "%xmm7", \
                      "%xmm8", \
                      "%xmm9", \
                      "memory"); \
__asm__ __volatile__ ("movlps %%xmm2, %0 \n\t" \
                      "movhps %%xmm2, %1 \n\t" \
                      : \
                      "=m" ((cc)->c[0]), \
                      "=m" ((cc)->c[1])\
                      : \
                      : \
                      "%xmm2", \
                      "memory"); \
__asm__ __volatile__ ("movlps %0, %%xmm0 \n\t" \
                      "unpcklps %%xmm0, %%xmm0 \n\t" \
                      "movaps %%xmm0, %%xmm12 \n\t" \
                      "movlps %1, %%xmm1 \n\t" \
                      "unpcklps %%xmm1, %%xmm1 \n\t" \
                      "movlhps %%xmm1, %%xmm12 \n\t" \
                      "movaps %%xmm1, %%xmm13 \n\t" \
                      "movhlps %%xmm0, %%xmm13 \n\t" \
                      "mulps %%xmm15, %%xmm12 \n\t" \
                      "mulps %%xmm15, %%xmm13 \n\t" \
                      "xorps %2, %%xmm13 \n\t" \
                      "shufps $0xb1, %%xmm13, %%xmm13 \n\t" \
                      "addps %%xmm13, %%xmm12 \n\t" \
                      "movlps %3, %%xmm11 \n\t" \
                      "movaps %%xmm11, %%xmm8 \n\t" \
                      "unpcklps %%xmm8, %%xmm8 \n\t" \
                      "movhlps %%xmm11, %%xmm11 \n\t" \
                      "mulps %%xmm8, %%xmm11 \n\t" \
                      "xorps %4, %%xmm11 \n\t" \
                      "shufps $0xb4, %%xmm11, %%xmm11 \n\t" \
                      "addps %%xmm11, %%xmm12 \n\t" \
                      "movhlps %%xmm12, %%xmm5 \n\t" \
                      "addps %%xmm5, %%xmm12 \n\t" \
                      : \
                      : \
                      "m" ((aa)->e[2][0]), \
                      "m" ((aa)->e[2][1]), \
                      "m" (_sse_sgn24), \
                      "m" ((aa)->e[2][2]), \
                      "m" (_sse_sgn4)\
                      : \
                      "%xmm0", \
                      "%xmm1", \
                      "%xmm11", \
                      "%xmm12", \
                      "%xmm13", \
                      "%xmm15", \
                      "%xmm5", \
                      "%xmm8", \
                      "memory"); \
__asm__ __volatile__ ("movlps %%xmm12, %0 \n\t" \
                      : \
                      "=m" ((cc)->c[2])\
                      : \
                      : \
                      "%xmm12", \
                      "memory"); \
}
