#define _inline_sse_su3_projector(aa,bb,cc) \
{ \
__asm__ __volatile__ ("movlps %0, %%xmm0 \n\t" \
                      "movhps %1, %%xmm0 \n\t" \
                      "movaps %%xmm0, %%xmm1 \n\t" \
                      "shufps $0xb1, %%xmm1, %%xmm1 \n\t" \
                      "xorps %2, %%xmm1 \n\t" \
                      "movlps %3, %%xmm2 \n\t" \
                      "movaps %%xmm2, %%xmm10 \n\t" \
                      "unpcklps %%xmm2, %%xmm2 \n\t" \
                      "movhlps %%xmm2, %%xmm3 \n\t" \
                      "movlhps %%xmm2, %%xmm2 \n\t" \
                      "movlhps %%xmm3, %%xmm3 \n\t" \
                      "movlps %4, %%xmm4 \n\t" \
                      "movlhps %%xmm4, %%xmm10 \n\t" \
                      "unpcklps %%xmm4, %%xmm4 \n\t" \
                      "movhlps %%xmm4, %%xmm5 \n\t" \
                      "movlhps %%xmm4, %%xmm4 \n\t" \
                      "movlhps %%xmm5, %%xmm5 \n\t" \
                      "movlps %5, %%xmm6" \
                      : \
                      : \
                      "m" ((bb)->c[0]), \
                      "m" ((bb)->c[1]), \
                      "m" (_sse_sgn24), \
                      "m" ((aa)->c[0]), \
                      "m" ((aa)->c[1]), \
                      "m" ((aa)->c[2])\
                      : \
                      "%xmm0", \
                      "%xmm1", \
                      "%xmm10", \
                      "%xmm2", \
                      "%xmm3", \
                      "%xmm4", \
                      "%xmm5", \
                      "%xmm6", \
                      "memory"); \
__asm__ __volatile__ ("unpcklps %%xmm6, %%xmm6 \n\t" \
                      "movaps %%xmm6, %%xmm8 \n\t" \
                      "movhlps %%xmm6, %%xmm7 \n\t" \
                      "movlhps %%xmm6, %%xmm6 \n\t" \
                      "movlhps %%xmm7, %%xmm7 \n\t" \
                      "mulps %%xmm0, %%xmm2 \n\t" \
                      "mulps %%xmm0, %%xmm4 \n\t" \
                      "mulps %%xmm0, %%xmm6 \n\t" \
                      "mulps %%xmm1, %%xmm3 \n\t" \
                      "mulps %%xmm1, %%xmm5 \n\t" \
                      "mulps %%xmm1, %%xmm7 \n\t" \
                      "addps %%xmm3, %%xmm2 \n\t" \
                      "xorps %0, %%xmm2 \n\t" \
                      : \
                      : \
                      "m" (_sse_sgn24)\
                      : \
                      "%xmm0", \
                      "%xmm1", \
                      "%xmm2", \
                      "%xmm3", \
                      "%xmm4", \
                      "%xmm5", \
                      "%xmm6", \
                      "%xmm7", \
                      "%xmm8", \
                      "memory"); \
__asm__ __volatile__ ("movlps %%xmm2, %0 \n\t" \
                      "movhps %%xmm2, %1 \n\t" \
                      "addps %%xmm5, %%xmm4 \n\t" \
                      : \
                      "=m" ((cc)->e[0][0]), \
                      "=m" ((cc)->e[0][1])\
                      : \
                      : \
                      "%xmm2", \
                      "%xmm4", \
                      "%xmm5", \
                      "memory"); \
__asm__ __volatile__ ("xorps %0, %%xmm4 \n\t" \
                      : \
                      : \
                      "m" (_sse_sgn24)\
                      : \
                      "%xmm4", \
                      "memory"); \
__asm__ __volatile__ ("movlps %%xmm4, %0 \n\t" \
                      "movhps %%xmm4, %1 \n\t" \
                      "addps %%xmm7, %%xmm6 \n\t" \
                      : \
                      "=m" ((cc)->e[1][0]), \
                      "=m" ((cc)->e[1][1])\
                      : \
                      : \
                      "%xmm4", \
                      "%xmm6", \
                      "%xmm7", \
                      "memory"); \
__asm__ __volatile__ ("xorps %0, %%xmm6 \n\t" \
                      : \
                      : \
                      "m" (_sse_sgn24)\
                      : \
                      "%xmm6", \
                      "memory"); \
__asm__ __volatile__ ("movlps %%xmm6, %0 \n\t" \
                      "movhps %%xmm6, %1 \n\t" \
                      "movaps %%xmm10, %%xmm1 \n\t" \
                      "shufps $0xb1, %%xmm1, %%xmm1 \n\t" \
                      : \
                      "=m" ((cc)->e[2][0]), \
                      "=m" ((cc)->e[2][1])\
                      : \
                      : \
                      "%xmm1", \
                      "%xmm10", \
                      "%xmm6", \
                      "memory"); \
__asm__ __volatile__ ("xorps %0, %%xmm1 \n\t" \
                      "movlps %1, %%xmm2 \n\t" \
                      "movaps %%xmm2, %%xmm11 \n\t" \
                      "unpcklps %%xmm2, %%xmm2 \n\t" \
                      "movhlps %%xmm2, %%xmm3 \n\t" \
                      "movlhps %%xmm2, %%xmm2 \n\t" \
                      "movlhps %%xmm3, %%xmm3 \n\t" \
                      "mulps %%xmm10, %%xmm2 \n\t" \
                      "mulps %%xmm1, %%xmm3 \n\t" \
                      "addps %%xmm3, %%xmm2 \n\t" \
                      : \
                      : \
                      "m" (_sse_sgn24), \
                      "m" ((bb)->c[2])\
                      : \
                      "%xmm1", \
                      "%xmm10", \
                      "%xmm11", \
                      "%xmm2", \
                      "%xmm3", \
                      "memory"); \
__asm__ __volatile__ ("movlps %%xmm2, %0 \n\t" \
                      "movhps %%xmm2, %1 \n\t" \
                      "movlhps %%xmm11, %%xmm11 \n\t" \
                      "mulps %%xmm11, %%xmm8 \n\t" \
                      : \
                      "=m" ((cc)->e[0][2]), \
                      "=m" ((cc)->e[1][2])\
                      : \
                      : \
                      "%xmm11", \
                      "%xmm2", \
                      "%xmm8", \
                      "memory"); \
__asm__ __volatile__ ("xorps %0, %%xmm8 \n\t" \
                      "shufps $0xb4, %%xmm8, %%xmm8 \n\t" \
                      "movhlps %%xmm8, %%xmm5 \n\t" \
                      "addps %%xmm5, %%xmm8 \n\t" \
                      : \
                      : \
                      "m" (_sse_sgn2)\
                      : \
                      "%xmm5", \
                      "%xmm8", \
                      "memory"); \
__asm__ __volatile__ ("movlps %%xmm8, %0 \n\t" \
                      : \
                      "=m" ((cc)->e[2][2])\
                      : \
                      : \
                      "%xmm8", \
                      "memory"); \
}
