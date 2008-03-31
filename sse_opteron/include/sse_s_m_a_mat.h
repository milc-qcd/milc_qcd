#define _inline_sse_scalar_mult_add_su3_matrix(aa,bb,cc,dd) \
{ \
__asm__ __volatile__ ("movss %0, %%xmm0 \n\t" \
                      "unpcklps %%xmm0, %%xmm0 \n\t" \
                      "movlhps %%xmm0, %%xmm0 \n\t" \
                      "movlps %1, %%xmm2 \n\t" \
                      "movhps %2, %%xmm2 \n\t" \
                      "mulps %%xmm0, %%xmm2 \n\t" \
                      "movlps %3, %%xmm1 \n\t" \
                      "movhps %4, %%xmm1 \n\t" \
                      "addps %%xmm2, %%xmm1 \n\t" \
                      : \
                      : \
                      "m" (cc), \
                      "m" ((bb)->e[0][0]), \
                      "m" ((bb)->e[0][1]), \
                      "m" ((aa)->e[0][0]), \
                      "m" ((aa)->e[0][1])\
                      : \
                      "%xmm0", \
                      "%xmm1", \
                      "%xmm2", \
                      "memory"); \
__asm__ __volatile__ ("movlps %%xmm1, %0 \n\t" \
                      "movhps %%xmm1, %1 \n\t" \
                      : \
                      "=m" ((dd)->e[0][0]), \
                      "=m" ((dd)->e[0][1])\
                      : \
                      : \
                      "%xmm1", \
                      "memory"); \
__asm__ __volatile__ ("movlps %0, %%xmm4 \n\t" \
                      "movhps %1, %%xmm4 \n\t" \
                      "mulps %%xmm0, %%xmm4 \n\t" \
                      "movlps %2, %%xmm3 \n\t" \
                      "movhps %3, %%xmm3 \n\t" \
                      "addps %%xmm4, %%xmm3 \n\t" \
                      : \
                      : \
                      "m" ((bb)->e[0][2]), \
                      "m" ((bb)->e[1][0]), \
                      "m" ((aa)->e[0][2]), \
                      "m" ((aa)->e[1][0])\
                      : \
                      "%xmm0", \
                      "%xmm3", \
                      "%xmm4", \
                      "memory"); \
__asm__ __volatile__ ("movlps %%xmm3, %0 \n\t" \
                      "movhps %%xmm3, %1 \n\t" \
                      : \
                      "=m" ((dd)->e[0][2]), \
                      "=m" ((dd)->e[1][0])\
                      : \
                      : \
                      "%xmm3", \
                      "memory"); \
__asm__ __volatile__ ("movlps %0, %%xmm6 \n\t" \
                      "movhps %1, %%xmm6 \n\t" \
                      "mulps %%xmm0, %%xmm6 \n\t" \
                      "movlps %2, %%xmm5 \n\t" \
                      "movhps %3, %%xmm5 \n\t" \
                      "addps %%xmm6, %%xmm5 \n\t" \
                      : \
                      : \
                      "m" ((bb)->e[1][1]), \
                      "m" ((bb)->e[1][2]), \
                      "m" ((aa)->e[1][1]), \
                      "m" ((aa)->e[1][2])\
                      : \
                      "%xmm0", \
                      "%xmm5", \
                      "%xmm6", \
                      "memory"); \
__asm__ __volatile__ ("movlps %%xmm5, %0 \n\t" \
                      "movhps %%xmm5, %1 \n\t" \
                      : \
                      "=m" ((dd)->e[1][1]), \
                      "=m" ((dd)->e[1][2])\
                      : \
                      : \
                      "%xmm5", \
                      "memory"); \
__asm__ __volatile__ ("movlps %0, %%xmm8 \n\t" \
                      "movhps %1, %%xmm8 \n\t" \
                      "mulps %%xmm0, %%xmm8 \n\t" \
                      "movlps %2, %%xmm7 \n\t" \
                      "movhps %3, %%xmm7 \n\t" \
                      "addps %%xmm8, %%xmm7 \n\t" \
                      : \
                      : \
                      "m" ((bb)->e[2][0]), \
                      "m" ((bb)->e[2][1]), \
                      "m" ((aa)->e[2][0]), \
                      "m" ((aa)->e[2][1])\
                      : \
                      "%xmm0", \
                      "%xmm7", \
                      "%xmm8", \
                      "memory"); \
__asm__ __volatile__ ("movlps %%xmm7, %0 \n\t" \
                      "movhps %%xmm7, %1 \n\t" \
                      : \
                      "=m" ((dd)->e[2][0]), \
                      "=m" ((dd)->e[2][1])\
                      : \
                      : \
                      "%xmm7", \
                      "memory"); \
__asm__ __volatile__ ("movlps %0, %%xmm10 \n\t" \
                      "mulps %%xmm0, %%xmm10 \n\t" \
                      "movlps %1, %%xmm9 \n\t" \
                      "addps %%xmm10, %%xmm9 \n\t" \
                      : \
                      : \
                      "m" ((bb)->e[2][2]), \
                      "m" ((aa)->e[2][2])\
                      : \
                      "%xmm0", \
                      "%xmm10", \
                      "%xmm9", \
                      "memory"); \
__asm__ __volatile__ ("movlps %%xmm9, %0 \n\t" \
                      : \
                      "=m" ((dd)->e[2][2])\
                      : \
                      : \
                      "%xmm9", \
                      "memory"); \
}
