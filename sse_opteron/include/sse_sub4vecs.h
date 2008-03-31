#define _inline_sse_sub_four_su3_vecs(aa,bb0,bb1,bb2,bb3) \
{ \
__asm__ __volatile__ ("movlps %0, %%xmm0 \n\t" \
                      "movhps %1, %%xmm0 \n\t" \
                      "movlps %2, %%xmm2 \n\t" \
                      "movhps %3, %%xmm2 \n\t" \
                      "subps %%xmm2, %%xmm0 \n\t" \
                      "movlps %4, %%xmm1 \n\t" \
                      "movlhps %%xmm1, %%xmm1 \n\t" \
                      "movlps %5, %%xmm3" \
                      : \
                      : \
                      "m" ((aa)->c[0]), \
                      "m" ((aa)->c[1]), \
                      "m" ((bb0)->c[0]), \
                      "m" ((bb0)->c[1]), \
                      "m" ((aa)->c[2]), \
                      "m" ((bb0)->c[2])\
                      : \
                      "%xmm0", \
                      "%xmm1", \
                      "%xmm2", \
                      "%xmm3", \
                      "memory"); \
__asm__ __volatile__ ("movlhps %%xmm3, %%xmm3 \n\t" \
                      "subps %%xmm3, %%xmm1 \n\t" \
                      "movlps %0, %%xmm4 \n\t" \
                      "movhps %1, %%xmm4 \n\t" \
                      "subps %%xmm4, %%xmm0 \n\t" \
                      "movlps %2, %%xmm5 \n\t" \
                      "movlhps %%xmm5, %%xmm5 \n\t" \
                      "subps %%xmm5, %%xmm1 \n\t" \
                      "movlps %3, %%xmm6 \n\t" \
                      "movhps %4, %%xmm6 \n\t" \
                      "subps %%xmm6, %%xmm0 \n\t" \
                      "movlps %5, %%xmm7" \
                      : \
                      : \
                      "m" ((bb1)->c[0]), \
                      "m" ((bb1)->c[1]), \
                      "m" ((bb1)->c[2]), \
                      "m" ((bb2)->c[0]), \
                      "m" ((bb2)->c[1]), \
                      "m" ((bb2)->c[2])\
                      : \
                      "%xmm0", \
                      "%xmm1", \
                      "%xmm3", \
                      "%xmm4", \
                      "%xmm5", \
                      "%xmm6", \
                      "%xmm7", \
                      "memory"); \
__asm__ __volatile__ ("movlhps %%xmm7, %%xmm7 \n\t" \
                      "subps %%xmm7, %%xmm1 \n\t" \
                      "movlps %0, %%xmm8 \n\t" \
                      "movhps %1, %%xmm8 \n\t" \
                      "subps %%xmm8, %%xmm0 \n\t" \
                      "movlps %2, %%xmm9 \n\t" \
                      "movlhps %%xmm9, %%xmm9 \n\t" \
                      "subps %%xmm9, %%xmm1 \n\t" \
                      : \
                      : \
                      "m" ((bb3)->c[0]), \
                      "m" ((bb3)->c[1]), \
                      "m" ((bb3)->c[2])\
                      : \
                      "%xmm0", \
                      "%xmm1", \
                      "%xmm7", \
                      "%xmm8", \
                      "%xmm9", \
                      "memory"); \
__asm__ __volatile__ ("movlps %%xmm0, %0 \n\t" \
                      "movhps %%xmm0, %1 \n\t" \
                      "movlps %%xmm1, %2 \n\t" \
                      : \
                      "=m" ((aa)->c[0]), \
                      "=m" ((aa)->c[1]), \
                      "=m" ((aa)->c[2])\
                      : \
                      : \
                      "%xmm0", \
                      "%xmm1", \
                      "memory"); \
}
