#define _inline_sse_sub_four_su3_vecs(aa,bb0,bb1,bb2,bb3) \
{ \
__asm__ __volatile__ ("movups %0, %%xmm0 \n\t" \
                      "movlps %1, %%xmm1 \n\t" \
                      "shufps $0x44, %%xmm1, %%xmm1 \n\t" \
                      "movups %2, %%xmm2 \n\t" \
                      "movlps %3, %%xmm3 \n\t" \
                      "shufps $0x44, %%xmm3, %%xmm3 \n\t" \
                      "subps %%xmm2, %%xmm0 \n\t" \
                      "subps %%xmm3, %%xmm1 \n\t" \
                      "movups %4, %%xmm2 \n\t" \
                      "movlps %5, %%xmm3" \
                      : \
                      : \
                      "m" ((aa)->c[0]), \
                      "m" ((aa)->c[2]), \
                      "m" ((bb0)->c[0]), \
                      "m" ((bb0)->c[2]), \
                      "m" ((bb1)->c[0]), \
                      "m" ((bb1)->c[2])\
                      : \
                      "%xmm0", \
                      "%xmm1", \
                      "%xmm2", \
                      "%xmm3", \
                      "memory"); \
__asm__ __volatile__ ("shufps $0x44, %%xmm3, %%xmm3 \n\t" \
                      "subps %%xmm2, %%xmm0 \n\t" \
                      "subps %%xmm3, %%xmm1 \n\t" \
                      "movups %0, %%xmm2 \n\t" \
                      "movlps %1, %%xmm3 \n\t" \
                      "shufps $0x44, %%xmm3, %%xmm3 \n\t" \
                      "subps %%xmm2, %%xmm0 \n\t" \
                      "subps %%xmm3, %%xmm1 \n\t" \
                      "movups %2, %%xmm2 \n\t" \
                      "movlps %3, %%xmm3 \n\t" \
                      "shufps $0x44, %%xmm3, %%xmm3 \n\t" \
                      "subps %%xmm2, %%xmm0 \n\t" \
                      "subps %%xmm3, %%xmm1 \n\t" \
                      : \
                      : \
                      "m" ((bb2)->c[0]), \
                      "m" ((bb2)->c[2]), \
                      "m" ((bb3)->c[0]), \
                      "m" ((bb3)->c[2])\
                      : \
                      "%xmm0", \
                      "%xmm1", \
                      "%xmm2", \
                      "%xmm3", \
                      "memory"); \
__asm__ __volatile__ ("movups %%xmm0, %0 \n\t" \
                      "movlps %%xmm1, %1 \n\t" \
                      : \
                      "=m" ((aa)->c[0]), \
                      "=m" ((aa)->c[2])\
                      : \
                      : \
                      "%xmm0", \
                      "%xmm1", \
                      "memory"); \
}
