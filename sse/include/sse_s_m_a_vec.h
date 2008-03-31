#define _inline_sse_scalar_mult_add_su3_vector(aa,bb,cc,dd) \
{ \
__asm__ __volatile__ ("movss %0, %%xmm4 \n\t" \
                      "shufps $0x00, %%xmm4, %%xmm4 \n\t" \
                      "movups %1, %%xmm0 \n\t" \
                      "movlps %2, %%xmm1 \n\t" \
                      "shufps $0x44, %%xmm1, %%xmm1 \n\t" \
                      "movups %3, %%xmm2 \n\t" \
                      "movlps %4, %%xmm3 \n\t" \
                      "shufps $0x44, %%xmm3, %%xmm3 \n\t" \
                      "mulps %%xmm4, %%xmm2 \n\t" \
                      "mulps %%xmm4, %%xmm3 \n\t" \
                      "addps %%xmm2, %%xmm0 \n\t" \
                      "addps %%xmm3, %%xmm1 \n\t" \
                      : \
                      : \
                      "m" ((cc)), \
                      "m" ((aa)->c[0]), \
                      "m" ((aa)->c[2]), \
                      "m" ((bb)->c[0]), \
                      "m" ((bb)->c[2])\
                      : \
                      "%xmm0", \
                      "%xmm1", \
                      "%xmm2", \
                      "%xmm3", \
                      "%xmm4", \
                      "memory"); \
__asm__ __volatile__ ("movups %%xmm0, %0 \n\t" \
                      "movlps %%xmm1, %1 \n\t" \
                      : \
                      "=m" ((dd)->c[0]), \
                      "=m" ((dd)->c[2])\
                      : \
                      : \
                      "%xmm0", \
                      "%xmm1", \
                      "memory"); \
}
