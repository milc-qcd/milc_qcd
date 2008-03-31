#define _inline_sse_add_su3_vector(aa,bb,cc) \
{ \
__asm__ __volatile__ ("movups %0, %%xmm0 \n\t" \
                      "movlps %1, %%xmm1 \n\t" \
                      "shufps $0x44, %%xmm1, %%xmm1 \n\t" \
                      "movups %2, %%xmm2 \n\t" \
                      "movlps %3, %%xmm3 \n\t" \
                      "shufps $0x44, %%xmm3, %%xmm3 \n\t" \
                      "addps %%xmm2, %%xmm0 \n\t" \
                      "addps %%xmm3, %%xmm1 \n\t" \
                      : \
                      : \
                      "m" ((aa)->c[0]), \
                      "m" ((aa)->c[2]), \
                      "m" ((bb)->c[0]), \
                      "m" ((bb)->c[2])\
                      : \
                      "%xmm0", \
                      "%xmm1", \
                      "%xmm2", \
                      "%xmm3", \
                      "memory"); \
__asm__ __volatile__ ("movups %%xmm0, %0 \n\t" \
                      "movlps %%xmm1, %1 \n\t" \
                      : \
                      "=m" ((cc)->c[0]), \
                      "=m" ((cc)->c[2])\
                      : \
                      : \
                      "%xmm0", \
                      "%xmm1", \
                      "memory"); \
}
