#define _inline_sse_add_su3_vector(aa,bb,cc) \
{ \
__asm__ __volatile__ ("movlps %0, %%xmm0 \n\t" \
                      "movhps %1, %%xmm0 \n\t" \
                      "movlps %2, %%xmm2 \n\t" \
                      "movhps %3, %%xmm2 \n\t" \
                      "addps %%xmm2, %%xmm0 \n\t" \
                      : \
                      : \
                      "m" ((aa)->c[0]), \
                      "m" ((aa)->c[1]), \
                      "m" ((bb)->c[0]), \
                      "m" ((bb)->c[1])\
                      : \
                      "%xmm0", \
                      "%xmm2", \
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
__asm__ __volatile__ ("movlps %0, %%xmm3 \n\t" \
                      "movlhps %%xmm3, %%xmm3 \n\t" \
                      "movlps %1, %%xmm1 \n\t" \
                      "movlhps %%xmm1, %%xmm1 \n\t" \
                      "addps %%xmm3, %%xmm1 \n\t" \
                      : \
                      : \
                      "m" ((bb)->c[2]), \
                      "m" ((aa)->c[2])\
                      : \
                      "%xmm1", \
                      "%xmm3", \
                      "memory"); \
__asm__ __volatile__ ("movlps %%xmm1, %0 \n\t" \
                      : \
                      "=m" ((cc)->c[2])\
                      : \
                      : \
                      "%xmm1", \
                      "memory"); \
}
