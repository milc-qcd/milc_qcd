#define _inline_sse_su3_projector(aa,bb,cc) \
{ \
__asm__ __volatile__ ("movlps %0, %%xmm0 \n\t" \
                      "movhps %1, %%xmm0 \n\t" \
                      "movaps %%xmm0, %%xmm1 \n\t" \
                      "shufps $0xb1, %%xmm1, %%xmm1 \n\t" \
                      "xorps %2, %%xmm1 \n\t" \
                      "movss %3, %%xmm2 \n\t" \
                      "shufps $0x00, %%xmm2, %%xmm2 \n\t" \
                      "movss %4, %%xmm3 \n\t" \
                      "shufps $0x00, %%xmm3, %%xmm3 \n\t" \
                      "mulps %%xmm0, %%xmm2 \n\t" \
                      "mulps %%xmm1, %%xmm3 \n\t" \
                      "addps %%xmm3, %%xmm2 \n\t" \
                      "xorps %5, %%xmm2" \
                      : \
                      : \
                      "m" ((bb)->c[0]), \
                      "m" ((bb)->c[1]), \
                      "m" (_sse_sgn24), \
                      "m" ((aa)->c[0].real), \
                      "m" ((aa)->c[0].imag), \
                      "m" (_sse_sgn24)\
                      : \
                      "%xmm0", \
                      "%xmm1", \
                      "%xmm2", \
                      "%xmm3", \
                      "memory"); \
__asm__ __volatile__ ("movups %%xmm2, %0 \n\t" \
                      : \
                      "=m" ((cc)->e[0][0])\
                      : \
                      : \
                      "%xmm2", \
                      "memory"); \
__asm__ __volatile__ ("movss %0, %%xmm2 \n\t" \
                      "shufps $0x00, %%xmm2, %%xmm2 \n\t" \
                      "movss %1, %%xmm3 \n\t" \
                      "shufps $0x00, %%xmm3, %%xmm3 \n\t" \
                      "mulps %%xmm0, %%xmm2 \n\t" \
                      "mulps %%xmm1, %%xmm3 \n\t" \
                      "addps %%xmm3, %%xmm2 \n\t" \
                      "xorps %2, %%xmm2 \n\t" \
                      : \
                      : \
                      "m" ((aa)->c[1].real), \
                      "m" ((aa)->c[1].imag), \
                      "m" (_sse_sgn24)\
                      : \
                      "%xmm0", \
                      "%xmm1", \
                      "%xmm2", \
                      "%xmm3", \
                      "memory"); \
__asm__ __volatile__ ("movups %%xmm2, %0 \n\t" \
                      : \
                      "=m" ((cc)->e[1][0])\
                      : \
                      : \
                      "%xmm2", \
                      "memory"); \
__asm__ __volatile__ ("movss %0, %%xmm2 \n\t" \
                      "shufps $0x00, %%xmm2, %%xmm2 \n\t" \
                      "movss %1, %%xmm3 \n\t" \
                      "shufps $0x00, %%xmm3, %%xmm3 \n\t" \
                      "mulps %%xmm0, %%xmm2 \n\t" \
                      "mulps %%xmm1, %%xmm3 \n\t" \
                      "addps %%xmm3, %%xmm2 \n\t" \
                      "xorps %2, %%xmm2 \n\t" \
                      : \
                      : \
                      "m" ((aa)->c[2].real), \
                      "m" ((aa)->c[2].imag), \
                      "m" (_sse_sgn24)\
                      : \
                      "%xmm0", \
                      "%xmm1", \
                      "%xmm2", \
                      "%xmm3", \
                      "memory"); \
__asm__ __volatile__ ("movups %%xmm2, %0 \n\t" \
                      : \
                      "=m" ((cc)->e[2][0])\
                      : \
                      : \
                      "%xmm2", \
                      "memory"); \
__asm__ __volatile__ ("movlps %0, %%xmm0 \n\t" \
                      "movhps %1, %%xmm0 \n\t" \
                      "movaps %%xmm0, %%xmm1 \n\t" \
                      "shufps $0xb1, %%xmm1, %%xmm1 \n\t" \
                      "xorps %2, %%xmm1 \n\t" \
                      "movss %3, %%xmm2 \n\t" \
                      "shufps $0x00, %%xmm2, %%xmm2 \n\t" \
                      "movss %4, %%xmm3 \n\t" \
                      "shufps $0x00, %%xmm3, %%xmm3 \n\t" \
                      "mulps %%xmm0, %%xmm2 \n\t" \
                      "mulps %%xmm1, %%xmm3 \n\t" \
                      "addps %%xmm3, %%xmm2 \n\t" \
                      : \
                      : \
                      "m" ((aa)->c[0]), \
                      "m" ((aa)->c[1]), \
                      "m" (_sse_sgn24), \
                      "m" ((bb)->c[2].real), \
                      "m" ((bb)->c[2].imag)\
                      : \
                      "%xmm0", \
                      "%xmm1", \
                      "%xmm2", \
                      "%xmm3", \
                      "memory"); \
__asm__ __volatile__ ("movlps %%xmm2, %0 \n\t" \
                      "movhps %%xmm2, %1 \n\t" \
                      : \
                      "=m" ((cc)->e[0][2]), \
                      "=m" ((cc)->e[1][2])\
                      : \
                      : \
                      "%xmm2", \
                      "memory"); \
__asm__ __volatile__ ("movlps %0, %%xmm0 \n\t" \
                      "shufps $0x14, %%xmm0, %%xmm0 \n\t" \
                      "movlps %1, %%xmm2 \n\t" \
                      "shufps $0x44, %%xmm2, %%xmm2 \n\t" \
                      "xorps %2, %%xmm2 \n\t" \
                      "mulps %%xmm0, %%xmm2 \n\t" \
                      "movaps %%xmm2, %%xmm1 \n\t" \
                      "shufps $0xd4, %%xmm1, %%xmm1 \n\t" \
                      "shufps $0x8c, %%xmm2, %%xmm2 \n\t" \
                      "addps %%xmm1, %%xmm2 \n\t" \
                      : \
                      : \
                      "m" ((aa)->c[2]), \
                      "m" ((bb)->c[2]), \
                      "m" (_sse_sgn4)\
                      : \
                      "%xmm0", \
                      "%xmm1", \
                      "%xmm2", \
                      "memory"); \
__asm__ __volatile__ ("movhps %%xmm2, %0 \n\t" \
                      : \
                      "=m" ((cc)->e[2][2])\
                      : \
                      : \
                      "%xmm2", \
                      "memory"); \
}
