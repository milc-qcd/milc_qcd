#define _inline_sse_su3_projector(aa,bb,cc) \
{ \
__asm__ __volatile__ ("movsd %0, %%xmm3 \n\t" \
                      "unpcklpd %%xmm3, %%xmm3 \n\t" \
                      "movsd %1, %%xmm4 \n\t" \
                      "unpcklpd %%xmm4, %%xmm4 \n\t" \
                      "movsd %2, %%xmm5 \n\t" \
                      "unpcklpd %%xmm5, %%xmm5 \n\t" \
                      "movsd %3, %%xmm6 \n\t" \
                      "unpcklpd %%xmm6, %%xmm6 \n\t" \
                      "movsd %4, %%xmm7 \n\t" \
                      "unpcklpd %%xmm7, %%xmm7 \n\t" \
                      "movupd %5, %%xmm0" \
                      : \
                      : \
                      "m" ((aa)->c[0].real), \
                      "m" ((aa)->c[0].imag), \
                      "m" ((aa)->c[1].real), \
                      "m" ((aa)->c[1].imag), \
                      "m" ((aa)->c[2].real), \
                      "m" ((bb)->c[0])\
                      : \
                      "%xmm0", \
                      "%xmm3", \
                      "%xmm4", \
                      "%xmm5", \
                      "%xmm6", \
                      "%xmm7", \
                      "memory"); \
__asm__ __volatile__ ("movupd %%xmm3, %%xmm1 \n\t" \
                      "movupd %%xmm4, %%xmm2 \n\t" \
                      "mulpd %%xmm0, %%xmm1 \n\t" \
                      "mulpd %%xmm0, %%xmm2 \n\t" \
                      "xorpd %0, %%xmm1 \n\t" \
                      "shufpd $0x1, %%xmm2, %%xmm2 \n\t" \
                      "addpd %%xmm1, %%xmm2 \n\t" \
                      : \
                      : \
                      "m" (_sse_sgn4)\
                      : \
                      "%xmm0", \
                      "%xmm1", \
                      "%xmm2", \
                      "%xmm3", \
                      "%xmm4", \
                      "memory"); \
__asm__ __volatile__ ("movupd %%xmm2, %0 \n\t" \
                      "movupd %%xmm5, %%xmm1 \n\t" \
                      "movupd %%xmm6, %%xmm2 \n\t" \
                      "mulpd %%xmm0, %%xmm2 \n\t" \
                      "mulpd %%xmm0, %%xmm1 \n\t" \
                      : \
                      "=m" ((cc)->e[0][0])\
                      : \
                      : \
                      "%xmm0", \
                      "%xmm1", \
                      "%xmm2", \
                      "%xmm5", \
                      "%xmm6", \
                      "memory"); \
__asm__ __volatile__ ("xorpd %0, %%xmm1 \n\t" \
                      "shufpd $0x1, %%xmm2, %%xmm2 \n\t" \
                      "addpd %%xmm1, %%xmm2 \n\t" \
                      : \
                      : \
                      "m" (_sse_sgn4)\
                      : \
                      "%xmm1", \
                      "%xmm2", \
                      "memory"); \
__asm__ __volatile__ ("movupd %%xmm2, %0 \n\t" \
                      "movupd %%xmm7, %%xmm1 \n\t" \
                      : \
                      "=m" ((cc)->e[1][0])\
                      : \
                      : \
                      "%xmm1", \
                      "%xmm2", \
                      "%xmm7", \
                      "memory"); \
__asm__ __volatile__ ("movsd %0, %%xmm2 \n\t" \
                      "unpcklpd %%xmm2, %%xmm2 \n\t" \
                      "mulpd %%xmm0, %%xmm2 \n\t" \
                      "mulpd %%xmm0, %%xmm1 \n\t" \
                      "xorpd %1, %%xmm1 \n\t" \
                      "shufpd $0x1, %%xmm2, %%xmm2 \n\t" \
                      "addpd %%xmm1, %%xmm2 \n\t" \
                      : \
                      : \
                      "m" ((aa)->c[2].imag), \
                      "m" (_sse_sgn4)\
                      : \
                      "%xmm0", \
                      "%xmm1", \
                      "%xmm2", \
                      "memory"); \
__asm__ __volatile__ ("movupd %%xmm2, %0 \n\t" \
                      : \
                      "=m" ((cc)->e[2][0])\
                      : \
                      : \
                      "%xmm2", \
                      "memory"); \
__asm__ __volatile__ ("movupd %0, %%xmm0 \n\t" \
                      "movupd %%xmm3, %%xmm1 \n\t" \
                      "movupd %%xmm4, %%xmm2 \n\t" \
                      "mulpd %%xmm0, %%xmm1 \n\t" \
                      "mulpd %%xmm0, %%xmm2 \n\t" \
                      "xorpd %1, %%xmm1 \n\t" \
                      "shufpd $0x1, %%xmm2, %%xmm2 \n\t" \
                      "addpd %%xmm1, %%xmm2 \n\t" \
                      : \
                      : \
                      "m" ((bb)->c[1]), \
                      "m" (_sse_sgn4)\
                      : \
                      "%xmm0", \
                      "%xmm1", \
                      "%xmm2", \
                      "%xmm3", \
                      "%xmm4", \
                      "memory"); \
__asm__ __volatile__ ("movupd %%xmm2, %0 \n\t" \
                      "movupd %%xmm5, %%xmm1 \n\t" \
                      "movupd %%xmm6, %%xmm2 \n\t" \
                      "mulpd %%xmm0, %%xmm2 \n\t" \
                      "mulpd %%xmm0, %%xmm1 \n\t" \
                      : \
                      "=m" ((cc)->e[0][1])\
                      : \
                      : \
                      "%xmm0", \
                      "%xmm1", \
                      "%xmm2", \
                      "%xmm5", \
                      "%xmm6", \
                      "memory"); \
__asm__ __volatile__ ("xorpd %0, %%xmm1 \n\t" \
                      "shufpd $0x1, %%xmm2, %%xmm2 \n\t" \
                      "addpd %%xmm1, %%xmm2 \n\t" \
                      : \
                      : \
                      "m" (_sse_sgn4)\
                      : \
                      "%xmm1", \
                      "%xmm2", \
                      "memory"); \
__asm__ __volatile__ ("movupd %%xmm2, %0 \n\t" \
                      "movupd %%xmm7, %%xmm1 \n\t" \
                      : \
                      "=m" ((cc)->e[1][1])\
                      : \
                      : \
                      "%xmm1", \
                      "%xmm2", \
                      "%xmm7", \
                      "memory"); \
__asm__ __volatile__ ("movsd %0, %%xmm2 \n\t" \
                      "unpcklpd %%xmm2, %%xmm2 \n\t" \
                      "mulpd %%xmm0, %%xmm2 \n\t" \
                      "mulpd %%xmm0, %%xmm1 \n\t" \
                      "xorpd %1, %%xmm1 \n\t" \
                      "shufpd $0x1, %%xmm2, %%xmm2 \n\t" \
                      "addpd %%xmm1, %%xmm2 \n\t" \
                      : \
                      : \
                      "m" ((aa)->c[2].imag), \
                      "m" (_sse_sgn4)\
                      : \
                      "%xmm0", \
                      "%xmm1", \
                      "%xmm2", \
                      "memory"); \
__asm__ __volatile__ ("movupd %%xmm2, %0 \n\t" \
                      : \
                      "=m" ((cc)->e[2][1])\
                      : \
                      : \
                      "%xmm2", \
                      "memory"); \
__asm__ __volatile__ ("movupd %0, %%xmm0 \n\t" \
                      "movupd %%xmm3, %%xmm1 \n\t" \
                      "movupd %%xmm4, %%xmm2 \n\t" \
                      "mulpd %%xmm0, %%xmm1 \n\t" \
                      "mulpd %%xmm0, %%xmm2 \n\t" \
                      "xorpd %1, %%xmm1 \n\t" \
                      "shufpd $0x1, %%xmm2, %%xmm2 \n\t" \
                      "addpd %%xmm1, %%xmm2 \n\t" \
                      : \
                      : \
                      "m" ((bb)->c[2]), \
                      "m" (_sse_sgn4)\
                      : \
                      "%xmm0", \
                      "%xmm1", \
                      "%xmm2", \
                      "%xmm3", \
                      "%xmm4", \
                      "memory"); \
__asm__ __volatile__ ("movupd %%xmm2, %0 \n\t" \
                      "movupd %%xmm5, %%xmm1 \n\t" \
                      "movupd %%xmm6, %%xmm2 \n\t" \
                      "mulpd %%xmm0, %%xmm1 \n\t" \
                      "mulpd %%xmm0, %%xmm2 \n\t" \
                      : \
                      "=m" ((cc)->e[0][2])\
                      : \
                      : \
                      "%xmm0", \
                      "%xmm1", \
                      "%xmm2", \
                      "%xmm5", \
                      "%xmm6", \
                      "memory"); \
__asm__ __volatile__ ("xorpd %0, %%xmm1 \n\t" \
                      "shufpd $0x1, %%xmm2, %%xmm2 \n\t" \
                      "addpd %%xmm1, %%xmm2 \n\t" \
                      : \
                      : \
                      "m" (_sse_sgn4)\
                      : \
                      "%xmm1", \
                      "%xmm2", \
                      "memory"); \
__asm__ __volatile__ ("movupd %%xmm2, %0 \n\t" \
                      "movupd %%xmm7, %%xmm1 \n\t" \
                      : \
                      "=m" ((cc)->e[1][2])\
                      : \
                      : \
                      "%xmm1", \
                      "%xmm2", \
                      "%xmm7", \
                      "memory"); \
__asm__ __volatile__ ("movsd %0, %%xmm2 \n\t" \
                      "unpcklpd %%xmm2, %%xmm2 \n\t" \
                      "mulpd %%xmm0, %%xmm1 \n\t" \
                      "mulpd %%xmm0, %%xmm2 \n\t" \
                      "xorpd %1, %%xmm1 \n\t" \
                      "shufpd $0x1, %%xmm2, %%xmm2 \n\t" \
                      "addpd %%xmm1, %%xmm2 \n\t" \
                      : \
                      : \
                      "m" ((aa)->c[2].imag), \
                      "m" (_sse_sgn4)\
                      : \
                      "%xmm0", \
                      "%xmm1", \
                      "%xmm2", \
                      "memory"); \
__asm__ __volatile__ ("movupd %%xmm2, %0 \n\t" \
                      : \
                      "=m" ((cc)->e[2][2])\
                      : \
                      : \
                      "%xmm2", \
                      "memory"); \
}
