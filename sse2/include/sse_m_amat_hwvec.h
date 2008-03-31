#define _inline_sse_mult_adj_su3_mat_hwvec(aa,bb,cc) \
{ \
__asm__ __volatile__ ("movupd %0, %%xmm0 \n\t" \
                      "movupd %1, %%xmm1 \n\t" \
                      "movupd %2, %%xmm2 \n\t" \
                      "movsd %3, %%xmm3 \n\t" \
                      "movsd %4, %%xmm6 \n\t" \
                      "movsd %5, %%xmm4" \
                      : \
                      : \
                      "m" ((bb)->h[0].c[0]), \
                      "m" ((bb)->h[0].c[1]), \
                      "m" ((bb)->h[0].c[2]), \
                      "m" ((aa)->e[0][0].real), \
                      "m" ((aa)->e[1][0].real), \
                      "m" ((aa)->e[0][1].real)\
                      : \
                      "%xmm0", \
                      "%xmm1", \
                      "%xmm2", \
                      "%xmm3", \
                      "%xmm4", \
                      "%xmm6", \
                      "memory"); \
__asm__ __volatile__ ("movsd %0, %%xmm7 \n\t" \
                      "movsd %1, %%xmm5 \n\t" \
                      "unpcklpd %%xmm3, %%xmm3 \n\t" \
                      "unpcklpd %%xmm6, %%xmm6 \n\t" \
                      "unpcklpd %%xmm4, %%xmm4 \n\t" \
                      "mulpd %%xmm0, %%xmm3 \n\t" \
                      "unpcklpd %%xmm7, %%xmm7 \n\t" \
                      "mulpd %%xmm1, %%xmm6 \n\t" \
                      "unpcklpd %%xmm5, %%xmm5 \n\t" \
                      "mulpd %%xmm0, %%xmm4 \n\t" \
                      "addpd %%xmm6, %%xmm3 \n\t" \
                      "mulpd %%xmm2, %%xmm7 \n\t" \
                      "mulpd %%xmm0, %%xmm5 \n\t" \
                      "addpd %%xmm7, %%xmm4 \n\t" \
                      "movsd %2, %%xmm6 \n\t" \
                      "movsd %3, %%xmm7 \n\t" \
                      "unpcklpd %%xmm6, %%xmm6 \n\t" \
                      "unpcklpd %%xmm7, %%xmm7 \n\t" \
                      "mulpd %%xmm1, %%xmm6 \n\t" \
                      "mulpd %%xmm2, %%xmm7 \n\t" \
                      "addpd %%xmm6, %%xmm5 \n\t" \
                      "addpd %%xmm7, %%xmm3 \n\t" \
                      "movsd %4, %%xmm6 \n\t" \
                      "movsd %5, %%xmm7" \
                      : \
                      : \
                      "m" ((aa)->e[2][1].real), \
                      "m" ((aa)->e[0][2].real), \
                      "m" ((aa)->e[1][2].real), \
                      "m" ((aa)->e[2][0].real), \
                      "m" ((aa)->e[1][1].real), \
                      "m" ((aa)->e[2][2].real)\
                      : \
                      "%xmm0", \
                      "%xmm1", \
                      "%xmm2", \
                      "%xmm3", \
                      "%xmm4", \
                      "%xmm5", \
                      "%xmm6", \
                      "%xmm7", \
                      "memory"); \
__asm__ __volatile__ ("unpcklpd %%xmm6, %%xmm6 \n\t" \
                      "unpcklpd %%xmm7, %%xmm7 \n\t" \
                      "mulpd %%xmm1, %%xmm6 \n\t" \
                      "mulpd %%xmm2, %%xmm7 \n\t" \
                      "addpd %%xmm6, %%xmm4 \n\t" \
                      "addpd %%xmm7, %%xmm5 \n\t" \
                      "movsd %0, %%xmm6 \n\t" \
                      "movsd %1, %%xmm7 \n\t" \
                      "shufpd $0x1, %%xmm0, %%xmm0 \n\t" \
                      "shufpd $0x1, %%xmm1, %%xmm1 \n\t" \
                      "shufpd $0x1, %%xmm2, %%xmm2 \n\t" \
                      "unpcklpd %%xmm6, %%xmm6 \n\t" \
                      "unpcklpd %%xmm7, %%xmm7 \n\t" \
                      "xorpd %2, %%xmm0 \n\t" \
                      "xorpd %3, %%xmm1 \n\t" \
                      "xorpd %4, %%xmm2 \n\t" \
                      "mulpd %%xmm0, %%xmm6 \n\t" \
                      "mulpd %%xmm1, %%xmm7 \n\t" \
                      "addpd %%xmm6, %%xmm3 \n\t" \
                      "addpd %%xmm7, %%xmm4 \n\t" \
                      "movsd %5, %%xmm6" \
                      : \
                      : \
                      "m" ((aa)->e[0][0].imag), \
                      "m" ((aa)->e[1][1].imag), \
                      "m" (_sse_sgn4), \
                      "m" (_sse_sgn4), \
                      "m" (_sse_sgn4), \
                      "m" ((aa)->e[2][2].imag)\
                      : \
                      "%xmm0", \
                      "%xmm1", \
                      "%xmm2", \
                      "%xmm3", \
                      "%xmm4", \
                      "%xmm5", \
                      "%xmm6", \
                      "%xmm7", \
                      "memory"); \
__asm__ __volatile__ ("movsd %0, %%xmm7 \n\t" \
                      "unpcklpd %%xmm6, %%xmm6 \n\t" \
                      "unpcklpd %%xmm7, %%xmm7 \n\t" \
                      "mulpd %%xmm2, %%xmm6 \n\t" \
                      "mulpd %%xmm0, %%xmm7 \n\t" \
                      "addpd %%xmm6, %%xmm5 \n\t" \
                      "addpd %%xmm7, %%xmm4 \n\t" \
                      "movsd %1, %%xmm6 \n\t" \
                      "movsd %2, %%xmm7 \n\t" \
                      "unpcklpd %%xmm6, %%xmm6 \n\t" \
                      "unpcklpd %%xmm7, %%xmm7 \n\t" \
                      "mulpd %%xmm1, %%xmm6 \n\t" \
                      "mulpd %%xmm0, %%xmm7 \n\t" \
                      "addpd %%xmm6, %%xmm3 \n\t" \
                      "addpd %%xmm7, %%xmm5 \n\t" \
                      "movsd %3, %%xmm0 \n\t" \
                      "movsd %4, %%xmm6 \n\t" \
                      "movsd %5, %%xmm7" \
                      : \
                      : \
                      "m" ((aa)->e[0][1].imag), \
                      "m" ((aa)->e[1][0].imag), \
                      "m" ((aa)->e[0][2].imag), \
                      "m" ((aa)->e[2][0].imag), \
                      "m" ((aa)->e[1][2].imag), \
                      "m" ((aa)->e[2][1].imag)\
                      : \
                      "%xmm0", \
                      "%xmm1", \
                      "%xmm2", \
                      "%xmm3", \
                      "%xmm4", \
                      "%xmm5", \
                      "%xmm6", \
                      "%xmm7", \
                      "memory"); \
__asm__ __volatile__ ("unpcklpd %%xmm0, %%xmm0 \n\t" \
                      "unpcklpd %%xmm6, %%xmm6 \n\t" \
                      "unpcklpd %%xmm7, %%xmm7 \n\t" \
                      "mulpd %%xmm2, %%xmm0 \n\t" \
                      "mulpd %%xmm1, %%xmm6 \n\t" \
                      "mulpd %%xmm2, %%xmm7 \n\t" \
                      "addpd %%xmm0, %%xmm3 \n\t" \
                      "addpd %%xmm7, %%xmm4 \n\t" \
                      "addpd %%xmm6, %%xmm5 \n\t" \
                      "movupd %%xmm3, %0 \n\t" \
                      "movupd %%xmm4, %1 \n\t" \
                      "movupd %%xmm5, %2 \n\t" \
                      : \
                      "=m" ((cc)->h[0].c[0]), \
                      "=m" ((cc)->h[0].c[1]), \
                      "=m" ((cc)->h[0].c[2])\
                      : \
                      : \
                      "%xmm0", \
                      "%xmm1", \
                      "%xmm2", \
                      "%xmm3", \
                      "%xmm4", \
                      "%xmm5", \
                      "%xmm6", \
                      "%xmm7", \
                      "memory"); \
__asm__ __volatile__ ("movupd %0, %%xmm0 \n\t" \
                      "movupd %1, %%xmm1 \n\t" \
                      "movupd %2, %%xmm2 \n\t" \
                      "movsd %3, %%xmm3 \n\t" \
                      "movsd %4, %%xmm6 \n\t" \
                      "movsd %5, %%xmm4" \
                      : \
                      : \
                      "m" ((bb)->h[1].c[0]), \
                      "m" ((bb)->h[1].c[1]), \
                      "m" ((bb)->h[1].c[2]), \
                      "m" ((aa)->e[0][0].real), \
                      "m" ((aa)->e[1][0].real), \
                      "m" ((aa)->e[0][1].real)\
                      : \
                      "%xmm0", \
                      "%xmm1", \
                      "%xmm2", \
                      "%xmm3", \
                      "%xmm4", \
                      "%xmm6", \
                      "memory"); \
__asm__ __volatile__ ("movsd %0, %%xmm7 \n\t" \
                      "movsd %1, %%xmm5 \n\t" \
                      "unpcklpd %%xmm3, %%xmm3 \n\t" \
                      "unpcklpd %%xmm6, %%xmm6 \n\t" \
                      "unpcklpd %%xmm4, %%xmm4 \n\t" \
                      "mulpd %%xmm0, %%xmm3 \n\t" \
                      "unpcklpd %%xmm7, %%xmm7 \n\t" \
                      "mulpd %%xmm1, %%xmm6 \n\t" \
                      "unpcklpd %%xmm5, %%xmm5 \n\t" \
                      "mulpd %%xmm0, %%xmm4 \n\t" \
                      "addpd %%xmm6, %%xmm3 \n\t" \
                      "mulpd %%xmm2, %%xmm7 \n\t" \
                      "mulpd %%xmm0, %%xmm5 \n\t" \
                      "addpd %%xmm7, %%xmm4 \n\t" \
                      "movsd %2, %%xmm6 \n\t" \
                      "movsd %3, %%xmm7 \n\t" \
                      "unpcklpd %%xmm6, %%xmm6 \n\t" \
                      "unpcklpd %%xmm7, %%xmm7 \n\t" \
                      "mulpd %%xmm1, %%xmm6 \n\t" \
                      "mulpd %%xmm2, %%xmm7 \n\t" \
                      "addpd %%xmm6, %%xmm5 \n\t" \
                      "addpd %%xmm7, %%xmm3 \n\t" \
                      "movsd %4, %%xmm6 \n\t" \
                      "movsd %5, %%xmm7" \
                      : \
                      : \
                      "m" ((aa)->e[2][1].real), \
                      "m" ((aa)->e[0][2].real), \
                      "m" ((aa)->e[1][2].real), \
                      "m" ((aa)->e[2][0].real), \
                      "m" ((aa)->e[1][1].real), \
                      "m" ((aa)->e[2][2].real)\
                      : \
                      "%xmm0", \
                      "%xmm1", \
                      "%xmm2", \
                      "%xmm3", \
                      "%xmm4", \
                      "%xmm5", \
                      "%xmm6", \
                      "%xmm7", \
                      "memory"); \
__asm__ __volatile__ ("unpcklpd %%xmm6, %%xmm6 \n\t" \
                      "unpcklpd %%xmm7, %%xmm7 \n\t" \
                      "mulpd %%xmm1, %%xmm6 \n\t" \
                      "mulpd %%xmm2, %%xmm7 \n\t" \
                      "addpd %%xmm6, %%xmm4 \n\t" \
                      "addpd %%xmm7, %%xmm5 \n\t" \
                      "movsd %0, %%xmm6 \n\t" \
                      "movsd %1, %%xmm7 \n\t" \
                      "shufpd $0x1, %%xmm0, %%xmm0 \n\t" \
                      "shufpd $0x1, %%xmm1, %%xmm1 \n\t" \
                      "shufpd $0x1, %%xmm2, %%xmm2 \n\t" \
                      "unpcklpd %%xmm6, %%xmm6 \n\t" \
                      "unpcklpd %%xmm7, %%xmm7 \n\t" \
                      "xorpd %2, %%xmm0 \n\t" \
                      "xorpd %3, %%xmm1 \n\t" \
                      "xorpd %4, %%xmm2 \n\t" \
                      "mulpd %%xmm0, %%xmm6 \n\t" \
                      "mulpd %%xmm1, %%xmm7 \n\t" \
                      "addpd %%xmm6, %%xmm3 \n\t" \
                      "addpd %%xmm7, %%xmm4 \n\t" \
                      "movsd %5, %%xmm6" \
                      : \
                      : \
                      "m" ((aa)->e[0][0].imag), \
                      "m" ((aa)->e[1][1].imag), \
                      "m" (_sse_sgn4), \
                      "m" (_sse_sgn4), \
                      "m" (_sse_sgn4), \
                      "m" ((aa)->e[2][2].imag)\
                      : \
                      "%xmm0", \
                      "%xmm1", \
                      "%xmm2", \
                      "%xmm3", \
                      "%xmm4", \
                      "%xmm5", \
                      "%xmm6", \
                      "%xmm7", \
                      "memory"); \
__asm__ __volatile__ ("movsd %0, %%xmm7 \n\t" \
                      "unpcklpd %%xmm6, %%xmm6 \n\t" \
                      "unpcklpd %%xmm7, %%xmm7 \n\t" \
                      "mulpd %%xmm2, %%xmm6 \n\t" \
                      "mulpd %%xmm0, %%xmm7 \n\t" \
                      "addpd %%xmm6, %%xmm5 \n\t" \
                      "addpd %%xmm7, %%xmm4 \n\t" \
                      "movsd %1, %%xmm6 \n\t" \
                      "movsd %2, %%xmm7 \n\t" \
                      "unpcklpd %%xmm6, %%xmm6 \n\t" \
                      "unpcklpd %%xmm7, %%xmm7 \n\t" \
                      "mulpd %%xmm1, %%xmm6 \n\t" \
                      "mulpd %%xmm0, %%xmm7 \n\t" \
                      "addpd %%xmm6, %%xmm3 \n\t" \
                      "addpd %%xmm7, %%xmm5 \n\t" \
                      "movsd %3, %%xmm0 \n\t" \
                      "movsd %4, %%xmm6 \n\t" \
                      "movsd %5, %%xmm7" \
                      : \
                      : \
                      "m" ((aa)->e[0][1].imag), \
                      "m" ((aa)->e[1][0].imag), \
                      "m" ((aa)->e[0][2].imag), \
                      "m" ((aa)->e[2][0].imag), \
                      "m" ((aa)->e[1][2].imag), \
                      "m" ((aa)->e[2][1].imag)\
                      : \
                      "%xmm0", \
                      "%xmm1", \
                      "%xmm2", \
                      "%xmm3", \
                      "%xmm4", \
                      "%xmm5", \
                      "%xmm6", \
                      "%xmm7", \
                      "memory"); \
__asm__ __volatile__ ("unpcklpd %%xmm0, %%xmm0 \n\t" \
                      "unpcklpd %%xmm6, %%xmm6 \n\t" \
                      "unpcklpd %%xmm7, %%xmm7 \n\t" \
                      "mulpd %%xmm2, %%xmm0 \n\t" \
                      "mulpd %%xmm1, %%xmm6 \n\t" \
                      "mulpd %%xmm2, %%xmm7 \n\t" \
                      "addpd %%xmm0, %%xmm3 \n\t" \
                      "addpd %%xmm7, %%xmm4 \n\t" \
                      "addpd %%xmm6, %%xmm5 \n\t" \
                      "movupd %%xmm3, %0 \n\t" \
                      "movupd %%xmm4, %1 \n\t" \
                      "movupd %%xmm5, %2 \n\t" \
                      : \
                      "=m" ((cc)->h[1].c[0]), \
                      "=m" ((cc)->h[1].c[1]), \
                      "=m" ((cc)->h[1].c[2])\
                      : \
                      : \
                      "%xmm0", \
                      "%xmm1", \
                      "%xmm2", \
                      "%xmm3", \
                      "%xmm4", \
                      "%xmm5", \
                      "%xmm6", \
                      "%xmm7", \
                      "memory"); \
}
