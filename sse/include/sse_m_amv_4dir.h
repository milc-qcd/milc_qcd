#define _inline_sse_mult_adj_su3_mat_vec_4dir(aa,bb,cc) \
{ \
__asm__ __volatile__ ("movups %0, %%xmm0 \n\t" \
                      "movaps %%xmm0, %%xmm1 \n\t" \
                      "shufps $0xB1, %%xmm1, %%xmm1 \n\t" \
                      "movups %1, %%xmm2 \n\t" \
                      "shufps $0xEB, %%xmm2, %%xmm2 \n\t" \
                      "movlps %2, %%xmm3 \n\t" \
                      "movhps %3, %%xmm3 \n\t" \
                      "movhps %4, %%xmm4 \n\t" \
                      "shufps $0xEE, %%xmm4, %%xmm4 \n\t" \
                      "movaps %%xmm3, %%xmm5 \n\t" \
                      "mulps %%xmm0, %%xmm3 \n\t" \
                      "mulps %%xmm1, %%xmm5 \n\t" \
                      "mulps %%xmm2, %%xmm4 \n\t" \
                      "movaps %%xmm5, %%xmm6 \n\t" \
                      "shufps $0x4E, %%xmm3, %%xmm6 \n\t" \
                      "shufps $0xE4, %%xmm3, %%xmm5 \n\t" \
                      "addps %%xmm5, %%xmm4 \n\t" \
                      "addps %%xmm6, %%xmm4 \n\t" \
                      "movlps %5, %%xmm3" \
                      : \
                      : \
                      "m" ((bb)->c[0]), \
                      "m" ((bb)->c[1]), \
                      "m" ((aa)[0].e[0][0]), \
                      "m" ((aa)[0].e[1][0]), \
                      "m" ((aa)[0].e[2][0]), \
                      "m" ((aa)[0].e[0][1])\
                      : \
                      "%xmm0", \
                      "%xmm1", \
                      "%xmm2", \
                      "%xmm3", \
                      "%xmm4", \
                      "%xmm5", \
                      "%xmm6", \
                      "memory"); \
__asm__ __volatile__ ("movhps %0, %%xmm3 \n\t" \
                      "movhps %1, %%xmm7 \n\t" \
                      "shufps $0xEE, %%xmm7, %%xmm7 \n\t" \
                      "movaps %%xmm3, %%xmm5 \n\t" \
                      "mulps %%xmm0, %%xmm3 \n\t" \
                      "mulps %%xmm1, %%xmm5 \n\t" \
                      "mulps %%xmm2, %%xmm7 \n\t" \
                      "movaps %%xmm5, %%xmm6 \n\t" \
                      "shufps $0x4E, %%xmm3, %%xmm6 \n\t" \
                      "shufps $0xE4, %%xmm3, %%xmm5 \n\t" \
                      "addps %%xmm5, %%xmm7 \n\t" \
                      "addps %%xmm6, %%xmm7 \n\t" \
                      "movaps %%xmm4, %%xmm5 \n\t" \
                      "shufps $0x22, %%xmm7, %%xmm5 \n\t" \
                      "shufps $0x77, %%xmm7, %%xmm4 \n\t" \
                      "xorps %2, %%xmm4 \n\t" \
                      "addps %%xmm4, %%xmm5 \n\t" \
                      : \
                      : \
                      "m" ((aa)[0].e[1][1]), \
                      "m" ((aa)[0].e[2][1]), \
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
                      "memory"); \
__asm__ __volatile__ ("movups %%xmm5, %0 \n\t" \
                      : \
                      "=m" ((cc)[0].c[0])\
                      : \
                      : \
                      "%xmm5", \
                      "memory"); \
__asm__ __volatile__ ("movlps %0, %%xmm3 \n\t" \
                      "movhps %1, %%xmm3 \n\t" \
                      "movhps %2, %%xmm4 \n\t" \
                      "shufps $0xEE, %%xmm4, %%xmm4 \n\t" \
                      "movaps %%xmm3, %%xmm5 \n\t" \
                      "mulps %%xmm0, %%xmm3 \n\t" \
                      "mulps %%xmm1, %%xmm5 \n\t" \
                      "mulps %%xmm2, %%xmm4 \n\t" \
                      "movaps %%xmm5, %%xmm6 \n\t" \
                      "shufps $0x4E, %%xmm3, %%xmm6 \n\t" \
                      "shufps $0xE4, %%xmm3, %%xmm5 \n\t" \
                      "addps %%xmm5, %%xmm4 \n\t" \
                      "addps %%xmm6, %%xmm4 \n\t" \
                      "movaps %%xmm4, %%xmm5 \n\t" \
                      "shufps $0x77, %%xmm4, %%xmm4 \n\t" \
                      "shufps $0x22, %%xmm5, %%xmm5 \n\t" \
                      "xorps %3, %%xmm4 \n\t" \
                      "addps %%xmm4, %%xmm5 \n\t" \
                      : \
                      : \
                      "m" ((aa)[0].e[0][2]), \
                      "m" ((aa)[0].e[1][2]), \
                      "m" ((aa)[0].e[2][2]), \
                      "m" (_sse_sgn24)\
                      : \
                      "%xmm0", \
                      "%xmm1", \
                      "%xmm2", \
                      "%xmm3", \
                      "%xmm4", \
                      "%xmm5", \
                      "%xmm6", \
                      "memory"); \
__asm__ __volatile__ ("movhps %%xmm5, %0 \n\t" \
                      : \
                      "=m" ((cc)[0].c[2])\
                      : \
                      : \
                      "%xmm5", \
                      "memory"); \
__asm__ __volatile__ ("movlps %0, %%xmm3 \n\t" \
                      "movhps %1, %%xmm3 \n\t" \
                      "movhps %2, %%xmm4 \n\t" \
                      "shufps $0xEE, %%xmm4, %%xmm4 \n\t" \
                      "movaps %%xmm3, %%xmm5 \n\t" \
                      "mulps %%xmm0, %%xmm3 \n\t" \
                      "mulps %%xmm1, %%xmm5 \n\t" \
                      "mulps %%xmm2, %%xmm4 \n\t" \
                      "movaps %%xmm5, %%xmm6 \n\t" \
                      "shufps $0x4E, %%xmm3, %%xmm6 \n\t" \
                      "shufps $0xE4, %%xmm3, %%xmm5 \n\t" \
                      "addps %%xmm5, %%xmm4 \n\t" \
                      "addps %%xmm6, %%xmm4 \n\t" \
                      "movlps %3, %%xmm3 \n\t" \
                      "movhps %4, %%xmm3 \n\t" \
                      "movhps %5, %%xmm7" \
                      : \
                      : \
                      "m" ((aa)[1].e[0][0]), \
                      "m" ((aa)[1].e[1][0]), \
                      "m" ((aa)[1].e[2][0]), \
                      "m" ((aa)[1].e[0][1]), \
                      "m" ((aa)[1].e[1][1]), \
                      "m" ((aa)[1].e[2][1])\
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
__asm__ __volatile__ ("shufps $0xEE, %%xmm7, %%xmm7 \n\t" \
                      "movaps %%xmm3, %%xmm5 \n\t" \
                      "mulps %%xmm0, %%xmm3 \n\t" \
                      "mulps %%xmm1, %%xmm5 \n\t" \
                      "mulps %%xmm2, %%xmm7 \n\t" \
                      "movaps %%xmm5, %%xmm6 \n\t" \
                      "shufps $0x4E, %%xmm3, %%xmm6 \n\t" \
                      "shufps $0xE4, %%xmm3, %%xmm5 \n\t" \
                      "addps %%xmm5, %%xmm7 \n\t" \
                      "addps %%xmm6, %%xmm7 \n\t" \
                      "movaps %%xmm4, %%xmm5 \n\t" \
                      "shufps $0x22, %%xmm7, %%xmm5 \n\t" \
                      "shufps $0x77, %%xmm7, %%xmm4 \n\t" \
                      "xorps %0, %%xmm4 \n\t" \
                      "addps %%xmm4, %%xmm5 \n\t" \
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
                      "memory"); \
__asm__ __volatile__ ("movups %%xmm5, %0 \n\t" \
                      : \
                      "=m" ((cc)[1].c[0])\
                      : \
                      : \
                      "%xmm5", \
                      "memory"); \
__asm__ __volatile__ ("movlps %0, %%xmm3 \n\t" \
                      "movhps %1, %%xmm3 \n\t" \
                      "movhps %2, %%xmm4 \n\t" \
                      "shufps $0xEE, %%xmm4, %%xmm4 \n\t" \
                      "movaps %%xmm3, %%xmm5 \n\t" \
                      "mulps %%xmm0, %%xmm3 \n\t" \
                      "mulps %%xmm1, %%xmm5 \n\t" \
                      "mulps %%xmm2, %%xmm4 \n\t" \
                      "movaps %%xmm5, %%xmm6 \n\t" \
                      "shufps $0x4E, %%xmm3, %%xmm6 \n\t" \
                      "shufps $0xE4, %%xmm3, %%xmm5 \n\t" \
                      "addps %%xmm5, %%xmm4 \n\t" \
                      "addps %%xmm6, %%xmm4 \n\t" \
                      "movaps %%xmm4, %%xmm5 \n\t" \
                      "shufps $0x77, %%xmm4, %%xmm4 \n\t" \
                      "shufps $0x22, %%xmm5, %%xmm5 \n\t" \
                      "xorps %3, %%xmm4 \n\t" \
                      "addps %%xmm4, %%xmm5 \n\t" \
                      : \
                      : \
                      "m" ((aa)[1].e[0][2]), \
                      "m" ((aa)[1].e[1][2]), \
                      "m" ((aa)[1].e[2][2]), \
                      "m" (_sse_sgn24)\
                      : \
                      "%xmm0", \
                      "%xmm1", \
                      "%xmm2", \
                      "%xmm3", \
                      "%xmm4", \
                      "%xmm5", \
                      "%xmm6", \
                      "memory"); \
__asm__ __volatile__ ("movhps %%xmm5, %0 \n\t" \
                      : \
                      "=m" ((cc)[1].c[2])\
                      : \
                      : \
                      "%xmm5", \
                      "memory"); \
__asm__ __volatile__ ("movlps %0, %%xmm3 \n\t" \
                      "movhps %1, %%xmm3 \n\t" \
                      "movhps %2, %%xmm4 \n\t" \
                      "shufps $0xEE, %%xmm4, %%xmm4 \n\t" \
                      "movaps %%xmm3, %%xmm5 \n\t" \
                      "mulps %%xmm0, %%xmm3 \n\t" \
                      "mulps %%xmm1, %%xmm5 \n\t" \
                      "mulps %%xmm2, %%xmm4 \n\t" \
                      "movaps %%xmm5, %%xmm6 \n\t" \
                      "shufps $0x4E, %%xmm3, %%xmm6 \n\t" \
                      "shufps $0xE4, %%xmm3, %%xmm5 \n\t" \
                      "addps %%xmm5, %%xmm4 \n\t" \
                      "addps %%xmm6, %%xmm4 \n\t" \
                      "movlps %3, %%xmm3 \n\t" \
                      "movhps %4, %%xmm3 \n\t" \
                      "movhps %5, %%xmm7" \
                      : \
                      : \
                      "m" ((aa)[2].e[0][0]), \
                      "m" ((aa)[2].e[1][0]), \
                      "m" ((aa)[2].e[2][0]), \
                      "m" ((aa)[2].e[0][1]), \
                      "m" ((aa)[2].e[1][1]), \
                      "m" ((aa)[2].e[2][1])\
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
__asm__ __volatile__ ("shufps $0xEE, %%xmm7, %%xmm7 \n\t" \
                      "movaps %%xmm3, %%xmm5 \n\t" \
                      "mulps %%xmm0, %%xmm3 \n\t" \
                      "mulps %%xmm1, %%xmm5 \n\t" \
                      "mulps %%xmm2, %%xmm7 \n\t" \
                      "movaps %%xmm5, %%xmm6 \n\t" \
                      "shufps $0x4E, %%xmm3, %%xmm6 \n\t" \
                      "shufps $0xE4, %%xmm3, %%xmm5 \n\t" \
                      "addps %%xmm5, %%xmm7 \n\t" \
                      "addps %%xmm6, %%xmm7 \n\t" \
                      "movaps %%xmm4, %%xmm5 \n\t" \
                      "shufps $0x22, %%xmm7, %%xmm5 \n\t" \
                      "shufps $0x77, %%xmm7, %%xmm4 \n\t" \
                      "xorps %0, %%xmm4 \n\t" \
                      "addps %%xmm4, %%xmm5 \n\t" \
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
                      "memory"); \
__asm__ __volatile__ ("movups %%xmm5, %0 \n\t" \
                      : \
                      "=m" ((cc)[2].c[0])\
                      : \
                      : \
                      "%xmm5", \
                      "memory"); \
__asm__ __volatile__ ("movlps %0, %%xmm3 \n\t" \
                      "movhps %1, %%xmm3 \n\t" \
                      "movhps %2, %%xmm4 \n\t" \
                      "shufps $0xEE, %%xmm4, %%xmm4 \n\t" \
                      "movaps %%xmm3, %%xmm5 \n\t" \
                      "mulps %%xmm0, %%xmm3 \n\t" \
                      "mulps %%xmm1, %%xmm5 \n\t" \
                      "mulps %%xmm2, %%xmm4 \n\t" \
                      "movaps %%xmm5, %%xmm6 \n\t" \
                      "shufps $0x4E, %%xmm3, %%xmm6 \n\t" \
                      "shufps $0xE4, %%xmm3, %%xmm5 \n\t" \
                      "addps %%xmm5, %%xmm4 \n\t" \
                      "addps %%xmm6, %%xmm4 \n\t" \
                      "movaps %%xmm4, %%xmm5 \n\t" \
                      "shufps $0x77, %%xmm4, %%xmm4 \n\t" \
                      "shufps $0x22, %%xmm5, %%xmm5 \n\t" \
                      "xorps %3, %%xmm4 \n\t" \
                      "addps %%xmm4, %%xmm5 \n\t" \
                      : \
                      : \
                      "m" ((aa)[2].e[0][2]), \
                      "m" ((aa)[2].e[1][2]), \
                      "m" ((aa)[2].e[2][2]), \
                      "m" (_sse_sgn24)\
                      : \
                      "%xmm0", \
                      "%xmm1", \
                      "%xmm2", \
                      "%xmm3", \
                      "%xmm4", \
                      "%xmm5", \
                      "%xmm6", \
                      "memory"); \
__asm__ __volatile__ ("movhps %%xmm5, %0 \n\t" \
                      : \
                      "=m" ((cc)[2].c[2])\
                      : \
                      : \
                      "%xmm5", \
                      "memory"); \
__asm__ __volatile__ ("movlps %0, %%xmm3 \n\t" \
                      "movhps %1, %%xmm3 \n\t" \
                      "movhps %2, %%xmm4 \n\t" \
                      "shufps $0xEE, %%xmm4, %%xmm4 \n\t" \
                      "movaps %%xmm3, %%xmm5 \n\t" \
                      "mulps %%xmm0, %%xmm3 \n\t" \
                      "mulps %%xmm1, %%xmm5 \n\t" \
                      "mulps %%xmm2, %%xmm4 \n\t" \
                      "movaps %%xmm5, %%xmm6 \n\t" \
                      "shufps $0x4E, %%xmm3, %%xmm6 \n\t" \
                      "shufps $0xE4, %%xmm3, %%xmm5 \n\t" \
                      "addps %%xmm5, %%xmm4 \n\t" \
                      "addps %%xmm6, %%xmm4 \n\t" \
                      "movlps %3, %%xmm3 \n\t" \
                      "movhps %4, %%xmm3 \n\t" \
                      "movhps %5, %%xmm7" \
                      : \
                      : \
                      "m" ((aa)[3].e[0][0]), \
                      "m" ((aa)[3].e[1][0]), \
                      "m" ((aa)[3].e[2][0]), \
                      "m" ((aa)[3].e[0][1]), \
                      "m" ((aa)[3].e[1][1]), \
                      "m" ((aa)[3].e[2][1])\
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
__asm__ __volatile__ ("shufps $0xEE, %%xmm7, %%xmm7 \n\t" \
                      "movaps %%xmm3, %%xmm5 \n\t" \
                      "mulps %%xmm0, %%xmm3 \n\t" \
                      "mulps %%xmm1, %%xmm5 \n\t" \
                      "mulps %%xmm2, %%xmm7 \n\t" \
                      "movaps %%xmm5, %%xmm6 \n\t" \
                      "shufps $0x4E, %%xmm3, %%xmm6 \n\t" \
                      "shufps $0xE4, %%xmm3, %%xmm5 \n\t" \
                      "addps %%xmm5, %%xmm7 \n\t" \
                      "addps %%xmm6, %%xmm7 \n\t" \
                      "movaps %%xmm4, %%xmm5 \n\t" \
                      "shufps $0x22, %%xmm7, %%xmm5 \n\t" \
                      "shufps $0x77, %%xmm7, %%xmm4 \n\t" \
                      "xorps %0, %%xmm4 \n\t" \
                      "addps %%xmm4, %%xmm5 \n\t" \
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
                      "memory"); \
__asm__ __volatile__ ("movups %%xmm5, %0 \n\t" \
                      : \
                      "=m" ((cc)[3].c[0])\
                      : \
                      : \
                      "%xmm5", \
                      "memory"); \
__asm__ __volatile__ ("movlps %0, %%xmm3 \n\t" \
                      "movhps %1, %%xmm3 \n\t" \
                      "movhps %2, %%xmm4 \n\t" \
                      "shufps $0xEE, %%xmm4, %%xmm4 \n\t" \
                      "movaps %%xmm3, %%xmm5 \n\t" \
                      "mulps %%xmm0, %%xmm3 \n\t" \
                      "mulps %%xmm1, %%xmm5 \n\t" \
                      "mulps %%xmm2, %%xmm4 \n\t" \
                      "movaps %%xmm5, %%xmm6 \n\t" \
                      "shufps $0x4E, %%xmm3, %%xmm6 \n\t" \
                      "shufps $0xE4, %%xmm3, %%xmm5 \n\t" \
                      "addps %%xmm5, %%xmm4 \n\t" \
                      "addps %%xmm6, %%xmm4 \n\t" \
                      "movaps %%xmm4, %%xmm5 \n\t" \
                      "shufps $0x77, %%xmm4, %%xmm4 \n\t" \
                      "shufps $0x22, %%xmm5, %%xmm5 \n\t" \
                      "xorps %3, %%xmm4 \n\t" \
                      "addps %%xmm4, %%xmm5 \n\t" \
                      : \
                      : \
                      "m" ((aa)[3].e[0][2]), \
                      "m" ((aa)[3].e[1][2]), \
                      "m" ((aa)[3].e[2][2]), \
                      "m" (_sse_sgn24)\
                      : \
                      "%xmm0", \
                      "%xmm1", \
                      "%xmm2", \
                      "%xmm3", \
                      "%xmm4", \
                      "%xmm5", \
                      "%xmm6", \
                      "memory"); \
__asm__ __volatile__ ("movhps %%xmm5, %0 \n\t" \
                      : \
                      "=m" ((cc)[3].c[2])\
                      : \
                      : \
                      "%xmm5", \
                      "memory"); \
}
