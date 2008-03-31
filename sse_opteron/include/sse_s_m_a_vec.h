#define _inline_sse_scalar_mult_add_su3_vector(aa,bb,cc,dd) \
{ \
__asm__ __volatile__ ("movss %0, %%xmm0 \n\t" \
                      "cvtss2sd %%xmm0, %%xmm6 \n\t" \
                      "unpcklpd %%xmm6, %%xmm6 \n\t" \
                      "cvtps2pd %1, %%xmm3 \n\t" \
                      "mulpd %%xmm6, %%xmm3 \n\t" \
                      "cvtps2pd %2, %%xmm0 \n\t" \
                      "addpd %%xmm3, %%xmm0 \n\t" \
                      "cvtps2pd %3, %%xmm4 \n\t" \
                      "mulpd %%xmm6, %%xmm4 \n\t" \
                      "cvtps2pd %4, %%xmm1 \n\t" \
                      "addpd %%xmm4, %%xmm1 \n\t" \
                      "cvtps2pd %5, %%xmm5" \
                      : \
                      : \
                      "m" (cc), \
                      "m" ((bb)->c[0]), \
                      "m" ((aa)->c[0]), \
                      "m" ((bb)->c[1]), \
                      "m" ((aa)->c[1]), \
                      "m" ((bb)->c[2])\
                      : \
                      "%xmm0", \
                      "%xmm1", \
                      "%xmm3", \
                      "%xmm4", \
                      "%xmm5", \
                      "%xmm6", \
                      "memory"); \
__asm__ __volatile__ ("mulpd %%xmm6, %%xmm5 \n\t" \
                      "cvtps2pd %0, %%xmm2 \n\t" \
                      "addpd %%xmm5, %%xmm2 \n\t" \
                      "cvtpd2ps %%xmm0, %%xmm3 \n\t" \
                      "cvtpd2ps %%xmm1, %%xmm4 \n\t" \
                      "cvtpd2ps %%xmm2, %%xmm5 \n\t" \
                      : \
                      : \
                      "m" ((aa)->c[2])\
                      : \
                      "%xmm0", \
                      "%xmm1", \
                      "%xmm2", \
                      "%xmm3", \
                      "%xmm4", \
                      "%xmm5", \
                      "%xmm6", \
                      "memory"); \
__asm__ __volatile__ ("movlps %%xmm3, %0 \n\t" \
                      "movlps %%xmm4, %1 \n\t" \
                      "movlps %%xmm5, %2 \n\t" \
                      : \
                      "=m" ((dd)->c[0]), \
                      "=m" ((dd)->c[1]), \
                      "=m" ((dd)->c[2])\
                      : \
                      : \
                      "%xmm3", \
                      "%xmm4", \
                      "%xmm5", \
                      "memory"); \
}
