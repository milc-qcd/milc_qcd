/*******************  m_mat_nn.c  (in su3.a) ****************************
*									*
* void mult_su3_nn( su3_matrix *a,*b,*c )				*
* matrix multiply, no adjoints 						*
* C  <-  A*B								*
*/
#define  _inline_C_mult_su3_nn( aa, bb, cc ) { \
  register su3_matrix *aaa, *bbb, *ccc; \
  register int iii,jjj; \
  register float t,ar,ai,br,bi,cr,ci; \
    aaa = (aa) ; bbb = (bb) ; ccc = (cc) ; \
    for(iii=0;iii<3;iii++)for(jjj=0;jjj<3;jjj++){ \
	ar=aaa->e[iii][0].real; ai=aaa->e[iii][0].imag; \
	br=bbb->e[0][jjj].real; bi=bbb->e[0][jjj].imag; \
	cr=ar*br; t=ai*bi; cr -= t; \
	ci=ar*bi; t=ai*br; ci += t; \
	ar=aaa->e[iii][1].real; ai=aaa->e[iii][1].imag; \
	br=bbb->e[1][jjj].real; bi=bbb->e[1][jjj].imag; \
	t=ar*br; cr += t; t=ai*bi; cr -= t; \
	t=ar*bi; ci += t; t=ai*br; ci += t; \
	ar=aaa->e[iii][2].real; ai=aaa->e[iii][2].imag; \
	br=bbb->e[2][jjj].real; bi=bbb->e[2][jjj].imag; \
	t=ar*br; cr += t; t=ai*bi; cr -= t; \
	t=ar*bi; ci += t; t=ai*br; ci += t; \
	ccc->e[iii][jjj].real=cr; \
	ccc->e[iii][jjj].imag=ci; \
    } \
}
