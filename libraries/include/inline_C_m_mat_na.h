/****************  m_mat_na.c  (in su3.a) *******************************
*									*
* void mult_su3_na( su3_matrix *a,*b,*c )				*
* matrix multiply, second matrix is adjoint 				*
* C  <-  A*B_adjoint							*
*/
#define  _inline_C_mult_su3_na( aa, bb, cc ) { \
  register su3_matrix *aaa, *bbb, *ccc; \
  register int iii,jjj; \
  register Real t,ar,ai,br,bi,cr,ci; \
    aaa = (aa) ; bbb = (bb) ; ccc = (cc) ; \
    for(iii=0;iii<3;iii++)for(jjj=0;jjj<3;jjj++){ \
	ar=aaa->e[iii][0].real; ai=aaa->e[iii][0].imag; \
	br=bbb->e[jjj][0].real; bi=bbb->e[jjj][0].imag; \
	cr=ar*br; t=ai*bi; cr += t; \
	ci=ai*br; t=ar*bi; ci -= t; \
	ar=aaa->e[iii][1].real; ai=aaa->e[iii][1].imag; \
	br=bbb->e[jjj][1].real; bi=bbb->e[jjj][1].imag; \
	t=ar*br; cr += t; t=ai*bi; cr += t; \
	t=ar*bi; ci -= t; t=ai*br; ci += t; \
	ar=aaa->e[iii][2].real; ai=aaa->e[iii][2].imag; \
	br=bbb->e[jjj][2].real; bi=bbb->e[jjj][2].imag; \
	t=ar*br; cr += t; t=ai*bi; cr += t; \
	t=ar*bi; ci -= t; t=ai*br; ci += t; \
	ccc->e[iii][jjj].real=cr; \
	ccc->e[iii][jjj].imag=ci; \
    } \
}
