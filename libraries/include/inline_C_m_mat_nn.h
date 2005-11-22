/*******************  m_mat_nn.c  (in su3.a) ****************************
*									*
* void mult_su3_nn( su3_matrix *a,*b,*c )				*
* matrix multiply, no adjoints 						*
* C  <-  A*B								*
*/

/*
for(i=0; i<3; i++) {
  for(j=0; j<3; j++) {
    c[i][j] = 0;
  }
  for(k=0; k<3; k++) {
    for(j=0; j<3; j++) {
      c[i][j] += a[i][k]*b[k][j];
    }
  }
}
*/

#if 1

#ifndef __inline_defs
#define __inline_defs

#define er(x,i,j) (((su3_matrix *)(x))->e[i][j].real)
#define ei(x,i,j) (((su3_matrix *)(x))->e[i][j].imag)

#define cop(rr,ri,op,ar,ai,br,bi) \
  rr op ar*br - ai*bi; \
  ri op ar*bi + ai*br

#define copr(rr,ri,op,ar,ai,br,bi) \
  rr op ar*br; \
  ri op ar*bi

#define copi(rr,ri,op,ar,ai,br,bi) \
  rr -= ai*bi; \
  ri += ai*br

#define copna(rr,ri,op,ar,ai,br,bi) \
  rr op ar*br + ai*bi; \
  ri op ai*br - ar*bi

#define copnar(rr,ri,op,ar,ai,br,bi) \
  rr op ar*br; \
  ri op ai*br

#define copnai(rr,ri,op,ar,ai,br,bi) \
  rr += ai*bi; \
  ri -= ar*bi

#endif

#define colopnn2(op,aa,bb,i,k) \
  aikr = er(aa,i,k); \
  aiki = ei(aa,i,k); \
  copr(c0r,c0i,op,aikr,aiki,er(bb,k,0),ei(bb,k,0)); \
  copr(c1r,c1i,op,aikr,aiki,er(bb,k,1),ei(bb,k,1)); \
  copr(c2r,c2i,op,aikr,aiki,er(bb,k,2),ei(bb,k,2)); \
  copi(c0r,c0i,op,aikr,aiki,er(bb,k,0),ei(bb,k,0)); \
  copi(c1r,c1i,op,aikr,aiki,er(bb,k,1),ei(bb,k,1)); \
  copi(c2r,c2i,op,aikr,aiki,er(bb,k,2),ei(bb,k,2))

#define colopnn(op,aa,bb,i,k) \
  aikr = er(aa,i,k); \
  aiki = ei(aa,i,k); \
  cop(c0r,c0i,op,aikr,aiki,er(bb,k,0),ei(bb,k,0)); \
  cop(c1r,c1i,op,aikr,aiki,er(bb,k,1),ei(bb,k,1)); \
  cop(c2r,c2i,op,aikr,aiki,er(bb,k,2),ei(bb,k,2))

#define _inline_C_mult_su3_nn( aa, bb, cc ) { \
  int iii; \
  for(iii=0; iii<3; iii++) { \
    Real c0r, c0i, c1r, c1i, c2r, c2i; \
    Real aikr, aiki; \
    colopnn2(=,aa,bb,iii,0); \
    colopnn2(+=,aa,bb,iii,1); \
    colopnn2(+=,aa,bb,iii,2); \
    er(cc,iii,0) = c0r; \
    ei(cc,iii,0) = c0i; \
    er(cc,iii,1) = c1r; \
    ei(cc,iii,1) = c1i; \
    er(cc,iii,2) = c2r; \
    ei(cc,iii,2) = c2i; \
  } \
}

#else
#define  _inline_C_mult_su3_nn( aa, bb, cc ) { \
  register su3_matrix *aaa, *bbb, *ccc; \
  register int iii,jjj; \
  register Real t,ar,ai,br,bi,cr,ci; \
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
#endif
