// scalar_mult_add_su3_matrix( su3_matrix *a, su3_matrix *b, float s,
//	 su3_matrix *c)
// c <- a + s*b
// file s_m_a_mat.m4, i860 version of s_m_a_mat.c
//
    define(A,r16)	// address of source 1
    define(B,r17)	// address of source 2
    define(C,r18)	// address of result
    define(S,f8)	// scalar

    define(a0,f10)	// complex number = register pair
    define(a0r,f10)	// real part
    define(a0i,f11)	// imag part

    define(b0,f16)
    define(b0r,f16)
    define(b0i,f17)

    define(c0,f22)
    define(c0r,f22)
    define(c0i,f23)

	.text
	.align	8
_scalar_mult_add_su3_matrix:

//.if FLOATOPTION=X167
      .align  8
	d.pfadd.ss f0,f0,f0;     fld.d 0(B),b0
//.else
	// convert double argument to single
//	.align	8
//	d.pfmov.ds S,f0;	fld.d 0(B),b0
//	d.pfmov.ds f0,f0;	nop
//	d.pfmov.ds f0,f0;	nop
//	d.pfmov.ds f0,S;	nop
//.endif

	// start multiplying first component, third op. is dummy
	d.pfmul.ss S,b0r,f0;	fld.d 0(A),a0
	d.pfmul.ss S,b0i,f0;	adds -8,C,C
	d.pfmul.ss f0,f0,f0;	nop
	// First component of A into adder
	d.pfadd.ss a0r,f0,f0;	fld.d 8(B)++,b0
	d.pfadd.ss a0i,f0,f0;	adds -1,r0,r20
	d.pfadd.ss f0,f0,f0;	or 6,r0,r21
	// multiply second component, add first
	d.m12apm.ss S,b0r,f0;	fld.d 8(A)++,a0
	d.m12apm.ss S,b0i,f0;	bla r20,r21,LOOP
	d.m12apm.ss f0,f0,f0;	nop
	// second component of A into adder
LOOP:	d.pfadd.ss a0r,f0,c0r;	fld.d 8(B)++,b0
	d.pfadd.ss a0i,f0,c0i;	nop
	d.pfadd.ss f0,f0,f0;	fst.d c0,8(C)++
	// multiply third component, add second
	d.m12apm.ss S,b0r,f0;	fld.d 8(A)++,a0
	d.m12apm.ss S,b0i,f0;	bla r20,r21,LOOP
	d.m12apm.ss f0,f0,f0;	nop
	// third component of A into adder
	d.pfadd.ss a0r,f0,c0r;	nop
	d.pfadd.ss a0i,f0,c0i;	nop
	d.pfadd.ss f0,f0,f0;	fst.d c0,8(C)++
	// empty pipes
	d.m12apm.ss f0,f0,f0;	nop
	d.m12apm.ss f0,f0,f0;	nop
	d.m12apm.ss f0,f0,f0;	nop
	// 
	pfadd.ss f0,f0,c0r;	nop
	pfadd.ss f0,f0,c0i;	nop
    				bri	r1
				fst.d c0,8(C)++

.globl	_scalar_mult_add_su3_matrix
