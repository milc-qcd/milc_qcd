// scalar_mult_add_su3_vector( su3_vector *a, su3_vector *b, float s,
//	 su3_vector *c)
// c <- a + s*b
// file s_m_a_vec.m4, i860 version of s_m_a_vec.c
//
    define(A,r16)	// address of source 1
    define(B,r17)	// address of source 2
    define(C,r18)	// address of result
    define(S,f8)	// scalar

    define(a0,f10)	// complex number = register pair
    define(a0r,f10)	// real part
    define(a0i,f11)	// imag part
    define(a1,f12)
    define(a1r,f12)
    define(a1i,f13)
    define(a2,f14)
    define(a2r,f14)
    define(a2i,f15)

    define(b0,f16)
    define(b0r,f16)
    define(b0i,f17)
    define(b1,f18)
    define(b1r,f18)
    define(b1i,f19)
    define(b2,f20)
    define(b2r,f20)
    define(b2i,f21)

    define(c0,f22)
    define(c0r,f22)
    define(c0i,f23)
    define(c1,f24)
    define(c1r,f24)
    define(c1i,f25)
    define(c2,f26)
    define(c2r,f26)
    define(c2i,f27)

	.text
	.align	8
_scalar_mult_add_su3_vector:

//.if FLOATOPTION=X167
      .align  8
	d.pfadd.ss f0,f0,f0;     fld.d 0(B),b0
//.else
	// convert double argument to single
//	.align	8
//	d.pfmov.ds f8,f0;	fld.d 0(B),b0
//	d.pfmov.ds f0,f0;	nop
//	d.pfmov.ds f0,f0;	nop
//	d.pfmov.ds f0,f8;	nop
//.endif

	// start multiplying first component, third op. is dummy
	d.pfmul.ss S,b0r,f0;	fld.d 0(A),a0
	d.pfmul.ss S,b0i,f0;	adds -8,C,C
	d.pfmul.ss f0,f0,f0;	nop
	// First component of A into adder
	d.pfadd.ss a0r,f0,f0;	fld.d 8(B)++,b1
	d.pfadd.ss a0i,f0,f0;	nop
	d.pfadd.ss f0,f0,f0;	nop
	// multiply second component, add first
	d.m12apm.ss S,b1r,f0;	fld.d 8(A)++,a1
	d.m12apm.ss S,b1i,f0;	nop
	d.m12apm.ss f0,f0,f0;	nop
	// second component of A into adder
	d.pfadd.ss a1r,f0,c0r;	fld.d 8(B)++,b2
	d.pfadd.ss a1i,f0,c0i;	nop
	d.pfadd.ss f0,f0,f0;	fst.d c0,8(C)++
	// multiply third component, add second
	d.m12apm.ss S,b2r,f0;	fld.d 8(A)++,a2
	d.m12apm.ss S,b2i,f0;	nop
	d.m12apm.ss f0,f0,f0;	nop
	// third component of A into adder
	d.pfadd.ss a2r,f0,c1r;	nop
	d.pfadd.ss a2i,f0,c1i;	nop
	d.pfadd.ss f0,f0,f0;	fst.d c1,8(C)++
	// empty pipes
	d.m12apm.ss f0,f0,f0;	nop
	d.m12apm.ss f0,f0,f0;	nop
	d.m12apm.ss f0,f0,f0;	nop
	// 
	pfadd.ss f0,f0,c2r;	nop
	pfadd.ss f0,f0,c2i;	nop

    				bri	r1
				fst.d c2,8(C)++

.globl	_scalar_mult_add_su3_vector
