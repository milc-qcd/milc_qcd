// float wvec_rdot( wilson_vector *a,*b );
// return real part of dot product of A and B
// file wvec_rdot.m4, i860 assembler version of wvec_rdot.c
//
    define(A,r16)	// address of source 1
    define(B,r17)	// address of source 2
    define(sum,f8)	// return value

    define(a0,f8)	// complex number = register pair
    define(a0r,f8)	// real part
    define(a0i,f9)	// imag part
    define(a1,f10)
    define(a1r,f10)
    define(a1i,f11)
    define(a2,f12)
    define(a2r,f12)
    define(a2i,f13)

    define(b0,f14)	// complex number = register pair
    define(b0r,f14)	// real part
    define(b0i,f15)	// imag part
    define(b1,f16)
    define(b1r,f16)
    define(b1i,f17)
    define(b2,f18)
    define(b2r,f18)
    define(b2i,f19)

    define(c0,f20)	// complex number = register pair
    define(c0r,f20)	// real part
    define(c0i,f21)	// imag part
    define(c1,f22)
    define(c1r,f22)
    define(c1i,f23)
    define(c2,f24)
    define(c2r,f24)
    define(c2i,f25)

    define(d0,f26)	// complex number = register pair
    define(d0r,f26)	// real part
    define(d0i,f27)	// imag part
    define(d1,f28)
    define(d1r,f28)
    define(d1i,f29)
    define(d2,f30)
    define(d2r,f30)
    define(d2i,f31)

	.text
	.align	8
_wvec_rdot:
	// enter zeroes into adder pipeline
	// (one extra, because need time for	fld's)
	d.pfadd.ss	f0,f0,f0;	fld.d	0(A),a0
	d.pfadd.ss	f0,f0,f0;	fld.d	0(B),b0
	d.pfadd.ss	f0,f0,f0;	fld.d	8(B),b1
	d.pfadd.ss	f0,f0,f0;	fld.d	8(A),a1
	d.pfmul.ss a0r,b0r,f0;	fld.d	16(A),a2
	d.pfmul.ss a0i,b0i,f0;	fld.d	16(B),b2
	d.pfmul.ss a1r,b1r,f0;	fld.d	24(B),d0
	d.m12apm.ss a1i,b1i,f0;	fld.d	24(A),c0
	d.m12apm.ss a2r,b2r,f0;	fld.d	32(A),c1
	d.m12apm.ss a2i,b2i,f0;	fld.d	32(B),d1
	d.m12apm.ss c0r,d0r,f0;	fld.d	40(B),d2
	d.m12apm.ss c0i,d0i,f0;	fld.d	40(A),c2
	d.m12apm.ss c1r,d1r,f0;	fld.d	48(A),a0
	d.m12apm.ss c1i,d1i,f0;	fld.d	48(B),b0
	d.m12apm.ss c2r,d2r,f0;	fld.d	56(B),b1
	d.m12apm.ss c2i,d2i,f0;	fld.d	56(A),a1
	d.m12apm.ss a0r,b0r,f0;	fld.d	64(A),a2
	d.m12apm.ss a0i,b0i,f0;	fld.d	64(B),b2
	d.m12apm.ss a1r,b1r,f0;	fld.d	72(B),d0
	d.m12apm.ss a1i,b1i,f0;	fld.d	72(A),c0
	d.m12apm.ss a2r,b2r,f0;	fld.d	80(A),c1
	d.m12apm.ss a2i,b2i,f0;	fld.d	80(B),d1
	m12apm.ss c0r,d0r,f0;	fld.d	88(B),d2
	m12apm.ss c0i,d0i,f0;	fld.d	88(A),c2
	m12apm.ss c1r,d1r,f0
	m12apm.ss c1i,d1i,f0
	m12apm.ss c2r,d2r,f0
	m12apm.ss c2i,d2i,f0

	//empty multiplier pipe
	m12apm.ss	f0,f0,f0
	m12apm.ss	f0,f0,f0
	m12apm.ss	f0,f0,f0
	//add three numbers in adder pipe
	pfadd.ss f0,f0,a0r
	pfadd.ss f0,f0,a1r
	pfadd.ss a0r,a1r,a2r
	pfadd.ss f0,f0,f0
	pfadd.ss f0,f0,f0
	pfadd.ss f0,f0,a0r
//.if FLOATOPTION=X167
	pfadd.ss a2r,a0r,f0
	pfadd.ss f0,f0,f0
	pfadd.ss f0,f0,f0
	pfadd.ss f0,f0,sum
//.else
//	pfadd.sd a2r,a0r,f0
//	pfadd.sd f0,f0,f0
//	pfadd.sd f0,f0,f0
//	pfadd.sd f0,f0,sum
//.endif
    bri	r1
    nop

.globl	_wvec_rdot
