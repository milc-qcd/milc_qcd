// float magsq_wvec( wilson_vector *a );
// return squared magnitude  of A
// file msq_wvec.m4, i860 assembler version of msq_wvec.c
//
    define(A,r16)	// address of source 
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

    define(c0,f20)	// complex number = register pair
    define(c0r,f20)	// real part
    define(c0i,f21)	// imag part
    define(c1,f22)
    define(c1r,f22)
    define(c1i,f23)
    define(c2,f24)
    define(c2r,f24)
    define(c2i,f25)

	.text
	.align	8
_magsq_wvec:
	// enter zeroes into adder pipeline
	d.pfadd.ss	f0,f0,f0;	fld.d	0(A),a0
	d.pfadd.ss	f0,f0,f0;	nop
	d.pfadd.ss	f0,f0,f0;	fld.d	8(A),a1
	d.pfmul.ss a0r,a0r,f0;	nop
	d.pfmul.ss a0i,a0i,f0;	fld.d	16(A),a2
	d.pfmul.ss a1r,a1r,f0;	nop
	d.m12apm.ss a1i,a1i,f0;	fld.d	24(A),c0
	d.m12apm.ss a2r,a2r,f0;	nop
	d.m12apm.ss a2i,a2i,f0;	fld.d	32(A),c1
	d.m12apm.ss c0r,c0r,f0;	nop
	d.m12apm.ss c0i,c0i,f0;	fld.d	40(A),c2
	d.m12apm.ss c1r,c1r,f0;	nop
	d.m12apm.ss c1i,c1i,f0;	fld.d	48(A),a0
	d.m12apm.ss c2r,c2r,f0;	nop
	d.m12apm.ss c2i,c2i,f0;	fld.d	56(A),a1
	d.m12apm.ss a0r,a0r,f0;	nop
	d.m12apm.ss a0i,a0i,f0;	fld.d	64(A),a2
	d.m12apm.ss a1r,a1r,f0;	nop
	d.m12apm.ss a1i,a1i,f0;	fld.d	72(A),c0
	d.m12apm.ss a2r,a2r,f0;	nop
	d.m12apm.ss a2i,a2i,f0;	fld.d	80(A),c1
	m12apm.ss c0r,c0r,f0;	nop
	m12apm.ss c0i,c0i,f0;	fld.d	88(A),c2
	m12apm.ss c1r,c1r,f0
	m12apm.ss c1i,c1i,f0
	m12apm.ss c2r,c2r,f0
	m12apm.ss c2i,c2i,f0

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

.globl	_magsq_wvec
