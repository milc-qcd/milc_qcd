// complex wvec_dot(a,b) wilson_vector *a,*b;
// 
// file wvec_dot.m4, i860 version of wvec_dot.c
//
// Use two thirds of three cycle accumulator pipe - one for real an
// one for imaginary part.  Every third operation is a dummy to push pipe.
    define(A,r17)	// address of source 1
    define(B,r18)	// address of source 2
    define(S,r16)	// address of complex scalar result

    define(zero,f0)	// f0 always is zero

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

	.text
	.align	8
_wvec_dot:

	.align 8
    // clear adder pipe
	d.pfadd.ss zero,zero,zero;
					fld.d 0(A),a0
	d.pfadd.ss zero,zero,zero;	fld.d 0(B),b0
	d.pfadd.ss zero,zero,zero;	nop
    // spin component 0 starts into multiplier
	d.pfmul.ss a0r,b0r,zero;	nop
	d.pfmul.ss a0r,b0i,zero;	fld.d 8(B)++,b1
	d.pfmul.ss f0,f0,f0;		adds -1,r0,r20	// dummy operation
	d.m12apm.ss a0i,b0i,zero;	fld.d 8(A)++,a1
	d.m12apm.ss a0i,b0r,zero;	or 2,r0,r21
	d.m12apm.ss f0,f0,f0;		fld.d 8(A)++,a2
	d.m12apm.ss a1r,b1r,zero;	bla r20,r21,LOOP
	d.m12asm.ss a1r,b1i,zero;	fld.d 8(B)++,b2
LOOP:	d.m12apm.ss f0,f0,f0;		nop
	d.m12apm.ss a1i,b1i,zero;	fld.d 8(B)++,b0
	d.m12apm.ss a1i,b1r,zero;	nop
	d.m12apm.ss f0,f0,f0;		fld.d 8(B)++,b1
	d.m12apm.ss a2r,b2r,zero;	nop
	d.m12asm.ss a2r,b2i,zero;	fld.d 8(A)++,a0
	d.m12apm.ss f0,f0,f0;		nop
	d.m12apm.ss a2i,b2i,zero;	fld.d 8(A)++,a1
	d.m12apm.ss a2i,b2r,zero;	nop
	d.m12apm.ss f0,f0,f0;		fld.d 8(A)++,a2
    // next spin component starts into multiplier
	d.m12apm.ss a0r,b0r,zero;	nop
	d.m12asm.ss a0r,b0i,zero;	fld.d 8(B)++,b2
	d.m12apm.ss f0,f0,f0;		nop
	d.m12apm.ss a0i,b0i,zero;	nop
	d.m12apm.ss a0i,b0r,zero;	nop
	d.m12apm.ss f0,f0,f0;		nop
	d.m12apm.ss a1r,b1r,zero;	bla r20,r21,LOOP
	d.m12asm.ss a1r,b1i,zero;	nop
	d.m12apm.ss f0,f0,f0;		nop
	d.m12apm.ss a1i,b1i,zero;	nop
	d.m12apm.ss a1i,b1r,zero;	nop
	d.m12apm.ss f0,f0,f0;		nop
	d.m12apm.ss a2r,b2r,zero;	nop
	d.m12asm.ss a2r,b2i,zero;	nop
	d.m12apm.ss f0,f0,f0;		nop
	d.m12apm.ss a2i,b2i,zero;	nop
	d.m12apm.ss a2i,b2r,zero;	nop
	d.m12apm.ss f0,f0,f0;		nop
   // empty pipelines
	d.m12apm.ss f0,f0,f0;		nop
	d.m12asm.ss f0,f0,f0;		nop
	d.m12apm.ss f0,f0,f0;		nop
	pfadd.ss f0,f0,f8;		nop	// real part of answer
	pfadd.ss f0,f0,f9;		fst.l f8,0(S)
					bri r1
					fst.l f9,4(S)

.globl _wvec_dot
