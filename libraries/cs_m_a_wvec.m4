// c_scalar_mult_add_wvec( wilson_vector *a, wilson_vector *b, complex *sp,
//	 wilson_vector *c)
// c <- a + s*b
// file cs_m_a_wvec.m4, i860 version of cs_m_a_wvec.c
//
    define(A,r16)	// address of source 1
    define(B,r17)	// address of source 2
    define(S,r18)	// address of complex scalar
    define(C,r19)	// address of result

    define(zero,f0)	// f0 always is zero
    define(s_real,f8)	// real part of scalar
    define(s_imag,f9)	// real part of scalar

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
_c_scalar_mult_add_wvec:

	// load complex scalar
				fld.l 0(S),s_real
				fld.l 4(S),s_imag
				fld.d 0(B),b0
				fld.d 8(B)++,b1
				fld.d 8(B)++,b2

	// start multiplying
    // spin component 0 starts into multiplier
	.align 8
	d.pfmul.ss s_real,b0r,zero;
					fld.d 0(A),a0
	d.pfmul.ss s_real,b1r,zero;	fld.d 8(A)++,a1
	d.pfmul.ss s_real,b2r,zero;	fld.d 8(A)++,a2
	d.pfadd.ss a0r,zero,zero;	adds -8,C,C
	d.pfadd.ss a1r,zero,zero;	nop;
	d.pfadd.ss a2r,zero,zero;	nop;
	d.m12apm.ss s_imag,b0i,zero;	nop
	d.m12apm.ss s_imag,b1i,zero;	nop
	d.m12apm.ss s_imag,b2i,zero;	nop
	d.m12asm.ss s_real,b0i,zero;	nop
	d.m12asm.ss s_real,b1i,zero;	nop
	d.m12asm.ss s_real,b2i,zero;	adds -1,r0,r20
	d.pfadd.ss a0i,zero,c0r;	or 2,r0,r21
	d.pfadd.ss a1i,zero,c1r;	bla r20,r21,LOOP
	d.pfadd.ss a2i,zero,c2r;	nop
LOOP:	d.m12apm.ss s_imag,b0r,zero;	nop
	d.m12apm.ss s_imag,b1r,zero;	fld.d 8(B)++,b0
	d.m12apm.ss s_imag,b2r,zero;	fld.d 8(B)++,b1
// next spin component starts down pipes
	d.m12apm.ss s_real,b0r,zero;	fld.d 8(B)++,b2
	d.m12apm.ss s_real,b1r,zero;	fld.d 8(A)++,a0
	d.m12apm.ss s_real,b2r,zero;	fld.d 8(A)++,a1
	d.pfadd.ss a0r,zero,c0i;	fld.d 8(A)++,a2
	d.pfadd.ss a1r,zero,c1i;	fst.d c0,8(C)++
	d.pfadd.ss a2r,zero,c2i;	nop
	d.m12apm.ss s_imag,b0i,zero;	fst.d c1,8(C)++
	d.m12apm.ss s_imag,b1i,zero;	nop
	d.m12apm.ss s_imag,b2i,zero;	fst.d c2,8(C)++
	d.m12asm.ss s_real,b0i,zero;	nop
	d.m12asm.ss s_real,b1i,zero;	nop
	d.m12asm.ss s_real,b2i,zero;	nop
	d.pfadd.ss a0i,zero,c0r;	nop
	d.pfadd.ss a1i,zero,c1r;	bla r20,r21,LOOP
	d.pfadd.ss a2i,zero,c2r;	nop	// last instruction in loop
	d.m12apm.ss s_imag,b0r,zero;	nop
	d.m12apm.ss s_imag,b1r,zero;	nop
	d.m12apm.ss s_imag,b2r,zero;	nop
//empty pipes
	d.m12apm.ss s_real,b0r,zero;	nop
	d.m12apm.ss s_real,b1r,zero;	nop
	d.m12apm.ss s_real,b2r,zero;	nop
	d.pfadd.ss zero,zero,c0i;	nop
	pfadd.ss zero,zero,c1i;	fst.d c0,8(C)++
	pfadd.ss zero,zero,c2i;	fst.d c1,8(C)++
        			bri     r1
				fst.d c2,8(C)++

.globl	_c_scalar_mult_add_wvec
