// scalar_mult_add_wvec( wilson_vector *a, wilson_vector *b, float s,
//	 wilson_vector *c)
// c <- a + s*b
// file s_m_a_wvec.m4, i860 version of s_m_a_wvec.c
//
    define(A,r16)	// address of source 1
    define(B,r17)	// address of source 2
    define(C,r18)	// address of result
    //S = f8,f9	// scalar, passed as double unless X167 option is used

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

    define(d0,f28)
    define(d0r,f28)
    define(d0i,f29)
    define(d1,f30)
    define(d1r,f30)
    define(d1i,f31)
    define(d2,f8)
    define(d2r,f8)
    define(d2i,f9)

    define(e0,f2)
    define(e0r,f2)
    define(e0i,f3)
    define(e1,f4)
    define(e1r,f4)
    define(e1i,f5)
    define(e2,f6)
    define(e2r,f6)
    define(e2i,f7)

	.text
	.align	8
_scalar_mult_add_wvec:

	// save f2-f7 on stack
        addu    -32,sp,sp
        fst.d   f2,0(sp)
        fst.d   f4,8(sp)
        fst.d   f6,16(sp)
//.if FLOATOPTION=X167
        // clear KR and KI - double precision (prevents traps)
				fld.d 0(B),b0
	.align 8
	d.r2pt.dd f0,f0,f0	
				fld.d 8(B),b1
	d.r2pt.dd f0,f0,f0;	fld.d 16(B),b2
//.else
	// convert double to single
        // clear KR and KI - double precision (prevents traps)
//	.align 8
//	d.r2pt.dd f0,f0,f0
//	d.i2pt.dd f0,f0,f0
//	d.pfmov.ds f8,f0; 	fld.d 0(B),b0
//	d.pfmov.ds f0,f0;	fld.d 8(B),b1
//	d.pfmov.ds f0,f0; 	nop
//	d.pfmov.ds f0,f8;	fld.d 16(B),b2
//.endif
	// load KR with scalar, single precision as it comes out of
	// adder pipe
	d.r2pt.ss f8,f0,f0;	nop
	// start multiplying
	d.r2p1.ss f0,b0r,f0;	fld.d 0(A),a0
	d.r2p1.ss f0,b0i,f0;	fld.d 8(A),a1
	d.r2p1.ss f0,b1r,f0;	fld.d 16(A),a2
	d.r2p1.ss a0r,b1i,f0;	fld.d 24(B),b0
	d.r2p1.ss a0i,b2r,f0;	nop
	d.r2p1.ss a1r,b2i,f0;	fld.d 32(B),b1
    // spin component 1 starts into multiplier
	d.r2p1.ss a1i,b0r,c0r;	fld.d 40(B),b2
	d.r2p1.ss a2r,b0i,c0i;	fld.d 24(A),a0
	d.r2p1.ss a2i,b1r,c1r;	fld.d 32(A),a1
	d.r2p1.ss a0r,b1i,c1i;	fld.d 48(B),b0
	d.r2p1.ss a0i,b2r,c2r;	fld.d 40(A),a2
	d.r2p1.ss a1r,b2i,c2i;	fld.d 56(B),b1
    // spin component 2 starts into multiplier
	d.r2p1.ss a1i,b0r,d0r;	fld.d 64(B),b2
	d.r2p1.ss a2r,b0i,d0i;	fld.d 48(A),a0
	d.r2p1.ss a2i,b1r,d1r;	fld.d 56(A),a1
	d.r2p1.ss a0r,b1i,d1i;	fld.d 72(B),b0
	d.r2p1.ss a0i,b2r,d2r;	fld.d 64(A),a2
	d.r2p1.ss a1r,b2i,d2i;	fld.d 80(B),b1
    // spin component 3 starts into multiplier
	d.r2p1.ss a1i,b0r,e0r;	fld.d 88(B),b2
	d.r2p1.ss a2r,b0i,e0i;	fld.d 72(A),a0
	d.r2p1.ss a2i,b1r,e1r;	fld.d 80(A),a1
	d.r2p1.ss a0r,b1i,e1i;	fld.d 88(A),a2
	d.r2p1.ss a0i,b2r,e2r;	fst.d c0,0(C)
	d.r2p1.ss a1r,b2i,e2i;	fst.d c1,8(C)
	d.r2p1.ss a1i,f0,a0r;	fst.d c2,16(C)
	d.r2p1.ss a2r,f0,a0i;	fst.d d0,24(C)
	d.r2p1.ss a2i,f0,a1r;	fst.d d1,32(C)
	d.r2pt.ss f0,f0,a1i;	fst.d d2,40(C)	//clear KR register
	r2p1.ss f0,f0,a2r;	fst.d e0,48(C)
	r2p1.ss f0,f0,a2i;	fst.d e1,56(C)
				fst.d e2,64(C)
				fst.d a0,72(C)
				fst.d a1,80(C)
				fst.d a2,88(C)
	//restore stack and exit
        fld.d   16(sp),f6
        fld.d   8(sp),f4
        fld.d   0(sp),f2
        bri     r1
        addu    32,sp,sp

.globl	_scalar_mult_add_wvec
