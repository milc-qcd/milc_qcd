// scalar_mult_sum_su3_vector( su3_vector *a, su3_vector *b, float s)
// a <- a + s*b
// file s_m_sum_vec.m4, i860 version of s_m_sum_vec.c
//
    define(A,r16)	// address of source 1
    define(B,r17)	// address of source 2
    //S = f8,f9	// scalar, passed as double

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
_scalar_mult_sum_su3_vector:

//.if FLOATOPTION=X167
        // clear KR and KI - double precision (prevents traps)
				fld.d 0(B),b0
	.align 8
	d.r2pt.ss f0,f0,f0	
				fld.d 8(B),b1
	d.r2pt.ss f0,f0,f0;	fld.d 16(B),b2
//.else
	// convert double to single
        // clear KR and KI - double precision (prevents traps
//	.align 8
//	d.r2pt.dd f0,f0,f0
//	d.i2pt.dd f0,f0,f0
//	d.pfmov.ds f8,f0;	fld.d 0(B),b0
//	d.pfmov.ds f0,f0;	fld.d 8(B),b1
//	d.pfmov.ds f0,f0;	nop
//	d.pfmov.ds f0,f8;	fld.d 16(B),b2
//.endif
	// load KR with scalar, single precision as it comes out of
	// adder pipe
	d.r2pt.ss f8,f0,f0;	nop
	// start multiplying
	d.r2p1.ss f0,b0r,f0;	fld.d 0(A),a0
	d.r2p1.ss f0,b0i,f0;	nop
	d.r2p1.ss f0,b1r,f0;	fld.d 8(A),a1
	d.r2p1.ss a0r,b1i,f0;	nop
	d.r2p1.ss a0i,b2r,f0;	fld.d 16(A),a2
	d.r2p1.ss a1r,b2i,f0;	nop
	d.r2p1.ss a1i,f0,c0r;	nop
	d.r2p1.ss a2r,f0,c0i;	nop
	d.r2p1.ss a2i,f0,c1r;	fst.d c0,0(A)
	r2pt.ss f0,f0,c1i;	nop	// clear KR register
	r2p1.ss f0,f0,c2r;	fst.d c1,8(A)
	r2p1.ss f0,f0,c2i
    bri	r1
				fst.d c2,16(A)

.globl	_scalar_mult_sum_su3_vector
