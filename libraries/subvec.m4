// subvec( su3_vector *a, su3_vector *b, su3_vector *c)
// c <- a+b
// file subvec.m4, i860 assembler version of subvec.c
//
    define(A,r16)	// address of source 1
    define(B,r17)	// address of source 2
    define(C,r18)	// address of result

    define(a0,f8)	// complex number = register pair
    define(a0r,f8)	// real part
    define(a0i,f9)	// imag part
    define(a1,f10)
    define(a1r,f10)
    define(a1i,f11)
    define(a2,f12)
    define(a2r,f12)
    define(a2i,f13)

    define(b0,f14)
    define(b0r,f14)
    define(b0i,f15)
    define(b1,f16)
    define(b1r,f16)
    define(b1i,f17)
    define(b2,f18)
    define(b2r,f18)
    define(b2i,f19)

    define(c0,f20)
    define(c0r,f20)
    define(c0i,f21)
    define(c1,f22)
    define(c1r,f22)
    define(c1i,f23)
    define(c2,f24)
    define(c2r,f24)
    define(c2i,f25)

	.text
	.align	8
_sub_su3_vector:
				fld.d	0(A),a0
				fld.d	0(B),b0
				fld.d	8(A),a1
				fld.d	8(B),b1
	pfsub.ss a0r,b0r,f0
				fld.d	16(A),a2
	pfsub.ss a0i,b0i,f0
				fld.d	16(B),b2
	pfsub.ss a1r,b1r,f0
.align 8
	d.pfsub.ss a1i,b1i,c0r
	d.pfsub.ss a2r,b2r,c0i
	d.pfsub.ss a2i,b2i,c1r;	fst.d	c0,0(C)
	pfsub.ss f0,f0,c1i;	nop
	pfsub.ss f0,f0,c2r;	fst.d	c1,8(C)
	pfsub.ss f0,f0,c2i
    bri	r1
				fst.d	c2,16(C)

.globl	_sub_su3_vector
