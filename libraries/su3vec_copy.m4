// su3vec_copy( su3_vector *a, su3_vector *b )
// b <- a
// file su3vec_copy.m4, i860 assembler version of su3vec_copy.c
//
    define(A,r16)	// address of source 
    define(B,r17)	// address of result

    define(a0,f8)	// complex number = register pair
    define(a0r,f8)	// real part
    define(a0i,f9)	// imag part
    define(a1,f10)
    define(a1r,f10)
    define(a1i,f11)
    define(a2,f12)
    define(a2r,f12)
    define(a2i,f13)

	.text
	.align	8
_su3vec_copy:
	fld.d	0(A),a0
	fld.d	8(A),a1
	fld.d	16(A),a2
	fst.d	a0,0(B)
	fst.d	a1,8(B)
bri r1
	fst.d	a2,16(B)
.globl	_su3vec_copy
