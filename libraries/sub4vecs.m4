// sub_four_su3_vecs( su3_vector *a,*b1,*b2,*b3,*b4 )
// a <- a - b1 -b2 -b3 -b4
// file sub4vecs.m4, i860 assembler version of sub4vecs.c
//
    define(A,r16)	// address of source/dest
    define(B1,r17)	// address of first vector to be subtracted
    define(B2,r18)	// address of second vector to be subtracted
    define(B3,r19)	// address of third vector to be subtracted
    define(B4,r20)	// address of fourth vector to be subtracted

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

	.text
	.align	8
_sub_four_su3_vecs:
				pfld.d	0(A),f0
				pfld.d	8(A),f0
				pfld.d	16(A),f0
				pfld.d 0(B1),a0
				pfld.d 8(B1),a1
				pfld.d 16(B1),a2
				pfld.d	0(B2),b0
	.align 8
	d.pfadd.ss f0,f0,f0;	// dummy - start dual intructions
				pfld.d	8(B2),b1
	d.pfsub.ss a0r,b0r,f0;	nop
	d.pfsub.ss a0i,b0i,f0;	pfld.d	16(B2),b2
	d.pfsub.ss a1r,b1r,f0;	nop
	d.pfsub.ss a1i,b1i,a0r;	pfld.d	0(B3),b0
	d.pfsub.ss a2r,b2r,a0i;	nop
	d.pfsub.ss a2i,b2i,a1r;	pfld.d	8(B3),b1
	d.pfsub.ss a0r,b0r,a1i;	nop
	d.pfsub.ss a0i,b0i,a2r;	pfld.d	16(B3),b2
	d.pfsub.ss a1r,b1r,a2i;	nop
	d.pfsub.ss a1i,b1i,a0r;	pfld.d	0(B4),b0
	d.pfsub.ss a2r,b2r,a0i;	nop
	d.pfsub.ss a2i,b2i,a1r;	pfld.d	8(B4),b1
	d.pfsub.ss a0r,b0r,a1i;	nop
	d.pfsub.ss a0i,b0i,a2r;	pfld.d	16(B4),b2
	d.pfsub.ss a1r,b1r,a2i;	nop
	d.pfsub.ss a1i,b1i,a0r;	pfld.d	0(B4),b0
	d.pfsub.ss a2r,b2r,a0i;	nop
	d.pfsub.ss a2i,b2i,a1r;	pfld.d	8(B4),b1
	d.pfsub.ss a0r,b0r,a1i;	nop
	d.pfsub.ss a0i,b0i,a2r;	pfld.d	16(B4),b2
	d.pfsub.ss a1r,b1r,a2i;	nop
	d.pfsub.ss a1i,b1i,a0r;	nop
	d.pfsub.ss a2r,b2r,a0i;	nop
	d.pfsub.ss a2i,b2i,a1r;	fst.d	a0,0(A)
	pfsub.ss f0,f0,a1i;	nop
	pfsub.ss f0,f0,a2r;	fst.d	a1,8(A)
	pfsub.ss f0,f0,a2i
    bri	r1
				fst.d	a2,16(A)

.globl	_sub_four_su3_vecs
