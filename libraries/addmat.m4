// addmat( su3_matrix *a, su3_matrix *b, su3_matrix *c)
// c <- a+b
// file addmat.m4, i860 assembler version of addmat.c
//
    define(A,r16)	// address of source 1
    define(B,r17)	// address of source 2
    define(C,r18)	// address of result

	.text
	.align	8
_add_su3_matrix:
				fld.d	0(A),f8
				fld.d	8(A),f10
				fld.d	16(A),f12
				fld.d	24(A),f14
				fld.d	32(A),f16
				fld.d	40(A),f18
				fld.d	48(A),f20
				fld.d	56(A),f22
				fld.d	64(A),f24
				fld.d	0(B),f26
				fld.d	8(B),f28
.align 8
	d.pfadd.ss f8,f26,f0;	nop
	d.pfadd.ss f9,f27,f0;	fld.d	16(B),f30
	d.pfadd.ss f10,f28,f0;	nop
	d.pfadd.ss f11,f29,f8;	fld.d	24(B),f26
	d.pfadd.ss f12,f30,f9;	nop
	d.pfadd.ss f13,f31,f10;	fld.d	32(B),f28
	d.pfadd.ss f14,f26,f11;	nop
	d.pfadd.ss f15,f27,f12;	fld.d	40(B),f30
	d.pfadd.ss f16,f28,f13;	nop
	d.pfadd.ss f17,f29,f14;	fld.d	48(B),f26
	d.pfadd.ss f18,f30,f15;	nop
	d.pfadd.ss f19,f31,f16;	fld.d	56(B),f28
	d.pfadd.ss f20,f26,f17;	nop
	d.pfadd.ss f21,f27,f18;	fld.d	64(B),f30
	d.pfadd.ss f22,f28,f19;	fst.d	f8,0(C)
	d.pfadd.ss f23,f29,f20;	fst.d	f10,8(C)
	d.pfadd.ss f24,f30,f21;	fst.d	f12,16(C)
	d.pfadd.ss f25,f31,f22;	fst.d	f14,24(C)
	d.pfadd.ss f0,f0,f23;	fst.d	f16,32(C)
	pfadd.ss f0,f0,f24;	fst.d	f18,40(C)
	pfadd.ss f0,f0,f25;	fst.d	f20,48(C)
				fst.d	f22,56(C)
    bri	r1
				fst.d	f24,64(C)

.globl	_add_su3_matrix
