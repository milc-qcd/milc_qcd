//************************* wp_shrink8.m4 **********************
//i860 assembler version of wp_shrink8.c
//
// wp_shrink_8dir(a,b,sign)
//
// Shrink a wilson vector in eight directions, producing eight
//  half_wilson_vectors.
// void wp_shrink_8dir(a,b,sign)
// wilson_vector *a; half_wilson_vector *b;
// int sign;	
// B <- (1 +- gamma)A,, projection	, for all directions, up and down
//  argument "sign" is sign of gamma matrix.
//  See wp_shrink.c for definitions of gamma matrices and eigenvectors.	
//
    define(A,r16)	// address of source
    define(B,r17)	// address of dest. - eight half vectors
    define(SIGN,r18)	// if sign = -1, reverse "UP" and "DOWN"
    define(INC,r19)	// +- 48 = amount to increment B
    define(SIZE,48)	// size of half_wilson_vector
    define(SIZE7,336)	// 7 times size of half_wilson_vector
// Use registers f8-f31 for source vector, f2-f7 for results

	.text
	.align	8
_wp_shrink_8dir:
        adds    1,r0,r30	// 1, for testing SIGN

        // save f2-f7 on stack
        addu    -32,sp,sp
        fst.d   f2,0(sp)
        fst.d   f4,8(sp)
        fst.d   f6,16(sp)

	adds	SIZE,r0,INC	// set properly if sign=1
        bte     r30,SIGN,.L4     //if sign = -1, execute next instructions

	subs	r0,INC,INC	// increment if sign = -1
	adds	SIZE7,B,B	// start B out at last half_vector
.L4:

// XUP direction
						fld.d	0(A),f8
						fld.d	8(A),f10
						fld.d	72(A),f26
						fld.d	80(A),f28
						fld.d	88(A),f30
.align 8
	d.pfsub.ss	f8,f27,f0	;	fld.d	16(A),f12
	d.pfadd.ss	f9,f26,f0	;	fld.d	24(A),f14
	d.pfsub.ss	f10,f29,f0	;	fld.d	32(A),f16
	d.pfadd.ss	f11,f28,f2	;	fld.d	40(A),f18
	d.pfsub.ss	f12,f31,f3	;	fld.d	48(A),f20
	d.pfadd.ss	f13,f30,f4	;	fld.d	56(A),f22
	d.pfsub.ss	f14,f21,f5	;	fld.d	64(A),f24
	d.pfadd.ss	f15,f20,f6	;	fst.d	f2,0(B)
	d.pfsub.ss	f16,f23,f7	;	nop
	d.pfadd.ss	f17,f22,f2	;	fst.d	f4,8(B)
	d.pfsub.ss	f18,f25,f3	;	nop
	d.pfadd.ss	f19,f24,f4	;	fst.d	f6,16(B)
// YUP direction starts into adder
	d.pfsub.ss	f8,f26,f5	;	nop
	d.pfsub.ss	f9,f27,f6	;	fst.d	f2,24(B)
	d.pfsub.ss	f10,f28,f7	;	nop
	d.pfsub.ss	f11,f29,f2	;	fst.d	f4,32(B)
	d.pfsub.ss	f12,f30,f3	;	nop
	d.pfsub.ss	f13,f31,f4	;	fst.d	f6,40(B)
	d.pfadd.ss	f14,f20,f5	;	adds 	INC,B,B
	d.pfadd.ss	f15,f21,f6	;	fst.d	f2,0(B)
	d.pfadd.ss	f16,f22,f7	;	nop
	d.pfadd.ss	f17,f23,f2	;	fst.d	f4,8(B)
	d.pfadd.ss	f18,f24,f3	;	nop
	d.pfadd.ss	f19,f25,f4	;	fst.d	f6,16(B)
// ZUP direction
	d.pfsub.ss	f8,f21,f5	;	nop
	d.pfadd.ss	f9,f20,f6	;	fst.d	f2,24(B)
	d.pfsub.ss	f10,f23,f7	;	nop
	d.pfadd.ss	f11,f22,f2	;	fst.d	f4,32(B)
	d.pfsub.ss	f12,f25,f3	;	nop
	d.pfadd.ss	f13,f24,f4	;	fst.d	f6,40(B)
	d.pfadd.ss	f14,f27,f5	;	adds	INC,B,B
	d.pfsub.ss	f15,f26,f6	;	fst.d	f2,0(B)
	d.pfadd.ss	f16,f29,f7	;	nop
	d.pfsub.ss	f17,f28,f2	;	fst.d	f4,8(B)
	d.pfadd.ss	f18,f31,f3	;	nop
	d.pfsub.ss	f19,f30,f4	;	fst.d	f6,16(B)
// TUP direction
	d.pfadd.ss	f8,f20,f5	;	nop
	d.pfadd.ss	f9,f21,f6	;	fst.d	f2,24(B)
	d.pfadd.ss	f10,f22,f7	;	nop
	d.pfadd.ss	f11,f23,f2	;	fst.d	f4,32(B)
	d.pfadd.ss	f12,f24,f3	;	nop
	d.pfadd.ss	f13,f25,f4	;	fst.d	f6,40(B)
	d.pfadd.ss	f14,f26,f5	;	adds	INC,B,B
	d.pfadd.ss	f15,f27,f6	;	fst.d	f2,0(B)
	d.pfadd.ss	f16,f28,f7	;	nop
	d.pfadd.ss	f17,f29,f2	;	fst.d	f4,8(B)
	d.pfadd.ss	f18,f30,f3	;	nop
	d.pfadd.ss	f19,f31,f4	;	fst.d	f6,16(B)
// TDOWN direction
	d.pfsub.ss	f8,f20,f5	;	nop
	d.pfsub.ss	f9,f21,f6	;	fst.d	f2,24(B)
	d.pfsub.ss	f10,f22,f7	;	nop
	d.pfsub.ss	f11,f23,f2	;	fst.d	f4,32(B)
	d.pfsub.ss	f12,f24,f3	;	nop
	d.pfsub.ss	f13,f25,f4	;	fst.d	f6,40(B)
	d.pfsub.ss	f14,f26,f5	;	adds	INC,B,B
	d.pfsub.ss	f15,f27,f6	;	fst.d	f2,0(B)
	d.pfsub.ss	f16,f28,f7	;	nop
	d.pfsub.ss	f17,f29,f2	;	fst.d	f4,8(B)
	d.pfsub.ss	f18,f30,f3	;	nop
	d.pfsub.ss	f19,f31,f4	;	fst.d	f6,16(B)
// ZDOWN direction
	d.pfadd.ss	f8,f21,f5	;	nop
	d.pfsub.ss	f9,f20,f6	;	fst.d	f2,24(B)
	d.pfadd.ss	f10,f23,f7	;	nop
	d.pfsub.ss	f11,f22,f2	;	fst.d	f4,32(B)
	d.pfadd.ss	f12,f25,f3	;	nop
	d.pfsub.ss	f13,f24,f4	;	fst.d	f6,40(B)
	d.pfsub.ss	f14,f27,f5	;	adds	INC,B,B
	d.pfadd.ss	f15,f26,f6	;	fst.d	f2,0(B)
	d.pfsub.ss	f16,f29,f7	;	nop
	d.pfadd.ss	f17,f28,f2	;	fst.d	f4,8(B)
	d.pfsub.ss	f18,f31,f3	;	nop
	d.pfadd.ss	f19,f30,f4	;	fst.d	f6,16(B)
// YDOWN direction
	d.pfadd.ss	f8,f26,f5	;	nop
	d.pfadd.ss	f9,f27,f6	;	fst.d	f2,24(B)
	d.pfadd.ss	f10,f28,f7	;	nop
	d.pfadd.ss	f11,f29,f2	;	fst.d	f4,32(B)
	d.pfadd.ss	f12,f30,f3	;	nop
	d.pfadd.ss	f13,f31,f4	;	fst.d	f6,40(B)
	d.pfsub.ss	f14,f20,f5	;	adds	INC,B,B
	d.pfsub.ss	f15,f21,f6	;	fst.d	f2,0(B)
	d.pfsub.ss	f16,f22,f7	;	nop
	d.pfsub.ss	f17,f23,f2	;	fst.d	f4,8(B)
	d.pfsub.ss	f18,f24,f3	;	nop
	d.pfsub.ss	f19,f25,f4	;	fst.d	f6,16(B)
//XDOWN direction
	d.pfadd.ss	f8,f27,f5	;	nop
	d.pfsub.ss	f9,f26,f6	;	fst.d   f2,24(B)
	d.pfadd.ss	f10,f29,f7	;	nop
	d.pfsub.ss	f11,f28,f2	;	fst.d   f4,32(B)
	d.pfadd.ss	f12,f31,f3	;	nop
	d.pfsub.ss	f13,f30,f4	;	fst.d   f6,40(B)
	d.pfadd.ss	f14,f21,f5	;	adds	INC,B,B
	d.pfsub.ss	f15,f20,f6	;	fst.d	f2,0(B)
	d.pfadd.ss	f16,f23,f7	;	nop
	d.pfsub.ss	f17,f22,f2	;	fst.d	f4,8(B)
	d.pfadd.ss	f18,f25,f3	;	nop
	d.pfsub.ss	f19,f24,f4	;	fst.d	f6,16(B)

	pfadd.ss	f0,f0,f5	;	nop
	pfadd.ss	f0,f0,f6	;	fst.d	f2,24(B)
	pfadd.ss	f0,f0,f7
						fst.d	f4,32(B)
						fst.d	f6,40(B)
        //restore stack and exit
        fld.d   0(sp),f2
        fld.d   8(sp),f4
        fld.d   16(sp),f6
        bri     r1
        addu    32,sp,sp

	.globl _wp_shrink_8dir
