//************************* wp_shrink4.m4 **********************
//i860 assembler version of wp_shrink4.c
//
// wp_shrink_4dir(a,b1,b2,b3,b4,sign)
//
// Shrink a wilson vector in four directions, producing four
//  half_wilson_vectors.
// void wp_shrink_4dir(a,b1,b2,b3,b4,sign)
// wilson_vector *a; half_wilson_vector *b1,*b2,*b3,*b4;
// int sign;	
// B1 <- (1 +- gamma_x)A,, projection	
//  argument "sign" is sign of gamma matrix.
//  See wp_shrink.c for definitions of gamma matrices and eigenvectors.	
//
    define(A,r16)	// address of source
    define(B1,r17)	//source shrunk with gamma_x
    define(B2,r18)	//source shrunk with gamma_y
    define(B3,r19)	//source shrunk with gamma_z
    define(B4,r20)	//source shrunk with gamma_t
    define(SIGN,r21)
// Use registers f8-f31 for source vector, f2-f7 for results

	.text
	.align	8
_wp_shrink_4dir:

        // save f2-f7 on stack
        addu    -32,sp,sp
        fst.d   f2,0(sp)
        fst.d   f4,8(sp)
        fst.d   f6,16(sp)

        adds    -1,r0,r30
        bte     r30,r21,.L4     //if isign== -1, use second code half

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
	d.pfadd.ss	f15,f20,f6	;	fst.d	f2,0(B1)
	d.pfsub.ss	f16,f23,f7	;	nop
	d.pfadd.ss	f17,f22,f2	;	fst.d	f4,8(B1)
	d.pfsub.ss	f18,f25,f3	;	nop
	d.pfadd.ss	f19,f24,f4	;	fst.d	f6,16(B1)
// YUP direction starts into adder
	d.pfsub.ss	f8,f26,f5	;	nop
	d.pfsub.ss	f9,f27,f6	;	fst.d	f2,24(B1)
	d.pfsub.ss	f10,f28,f7	;	nop
	d.pfsub.ss	f11,f29,f2	;	fst.d	f4,32(B1)
	d.pfsub.ss	f12,f30,f3	;	nop
	d.pfsub.ss	f13,f31,f4	;	fst.d	f6,40(B1)
	d.pfadd.ss	f14,f20,f5	;	nop
	d.pfadd.ss	f15,f21,f6	;	fst.d	f2,0(B2)
	d.pfadd.ss	f16,f22,f7	;	nop
	d.pfadd.ss	f17,f23,f2	;	fst.d	f4,8(B2)
	d.pfadd.ss	f18,f24,f3	;	nop
	d.pfadd.ss	f19,f25,f4	;	fst.d	f6,16(B2)
// ZUP direction
	d.pfsub.ss	f8,f21,f5	;	nop
	d.pfadd.ss	f9,f20,f6	;	fst.d	f2,24(B2)
	d.pfsub.ss	f10,f23,f7	;	nop
	d.pfadd.ss	f11,f22,f2	;	fst.d	f4,32(B2)
	d.pfsub.ss	f12,f25,f3	;	nop
	d.pfadd.ss	f13,f24,f4	;	fst.d	f6,40(B2)
	d.pfadd.ss	f14,f27,f5	;	nop
	d.pfsub.ss	f15,f26,f6	;	fst.d	f2,0(B3)
	d.pfadd.ss	f16,f29,f7	;	nop
	d.pfsub.ss	f17,f28,f2	;	fst.d	f4,8(B3)
	d.pfadd.ss	f18,f31,f3	;	nop
	d.pfsub.ss	f19,f30,f4	;	fst.d	f6,16(B3)
// TUP direction
	d.pfadd.ss	f8,f20,f5	;	nop
	d.pfadd.ss	f9,f21,f6	;	fst.d	f2,24(B3)
	d.pfadd.ss	f10,f22,f7	;	nop
	d.pfadd.ss	f11,f23,f2	;	fst.d	f4,32(B3)
	d.pfadd.ss	f12,f24,f3	;	nop
	d.pfadd.ss	f13,f25,f4	;	fst.d	f6,40(B3)
	d.pfadd.ss	f14,f26,f5	;	nop
	d.pfadd.ss	f15,f27,f6	;	fst.d	f2,0(B4)
	d.pfadd.ss	f16,f28,f7	;	nop
	d.pfadd.ss	f17,f29,f2	;	fst.d	f4,8(B4)
	d.pfadd.ss	f18,f30,f3	;	nop
	d.pfadd.ss	f19,f31,f4	;	fst.d	f6,16(B4)
	pfadd.ss	f0,f0,f5	;	nop
	pfadd.ss	f0,f0,f6	;	fst.d	f2,24(B4)
	pfadd.ss	f0,f0,f7
						fst.d	f4,32(B4)
						fst.d	f6,40(B4)
        //restore stack and exit
        fld.d   16(sp),f6
        fld.d   8(sp),f4
        fld.d   0(sp),f2
        bri     r1
        addu    32,sp,sp

.L4:
//XDOWN direction
						fld.d	0(A),f8
						fld.d	8(A),f10
						fld.d	72(A),f26
						fld.d	80(A),f28
						fld.d	88(A),f30
.align 8
	d.pfadd.ss	f8,f27,f0	;	fld.d	16(A),f12
	d.pfsub.ss	f9,f26,f0	;	fld.d	24(A),f14
	d.pfadd.ss	f10,f29,f0	;	fld.d	48(A),f20
	d.pfsub.ss	f11,f28,f2	;	fld.d	32(A),f16
	d.pfadd.ss	f12,f31,f3	;	fld.d	56(A),f22
	d.pfsub.ss	f13,f30,f4	;	fld.d	40(A),f18
	d.pfadd.ss	f14,f21,f5	;	fld.d	64(A),f24
	d.pfsub.ss	f15,f20,f6	;	fst.d	f2,0(B1)
	d.pfadd.ss	f16,f23,f7	;	nop
	d.pfsub.ss	f17,f22,f2	;	fst.d	f4,8(B1)
	d.pfadd.ss	f18,f25,f3	;	nop
	d.pfsub.ss	f19,f24,f4	;	fst.d	f6,16(B1)
// YDOWN direction
	d.pfadd.ss	f8,f26,f5	;	nop
	d.pfadd.ss	f9,f27,f6	;	fst.d	f2,24(B1)
	d.pfadd.ss	f10,f28,f7	;	nop
	d.pfadd.ss	f11,f29,f2	;	fst.d	f4,32(B1)
	d.pfadd.ss	f12,f30,f3	;	nop
	d.pfadd.ss	f13,f31,f4	;	fst.d	f6,40(B1)
	d.pfsub.ss	f14,f20,f5	;	nop
	d.pfsub.ss	f15,f21,f6	;	fst.d	f2,0(B2)
	d.pfsub.ss	f16,f22,f7	;	nop
	d.pfsub.ss	f17,f23,f2	;	fst.d	f4,8(B2)
	d.pfsub.ss	f18,f24,f3	;	nop
	d.pfsub.ss	f19,f25,f4	;	fst.d	f6,16(B2)
// ZDOWN direction
	d.pfadd.ss	f8,f21,f5	;	nop
	d.pfsub.ss	f9,f20,f6	;	fst.d	f2,24(B2)
	d.pfadd.ss	f10,f23,f7	;	nop
	d.pfsub.ss	f11,f22,f2	;	fst.d	f4,32(B2)
	d.pfadd.ss	f12,f25,f3	;	nop
	d.pfsub.ss	f13,f24,f4	;	fst.d	f6,40(B2)
	d.pfsub.ss	f14,f27,f5	;	nop
	d.pfadd.ss	f15,f26,f6	;	fst.d	f2,0(B3)
	d.pfsub.ss	f16,f29,f7	;	nop
	d.pfadd.ss	f17,f28,f2	;	fst.d	f4,8(B3)
	d.pfsub.ss	f18,f31,f3	;	nop
	d.pfadd.ss	f19,f30,f4	;	fst.d	f6,16(B3)
// TDOWN direction
	d.pfsub.ss	f8,f20,f5	;	nop
	d.pfsub.ss	f9,f21,f6	;	fst.d	f2,24(B3)
	d.pfsub.ss	f10,f22,f7	;	nop
	d.pfsub.ss	f11,f23,f2	;	fst.d	f4,32(B3)
	d.pfsub.ss	f12,f24,f3	;	nop
	d.pfsub.ss	f13,f25,f4	;	fst.d	f6,40(B3)
	d.pfsub.ss	f14,f26,f5	;	nop
	d.pfsub.ss	f15,f27,f6	;	fst.d	f2,0(B4)
	d.pfsub.ss	f16,f28,f7	;	nop
	d.pfsub.ss	f17,f29,f2	;	fst.d	f4,8(B4)
	d.pfsub.ss	f18,f30,f3	;	nop
	d.pfsub.ss	f19,f31,f4	;	fst.d	f6,16(B4)
	pfadd.ss	f0,f0,f5	;	nop
	pfadd.ss	f0,f0,f6	;	fst.d	f2,24(B4)
	pfadd.ss	f0,f0,f7
						fst.d	f4,32(B4)
						fst.d	f6,40(B4)
        //restore stack and exit
        fld.d   0(sp),f2
        fld.d   8(sp),f4
        fld.d   16(sp),f6
        bri     r1
        addu    32,sp,sp

	.globl _wp_shrink_4dir
