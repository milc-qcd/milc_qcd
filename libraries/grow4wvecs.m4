//*********************** grow4wvecs.m4 *******************
// i860 assembler version of grow4wvecs.c
//
//void grow_add_four_wvecs( dest, src1, src2, src3, src4, sign, sum)
//half_wilson_vector *src1,*src2,*src3,*src4;
//wilson_vector *dest;
//int sign, sum;

//register usage:
//	f2-3 = source color 0 (all spins)
//	f4-5 = source color 1 (all spins)
//	f6-7 = source color 2 (all spins)
//	f8-9   = dest spin 0 color 0
//	f10-11 = dest spin 0 color 1
//	f12-13 = dest spin 0 color 2
//	f14-15 = dest spin 1 color 0
//	f16-17 = dest spin 1 color 1
//	f18-19 = dest spin 1 color 2
//	f20-21 = dest spin 2 color 0
//	f22-23 = dest spin 2 color 1
//	f24-25 = dest spin 2 color 2
//	f26-27 = dest spin 3 color 0
//	f28-29 = dest spin 3 color 1
//	f30-31 = dest spin 3 color 2

	.text
	.align	8
_grow_add_four_wvecs:
//	    .bf
	// save f2-f7 on stack
	addu	-32,sp,sp
	fst.d	f2,0(sp)
	fst.d	f4,8(sp)
	fst.d	f6,16(sp)

	bte	r0,r22,NOSUM	// if sum==0, don't add in current dest

	adds	-1,r0,r30
	bte	r30,r21,.L4	//if isign== -1, use second code half

// XUP direction - read in and add to dest
				fld.d	0(r17),f2	//src spin 0 color 0
				fld.d	0(r16),f8	//dest spin 0 color 0
				fld.d	8(r17),f4	//src spin 0 color 1
				fld.d	8(r16),f10	//dest spin 0 color 1
.align 8
	d.pfadd.ss f2,f8,f0;	fld.d	16(r17),f6	//src spin 0 color 2
	d.pfadd.ss f3,f9,f0;	fld.d	16(r16),f12	//dest spin 0 color 2
	d.pfadd.ss f4,f10,f0;	fld.d	72(r16),f26	//dest spin 3 color 0
	d.pfadd.ss f5,f11,f8;	fld.d	80(r16),f28	//dest spin 3 color 1
	d.pfadd.ss f6,f12,f9;	fld.d	88(r16),f30	//dest spin 3 color 2
	d.pfadd.ss f7,f13,f10;	nop
	d.pfadd.ss f26,f3,f11;	fld.d	24(r16),f14	//dest spin 1 color 0
	d.pfsub.ss f27,f2,f12;	nop
	d.pfadd.ss f28,f5,f13;	fld.d	24(r17),f2	//src spin 1 color 0
	d.pfsub.ss f29,f4,f26;	nop
	d.pfadd.ss f30,f7,f27;	fld.d	32(r17),f4	//src spin 1 color 1
	d.pfsub.ss f31,f6,f28;	fld.d	32(r16),f16	//dest spin 1 color 1

	d.pfadd.ss f2,f14,f29;	fld.d	40(r17),f6	//src spin 1 color 2
	d.pfadd.ss f3,f15,f30;	fld.d	40(r16),f18	//dest spin 1 color 2
	d.pfadd.ss f4,f16,f31;	fld.d	48(r16),f20	//dest spin 2 color 0
	d.pfadd.ss f5,f17,f14;	fld.d	56(r16),f22	//dest spin 2 color 1
	pfadd.ss f6,f18,f15;	fld.d	64(r16),f24	//dest spin 2 color 2
	pfadd.ss f7,f19,f16;	nop
.align 8	// had 2 single instructions, should be aligned already
	d.pfadd.ss f20,f3,f17
	d.pfsub.ss f21,f2,f18
	d.pfadd.ss f22,f5,f19;	fld.d	0(r18),f2	//src spin 0 color 0
	d.pfsub.ss f23,f4,f20;	nop
	d.pfadd.ss f24,f7,f21;	fld.d	8(r18),f4	//src spin 0 color 1
	d.pfsub.ss f25,f6,f22;	nop
// YUP direction
	pfadd.ss f2,f8,f23;	fld.d	16(r18),f6	//src spin 0 color 2
	pfadd.ss f3,f9,f24;	nop
	pfadd.ss f4,f10,f25	
	pfadd.ss f5,f11,f8
	pfadd.ss f6,f12,f9
	pfadd.ss f7,f13,f10
.align 8	// had 6 single instructions, should be aligned already
	d.pfsub.ss f26,f2,f11
	d.pfsub.ss f27,f3,f12
	d.pfsub.ss f28,f4,f13;	fld.d	24(r18),f2	//src spin 1 color 0
	d.pfsub.ss f29,f5,f26;	nop
	d.pfsub.ss f30,f6,f27;	fld.d	32(r18),f4	//src spin 1 color 1
	d.pfsub.ss f31,f7,f28;	nop

	pfadd.ss f2,f14,f29;	fld.d	40(r18),f6	//src spin 1 color 2
	pfadd.ss f3,f15,f30;	nop
	pfadd.ss f4,f16,f31
	pfadd.ss f5,f17,f14
	pfadd.ss f6,f18,f15
	pfadd.ss f7,f19,f16
.align 8	// had 6 single instructions, should be aligned already
	d.pfadd.ss f20,f2,f17
	d.pfadd.ss f21,f3,f18
	d.pfadd.ss f22,f4,f19;	fld.d	0(r19),f2	//src spin 0 color 0
	d.pfadd.ss f23,f5,f20;	nop
	d.pfadd.ss f24,f6,f21;	fld.d	8(r19),f4	//src spin 0 color 1
	d.pfadd.ss f25,f7,f22;	nop
// ZUP direction
	pfadd.ss f2,f8,f23;	fld.d	16(r19),f6	//src spin 0 color 2
	pfadd.ss f3,f9,f24;	nop
	pfadd.ss f4,f10,f25
	pfadd.ss f5,f11,f8
	pfadd.ss f6,f12,f9
	pfadd.ss f7,f13,f10
.align 8	// had 6 single instructions, should be aligned already
	d.pfadd.ss f20,f3,f11
	d.pfsub.ss f21,f2,f12
	d.pfadd.ss f22,f5,f13;	fld.d	24(r19),f2	//src spin 1 color 0
	d.pfsub.ss f23,f4,f20;	nop
	d.pfadd.ss f24,f7,f21;	fld.d	32(r19),f4	//src spin 1 color 1
	d.pfsub.ss f25,f6,f22;	nop

	pfadd.ss f2,f14,f23;	fld.d	40(r19),f6	//src spin 1 color 2
	pfadd.ss f3,f15,f24;	nop
	pfadd.ss f4,f16,f25
	pfadd.ss f5,f17,f14
	pfadd.ss f6,f18,f15
	pfadd.ss f7,f19,f16
.align 8	// had 6 single instructions, should be aligned already
	d.pfsub.ss f26,f3,f17
	d.pfadd.ss f27,f2,f18
	d.pfsub.ss f28,f5,f19;	fld.d	0(r20),f2	//src spin 0 color 0
	d.pfadd.ss f29,f4,f26;	nop
	d.pfsub.ss f30,f7,f27;	fld.d	8(r20),f4	//src spin 0 color 1
	d.pfadd.ss f31,f6,f28;	nop
// TUP direction
	pfadd.ss f2,f8,f29;	fld.d	16(r20),f6	//src spin 0 color 2
	pfadd.ss f3,f9,f30;	nop
	pfadd.ss f4,f10,f31
	pfadd.ss f5,f11,f8
	pfadd.ss f6,f12,f9
	pfadd.ss f7,f13,f10
.align 8	// had 6 single instructions, should be aligned already
	d.pfadd.ss f20,f2,f11
	d.pfadd.ss f21,f3,f12
	d.pfadd.ss f22,f4,f13;	fld.d	24(r20),f2	//src spin 1 color 0
	d.pfadd.ss f23,f5,f20;	nop
	d.pfadd.ss f24,f6,f21;	fld.d	32(r20),f4	//src spin 1 color 1
	d.pfadd.ss f25,f7,f22;	nop

	d.pfadd.ss f2,f14,f23;	fld.d	40(r20),f6	//src spin 1 color 2
	d.pfadd.ss f3,f15,f24;	fst.d	f8,0(r16)	//dest spin 0 color 0
	d.pfadd.ss f4,f16,f25;	fst.d	f10,8(r16)	//dest spin 0 color 1
	d.pfadd.ss f5,f17,f14;	fst.d	f12,16(r16)	//dest spin 0 color 2
	d.pfadd.ss f6,f18,f15;	fst.d	f20,48(r16)	//dest spin 2 color 0
	d.pfadd.ss f7,f19,f16;	fst.d	f22,56(r16)	//dest spin 2 color 1
	d.pfadd.ss f26,f2,f17;	fst.d	f24,64(r16)	//dest spin 2 color 2
	d.pfadd.ss f27,f3,f18;	fst.d	f14,24(r16)	//dest spin 1 color 0
	d.pfadd.ss f28,f4,f19;	fst.d	f16,32(r16)	//dest spin 1 color 1
	d.pfadd.ss f29,f5,f26;	nop
	d.pfadd.ss f30,f6,f27;	fst.d	f18,40(r16)	//dest spin 1 color 2
	d.pfadd.ss f31,f7,f28;	nop
	d.pfadd.ss f0,f0,f29;	fst.d	f26,72(r16)	//dest spin 3 color 0
	pfadd.ss f0,f0,f30;	nop
	pfadd.ss f0,f0,f31;	fst.d	f28,80(r16)	//dest spin 3 color 1
				fst.d	f30,88(r16)	//dest spin 3 color 2

	//restore stack and exit
	fld.d	16(sp),f6
	fld.d	8(sp),f4
	fld.d	0(sp),f2
	bri	r1
	addu	32,sp,sp


// sign = -1
.L4:
// XDOWN direction:
				fld.d	0(r17),f2	//src spin 0 color 0
				fld.d	0(r16),f8	//dest spin 0 color 0
				fld.d	8(r17),f4	//src spin 0 color 1
				fld.d	8(r16),f10	//dest spin 0 color 1
.align 8
	d.pfadd.ss f2,f8,f0;	fld.d	16(r17),f6	//src spin 0 color 2
	d.pfadd.ss f3,f9,f0;	fld.d	16(r16),f12	//dest spin 0 color 2
	d.pfadd.ss f4,f10,f0;	fld.d	72(r16),f26	//dest spin 3 color 0
	d.pfadd.ss f5,f11,f8;	fld.d	80(r16),f28	//dest spin 3 color 1
	d.pfadd.ss f6,f12,f9;	fld.d	88(r16),f30	//dest spin 3 color 2
	d.pfadd.ss f7,f13,f10;	nop
	d.pfsub.ss f26,f3,f11;	fld.d	24(r16),f14	//dest spin 1 color 0
	d.pfadd.ss f27,f2,f12;	nop
	d.pfsub.ss f28,f5,f13;	fld.d	24(r17),f2	//src spin 1 color 0
	d.pfadd.ss f29,f4,f26;	nop
	d.pfsub.ss f30,f7,f27;	fld.d	32(r17),f4	//src spin 1 color 1
	d.pfadd.ss f31,f6,f28;	fld.d	32(r16),f16	//dest spin 1 color 1
				
	d.pfadd.ss f2,f14,f29;	fld.d	40(r17),f6	//src spin 1 color 2
	d.pfadd.ss f3,f15,f30;	fld.d	40(r16),f18	//dest spin 1 color 2
	d.pfadd.ss f4,f16,f31;	fld.d	48(r16),f20	//dest spin 2 color 0
	d.pfadd.ss f5,f17,f14;	fld.d	56(r16),f22	//dest spin 2 color 1
	pfadd.ss f6,f18,f15;	fld.d	64(r16),f24	//dest spin 2 color 2
	pfadd.ss f7,f19,f16;	nop
.align 8	// had 2 single instructions, should be aligned already
	d.pfsub.ss f20,f3,f17
	d.pfadd.ss f21,f2,f18
	d.pfsub.ss f22,f5,f19;	fld.d	0(r18),f2	//src spin 0 color 0
	d.pfadd.ss f23,f4,f20;	nop
	d.pfsub.ss f24,f7,f21;	fld.d	8(r18),f4	//src spin 0 color 1
	d.pfadd.ss f25,f6,f22;	nop
// YDOWN direction:
	pfadd.ss f2,f8,f23;	fld.d	16(r18),f6	//src spin 0 color 2
	pfadd.ss f3,f9,f24;	nop
	pfadd.ss f4,f10,f25
	pfadd.ss f5,f11,f8
	pfadd.ss f6,f12,f9
	pfadd.ss f7,f13,f10
.align 8	// had 6 single instructions, should be aligned already
	d.pfadd.ss f26,f2,f11
	d.pfadd.ss f27,f3,f12
	d.pfadd.ss f28,f4,f13;	fld.d	24(r18),f2	//src spin 1 color 0
	d.pfadd.ss f29,f5,f26;	nop
	d.pfadd.ss f30,f6,f27;	fld.d	32(r18),f4	//src spin 1 color 1
	d.pfadd.ss f31,f7,f28;	nop

	pfadd.ss f2,f14,f29;	fld.d	40(r18),f6	//src spin 1 color 2
	pfadd.ss f3,f15,f30;	nop
	pfadd.ss f4,f16,f31
	pfadd.ss f5,f17,f14
	pfadd.ss f6,f18,f15
	pfadd.ss f7,f19,f16
.align 8	// had 6 single instructions, should be aligned already
	d.pfsub.ss f20,f2,f17
	d.pfsub.ss f21,f3,f18
	d.pfsub.ss f22,f4,f19;	fld.d	0(r19),f2	//src spin 0 color 0
	d.pfsub.ss f23,f5,f20;	nop
	d.pfsub.ss f24,f6,f21;	fld.d	8(r19),f4	//src spin 0 color 1
	d.pfsub.ss f25,f7,f22;	nop
// ZDOWN direction:
	pfadd.ss f2,f8,f23;	fld.d	16(r19),f6	//src spin 0 color 2
	pfadd.ss f3,f9,f24;	nop
	pfadd.ss f4,f10,f25
	pfadd.ss f5,f11,f8
	pfadd.ss f6,f12,f9
	pfadd.ss f7,f13,f10
.align 8	// had 6 single instructions, should be aligned already
	d.pfsub.ss f20,f3,f11
	d.pfadd.ss f21,f2,f12
	d.pfsub.ss f22,f5,f13;	fld.d	24(r19),f2	//src spin 1 color 0
	d.pfadd.ss f23,f4,f20;	nop
	d.pfsub.ss f24,f7,f21;	fld.d	32(r19),f4	//src spin 1 color 1
	d.pfadd.ss f25,f6,f22;	nop
				
	pfadd.ss f2,f14,f23;	fld.d	40(r19),f6	//src spin 1 color 2
	pfadd.ss f3,f15,f24;	nop
	pfadd.ss f4,f16,f25
	pfadd.ss f5,f17,f14
	pfadd.ss f6,f18,f15
	pfadd.ss f7,f19,f16
.align 8	// had 6 single instructions, should be aligned already
	d.pfadd.ss f26,f3,f17
	d.pfsub.ss f27,f2,f18
	d.pfadd.ss f28,f5,f19;	fld.d	0(r20),f2	//src spin 0 color 0
	d.pfsub.ss f29,f4,f26;	nop
	d.pfadd.ss f30,f7,f27;	fld.d	8(r20),f4	//src spin 0 color 1
	d.pfsub.ss f31,f6,f28;	nop
// TDOWN direction:
	pfadd.ss f2,f8,f29;	fld.d	16(r20),f6	//src spin 0 color 2
	pfadd.ss f3,f9,f30;	nop
	pfadd.ss f4,f10,f31
	pfadd.ss f5,f11,f8
	pfadd.ss f6,f12,f9
	pfadd.ss f7,f13,f10
.align 8	// had 6 single instructions, should be aligned already
	d.pfsub.ss f20,f2,f11
	d.pfsub.ss f21,f3,f12
	d.pfsub.ss f22,f4,f13;	fld.d	24(r20),f2	//src spin 1 color 0
	d.pfsub.ss f23,f5,f20;	nop
	d.pfsub.ss f24,f6,f21;	fld.d	32(r20),f4	//src spin 1 color 1
	d.pfsub.ss f25,f7,f22;	nop
				
	d.pfadd.ss f2,f14,f23;	fld.d	40(r20),f6	//src spin 1 color 2
	d.pfadd.ss f3,f15,f24;	fst.d	f8,0(r16)	//dest spin 0 color 0
	d.pfadd.ss f4,f16,f25;	fst.d	f10,8(r16)	//dest spin 0 color 1
	d.pfadd.ss f5,f17,f14;	fst.d	f12,16(r16)	//dest spin 0 color 2
	d.pfadd.ss f6,f18,f15;	fst.d	f20,48(r16)	//dest spin 2 color 0
	d.pfadd.ss f7,f19,f16;	fst.d	f22,56(r16)	//dest spin 2 color 1
	d.pfsub.ss f26,f2,f17;	fst.d	f24,64(r16)	//dest spin 2 color 2
	d.pfsub.ss f27,f3,f18;	fst.d	f14,24(r16)	//dest spin 1 color 0
	d.pfsub.ss f28,f4,f19;	fst.d	f16,32(r16)	//dest spin 1 color 1
	d.pfsub.ss f29,f5,f26;	nop
	d.pfsub.ss f30,f6,f27;	fst.d	f18,40(r16)	//dest spin 1 color 2
	d.pfsub.ss f31,f7,f28;	nop
	d.pfadd.ss f0,f0,f29;	fst.d	f26,72(r16)	//dest spin 3 color 0
	pfadd.ss f0,f0,f30;	nop
	pfadd.ss f0,f0,f31;	fst.d	f28,80(r16)	//dest spin 3 color 1
				fst.d	f30,88(r16)	//dest spin 3 color 2

	//restore stack and exit
	fld.d	16(sp),f6
	fld.d	8(sp),f4
	fld.d	0(sp),f2
	bri	r1
	addu	32,sp,sp

// Code for not summing in current dest
NOSUM:
	adds	-1,r0,r30
	bte	r30,r21,.L5	//if isign== -1, use second code half

// XUP direction - read in and add to YUP direction
				fld.d	0(r17),f26	//src spin 0 color 0
				fld.d	8(r17),f28	//src spin 0 color 1
				fld.d	16(r17),f30	//src spin 0 color 2
				fld.d	0(r18),f2	//src spin 0 color 0
				fld.d	8(r18),f4	//src spin 0 color 1
// YUP direction
.align 8
	pfadd.ss f2,f26,f0;	nop
	pfadd.ss f3,f27,f0;	fld.d	16(r18),f6	//src spin 0 color 2
	pfadd.ss f4,f28,f0;	nop
	d.pfadd.ss f5,f29,f8;	fld.d	24(r17),f20	//src spin 1 color 0
	d.pfadd.ss f6,f30,f9;	nop
	d.pfadd.ss f7,f31,f10;	fld.d	32(r17),f22	//src spin 1 color 1
	// The second, fourth and sixth adds should be -A-B,  will fix
	// the sign the next time they are used
	d.pfsub.ss f27,f2,f11;	nop
	d.pfadd.ss f26,f3,f12;	fld.d	40(r17),f24	//src spin 1 color 2
	d.pfsub.ss f29,f4,f13;	nop
	d.pfadd.ss f28,f5,f26;	fld.d	24(r18),f2	//src spin 1 color 0
	d.pfsub.ss f31,f6,f27;	nop
	d.pfadd.ss f30,f7,f28;	fld.d	32(r18),f4	//src spin 1 color 1

	pfadd.ss f2,f20,f29;	nop
	pfadd.ss f3,f21,f30;	fld.d	40(r18),f6	//src spin 1 color 2
	pfadd.ss f4,f22,f31	
	pfadd.ss f5,f23,f14	
	pfadd.ss f6,f24,f15
	pfadd.ss f7,f25,f16
.align 8	// had 6 single instructions, should be aligned already
	d.pfadd.ss f21,f2,f17
	d.pfsub.ss f3,f20,f18
	d.pfadd.ss f23,f4,f19;	fld.d	0(r19),f2	//src spin 0 color 0
	d.pfsub.ss f5,f22,f20;	nop
	d.pfadd.ss f25,f6,f21;	fld.d	8(r19),f4	//src spin 0 color 1
	d.pfsub.ss f7,f24,f22;	nop
// ZUP direction
	pfadd.ss f2,f8,f23;	fld.d	16(r19),f6	//src spin 0 color 2
	pfadd.ss f3,f9,f24;	nop
	pfadd.ss f4,f10,f25
	pfadd.ss f5,f11,f8
	pfadd.ss f6,f12,f9
	pfadd.ss f7,f13,f10
.align 8	// had 6 single instructions, should be aligned already
	d.pfadd.ss f20,f3,f11
	d.pfsub.ss f21,f2,f12
	d.pfadd.ss f22,f5,f13;	fld.d	24(r19),f2	//src spin 1 color 0
	d.pfsub.ss f23,f4,f20;	nop
	d.pfadd.ss f24,f7,f21;	fld.d	32(r19),f4	//src spin 1 color 1
	d.pfsub.ss f25,f6,f22;	nop

	pfadd.ss f2,f14,f23;	fld.d	40(r19),f6	//src spin 1 color 2
	pfadd.ss f3,f15,f24;	nop
	pfadd.ss f4,f16,f25
	pfadd.ss f5,f17,f14
	pfadd.ss f6,f18,f15
	pfadd.ss f7,f19,f16
.align 8	// had 6 single instructions, should be aligned already
	d.pfsub.ss f26,f3,f17
	d.pfsub.ss f2,f27,f18		// fixes sign
	d.pfsub.ss f28,f5,f19;	fld.d	0(r20),f2	//src spin 0 color 0
	d.pfsub.ss f4,f29,f26;	nop	// fixes sign
	d.pfsub.ss f30,f7,f27;	fld.d	8(r20),f4	//src spin 0 color 1
	d.pfsub.ss f6,f31,f28;	nop	// fixes sign
// TUP direction
	pfadd.ss f2,f8,f29;	fld.d	16(r20),f6	//src spin 0 color 2
	pfadd.ss f3,f9,f30;	nop
	pfadd.ss f4,f10,f31
	pfadd.ss f5,f11,f8
	pfadd.ss f6,f12,f9
	pfadd.ss f7,f13,f10
.align 8	// had 6 single instructions, should be aligned already
	d.pfadd.ss f20,f2,f11
	d.pfadd.ss f21,f3,f12
	d.pfadd.ss f22,f4,f13;	fld.d	24(r20),f2	//src spin 1 color 0
	d.pfadd.ss f23,f5,f20;	nop
	d.pfadd.ss f24,f6,f21;	fld.d	32(r20),f4	//src spin 1 color 1
	d.pfadd.ss f25,f7,f22;	nop

	d.pfadd.ss f2,f14,f23;	fld.d	40(r20),f6	//src spin 1 color 2
	d.pfadd.ss f3,f15,f24;	fst.d	f8,0(r16)	//dest spin 0 color 0
	d.pfadd.ss f4,f16,f25;	fst.d	f10,8(r16)	//dest spin 0 color 1
	d.pfadd.ss f5,f17,f14;	fst.d	f12,16(r16)	//dest spin 0 color 2
	d.pfadd.ss f6,f18,f15;	fst.d	f20,48(r16)	//dest spin 2 color 0
	d.pfadd.ss f7,f19,f16;	fst.d	f22,56(r16)	//dest spin 2 color 1
	d.pfadd.ss f26,f2,f17;	fst.d	f24,64(r16)	//dest spin 2 color 2
	d.pfadd.ss f27,f3,f18;	fst.d	f14,24(r16)	//dest spin 1 color 0
	d.pfadd.ss f28,f4,f19;	fst.d	f16,32(r16)	//dest spin 1 color 1
	d.pfadd.ss f29,f5,f26;	nop
	d.pfadd.ss f30,f6,f27;	fst.d	f18,40(r16)	//dest spin 1 color 2
	d.pfadd.ss f31,f7,f28;	nop
	d.pfadd.ss f0,f0,f29;	fst.d	f26,72(r16)	//dest spin 3 color 0
	pfadd.ss f0,f0,f30;	nop
	pfadd.ss f0,f0,f31;	fst.d	f28,80(r16)	//dest spin 3 color 1
				fst.d	f30,88(r16)	//dest spin 3 color 2

	//restore stack and exit
	fld.d	16(sp),f6
	fld.d	8(sp),f4
	fld.d	0(sp),f2
	bri	r1
	addu	32,sp,sp


// sign = -1
.L5:
// XDOWN direction:
				fld.d	0(r17),f26	//src spin 0 color 0
				fld.d	8(r17),f28	//src spin 0 color 1
				fld.d	16(r17),f30	//src spin 0 color 2
				fld.d	24(r17),f14	//src spin 1 color 0
				fld.d	32(r17),f16	//src spin 1 color 1
.align 8
	d.pfsub.ss f0,f15,f0;	fld.d	40(r17),f18	//src spin 1 color 2
	d.pfadd.ss f0,f14,f0;	nop
	d.pfsub.ss f0,f17,f0;	fld.d	0(r18),f2	//src spin 0 color 0
	d.pfadd.ss f0,f16,f20;	nop
	d.pfsub.ss f0,f19,f21;	fld.d	8(r18),f4	//src spin 0 color 1
	d.pfadd.ss f0,f18,f22;	nop
// YDOWN direction:
	pfadd.ss f2,f26,f23;	fld.d	16(r18),f6	//src spin 0 color 2
	pfadd.ss f3,f27,f24;	nop
	pfadd.ss f4,f28,f25
	pfadd.ss f5,f29,f8
	pfadd.ss f6,f30,f9
	pfadd.ss f7,f31,f10
.align 8	// had 6 single instructions, should be aligned already
	d.pfsub.ss f2,f27,f11
	d.pfadd.ss f26,f3,f12
	d.pfsub.ss f4,f29,f13;	fld.d	24(r18),f2	//src spin 1 color 0
	d.pfadd.ss f28,f5,f26;	nop
	d.pfsub.ss f6,f31,f27;	fld.d	32(r18),f4	//src spin 1 color 1
	d.pfadd.ss f30,f7,f28;	nop

	pfadd.ss f2,f14,f29;	fld.d	40(r18),f6	//src spin 1 color 2
	pfadd.ss f3,f15,f30;	nop
	pfadd.ss f4,f16,f31
	pfadd.ss f5,f17,f14
	pfadd.ss f6,f18,f15
	pfadd.ss f7,f19,f16
.align 8	// had 6 single instructions, should be aligned already
	d.pfsub.ss f20,f2,f17
	d.pfsub.ss f21,f3,f18
	d.pfsub.ss f22,f4,f19;	fld.d	0(r19),f2	//src spin 0 color 0
	d.pfsub.ss f23,f5,f20;	nop
	d.pfsub.ss f24,f6,f21;	fld.d	8(r19),f4	//src spin 0 color 1
	d.pfsub.ss f25,f7,f22;	nop
// ZDOWN direction:
	pfadd.ss f2,f8,f23;	fld.d	16(r19),f6	//src spin 0 color 2
	pfadd.ss f3,f9,f24;	nop
	pfadd.ss f4,f10,f25
	pfadd.ss f5,f11,f8
	pfadd.ss f6,f12,f9
	pfadd.ss f7,f13,f10
.align 8	// had 6 single instructions, should be aligned already
	d.pfsub.ss f20,f3,f11
	d.pfadd.ss f21,f2,f12
	d.pfsub.ss f22,f5,f13;	fld.d	24(r19),f2	//src spin 1 color 0
	d.pfadd.ss f23,f4,f20;	nop
	d.pfsub.ss f24,f7,f21;	fld.d	32(r19),f4	//src spin 1 color 1
	d.pfadd.ss f25,f6,f22;	nop
				
	pfadd.ss f2,f14,f23;	fld.d	40(r19),f6	//src spin 1 color 2
	pfadd.ss f3,f15,f24;	nop
	pfadd.ss f4,f16,f25
	pfadd.ss f5,f17,f14
	pfadd.ss f6,f18,f15
	pfadd.ss f7,f19,f16
.align 8	// had 6 single instructions, should be aligned already
	d.pfadd.ss f26,f3,f17
	d.pfsub.ss f27,f2,f18
	d.pfadd.ss f28,f5,f19;	fld.d	0(r20),f2	//src spin 0 color 0
	d.pfsub.ss f29,f4,f26;	nop
	d.pfadd.ss f30,f7,f27;	fld.d	8(r20),f4	//src spin 0 color 1
	d.pfsub.ss f31,f6,f28;	nop
// TDOWN direction:
	pfadd.ss f2,f8,f29;	fld.d	16(r20),f6	//src spin 0 color 2
	pfadd.ss f3,f9,f30;	nop
	pfadd.ss f4,f10,f31
	pfadd.ss f5,f11,f8
	pfadd.ss f6,f12,f9
	pfadd.ss f7,f13,f10
.align 8	// had 6 single instructions, should be aligned already
	d.pfsub.ss f20,f2,f11
	d.pfsub.ss f21,f3,f12
	d.pfsub.ss f22,f4,f13;	fld.d	24(r20),f2	//src spin 1 color 0
	d.pfsub.ss f23,f5,f20;	nop
	d.pfsub.ss f24,f6,f21;	fld.d	32(r20),f4	//src spin 1 color 1
	d.pfsub.ss f25,f7,f22;	nop
				
	d.pfadd.ss f2,f14,f23;	fld.d	40(r20),f6	//src spin 1 color 2
	d.pfadd.ss f3,f15,f24;	fst.d	f8,0(r16)	//dest spin 0 color 0
	d.pfadd.ss f4,f16,f25;	fst.d	f10,8(r16)	//dest spin 0 color 1
	d.pfadd.ss f5,f17,f14;	fst.d	f12,16(r16)	//dest spin 0 color 2
	d.pfadd.ss f6,f18,f15;	fst.d	f20,48(r16)	//dest spin 2 color 0
	d.pfadd.ss f7,f19,f16;	fst.d	f22,56(r16)	//dest spin 2 color 1
	d.pfsub.ss f26,f2,f17;	fst.d	f24,64(r16)	//dest spin 2 color 2
	d.pfsub.ss f27,f3,f18;	fst.d	f14,24(r16)	//dest spin 1 color 0
	d.pfsub.ss f28,f4,f19;	fst.d	f16,32(r16)	//dest spin 1 color 1
	d.pfsub.ss f29,f5,f26;	nop
	d.pfsub.ss f30,f6,f27;	fst.d	f18,40(r16)	//dest spin 1 color 2
	d.pfsub.ss f31,f7,f28;	nop
	d.pfadd.ss f0,f0,f29;	fst.d	f26,72(r16)	//dest spin 3 color 0
	pfadd.ss f0,f0,f30;	nop
	pfadd.ss f0,f0,f31;	fst.d	f28,80(r16)	//dest spin 3 color 1
				fst.d	f30,88(r16)	//dest spin 3 color 2

	//restore stack and exit
	fld.d	16(sp),f6
	fld.d	8(sp),f4
	fld.d	0(sp),f2
	bri	r1
	addu	32,sp,sp

//_dest	r16	local
//_src1	r17	local
//_src2	r18	local
//_src3	r19	local
//_src4	r20	local
//_sign	r21	local
//_sum	r22	local

	.text
	.data
	.globl	_grow_add_four_wvecs

	.text
