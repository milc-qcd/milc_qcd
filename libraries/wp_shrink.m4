//************************ wp_shrink.m4 **********************
//i860 assembler version of wp_shrink.c
//
// ccom  -OLM -X22 -X70 -X74 -X80 -X83 -X85 -X104 -X247 -X254 -X266 -X278 -X325
//	 -X326 -X350 -X370 -X383 -X393 -X412 -X424 -X494 -X501 -X523 -X524
//	 -X525 -X543 -X616 -X618 -X630
    define(SRC,r16)
    define(DEST,r17)
    define(DIR,r18)
    define(SIGN,r19)

	.text
	.align	4
	.align	8
_wp_shrink:
//	    .bf

	adds	-16,sp,sp
	st.l	r1,12(sp)
	// if sign = -1, dir <- 7-dir
	adds	-1,r0,r30
	btne	r30,SIGN,.L4
	subs	7,DIR,DIR
.L4:

	subu	7,DIR,r0
	bnc	.L55
	shl	2,DIR,DIR
	orh	h%.L153,r0,r30
	or	l%.L153,r30,r30
	adds	r30,DIR,r30
	ld.l	0(r30),r30
	bri	r30
	 nop

.L7:	//case XUP
	fld.d	0(SRC),f8
	fld.d	72(SRC),f14
	fld.d	8(SRC),f10
	fld.d	80(SRC),f16
.align 8
	d.pfsub.ss	f8,f15,f0	;	fld.d	16(SRC),f12
	d.pfadd.ss	f9,f14,f0	;	fld.d	88(SRC),f18
	d.pfsub.ss	f10,f17,f0	;	fld.d	24(SRC),f8
	d.pfadd.ss	f11,f16,f20	;	fld.d	48(SRC),f14
	d.pfsub.ss	f12,f19,f21	;	fld.d	32(SRC),f10
	d.pfadd.ss	f13,f18,f22	;	fld.d	56(SRC),f16
	d.pfsub.ss	f8,f15,f23	;	fld.d	40(SRC),f12
	d.pfadd.ss	f9,f14,f24	;	fld.d	64(SRC),f18
	d.pfsub.ss	f10,f17,f25	;	fst.d	f20,0(DEST)
	d.pfadd.ss	f11,f16,f26	;	fst.d	f22,8(DEST)
	d.pfsub.ss	f12,f19,f27	;	fst.d	f24,16(DEST)
	d.pfadd.ss	f13,f18,f28	;	fst.d	f26,24(DEST)
	pfadd.ss	f0,f0,f29	;	nop
	pfadd.ss	f0,f0,f30	;	fst.d	f28,32(DEST)
	pfadd.ss	f0,f0,f31
	fst.d	f30,40(DEST)
	ld.l	12(sp),r1
	adds	16,sp,sp
	bri	r1
	 nop
.L13:	//case XDOWN

	fld.d	0(SRC),f8
	fld.d	72(SRC),f14
	fld.d	8(SRC),f10
	fld.d	80(SRC),f16
.align 8
	d.pfadd.ss	f8,f15,f0	;	fld.d	16(SRC),f12
	d.pfsub.ss	f9,f14,f0	;	fld.d	88(SRC),f18
	d.pfadd.ss	f10,f17,f0	;	fld.d	24(SRC),f8
	d.pfsub.ss	f11,f16,f20	;	fld.d	48(SRC),f14
	d.pfadd.ss	f12,f19,f21	;	fld.d	32(SRC),f10
	d.pfsub.ss	f13,f18,f22	;	fld.d	56(SRC),f16
	d.pfadd.ss	f8,f15,f23	;	fld.d	40(SRC),f12
	d.pfsub.ss	f9,f14,f24	;	fld.d	64(SRC),f18
	d.pfadd.ss	f10,f17,f25	;	fst.d	f20,0(DEST)
	d.pfsub.ss	f11,f16,f26	;	fst.d	f22,8(DEST)
	d.pfadd.ss	f12,f19,f27	;	fst.d	f24,16(DEST)
	d.pfsub.ss	f13,f18,f28	;	fst.d	f26,24(DEST)
	pfadd.ss	f0,f0,f29	;	nop
	pfadd.ss	f0,f0,f30	;	fst.d	f28,32(DEST)
	pfadd.ss	f0,f0,f31
	fst.d	f30,40(DEST)

	ld.l	12(sp),r1
	adds	16,sp,sp
	bri	r1
	 nop
.L19:	//case YUP

	fld.d	0(SRC),f8
	fld.d	72(SRC),f14
	fld.d	8(SRC),f10
	fld.d	80(SRC),f16
.align 8
	d.pfsub.ss	f8,f14,f0	;	fld.d	16(SRC),f12
	d.pfsub.ss	f9,f15,f0	;	fld.d	88(SRC),f18
	d.pfsub.ss	f10,f16,f0	;	fld.d	24(SRC),f8
	d.pfsub.ss	f11,f17,f20	;	fld.d	48(SRC),f14
	d.pfsub.ss	f12,f18,f21	;	fld.d	32(SRC),f10
	d.pfsub.ss	f13,f19,f22	;	fld.d	56(SRC),f16
	d.pfadd.ss	f8,f14,f23	;	fld.d	40(SRC),f12
	d.pfadd.ss	f9,f15,f24	;	fld.d	64(SRC),f18
	d.pfadd.ss	f10,f16,f25	;	fst.d	f20,0(DEST)
	d.pfadd.ss	f11,f17,f26	;	fst.d	f22,8(DEST)
	d.pfadd.ss	f12,f18,f27	;	fst.d	f24,16(DEST)
	d.pfadd.ss	f13,f19,f28	;	fst.d	f26,24(DEST)
	pfadd.ss	f0,f0,f29	;	nop
	pfadd.ss	f0,f0,f30	;	fst.d	f28,32(DEST)
	pfadd.ss	f0,f0,f31
	fst.d	f30,40(DEST)

	ld.l	12(sp),r1
	adds	16,sp,sp
	bri	r1
	 nop
.L25:	//case YDOWN
	fld.d	0(SRC),f8
	fld.d	72(SRC),f14
	fld.d	8(SRC),f10
	fld.d	80(SRC),f16
.align 8
	d.pfadd.ss	f8,f14,f0	;	fld.d	16(SRC),f12
	d.pfadd.ss	f9,f15,f0	;	fld.d	88(SRC),f18
	d.pfadd.ss	f10,f16,f0	;	fld.d	24(SRC),f8
	d.pfadd.ss	f11,f17,f20	;	fld.d	48(SRC),f14
	d.pfadd.ss	f12,f18,f21	;	fld.d	32(SRC),f10
	d.pfadd.ss	f13,f19,f22	;	fld.d	56(SRC),f16
	d.pfsub.ss	f8,f14,f23	;	fld.d	40(SRC),f12
	d.pfsub.ss	f9,f15,f24	;	fld.d	64(SRC),f18
	d.pfsub.ss	f10,f16,f25	;	fst.d	f20,0(DEST)
	d.pfsub.ss	f11,f17,f26	;	fst.d	f22,8(DEST)
	d.pfsub.ss	f12,f18,f27	;	fst.d	f24,16(DEST)
	d.pfsub.ss	f13,f19,f28	;	fst.d	f26,24(DEST)
	pfadd.ss	f0,f0,f29	;	nop
	pfadd.ss	f0,f0,f30	;	fst.d	f28,32(DEST)
	pfadd.ss	f0,f0,f31
	fst.d	f30,40(DEST)

	ld.l	12(sp),r1
	adds	16,sp,sp
	bri	r1
	 nop
.L31:	//case ZUP
	fld.d	0(SRC),f8
	fld.d	48(SRC),f14
	fld.d	8(SRC),f10
	fld.d	56(SRC),f16
.align 8
	d.pfsub.ss	f8,f15,f0	;	fld.d	16(SRC),f12
	d.pfadd.ss	f9,f14,f0	;	fld.d	64(SRC),f18
	d.pfsub.ss	f10,f17,f0	;	fld.d	24(SRC),f8
	d.pfadd.ss	f11,f16,f20	;	fld.d	72(SRC),f14
	d.pfsub.ss	f12,f19,f21	;	fld.d	32(SRC),f10
	d.pfadd.ss	f13,f18,f22	;	fld.d	80(SRC),f16
	d.pfadd.ss	f8,f15,f23	;	fld.d	40(SRC),f12
	d.pfsub.ss	f9,f14,f24	;	fld.d	88(SRC),f18
	d.pfadd.ss	f10,f17,f25	;	fst.d	f20,0(DEST)
	d.pfsub.ss	f11,f16,f26	;	fst.d	f22,8(DEST)
	d.pfadd.ss	f12,f19,f27	;	fst.d	f24,16(DEST)
	d.pfsub.ss	f13,f18,f28	;	fst.d	f26,24(DEST)
	pfadd.ss	f0,f0,f29	;	nop
	pfadd.ss	f0,f0,f30	;	fst.d	f28,32(DEST)
	pfadd.ss	f0,f0,f31
	fst.d	f30,40(DEST)

	ld.l	12(sp),r1
	adds	16,sp,sp
	bri	r1
	 nop
.L37:	//case ZDOWN
	fld.d	0(SRC),f8
	fld.d	48(SRC),f14
	fld.d	8(SRC),f10
	fld.d	56(SRC),f16
.align 8
	d.pfadd.ss	f8,f15,f0	;	fld.d	16(SRC),f12
	d.pfsub.ss	f9,f14,f0	;	fld.d	64(SRC),f18
	d.pfadd.ss	f10,f17,f0	;	fld.d	24(SRC),f8
	d.pfsub.ss	f11,f16,f20	;	fld.d	72(SRC),f14
	d.pfadd.ss	f12,f19,f21	;	fld.d	32(SRC),f10
	d.pfsub.ss	f13,f18,f22	;	fld.d	80(SRC),f16
	d.pfsub.ss	f8,f15,f23	;	fld.d	40(SRC),f12
	d.pfadd.ss	f9,f14,f24	;	fld.d	88(SRC),f18
	d.pfsub.ss	f10,f17,f25	;	fst.d	f20,0(DEST)
	d.pfadd.ss	f11,f16,f26	;	fst.d	f22,8(DEST)
	d.pfsub.ss	f12,f19,f27	;	fst.d	f24,16(DEST)
	d.pfadd.ss	f13,f18,f28	;	fst.d	f26,24(DEST)
	pfadd.ss	f0,f0,f29	;	nop
	pfadd.ss	f0,f0,f30	;	fst.d	f28,32(DEST)
	pfadd.ss	f0,f0,f31
	fst.d	f30,40(DEST)

	ld.l	12(sp),r1
	adds	16,sp,sp
	bri	r1
	 nop
.L43:	//case TUP
	fld.d	0(SRC),f8
	fld.d	48(SRC),f14
	fld.d	8(SRC),f10
	fld.d	56(SRC),f16
.align 8
	d.pfadd.ss	f8,f14,f0	;	fld.d	16(SRC),f12
	d.pfadd.ss	f9,f15,f0	;	fld.d	64(SRC),f18
	d.pfadd.ss	f10,f16,f0	;	fld.d	24(SRC),f8
	d.pfadd.ss	f11,f17,f20	;	fld.d	72(SRC),f14
	d.pfadd.ss	f12,f18,f21	;	fld.d	32(SRC),f10
	d.pfadd.ss	f13,f19,f22	;	fld.d	80(SRC),f16
	d.pfadd.ss	f8,f14,f23	;	fld.d	40(SRC),f12
	d.pfadd.ss	f9,f15,f24	;	fld.d	88(SRC),f18
	d.pfadd.ss	f10,f16,f25	;	fst.d	f20,0(DEST)
	d.pfadd.ss	f11,f17,f26	;	fst.d	f22,8(DEST)
	d.pfadd.ss	f12,f18,f27	;	fst.d	f24,16(DEST)
	d.pfadd.ss	f13,f19,f28	;	fst.d	f26,24(DEST)
	pfadd.ss	f0,f0,f29	;	nop
	pfadd.ss	f0,f0,f30	;	fst.d	f28,32(DEST)
	pfadd.ss	f0,f0,f31
	fst.d	f30,40(DEST)

	ld.l	12(sp),r1
	adds	16,sp,sp
	bri	r1
	 nop
.L49:	//case TDOWN
	fld.d	0(SRC),f8
	fld.d	48(SRC),f14
	fld.d	8(SRC),f10
	fld.d	56(SRC),f16
.align 8
	d.pfsub.ss	f8,f14,f0	;	fld.d	16(SRC),f12
	d.pfsub.ss	f9,f15,f0	;	fld.d	64(SRC),f18
	d.pfsub.ss	f10,f16,f0	;	fld.d	24(SRC),f8
	d.pfsub.ss	f11,f17,f20	;	fld.d	72(SRC),f14
	d.pfsub.ss	f12,f18,f21	;	fld.d	32(SRC),f10
	d.pfsub.ss	f13,f19,f22	;	fld.d	80(SRC),f16
	d.pfsub.ss	f8,f14,f23	;	fld.d	40(SRC),f12
	d.pfsub.ss	f9,f15,f24	;	fld.d	88(SRC),f18
	d.pfsub.ss	f10,f16,f25	;	fst.d	f20,0(DEST)
	d.pfsub.ss	f11,f17,f26	;	fst.d	f22,8(DEST)
	d.pfsub.ss	f12,f18,f27	;	fst.d	f24,16(DEST)
	d.pfsub.ss	f13,f19,f28	;	fst.d	f26,24(DEST)
	pfadd.ss	f0,f0,f29	;	nop
	pfadd.ss	f0,f0,f30	;	fst.d	f28,32(DEST)
	pfadd.ss	f0,f0,f31
	fst.d	f30,40(DEST)

	ld.l	12(sp),r1
	adds	16,sp,sp
	bri	r1
	 nop
.L55:
	orh	h%.L164,r0,DIR
	call	_printf
	 or	l%.L164,DIR,DIR


	ld.l	12(sp),r1
	adds	16,sp,sp
	bri	r1
	 nop
	.align	4
	.data
//_i	r16	local
//.L152	.L164	static
.L56:
.L153:	.long	.L7
	.long	.L19
	.long	.L31
	.long	.L43
	.long	.L49
	.long	.L37
	.long	.L25
	.long	.L13
.L164:	.byte	66,65,68,32
	.byte	67,65,76,76
	.byte	32,84,79,32
	.byte	87,80,95,83
	.byte	72,82,73,78
	.byte	75,40,41,10
	.byte	0

//_src	r16	local
//_dest	r17	local
//_dir	r18	local
//_sign	r19	local

	.globl	_wp_shrink
