//**************** wp_grow_a.m4 **************************
//i860 assembler language version of wp_grow_a.c
//
//void wp_grow_add( src, dest, dir, sign)
//half_wilson_vector *src;
//wilson_vector *dest;
//int dir,sign;

	.text
	.align	8
_wp_grow_add:
//	    .bf
	adds	-16,sp,sp
	st.l	r1,12(sp)
	adds	-1,r0,r30
	btne	r30,r19,.L4	//if isign== -1, dir = OPPDIR(dir)
	subs	7,r18,r18
.L4:

	subu	7,r18,r0	//check for dir out of range
	shl	2,r18,r18
	bnc	.L55
	orh	h%.L153,r0,r30
	or	l%.L153,r30,r30
	adds	r30,r18,r30
	ld.l	0(r30),r30
	bri	r30		// switch(dir)
	 nop

//	case XUP
//	for(i=0;i<3;i++){
//	    CSUM( dest->d[0].c[i], src->h[0].c[i]);
//	    CSUM( dest->d[1].c[i], src->h[1].c[i]);
//	    CSUM_TMI( dest->d[2].c[i], src->h[1].c[i] );
//	    CSUM_TMI( dest->d[3].c[i], src->h[0].c[i] );
//	}
//	break;
.L7:
				fld.d	0(r16),f8	//src spin 0 color 0
				fld.d	0(r17),f14	//dest spin 0 color 0
.align 8
	d.pfadd.ss f0,f0,f0;	fld.d	8(r16),f10	//src spin 0 color 1
	d.pfadd.ss f0,f0,f0;	fld.d	8(r17),f16	//dest spin 0 color 1
	d.pfadd.ss f8,f14,f0;	fld.d	16(r16),f12	//src spin 0 color 2
	d.pfadd.ss f9,f15,f0;	fld.d	16(r17),f18	//dest spin 0 color 2
	d.pfadd.ss f10,f16,f0;	fld.d	72(r17),f20	//dest spin 3 color 0
	d.pfadd.ss f11,f17,f26;	fld.d	80(r17),f22	//dest spin 3 color 1
	d.pfadd.ss f12,f18,f27;	fld.d	88(r17),f24	//dest spin 3 color 2
	d.pfadd.ss f13,f19,f28;	nop
	d.pfadd.ss f20,f9,f29;	fst.d	f26,0(r17)	//dest spin 0 color 0
	d.pfsub.ss f21,f8,f30;	nop
	d.pfadd.ss f22,f11,f31;	fst.d	f28,8(r17)	//dest spin 0 color 1
	d.pfsub.ss f23,f10,f26;	nop
	d.pfadd.ss f24,f13,f27;	fst.d	f30,16(r17)	//dest spin 0 color 2
	d.pfsub.ss f25,f12,f28;	nop
	d.pfadd.ss f0,f0,f29;	fst.d	f26,72(r17)	//dest spin 3 color 0
	pfadd.ss f0,f0,f30;	nop
	pfadd.ss f0,f0,f31;	fst.d	f28,80(r17)	//dest spin 3 color 1
				fst.d	f30,88(r17)	//dest spin 3 color 2

				fld.d	24(r16),f8	//src spin 1 color 0
				fld.d	24(r17),f14	//dest spin 1 color 0
.align 8
	d.pfadd.ss f0,f0,f0;	fld.d	32(r16),f10	//src spin 1 color 1
	d.pfadd.ss f0,f0,f0;	fld.d	32(r17),f16	//dest spin 1 color 1
	d.pfadd.ss f8,f14,f0;	fld.d	40(r16),f12	//src spin 1 color 2
	d.pfadd.ss f9,f15,f0;	fld.d	40(r17),f18	//dest spin 1 color 2
	d.pfadd.ss f10,f16,f0;	fld.d	48(r17),f20	//dest spin 2 color 0
	d.pfadd.ss f11,f17,f26;	fld.d	56(r17),f22	//dest spin 2 color 1
	d.pfadd.ss f12,f18,f27;	fld.d	64(r17),f24	//dest spin 2 color 2
	d.pfadd.ss f13,f19,f28;	nop
	d.pfadd.ss f20,f9,f29;	fst.d	f26,24(r17)	//dest spin 1 color 0
	d.pfsub.ss f21,f8,f30;	nop
	d.pfadd.ss f22,f11,f31;	fst.d	f28,32(r17)	//dest spin 1 color 1
	d.pfsub.ss f23,f10,f26;	nop
	d.pfadd.ss f24,f13,f27;	fst.d	f30,40(r17)	//dest spin 1 color 2
	d.pfsub.ss f25,f12,f28;	nop
	d.pfadd.ss f0,f0,f29;	fst.d	f26,48(r17)	//dest spin 2 color 0
	pfadd.ss f0,f0,f30;	nop
	pfadd.ss f0,f0,f31;	fst.d	f28,56(r17)	//dest spin 2 color 1
				fst.d	f30,64(r17)	//dest spin 2 color 2

	ld.l	12(sp),r1
	adds	16,sp,sp
	bri	r1
	 nop

//    case XDOWN:
//	for(i=0;i<3;i++){
//	    CSUM( dest->d[0].c[i], src->h[0].c[i]);
//	    CSUM( dest->d[1].c[i], src->h[1].c[i]);
//	    CSUM_TPI( dest->d[2].c[i], src->h[1].c[i] );
//	    CSUM_TPI( dest->d[3].c[i], src->h[0].c[i] );
//	}
//	break;
.L13:
				fld.d	0(r16),f8	//src spin 0 color 0
				fld.d	0(r17),f14	//dest spin 0 color 0
.align 8
	d.pfadd.ss f0,f0,f0;	fld.d	8(r16),f10	//src spin 0 color 1
	d.pfadd.ss f0,f0,f0;	fld.d	8(r17),f16	//dest spin 0 color 1
	d.pfadd.ss f8,f14,f0;	fld.d	16(r16),f12	//src spin 0 color 2
	d.pfadd.ss f9,f15,f0;	fld.d	16(r17),f18	//dest spin 0 color 2
	d.pfadd.ss f10,f16,f0;	fld.d	72(r17),f20	//dest spin 3 color 0
	d.pfadd.ss f11,f17,f26;	fld.d	80(r17),f22	//dest spin 3 color 1
	d.pfadd.ss f12,f18,f27;	fld.d	88(r17),f24	//dest spin 3 color 2
	d.pfadd.ss f13,f19,f28;	nop
	d.pfsub.ss f20,f9,f29;	fst.d	f26,0(r17)	//dest spin 0 color 0
	d.pfadd.ss f21,f8,f30;	nop
	d.pfsub.ss f22,f11,f31;	fst.d	f28,8(r17)	//dest spin 0 color 1
	d.pfadd.ss f23,f10,f26;	nop
	d.pfsub.ss f24,f13,f27;	fst.d	f30,16(r17)	//dest spin 0 color 2
	d.pfadd.ss f25,f12,f28;	nop
	d.pfadd.ss f0,f0,f29;	fst.d	f26,72(r17)	//dest spin 3 color 0
	pfadd.ss f0,f0,f30;	nop
	pfadd.ss f0,f0,f31;	fst.d	f28,80(r17)	//dest spin 3 color 1
				fst.d	f30,88(r17)	//dest spin 3 color 2

				fld.d	24(r16),f8	//src spin 1 color 0
				fld.d	24(r17),f14	//dest spin 1 color 0
.align 8
	d.pfadd.ss f0,f0,f0;	fld.d	32(r16),f10	//src spin 1 color 1
	d.pfadd.ss f0,f0,f0;	fld.d	32(r17),f16	//dest spin 1 color 1
	d.pfadd.ss f8,f14,f0;	fld.d	40(r16),f12	//src spin 1 color 2
	d.pfadd.ss f9,f15,f0;	fld.d	40(r17),f18	//dest spin 1 color 2
	d.pfadd.ss f10,f16,f0;	fld.d	48(r17),f20	//dest spin 2 color 0
	d.pfadd.ss f11,f17,f26;	fld.d	56(r17),f22	//dest spin 2 color 1
	d.pfadd.ss f12,f18,f27;	fld.d	64(r17),f24	//dest spin 2 color 2
	d.pfadd.ss f13,f19,f28;	nop
	d.pfsub.ss f20,f9,f29;	fst.d	f26,24(r17)	//dest spin 1 color 0
	d.pfadd.ss f21,f8,f30;	nop
	d.pfsub.ss f22,f11,f31;	fst.d	f28,32(r17)	//dest spin 1 color 1
	d.pfadd.ss f23,f10,f26;	nop
	d.pfsub.ss f24,f13,f27;	fst.d	f30,40(r17)	//dest spin 1 color 2
	d.pfadd.ss f25,f12,f28;	nop
	d.pfsub.ss f0,f0,f29;	fst.d	f26,48(r17)	//dest spin 2 color 0
	pfadd.ss f0,f0,f30;	nop
	pfadd.ss f0,f0,f31;	fst.d	f28,56(r17)	//dest spin 2 color 1
				fst.d	f30,64(r17)	//dest spin 2 color 2

	ld.l	12(sp),r1
	adds	16,sp,sp
	bri	r1
	 nop

//    case YUP:
//	for(i=0;i<3;i++){
//	    CSUM( dest->d[0].c[i], src->h[0].c[i]);
//	    CSUM( dest->d[1].c[i], src->h[1].c[i]);
//	    CSUM( dest->d[2].c[i], src->h[1].c[i]);
//	    CSUB( dest->d[3].c[i], src->h[0].c[i], dest->d[3].c[i] );
//	}
//	break;
.L19:
				fld.d	0(r16),f8	//src spin 0 color 0
				fld.d	0(r17),f14	//dest spin 0 color 0
.align 8
	d.pfadd.ss f0,f0,f0;	fld.d	8(r16),f10	//src spin 0 color 1
	d.pfadd.ss f0,f0,f0;	fld.d	8(r17),f16	//dest spin 0 color 1
	d.pfadd.ss f8,f14,f0;	fld.d	16(r16),f12	//src spin 0 color 2
	d.pfadd.ss f9,f15,f0;	fld.d	16(r17),f18	//dest spin 0 color 2
	d.pfadd.ss f10,f16,f0;	fld.d	72(r17),f20	//dest spin 3 color 0
	d.pfadd.ss f11,f17,f26;	fld.d	80(r17),f22	//dest spin 3 color 1
	d.pfadd.ss f12,f18,f27;	fld.d	88(r17),f24	//dest spin 3 color 2
	d.pfadd.ss f13,f19,f28;	nop
	d.pfsub.ss f20,f8,f29;	fst.d	f26,0(r17)	//dest spin 0 color 0
	d.pfsub.ss f21,f9,f30;	nop
	d.pfsub.ss f22,f10,f31;	fst.d	f28,8(r17)	//dest spin 0 color 1
	d.pfsub.ss f23,f11,f26;	nop
	d.pfsub.ss f24,f12,f27;	fst.d	f30,16(r17)	//dest spin 0 color 2
	d.pfsub.ss f25,f13,f28;	nop
	d.pfadd.ss f0,f0,f29;	fst.d	f26,72(r17)	//dest spin 3 color 0
	pfadd.ss f0,f0,f30;	nop
	pfadd.ss f0,f0,f31;	fst.d	f28,80(r17)	//dest spin 3 color 1
				fst.d	f30,88(r17)	//dest spin 3 color 2

				fld.d	24(r16),f8	//src spin 1 color 0
				fld.d	24(r17),f14	//dest spin 1 color 0
.align 8
	d.pfadd.ss f0,f0,f0;	fld.d	32(r16),f10	//src spin 1 color 1
	d.pfadd.ss f0,f0,f0;	fld.d	32(r17),f16	//dest spin 1 color 1
	d.pfadd.ss f8,f14,f0;	fld.d	40(r16),f12	//src spin 1 color 2
	d.pfadd.ss f9,f15,f0;	fld.d	40(r17),f18	//dest spin 1 color 2
	d.pfadd.ss f10,f16,f0;	fld.d	48(r17),f20	//dest spin 2 color 0
	d.pfadd.ss f11,f17,f26;	fld.d	56(r17),f22	//dest spin 2 color 1
	d.pfadd.ss f12,f18,f27;	fld.d	64(r17),f24	//dest spin 2 color 2
	d.pfadd.ss f13,f19,f28;	nop
	d.pfadd.ss f20,f8,f29;	fst.d	f26,24(r17)	//dest spin 1 color 0
	d.pfadd.ss f21,f9,f30;	nop
	d.pfadd.ss f22,f10,f31;	fst.d	f28,32(r17)	//dest spin 1 color 1
	d.pfadd.ss f23,f11,f26;	nop
	d.pfadd.ss f24,f12,f27;	fst.d	f30,40(r17)	//dest spin 1 color 2
	d.pfadd.ss f25,f13,f28;	nop
	d.pfadd.ss f0,f0,f29;	fst.d	f26,48(r17)	//dest spin 2 color 0
	pfadd.ss f0,f0,f30;	nop
	pfadd.ss f0,f0,f31;	fst.d	f28,56(r17)	//dest spin 2 color 1
				fst.d	f30,64(r17)	//dest spin 2 color 2
	ld.l	12(sp),r1
	adds	16,sp,sp
	bri	r1
	 nop

//    case YDOWN:
//	for(i=0;i<3;i++){
//	    CSUM( dest->d[0].c[i], src->h[0].c[i]);
//	    CSUM( dest->d[1].c[i], src->h[1].c[i]);
//	    CSUB( dest->d[2].c[i], src->h[1].c[i], dest->d[2].c[i] );
//	    CSUM( dest->d[3].c[i], src->h[0].c[i]);
//	}
//	break;
.L25:
				fld.d	0(r16),f8	//src spin 0 color 0
				fld.d	0(r17),f14	//dest spin 0 color 0
.align 8
	d.pfadd.ss f0,f0,f0;	fld.d	8(r16),f10	//src spin 0 color 1
	d.pfadd.ss f0,f0,f0;	fld.d	8(r17),f16	//dest spin 0 color 1
	d.pfadd.ss f8,f14,f0;	fld.d	16(r16),f12	//src spin 0 color 2
	d.pfadd.ss f9,f15,f0;	fld.d	16(r17),f18	//dest spin 0 color 2
	d.pfadd.ss f10,f16,f0;	fld.d	72(r17),f20	//dest spin 3 color 0
	d.pfadd.ss f11,f17,f26;	fld.d	80(r17),f22	//dest spin 3 color 1
	d.pfadd.ss f12,f18,f27;	fld.d	88(r17),f24	//dest spin 3 color 2
	d.pfadd.ss f13,f19,f28;	nop
	d.pfadd.ss f20,f8,f29;	fst.d	f26,0(r17)	//dest spin 0 color 0
	d.pfadd.ss f21,f9,f30;	nop
	d.pfadd.ss f22,f10,f31;	fst.d	f28,8(r17)	//dest spin 0 color 1
	d.pfadd.ss f23,f11,f26;	nop
	d.pfadd.ss f24,f12,f27;	fst.d	f30,16(r17)	//dest spin 0 color 2
	d.pfadd.ss f25,f13,f28;	nop
	d.pfadd.ss f0,f0,f29;	fst.d	f26,72(r17)	//dest spin 3 color 0
	pfadd.ss f0,f0,f30;	nop
	pfadd.ss f0,f0,f31;	fst.d	f28,80(r17)	//dest spin 3 color 1
				fst.d	f30,88(r17)	//dest spin 3 color 2

				fld.d	24(r16),f8	//src spin 1 color 0
				fld.d	24(r17),f14	//dest spin 1 color 0
.align 8
	d.pfadd.ss f0,f0,f0;	fld.d	32(r16),f10	//src spin 1 color 1
	d.pfadd.ss f0,f0,f0;	fld.d	32(r17),f16	//dest spin 1 color 1
	d.pfadd.ss f8,f14,f0;	fld.d	40(r16),f12	//src spin 1 color 2
	d.pfadd.ss f9,f15,f0;	fld.d	40(r17),f18	//dest spin 1 color 2
	d.pfadd.ss f10,f16,f0;	fld.d	48(r17),f20	//dest spin 2 color 0
	d.pfadd.ss f11,f17,f26;	fld.d	56(r17),f22	//dest spin 2 color 1
	d.pfadd.ss f12,f18,f27;	fld.d	64(r17),f24	//dest spin 2 color 2
	d.pfadd.ss f13,f19,f28;	nop
	d.pfsub.ss f20,f8,f29;	fst.d	f26,24(r17)	//dest spin 1 color 0
	d.pfsub.ss f21,f9,f30;	nop
	d.pfsub.ss f22,f10,f31;	fst.d	f28,32(r17)	//dest spin 1 color 1
	d.pfsub.ss f23,f11,f26;	nop
	d.pfsub.ss f24,f12,f27;	fst.d	f30,40(r17)	//dest spin 1 color 2
	d.pfsub.ss f25,f13,f28;	nop
	d.pfadd.ss f0,f0,f29;	fst.d	f26,48(r17)	//dest spin 2 color 0
	pfadd.ss f0,f0,f30;	nop
	pfadd.ss f0,f0,f31;	fst.d	f28,56(r17)	//dest spin 2 color 1
				fst.d	f30,64(r17)	//dest spin 2 color 2
	ld.l	12(sp),r1
	adds	16,sp,sp
	bri	r1
	 nop

//    case ZUP:
//	for(i=0;i<3;i++){
//	    CSUM( dest->d[0].c[i], src->h[0].c[i]);
//	    CSUM( dest->d[1].c[i], src->h[1].c[i]);
//	    CSUM_TMI( dest->d[2].c[i], src->h[0].c[i] );
//	    CSUM_TPI( dest->d[3].c[i], src->h[1].c[i] );
//	}
//	break;
.L31:
				fld.d	0(r16),f8	//src spin 0 color 0
				fld.d	0(r17),f14	//dest spin 0 color 0
.align 8
	d.pfadd.ss f0,f0,f0;	fld.d	8(r16),f10	//src spin 0 color 1
	d.pfadd.ss f0,f0,f0;	fld.d	8(r17),f16	//dest spin 0 color 1
	d.pfadd.ss f8,f14,f0;	fld.d	16(r16),f12	//src spin 0 color 2
	d.pfadd.ss f9,f15,f0;	fld.d	16(r17),f18	//dest spin 0 color 2
	d.pfadd.ss f10,f16,f0;	fld.d	48(r17),f20	//dest spin 2 color 0
	d.pfadd.ss f11,f17,f26;	fld.d	56(r17),f22	//dest spin 2 color 1
	d.pfadd.ss f12,f18,f27;	fld.d	64(r17),f24	//dest spin 2 color 2
	d.pfadd.ss f13,f19,f28;	nop
	d.pfadd.ss f20,f9,f29;	fst.d	f26,0(r17)	//dest spin 0 color 0
	d.pfsub.ss f21,f8,f30;	nop
	d.pfadd.ss f22,f11,f31;	fst.d	f28,8(r17)	//dest spin 0 color 1
	d.pfsub.ss f23,f10,f26;	nop
	d.pfadd.ss f24,f13,f27;	fst.d	f30,16(r17)	//dest spin 0 color 2
	d.pfsub.ss f25,f12,f28;	nop
	d.pfadd.ss f0,f0,f29;	fst.d	f26,48(r17)	//dest spin 2 color 0
	pfadd.ss f0,f0,f30;	nop
	pfadd.ss f0,f0,f31;	fst.d	f28,56(r17)	//dest spin 2 color 1
				fst.d	f30,64(r17)	//dest spin 2 color 2

				fld.d	24(r16),f8	//src spin 1 color 0
				fld.d	24(r17),f14	//dest spin 1 color 0
.align 8
	d.pfadd.ss f0,f0,f0;	fld.d	32(r16),f10	//src spin 1 color 1
	d.pfadd.ss f0,f0,f0;	fld.d	32(r17),f16	//dest spin 1 color 1
	d.pfadd.ss f8,f14,f0;	fld.d	40(r16),f12	//src spin 1 color 2
	d.pfadd.ss f9,f15,f0;	fld.d	40(r17),f18	//dest spin 1 color 2
	d.pfadd.ss f10,f16,f0;	fld.d	72(r17),f20	//dest spin 3 color 0
	d.pfadd.ss f11,f17,f26;	fld.d	80(r17),f22	//dest spin 3 color 1
	d.pfadd.ss f12,f18,f27;	fld.d	88(r17),f24	//dest spin 3 color 2
	d.pfadd.ss f13,f19,f28;	nop
	d.pfsub.ss f20,f9,f29;	fst.d	f26,24(r17)	//dest spin 1 color 0
	d.pfadd.ss f21,f8,f30;	nop
	d.pfsub.ss f22,f11,f31;	fst.d	f28,32(r17)	//dest spin 1 color 1
	d.pfadd.ss f23,f10,f26;	nop
	d.pfsub.ss f24,f13,f27;	fst.d	f30,40(r17)	//dest spin 1 color 2
	d.pfadd.ss f25,f12,f28;	nop
	d.pfadd.ss f0,f0,f29;	fst.d	f26,72(r17)	//dest spin 3 color 0
	pfadd.ss f0,f0,f30;	nop
	pfadd.ss f0,f0,f31;	fst.d	f28,80(r17)	//dest spin 3 color 1
				fst.d	f30,88(r17)	//dest spin 3 color 2
	ld.l	12(sp),r1
	adds	16,sp,sp
	bri	r1
	 nop

//    case ZDOWN:
//	for(i=0;i<3;i++){
//	    CSUM( dest->d[0].c[i], src->h[0].c[i]);
//	    CSUM( dest->d[1].c[i], src->h[1].c[i]);
//	    CSUM_TPI( dest->d[2].c[i], src->h[0].c[i] );
//	    CSUM_TMI( dest->d[3].c[i], src->h[1].c[i] );
//	}
//	break;
.L37:
				fld.d	0(r16),f8	//src spin 0 color 0
				fld.d	0(r17),f14	//dest spin 0 color 0
.align 8
	d.pfadd.ss f0,f0,f0;	fld.d	8(r16),f10	//src spin 0 color 1
	d.pfadd.ss f0,f0,f0;	fld.d	8(r17),f16	//dest spin 0 color 1
	d.pfadd.ss f8,f14,f0;	fld.d	16(r16),f12	//src spin 0 color 2
	d.pfadd.ss f9,f15,f0;	fld.d	16(r17),f18	//dest spin 0 color 2
	d.pfadd.ss f10,f16,f0;	fld.d	48(r17),f20	//dest spin 2 color 0
	d.pfadd.ss f11,f17,f26;	fld.d	56(r17),f22	//dest spin 2 color 1
	d.pfadd.ss f12,f18,f27;	fld.d	64(r17),f24	//dest spin 2 color 2
	d.pfadd.ss f13,f19,f28;	nop
	d.pfsub.ss f20,f9,f29;	fst.d	f26,0(r17)	//dest spin 0 color 0
	d.pfadd.ss f21,f8,f30;	nop
	d.pfsub.ss f22,f11,f31;	fst.d	f28,8(r17)	//dest spin 0 color 1
	d.pfadd.ss f23,f10,f26;	nop
	d.pfsub.ss f24,f13,f27;	fst.d	f30,16(r17)	//dest spin 0 color 2
	d.pfadd.ss f25,f12,f28;	nop
	d.pfadd.ss f0,f0,f29;	fst.d	f26,48(r17)	//dest spin 2 color 0
	pfadd.ss f0,f0,f30;	nop
	pfadd.ss f0,f0,f31;	fst.d	f28,56(r17)	//dest spin 2 color 1
				fst.d	f30,64(r17)	//dest spin 2 color 2

				fld.d	24(r16),f8	//src spin 1 color 0
				fld.d	24(r17),f14	//dest spin 1 color 0
.align 8
	d.pfadd.ss f0,f0,f0;	fld.d	32(r16),f10	//src spin 1 color 1
	d.pfadd.ss f0,f0,f0;	fld.d	32(r17),f16	//dest spin 1 color 1
	d.pfadd.ss f8,f14,f0;	fld.d	40(r16),f12	//src spin 1 color 2
	d.pfadd.ss f9,f15,f0;	fld.d	40(r17),f18	//dest spin 1 color 2
	d.pfadd.ss f10,f16,f0;	fld.d	72(r17),f20	//dest spin 3 color 0
	d.pfadd.ss f11,f17,f26;	fld.d	80(r17),f22	//dest spin 3 color 1
	d.pfadd.ss f12,f18,f27;	fld.d	88(r17),f24	//dest spin 3 color 2
	d.pfadd.ss f13,f19,f28;	nop
	d.pfadd.ss f20,f9,f29;	fst.d	f26,24(r17)	//dest spin 1 color 0
	d.pfsub.ss f21,f8,f30;	nop
	d.pfadd.ss f22,f11,f31;	fst.d	f28,32(r17)	//dest spin 1 color 1
	d.pfsub.ss f23,f10,f26;	nop
	d.pfadd.ss f24,f13,f27;	fst.d	f30,40(r17)	//dest spin 1 color 2
	d.pfsub.ss f25,f12,f28;	nop
	d.pfadd.ss f0,f0,f29;	fst.d	f26,72(r17)	//dest spin 3 color 0
	pfadd.ss f0,f0,f30;	nop
	pfadd.ss f0,f0,f31;	fst.d	f28,80(r17)	//dest spin 3 color 1
				fst.d	f30,88(r17)	//dest spin 3 color 2
	ld.l	12(sp),r1
	adds	16,sp,sp
	bri	r1
	 nop

//	line	148
//    case TUP:
//	for(i=0;i<3;i++){
//	    CSUM( dest->d[0].c[i], src->h[0].c[i]);
//	    CSUM( dest->d[1].c[i], src->h[1].c[i]);
//	    CSUM( dest->d[2].c[i], src->h[0].c[i]);
//	    CSUM( dest->d[3].c[i], src->h[1].c[i]);
//	}
.L43:
//	break;
				fld.d	0(r16),f8	//src spin 0 color 0
				fld.d	0(r17),f14	//dest spin 0 color 0
.align 8
	d.pfadd.ss f0,f0,f0;	fld.d	8(r16),f10	//src spin 0 color 1
	d.pfadd.ss f0,f0,f0;	fld.d	8(r17),f16	//dest spin 0 color 1
	d.pfadd.ss f8,f14,f0;	fld.d	16(r16),f12	//src spin 0 color 2
	d.pfadd.ss f9,f15,f0;	fld.d	16(r17),f18	//dest spin 0 color 2
	d.pfadd.ss f10,f16,f0;	fld.d	48(r17),f20	//dest spin 2 color 0
	d.pfadd.ss f11,f17,f26;	fld.d	56(r17),f22	//dest spin 2 color 1
	d.pfadd.ss f12,f18,f27;	fld.d	64(r17),f24	//dest spin 2 color 2
	d.pfadd.ss f13,f19,f28;	nop
	d.pfadd.ss f20,f8,f29;	fst.d	f26,0(r17)	//dest spin 0 color 0
	d.pfadd.ss f21,f9,f30;	nop
	d.pfadd.ss f22,f10,f31;	fst.d	f28,8(r17)	//dest spin 0 color 1
	d.pfadd.ss f23,f11,f26;	nop
	d.pfadd.ss f24,f12,f27;	fst.d	f30,16(r17)	//dest spin 0 color 2
	d.pfadd.ss f25,f13,f28;	nop
	d.pfadd.ss f0,f0,f29;	fst.d	f26,48(r17)	//dest spin 2 color 0
	pfadd.ss f0,f0,f30;	nop
	pfadd.ss f0,f0,f31;	fst.d	f28,56(r17)	//dest spin 2 color 1
				fst.d	f30,64(r17)	//dest spin 2 color 2

				fld.d	24(r16),f8	//src spin 1 color 0
				fld.d	24(r17),f14	//dest spin 1 color 0
.align 8
	d.pfadd.ss f0,f0,f0;	fld.d	32(r16),f10	//src spin 1 color 1
	d.pfadd.ss f0,f0,f0;	fld.d	32(r17),f16	//dest spin 1 color 1
	d.pfadd.ss f8,f14,f0;	fld.d	40(r16),f12	//src spin 1 color 2
	d.pfadd.ss f9,f15,f0;	fld.d	40(r17),f18	//dest spin 1 color 2
	d.pfadd.ss f10,f16,f0;	fld.d	72(r17),f20	//dest spin 3 color 0
	d.pfadd.ss f11,f17,f26;	fld.d	80(r17),f22	//dest spin 3 color 1
	d.pfadd.ss f12,f18,f27;	fld.d	88(r17),f24	//dest spin 3 color 2
	d.pfadd.ss f13,f19,f28;	nop
	d.pfadd.ss f20,f8,f29;	fst.d	f26,24(r17)	//dest spin 1 color 0
	d.pfadd.ss f21,f9,f30;	nop
	d.pfadd.ss f22,f10,f31;	fst.d	f28,32(r17)	//dest spin 1 color 1
	d.pfadd.ss f23,f11,f26;	nop
	d.pfadd.ss f24,f12,f27;	fst.d	f30,40(r17)	//dest spin 1 color 2
	d.pfadd.ss f25,f13,f28;	nop
	d.pfadd.ss f0,f0,f29;	fst.d	f26,72(r17)	//dest spin 3 color 0
	pfadd.ss f0,f0,f30;	nop
	pfadd.ss f0,f0,f31;	fst.d	f28,80(r17)	//dest spin 3 color 1
				fst.d	f30,88(r17)	//dest spin 3 color 2
	ld.l	12(sp),r1
	adds	16,sp,sp
	bri	r1
	 nop

//	line	156
//    case TDOWN:
//	for(i=0;i<3;i++){
//	    CSUM( dest->d[0].c[i], src->h[0].c[i]);
//	    CSUM( dest->d[1].c[i], src->h[1].c[i]);
//	    CSUB( dest->d[2].c[i], src->h[0].c[i], dest->d[2].c[i] );
//	    CSUB( dest->d[3].c[i], src->h[1].c[i], dest->d[3].c[i] );
//	}
.L49:
//	break;
				fld.d	0(r16),f8	//src spin 0 color 0
				fld.d	0(r17),f14	//dest spin 0 color 0
.align 8
	d.pfadd.ss f0,f0,f0;	fld.d	8(r16),f10	//src spin 0 color 1
	d.pfadd.ss f0,f0,f0;	fld.d	8(r17),f16	//dest spin 0 color 1
	d.pfadd.ss f8,f14,f0;	fld.d	16(r16),f12	//src spin 0 color 2
	d.pfadd.ss f9,f15,f0;	fld.d	16(r17),f18	//dest spin 0 color 2
	d.pfadd.ss f10,f16,f0;	fld.d	48(r17),f20	//dest spin 2 color 0
	d.pfadd.ss f11,f17,f26;	fld.d	56(r17),f22	//dest spin 2 color 1
	d.pfadd.ss f12,f18,f27;	fld.d	64(r17),f24	//dest spin 2 color 2
	d.pfadd.ss f13,f19,f28;	nop
	d.pfsub.ss f20,f8,f29;	fst.d	f26,0(r17)	//dest spin 0 color 0
	d.pfsub.ss f21,f9,f30;	nop
	d.pfsub.ss f22,f10,f31;	fst.d	f28,8(r17)	//dest spin 0 color 1
	d.pfsub.ss f23,f11,f26;	nop
	d.pfsub.ss f24,f12,f27;	fst.d	f30,16(r17)	//dest spin 0 color 2
	d.pfsub.ss f25,f13,f28;	nop
	d.pfadd.ss f0,f0,f29;	fst.d	f26,48(r17)	//dest spin 2 color 0
	pfadd.ss f0,f0,f30;	nop
	pfadd.ss f0,f0,f31;	fst.d	f28,56(r17)	//dest spin 2 color 1
				fst.d	f30,64(r17)	//dest spin 2 color 2

				fld.d	24(r16),f8	//src spin 1 color 0
				fld.d	24(r17),f14	//dest spin 1 color 0
.align 8
	d.pfadd.ss f0,f0,f0;	fld.d	32(r16),f10	//src spin 1 color 1
	d.pfadd.ss f0,f0,f0;	fld.d	32(r17),f16	//dest spin 1 color 1
	d.pfadd.ss f8,f14,f0;	fld.d	40(r16),f12	//src spin 1 color 2
	d.pfadd.ss f9,f15,f0;	fld.d	40(r17),f18	//dest spin 1 color 2
	d.pfadd.ss f10,f16,f0;	fld.d	72(r17),f20	//dest spin 3 color 0
	d.pfadd.ss f11,f17,f26;	fld.d	80(r17),f22	//dest spin 3 color 1
	d.pfadd.ss f12,f18,f27;	fld.d	88(r17),f24	//dest spin 3 color 2
	d.pfadd.ss f13,f19,f28;	nop
	d.pfsub.ss f20,f8,f29;	fst.d	f26,24(r17)	//dest spin 1 color 0
	d.pfsub.ss f21,f9,f30;	nop
	d.pfsub.ss f22,f10,f31;	fst.d	f28,32(r17)	//dest spin 1 color 1
	d.pfsub.ss f23,f11,f26;	nop
	d.pfsub.ss f24,f12,f27;	fst.d	f30,40(r17)	//dest spin 1 color 2
	d.pfsub.ss f25,f13,f28;	nop
	d.pfadd.ss f0,f0,f29;	fst.d	f26,72(r17)	//dest spin 3 color 0
	pfadd.ss f0,f0,f30;	nop
	pfadd.ss f0,f0,f31;	fst.d	f28,80(r17)	//dest spin 3 color 1
				fst.d	f30,88(r17)	//dest spin 3 color 2
	ld.l	12(sp),r1
	adds	16,sp,sp
	bri	r1
	 nop

.L55:
//	line	164
	orh	h%.L164,r0,r16
	call	_printf
	 or	l%.L164,r16,r16
//	    .ef
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
	.byte	87,80,95,71
	.byte	82,79,87,40
	.byte	41,10
	.byte	0

//_src	r16	local
//_dest	r17	local
//_dir	r18	local
//_sign	r19	local

	.text
	.data
	.globl	_wp_grow_add

	.text
