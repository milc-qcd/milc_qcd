;
; su3_projector( su3_vector *a, su3_vector *b, su3_matrix *c)
;
; C = a (outer) b* 
;
; Modified by Steve Whalen, Cray, Inc.

	bits		64

global su3_projector
su3_projector:

	; First multiply all of vector a by b[0] and b[1]
	movlps		xmm0,[rsi]			;				<(bb)->c[0]>
	movhps		xmm0,[rsi+8]			;				<(bb)->c[1]>
	movaps		xmm1,xmm0
	shufps		xmm1,xmm1,0xb1
	xorps		xmm1,[negate]			;				<_sse_sgn24>

	movlps		xmm2,[rdi]			;				<(aa)->c[0]>
	movaps		xmm10,xmm2	; for later use
	unpcklps	xmm2,xmm2
	movhlps		xmm3,xmm2
	movlhps		xmm2,xmm2
	movlhps		xmm3,xmm3
	movlps		xmm4,[rdi+8]			;				<(aa)->c[1]>
	movlhps		xmm10,xmm4	; for later use
	unpcklps	xmm4,xmm4
	movhlps		xmm5,xmm4
	movlhps		xmm4,xmm4
	movlhps		xmm5,xmm5
	movlps		xmm6,[rdi+16]			;				<(aa)->c[2]>
	unpcklps	xmm6,xmm6
	movaps		xmm8,xmm6	; for later use
	movhlps		xmm7,xmm6
	movlhps		xmm6,xmm6
	movlhps		xmm7,xmm7

	mulps		xmm2,xmm0
	mulps		xmm4,xmm0
	mulps		xmm6,xmm0

	mulps		xmm3,xmm1
	mulps		xmm5,xmm1
	mulps		xmm7,xmm1

	addps		xmm2,xmm3
	xorps		xmm2,[negate]			;				<_sse_sgn24>
	movlps		[rdx],xmm2			;				<(cc)->e[0][0]>
	movhps		[rdx+8],xmm2			;				<(cc)->e[0][1]>
	addps		xmm4,xmm5
	xorps		xmm4,[negate]			;				<_sse_sgn24>
	movlps		[rdx+24],xmm4			;				<(cc)->e[1][0]>
	movhps		[rdx+32],xmm4			;				<(cc)->e[1][1]>
	addps		xmm6,xmm7
	xorps		xmm6,[negate]			;				<_sse_sgn24>
	movlps		[rdx+48],xmm6			;				<(cc)->e[2][0]>
	movhps		[rdx+56],xmm6			;				<(cc)->e[2][1]>

	; Next, multiply b[2] by a[0] and a[1]
	movaps		xmm1,xmm10
	shufps		xmm1,xmm1,0xb1
	xorps		xmm1,[negate]			;				<_sse_sgn24>

	movlps		xmm2,[rsi+16]			;				<(bb)->c[2]>
	movaps		xmm11,xmm2	; for later use
	unpcklps	xmm2,xmm2
	movhlps		xmm3,xmm2
	movlhps		xmm2,xmm2
	movlhps		xmm3,xmm3
	mulps		xmm2,xmm10
	mulps		xmm3,xmm1
	addps		xmm2,xmm3
	movlps		[rdx+16],xmm2			;				<(cc)->e[0][2]>
	movhps		[rdx+40],xmm2			;				<(cc)->e[1][2]>

	; Finally, do a[2] X b[2]
	movlhps		xmm11,xmm11		; b2.i b2.r b2.i b2.r
	mulps		xmm8,xmm11
	xorps		xmm8,[neg2]		; <_sse_sgn2>
	shufps		xmm8,xmm8,0xb4

	; need to add high bits to low bits in xmm8
	movhlps         xmm5,xmm8
	addps           xmm8,xmm5
	movlps		[rdx+64],xmm8		;				<(cc)->e[2][2]>

here:	ret
	
	align		16
negate:	dd		0x00000000
	dd		0x80000000
	dd		0x00000000
	dd		0x80000000

neg2:	dd		0x00000000
	dd		0x80000000
	dd		0x00000000
	dd		0x00000000

	
