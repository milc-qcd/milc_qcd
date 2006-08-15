;
; scalar_mult_add_su3_vector( su3_vector *a, su3_vector *b, float s, su3_vector *c)
;
; c <- a + s*b
;
; Modified by Steve Whalen, Cray, Inc.

	bits		64

global scalar_mult_add_su3_vector
scalar_mult_add_su3_vector:

	;!inline	movss xmm0,<cc>

	cvtss2sd	xmm6,xmm0			; float s		<(cc)>
	unpcklpd	xmm6,xmm6
	cvtps2pd	xmm3,[rsi]			;			<(bb)->c[0]>
	mulpd		xmm3,xmm6
	cvtps2pd	xmm0,[rdi]			;			<(aa)->c[0]>
	addpd		xmm0,xmm3
	cvtps2pd	xmm4,[rsi+8]			;			<(bb)->c[1]>
	mulpd		xmm4,xmm6
	cvtps2pd	xmm1,[rdi+8]			;			<(aa)->c[1]>
	addpd		xmm1,xmm4
	cvtps2pd	xmm5,[rsi+16]			;			<(bb)->c[2]>
	mulpd		xmm5,xmm6
	cvtps2pd	xmm2,[rdi+16]			;			<(aa)->c[2]>
	addpd		xmm2,xmm5
	cvtpd2ps	xmm3,xmm0
	cvtpd2ps	xmm4,xmm1
	cvtpd2ps	xmm5,xmm2
	movlps		[rdx],xmm3			;			<(dd)->c[0]>
	movlps		[rdx+8],xmm4			;			<(dd)->c[1]>
	movlps		[rdx+16],xmm5			;			<(dd)->c[2]>

here:	ret
	
