;
; add_su3_vector( su3_vector *a, su3_vector *b, su3_vector *c)
; 
; c <- a + b
;
; Modified by Steve Whalen, Cray, Inc.

	bits		64

global add_su3_vector
add_su3_vector:

	movlps		xmm0,[rdi]			;			<(aa)->c[0]>
	movhps		xmm0,[rdi+8]			;			<(aa)->c[1]>
	movlps		xmm2,[rsi]			;			<(bb)->c[0]>
	movhps		xmm2,[rsi+8]			;			<(bb)->c[1]>
	addps		xmm0,xmm2
	movlps		[rdx],xmm0			;			<(cc)->c[0]>
	movhps		[rdx+8],xmm0			;			<(cc)->c[1]>

	movlps		xmm3,[rsi+16]			;			<(bb)->c[2]>
	movlhps		xmm3,xmm3
	movlps		xmm1,[rdi+16]			;			<(aa)->c[2]>
	movlhps		xmm1,xmm1
	addps		xmm1,xmm3
	
	movlps		[rdx+16],xmm1			;			<(cc)->c[2]>

here:	ret
	
