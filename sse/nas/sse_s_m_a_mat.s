;
; scalar_mult_add_su3_matrix( su3_matrix *a, su3_matrix *b, float s, su3_matrix *c)
; 
; C <- A + s*B
;
; Modified by Steve Whalen, Cray, Inc.

	bits		64

global scalar_mult_add_su3_matrix
scalar_mult_add_su3_matrix:

	; get scalar value, s
	;!inline	movss xmm0,<cc>
	unpcklps	xmm0,xmm0
	movlhps		xmm0,xmm0

	; process (1,1) and (1,2)
	movlps		xmm2,[rsi]			;			<(bb)->e[0][0]>
	movhps		xmm2,[rsi+8]			;			<(bb)->e[0][1]>
	mulps		xmm2,xmm0
	movlps		xmm1,[rdi]			;			<(aa)->e[0][0]>
	movhps		xmm1,[rdi+8]			;			<(aa)->e[0][1]>
	addps		xmm1,xmm2
	movlps		[rdx],xmm1			;			<(dd)->e[0][0]>
	movhps		[rdx+8],xmm1			;			<(dd)->e[0][1]>

	; process (1,3) and (2,1)
	movlps		xmm4,[rsi+16]			;			<(bb)->e[0][2]>
	movhps		xmm4,[rsi+24]			;			<(bb)->e[1][0]>
	mulps		xmm4,xmm0
	movlps		xmm3,[rdi+16]			;			<(aa)->e[0][2]>
	movhps		xmm3,[rdi+24]			;			<(aa)->e[1][0]>
	addps		xmm3,xmm4
	movlps		[rdx+16],xmm3			;			<(dd)->e[0][2]>
	movhps		[rdx+24],xmm3			;			<(dd)->e[1][0]>
	
	; process (2,2) and (2,3)
	movlps		xmm6,[rsi+32]			;			<(bb)->e[1][1]>
	movhps		xmm6,[rsi+40]			;			<(bb)->e[1][2]>
	mulps		xmm6,xmm0
	movlps		xmm5,[rdi+32]			;			<(aa)->e[1][1]>
	movhps		xmm5,[rdi+40]			;			<(aa)->e[1][2]>
	addps		xmm5,xmm6
	movlps		[rdx+32],xmm5			;			<(dd)->e[1][1]>
	movhps		[rdx+40],xmm5			;			<(dd)->e[1][2]>
	
	; process (3,1) and (3,2)
	movlps		xmm8,[rsi+48]			;			<(bb)->e[2][0]>
	movhps		xmm8,[rsi+56]			;			<(bb)->e[2][1]>
	mulps		xmm8,xmm0
	movlps		xmm7,[rdi+48]			;			<(aa)->e[2][0]>
	movhps		xmm7,[rdi+56]			;			<(aa)->e[2][1]>
	addps		xmm7,xmm8
	movlps		[rdx+48],xmm7			;			<(dd)->e[2][0]>
	movhps		[rdx+56],xmm7			;			<(dd)->e[2][1]>
	
	; process (3,3)
	movlps		xmm10,[rsi+64]			;			<(bb)->e[2][2]>
	mulps		xmm10,xmm0
	movlps		xmm9,[rdi+64]			;			<(aa)->e[2][2]>
	addps		xmm9,xmm10
	movlps		[rdx+64],xmm9			;			<(dd)->e[2][2]>
	
here:	ret
	
