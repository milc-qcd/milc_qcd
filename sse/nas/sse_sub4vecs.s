;
; sub_four_su3_vecs( su3_vector *a, su3_vector *b0, su3_vector *b1,
;                    su3_vector *b2, su3_vector *b3 )
;
; a <- a - b0 - b1 - b2 - b3
;
; Modified by Steve Whalen, Cray, Inc.

	bits		64

global sub_four_su3_vecs
sub_four_su3_vecs:

	movlps		xmm0,[rdi]			;			<(aa)->c[0]>
	movhps		xmm0,[rdi+8]			;			<(aa)->c[1]>
	movlps		xmm2,[rsi]			;			<(bb0)->c[0]>
	movhps		xmm2,[rsi+8]			;			<(bb0)->c[1]>
	subps		xmm0,xmm2
	movlps		xmm1,[rdi+16]			;			<(aa)->c[2]>
	movlhps		xmm1,xmm1
	movlps		xmm3,[rsi+16]			;			<(bb0)->c[2]>
	movlhps		xmm3,xmm3
	subps		xmm1,xmm3
	
	movlps		xmm4,[rdx]			;			<(bb1)->c[0]>
	movhps		xmm4,[rdx+8]			;			<(bb1)->c[1]>
	subps		xmm0,xmm4
	movlps		xmm5,[rdx+16]			;			<(bb1)->c[2]>
	movlhps		xmm5,xmm5
	subps		xmm1,xmm5

	movlps		xmm6,[rcx]			;			<(bb2)->c[0]>
	movhps		xmm6,[rcx+8]			;			<(bb2)->c[1]>
	subps		xmm0,xmm6
	movlps		xmm7,[rcx+16]			;			<(bb2)->c[2]>
	movlhps		xmm7,xmm7
	subps		xmm1,xmm7

	movlps		xmm8,[r8]			;			<(bb3)->c[0]>
	movhps		xmm8,[r8+8]			;			<(bb3)->c[1]>
	subps		xmm0,xmm8
	movlps		xmm9,[r8+16]			;			<(bb3)->c[2]>
	movlhps		xmm9,xmm9
	subps		xmm1,xmm9

	movlps		[rdi],xmm0			;			<(aa)->c[0]>
	movhps		[rdi+8],xmm0			;			<(aa)->c[1]>
	movlps		[rdi+16],xmm1			;			<(aa)->c[2]>

here:	ret
