;
; sub_four_su3_vecs( su3_vector *a, su3_vector *b0, su3_vector *b1,
;                                       su3_vector *b2, su3_vector *b3 )
;

global sub_four_su3_vecs
sub_four_su3_vecs:
	push		ebp
	mov		ebp,esp
	push		eax
	push		ebx
	push		ecx
	mov		eax,[ebp+8]			; su3_vector *a
	mov		ebx,[ebp+12]			; su3_vector *b

	movups		xmm0,[eax]			;			<(aa)->c[0]>
	movlps		xmm1,[eax+16]			;			<(aa)->c[2]>
	shufps		xmm1,xmm1,0x44
	movups		xmm2,[ebx]			;			<(bb0)->c[0]>
	movlps		xmm3,[ebx+16]			;			<(bb0)->c[2]>
	shufps		xmm3,xmm3,0x44
	subps		xmm0,xmm2
	subps		xmm1,xmm3
	
	mov		ebx,[ebp+16]			; su3_vector *b
	movups		xmm2,[ebx]			;			<(bb1)->c[0]>
	movlps		xmm3,[ebx+16]			;			<(bb1)->c[2]>
	shufps		xmm3,xmm3,0x44
	subps		xmm0,xmm2
	subps		xmm1,xmm3

	mov		ebx,[ebp+20]			; su3_vector *b
	movups		xmm2,[ebx]			;			<(bb2)->c[0]>
	movlps		xmm3,[ebx+16]			;			<(bb2)->c[2]>
	shufps		xmm3,xmm3,0x44
	subps		xmm0,xmm2
	subps		xmm1,xmm3

	mov		ebx,[ebp+24]			; su3_vector *b
	movups		xmm2,[ebx]			;			<(bb3)->c[0]>
	movlps		xmm3,[ebx+16]			;			<(bb3)->c[2]>
	shufps		xmm3,xmm3,0x44
	subps		xmm0,xmm2
	subps		xmm1,xmm3

	movups		[eax],xmm0			;			<(aa)->c[0]>
	movlps		[eax+16],xmm1			;			<(aa)->c[2]>

here:	pop	ecx
	pop	ebx
	pop	eax
	mov	esp,ebp
	pop	ebp
	ret
	
	
