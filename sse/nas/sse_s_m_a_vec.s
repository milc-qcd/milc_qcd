;
; scalar_mult_add_su3_vector( su3_vector *a, su3_vector *b, float s, su3_vector *c)
; 
;

global scalar_mult_add_su3_vector
scalar_mult_add_su3_vector:
	push		ebp
	mov		ebp,esp
	push		eax
	push		ebx
	push		ecx
	mov		eax,[ebp+8]			; su3_vector *a
	mov		ebx,[ebp+12]			; su3_vector *b
	mov		ecx,[ebp+20]			; su3_vector *c
	
	movss		xmm4,[ebp+16]			; float s		<(cc)>
	shufps		xmm4,xmm4,0x00
	movups		xmm0,[eax]			;			<(aa)->c[0]>
	movlps		xmm1,[eax+16]			;			<(aa)->c[2]>
	shufps		xmm1,xmm1,0x44			; make sure high floats don't cause fp errors
	movups		xmm2,[ebx]			;			<(bb)->c[0]>
	movlps		xmm3,[ebx+16]			;			<(bb)->c[2]>
	shufps		xmm3,xmm3,0x44			; make sure high floats don't cause fp errors
	mulps		xmm2,xmm4
	mulps		xmm3,xmm4
	addps		xmm0,xmm2
	addps		xmm1,xmm3
	movups		[ecx],xmm0			;			<(dd)->c[0]>
	movlps		[ecx+16],xmm1			;			<(dd)->c[2]>

here:	pop	ecx
	pop	ebx
	pop	eax
	mov	esp,ebp
	pop	ebp
	ret
	
	
