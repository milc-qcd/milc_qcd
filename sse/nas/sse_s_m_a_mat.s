;
; scalar_mult_add_su3_matrix( su3_matrix *a, su3_matrix *b, float s, su3_matrix *c)
; 
; C <- A + s*B

global scalar_mult_add_su3_matrix
scalar_mult_add_su3_matrix:
	push		ebp
	mov		ebp,esp
	push		eax
	push		ebx
	push		ecx
	mov		eax,[ebp+8]			; su3_matrix *a
	mov		ebx,[ebp+12]			; su3_matrix *b
	mov		ecx,[ebp+20]			; su3_matrix *c
	
	; get scalar value, s
	movss		xmm4,[ebp+16]			; float s		<(cc)>
	shufps		xmm4,xmm4,0x00

	; process (1,1) and (1,2)
	movups		xmm0,[eax]			;			<(aa)->e[0][0]>
	movups		xmm1,[ebx]			;			<(bb)->e[0][0]>
	mulps		xmm1,xmm4
	addps		xmm0,xmm1
	movups		[ecx],xmm0			;			<(dd)->e[0][0]>

	; process (1,3) and (2,1)
	movups		xmm0,[eax+16]			;			<(aa)->e[0][2]>
	movups		xmm1,[ebx+16]			;			<(bb)->e[0][2]>
	mulps		xmm1,xmm4
	addps		xmm0,xmm1
	movups		[ecx+16],xmm0			;			<(dd)->e[0][2]>
	
	; process (2,2) and (2,3)
	movups		xmm0,[eax+32]			;			<(aa)->e[1][1]>
	movups		xmm1,[ebx+32]			;			<(bb)->e[1][1]>
	mulps		xmm1,xmm4
	addps		xmm0,xmm1
	movups		[ecx+32],xmm0			;			<(dd)->e[1][1]>
	
	; process (3,1) and (3,2)
	movups		xmm0,[eax+48]			;			<(aa)->e[2][0]>
	movups		xmm1,[ebx+48]			;			<(bb)->e[2][0]>
	mulps		xmm1,xmm4
	addps		xmm0,xmm1
	movups		[ecx+48],xmm0			;			<(dd)->e[2][0]>
	
	; process (3,3)
	movlps		xmm0,[eax+64]			;			<(aa)->e[2][2]>
	movlps		xmm1,[ebx+64]			;			<(bb)->e[2][2]>
	mulps		xmm1,xmm4
	addps		xmm0,xmm1
	movlps		[ecx+64],xmm0			;			<(dd)->e[2][2]>
	
here:	pop	ecx
	pop	ebx
	pop	eax
	mov	esp,ebp
	pop	ebp
	ret
	
	
