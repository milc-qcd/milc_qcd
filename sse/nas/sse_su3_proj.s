;
; su3_projector( su3_vector *a, su3_vector *b, su3_matrix *c)
;
; C = a (outer) b* 

global su3_projector
su3_projector:
	push		ebp
	mov		ebp,esp
	push		eax
	push		ebx
	push		ecx
	mov		eax,[ebp+8]			; su3_vector *a
	mov		ebx,[ebp+12]			; su3_vector *b
	mov		ecx,[ebp+16]			; su3_matrix *c

	; First multiply all of vector a by b[0] and b[1]
	movlps		xmm0,[ebx]			;				<(bb)->c[0]>
	movhps		xmm0,[ebx+8]			;				<(bb)->c[1]>
	movaps		xmm1,xmm0
	shufps		xmm1,xmm1,0xb1
	xorps		xmm1,[negate]			;				<_sse_sgn24>

	movss		xmm2,[eax]			;				<(aa)->c[0].real>
	shufps		xmm2,xmm2,0x00
	movss		xmm3,[eax+4]			;				<(aa)->c[0].imag>
	shufps		xmm3,xmm3,0x00
	mulps		xmm2,xmm0
	mulps		xmm3,xmm1
	addps		xmm2,xmm3
	xorps		xmm2,[negate]			;				<_sse_sgn24>
	movups		[ecx],xmm2			;				<(cc)->e[0][0]>

	movss		xmm2,[eax+8]			;				<(aa)->c[1].real>
	shufps		xmm2,xmm2,0x00
	movss		xmm3,[eax+12]			;				<(aa)->c[1].imag>
	shufps		xmm3,xmm3,0x00
	mulps		xmm2,xmm0
	mulps		xmm3,xmm1
	addps		xmm2,xmm3
	xorps		xmm2,[negate]			;				<_sse_sgn24>
	movups		[ecx+24],xmm2			;				<(cc)->e[1][0]>
	
	movss		xmm2,[eax+16]			;				<(aa)->c[2].real>
	shufps		xmm2,xmm2,0x00
	movss		xmm3,[eax+20]			;				<(aa)->c[2].imag>
	shufps		xmm3,xmm3,0x00
	mulps		xmm2,xmm0
	mulps		xmm3,xmm1
	addps		xmm2,xmm3
	xorps		xmm2,[negate]			;				<_sse_sgn24>
	movups		[ecx+48],xmm2			;				<(cc)->e[2][0]>

	; Next, multiply b[2] by a[0] and a[1]
	movlps		xmm0,[eax]			;				<(aa)->c[0]>
	movhps		xmm0,[eax+8]			;				<(aa)->c[1]>
	movaps		xmm1,xmm0
	shufps		xmm1,xmm1,0xb1
	xorps		xmm1,[negate]			;				<_sse_sgn24>

	movss		xmm2,[ebx+16]			;				<(bb)->c[2].real>
	shufps		xmm2,xmm2,0x00
	movss		xmm3,[ebx+20]			;				<(bb)->c[2].imag>
	shufps		xmm3,xmm3,0x00
	mulps		xmm2,xmm0
	mulps		xmm3,xmm1
	addps		xmm2,xmm3
	movlps		[ecx+16],xmm2			;				<(cc)->e[0][2]>
	movhps		[ecx+40],xmm2			;				<(cc)->e[1][2]>

	; Finally, do a[2] X b[2]
	movlps		xmm0,[eax+16]			;				<(aa)->c[2]>
	shufps		xmm0,xmm0,0x14
	movlps		xmm2,[ebx+16]			;				<(bb)->c[2]>
	shufps		xmm2,xmm2,0x44
	xorps		xmm2,[neg2]			;				<_sse_sgn4>
	mulps		xmm2,xmm0
	movaps		xmm1,xmm2
	shufps		xmm1,xmm1,0xd4
	shufps		xmm2,xmm2,0x8c
	addps		xmm2,xmm1
	movhps		[ecx+64],xmm2			;				<(cc)->e[2][2]>

here:	pop	ecx
	pop	ebx
	pop	eax
	mov	esp,ebp
	pop	ebp
	ret
	
	align		16
negate:	dd		0x00000000
	dd		0x80000000
	dd		0x00000000
	dd		0x80000000

neg2:	dd		0x00000000
	dd		0x00000000
	dd		0x00000000
	dd		0x80000000

	
