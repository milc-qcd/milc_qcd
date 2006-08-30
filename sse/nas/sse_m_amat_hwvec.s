;
; mult_adj_su3_mat_hwvec( su3_matrix *a, half_wilson_vector *b, *c)
;

global mult_adj_su3_mat_hwvec
mult_adj_su3_mat_hwvec:
	push		ebp
	mov		ebp,esp
	push		eax
	push		ebx
	push		ecx
	mov		eax,[ebp+8]			; su3_matrix *a
	mov		ebx,[ebp+12]			; su3_vector *b[2]
	mov		ecx,[ebp+16]			; su3_vector *c[2]

	; bring in half wilson vector
	movlps		xmm0,[ebx]			;				<(bb)->h[0].c[0]>
	movlps		xmm1,[ebx+8]			;				<(bb)->h[0].c[1]>
	movlps		xmm2,[ebx+16]			;				<(bb)->h[0].c[2]>
	movhps		xmm0,[ebx+24]			; b[1]1 b[0]1			<(bb)->h[1].c[0]>
	movhps		xmm1,[ebx+32]			; b[1]2 b[0]2			<(bb)->h[1].c[1]>
	movhps		xmm2,[ebx+40]			; b[1]3 b[0]3			<(bb)->h[1].c[2]>

	; multiply vectors by real components of a
	movss		xmm3,[eax]			; a(1,1)			<(aa)->e[0][0].real>
	movss		xmm6,[eax+24]			; a(1,2)			<(aa)->e[1][0].real>
	movss		xmm4,[eax+8]			; a(2,1)			<(aa)->e[0][1].real>
	movss		xmm7,[eax+56]			; a(2,3)			<(aa)->e[2][1].real>
	movss		xmm5,[eax+16]			; a(3,1)			<(aa)->e[0][2].real>
	shufps		xmm3,xmm3,0x00
	shufps		xmm6,xmm6,0x00
	shufps		xmm4,xmm4,0x00
	mulps		xmm3,xmm0
	shufps		xmm7,xmm7,0x00
	mulps		xmm6,xmm1
	shufps		xmm5,xmm5,0x00
	mulps		xmm4,xmm0
	addps		xmm3,xmm6
	mulps		xmm7,xmm2
	mulps		xmm5,xmm0
	addps		xmm4,xmm7
	movss		xmm6,[eax+40]			; a(3,2)			<(aa)->e[1][2].real>
	movss		xmm7,[eax+48]			; a(1,3)			<(aa)->e[2][0].real>
	shufps		xmm6,xmm6,0x00
	shufps		xmm7,xmm7,0x00
	mulps		xmm6,xmm1
	mulps		xmm7,xmm2
	addps		xmm5,xmm6
	addps		xmm3,xmm7
	movss		xmm6,[eax+32]			; a(2,2)			<(aa)->e[1][1].real>
	movss		xmm7,[eax+64]			; a(3,3)			<(aa)->e[2][2].real>
	shufps		xmm6,xmm6,0x00
	shufps		xmm7,xmm7,0x00
	mulps		xmm6,xmm1
	mulps		xmm7,xmm2
	addps		xmm4,xmm6
	addps		xmm5,xmm7
	
	; multiply vectors by imaginary components of a
	movss		xmm6,[eax+4]			; a(1,1)			<(aa)->e[0][0].imag>
	movss		xmm7,[eax+36]			; a(2,2)			<(aa)->e[1][1].imag>
	shufps		xmm0,xmm0,0xb1
	shufps		xmm1,xmm1,0xb1
	shufps		xmm2,xmm2,0xb1
	shufps		xmm6,xmm6,0x00
	shufps		xmm7,xmm7,0x00
	xorps		xmm0,[negate]			;				<_sse_sgn24>
	xorps		xmm1,[negate]			; 				<_sse_sgn24>
	xorps		xmm2,[negate]			;				<_sse_sgn24> 
	mulps		xmm6,xmm0
	mulps		xmm7,xmm1
	addps		xmm3,xmm6
	addps		xmm4,xmm7
	movss		xmm6,[eax+68]			; a(3,3)			<(aa)->e[2][2].imag>
	movss		xmm7,[eax+12]			; a(2,1)			<(aa)->e[0][1].imag>
	shufps		xmm6,xmm6,0x00
	shufps		xmm7,xmm7,0x00
	mulps		xmm6,xmm2
	mulps		xmm7,xmm0
	addps		xmm5,xmm6
	addps		xmm4,xmm7
	movss		xmm6,[eax+28]			; a(1,2)			<(aa)->e[1][0].imag>
	movss		xmm7,[eax+20]			; a(3,1)			<(aa)->e[0][2].imag>
	shufps		xmm6,xmm6,0x00
	shufps		xmm7,xmm7,0x00
	mulps		xmm6,xmm1
	mulps		xmm7,xmm0
	addps		xmm3,xmm6
	addps		xmm5,xmm7
	movss		xmm0,[eax+52]			; a(1,3)			<(aa)->e[2][0].imag>
	movss		xmm6,[eax+44]			; a(3,2)			<(aa)->e[1][2].imag>
	movss		xmm7,[eax+60]			; a(2,3)			<(aa)->e[2][1].imag>
	shufps		xmm0,xmm0,0x00
	shufps		xmm6,xmm6,0x00
	shufps		xmm7,xmm7,0x00
	mulps		xmm0,xmm2
	mulps		xmm6,xmm1
	mulps		xmm7,xmm2
	addps		xmm3,xmm0
	addps		xmm5,xmm6
	addps		xmm4,xmm7

	; move results into output
	movlps		[ecx],xmm3			; 				<(cc)->h[0].c[0]>
	movlps		[ecx+8],xmm4			;				<(cc)->h[0].c[1]> 
	movlps		[ecx+16],xmm5			; 				<(cc)->h[0].c[2]>
	movhps		[ecx+24],xmm3			; 				<(cc)->h[1].c[0]>
	movhps		[ecx+32],xmm4			;				<(cc)->h[1].c[1]> 
	movhps		[ecx+40],xmm5			;				<(cc)->h[1].c[2]> 
		
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
	
