;
; mult_su3_mat_vec( su3_matrix *a, su3_vector *b, su3_vector *c)
;

global mult_su3_mat_vec
mult_su3_mat_vec:
	push		ebp
	mov		ebp,esp
	push		eax
	push		ebx
	push		ecx
	mov		eax,[ebp+8]			; su3_matrix *a
	mov		ebx,[ebp+12]			; su3_vector *b
	mov		ecx,[ebp+16]			; su3_vector *c

	;  bring in real and imaginary b vector
	movlps		xmm0,[ebx]			; x,x,b0i,b0r		<(bb)->c[0]>
	movlps		xmm1,[ebx+8]			; x,x,b1i,b1r		<(bb)->c[1]>
	movlps		xmm2,[ebx+16]			; x,x,b2i,b2r		<(bb)->c[2]>
	shufps		xmm0,xmm0,0x44			; b0i,b0r,b0i,b0r
	shufps		xmm1,xmm1,0x44			; b1i,b1r,b1i,b1r
	shufps		xmm2,xmm2,0x44			; b2i,b2r,b2i,b2r

	; bring in real components of first two rows of matrix a
	movss		xmm3,[eax]			; x,x,x,c00r		<(aa)->e[0][0].real>
	movss		xmm7,[eax+24]			; x,x,x,c10r		<(aa)->e[1][0].real>
	shufps		xmm3,xmm7,0x00			; c10r,c10r,c00r,c00r
	movss		xmm4,[eax+8]			; x,x,x,c01r		<(aa)->e[0][1].real>
	movss		xmm7,[eax+32]			; x,x,x,c11r		<(aa)->e[1][1].real>
	shufps		xmm4,xmm7,0x00			; c11r,c11r,c01r,c01r
	mulps		xmm3,xmm0
	mulps		xmm4,xmm1
	addps		xmm3,xmm4
	movss		xmm5,[eax+16]			; x,x,x,c02r		<(aa)->e[0][2].real>
	movss		xmm7,[eax+40]			; x,x,x,c12r		<(aa)->e[1][2].real>
	shufps		xmm5,xmm7,0x00			; c12r,c12r,c02r,c02r
	mulps		xmm5,xmm2
	addps		xmm3,xmm5

	; special handling of the 3rd row of matrix a
	shufps		xmm1,xmm0,0x44			; b0i,b0r,b1i,b1r
	movss		xmm7,[eax+48]			; x,x,x,c20r		<(aa)->e[2][0].real>
	movss		xmm6,[eax+56]			; x,x,x,c21r		<(aa)->e[2][1].real>
	shufps		xmm6,xmm7,0x00			; c20r,c20r,c21r,c21r
	mulps		xmm6,xmm1

	; shuffle b vector for imaginary components of matrix a
	shufps		xmm0,xmm0,0xB1			; b0r,b0i,b0r,b0i
	xorps		xmm0,[negate]			; b0r,-b0i,b0r,-b0i	<_sse_sgn13>
	shufps		xmm1,xmm1,0x11			; b1r,b1i,b1r,b1i
	xorps		xmm1,[negate]			; b1r,-b1i,b1r,-b1i	<_sse_sgn13>
	shufps		xmm2,xmm2,0xB1			; b2r,b2i,b2r,b2i
	xorps		xmm2,[negate]			; b2r,-b2i,b2r,-b2i	<_sse_sgn13>

	; bring in imaginary components of first two rows of matrix b
	movss		xmm4,[eax+4]			; x,x,x,c00i		<(aa)->e[0][0].imag>
	movss		xmm7,[eax+28]			; x,x,x,c10i		<(aa)->e[1][0].imag>
	shufps		xmm4,xmm7,0x00			; c10i,c10i,c00i,c00i
	mulps		xmm4,xmm0
	addps		xmm3,xmm4
	movss		xmm5,[eax+12]			; x,x,x,c01i		<(aa)->e[0][1].imag>
	movss		xmm7,[eax+36]			; x,x,x,c11i		<(aa)->e[1][1].imag>
	shufps		xmm5,xmm7,0x00			; c11i,c11i,c01i,c01i
	mulps		xmm5,xmm1
	addps		xmm3,xmm5
	movss		xmm5,[eax+20]			; x,x,x,c02i		<(aa)->e[0][2].imag>
	movss		xmm7,[eax+44]			; x,x,x,c12i		<(aa)->e[1][2].imag>
	shufps		xmm5,xmm7,0x00			; c12i,c12i,c02i,c02i
	mulps		xmm5,xmm2
	addps		xmm3,xmm5			; d1i,d1r,d0i,d0r
	movups		[ecx],xmm3			; store result		<(cc)->c[0]>

	; more special handling of the 3rd row of matrix a
	shufps		xmm1,xmm0,0x44			; b0r,-b0i,b1r,-b1i
	movss		xmm7,[eax+52]			; x,x,x,c20i		<(aa)->e[2][0].imag>
	movss		xmm5,[eax+60]			; x,x,x,c21i		<(aa)->e[2][1].imag>
	shufps		xmm5,xmm7,0x00			; c20i,c20i,c21i,c21i
	mulps		xmm5,xmm1
	addps		xmm6,xmm5
	shufps		xmm2,xmm2,0xB4			; -b2i,b2r,b2r,-b2i
	xorps		xmm2,[neg2]			; b2i,b2r,b2r,-b2i	<_sse_sgn4>
	movlps		xmm7,[eax+64]			; x,x,c22i,c22r		<(aa)->e[2][2]>
	shufps		xmm7,xmm7,0x05			; c22r,c22r,c22i,c22i
	mulps		xmm7,xmm2
	addps		xmm6,xmm7
	movaps		xmm7,xmm6
	shufps		xmm7,xmm7,0xEE
	addps		xmm6,xmm7
	movlps		[ecx+16],xmm6			;			<(cc)->c[2]>

	
	; *******************************************************************	

here:	pop	ecx
	pop	ebx
	pop	eax
	mov	esp,ebp
	pop	ebp
	ret
	
	align		16
negate:	dd		0x80000000
	dd		0x00000000
	dd		0x80000000
	dd		0x00000000
	
	align		16
neg2:   dd		0x00000000
	dd		0x00000000
	dd		0x00000000
	dd		0x80000000
	
