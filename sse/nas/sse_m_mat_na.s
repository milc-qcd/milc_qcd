;
; mult_su3_na( su3_matrix *a, su3_matrix *b, su3_matrix *c)
;

global mult_su3_na
mult_su3_na:
	push		ebp
	mov		ebp,esp
	push		eax
	push		ebx
	push		ecx
	mov		eax,[ebp+8]			; su3_matrix *a
	mov		ebx,[ebp+12]			; su3_vector *b
	mov		ecx,[ebp+16]			; su3_vector *c

	; bring in "half wilson vector" (first two columns)
	movlps		xmm0,[ebx]			;			<(bb)->e[0][0]>
	movlps		xmm1,[ebx+8]			;			<(bb)->e[0][1]>
	movlps		xmm2,[ebx+16]			;			<(bb)->e[0][2]>
	movhps		xmm0,[ebx+24]			; b[1]1 b[0]1		<(bb)->e[1][0]>
	movhps		xmm1,[ebx+32]			; b[1]2 b[0]2		<(bb)->e[1][1]>
	movhps		xmm2,[ebx+40]			; b[1]3 b[0]3		<(bb)->e[1][2]>

	; multiply vectors by real components of a
	movss		xmm3,[eax]			; a(1,1)		<(aa)->e[0][0].real>
	movss		xmm6,[eax+8]			; a(1,2)		<(aa)->e[0][1].real>
	movss		xmm4,[eax+24]			; a(2,1)		<(aa)->e[1][0].real>
	movss		xmm7,[eax+40]			; a(2,3)		<(aa)->e[1][2].real>
	movss		xmm5,[eax+48]			; a(3,1)		<(aa)->e[2][0].real>
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
	movss		xmm6,[eax+56]			; a(3,2)		<(aa)->e[2][1].real>
	movss		xmm7,[eax+16]			; a(1,3)		<(aa)->e[0][2].real>
	shufps		xmm6,xmm6,0x00
	shufps		xmm7,xmm7,0x00
	mulps		xmm6,xmm1
	mulps		xmm7,xmm2
	addps		xmm5,xmm6
	addps		xmm3,xmm7
	movss		xmm6,[eax+32]			; a(2,2)		<(aa)->e[1][1].real>
	movss		xmm7,[eax+64]			; a(3,3)		<(aa)->e[2][2].real>
	shufps		xmm6,xmm6,0x00
	shufps		xmm7,xmm7,0x00
	mulps		xmm6,xmm1
	mulps		xmm7,xmm2
	addps		xmm4,xmm6
	addps		xmm5,xmm7
	
	; multiply vectors by imaginary components of a
	movss		xmm6,[eax+4]			; a(1,1)		<(aa)->e[0][0].imag>
	movss		xmm7,[eax+36]			; a(2,2)		<(aa)->e[1][1].imag>
	shufps		xmm0,xmm0,0xb1
	shufps		xmm1,xmm1,0xb1
	shufps		xmm2,xmm2,0xb1
	shufps		xmm6,xmm6,0x00
	shufps		xmm7,xmm7,0x00
	xorps		xmm0,[negate]			;			<_sse_sgn24>
	xorps		xmm1,[negate]			;			<_sse_sgn24>
	xorps		xmm2,[negate]			;			<_sse_sgn24>
	mulps		xmm6,xmm0
	mulps		xmm7,xmm1
	addps		xmm3,xmm6
	addps		xmm4,xmm7
	movss		xmm6,[eax+68]			; a(3,3)		<(aa)->e[2][2].imag>
	movss		xmm7,[eax+28]			; a(2,1)		<(aa)->e[1][0].imag>
	shufps		xmm6,xmm6,0x00
	shufps		xmm7,xmm7,0x00
	mulps		xmm6,xmm2
	mulps		xmm7,xmm0
	addps		xmm5,xmm6
	addps		xmm4,xmm7
	movss		xmm6,[eax+12]			; a(1,2)		<(aa)->e[0][1].imag>
	movss		xmm7,[eax+52]			; a(3,1)		<(aa)->e[2][0].imag>
	shufps		xmm6,xmm6,0x00
	shufps		xmm7,xmm7,0x00
	mulps		xmm6,xmm1
	mulps		xmm7,xmm0
	addps		xmm3,xmm6
	addps		xmm5,xmm7
	movss		xmm0,[eax+20]			; a(1,3)		<(aa)->e[0][2].imag>
	movss		xmm6,[eax+60]			; a(3,2)		<(aa)->e[2][1].imag>
	movss		xmm7,[eax+44]			; a(2,3)		<(aa)->e[1][2].imag>
	shufps		xmm0,xmm0,0x00
	shufps		xmm6,xmm6,0x00
	shufps		xmm7,xmm7,0x00
	mulps		xmm0,xmm2
	mulps		xmm6,xmm1
	mulps		xmm7,xmm2
	addps		xmm3,xmm0
	addps		xmm5,xmm6
	addps		xmm4,xmm7
	xorps		xmm3,[negate]			;			<_sse_sgn24>
	xorps		xmm4,[negate]			;			<_sse_sgn24>
	xorps		xmm5,[negate]			;			<_sse_sgn24>

	; move results into output
	movlps		[ecx],xmm3			;			<(cc)->e[0][0]>
	movlps		[ecx+24],xmm4			;			<(cc)->e[1][0]>
	movlps		[ecx+48],xmm5			;			<(cc)->e[2][0]>
	movhps		[ecx+8],xmm3			;			<(cc)->e[0][1]>
	movhps		[ecx+32],xmm4			;			<(cc)->e[1][1]>
	movhps		[ecx+56],xmm5			;			<(cc)->e[2][1]>

	;  Now handle 3rd column
	movlps		xmm0,[ebx+48]			; x,x,b0i,b0r		<(bb)->e[2][0]>
	movlps		xmm1,[ebx+56]			; x,x,b1i,b1r		<(bb)->e[2][1]>
	movlps		xmm2,[ebx+64]			; x,x,b2i,b2r		<(bb)->e[2][2]>
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
	xorps		xmm0,[negate]			; b0r,-b0i,b0r,-b0i	<_sse_sgn24>
	shufps		xmm1,xmm1,0x11			; b1r,b1i,b1r,b1i
	xorps		xmm1,[negate]			; b1r,-b1i,b1r,-b1i	<_sse_sgn24>
	shufps		xmm2,xmm2,0xB1			; b2r,b2i,b2r,b2i
	xorps		xmm2,[negate]			; b2r,-b2i,b2r,-b2i	<_sse_sgn24>

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
	xorps		xmm3,[negate]			;			<_sse_sgn24>
	movlps		[ecx+16],xmm3			; store result		<(cc)->e[0][2]>
	movhps		[ecx+40],xmm3			;			<(cc)->e[1][2]>

	; more special handling of the 3rd row of matrix a
	shufps		xmm1,xmm0,0x44			; b0r,-b0i,b1r,-b1i
	movss		xmm7,[eax+52]			; x,x,x,c20i		<(aa)->e[2][0].imag>
	movss		xmm5,[eax+60]			; x,x,x,c21i		<(aa)->e[2][1].imag>
	shufps		xmm5,xmm7,0x00			; c20i,c20i,c21i,c21i
	mulps		xmm5,xmm1
	addps		xmm6,xmm5
	shufps		xmm2,xmm2,0xB4			; -b2i,b2r,b2r,-b2i
	xorps		xmm2,[neg2]			; b2i,b2r,b2r,-b2i	<_sse_sgn3>
	movlps		xmm7,[eax+64]			; x,x,c22i,c22r		<(aa)->e[2][2]>
	shufps		xmm7,xmm7,0x05			; c22r,c22r,c22i,c22i
	mulps		xmm7,xmm2
	addps		xmm6,xmm7
	movaps		xmm7,xmm6
	shufps		xmm7,xmm7,0xEE
	addps		xmm6,xmm7
	xorps		xmm6,[negate]			;			<_sse_sgn24>
	movlps		[ecx+64],xmm6			;			<(cc)->e[2][2]>

	
	; *******************************************************************	

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
	
        align           16
neg2:   dd              0x00000000
        dd              0x00000000
        dd              0x80000000
        dd              0x00000000
