;
; mult_su3_mat_vec_sum_4dir( su3_matrix *a, su3_vector *b0, su3_vector *b1, 
;                                su3_vector *b2, su3_vector *b3, su3_vector *c)
; 
; Multiply an array of 4 su4 matrices, a[4], by an array of vectors, b[4], summing
; the result into an su3 vector c.
;

global mult_su3_mat_vec_sum_4dir
mult_su3_mat_vec_sum_4dir:
	push		ebp
	mov		ebp,esp
	push		eax
	push		ebx
	push		ecx
	mov		eax,[ebp+8]			; su3_matrix *a
	mov		ebx,[ebp+12]			; su3_vector *b
	mov		ecx,[ebp+28]			; su3_vector *c

	;  bring in real and imaginary b[0] vector
	movlps		xmm0,[ebx]			; x,x,b0i,b0r		<(bb0)->c[0]>
	movlps		xmm1,[ebx+8]			; x,x,b1i,b1r		<(bb0)->c[1]>
	movlps		xmm2,[ebx+16]			; x,x,b2i,b2r		<(bb0)->c[2]>
	shufps		xmm0,xmm0,0x44			; b0i,b0r,b0i,b0r
	shufps		xmm1,xmm1,0x44			; b1i,b1r,b1i,b1r
	shufps		xmm2,xmm2,0x44			; b2i,b2r,b2i,b2r

	; bring in real components of first two rows of matrix a[0]
	movss		xmm3,[eax]			; x,x,x,c00r		<(aa)[0].e[0][0].real>
	movss		xmm7,[eax+24]			; x,x,x,c10r		<(aa)[0].e[1][0].real>
	shufps		xmm3,xmm7,0x00			; c10r,c10r,c00r,c00r
	movss		xmm4,[eax+8]			; x,x,x,c01r		<(aa)[0].e[0][1].real>
	movss		xmm7,[eax+32]			; x,x,x,c11r		<(aa)[0].e[1][1].real>
	shufps		xmm4,xmm7,0x00			; c11r,c11r,c01r,c01r
	mulps		xmm3,xmm0
	mulps		xmm4,xmm1
	addps		xmm3,xmm4
	movss		xmm5,[eax+16]			; x,x,x,c02r		<(aa)[0].e[0][2].real>
	movss		xmm7,[eax+40]			; x,x,x,c12r		<(aa)[0].e[1][2].real>
	shufps		xmm5,xmm7,0x00			; c12r,c12r,c02r,c02r
	mulps		xmm5,xmm2
	addps		xmm3,xmm5
	
	; special handling of the 3rd row of matrix a[0]
	shufps		xmm1,xmm0,0x44			; b0i,b0r,b1i,b1r
	movss		xmm7,[eax+48]			; x,x,x,c20r		<(aa)[0].e[2][0].real>
	movss		xmm6,[eax+56]			; x,x,x,c21r		<(aa)[0].e[2][1].real>
	shufps		xmm6,xmm7,0x00			; c20r,c20r,c21r,c21r
	mulps		xmm6,xmm1

	; shuffle b vector for imaginary components of matrix a[0]
	shufps		xmm0,xmm0,0xB1			; b0r,b0i,b0r,b0i
	xorps		xmm0,[negate]			; b0r,-b0i,b0r,-b0i	<_sse_sgn13>
	shufps		xmm1,xmm1,0x11			; b1r,b1i,b1r,b1i
	xorps		xmm1,[negate]			; b1r,-b1i,b1r,-b1i	<_sse_sgn13>
	shufps		xmm2,xmm2,0xB1			; b2r,b2i,b2r,b2i
	xorps		xmm2,[negate]			; b2r,-b2i,b2r,-b2i	<_sse_sgn13>

	; bring in imaginary components of first two rows of matrix b[0]
	movss		xmm4,[eax+4]			; x,x,x,c00i		<(aa)[0].e[0][0].imag>
	movss		xmm7,[eax+28]			; x,x,x,c10i		<(aa)[0].e[1][0].imag>
	shufps		xmm4,xmm7,0x00			; c10i,c10i,c00i,c00i
	mulps		xmm4,xmm0
	addps		xmm3,xmm4
	movss		xmm5,[eax+12]			; x,x,x,c01i		<(aa)[0].e[0][1].imag>
	movss		xmm7,[eax+36]			; x,x,x,c11i		<(aa)[0].e[1][1].imag>
	shufps		xmm5,xmm7,0x00			; c11i,c11i,c01i,c01i
	mulps		xmm5,xmm1
	addps		xmm3,xmm5
	movss		xmm5,[eax+20]			; x,x,x,c02i		<(aa)[0].e[0][2].imag>
	movss		xmm7,[eax+44]			; x,x,x,c12i		<(aa)[0].e[1][2].imag>
	shufps		xmm5,xmm7,0x00			; c12i,c12i,c02i,c02i
	mulps		xmm5,xmm2
	addps		xmm3,xmm5			; d1i,d1r,d0i,d0r
	movups		[ecx],xmm3			; store result		<(cc)->c[0]>

	; more special handling of the 3rd row of matrix a[0]
	shufps		xmm1,xmm0,0x44			; b0r,-b0i,b1r,-b1i
	movss		xmm7,[eax+52]			; x,x,x,c20i		<(aa)[0].e[2][0].imag>
	movss		xmm5,[eax+60]			; x,x,x,c21i		<(aa)[0].e[2][1].imag>
	shufps		xmm5,xmm7,0x00			; c20i,c20i,c21i,c21i
	mulps		xmm5,xmm1
	addps		xmm6,xmm5
	shufps		xmm2,xmm2,0xB4			; -b2i,b2r,b2r,-b2i
	xorps		xmm2,[neg2]			; b2i,b2r,b2r,-b2i	<_sse_sgn4>
	movlps		xmm7,[eax+64]			; x,x,c22i,c22r		<(aa)[0].e[2][2]>
	shufps		xmm7,xmm7,0x05			; c22r,c22r,c22i,c22i
	mulps		xmm7,xmm2
	addps		xmm6,xmm7
	movaps		xmm7,xmm6
	shufps		xmm7,xmm7,0xEE
	addps		xmm6,xmm7
	movlps		[ecx+16],xmm6			;			<(cc)->c[2]>
	
	; *******************************************************************	

        mov             ebx,[ebp+16]                    ; su3_vector *b
        add             eax,72                          ; su3_matrix *a[1]
	;  bring in real and imaginary b[1] vector
	movlps		xmm0,[ebx]			; x,x,b0i,b0r		<(bb1)->c[0]>
	movlps		xmm1,[ebx+8]			; x,x,b1i,b1r		<(bb1)->c[1]>
	movlps		xmm2,[ebx+16]			; x,x,b2i,b2r		<(bb1)->c[2]>
	shufps		xmm0,xmm0,0x44			; b0i,b0r,b0i,b0r
	shufps		xmm1,xmm1,0x44			; b1i,b1r,b1i,b1r
	shufps		xmm2,xmm2,0x44			; b2i,b2r,b2i,b2r

	; bring in real components of first two rows of matrix a[1]
	movss		xmm3,[eax]			; x,x,x,c00r		<(aa)[1].e[0][0].real>
	movss		xmm7,[eax+24]			; x,x,x,c10r		<(aa)[1].e[1][0].real>
	shufps		xmm3,xmm7,0x00			; c10r,c10r,c00r,c00r
	movss		xmm4,[eax+8]			; x,x,x,c01r		<(aa)[1].e[0][1].real>
	movss		xmm7,[eax+32]			; x,x,x,c11r		<(aa)[1].e[1][1].real>
	shufps		xmm4,xmm7,0x00			; c11r,c11r,c01r,c01r
	mulps		xmm3,xmm0
	mulps		xmm4,xmm1
	addps		xmm3,xmm4
	movss		xmm5,[eax+16]			; x,x,x,c02r		<(aa)[1].e[0][2].real>
	movss		xmm7,[eax+40]			; x,x,x,c12r		<(aa)[1].e[1][2].real>
	shufps		xmm5,xmm7,0x00			; c12r,c12r,c02r,c02r
	mulps		xmm5,xmm2
	addps		xmm3,xmm5

	; special handling of the 3rd row of matrix a[1]
	shufps		xmm1,xmm0,0x44			; b0i,b0r,b1i,b1r
	movss		xmm7,[eax+48]			; x,x,x,c20r		<(aa)[1].e[2][0].real>
	movss		xmm6,[eax+56]			; x,x,x,c21r		<(aa)[1].e[2][1].real>
	shufps		xmm6,xmm7,0x00			; c20r,c20r,c21r,c21r
	mulps		xmm6,xmm1

	; shuffle b vector for imaginary components of matrix a[1]
	shufps		xmm0,xmm0,0xB1			; b0r,b0i,b0r,b0i
	xorps		xmm0,[negate]			; b0r,-b0i,b0r,-b0i	<_sse_sgn13>
	shufps		xmm1,xmm1,0x11			; b1r,b1i,b1r,b1i
	xorps		xmm1,[negate]			; b1r,-b1i,b1r,-b1i	<_sse_sgn13>
	shufps		xmm2,xmm2,0xB1			; b2r,b2i,b2r,b2i
	xorps		xmm2,[negate]			; b2r,-b2i,b2r,-b2i	<_sse_sgn13>

	; bring in imaginary components of first two rows of matrix b[1]
	movss		xmm4,[eax+4]			; x,x,x,c00i		<(aa)[1].e[0][0].imag>
	movss		xmm7,[eax+28]			; x,x,x,c10i		<(aa)[1].e[1][0].imag>
	shufps		xmm4,xmm7,0x00			; c10i,c10i,c00i,c00i
	mulps		xmm4,xmm0
	addps		xmm3,xmm4
	movss		xmm5,[eax+12]			; x,x,x,c01i		<(aa)[1].e[0][1].imag>
	movss		xmm7,[eax+36]			; x,x,x,c11i		<(aa)[1].e[1][1].imag>
	shufps		xmm5,xmm7,0x00			; c11i,c11i,c01i,c01i
	mulps		xmm5,xmm1
	addps		xmm3,xmm5
	movss		xmm5,[eax+20]			; x,x,x,c02i		<(aa)[1].e[0][2].imag>
	movss		xmm7,[eax+44]			; x,x,x,c12i		<(aa)[1].e[1][2].imag>
	shufps		xmm5,xmm7,0x00			; c12i,c12i,c02i,c02i
	mulps		xmm5,xmm2
	addps		xmm3,xmm5			; d1i,d1r,d0i,d0r
	movups		xmm7,[ecx]			;			<(cc)->c[0]>
	addps		xmm7,xmm3
	movups		[ecx],xmm7			; store result		<(cc)->c[0]>

	; more special handling of the 3rd row of matrix a[1]
	shufps		xmm1,xmm0,0x44			; b0r,-b0i,b1r,-b1i
	movss		xmm7,[eax+52]			; x,x,x,c20i		<(aa)[1].e[2][0].imag>
	movss		xmm5,[eax+60]			; x,x,x,c21i		<(aa)[1].e[2][1].imag>
	shufps		xmm5,xmm7,0x00			; c20i,c20i,c21i,c21i
	mulps		xmm5,xmm1
	addps		xmm6,xmm5
	shufps		xmm2,xmm2,0xB4			; -b2i,b2r,b2r,-b2i
	xorps		xmm2,[neg2]			; b2i,b2r,b2r,-b2i	<_sse_sgn4>
	movlps		xmm7,[eax+64]			; x,x,c22i,c22r		<(aa)[1].e[2][2]>
	shufps		xmm7,xmm7,0x05			; c22r,c22r,c22i,c22i
	mulps		xmm7,xmm2
	addps		xmm6,xmm7
	movaps		xmm7,xmm6
	shufps		xmm7,xmm7,0xEE
	addps		xmm6,xmm7
	movlps		xmm7,[ecx+16]			;			<(cc)->c[2]>
	addps		xmm7,xmm6
	movlps		[ecx+16],xmm7			;			<(cc)->c[2]>

	
	; *******************************************************************	


        mov             ebx,[ebp+20]                    ; su3_vector *b
        add             eax,72                          ; su3_matrix *a[1]
	;  bring in real and imaginary b[2] vector
	movlps		xmm0,[ebx]			; x,x,b0i,b0r		<(bb2)->c[0]>
	movlps		xmm1,[ebx+8]			; x,x,b1i,b1r		<(bb2)->c[1]>
	movlps		xmm2,[ebx+16]			; x,x,b2i,b2r		<(bb2)->c[2]>
	shufps		xmm0,xmm0,0x44			; b0i,b0r,b0i,b0r
	shufps		xmm1,xmm1,0x44			; b1i,b1r,b1i,b1r
	shufps		xmm2,xmm2,0x44			; b2i,b2r,b2i,b2r

	; bring in real components of first two rows of matrix a[2]
	movss		xmm3,[eax]			; x,x,x,c00r		<(aa)[2].e[0][0].real>
	movss		xmm7,[eax+24]			; x,x,x,c10r		<(aa)[2].e[1][0].real>
	shufps		xmm3,xmm7,0x00			; c10r,c10r,c00r,c00r
	movss		xmm4,[eax+8]			; x,x,x,c01r		<(aa)[2].e[0][1].real>
	movss		xmm7,[eax+32]			; x,x,x,c11r		<(aa)[2].e[1][1].real>
	shufps		xmm4,xmm7,0x00			; c11r,c11r,c01r,c01r
	mulps		xmm3,xmm0
	mulps		xmm4,xmm1
	addps		xmm3,xmm4
	movss		xmm5,[eax+16]			; x,x,x,c02r		<(aa)[2].e[0][2].real>
	movss		xmm7,[eax+40]			; x,x,x,c12r		<(aa)[2].e[1][2].real>
	shufps		xmm5,xmm7,0x00			; c12r,c12r,c02r,c02r
	mulps		xmm5,xmm2
	addps		xmm3,xmm5

	; special handling of the 3rd row of matrix a[2]
	shufps		xmm1,xmm0,0x44			; b0i,b0r,b1i,b1r
	movss		xmm7,[eax+48]			; x,x,x,c20r		<(aa)[2].e[2][0].real>
	movss		xmm6,[eax+56]			; x,x,x,c21r		<(aa)[2].e[2][1].real>
	shufps		xmm6,xmm7,0x00			; c20r,c20r,c21r,c21r
	mulps		xmm6,xmm1

	; shuffle b vector for imaginary components of matrix a[2]
	shufps		xmm0,xmm0,0xB1			; b0r,b0i,b0r,b0i
	xorps		xmm0,[negate]			; b0r,-b0i,b0r,-b0i	<_sse_sgn13>
	shufps		xmm1,xmm1,0x11			; b1r,b1i,b1r,b1i
	xorps		xmm1,[negate]			; b1r,-b1i,b1r,-b1i	<_sse_sgn13>
	shufps		xmm2,xmm2,0xB1			; b2r,b2i,b2r,b2i
	xorps		xmm2,[negate]			; b2r,-b2i,b2r,-b2i	<_sse_sgn13>

	; bring in imaginary components of first two rows of matrix b[2]
	movss		xmm4,[eax+4]			; x,x,x,c00i		<(aa)[2].e[0][0].imag>
	movss		xmm7,[eax+28]			; x,x,x,c10i		<(aa)[2].e[1][0].imag>
	shufps		xmm4,xmm7,0x00			; c10i,c10i,c00i,c00i
	mulps		xmm4,xmm0
	addps		xmm3,xmm4
	movss		xmm5,[eax+12]			; x,x,x,c01i		<(aa)[2].e[0][1].imag>
	movss		xmm7,[eax+36]			; x,x,x,c11i		<(aa)[2].e[1][1].imag>
	shufps		xmm5,xmm7,0x00			; c11i,c11i,c01i,c01i
	mulps		xmm5,xmm1
	addps		xmm3,xmm5
	movss		xmm5,[eax+20]			; x,x,x,c02i		<(aa)[2].e[0][2].imag>
	movss		xmm7,[eax+44]			; x,x,x,c12i		<(aa)[2].e[1][2].imag>
	shufps		xmm5,xmm7,0x00			; c12i,c12i,c02i,c02i
	mulps		xmm5,xmm2
	addps		xmm3,xmm5			; d1i,d1r,d0i,d0r
	movups		xmm7,[ecx]			;			<(cc)->c[0]>
	addps		xmm7,xmm3
	movups		[ecx],xmm7			; store result		<(cc)->c[0]>

	; more special handling of the 3rd row of matrix a[2]
	shufps		xmm1,xmm0,0x44			; b0r,-b0i,b1r,-b1i
	movss		xmm7,[eax+52]			; x,x,x,c20i		<(aa)[2].e[2][0].imag>
	movss		xmm5,[eax+60]			; x,x,x,c21i		<(aa)[2].e[2][1].imag>
	shufps		xmm5,xmm7,0x00			; c20i,c20i,c21i,c21i
	mulps		xmm5,xmm1
	addps		xmm6,xmm5
	shufps		xmm2,xmm2,0xB4			; -b2i,b2r,b2r,-b2i
	xorps		xmm2,[neg2]			; b2i,b2r,b2r,-b2i	<_sse_sgn4>
	movlps		xmm7,[eax+64]			; x,x,c22i,c22r		<(aa)[2].e[2][2]>
	shufps		xmm7,xmm7,0x05			; c22r,c22r,c22i,c22i
	mulps		xmm7,xmm2
	addps		xmm6,xmm7
	movaps		xmm7,xmm6
	shufps		xmm7,xmm7,0xEE
	addps		xmm6,xmm7
	movlps		xmm7,[ecx+16]			;			<(cc)->c[2]>
	addps		xmm7,xmm6
	movlps		[ecx+16],xmm7			;			<(cc)->c[2]>

	
	; *******************************************************************	

	

        mov             ebx,[ebp+24]                    ; su3_vector *b
        add             eax,72                          ; su3_matrix *a[1]
	;  bring in real and imaginary b[3] vector
	movlps		xmm0,[ebx]			; x,x,b0i,b0r		<(bb3)->c[0]>
	movlps		xmm1,[ebx+8]			; x,x,b1i,b1r		<(bb3)->c[1]>
	movlps		xmm2,[ebx+16]			; x,x,b2i,b2r		<(bb3)->c[2]>
	shufps		xmm0,xmm0,0x44			; b0i,b0r,b0i,b0r
	shufps		xmm1,xmm1,0x44			; b1i,b1r,b1i,b1r
	shufps		xmm2,xmm2,0x44			; b2i,b2r,b2i,b2r

	; bring in real components of first two rows of matrix a[3]
	movss		xmm3,[eax]			; x,x,x,c00r		<(aa)[3].e[0][0].real>
	movss		xmm7,[eax+24]			; x,x,x,c10r		<(aa)[3].e[1][0].real>
	shufps		xmm3,xmm7,0x00			; c10r,c10r,c00r,c00r
	movss		xmm4,[eax+8]			; x,x,x,c01r		<(aa)[3].e[0][1].real>
	movss		xmm7,[eax+32]			; x,x,x,c11r		<(aa)[3].e[1][1].real>
	shufps		xmm4,xmm7,0x00			; c11r,c11r,c01r,c01r
	mulps		xmm3,xmm0
	mulps		xmm4,xmm1
	addps		xmm3,xmm4
	movss		xmm5,[eax+16]			; x,x,x,c02r		<(aa)[3].e[0][2].real>
	movss		xmm7,[eax+40]			; x,x,x,c12r		<(aa)[3].e[1][2].real>
	shufps		xmm5,xmm7,0x00			; c12r,c12r,c02r,c02r
	mulps		xmm5,xmm2
	addps		xmm3,xmm5

	; special handling of the 3rd row of matrix a[3]
	shufps		xmm1,xmm0,0x44			; b0i,b0r,b1i,b1r
	movss		xmm7,[eax+48]			; x,x,x,c20r		<(aa)[3].e[2][0].real>
	movss		xmm6,[eax+56]			; x,x,x,c21r		<(aa)[3].e[2][1].real>
	shufps		xmm6,xmm7,0x00			; c20r,c20r,c21r,c21r
	mulps		xmm6,xmm1

	; shuffle b vector for imaginary components of matrix a[3]
	shufps		xmm0,xmm0,0xB1			; b0r,b0i,b0r,b0i
	xorps		xmm0,[negate]			; b0r,-b0i,b0r,-b0i	<_sse_sgn13>
	shufps		xmm1,xmm1,0x11			; b1r,b1i,b1r,b1i
	xorps		xmm1,[negate]			; b1r,-b1i,b1r,-b1i	<_sse_sgn13>
	shufps		xmm2,xmm2,0xB1			; b2r,b2i,b2r,b2i
	xorps		xmm2,[negate]			; b2r,-b2i,b2r,-b2i	<_sse_sgn13>

	; bring in imaginary components of first two rows of matrix b[3]
	movss		xmm4,[eax+4]			; x,x,x,c00i		<(aa)[3].e[0][0].imag>
	movss		xmm7,[eax+28]			; x,x,x,c10i		<(aa)[3].e[1][0].imag>
	shufps		xmm4,xmm7,0x00			; c10i,c10i,c00i,c00i
	mulps		xmm4,xmm0
	addps		xmm3,xmm4
	movss		xmm5,[eax+12]			; x,x,x,c01i		<(aa)[3].e[0][1].imag>
	movss		xmm7,[eax+36]			; x,x,x,c11i		<(aa)[3].e[1][1].imag>
	shufps		xmm5,xmm7,0x00			; c11i,c11i,c01i,c01i
	mulps		xmm5,xmm1
	addps		xmm3,xmm5
	movss		xmm5,[eax+20]			; x,x,x,c02i		<(aa)[3].e[0][2].imag>
	movss		xmm7,[eax+44]			; x,x,x,c12i		<(aa)[3].e[1][2].imag>
	shufps		xmm5,xmm7,0x00			; c12i,c12i,c02i,c02i
	mulps		xmm5,xmm2
	addps		xmm3,xmm5			; d1i,d1r,d0i,d0r
	movups		xmm7,[ecx]			;			<(cc)->c[0]>
	addps		xmm7,xmm3
	movups		[ecx],xmm7			; store result		<(cc)->c[0]>

	; more special handling of the 3rd row of matrix a[3]
	shufps		xmm1,xmm0,0x44			; b0r,-b0i,b1r,-b1i
	movss		xmm7,[eax+52]			; x,x,x,c20i		<(aa)[3].e[2][0].imag>
	movss		xmm5,[eax+60]			; x,x,x,c21i		<(aa)[3].e[2][1].imag>
	shufps		xmm5,xmm7,0x00			; c20i,c20i,c21i,c21i
	mulps		xmm5,xmm1
	addps		xmm6,xmm5
	shufps		xmm2,xmm2,0xB4			; -b2i,b2r,b2r,-b2i
	xorps		xmm2,[neg2]			; b2i,b2r,b2r,-b2i	<_sse_sgn4>
	movlps		xmm7,[eax+64]			; x,x,c22i,c22r		<(aa)[3].e[2][2]>
	shufps		xmm7,xmm7,0x05			; c22r,c22r,c22i,c22i
	mulps		xmm7,xmm2
	addps		xmm6,xmm7
	movaps		xmm7,xmm6
	shufps		xmm7,xmm7,0xEE
	addps		xmm6,xmm7
	movlps		xmm7,[ecx+16]			;			<(cc)->c[2]>
	addps		xmm7,xmm6
	movlps		[ecx+16],xmm7			;			<(cc)->c[2]>
	
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
	
