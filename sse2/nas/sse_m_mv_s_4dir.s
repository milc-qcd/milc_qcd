;
; mult_su3_mat_vec_sum_4dir( su3_matrix *a, su3_vector *b0, su3_vector *b1,
;                                su3_vector *b2, su3_vector *b3, su3_vector *c)
;
; Multiply an array of 4 su4 matrices, a[4], by an array of vectors, b[4], summing
; the result into an su3 vector c.
;
 
global mult_su3_mat_vec_sum_4dir
mult_su3_mat_vec_sum_4dir:
	push	ebp
	mov	ebp, esp
	push	eax
	push	ebx
	push	ecx
	mov	eax, [ebp+8]		; su3_matrix *a
	mov	ebx, [ebp+12]		; su3_vector *b
	mov	ecx, [ebp+28]		; su3_vector *c

	movupd	xmm0,[ebx]		; <(bb0)->c[0]>
	movupd	xmm1,[ebx+16]		; <(bb0)->c[1]>
	movupd	xmm2,[ebx+32]		; <(bb0)->c[2]>

	movsd	xmm3,[eax]		; xmm3 = x,a00r	<aa[0].e[0][0].real>
	movsd	xmm6,[eax+16]		; xmm6 = x,a01r	<aa[0].e[0][1].real>
	movsd	xmm4,[eax+48]		; xmm4 = x,a10r	<aa[0].e[1][0].real>
	movsd	xmm7,[eax+80]		; xmm7 = x,a12r	<aa[0].e[1][2].real>
	movsd	xmm5,[eax+96]		; xmm5 = x,a20r	<aa[0].e[2][0].real>
	unpcklpd xmm3,xmm3		; xmm3 = a00r,a00r
	unpcklpd xmm6,xmm6		; xmm6 = a01r,a01r
	unpcklpd xmm4,xmm4		; xmm4 = a10r,a10r
	mulpd	xmm3,xmm0		; xmm3 = a00r*b0r,a00r*b0i
	unpcklpd xmm7,xmm7		; xmm7 = a12r,a12r
	mulpd 	xmm6,xmm1		; xmm6 = a01r*b1r,a01r*b1i
	unpcklpd xmm5,xmm5		; xmm5 = a20r,a20r
	mulpd 	xmm4,xmm0		; xmm4 = a10r*b0r,a10r*b0i
	addpd 	xmm3,xmm6		; xmm3 = a00r*b0r+a01r*b1r,a00r*b0i+a01r*b1i
	mulpd 	xmm7,xmm2		; xmm7 = a12r*b2r,a12r*b2i
	mulpd 	xmm5,xmm0		; xmm5 = a20r*b0r,a20r*b0i
	addpd 	xmm4,xmm7		; xmm4 = a10r*b0r+a12r*b2r,a10r*b0i+a12r*b2i
	movsd 	xmm6,[eax+112]		; xmm6 = x,a21r	<aa[0].e[2][1].real>
	movsd 	xmm7,[eax+32]		; xmm7 = x,a02r <aa[0].e[0][2].real>
	unpcklpd xmm6,xmm6		; xmm6 = a21r,a21r
	unpcklpd xmm7,xmm7		; xmm7 = a02r,a02r
	mulpd 	xmm6,xmm1		; xmm6 = a21r*b1r,a21r*b1i
	mulpd 	xmm7,xmm2		; xmm7 = a02r*b2r,a02r*b2i
	addpd 	xmm5,xmm6		; xmm5 = a20r*b0r+a21r*b1r,a20r*b0i+a21r*b1i
	addpd 	xmm3,xmm7		; xmm3 = a00r*b0r + a01r*b1r + a02r*b2r , a00r*b0i + a01r*b1i + a02r*b2i
	movsd 	xmm6,[eax+64]		; xmm6 = x,a11r	<aa[0].e[1][1].real>
	movsd 	xmm7,[eax+128]		; xmm7 = x,a22r	<aa[0].e[2][2].real>
	unpcklpd xmm6,xmm6		; xmm6 = a11r,a11r
	unpcklpd xmm7,xmm7		; xmm7 = a22r,a22r
	mulpd 	xmm6,xmm1   		; xmm6 = a11r*b1r,a11r*b1i
	mulpd 	xmm7,xmm2   		; xmm7 = a22r*b2r,a22r*b2i
	addpd 	xmm4,xmm6		; xmm4 = a10r*b0r + a12r*b2r + a11r*b1r , a10r*b0i + a12r*b2i + a11r*b1i
	addpd 	xmm5,xmm7		; xmm5 = a20r*b0r + a21r*b1r + a22r*b2r , a20r*b0i + a21r*b1i + a22r*b2i
	;	
	; IMAGINARY
	;
	movsd	xmm6,[eax+8]		; xmm6 = x,a00i	<aa[0].e[0][0].imag>
	movsd	xmm7,[eax+72]		; xmm7 = x,a11i	<aa[0].e[1][1].imag>
	shufpd	xmm0,xmm0,0x1		; xmm0 = b0i,b0r
	shufpd	xmm1,xmm1,0x1		; xmm1 = b1i,b1r
	shufpd	xmm2,xmm2,0x1		; xmm2 = b2i,b2r
	unpcklpd xmm6,xmm6		; xmm6 = a00i,a00i
	unpcklpd xmm7,xmm7		; xmm7 = a11i,a11i
	xorpd	xmm0,[negate]		; xmm0 = -b0i,b0r	<_sse_sgn2>
	xorpd	xmm1,[negate]		; xmm1 = -b1i,b1r	<_sse_sgn2>
	xorpd	xmm2,[negate]		; xmm2 = -b2i,b2r	<_sse_sgn2>
	mulpd	xmm6,xmm0		; xmm6 = -a00i*b0i,a00i*b0r
	mulpd	xmm7,xmm1		; xmm7 = -a11i*b1i,a11i*b1r
	addpd 	xmm3,xmm6		; xmm3 = a00r*b0r + a01r*b1r + a02r*b2r - a00i*b0i , a00r*b0i + a01r*b1i + a02r*b2i + a00i*b0r
	addpd 	xmm4,xmm7		; xmm4 = a10r*b0r + a12r*b2r + a11r*b1r - a11i*b1i , a10r*b0i + a12r*b2i + a11r*b1i + a11i*b1r
	movsd	xmm6,[eax+136]		; xmm6 = x,a22i	<aa[0].e[2][2].imag>
	movsd	xmm7,[eax+56]		; xmm7 = x,a10i	<aa[0].e[1][0].imag>
	unpcklpd xmm6,xmm6		; xmm6 = a22i,a22i
	unpcklpd xmm7,xmm7		; xmm7 = a10i,a10i
	mulpd 	xmm6,xmm2		; xmm6 = -a22i*b2i,a22i*b2r
	mulpd	xmm7,xmm0		; xmm7 = -a10i*b0i,a10i*b0r
	addpd 	xmm5,xmm6		; xmm5 = a20r*b0r + a21r*b1r + a22r*b2r - a22i*b2i , a20r*b0i + a21r*b1i + a22r*b2i + a22i*b2r
	addpd 	xmm4,xmm7		; xmm4 = a10r*b0r + a12r*b2r + a11r*b1r - a11i*b1i - a10i*b0i, a10r*b0i + a12r*b2i + a11r*b1i + a11i*b1r + a10i*b0r
	movsd	xmm6,[eax+24]		; xmm6 = x,a01i	<aa[0].e[0][1].imag>
	movsd	xmm7,[eax+104]		; xmm7 = x,a20i	<aa[0].e[2][0].imag>
	unpcklpd xmm6,xmm6		; xmm6 = a01i,a01i
	unpcklpd xmm7,xmm7		; xmm7 = a20i,a20i
	mulpd 	xmm6,xmm1		; xmm6 = -a01i*b1i,a01i*b1r
	mulpd 	xmm7,xmm0		; xmm7 = -a20i*b0i,a20i*b0r
	addpd 	xmm3,xmm6		; xmm3 = a00r*b0r + a01r*b1r + a02r*b2r - a00i*b0i - a01i*b1i , a00r*b0i + a01r*b1i + a02r*b2i + a00i*b0r + a01i*b1r
	addpd 	xmm5,xmm7		; xmm5 = a20r*b0r + a21r*b1r + a22r*b2r - a22i*b2i - a20i*b0i , a20r*b0i + a21r*b1i + a22r*b2i + a22i*b2r + a20i*b0r
	movsd 	xmm0,[eax+40]		; xmm0 = x,a02i	<aa[0].e[0][2].imag>
	movsd	xmm6,[eax+120]		; xmm6 = x,a21i	<aa[0].e[2][1].imag>
	movsd 	xmm7,[eax+88]		; xmm7 = x,a12i	<aa[0].e[1][2].imag>
	unpcklpd xmm0,xmm0		; xmm0 = a02i,a02i
	unpcklpd xmm6,xmm6		; xmm6 = a21i,a21i
	unpcklpd xmm7,xmm7		; xmm7 = a12i,a12i
	mulpd 	xmm0,xmm2		; xmm0 = -a02i*b2i,a02i*b2r
	mulpd	xmm6,xmm1		; xmm6 = -a21i*b1i,a21i*b1r
	mulpd	xmm7,xmm2		; xmm7 = -a12i*b2i,a12i*b2r
	addpd 	xmm3,xmm0		; xmm3 = a00r*b0r + a01r*b1r + a02r*b2r - a00i*b0i - a01i*b1i - a02i*b2i , a00r*b0i + a01r*b1i + a02r*b2i + a00i*b0r + a01i*b1r + a02i*b2r
		                        ;        a00r*b0r - a00i*b0i + a01r*b1r - a01i*b1i + a02r*b2r - a02i*b2i , a00r*b0i + a00i*b0r + a01r*b1i + a01i*b1r + a02r*b2i + a02i*b2r 
	addpd 	xmm4,xmm7		; xmm4 = a10r*b0r + a12r*b2r + a11r*b1r - a11i*b1i - a10i*b0i - a12i*b2i , a10r*b0i + a12r*b2i + a11r*b1i + a11i*b1r + a10i*b0r + a12i*b2r
					;        a10r*b0r - a10i*b0i + a11r*b1r - a11i*b1i + a12r*b2r - a12i*b2i , a10r*b0i + a10i*b0r + a11r*b1i + a11i*b1r + a12r*b2i + a12i*b2r
	addpd 	xmm5,xmm6		; xmm5 = a20r*b0r + a21r*b1r + a22r*b2r - a22i*b2i - a20i*b0i - a21i*b1i , a20r*b0i + a21r*b1i + a22r*b2i + a22i*b2r + a20i*b0r + a21i*b1r
        	                        ;        a20r*b0r - a20i*b0i + a21r*b1r - a21i*b1i + a22r*b2r - a22i*b2i , a20r*b0i + a20i*b0r + a21r*b1i + a21i*b1r + a22r*b2i + a22i*b2r 
	movupd	[ecx],xmm3		;	<(cc)->c[0]>
	movupd	[ecx+16],xmm4		;	<(cc)->c[1]>
	movupd	[ecx+32],xmm5		;	<(cc)->c[2]>
	;
	; second
	;
	mov	ebx,[ebp+16]		; su3_vector *b1
	add	eax,144			; su3_matrix *a[1]
	;
	movupd	xmm0,[ebx]		; <(bb1)->c[0]>
	movupd	xmm1,[ebx+16]		; <(bb1)->c[1]>
	movupd	xmm2,[ebx+32]		; <(bb1)->c[2]>

	movsd	xmm3,[eax]		; xmm3 = x,a00r	<aa[1].e[0][0].real>
	movsd	xmm6,[eax+16]		; xmm6 = x,a01r	<aa[1].e[0][1].real>
	movsd	xmm4,[eax+48]		; xmm4 = x,a10r	<aa[1].e[1][0].real>
	movsd	xmm7,[eax+80]		; xmm7 = x,a12r	<aa[1].e[1][2].real>
	movsd	xmm5,[eax+96]		; xmm5 = x,a20r	<aa[1].e[2][0].real>
	unpcklpd xmm3,xmm3		; xmm3 = a00r,a00r
	unpcklpd xmm6,xmm6		; xmm6 = a01r,a01r
	unpcklpd xmm4,xmm4		; xmm4 = a10r,a10r
	mulpd	xmm3,xmm0		; xmm3 = a00r*b0r,a00r*b0i
	unpcklpd xmm7,xmm7		; xmm7 = a12r,a12r
	mulpd 	xmm6,xmm1		; xmm6 = a01r*b1r,a01r*b1i
	unpcklpd xmm5,xmm5		; xmm5 = a20r,a20r
	mulpd 	xmm4,xmm0		; xmm4 = a10r*b0r,a10r*b0i
	addpd 	xmm3,xmm6		; xmm3 = a00r*b0r+a01r*b1r,a00r*b0i+a01r*b1i
	mulpd 	xmm7,xmm2		; xmm7 = a12r*b2r,a12r*b2i
	mulpd 	xmm5,xmm0		; xmm5 = a20r*b0r,a20r*b0i
	addpd 	xmm4,xmm7		; xmm4 = a10r*b0r+a12r*b2r,a10r*b0i+a12r*b2i
	movsd 	xmm6,[eax+112]		; xmm6 = x,a21r	<aa[1].e[2][1].real>
	movsd 	xmm7,[eax+32]		; xmm7 = x,a02r <aa[1].e[0][2].real>
	unpcklpd xmm6,xmm6		; xmm6 = a21r,a21r
	unpcklpd xmm7,xmm7		; xmm7 = a02r,a02r
	mulpd 	xmm6,xmm1		; xmm6 = a21r*b1r,a21r*b1i
	mulpd 	xmm7,xmm2		; xmm7 = a02r*b2r,a02r*b2i
	addpd 	xmm5,xmm6		; xmm5 = a20r*b0r+a21r*b1r,a20r*b0i+a21r*b1i
	addpd 	xmm3,xmm7		; xmm3 = a00r*b0r + a01r*b1r + a02r*b2r , a00r*b0i + a01r*b1i + a02r*b2i
	movsd 	xmm6,[eax+64]		; xmm6 = x,a11r	<aa[1].e[1][1].real>
	movsd 	xmm7,[eax+128]		; xmm7 = x,a22r	<aa[1].e[2][2].real>
	unpcklpd xmm6,xmm6		; xmm6 = a11r,a11r
	unpcklpd xmm7,xmm7		; xmm7 = a22r,a22r
	mulpd 	xmm6,xmm1   		; xmm6 = a11r*b1r,a11r*b1i
	mulpd 	xmm7,xmm2   		; xmm7 = a22r*b2r,a22r*b2i
	addpd 	xmm4,xmm6		; xmm4 = a10r*b0r + a12r*b2r + a11r*b1r , a10r*b0i + a12r*b2i + a11r*b1i
	addpd 	xmm5,xmm7		; xmm5 = a20r*b0r + a21r*b1r + a22r*b2r , a20r*b0i + a21r*b1i + a22r*b2i
	;	
	; IMAGINARY
	;
	movsd	xmm6,[eax+8]		; xmm6 = x,a00i	<aa[1].e[0][0].imag>
	movsd	xmm7,[eax+72]		; xmm7 = x,a11i	<aa[1].e[1][1].imag>
	shufpd	xmm0,xmm0,0x1		; xmm0 = b0i,b0r
	shufpd	xmm1,xmm1,0x1		; xmm1 = b1i,b1r
	shufpd	xmm2,xmm2,0x1		; xmm2 = b2i,b2r
	unpcklpd xmm6,xmm6		; xmm6 = a00i,a00i
	unpcklpd xmm7,xmm7		; xmm7 = a11i,a11i
	xorpd	xmm0,[negate]		; xmm0 = -b0i,b0r	<_sse_sgn2>
	xorpd	xmm1,[negate]		; xmm1 = -b1i,b1r	<_sse_sgn2>
	xorpd	xmm2,[negate]		; xmm2 = -b2i,b2r	<_sse_sgn2>
	mulpd	xmm6,xmm0		; xmm6 = -a00i*b0i,a00i*b0r
	mulpd	xmm7,xmm1		; xmm7 = -a11i*b1i,a11i*b1r
	addpd 	xmm3,xmm6		; xmm3 = a00r*b0r + a01r*b1r + a02r*b2r - a00i*b0i , a00r*b0i + a01r*b1i + a02r*b2i + a00i*b0r
	addpd 	xmm4,xmm7		; xmm4 = a10r*b0r + a12r*b2r + a11r*b1r - a11i*b1i , a10r*b0i + a12r*b2i + a11r*b1i + a11i*b1r
	movsd	xmm6,[eax+136]		; xmm6 = x,a22i	<aa[1].e[2][2].imag>
	movsd	xmm7,[eax+56]		; xmm7 = x,a10i	<aa[1].e[1][0].imag>
	unpcklpd xmm6,xmm6		; xmm6 = a22i,a22i
	unpcklpd xmm7,xmm7		; xmm7 = a10i,a10i
	mulpd 	xmm6,xmm2		; xmm6 = -a22i*b2i,a22i*b2r
	mulpd	xmm7,xmm0		; xmm7 = -a10i*b0i,a10i*b0r
	addpd 	xmm5,xmm6		; xmm5 = a20r*b0r + a21r*b1r + a22r*b2r - a22i*b2i , a20r*b0i + a21r*b1i + a22r*b2i + a22i*b2r
	addpd 	xmm4,xmm7		; xmm4 = a10r*b0r + a12r*b2r + a11r*b1r - a11i*b1i - a10i*b0i, a10r*b0i + a12r*b2i + a11r*b1i + a11i*b1r + a10i*b0r
	movsd	xmm6,[eax+24]		; xmm6 = x,a01i	<aa[1].e[0][1].imag>
	movsd	xmm7,[eax+104]		; xmm7 = x,a20i	<aa[1].e[2][0].imag>
	unpcklpd xmm6,xmm6		; xmm6 = a01i,a01i
	unpcklpd xmm7,xmm7		; xmm7 = a20i,a20i
	mulpd 	xmm6,xmm1		; xmm6 = -a01i*b1i,a01i*b1r
	mulpd 	xmm7,xmm0		; xmm7 = -a20i*b0i,a20i*b0r
	addpd 	xmm3,xmm6		; xmm3 = a00r*b0r + a01r*b1r + a02r*b2r - a00i*b0i - a01i*b1i , a00r*b0i + a01r*b1i + a02r*b2i + a00i*b0r + a01i*b1r
	addpd 	xmm5,xmm7		; xmm5 = a20r*b0r + a21r*b1r + a22r*b2r - a22i*b2i - a20i*b0i , a20r*b0i + a21r*b1i + a22r*b2i + a22i*b2r + a20i*b0r
	movsd 	xmm0,[eax+40]		; xmm0 = x,a02i	<aa[1].e[0][2].imag>
	movsd	xmm6,[eax+120]		; xmm6 = x,a21i	<aa[1].e[2][1].imag>
	movsd 	xmm7,[eax+88]		; xmm7 = x,a12i	<aa[1].e[1][2].imag>
	unpcklpd xmm0,xmm0		; xmm0 = a02i,a02i
	unpcklpd xmm6,xmm6		; xmm6 = a21i,a21i
	unpcklpd xmm7,xmm7		; xmm7 = a12i,a12i
	mulpd 	xmm0,xmm2		; xmm0 = -a02i*b2i,a02i*b2r
	mulpd	xmm6,xmm1		; xmm6 = -a21i*b1i,a21i*b1r
	mulpd	xmm7,xmm2		; xmm7 = -a12i*b2i,a12i*b2r
	addpd 	xmm3,xmm0		; xmm3 = a00r*b0r + a01r*b1r + a02r*b2r - a00i*b0i - a01i*b1i - a02i*b2i , a00r*b0i + a01r*b1i + a02r*b2i + a00i*b0r + a01i*b1r + a02i*b2r
		                        ;        a00r*b0r - a00i*b0i + a01r*b1r - a01i*b1i + a02r*b2r - a02i*b2i , a00r*b0i + a00i*b0r + a01r*b1i + a01i*b1r + a02r*b2i + a02i*b2r 
	addpd 	xmm4,xmm7		; xmm4 = a10r*b0r + a12r*b2r + a11r*b1r - a11i*b1i - a10i*b0i - a12i*b2i , a10r*b0i + a12r*b2i + a11r*b1i + a11i*b1r + a10i*b0r + a12i*b2r
					;        a10r*b0r - a10i*b0i + a11r*b1r - a11i*b1i + a12r*b2r - a12i*b2i , a10r*b0i + a10i*b0r + a11r*b1i + a11i*b1r + a12r*b2i + a12i*b2r
	addpd 	xmm5,xmm6		; xmm5 = a20r*b0r + a21r*b1r + a22r*b2r - a22i*b2i - a20i*b0i - a21i*b1i , a20r*b0i + a21r*b1i + a22r*b2i + a22i*b2r + a20i*b0r + a21i*b1r
        	                        ;        a20r*b0r - a20i*b0i + a21r*b1r - a21i*b1i + a22r*b2r - a22i*b2i , a20r*b0i + a20i*b0r + a21r*b1i + a21i*b1r + a22r*b2i + a22i*b2r 
	movupd	xmm0,[ecx]              ;       <(cc)->c[0]>
	movupd	xmm6,[ecx+16]		; 	<(cc)->c[1]>
	movupd	xmm7,[ecx+32]		; 	<(cc)->c[2]>
	addpd	xmm3,xmm0
	addpd	xmm4,xmm6
	addpd	xmm5,xmm7
	movupd	[ecx],xmm3		;	<(cc)->c[0]>
	movupd	[ecx+16],xmm4		;	<(cc)->c[1]>
	movupd	[ecx+32],xmm5		;	<(cc)->c[2]>
	;
	; third
	;
	mov	ebx,[ebp+20]		; su3_vector *b
	add	eax,144			; su3_matrix *a[2]
	;
	movupd	xmm0,[ebx]		; <(bb2)->c[0]>
	movupd	xmm1,[ebx+16]		; <(bb2)->c[1]>
	movupd	xmm2,[ebx+32]		; <(bb2)->c[2]>

	movsd	xmm3,[eax]		; xmm3 = x,a00r	<aa[2].e[0][0].real>
	movsd	xmm6,[eax+16]		; xmm6 = x,a01r	<aa[2].e[0][1].real>
	movsd	xmm4,[eax+48]		; xmm4 = x,a10r	<aa[2].e[1][0].real>
	movsd	xmm7,[eax+80]		; xmm7 = x,a12r	<aa[2].e[1][2].real>
	movsd	xmm5,[eax+96]		; xmm5 = x,a20r	<aa[2].e[2][0].real>
	unpcklpd xmm3,xmm3		; xmm3 = a00r,a00r
	unpcklpd xmm6,xmm6		; xmm6 = a01r,a01r
	unpcklpd xmm4,xmm4		; xmm4 = a10r,a10r
	mulpd	xmm3,xmm0		; xmm3 = a00r*b0r,a00r*b0i
	unpcklpd xmm7,xmm7		; xmm7 = a12r,a12r
	mulpd 	xmm6,xmm1		; xmm6 = a01r*b1r,a01r*b1i
	unpcklpd xmm5,xmm5		; xmm5 = a20r,a20r
	mulpd 	xmm4,xmm0		; xmm4 = a10r*b0r,a10r*b0i
	addpd 	xmm3,xmm6		; xmm3 = a00r*b0r+a01r*b1r,a00r*b0i+a01r*b1i
	mulpd 	xmm7,xmm2		; xmm7 = a12r*b2r,a12r*b2i
	mulpd 	xmm5,xmm0		; xmm5 = a20r*b0r,a20r*b0i
	addpd 	xmm4,xmm7		; xmm4 = a10r*b0r+a12r*b2r,a10r*b0i+a12r*b2i
	movsd 	xmm6,[eax+112]		; xmm6 = x,a21r	<aa[2].e[2][1].real>
	movsd 	xmm7,[eax+32]		; xmm7 = x,a02r <aa[2].e[0][2].real>
	unpcklpd xmm6,xmm6		; xmm6 = a21r,a21r
	unpcklpd xmm7,xmm7		; xmm7 = a02r,a02r
	mulpd 	xmm6,xmm1		; xmm6 = a21r*b1r,a21r*b1i
	mulpd 	xmm7,xmm2		; xmm7 = a02r*b2r,a02r*b2i
	addpd 	xmm5,xmm6		; xmm5 = a20r*b0r+a21r*b1r,a20r*b0i+a21r*b1i
	addpd 	xmm3,xmm7		; xmm3 = a00r*b0r + a01r*b1r + a02r*b2r , a00r*b0i + a01r*b1i + a02r*b2i
	movsd 	xmm6,[eax+64]		; xmm6 = x,a11r	<aa[2].e[1][1].real>
	movsd 	xmm7,[eax+128]		; xmm7 = x,a22r	<aa[2].e[2][2].real>
	unpcklpd xmm6,xmm6		; xmm6 = a11r,a11r
	unpcklpd xmm7,xmm7		; xmm7 = a22r,a22r
	mulpd 	xmm6,xmm1   		; xmm6 = a11r*b1r,a11r*b1i
	mulpd 	xmm7,xmm2   		; xmm7 = a22r*b2r,a22r*b2i
	addpd 	xmm4,xmm6		; xmm4 = a10r*b0r + a12r*b2r + a11r*b1r , a10r*b0i + a12r*b2i + a11r*b1i
	addpd 	xmm5,xmm7		; xmm5 = a20r*b0r + a21r*b1r + a22r*b2r , a20r*b0i + a21r*b1i + a22r*b2i
	;	
	; IMAGINARY
	;
	movsd	xmm6,[eax+8]		; xmm6 = x,a00i	<aa[2].e[0][0].imag>
	movsd	xmm7,[eax+72]		; xmm7 = x,a11i	<aa[2].e[1][1].imag>
	shufpd	xmm0,xmm0,0x1		; xmm0 = b0i,b0r
	shufpd	xmm1,xmm1,0x1		; xmm1 = b1i,b1r
	shufpd	xmm2,xmm2,0x1		; xmm2 = b2i,b2r
	unpcklpd xmm6,xmm6		; xmm6 = a00i,a00i
	unpcklpd xmm7,xmm7		; xmm7 = a11i,a11i
	xorpd	xmm0,[negate]		; xmm0 = -b0i,b0r	<_sse_sgn2>
	xorpd	xmm1,[negate]		; xmm1 = -b1i,b1r	<_sse_sgn2>
	xorpd	xmm2,[negate]		; xmm2 = -b2i,b2r	<_sse_sgn2>
	mulpd	xmm6,xmm0		; xmm6 = -a00i*b0i,a00i*b0r
	mulpd	xmm7,xmm1		; xmm7 = -a11i*b1i,a11i*b1r
	addpd 	xmm3,xmm6		; xmm3 = a00r*b0r + a01r*b1r + a02r*b2r - a00i*b0i , a00r*b0i + a01r*b1i + a02r*b2i + a00i*b0r
	addpd 	xmm4,xmm7		; xmm4 = a10r*b0r + a12r*b2r + a11r*b1r - a11i*b1i , a10r*b0i + a12r*b2i + a11r*b1i + a11i*b1r
	movsd	xmm6,[eax+136]		; xmm6 = x,a22i	<aa[2].e[2][2].imag>
	movsd	xmm7,[eax+56]		; xmm7 = x,a10i	<aa[2].e[1][0].imag>
	unpcklpd xmm6,xmm6		; xmm6 = a22i,a22i
	unpcklpd xmm7,xmm7		; xmm7 = a10i,a10i
	mulpd 	xmm6,xmm2		; xmm6 = -a22i*b2i,a22i*b2r
	mulpd	xmm7,xmm0		; xmm7 = -a10i*b0i,a10i*b0r
	addpd 	xmm5,xmm6		; xmm5 = a20r*b0r + a21r*b1r + a22r*b2r - a22i*b2i , a20r*b0i + a21r*b1i + a22r*b2i + a22i*b2r
	addpd 	xmm4,xmm7		; xmm4 = a10r*b0r + a12r*b2r + a11r*b1r - a11i*b1i - a10i*b0i, a10r*b0i + a12r*b2i + a11r*b1i + a11i*b1r + a10i*b0r
	movsd	xmm6,[eax+24]		; xmm6 = x,a01i	<aa[2].e[0][1].imag>
	movsd	xmm7,[eax+104]		; xmm7 = x,a20i	<aa[2].e[2][0].imag>
	unpcklpd xmm6,xmm6		; xmm6 = a01i,a01i
	unpcklpd xmm7,xmm7		; xmm7 = a20i,a20i
	mulpd 	xmm6,xmm1		; xmm6 = -a01i*b1i,a01i*b1r
	mulpd 	xmm7,xmm0		; xmm7 = -a20i*b0i,a20i*b0r
	addpd 	xmm3,xmm6		; xmm3 = a00r*b0r + a01r*b1r + a02r*b2r - a00i*b0i - a01i*b1i , a00r*b0i + a01r*b1i + a02r*b2i + a00i*b0r + a01i*b1r
	addpd 	xmm5,xmm7		; xmm5 = a20r*b0r + a21r*b1r + a22r*b2r - a22i*b2i - a20i*b0i , a20r*b0i + a21r*b1i + a22r*b2i + a22i*b2r + a20i*b0r
	movsd 	xmm0,[eax+40]		; xmm0 = x,a02i	<aa[2].e[0][2].imag>
	movsd	xmm6,[eax+120]		; xmm6 = x,a21i	<aa[2].e[2][1].imag>
	movsd 	xmm7,[eax+88]		; xmm7 = x,a12i	<aa[2].e[1][2].imag>
	unpcklpd xmm0,xmm0		; xmm0 = a02i,a02i
	unpcklpd xmm6,xmm6		; xmm6 = a21i,a21i
	unpcklpd xmm7,xmm7		; xmm7 = a12i,a12i
	mulpd 	xmm0,xmm2		; xmm0 = -a02i*b2i,a02i*b2r
	mulpd	xmm6,xmm1		; xmm6 = -a21i*b1i,a21i*b1r
	mulpd	xmm7,xmm2		; xmm7 = -a12i*b2i,a12i*b2r
	addpd 	xmm3,xmm0		; xmm3 = a00r*b0r + a01r*b1r + a02r*b2r - a00i*b0i - a01i*b1i - a02i*b2i , a00r*b0i + a01r*b1i + a02r*b2i + a00i*b0r + a01i*b1r + a02i*b2r
		                        ;        a00r*b0r - a00i*b0i + a01r*b1r - a01i*b1i + a02r*b2r - a02i*b2i , a00r*b0i + a00i*b0r + a01r*b1i + a01i*b1r + a02r*b2i + a02i*b2r 
	addpd 	xmm4,xmm7		; xmm4 = a10r*b0r + a12r*b2r + a11r*b1r - a11i*b1i - a10i*b0i - a12i*b2i , a10r*b0i + a12r*b2i + a11r*b1i + a11i*b1r + a10i*b0r + a12i*b2r
					;        a10r*b0r - a10i*b0i + a11r*b1r - a11i*b1i + a12r*b2r - a12i*b2i , a10r*b0i + a10i*b0r + a11r*b1i + a11i*b1r + a12r*b2i + a12i*b2r
	addpd 	xmm5,xmm6		; xmm5 = a20r*b0r + a21r*b1r + a22r*b2r - a22i*b2i - a20i*b0i - a21i*b1i , a20r*b0i + a21r*b1i + a22r*b2i + a22i*b2r + a20i*b0r + a21i*b1r
        	                        ;        a20r*b0r - a20i*b0i + a21r*b1r - a21i*b1i + a22r*b2r - a22i*b2i , a20r*b0i + a20i*b0r + a21r*b1i + a21i*b1r + a22r*b2i + a22i*b2r 
        movupd  xmm0,[ecx]		;	<(cc)->c[0]>
        movupd  xmm6,[ecx+16]		; 	<(cc)->c[1]>
        movupd  xmm7,[ecx+32]		; 	<(cc)->c[2]>
        addpd   xmm3,xmm0
        addpd   xmm4,xmm6
        addpd   xmm5,xmm7
        movupd  [ecx],xmm3              ;       <(cc)->c[0]>
        movupd  [ecx+16],xmm4           ;       <(cc)->c[1]>
        movupd  [ecx+32],xmm5           ;       <(cc)->c[2]>
	;
	; fourth
	;
	mov	ebx,[ebp+24]		; su3_vector *b
	add	eax,144			; su3_matrix *a[1]
	;
	movupd	xmm0,[ebx]		; <(bb3)->c[0]>
	movupd	xmm1,[ebx+16]		; <(bb3)->c[1]>
	movupd	xmm2,[ebx+32]		; <(bb3)->c[2]>

	movsd	xmm3,[eax]		; xmm3 = x,a00r	<aa[3].e[0][0].real>
	movsd	xmm6,[eax+16]		; xmm6 = x,a01r	<aa[3].e[0][1].real>
	movsd	xmm4,[eax+48]		; xmm4 = x,a10r	<aa[3].e[1][0].real>
	movsd	xmm7,[eax+80]		; xmm7 = x,a12r	<aa[3].e[1][2].real>
	movsd	xmm5,[eax+96]		; xmm5 = x,a20r	<aa[3].e[2][0].real>
	unpcklpd xmm3,xmm3		; xmm3 = a00r,a00r
	unpcklpd xmm6,xmm6		; xmm6 = a01r,a01r
	unpcklpd xmm4,xmm4		; xmm4 = a10r,a10r
	mulpd	xmm3,xmm0		; xmm3 = a00r*b0r,a00r*b0i
	unpcklpd xmm7,xmm7		; xmm7 = a12r,a12r
	mulpd 	xmm6,xmm1		; xmm6 = a01r*b1r,a01r*b1i
	unpcklpd xmm5,xmm5		; xmm5 = a20r,a20r
	mulpd 	xmm4,xmm0		; xmm4 = a10r*b0r,a10r*b0i
	addpd 	xmm3,xmm6		; xmm3 = a00r*b0r+a01r*b1r,a00r*b0i+a01r*b1i
	mulpd 	xmm7,xmm2		; xmm7 = a12r*b2r,a12r*b2i
	mulpd 	xmm5,xmm0		; xmm5 = a20r*b0r,a20r*b0i
	addpd 	xmm4,xmm7		; xmm4 = a10r*b0r+a12r*b2r,a10r*b0i+a12r*b2i
	movsd 	xmm6,[eax+112]		; xmm6 = x,a21r	<aa[3].e[2][1].real>
	movsd 	xmm7,[eax+32]		; xmm7 = x,a02r <aa[3].e[0][2].real>
	unpcklpd xmm6,xmm6		; xmm6 = a21r,a21r
	unpcklpd xmm7,xmm7		; xmm7 = a02r,a02r
	mulpd 	xmm6,xmm1		; xmm6 = a21r*b1r,a21r*b1i
	mulpd 	xmm7,xmm2		; xmm7 = a02r*b2r,a02r*b2i
	addpd 	xmm5,xmm6		; xmm5 = a20r*b0r+a21r*b1r,a20r*b0i+a21r*b1i
	addpd 	xmm3,xmm7		; xmm3 = a00r*b0r + a01r*b1r + a02r*b2r , a00r*b0i + a01r*b1i + a02r*b2i
	movsd 	xmm6,[eax+64]		; xmm6 = x,a11r	<aa[3].e[1][1].real>
	movsd 	xmm7,[eax+128]		; xmm7 = x,a22r	<aa[3].e[2][2].real>
	unpcklpd xmm6,xmm6		; xmm6 = a11r,a11r
	unpcklpd xmm7,xmm7		; xmm7 = a22r,a22r
	mulpd 	xmm6,xmm1   		; xmm6 = a11r*b1r,a11r*b1i
	mulpd 	xmm7,xmm2   		; xmm7 = a22r*b2r,a22r*b2i
	addpd 	xmm4,xmm6		; xmm4 = a10r*b0r + a12r*b2r + a11r*b1r , a10r*b0i + a12r*b2i + a11r*b1i
	addpd 	xmm5,xmm7		; xmm5 = a20r*b0r + a21r*b1r + a22r*b2r , a20r*b0i + a21r*b1i + a22r*b2i
	;	
	; IMAGINARY
	;
	movsd	xmm6,[eax+8]		; xmm6 = x,a00i	<aa[3].e[0][0].imag>
	movsd	xmm7,[eax+72]		; xmm7 = x,a11i	<aa[3].e[1][1].imag>
	shufpd	xmm0,xmm0,0x1		; xmm0 = b0i,b0r
	shufpd	xmm1,xmm1,0x1		; xmm1 = b1i,b1r
	shufpd	xmm2,xmm2,0x1		; xmm2 = b2i,b2r
	unpcklpd xmm6,xmm6		; xmm6 = a00i,a00i
	unpcklpd xmm7,xmm7		; xmm7 = a11i,a11i
	xorpd	xmm0,[negate]		; xmm0 = -b0i,b0r	<_sse_sgn2>
	xorpd	xmm1,[negate]		; xmm1 = -b1i,b1r	<_sse_sgn2>
	xorpd	xmm2,[negate]		; xmm2 = -b2i,b2r	<_sse_sgn2>
	mulpd	xmm6,xmm0		; xmm6 = -a00i*b0i,a00i*b0r
	mulpd	xmm7,xmm1		; xmm7 = -a11i*b1i,a11i*b1r
	addpd 	xmm3,xmm6		; xmm3 = a00r*b0r + a01r*b1r + a02r*b2r - a00i*b0i , a00r*b0i + a01r*b1i + a02r*b2i + a00i*b0r
	addpd 	xmm4,xmm7		; xmm4 = a10r*b0r + a12r*b2r + a11r*b1r - a11i*b1i , a10r*b0i + a12r*b2i + a11r*b1i + a11i*b1r
	movsd	xmm6,[eax+136]		; xmm6 = x,a22i	<aa[3].e[2][2].imag>
	movsd	xmm7,[eax+56]		; xmm7 = x,a10i	<aa[3].e[1][0].imag>
	unpcklpd xmm6,xmm6		; xmm6 = a22i,a22i
	unpcklpd xmm7,xmm7		; xmm7 = a10i,a10i
	mulpd 	xmm6,xmm2		; xmm6 = -a22i*b2i,a22i*b2r
	mulpd	xmm7,xmm0		; xmm7 = -a10i*b0i,a10i*b0r
	addpd 	xmm5,xmm6		; xmm5 = a20r*b0r + a21r*b1r + a22r*b2r - a22i*b2i , a20r*b0i + a21r*b1i + a22r*b2i + a22i*b2r
	addpd 	xmm4,xmm7		; xmm4 = a10r*b0r + a12r*b2r + a11r*b1r - a11i*b1i - a10i*b0i, a10r*b0i + a12r*b2i + a11r*b1i + a11i*b1r + a10i*b0r
	movsd	xmm6,[eax+24]		; xmm6 = x,a01i	<aa[3].e[0][1].imag>
	movsd	xmm7,[eax+104]		; xmm7 = x,a20i	<aa[3].e[2][0].imag>
	unpcklpd xmm6,xmm6		; xmm6 = a01i,a01i
	unpcklpd xmm7,xmm7		; xmm7 = a20i,a20i
	mulpd 	xmm6,xmm1		; xmm6 = -a01i*b1i,a01i*b1r
	mulpd 	xmm7,xmm0		; xmm7 = -a20i*b0i,a20i*b0r
	addpd 	xmm3,xmm6		; xmm3 = a00r*b0r + a01r*b1r + a02r*b2r - a00i*b0i - a01i*b1i , a00r*b0i + a01r*b1i + a02r*b2i + a00i*b0r + a01i*b1r
	addpd 	xmm5,xmm7		; xmm5 = a20r*b0r + a21r*b1r + a22r*b2r - a22i*b2i - a20i*b0i , a20r*b0i + a21r*b1i + a22r*b2i + a22i*b2r + a20i*b0r
	movsd 	xmm0,[eax+40]		; xmm0 = x,a02i	<aa[3].e[0][2].imag>
	movsd	xmm6,[eax+120]		; xmm6 = x,a21i	<aa[3].e[2][1].imag>
	movsd 	xmm7,[eax+88]		; xmm7 = x,a12i	<aa[3].e[1][2].imag>
	unpcklpd xmm0,xmm0		; xmm0 = a02i,a02i
	unpcklpd xmm6,xmm6		; xmm6 = a21i,a21i
	unpcklpd xmm7,xmm7		; xmm7 = a12i,a12i
	mulpd 	xmm0,xmm2		; xmm0 = -a02i*b2i,a02i*b2r
	mulpd	xmm6,xmm1		; xmm6 = -a21i*b1i,a21i*b1r
	mulpd	xmm7,xmm2		; xmm7 = -a12i*b2i,a12i*b2r
	addpd 	xmm3,xmm0		; xmm3 = a00r*b0r + a01r*b1r + a02r*b2r - a00i*b0i - a01i*b1i - a02i*b2i , a00r*b0i + a01r*b1i + a02r*b2i + a00i*b0r + a01i*b1r + a02i*b2r
		                        ;        a00r*b0r - a00i*b0i + a01r*b1r - a01i*b1i + a02r*b2r - a02i*b2i , a00r*b0i + a00i*b0r + a01r*b1i + a01i*b1r + a02r*b2i + a02i*b2r 
	addpd 	xmm4,xmm7		; xmm4 = a10r*b0r + a12r*b2r + a11r*b1r - a11i*b1i - a10i*b0i - a12i*b2i , a10r*b0i + a12r*b2i + a11r*b1i + a11i*b1r + a10i*b0r + a12i*b2r
					;        a10r*b0r - a10i*b0i + a11r*b1r - a11i*b1i + a12r*b2r - a12i*b2i , a10r*b0i + a10i*b0r + a11r*b1i + a11i*b1r + a12r*b2i + a12i*b2r
	addpd 	xmm5,xmm6		; xmm5 = a20r*b0r + a21r*b1r + a22r*b2r - a22i*b2i - a20i*b0i - a21i*b1i , a20r*b0i + a21r*b1i + a22r*b2i + a22i*b2r + a20i*b0r + a21i*b1r
        	                        ;        a20r*b0r - a20i*b0i + a21r*b1r - a21i*b1i + a22r*b2r - a22i*b2i , a20r*b0i + a20i*b0r + a21r*b1i + a21i*b1r + a22r*b2i + a22i*b2r 
        movupd  xmm0,[ecx]		;	<(cc)->c[0]>
        movupd  xmm6,[ecx+16]		; 	<(cc)->c[1]>
        movupd  xmm7,[ecx+32]		;  	<(cc)->c[2]>
        addpd   xmm3,xmm0
        addpd   xmm4,xmm6
        addpd   xmm5,xmm7
        movupd  [ecx],xmm3              ;       <(cc)->c[0]>
        movupd  [ecx+16],xmm4           ;       <(cc)->c[1]>
        movupd  [ecx+32],xmm5           ;       <(cc)->c[2]>	

here:	pop	ecx
	pop	ebx
	pop	eax
	mov	esp, ebp
	pop	ebp
	ret

	align	16
negate: dd	0x00000000
	dd	0x80000000
	dd	0x00000000
	dd	0x00000000
