;
; mult_adj_su3_mat_vec( su3_matrix *a, su3_vector *b, su3_vector *c)
;
                                                                                                             
global mult_adj_su3_mat_vec
mult_adj_su3_mat_vec:
	push	ebp
	mov	ebp, esp
	push	eax
	push	ebx
	push	ecx
	mov	eax, [ebp+8]		; su3_matrix *a
	mov	ebx, [ebp+12]		; su3_vector *b
	mov	ecx, [ebp+16]		; su3_vector *c

	movupd	xmm0,[ebx]		; <(bb)->c[0]>
	movupd	xmm1,[ebx+16]		; <(bb)->c[1]>
	movupd	xmm2,[ebx+32]		; <(bb)->c[2]>

	movsd	xmm3,[eax]		; xmm3 = x,a00r	<(aa)->e[0][0].real>
	movsd	xmm6,[eax+48]		; xmm6 = x,a10r	<(aa)->e[1][0].real>
	movsd	xmm4,[eax+16]		; xmm4 = x,a01r	<(aa)->e[0][1].real>
	movsd	xmm7,[eax+112]		; xmm7 = x,a21r	<(aa)->e[2][1].real>
	movsd	xmm5,[eax+32]		; xmm5 = x,a02r	<(aa)->e[0][2].real>
	unpcklpd xmm3,xmm3		; xmm3 = a00r,a00r
	unpcklpd xmm6,xmm6		; xmm6 = a10r,a10r
	unpcklpd xmm4,xmm4		; xmm4 = a01r,a01r
	mulpd	xmm3,xmm0		; xmm3 = a00r*b0r,a00r*b0i
	unpcklpd xmm7,xmm7		; xmm7 = a21r,a21r
	mulpd 	xmm6,xmm1		; xmm6 = a10r*b1r,a10r*b1i
	unpcklpd xmm5,xmm5		; xmm5 = a02r,a02r
	mulpd 	xmm4,xmm0		; xmm4 = a01r*b0r,a01r*b0i
	addpd 	xmm3,xmm6		; xmm3 = a00r*b0r+a10r*b1r,a00r*b0i+a10r*b1i
	mulpd 	xmm7,xmm2		; xmm7 = a21r*b2r,a21r*b2i
	mulpd 	xmm5,xmm0		; xmm5 = a02r*b0r,a02r*b0i
	addpd 	xmm4,xmm7		; xmm4 = a01r*b0r+a21r*b2r,a01r*b0i+a21r*b2i
	movsd 	xmm6,[eax+80]		; xmm6 = x,a12r	<(aa)->e[1][2].real>
	movsd 	xmm7,[eax+96]		; xmm7 = x,a20r <(aa)->e[2][0].real>
	unpcklpd xmm6,xmm6		; xmm6 = a12r,a12r
	unpcklpd xmm7,xmm7		; xmm7 = a20r,a20r
	mulpd 	xmm6,xmm1		; xmm6 = a12r*b1r,a12r*b1i
	mulpd 	xmm7,xmm2		; xmm7 = a20r*b2r,a20r*b2i
	addpd 	xmm5,xmm6		; xmm5 = a02r*b0r+a12r*b1r,a02r*b0i+a12r*b1i
	addpd 	xmm3,xmm7		; xmm3 = a00r*b0r + a10r*b1r + a20r*b2r , a00r*b0i + a10r*b1i + a20r*b2i
	movsd 	xmm6,[eax+64]		; xmm6 = x,a11r	<(aa)->e[1][1].real>
	movsd 	xmm7,[eax+128]		; xmm7 = x,a22r	<(aa)->e[2][2].real>
	unpcklpd xmm6,xmm6		; xmm6 = a11r,a11r
	unpcklpd xmm7,xmm7		; xmm7 = a22r,a22r
	mulpd 	xmm6,xmm1   		; xmm6 = a11r*b1r,a11r*b1i
	mulpd 	xmm7,xmm2   		; xmm7 = a22r*b2r,a22r*b2i
	addpd 	xmm4,xmm6		; xmm4 = a01r*b0r + a21r*b2r + a11r*b1r , a01r*b0i + a21r*b2i + a11r*b1i
	addpd 	xmm5,xmm7		; xmm5 = a02r*b0r + a12r*b1r + a22r*b2r , a02r*b0i + a12r*b1i + a22r*b2i
	;	
	; IMAGINARY
	;
	movsd	xmm6,[eax+8]		; xmm6 = x,a00i	<(aa)->e[0][0].imag>
	movsd	xmm7,[eax+72]		; xmm7 = x,a11i	<(aa)->e[1][1].imag>
	shufpd	xmm0,xmm0,0x1		; xmm0 = b0i,b0r
	shufpd	xmm1,xmm1,0x1		; xmm1 = b1i,b1r
	shufpd	xmm2,xmm2,0x1		; xmm2 = b2i,b2r
	unpcklpd xmm6,xmm6		; xmm6 = a00i,a00i
	unpcklpd xmm7,xmm7		; xmm7 = a11i,a11i
	xorpd	xmm0,[negate]		; xmm0 = b0i,-b0r	<_sse_sgn4>
	xorpd	xmm1,[negate]		; xmm1 = b1i,-b1r	<_sse_sgn4>
	xorpd	xmm2,[negate]		; xmm2 = b2i,-b2r	<_sse_sgn4>
	mulpd	xmm6,xmm0		; xmm6 = a00i*b0i,-a00i*b0r
	mulpd	xmm7,xmm1		; xmm7 = a11i*b1i,-a11i*b1r
	addpd 	xmm3,xmm6		; xmm3 = a00r*b0r + a10r*b1r + a20r*b2r + a00i*b0i , a00r*b0i + a10r*b1i + a20r*b2i - a00i*b0r
	addpd 	xmm4,xmm7		; xmm4 = a01r*b0r + a21r*b2r + a11r*b1r + a11i*b1i , a01r*b0i + a21r*b2i + a11r*b1i - a11i*b1r
	movsd	xmm6,[eax+136]		; xmm6 = x,a22i	<(aa)->e[2][2].imag>
	movsd	xmm7,[eax+24]		; xmm7 = x,a01i	<(aa)->e[0][1].imag>
	unpcklpd xmm6,xmm6		; xmm6 = a22i,a22i
	unpcklpd xmm7,xmm7		; xmm7 = a01i,a01i
	mulpd 	xmm6,xmm2		; xmm6 = a22i*b2i,-a22i*b2r
	mulpd	xmm7,xmm0		; xmm7 = a01i*b0i,-a01i*b0r
	addpd 	xmm5,xmm6		; xmm5 = a02r*b0r + a12r*b1r + a22r*b2r + a22i*b2i , a02r*b0i + a12r*b1i + a22r*b2i - a22i*b2r
	addpd 	xmm4,xmm7		; xmm4 = a01r*b0r + a21r*b2r + a11r*b1r + a11i*b1i + a01i*b0i, a01r*b0i + a21r*b2i + a11r*b1i - a11i*b1r - a01i*b0r
	movsd	xmm6,[eax+56]		; xmm6 = x,a10i	<(aa)->e[1][0].imag>
	movsd	xmm7,[eax+40]		; xmm7 = x,a02i	<(aa)->e[0][2].imag>
	unpcklpd xmm6,xmm6		; xmm6 = a10i,a10i
	unpcklpd xmm7,xmm7		; xmm7 = a02i,a02i
	mulpd 	xmm6,xmm1		; xmm6 = a10i*b1i,-a10i*b1r
	mulpd 	xmm7,xmm0		; xmm7 = a02i*b0i,-a02i*b0r
	addpd 	xmm3,xmm6		; xmm3 = a00r*b0r + a10r*b1r + a20r*b2r + a00i*b0i + a10i*b1i , a00r*b0i + a10r*b1i + a20r*b2i - a00i*b0r - a10i*b1r
	addpd 	xmm5,xmm7		; xmm5 = a02r*b0r + a12r*b1r + a22r*b2r + a22i*b2i + a02i*b0i , a02r*b0i + a12r*b1i + a22r*b2i - a22i*b2r - a02i*b0r
	movsd 	xmm0,[eax+104]		; xmm0 = x,a20i	<(aa)->e[2][0].imag>
	movsd	xmm6,[eax+88]		; xmm6 = x,a12i	<(aa)->e[1][2].imag>
	movsd 	xmm7,[eax+120]		; xmm7 = x,a21i	<(aa)->e[2][1].imag>
	unpcklpd xmm0,xmm0		; xmm0 = a20i,a20i
	unpcklpd xmm6,xmm6		; xmm6 = a12i,a12i
	unpcklpd xmm7,xmm7		; xmm7 = a21i,a21i
	mulpd 	xmm0,xmm2		; xmm0 = a20i*b2i,-a20i*b2r
	mulpd	xmm6,xmm1		; xmm6 = a12i*b1i,-a12i*b1r
	mulpd	xmm7,xmm2		; xmm7 = a21i*b2i,-a21i*b2r
	addpd 	xmm3,xmm0		; xmm3 = a00r*b0r + a10r*b1r + a20r*b2r + a00i*b0i + a10i*b1i + a20i*b2i , a00r*b0i + a10r*b1i + a20r*b2i - a00i*b0r - a10i*b1r - a20i*b2r
		                        ;        a00r*b0r + a00i*b0i + a10r*b1r + a10i*b1i + a20r*b2r + a20i*b2i , a00r*b0i - a00i*b0r + a10r*b1i - a10i*b1r + a20r*b2i - a20i*b2r 
	addpd 	xmm4,xmm7		; xmm4 = a01r*b0r + a21r*b2r + a11r*b1r + a11i*b1i + a01i*b0i + a21i*b2i , a01r*b0i + a21r*b2i + a11r*b1i - a11i*b1r - a01i*b0r - a21i*b2r
					;        a01r*b0r - a01i*b0i + a11r*b1r + a11i*b1i + a21r*b2r + a21i*b2i , a01r*b0i - a01i*b0r + a11r*b1i - a11i*b1r + a21r*b2i - a21i*b2r
	addpd 	xmm5,xmm6		; xmm5 = a02r*b0r + a12r*b1r + a22r*b2r + a22i*b2i + a02i*b0i + a12i*b1i , a02r*b0i + a12r*b1i + a22r*b2i - a22i*b2r - a02i*b0r - a12i*b1r
        	                        ;        a02r*b0r + a02i*b0i + a12r*b1r + a12i*b1i + a22r*b2r + a22i*b2i , a02r*b0i - a02i*b0r + a12r*b1i - a12i*b1r + a22r*b2i - a22i*b2r 
	movupd	[ecx],xmm3		;	<(cc)->c[0]>
	movupd	[ecx+16],xmm4		;	<(cc)->c[1]>
	movupd	[ecx+32],xmm5		;	<(cc)->c[2]>

	;
here:	pop	ecx
	pop	ebx
	pop	eax
	mov	esp, ebp
	pop	ebp
	ret

	align	16
negate: dd	0x00000000
	dd	0x00000000
	dd	0x00000000
	dd	0x80000000
