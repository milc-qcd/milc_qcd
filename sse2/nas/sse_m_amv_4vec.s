;
; mult_adj_su3_mat_4vec( su3_matrix *a, su3_vector *b, su3_vector *c0,
;                                su3_vector *c1, su3_vector *c2, su3_vector *c3)
;
; Multiply the adjoint of each of four input matrices by an input vector,
; storing the resulting vectors in 4 separate destinations.
;
 
global mult_adj_su3_mat_4vec
mult_adj_su3_mat_4vec:
	push	ebp
	mov	ebp, esp
	push	eax
	push	ebx
	push	ecx
	mov	eax, [ebp+8]		; su3_matrix *a
	mov	ebx, [ebp+12]		; su3_vector *b
	mov	ecx, [ebp+16]		; su3_vector *c
	;
	; first
	;
	movupd	xmm0,[ebx]		; <(bb)->c[0]>
	movupd	xmm1,[ebx+16]		; <(bb)->c[1]>
	movupd	xmm2,[ebx+32]		; <(bb)->c[2]>

	movsd	xmm3,[eax]		; xmm3 = x,a00r	<aa[0].e[0][0].real>
	movsd	xmm6,[eax+48]		; xmm6 = x,a10r	<aa[0].e[1][0].real>
	movsd	xmm4,[eax+16]		; xmm4 = x,a01r	<aa[0].e[0][1].real>
	movsd	xmm7,[eax+112]		; xmm7 = x,a21r	<aa[0].e[2][1].real>
	movsd	xmm5,[eax+32]		; xmm5 = x,a02r	<aa[0].e[0][2].real>
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
	movsd 	xmm6,[eax+80]		; xmm6 = x,a12r	<aa[0].e[1][2].real>
	movsd 	xmm7,[eax+96]		; xmm7 = x,a20r <aa[0].e[2][0].real>
	unpcklpd xmm6,xmm6		; xmm6 = a12r,a12r
	unpcklpd xmm7,xmm7		; xmm7 = a20r,a20r
	mulpd 	xmm6,xmm1		; xmm6 = a12r*b1r,a12r*b1i
	mulpd 	xmm7,xmm2		; xmm7 = a20r*b2r,a20r*b2i
	addpd 	xmm5,xmm6		; xmm5 = a02r*b0r+a12r*b1r,a02r*b0i+a12r*b1i
	addpd 	xmm3,xmm7		; xmm3 = a00r*b0r + a10r*b1r + a20r*b2r , a00r*b0i + a10r*b1i + a20r*b2i
	movsd 	xmm6,[eax+64]		; xmm6 = x,a11r	<aa[0].e[1][1].real>
	movsd 	xmm7,[eax+128]		; xmm7 = x,a22r	<aa[0].e[2][2].real>
	unpcklpd xmm6,xmm6		; xmm6 = a11r,a11r
	unpcklpd xmm7,xmm7		; xmm7 = a22r,a22r
	mulpd 	xmm6,xmm1   		; xmm6 = a11r*b1r,a11r*b1i
	mulpd 	xmm7,xmm2   		; xmm7 = a22r*b2r,a22r*b2i
	addpd 	xmm4,xmm6		; xmm4 = a01r*b0r + a21r*b2r + a11r*b1r , a01r*b0i + a21r*b2i + a11r*b1i
	addpd 	xmm5,xmm7		; xmm5 = a02r*b0r + a12r*b1r + a22r*b2r , a02r*b0i + a12r*b1i + a22r*b2i
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
	xorpd	xmm0,[negate]		; xmm0 = b0i,-b0r	<_sse_sgn4>
	xorpd	xmm1,[negate]		; xmm1 = b1i,-b1r	<_sse_sgn4>
	xorpd	xmm2,[negate]		; xmm2 = b2i,-b2r	<_sse_sgn4>
	mulpd	xmm6,xmm0		; xmm6 = a00i*b0i,-a00i*b0r
	mulpd	xmm7,xmm1		; xmm7 = a11i*b1i,-a11i*b1r
	addpd 	xmm3,xmm6		; xmm3 = a00r*b0r + a10r*b1r + a20r*b2r + a00i*b0i , a00r*b0i + a10r*b1i + a20r*b2i - a00i*b0r
	addpd 	xmm4,xmm7		; xmm4 = a01r*b0r + a21r*b2r + a11r*b1r + a11i*b1i , a01r*b0i + a21r*b2i + a11r*b1i - a11i*b1r
	movsd	xmm6,[eax+136]		; xmm6 = x,a22i	<aa[0].e[2][2].imag>
	movsd	xmm7,[eax+24]		; xmm7 = x,a01i	<aa[0].e[0][1].imag>
	unpcklpd xmm6,xmm6		; xmm6 = a22i,a22i
	unpcklpd xmm7,xmm7		; xmm7 = a01i,a01i
	mulpd 	xmm6,xmm2		; xmm6 = a22i*b2i,-a22i*b2r
	mulpd	xmm7,xmm0		; xmm7 = a01i*b0i,-a01i*b0r
	addpd 	xmm5,xmm6		; xmm5 = a02r*b0r + a12r*b1r + a22r*b2r + a22i*b2i , a02r*b0i + a12r*b1i + a22r*b2i - a22i*b2r
	addpd 	xmm4,xmm7		; xmm4 = a01r*b0r + a21r*b2r + a11r*b1r + a11i*b1i + a01i*b0i, a01r*b0i + a21r*b2i + a11r*b1i - a11i*b1r - a01i*b0r
	movsd	xmm6,[eax+56]		; xmm6 = x,a10i	<aa[0].e[1][0].imag>
	movsd	xmm7,[eax+40]		; xmm7 = x,a02i	<aa[0].e[0][2].imag>
	unpcklpd xmm6,xmm6		; xmm6 = a10i,a10i
	unpcklpd xmm7,xmm7		; xmm7 = a02i,a02i
	mulpd 	xmm6,xmm1		; xmm6 = a10i*b1i,-a10i*b1r
	mulpd 	xmm7,xmm0		; xmm7 = a02i*b0i,-a02i*b0r
	addpd 	xmm3,xmm6		; xmm3 = a00r*b0r + a10r*b1r + a20r*b2r + a00i*b0i + a10i*b1i , a00r*b0i + a10r*b1i + a20r*b2i - a00i*b0r - a10i*b1r
	addpd 	xmm5,xmm7		; xmm5 = a02r*b0r + a12r*b1r + a22r*b2r + a22i*b2i + a02i*b0i , a02r*b0i + a12r*b1i + a22r*b2i - a22i*b2r - a02i*b0r
	movsd 	xmm0,[eax+104]		; xmm0 = x,a20i	<aa[0].e[2][0].imag>
	movsd	xmm6,[eax+88]		; xmm6 = x,a12i	<aa[0].e[1][2].imag>
	movsd 	xmm7,[eax+120]		; xmm7 = x,a21i	<aa[0].e[2][1].imag>
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
	movupd	[ecx],xmm3		;	<(cc0)->c[0]>
	movupd	[ecx+16],xmm4		;	<(cc0)->c[1]>
	movupd	[ecx+32],xmm5		;	<(cc0)->c[2]>
	;
	; restore original su3_vector b values
	;
	movupd  xmm0,[ebx]		; <(bb)->c[0]>
	shufpd  xmm1,xmm1,0x1
	shufpd  xmm2,xmm2,0x1
	xorpd	xmm1,[negate2]		; <_sse_sgn2>
	xorpd	xmm2,[negate2]		; <_sse_sgn2>
	;
	; second
	;
	mov	ecx, [ebp+20]		; su3_vector *c
	;
	movsd	xmm3,[eax+144]		; xmm3 = x,a00r	<aa[1].e[0][0].real>
	movsd	xmm6,[eax+192]		; xmm6 = x,a10r	<aa[1].e[1][0].real>
	movsd	xmm4,[eax+160]		; xmm4 = x,a01r	<aa[1].e[0][1].real>
	movsd	xmm7,[eax+256]		; xmm7 = x,a21r	<aa[1].e[2][1].real>
	movsd	xmm5,[eax+176]		; xmm5 = x,a02r	<aa[1].e[0][2].real>
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
	movsd 	xmm6,[eax+224]		; xmm6 = x,a12r	<aa[1].e[1][2].real>
	movsd 	xmm7,[eax+240]		; xmm7 = x,a20r <aa[1].e[2][0].real>
	unpcklpd xmm6,xmm6		; xmm6 = a12r,a12r
	unpcklpd xmm7,xmm7		; xmm7 = a20r,a20r
	mulpd 	xmm6,xmm1		; xmm6 = a12r*b1r,a12r*b1i
	mulpd 	xmm7,xmm2		; xmm7 = a20r*b2r,a20r*b2i
	addpd 	xmm5,xmm6		; xmm5 = a02r*b0r+a12r*b1r,a02r*b0i+a12r*b1i
	addpd 	xmm3,xmm7		; xmm3 = a00r*b0r + a10r*b1r + a20r*b2r , a00r*b0i + a10r*b1i + a20r*b2i
	movsd 	xmm6,[eax+208]		; xmm6 = x,a11r	<aa[1].e[1][1].real>
	movsd 	xmm7,[eax+272]		; xmm7 = x,a22r	<aa[1].e[2][2].real>
	unpcklpd xmm6,xmm6		; xmm6 = a11r,a11r
	unpcklpd xmm7,xmm7		; xmm7 = a22r,a22r
	mulpd 	xmm6,xmm1   		; xmm6 = a11r*b1r,a11r*b1i
	mulpd 	xmm7,xmm2   		; xmm7 = a22r*b2r,a22r*b2i
	addpd 	xmm4,xmm6		; xmm4 = a01r*b0r + a21r*b2r + a11r*b1r , a01r*b0i + a21r*b2i + a11r*b1i
	addpd 	xmm5,xmm7		; xmm5 = a02r*b0r + a12r*b1r + a22r*b2r , a02r*b0i + a12r*b1i + a22r*b2i
	;	
	; IMAGINARY
	;
	movsd	xmm6,[eax+152]		; xmm6 = x,a00i	<aa[1].e[0][0].imag>
	movsd	xmm7,[eax+216]		; xmm7 = x,a11i	<aa[1].e[1][1].imag>
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
	movsd	xmm6,[eax+280]		; xmm6 = x,a22i	<aa[1].e[2][2].imag>
	movsd	xmm7,[eax+168]		; xmm7 = x,a01i	<aa[1].e[0][1].imag>
	unpcklpd xmm6,xmm6		; xmm6 = a22i,a22i
	unpcklpd xmm7,xmm7		; xmm7 = a01i,a01i
	mulpd 	xmm6,xmm2		; xmm6 = a22i*b2i,-a22i*b2r
	mulpd	xmm7,xmm0		; xmm7 = a01i*b0i,-a01i*b0r
	addpd 	xmm5,xmm6		; xmm5 = a02r*b0r + a12r*b1r + a22r*b2r + a22i*b2i , a02r*b0i + a12r*b1i + a22r*b2i - a22i*b2r
	addpd 	xmm4,xmm7		; xmm4 = a01r*b0r + a21r*b2r + a11r*b1r + a11i*b1i + a01i*b0i, a01r*b0i + a21r*b2i + a11r*b1i - a11i*b1r - a01i*b0r
	movsd	xmm6,[eax+200]		; xmm6 = x,a10i	<aa[1].e[1][0].imag>
	movsd	xmm7,[eax+184]		; xmm7 = x,a02i	<aa[1].e[0][2].imag>
	unpcklpd xmm6,xmm6		; xmm6 = a10i,a10i
	unpcklpd xmm7,xmm7		; xmm7 = a02i,a02i
	mulpd 	xmm6,xmm1		; xmm6 = a10i*b1i,-a10i*b1r
	mulpd 	xmm7,xmm0		; xmm7 = a02i*b0i,-a02i*b0r
	addpd 	xmm3,xmm6		; xmm3 = a00r*b0r + a10r*b1r + a20r*b2r + a00i*b0i + a10i*b1i , a00r*b0i + a10r*b1i + a20r*b2i - a00i*b0r - a10i*b1r
	addpd 	xmm5,xmm7		; xmm5 = a02r*b0r + a12r*b1r + a22r*b2r + a22i*b2i + a02i*b0i , a02r*b0i + a12r*b1i + a22r*b2i - a22i*b2r - a02i*b0r
	movsd 	xmm0,[eax+248]		; xmm0 = x,a20i	<aa[1].e[2][0].imag>
	movsd	xmm6,[eax+232]		; xmm6 = x,a12i	<aa[1].e[1][2].imag>
	movsd 	xmm7,[eax+264]		; xmm7 = x,a21i	<aa[1].e[2][1].imag>
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
	movupd	[ecx],xmm3		;	<(cc1)->c[0]>
	movupd	[ecx+16],xmm4		;	<(cc1)->c[1]>
	movupd	[ecx+32],xmm5		;	<(cc1)->c[2]>
	;
	; restore original su3_vector b values
	;
	movupd  xmm0,[ebx]		; <(bb)->c[0]>
	shufpd  xmm1,xmm1,0x1
	shufpd  xmm2,xmm2,0x1
	xorpd	xmm1,[negate2]		; <_sse_sgn2>
	xorpd	xmm2,[negate2]          ; <_sse_sgn2>
	;
	; third
	;
	mov	ecx, [ebp+24]
	;
	movsd	xmm3,[eax+288]		; xmm3 = x,a00r	<aa[2].e[0][0].real>
	movsd	xmm6,[eax+336]		; xmm6 = x,a10r	<aa[2].e[1][0].real>
	movsd	xmm4,[eax+304]		; xmm4 = x,a01r	<aa[2].e[0][1].real>
	movsd	xmm7,[eax+400]		; xmm7 = x,a21r	<aa[2].e[2][1].real>
	movsd	xmm5,[eax+320]		; xmm5 = x,a02r	<aa[2].e[0][2].real>
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
	movsd 	xmm6,[eax+368]		; xmm6 = x,a12r	<aa[2].e[1][2].real>
	movsd 	xmm7,[eax+384]		; xmm7 = x,a20r <aa[2].e[2][0].real>
	unpcklpd xmm6,xmm6		; xmm6 = a12r,a12r
	unpcklpd xmm7,xmm7		; xmm7 = a20r,a20r
	mulpd 	xmm6,xmm1		; xmm6 = a12r*b1r,a12r*b1i
	mulpd 	xmm7,xmm2		; xmm7 = a20r*b2r,a20r*b2i
	addpd 	xmm5,xmm6		; xmm5 = a02r*b0r+a12r*b1r,a02r*b0i+a12r*b1i
	addpd 	xmm3,xmm7		; xmm3 = a00r*b0r + a10r*b1r + a20r*b2r , a00r*b0i + a10r*b1i + a20r*b2i
	movsd 	xmm6,[eax+352]		; xmm6 = x,a11r	<aa[2].e[1][1].real>
	movsd 	xmm7,[eax+416]		; xmm7 = x,a22r	<aa[2].e[2][2].real>
	unpcklpd xmm6,xmm6		; xmm6 = a11r,a11r
	unpcklpd xmm7,xmm7		; xmm7 = a22r,a22r
	mulpd 	xmm6,xmm1   		; xmm6 = a11r*b1r,a11r*b1i
	mulpd 	xmm7,xmm2   		; xmm7 = a22r*b2r,a22r*b2i
	addpd 	xmm4,xmm6		; xmm4 = a01r*b0r + a21r*b2r + a11r*b1r , a01r*b0i + a21r*b2i + a11r*b1i
	addpd 	xmm5,xmm7		; xmm5 = a02r*b0r + a12r*b1r + a22r*b2r , a02r*b0i + a12r*b1i + a22r*b2i
	;	
	; IMAGINARY
	;
	movsd	xmm6,[eax+296]		; xmm6 = x,a00i	<aa[2].e[0][0].imag>
	movsd	xmm7,[eax+360]		; xmm7 = x,a11i	<aa[2].e[1][1].imag>
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
	movsd	xmm6,[eax+424]		; xmm6 = x,a22i	<aa[2].e[2][2].imag>
	movsd	xmm7,[eax+312]		; xmm7 = x,a01i	<aa[2].e[0][1].imag>
	unpcklpd xmm6,xmm6		; xmm6 = a22i,a22i
	unpcklpd xmm7,xmm7		; xmm7 = a01i,a01i
	mulpd 	xmm6,xmm2		; xmm6 = a22i*b2i,-a22i*b2r
	mulpd	xmm7,xmm0		; xmm7 = a01i*b0i,-a01i*b0r
	addpd 	xmm5,xmm6		; xmm5 = a02r*b0r + a12r*b1r + a22r*b2r + a22i*b2i , a02r*b0i + a12r*b1i + a22r*b2i - a22i*b2r
	addpd 	xmm4,xmm7		; xmm4 = a01r*b0r + a21r*b2r + a11r*b1r + a11i*b1i + a01i*b0i, a01r*b0i + a21r*b2i + a11r*b1i - a11i*b1r - a01i*b0r
	movsd	xmm6,[eax+344]		; xmm6 = x,a10i	<aa[2].e[1][0].imag>
	movsd	xmm7,[eax+328]		; xmm7 = x,a02i	<aa[2].e[0][2].imag>
	unpcklpd xmm6,xmm6		; xmm6 = a10i,a10i
	unpcklpd xmm7,xmm7		; xmm7 = a02i,a02i
	mulpd 	xmm6,xmm1		; xmm6 = a10i*b1i,-a10i*b1r
	mulpd 	xmm7,xmm0		; xmm7 = a02i*b0i,-a02i*b0r
	addpd 	xmm3,xmm6		; xmm3 = a00r*b0r + a10r*b1r + a20r*b2r + a00i*b0i + a10i*b1i , a00r*b0i + a10r*b1i + a20r*b2i - a00i*b0r - a10i*b1r
	addpd 	xmm5,xmm7		; xmm5 = a02r*b0r + a12r*b1r + a22r*b2r + a22i*b2i + a02i*b0i , a02r*b0i + a12r*b1i + a22r*b2i - a22i*b2r - a02i*b0r
	movsd 	xmm0,[eax+392]		; xmm0 = x,a20i	<aa[2].e[2][0].imag>
	movsd	xmm6,[eax+376]		; xmm6 = x,a12i	<aa[2].e[1][2].imag>
	movsd 	xmm7,[eax+408]		; xmm7 = x,a21i	<aa[2].e[2][1].imag>
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
	movupd	[ecx],xmm3		;	<(cc2)->c[0]>
	movupd	[ecx+16],xmm4		;	<(cc2)->c[1]>
	movupd	[ecx+32],xmm5		;	<(cc2)->c[2]>
	;
	; restore original su3_vector b values
	;
	movupd  xmm0,[ebx]		; <(bb)->c[0]>
	shufpd  xmm1,xmm1,0x1
	shufpd  xmm2,xmm2,0x1
	xorpd	xmm1,[negate2]          ; <_sse_sgn2>
	xorpd	xmm2,[negate2]          ; <_sse_sgn2>
	;
	; fourth
	;
	mov	ecx, [ebp+28]
	;
	movsd	xmm3,[eax+432]		; xmm3 = x,a00r	<aa[3].e[0][0].real>
	movsd	xmm6,[eax+480]		; xmm6 = x,a10r	<aa[3].e[1][0].real>
	movsd	xmm4,[eax+448]		; xmm4 = x,a01r	<aa[3].e[0][1].real>
	movsd	xmm7,[eax+544]		; xmm7 = x,a21r	<aa[3].e[2][1].real>
	movsd	xmm5,[eax+464]		; xmm5 = x,a02r	<aa[3].e[0][2].real>
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
	movsd 	xmm6,[eax+512]		; xmm6 = x,a12r	<aa[3].e[1][2].real>
	movsd 	xmm7,[eax+528]		; xmm7 = x,a20r <aa[3].e[2][0].real>
	unpcklpd xmm6,xmm6		; xmm6 = a12r,a12r
	unpcklpd xmm7,xmm7		; xmm7 = a20r,a20r
	mulpd 	xmm6,xmm1		; xmm6 = a12r*b1r,a12r*b1i
	mulpd 	xmm7,xmm2		; xmm7 = a20r*b2r,a20r*b2i
	addpd 	xmm5,xmm6		; xmm5 = a02r*b0r+a12r*b1r,a02r*b0i+a12r*b1i
	addpd 	xmm3,xmm7		; xmm3 = a00r*b0r + a10r*b1r + a20r*b2r , a00r*b0i + a10r*b1i + a20r*b2i
	movsd 	xmm6,[eax+496]		; xmm6 = x,a11r	<aa[3].e[1][1].real>
	movsd 	xmm7,[eax+560]		; xmm7 = x,a22r	<aa[3].e[2][2].real>
	unpcklpd xmm6,xmm6		; xmm6 = a11r,a11r
	unpcklpd xmm7,xmm7		; xmm7 = a22r,a22r
	mulpd 	xmm6,xmm1   		; xmm6 = a11r*b1r,a11r*b1i
	mulpd 	xmm7,xmm2   		; xmm7 = a22r*b2r,a22r*b2i
	addpd 	xmm4,xmm6		; xmm4 = a01r*b0r + a21r*b2r + a11r*b1r , a01r*b0i + a21r*b2i + a11r*b1i
	addpd 	xmm5,xmm7		; xmm5 = a02r*b0r + a12r*b1r + a22r*b2r , a02r*b0i + a12r*b1i + a22r*b2i
	;	
	; IMAGINARY
	;
	movsd	xmm6,[eax+440]		; xmm6 = x,a00i	<aa[3].e[0][0].imag>
	movsd	xmm7,[eax+504]		; xmm7 = x,a11i	<aa[3].e[1][1].imag>
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
	movsd	xmm6,[eax+568]		; xmm6 = x,a22i	<aa[3].e[2][2].imag>
	movsd	xmm7,[eax+456]		; xmm7 = x,a01i	<aa[3].e[0][1].imag>
	unpcklpd xmm6,xmm6		; xmm6 = a22i,a22i
	unpcklpd xmm7,xmm7		; xmm7 = a01i,a01i
	mulpd 	xmm6,xmm2		; xmm6 = a22i*b2i,-a22i*b2r
	mulpd	xmm7,xmm0		; xmm7 = a01i*b0i,-a01i*b0r
	addpd 	xmm5,xmm6		; xmm5 = a02r*b0r + a12r*b1r + a22r*b2r + a22i*b2i , a02r*b0i + a12r*b1i + a22r*b2i - a22i*b2r
	addpd 	xmm4,xmm7		; xmm4 = a01r*b0r + a21r*b2r + a11r*b1r + a11i*b1i + a01i*b0i, a01r*b0i + a21r*b2i + a11r*b1i - a11i*b1r - a01i*b0r
	movsd	xmm6,[eax+488]		; xmm6 = x,a10i	<aa[3].e[1][0].imag>
	movsd	xmm7,[eax+472]		; xmm7 = x,a02i	<aa[3].e[0][2].imag>
	unpcklpd xmm6,xmm6		; xmm6 = a10i,a10i
	unpcklpd xmm7,xmm7		; xmm7 = a02i,a02i
	mulpd 	xmm6,xmm1		; xmm6 = a10i*b1i,-a10i*b1r
	mulpd 	xmm7,xmm0		; xmm7 = a02i*b0i,-a02i*b0r
	addpd 	xmm3,xmm6		; xmm3 = a00r*b0r + a10r*b1r + a20r*b2r + a00i*b0i + a10i*b1i , a00r*b0i + a10r*b1i + a20r*b2i - a00i*b0r - a10i*b1r
	addpd 	xmm5,xmm7		; xmm5 = a02r*b0r + a12r*b1r + a22r*b2r + a22i*b2i + a02i*b0i , a02r*b0i + a12r*b1i + a22r*b2i - a22i*b2r - a02i*b0r
	movsd 	xmm0,[eax+536]		; xmm0 = x,a20i	<aa[3].e[2][0].imag>
	movsd	xmm6,[eax+520]		; xmm6 = x,a12i	<aa[3].e[1][2].imag>
	movsd 	xmm7,[eax+552]		; xmm7 = x,a21i	<aa[3].e[2][1].imag>
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
	movupd	[ecx],xmm3		;	<(cc3)->c[0]>
	movupd	[ecx+16],xmm4		;	<(cc3)->c[1]>
	movupd	[ecx+32],xmm5		;	<(cc3)->c[2]>
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

negate2: dd 	0x00000000
	 dd     0x80000000
	 dd	0x00000000
	 dd	0x00000000
