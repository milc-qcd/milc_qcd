;
; mult_su3_na( su3_matrix *a, su3_matrix *b, su3_matrix *c)
;
; | a00 a01 a02 | | b00 b01 b02 | | c00 c01 c02 |
; | a10 a11 a12 |*| b10 b11 b12 |=| c10 c11 c12 |
; | a20 a21 a22 | | b20 b21 b22 | | c20 c21 c22 |

global mult_su3_na
mult_su3_na:
	push	ebp
	mov	ebp, esp
	push	eax
	push	ebx
	push	ecx
	mov	eax, [ebp+8]		; su3_matrix *a
	mov	ebx, [ebp+12]		; su3_matrix *b
	mov	ecx, [ebp+16]		; su3_matrix *c

	movupd	xmm0,[ebx]		; xmm0 = b00r,b00i <(bb)->e[0][0]>
	movupd	xmm1,[ebx+16]		; xmm1 = b01r,b01i <(bb)->e[0][1]>
	movupd	xmm2,[ebx+32]		; xmm2 = b02r,b02i <(bb)->e[0][2]>
	xorpd	xmm0,[negate]		; xmm0 = b00r,-b00i	<_sse_sgn4>
	xorpd	xmm1,[negate]		; xmm1 = b10r,-b10i	<_sse_sgn4>
	xorpd	xmm2,[negate]		; xmm2 = b20r,-b20i	<_sse_sgn4>
	;
	; REAL	Part I
	;
	movsd	xmm3,[eax]		; xmm3 = x,a00r	<(aa)->e[0][0].real>
	movsd	xmm6,[eax+16]		; xmm6 = x,a01r	<(aa)->e[0][1].real>
	movsd	xmm4,[eax+48]		; xmm4 = x,a10r	<(aa)->e[1][0].real>
	movsd	xmm7,[eax+80]		; xmm7 = x,a12r	<(aa)->e[1][2].real>
	movsd	xmm5,[eax+96]		; xmm5 = x,a20r	<(aa)->e[2][0].real>
	unpcklpd xmm3,xmm3		; xmm3 = a00r,a00r
	unpcklpd xmm6,xmm6		; xmm6 = a01r,a01r
	unpcklpd xmm4,xmm4		; xmm4 = a10r,a10r
	mulpd	xmm3,xmm0		; xmm3 = a00r*b00r,-a00r*b00i
	unpcklpd xmm7,xmm7		; xmm7 = a12r,a12r
	mulpd 	xmm6,xmm1		; xmm6 = a01r*b10r,-a01r*b10i
	unpcklpd xmm5,xmm5		; xmm5 = a20r,a20r
	mulpd 	xmm4,xmm0		; xmm4 = a10r*b00r,-a10r*b00i
	addpd 	xmm3,xmm6		; xmm3 = a00r*b00r+a01r*b10r,-a00r*b00i-a01r*b10i
	mulpd 	xmm7,xmm2		; xmm7 = a12r*b20r,-a12r*b20i
	mulpd 	xmm5,xmm0		; xmm5 = a20r*b00r,-a20r*b00i
	addpd 	xmm4,xmm7		; xmm4 = a10r*b00r+a12r*b20r,-a10r*b00i-a12r*b20i
	movsd 	xmm6,[eax+112]		; xmm6 = x,a21r	<(aa)->e[2][1].real>
	movsd 	xmm7,[eax+32]		; xmm7 = x,a02r <(aa)->e[0][2].real>
	unpcklpd xmm6,xmm6		; xmm6 = a21r,a21r
	unpcklpd xmm7,xmm7		; xmm7 = a02r,a02r
	mulpd 	xmm6,xmm1		; xmm6 = a21r*b10r,-a21r*b10i
	mulpd 	xmm7,xmm2		; xmm7 = a02r*b20r,-a02r*b20i
	addpd 	xmm5,xmm6		; xmm5 = a20r*b00r+a21r*b10r,a20r*b00i+a21r*b10i
	addpd 	xmm3,xmm7		; xmm3 = a00r*b00r + a01r*b10r + a02r*b20r , a00r*b00i + a01r*b10i + a02r*b20i
	movsd 	xmm6,[eax+64]		; xmm6 = x,a11r	<(aa)->e[1][1].real>
	movsd 	xmm7,[eax+128]		; xmm7 = x,a22r	<(aa)->e[2][2].real>
	unpcklpd xmm6,xmm6		; xmm6 = a11r,a11r
	unpcklpd xmm7,xmm7		; xmm7 = a22r,a22r
	mulpd 	xmm6,xmm1   		; xmm6 = a11r*b10r,a11r*b10i
	mulpd 	xmm7,xmm2   		; xmm7 = a22r*b20r,a22r*b20i
	addpd 	xmm4,xmm6		; xmm4 = a10r*b00r + a12r*b20r + a11r*b10r , a10r*b00i + a12r*b20i + a11r*b10i
	addpd 	xmm5,xmm7		; xmm5 = a20r*b00r + a21r*b10r + a22r*b20r , a20r*b00i + a21r*b10i + a22r*b20i
	;	
	; IMAGINARY	Part I
	;
	movsd	xmm6,[eax+8]		; xmm6 = x,a00i	<(aa)->e[0][0].imag>
	movsd	xmm7,[eax+72]		; xmm7 = x,a11i	<(aa)->e[1][1].imag>

	xorpd	xmm0,[negate]		; xmm0 = b00r,b00i	<_sse_sgn4>
	xorpd	xmm1,[negate]		; xmm1 = b10r,b10i	<_sse_sgn4>
	xorpd	xmm2,[negate]		; xmm2 = b20r,b20i	<_sse_sgn4>

	shufpd	xmm0,xmm0,0x1		; xmm0 = b00i,b00r
	shufpd	xmm1,xmm1,0x1		; xmm1 = b10i,b10r
	shufpd	xmm2,xmm2,0x1		; xmm2 = b20i,b20r
	unpcklpd xmm6,xmm6		; xmm6 = a00i,a00i
	unpcklpd xmm7,xmm7		; xmm7 = a11i,a11i
	mulpd	xmm6,xmm0		; xmm6 = -a00i*b00i,a00i*b00r
	mulpd	xmm7,xmm1		; xmm7 = -a11i*b10i,a11i*b10r
	addpd 	xmm3,xmm6		; xmm3 = a00r*b00r + a01r*b10r + a02r*b20r - a00i*b00i , a00r*b00i + a01r*b10i + a02r*b20i + a00i*b00r
	addpd 	xmm4,xmm7		; xmm4 = a10r*b00r + a12r*b20r + a11r*b10r - a11i*b10i , a10r*b00i + a12r*b20i + a11r*b10i + a11i*b10r
	movsd	xmm6,[eax+136]		; xmm6 = x,a22i	<(aa)->e[2][2].imag>
	movsd	xmm7,[eax+56]		; xmm7 = x,a10i	<(aa)->e[1][0].imag>
	unpcklpd xmm6,xmm6		; xmm6 = a22i,a22i
	unpcklpd xmm7,xmm7		; xmm7 = a10i,a10i
	mulpd 	xmm6,xmm2		; xmm6 = -a22i*b20i,a22i*b20r
	mulpd	xmm7,xmm0		; xmm7 = -a10i*b00i,a10i*b00r
	addpd 	xmm5,xmm6		; xmm5 = a20r*b00r + a21r*b10r + a22r*b20r - a22i*b20i , a20r*b00i + a21r*b10i + a22r*b20i + a22i*b20r
	addpd 	xmm4,xmm7		; xmm4 = a10r*b00r + a12r*b20r + a11r*b10r - a11i*b10i - a10i*b00i0, a10r*b00i + a12r*b20i + a11r*b10i + a11i*b10r + a10i*b00r
	movsd	xmm6,[eax+24]		; xmm6 = x,a01i	<(aa)->e[0][1].imag>
	movsd	xmm7,[eax+104]		; xmm7 = x,a20i	<(aa)->e[2][0].imag>
	unpcklpd xmm6,xmm6		; xmm6 = a01i,a01i
	unpcklpd xmm7,xmm7		; xmm7 = a20i,a20i
	mulpd 	xmm6,xmm1		; xmm6 = -a01i*b10i,a01i*b10r
	mulpd 	xmm7,xmm0		; xmm7 = -a20i*b00i,a20i*b00r
	addpd 	xmm3,xmm6		; xmm3 = a00r*b00r + a01r*b10r + a02r*b20r - a00i*b00i - a01i*b10i , a00r*b00i + a01r*b10i + a02r*b20i + a00i*b00r + a01i*b10r
	addpd 	xmm5,xmm7		; xmm5 = a20r*b00r + a21r*b10r + a22r*b20r - a22i*b20i - a20i*b00i , a20r*b00i + a21r*b10i + a22r*b20i + a22i*b20r + a20i*b00r
	movsd 	xmm0,[eax+40]		; xmm0 = x,a02i	<(aa)->e[0][2].imag>
	movsd	xmm6,[eax+120]		; xmm6 = x,a21i	<(aa)->e[2][1].imag>
	movsd 	xmm7,[eax+88]		; xmm7 = x,a12i	<(aa)->e[1][2].imag>
	unpcklpd xmm0,xmm0		; xmm0 = a02i,a02i
	unpcklpd xmm6,xmm6		; xmm6 = a21i,a21i
	unpcklpd xmm7,xmm7		; xmm7 = a12i,a12i
	mulpd 	xmm0,xmm2		; xmm0 = -a02i*b20i,a02i*b20r
	mulpd	xmm6,xmm1		; xmm6 = -a21i*b10i,a21i*b10r
	mulpd	xmm7,xmm2		; xmm7 = -a12i*b20i,a12i*b20r
	addpd 	xmm3,xmm0		; xmm3 = a00r*b00r + a01r*b10r + a02r*b20r - a00i*b00i - a01i*b10i - a02i*b20i , a00r*b00i + a01r*b10i + a02r*b20i + a00i*b00r + a01i*b10r + a02i*b20r
		                        ;        a00r*b00r - a00i*b00i + a01r*b10r - a01i*b10i + a02r*b20r - a02i*b20i , a00r*b00i + a00i*b00r + a01r*b10i + a01i*b10r + a02r*b20i + a02i*b20r 
	addpd 	xmm4,xmm7		; xmm4 = a10r*b00r + a12r*b20r + a11r*b10r - a11i*b10i - a10i*b00i - a12i*b20i , a10r*b00i + a12r*b20i + a11r*b10i + a11i*b10r + a10i*b00r + a12i*b20r
					;        a10r*b00r - a10i*b00i + a11r*b10r - a11i*b10i + a12r*b20r - a12i*b20i , a10r*b00i + a10i*b00r + a11r*b10i + a11i*b10r + a12r*b20i + a12i*b20r
	addpd 	xmm5,xmm6		; xmm5 = a20r*b00r + a21r*b10r + a22r*b20r - a22i*b20i - a20i*b00i - a21i*b10i , a20r*b00i + a21r*b10i + a22r*b20i + a22i*b20r + a20i*b00r + a21i*b10r
        	                        ;        a20r*b00r - a20i*b00i + a21r*b10r - a21i*b10i + a22r*b20r - a22i*b20i , a20r*b00i + a20i*b00r + a21r*b10i + a21i*b10r + a22r*b20i + a22i*b20r 
	movupd	[ecx],xmm3		;	<(cc)->e[0][0]>
	movupd	[ecx+48],xmm4		;	<(cc)->e[1][0]>
	movupd	[ecx+96],xmm5		;	<(cc)->e[2][0]>
	;
	; REAL	Part II
	;
	movupd	xmm0,[ebx+48]		; xmm0 = b10r,b10i <(bb)->e[1][0]>
	movupd	xmm1,[ebx+64]		; xmm1 = b11r,b11i <(bb)->e[1][1]>
	movupd	xmm2,[ebx+80]		; xmm2 = b12r,b12i <(bb)->e[1][2]>
	xorpd	xmm0,[negate]		; xmm0 = b10r,-b10i	<_sse_sgn4>
	xorpd	xmm1,[negate]		; xmm1 = b11r,-b11i	<_sse_sgn4>
	xorpd	xmm2,[negate]		; xmm2 = b12r,-b12i	<_sse_sgn4>

	movsd	xmm3,[eax]		; xmm3 = x,a00r	<(aa)->e[0][0].real>
	movsd	xmm6,[eax+16]		; xmm6 = x,a01r	<(aa)->e[0][1].real>
	movsd	xmm4,[eax+48]		; xmm4 = x,a10r	<(aa)->e[1][0].real>
	movsd	xmm7,[eax+80]		; xmm7 = x,a12r	<(aa)->e[1][2].real>
	movsd	xmm5,[eax+96]		; xmm5 = x,a20r	<(aa)->e[2][0].real>
	unpcklpd xmm3,xmm3		; xmm3 = a00r,a00r
	unpcklpd xmm6,xmm6		; xmm6 = a01r,a01r
	unpcklpd xmm4,xmm4		; xmm4 = a10r,a10r
	mulpd	xmm3,xmm0		; xmm3 = a00r*b01r,a00r*b01i
	unpcklpd xmm7,xmm7		; xmm7 = a12r,a12r
	mulpd 	xmm6,xmm1		; xmm6 = a01r*b11r,a01r*b11i
	unpcklpd xmm5,xmm5		; xmm5 = a20r,a20r
	mulpd 	xmm4,xmm0		; xmm4 = a10r*b01r,a10r*b01i
	addpd 	xmm3,xmm6		; xmm3 = a00r*b01r+a01r*b11r,a00r*b01i+a01r*b11i
	mulpd 	xmm7,xmm2		; xmm7 = a12r*b21r,a12r*b21i
	mulpd 	xmm5,xmm0		; xmm5 = a20r*b01r,a20r*b01i
	addpd 	xmm4,xmm7		; xmm4 = a10r*b01r+a12r*b21r,a10r*b01i+a12r*b21i
	movsd 	xmm6,[eax+112]		; xmm6 = x,a21r	<(aa)->e[2][1].real>
	movsd 	xmm7,[eax+32]		; xmm7 = x,a02r <(aa)->e[0][2].real>
	unpcklpd xmm6,xmm6		; xmm6 = a21r,a21r
	unpcklpd xmm7,xmm7		; xmm7 = a02r,a02r
	mulpd 	xmm6,xmm1		; xmm6 = a21r*b11r,a21r*b11i
	mulpd 	xmm7,xmm2		; xmm7 = a02r*b21r,a02r*b21i
	addpd 	xmm5,xmm6		; xmm5 = a20r*b01r+a21r*b11r,a20r*b01i+a21r*b11i
	addpd 	xmm3,xmm7		; xmm3 = a00r*b01r + a01r*b11r + a02r*b21r , a00r*b01i + a01r*b11i + a02r*b21i
	movsd 	xmm6,[eax+64]		; xmm6 = x,a11r	<(aa)->e[1][1].real>
	movsd 	xmm7,[eax+128]		; xmm7 = x,a22r	<(aa)->e[2][2].real>
	unpcklpd xmm6,xmm6		; xmm6 = a11r,a11r
	unpcklpd xmm7,xmm7		; xmm7 = a22r,a22r
	mulpd 	xmm6,xmm1   		; xmm6 = a11r*b11r,a11r*b11i
	mulpd 	xmm7,xmm2   		; xmm7 = a22r*b21r,a22r*b21i
	addpd 	xmm4,xmm6		; xmm4 = a10r*b01r + a12r*b21r + a11r*b11r , a10r*b01i + a12r*b21i + a11r*b11i
	addpd 	xmm5,xmm7		; xmm5 = a20r*b01r + a21r*b11r + a22r*b21r , a20r*b01i + a21r*b11i + a22r*b21i
	;	
	; IMAGINARY	Part II
	;
	movsd	xmm6,[eax+8]		; xmm6 = x,a00i	<(aa)->e[0][0].imag>
	movsd	xmm7,[eax+72]		; xmm7 = x,a11i	<(aa)->e[1][1].imag>

	xorpd	xmm0,[negate]		; xmm0 = b01r,b01i	<_sse_sgn4>
	xorpd	xmm1,[negate]		; xmm1 = b11r,b11i	<_sse_sgn4>
	xorpd	xmm2,[negate]		; xmm2 = b21r,b21i	<_sse_sgn4>

	shufpd	xmm0,xmm0,0x1		; xmm0 = b01i,b01r
	shufpd	xmm1,xmm1,0x1		; xmm1 = b11i,b11r
	shufpd	xmm2,xmm2,0x1		; xmm2 = b21i,b21r
	unpcklpd xmm6,xmm6		; xmm6 = a00i,a00i
	unpcklpd xmm7,xmm7		; xmm7 = a11i,a11i
	mulpd	xmm6,xmm0		; xmm6 = -a00i*b01i,a00i*b01r
	mulpd	xmm7,xmm1		; xmm7 = -a11i*b11i,a11i*b11r
	addpd 	xmm3,xmm6		; xmm3 = a00r*b01r + a01r*b11r + a02r*b21r - a00i*b01i , a00r*b01i + a01r*b11i + a02r*b21i + a00i*b01r
	addpd 	xmm4,xmm7		; xmm4 = a10r*b01r + a12r*b21r + a11r*b11r - a11i*b11i , a10r*b01i + a12r*b21i + a11r*b11i + a11i*b11r
	movsd	xmm6,[eax+136]		; xmm6 = x,a22i	<(aa)->e[2][2].imag>
	movsd	xmm7,[eax+56]		; xmm7 = x,a10i	<(aa)->e[1][0].imag>
	unpcklpd xmm6,xmm6		; xmm6 = a22i,a22i
	unpcklpd xmm7,xmm7		; xmm7 = a10i,a10i
	mulpd 	xmm6,xmm2		; xmm6 = -a22i*b21i,a22i*b21r
	mulpd	xmm7,xmm0		; xmm7 = -a10i*b01i,a10i*b01r
	addpd 	xmm5,xmm6		; xmm5 = a20r*b01r + a21r*b11r + a22r*b21r - a22i*b21i , a20r*b01i + a21r*b11i + a22r*b21i + a22i*b21r
	addpd 	xmm4,xmm7		; xmm4 = a10r*b01r + a12r*b21r + a11r*b11r - a11i*b11i - a10i*b01i , a10r*b01i + a12r*b21i + a11r*b11i + a11i*b11r + a10i*b01r
	movsd	xmm6,[eax+24]		; xmm6 = x,a01i	<(aa)->e[0][1].imag>
	movsd	xmm7,[eax+104]		; xmm7 = x,a20i	<(aa)->e[2][0].imag>
	unpcklpd xmm6,xmm6		; xmm6 = a01i,a01i
	unpcklpd xmm7,xmm7		; xmm7 = a20i,a20i
	mulpd 	xmm6,xmm1		; xmm6 = -a01i*b11i,a01i*b11r
	mulpd 	xmm7,xmm0		; xmm7 = -a20i*b01i,a20i*b01r
	addpd 	xmm3,xmm6		; xmm3 = a00r*b01r + a01r*b11r + a02r*b21r - a00i*b01i - a01i*b11i , a00r*b01i + a01r*b11i + a02r*b21i + a00i*b01r + a01i*b11r
	addpd 	xmm5,xmm7		; xmm5 = a20r*b01r + a21r*b11r + a22r*b21r - a22i*b21i - a20i*b01i , a20r*b01i + a21r*b11i + a22r*b21i + a22i*b21r + a20i*b01r
	movsd 	xmm0,[eax+40]		; xmm0 = x,a02i	<(aa)->e[0][2].imag>
	movsd	xmm6,[eax+120]		; xmm6 = x,a21i	<(aa)->e[2][1].imag>
	movsd 	xmm7,[eax+88]		; xmm7 = x,a12i	<(aa)->e[1][2].imag>
	unpcklpd xmm0,xmm0		; xmm0 = a02i,a02i
	unpcklpd xmm6,xmm6		; xmm6 = a21i,a21i
	unpcklpd xmm7,xmm7		; xmm7 = a12i,a12i
	mulpd 	xmm0,xmm2		; xmm0 = -a02i*b21i,a02i*b21r
	mulpd	xmm6,xmm1		; xmm6 = -a21i*b11i,a21i*b11r
	mulpd	xmm7,xmm2		; xmm7 = -a12i*b21i,a12i*b21r
	addpd 	xmm3,xmm0		; xmm3 = a00r*b01r + a01r*b11r + a02r*b21r - a00i*b01i - a01i*b11i - a02i*b21i , a00r*b01i + a01r*b11i + a02r*b21i + a00i*b01r + a01i*b11r + a02i*b21r
		                        ;        a00r*b01r - a00i*b01i + a01r*b11r - a01i*b11i + a02r*b21r - a02i*b21i , a00r*b01i + a00i*b01r + a01r*b11i + a01i*b11r + a02r*b21i + a02i*b21r 
	addpd 	xmm4,xmm7		; xmm4 = a10r*b01r + a12r*b21r + a11r*b11r - a11i*b11i - a10i*b01i - a12i*b21i , a10r*b01i + a12r*b21i + a11r*b11i + a11i*b11r + a10i*b01r + a12i*b21r
					;        a10r*b01r - a10i*b01i + a11r*b11r - a11i*b11i + a12r*b21r - a12i*b21i , a10r*b01i + a10i*b01r + a11r*b11i + a11i*b11r + a12r*b21i + a12i*b21r
	addpd 	xmm5,xmm6		; xmm5 = a20r*b01r + a21r*b11r + a22r*b21r - a22i*b21i - a20i*b01i - a21i*b11i , a20r*b01i + a21r*b11i + a22r*b21i + a22i*b21r + a20i*b01r + a21i*b11r
        	                        ;        a20r*b01r - a20i*b01i + a21r*b11r - a21i*b11i + a22r*b21r - a22i*b21i , a20r*b01i + a20i*b01r + a21r*b11i + a21i*b11r + a22r*b21i + a22i*b21r 
	movupd	[ecx+16],xmm3		;	<(cc)->e[0][1]>
	movupd	[ecx+64],xmm4		;	<(cc)->e[1][1]>
	movupd	[ecx+112],xmm5		;	<(cc)->e[2][1]>
	;
	; REAL	Part III
	;
	movupd	xmm0,[ebx+96]		; xmm0 = b20r,b20i <(bb)->e[2][0]>
	movupd	xmm1,[ebx+112]		; xmm1 = b21r,b21i <(bb)->e[2][1]>
	movupd	xmm2,[ebx+128]		; xmm2 = b22r,b22i <(bb)->e[2][2]>
	xorpd	xmm0,[negate]		; xmm0 = b20r,-b20i	<_sse_sgn4>
	xorpd	xmm1,[negate]		; xmm1 = b21r,-b21i	<_sse_sgn4>
	xorpd	xmm2,[negate]		; xmm2 = b22r,-b22i	<_sse_sgn4>

	movsd	xmm3,[eax]		; xmm3 = x,a00r	<(aa)->e[0][0].real>
	movsd	xmm6,[eax+16]		; xmm6 = x,a01r	<(aa)->e[0][1].real>
	movsd	xmm4,[eax+48]		; xmm4 = x,a10r	<(aa)->e[1][0].real>
	movsd	xmm7,[eax+80]		; xmm7 = x,a12r	<(aa)->e[1][2].real>
	movsd	xmm5,[eax+96]		; xmm5 = x,a20r	<(aa)->e[2][0].real>
	unpcklpd xmm3,xmm3		; xmm3 = a00r,a00r
	unpcklpd xmm6,xmm6		; xmm6 = a01r,a01r
	unpcklpd xmm4,xmm4		; xmm4 = a10r,a10r
	mulpd	xmm3,xmm0		; xmm3 = a00r*b02r,a00r*b02i
	unpcklpd xmm7,xmm7		; xmm7 = a12r,a12r
	mulpd 	xmm6,xmm1		; xmm6 = a01r*b12r,a01r*b12i
	unpcklpd xmm5,xmm5		; xmm5 = a20r,a20r
	mulpd 	xmm4,xmm0		; xmm4 = a10r*b02r,a10r*b02i
	addpd 	xmm3,xmm6		; xmm3 = a00r*b02r+a01r*b12r,a00r*b02i+a01r*b12i
	mulpd 	xmm7,xmm2		; xmm7 = a12r*b22r,a12r*b22i
	mulpd 	xmm5,xmm0		; xmm5 = a20r*b02r,a20r*b02i
	addpd 	xmm4,xmm7		; xmm4 = a10r*b02r+a12r*b22r,a10r*b02i+a12r*b22i
	movsd 	xmm6,[eax+112]		; xmm6 = x,a21r	<(aa)->e[2][1].real>
	movsd 	xmm7,[eax+32]		; xmm7 = x,a02r <(aa)->e[0][2].real>
	unpcklpd xmm6,xmm6		; xmm6 = a21r,a21r
	unpcklpd xmm7,xmm7		; xmm7 = a02r,a02r
	mulpd 	xmm6,xmm1		; xmm6 = a21r*b12r,a21r*b12i
	mulpd 	xmm7,xmm2		; xmm7 = a02r*b22r,a02r*b22i
	addpd 	xmm5,xmm6		; xmm5 = a20r*b02r+a21r*b12r,a20r*b02i+a21r*b12i
	addpd 	xmm3,xmm7		; xmm3 = a00r*b02r + a01r*b12r + a02r*b22r , a00r*b02i + a01r*b12i + a02r*b22i
	movsd 	xmm6,[eax+64]		; xmm6 = x,a11r	<(aa)->e[1][1].real>
	movsd 	xmm7,[eax+128]		; xmm7 = x,a22r	<(aa)->e[2][2].real>
	unpcklpd xmm6,xmm6		; xmm6 = a11r,a11r
	unpcklpd xmm7,xmm7		; xmm7 = a22r,a22r
	mulpd 	xmm6,xmm1   		; xmm6 = a11r*b12r,a11r*b12i
	mulpd 	xmm7,xmm2   		; xmm7 = a22r*b22r,a22r*b22i
	addpd 	xmm4,xmm6		; xmm4 = a10r*b02r + a12r*b22r + a11r*b12r , a10r*b02i + a12r*b22i + a11r*b12i
	addpd 	xmm5,xmm7		; xmm5 = a20r*b02r + a21r*b12r + a22r*b22r , a20r*b02i + a21r*b12i + a22r*b22i
	;	
	; IMAGINARY	Part III
	;
	movsd	xmm6,[eax+8]		; xmm6 = x,a00i	<(aa)->e[0][0].imag>
	movsd	xmm7,[eax+72]		; xmm7 = x,a11i	<(aa)->e[1][1].imag>

	xorpd	xmm0,[negate]		; xmm0 = b02r,b02i	<_sse_sgn4>
	xorpd	xmm1,[negate]		; xmm1 = b12r,b12i	<_sse_sgn4>
	xorpd	xmm2,[negate]		; xmm2 = b22r,b22i	<_sse_sgn4>

	shufpd	xmm0,xmm0,0x1		; xmm0 = b02i,b02r
	shufpd	xmm1,xmm1,0x1		; xmm1 = b12i,b12r
	shufpd	xmm2,xmm2,0x1		; xmm2 = b22i,b22r
	unpcklpd xmm6,xmm6		; xmm6 = a00i,a00i
	unpcklpd xmm7,xmm7		; xmm7 = a11i,a11i
	mulpd	xmm6,xmm0		; xmm6 = -a00i*b02i,a00i*b02r
	mulpd	xmm7,xmm1		; xmm7 = -a11i*b12i,a11i*b12r
	addpd 	xmm3,xmm6		; xmm3 = a00r*b02r + a01r*b12r + a02r*b22r - a00i*b02i , a00r*b02i + a01r*b12i + a02r*b22i + a00i*b02r
	addpd 	xmm4,xmm7		; xmm4 = a10r*b02r + a12r*b22r + a11r*b12r - a11i*b12i , a10r*b02i + a12r*b22i + a11r*b12i + a11i*b12r
	movsd	xmm6,[eax+136]		; xmm6 = x,a22i	<(aa)->e[2][2].imag>
	movsd	xmm7,[eax+56]		; xmm7 = x,a10i	<(aa)->e[1][0].imag>
	unpcklpd xmm6,xmm6		; xmm6 = a22i,a22i
	unpcklpd xmm7,xmm7		; xmm7 = a10i,a10i
	mulpd 	xmm6,xmm2		; xmm6 = -a22i*b22i,a22i*b22r
	mulpd	xmm7,xmm0		; xmm7 = -a10i*b02i,a10i*b02r
	addpd 	xmm5,xmm6		; xmm5 = a20r*b02r + a21r*b12r + a22r*b22r - a22i*b22i , a20r*b02i + a21r*b12i + a22r*b22i + a22i*b22r
	addpd 	xmm4,xmm7		; xmm4 = a10r*b02r + a12r*b22r + a11r*b12r - a11i*b12i - a10i*b02i , a10r*b02i + a12r*b22i + a11r*b12i + a11i*b12r + a10i*b02r
	movsd	xmm6,[eax+24]		; xmm6 = x,a01i	<(aa)->e[0][1].imag>
	movsd	xmm7,[eax+104]		; xmm7 = x,a20i	<(aa)->e[2][0].imag>
	unpcklpd xmm6,xmm6		; xmm6 = a01i,a01i
	unpcklpd xmm7,xmm7		; xmm7 = a20i,a20i
	mulpd 	xmm6,xmm1		; xmm6 = -a01i*b12i,a01i*b12r
	mulpd 	xmm7,xmm0		; xmm7 = -a20i*b02i,a20i*b02r
	addpd 	xmm3,xmm6		; xmm3 = a00r*b02r + a01r*b12r + a02r*b22r - a00i*b02i - a01i*b12i , a00r*b02i + a01r*b12i + a02r*b22i + a00i*b02r + a01i*b12r
	addpd 	xmm5,xmm7		; xmm5 = a20r*b02r + a21r*b12r + a22r*b22r - a22i*b22i - a20i*b02i , a20r*b02i + a21r*b12i + a22r*b22i + a22i*b22r + a20i*b02r
	movsd 	xmm0,[eax+40]		; xmm0 = x,a02i	<(aa)->e[0][2].imag>
	movsd	xmm6,[eax+120]		; xmm6 = x,a21i	<(aa)->e[2][1].imag>
	movsd 	xmm7,[eax+88]		; xmm7 = x,a12i	<(aa)->e[1][2].imag>
	unpcklpd xmm0,xmm0		; xmm0 = a02i,a02i
	unpcklpd xmm6,xmm6		; xmm6 = a21i,a21i
	unpcklpd xmm7,xmm7		; xmm7 = a12i,a12i
	mulpd 	xmm0,xmm2		; xmm0 = -a02i*b22i,a02i*b22r
	mulpd	xmm6,xmm1		; xmm6 = -a21i*b12i,a21i*b12r
	mulpd	xmm7,xmm2		; xmm7 = -a12i*b22i,a12i*b22r
	addpd 	xmm3,xmm0		; xmm3 = a00r*b02r + a01r*b12r + a02r*b22r - a00i*b02i - a01i*b12i - a02i*b22i , a00r*b02i + a01r*b12i + a02r*b22i + a00i*b02r + a01i*b12r + a02i*b22r
		                        ;        a00r*b02r - a00i*b02i + a01r*b12r - a01i*b12i + a02r*b22r - a02i*b22i , a00r*b02i + a00i*b02r + a01r*b12i + a01i*b12r + a02r*b22i + a02i*b22r 
	addpd 	xmm4,xmm7		; xmm4 = a10r*b02r + a12r*b22r + a11r*b12r - a11i*b12i - a10i*b02i - a12i*b22i , a10r*b02i + a12r*b22i + a11r*b12i + a11i*b12r + a10i*b02r + a12i*b22r
					;        a10r*b02r - a10i*b02i + a11r*b12r - a11i*b12i + a12r*b22r - a12i*b22i , a10r*b02i + a10i*b02r + a11r*b12i + a11i*b12r + a12r*b22i + a12i*b22r
	addpd 	xmm5,xmm6		; xmm5 = a20r*b02r + a21r*b12r + a22r*b22r - a22i*b22i - a20i*b02i - a21i*b12i , a20r*b02i + a21r*b12i + a22r*b22i + a22i*b22r + a20i*b02r + a21i*b12r
        	                        ;        a20r*b02r - a20i*b02i + a21r*b12r - a21i*b12i + a22r*b22r - a22i*b22i , a20r*b02i + a20i*b02r + a21r*b12i + a21i*b12r + a22r*b22i + a22i*b22r 
	movupd	[ecx+32],xmm3		;	<(cc)->e[0][2]>
	movupd	[ecx+80],xmm4		;	<(cc)->e[1][2]>
	movupd	[ecx+128],xmm5		;	<(cc)->e[2][2]>
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

