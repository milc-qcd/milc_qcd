;
; su3_projector( su3_vector *a, su3_vector *b, su3_matrix *c)
;
; | a00  |*| b00 b01 b02 | | c00 c01 c02 |
; | a10  |                =| c10 c11 c12 |
; | a20  |                 | c20 c21 c22 |
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

	; First load part of vector a into registers

	; a00
	movsd		xmm3,[eax]			; xmm3=x,a00r		<(aa)->c[0].real>
	unpcklpd	xmm3,xmm3			; xmm3=a00r,a00r
	movsd		xmm4,[eax+8]			; xmm4=x,a00i		<(aa)->c[0].imag>
	unpcklpd	xmm4,xmm4			; xmm4=a00i,a00i
	; a10
	movsd		xmm5,[eax+16]			; xmm5=x,a10r		<(aa)->c[1].real>
	unpcklpd	xmm5,xmm5			; xmm5=a10r,a10r
	movsd		xmm6,[eax+24]			; xmm6=x,a10i		<(aa)->c[1].imag>
	unpcklpd	xmm6,xmm6			; xmm6=a10i,a10i
	; a20r
        movsd           xmm7,[eax+32]                   ; xmm7=x,a20r           <(aa)->c[2].real>
        unpcklpd        xmm7,xmm7                       ; xmm7=a20r,a20r

	; First multiply all of vector a by b[0] and b[1]

	; a00 x b00
	movupd		xmm0,[ebx]			; xmm0=b00i,b00r	<(bb)->c[0]>
	movupd		xmm1,xmm3			; xmm1=a00r,a00r
	movupd		xmm2,xmm4			; xmm2=a00i,a00i
	mulpd		xmm1,xmm0			; xmm1=a00r*b00i,a00r*b00r
	mulpd		xmm2,xmm0			; xmm2=a00i*b00i,a00i*b00r
	xorpd		xmm1,[negate]			; xmm1=-a00r*b00i,a00r*b00r	<_sse_sgn4>
	shufpd		xmm2,xmm2,0x1			; xmm2=a00i*b00r,a00i*b00i
	addpd		xmm2,xmm1			; xmm2=-a00r*b00i+a00i*b00r,a00r*b00r+a00i*b00i
	movupd		[ecx],xmm2			;<(cc)->e[0][0]>

	; a10 x b00
	movupd		xmm1,xmm5			; xmm1=a10r,a10r
	movupd		xmm2,xmm6			; xmm2=a10i,a10i
	mulpd		xmm2,xmm0			; xmm2=a10i*b00i,a10i*b00r
	mulpd		xmm1,xmm0			; xmm1=a10r*b00i,a10r*b00r
	xorpd		xmm1,[negate]			; xmm1=-a10r*b00i,a10r*b00r	<_sse_sgn4>
	shufpd		xmm2,xmm2,0x1			; xmm2=a10i*b00r,a10i*b00i
	addpd		xmm2,xmm1			; xmm2=-a10r*b00i+a10i*b00r,a10r*b00r+a10i*b00i
	movupd		[ecx+48],xmm2			;<(cc)->e[1][0]>

	; a20 x b00
	movupd		xmm1,xmm7			; xmm1=a20r,a20r
	movsd		xmm2,[eax+40]			; xmm2=x,a20i		<(aa)->c[2].imag>
	unpcklpd	xmm2,xmm2			; xmm2=a20i,a20i
	mulpd		xmm2,xmm0			; xmm2=a20i*b00i,a20i*b00r
	mulpd		xmm1,xmm0			; xmm1=a20r*b00i,a20r*b00r
	xorpd		xmm1,[negate]			; xmm1=-a20r*b00i,a20r*b00r	<_sse_sgn4>
	shufpd		xmm2,xmm2,0x1			; xmm2=a20i*b00r,a20i*b00i
	addpd		xmm2,xmm1			; xmm2=-a20r*b00i+a20i*b00r,a20r*b00r+a20i*b00i
	movupd		[ecx+96],xmm2			;<(cc)->e[2][0]>
	
	; a00 x b01
	movupd		xmm0,[ebx+16]			; xmm0=b01i,b01r	<(bb)->c[1]>
	movupd		xmm1,xmm3			; xmm1=a00r,a00r
	movupd		xmm2,xmm4			; xmm2=a00i,a00i
	mulpd		xmm1,xmm0			; xmm1=a00r*b01i,a00r*b01r
	mulpd		xmm2,xmm0			; xmm2=a00i*b01i,a00i*b01r
	xorpd		xmm1,[negate]			; xmm1=-a00r*b01i,a00r*b01r	<_sse_sgn4>
	shufpd		xmm2,xmm2,0x1			; xmm2=a00i*b01r,a00i*b01i
	addpd		xmm2,xmm1			; xmm2=-a00r*b01i+a00i*b01r,a00r*b01r+a00i*b01i
	movupd		[ecx+16],xmm2			;<(cc)->e[0][1]>

	; a10 x b01
	movupd		xmm1,xmm5			; xmm1=a10r,a10r
        movupd          xmm2,xmm6                       ; xmm2=a10i,a10i
	mulpd		xmm2,xmm0			; xmm2=a10i*b01i,a10i*b01r
	mulpd		xmm1,xmm0			; xmm1=a10r*b01i,a10r*b01r
	xorpd		xmm1,[negate]			; xmm1=-a10r*b01i,a10r*b01r	<_sse_sgn4>
	shufpd		xmm2,xmm2,0x1			; xmm2=a10i*b01r,a10i*b01i
	addpd		xmm2,xmm1			; xmm2=-a10r*b01i+a10i*b01r,a10r*b01r+a10i*b01i
	movupd		[ecx+64],xmm2			;<(cc)->e[1][1]>

	; a20 x b01
        movupd          xmm1,xmm7                       ; xmm1=a20r,a20r
	movsd		xmm2,[eax+40]			; xmm2=x,a20i		<(aa)->c[2].imag>
	unpcklpd	xmm2,xmm2			; xmm2=a20i,a20i
	mulpd		xmm2,xmm0			; xmm2=a20i*b01i,a20i*b01r
	mulpd		xmm1,xmm0			; xmm1=a20r*b01i,a20r*b01r
	xorpd		xmm1,[negate]			; xmm1=-a20r*b01i,a20r*b01r	<_sse_sgn4>
	shufpd		xmm2,xmm2,0x1			; xmm2=a20i*b01r,a20i*b01i
	addpd		xmm2,xmm1			; xmm2=-a20r*b01i+a20i*b01r,a20r*b01r+a20i*b01i
	movupd		[ecx+112],xmm2			;<(cc)->e[2][1]>

	; Next, multiply b[2] by a[0] and a[1]

	movupd		xmm0,[ebx+32]			; xmm0=b02i,b02r	<(bb)->c[2]>
	
	; a00 x b02
	movupd		xmm1,xmm3			; xmm1=a00r,a00r
	movupd		xmm2,xmm4			; xmm2=a00i,a00i
	mulpd		xmm1,xmm0			; xmm1=a00r*b02i,a00r*b02r
	mulpd		xmm2,xmm0			; xmm2=a00i*b02i,a00i*b02r
	xorpd		xmm1,[negate]			; xmm1=-a00r*b02i,a00r*b02r	<_sse_sgn4>
	shufpd		xmm2,xmm2,0x1			; xmm2=a00i*b02r,a00i*b02i
	addpd		xmm2,xmm1			; xmm2=-a00r*b02i+a00i*b02r,a00r*b02r+a00i*b02i
	movupd		[ecx+32],xmm2			;<(cc)->e[0][2]>

	; a10 x b02
	movupd		xmm1,xmm5			; xmm1=a10r,a10r
	movupd		xmm2,xmm6			; xmm2=a10i,a10i
	mulpd		xmm1,xmm0			; xmm1=a10r*b02i,a10r*b02r
	mulpd		xmm2,xmm0			; xmm2=a10i*b02i,a10i*b02r
	xorpd		xmm1,[negate]			; xmm1=-a10r*b02i,a10r*b02r	<_sse_sgn4>
	shufpd		xmm2,xmm2,0x1			; xmm2=a10i*b02r,a10i*b02i
	addpd		xmm2,xmm1			; xmm2=-a10r*b02i+a10i*b02r,a10r*b02r+a10i*b02i
	movupd		[ecx+80],xmm2			;<(cc)->e[1][2]>
	
	; Finally, do a[2] X b[2]

	; a20 x b02
        movupd          xmm1,xmm7                       ; xmm1=a20r,a20r
	movsd		xmm2,[eax+40]			; xmm2=x,a20i		<(aa)->c[2].imag>
	unpcklpd	xmm2,xmm2			; xmm2=a20i,a20i
	mulpd		xmm1,xmm0			; xmm1=a20r*b02i,a20r*b02r
	mulpd		xmm2,xmm0			; xmm2=a20i*b02i,a20i*b02r
	xorpd		xmm1,[negate]			; xmm1=-a20r*b02i,a20r*b02r	<_sse_sgn4>
	shufpd		xmm2,xmm2,0x1			; xmm2=a20i*b02r,a20i*b02i
	addpd		xmm2,xmm1			; xmm2=-a20r*b02i+a20i*b02r,a20r*b02r+a20i*b02i
	movupd		[ecx+128],xmm2			;<(cc)->e[2][2]>

here:	pop	ecx
	pop	ebx
	pop	eax
	mov	esp,ebp
	pop	ebp
	ret
	
	align		16
negate:	dd		0x00000000
	dd		0x00000000
	dd		0x00000000
	dd		0x80000000
