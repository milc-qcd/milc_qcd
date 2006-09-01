;
; mult_adj_su3_mat_vec( su3_matrix *a, su3_vector *b, su3_vector *c)
;
; Opteron SSE code for 3x3 matrix-vector multiply
; c <- (A^H)b
; 
; Steve Whalen
; Cray, Inc.

	bits		64

global mult_adj_su3_mat_vec
mult_adj_su3_mat_vec:

	;; rows 1 and 2

        ; real part of b
	movlps		xmm15,[rsi]		;              b0.i  b0.r	<(bb)->c[0]>
        movaps          xmm2,xmm15
	unpcklps	xmm2,xmm2		;  b0.i  b0.i  b0.r  b0.r
	movhlps		xmm6,xmm2		;              b0.i  b0.i
	movlhps		xmm2,xmm2		;  b0.r  b0.r  b0.r  b0.r
	movhps		xmm15,[rsi+8]		;  b1.i  b1.r			<(bb)->c[1]>
        movaps          xmm3,xmm15
	unpckhps	xmm3,xmm3		;  b1.i  b1.i  b1.r  b1.r
	movhlps		xmm7,xmm3		;              b1.i  b1.i
	movlhps		xmm3,xmm3		;  b1.r  b1.r  b1.r  b1.r
	movhps		xmm11,[rsi+16]		;  b2.i  b2.r			<(bb)->c[2]>
        movaps          xmm4,xmm11
	unpckhps	xmm4,xmm4		;  b2.i  b2.i  b2.r  b2.r
	movhlps		xmm9,xmm4		;              b2.i  b2.i
	movlhps		xmm4,xmm4		;  b2.r  b2.r  b2.r  b2.r

	movlps		xmm5,[rdi]		;             a00.i a00.r	<(aa)->e[0][0]>
	movhps		xmm5,[rdi+8]		; a01.i a01.r			<(aa)->e[0][1]>
	mulps           xmm2,xmm5
	movlps		xmm8,[rdi+24]		;             a10.i a10.r	<(aa)->e[1][0]>
	movhps		xmm8,[rdi+32]		; a11.i a11.r			<(aa)->e[1][1]>
	mulps           xmm3,xmm8
	movlps		xmm10,[rdi+48]		;             a20.i a20.r	<(aa)->e[2][0]>
	movhps		xmm10,[rdi+56]		; a21.i a21.r			<(aa)->e[2][1]>
        mulps           xmm4,xmm10
        addps           xmm2,xmm3
        addps           xmm2,xmm4

        ; imaginary part of b
        movlhps		xmm6,xmm6		;  b0.i  b0.i  b0.i  b0.i
        movlhps		xmm7,xmm7		;  b1.i  b1.i  b1.i  b1.i
        movlhps		xmm9,xmm9		;  b2.i  b2.i  b2.i  b2.i
        mulps           xmm6,xmm5
        mulps           xmm7,xmm8
	; xmm8 can now be modified
        mulps           xmm9,xmm10
	; xmm10 can now be modified
        addps           xmm6,xmm7
        addps           xmm6,xmm9
        xorps           xmm2,[negate]		; <_sse_sgn24>
        shufps          xmm6,xmm6,0xb1
	addps           xmm2,xmm6
	movlps		[rdx],xmm2		; <(cc)->c[0]>
	movhps		[rdx+8],xmm2		; <(cc)->c[1]>


	; third row
	movlps		xmm0,[rdi+16]		;             a02.i a02.r	<(aa)->e[0][2]>
	unpcklps	xmm0,xmm0		; a02.i a02.i a02.r a02.r
	movaps		xmm12,xmm0
	movlps		xmm1,[rdi+40]		;             a12.i a12.r	<(aa)->e[1][2]>
	unpcklps	xmm1,xmm1		; a12.i a12.i a12.r a12.r
	movlhps		xmm12,xmm1		; a12.r a12.r a02.r a02.r
	movhlps		xmm1,xmm0		; a12.i a12.i a02.i a02.i

	mulps		xmm12,xmm15
	mulps		xmm1,xmm15
	shufps		xmm1,xmm1,0xb1
	xorps		xmm1,[negate]		; <_sse_sgn24>
	addps		xmm12,xmm1
	; xmm12 now contains a02*b0 and a12*b1

	movlps		xmm8,[rdi+64]		;             a22.i a22.r	<(aa)->e[2][2]>
	unpcklps	xmm8,xmm8		; a22.i a22.i a22.r a22.r
	movhlps		xmm11,xmm11		;  b2.i  b2.r  b2.i  b2.r
	mulps		xmm11,xmm8
	shufps		xmm11,xmm11,0xb4
	xorps		xmm11,[neg4]		; <_sse_sgn4>
	addps		xmm12,xmm11

	; need to add high bits to low bits in xmm12
	movhlps		xmm5,xmm12
	addps		xmm12,xmm5
	movlps		[rdx+16],xmm12		; <(cc)->c[2]>

	; Done!

	; *******************************************************************	

here:	ret
	
	align		16
negate:	dd		0x00000000
	dd		0x80000000
	dd		0x00000000
	dd		0x80000000

	align		16
neg4:	dd		0x00000000
	dd		0x00000000
	dd		0x00000000
	dd		0x80000000
