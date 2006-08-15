;
; mult_adj_su3_mat_4vec( su3_matrix *a[4], su3_vector *b,
;                        su3_vector *c0, su3_vector *c1,
;                        su3_vector *c2, su3_vector *c3)
;
; Multiply the adjoint of each of four input matrices by an input vector, 
; storing the resulting vectors in 4 separate destinations.
;
; Steve Whalen
; Cray, Inc.

	bits		64
	
global mult_adj_su3_mat_4vec
mult_adj_su3_mat_4vec:

	;; first product: c0 <- (A[0]^H)b

	;; rows 1 and 2

        ; real part of b
	movlps		xmm14,[rsi]		;              b0.i  b0.r	<(bb)->c[0]>
        movaps          xmm8,xmm14
	unpcklps	xmm8,xmm8		;  b0.i  b0.i  b0.r  b0.r
	movhlps		xmm11,xmm8		;              b0.i  b0.i
	movlhps		xmm8,xmm8		;  b0.r  b0.r  b0.r  b0.r
	movhps		xmm14,[rsi+8]		;  b1.i  b1.r			<(bb)->c[1]>
        movaps          xmm9,xmm14
	unpckhps	xmm9,xmm9		;  b1.i  b1.i  b1.r  b1.r
	movhlps		xmm12,xmm9		;              b1.i  b1.i
	movlhps		xmm9,xmm9		;  b1.r  b1.r  b1.r  b1.r
	movhps		xmm15,[rsi+16]		;  b2.i  b2.r			<(bb)->c[2]>
        movaps          xmm10,xmm15
	unpckhps	xmm10,xmm10		;  b2.i  b2.i  b2.r  b2.r
	movhlps		xmm13,xmm10		;              b2.i  b2.i
	movlhps		xmm10,xmm10		;  b2.r  b2.r  b2.r  b2.r

	movaps		xmm0,xmm8
	movlps		xmm3,[rdi]		;             a00.i a00.r	<(aa)[0].e[0][0]>
	movhps		xmm3,[rdi+8]		; a01.i a01.r			<(aa)[0].e[0][1]>
	mulps           xmm0,xmm3
	movaps		xmm1,xmm9
	movlps		xmm4,[rdi+24]		;             a10.i a10.r	<(aa)[0].e[1][0]>
	movhps		xmm4,[rdi+32]		; a11.i a11.r			<(aa)[0].e[1][1]>
	mulps           xmm1,xmm4
	movaps		xmm2,xmm10
	movlps		xmm5,[rdi+48]		;             a20.i a20.r	<(aa)[0].e[2][0]>
	movhps		xmm5,[rdi+56]		; a21.i a21.r			<(aa)[0].e[2][1]>
        mulps           xmm2,xmm5
        addps           xmm0,xmm1
        addps           xmm0,xmm2

        ; imaginary part of b
        movlhps		xmm11,xmm11		;  b0.i  b0.i  b0.i  b0.i
        movlhps		xmm12,xmm12		;  b1.i  b1.i  b1.i  b1.i
        movlhps		xmm13,xmm13		;  b2.i  b2.i  b2.i  b2.i
        mulps           xmm3,xmm11
        mulps           xmm4,xmm12
        mulps           xmm5,xmm13
        addps           xmm3,xmm4
        addps           xmm3,xmm5
        xorps           xmm0,[negate]		; <_sse_sgn24>
        shufps          xmm3,xmm3,0xb1
	addps           xmm0,xmm3
	movlps		[rdx],xmm0		; <(cc0)->c[0]>
	movhps		[rdx+8],xmm0		; <(cc0)->c[1]>


	; third row
	movlps		xmm1,[rdi+16]		;             a02.i a02.r	<(aa)[0].e[0][2]>
	unpcklps	xmm1,xmm1		; a02.i a02.i a02.r a02.r
	movaps		xmm6,xmm1
	movlps		xmm2,[rdi+40]		;             a12.i a12.r	<(aa)[0].e[1][2]>
	unpcklps	xmm2,xmm2		; a12.i a12.i a12.r a12.r
	movlhps		xmm6,xmm2		; a12.r a12.r a02.r a02.r
	movhlps		xmm2,xmm1		; a12.i a12.i a02.i a02.i

	mulps		xmm6,xmm14
	mulps		xmm2,xmm14
	shufps		xmm2,xmm2,0xb1
	xorps		xmm2,[negate]		; <_sse_sgn24>
	addps		xmm6,xmm2
	; xmm6 now contains a02*b0 and a12*b1

	movlps		xmm7,[rdi+64]		;             a22.i a22.r	<(aa)[0].e[2][2]>
	unpcklps	xmm7,xmm7		; a22.i a22.i a22.r a22.r
	movhlps		xmm15,xmm15		;  b2.i  b2.r  b2.i  b2.r
	mulps		xmm7,xmm15
	shufps		xmm7,xmm7,0xb4
	xorps		xmm7,[neg4]		; <_sse_sgn4>
	addps		xmm6,xmm7

	; need to add high bits to low bits in xmm6
	movhlps		xmm5,xmm6
	addps		xmm6,xmm5
	movlps		[rdx+16],xmm6		; <(cc0)->c[2]>

	; At this point, xmm8-xmm15 contain values for reuse:
	; xmm8  = [b0.r, b0.r, b0.r, b0.r]
	; xmm9  = [b1.r, b1.r, b1.r, b1.r]
	; xmm10 = [b2.r, b2.r, b2.r, b2.r]
	; xmm11 = [b0.i, b0.i, b0.i, b0.i]
	; xmm12 = [b1.i, b1.i, b1.i, b1.i]
	; xmm13 = [b2.i, b2.i, b2.i, b2.i]
	; xmm14 = [b1.i, b1.r, b0.i, b0.r]
	; xmm15 = [b2.i, b2.r, b2.i, b2.r]

	; *******************************************************************	
	; *******************************************************************	

	;; second product: c1 <- (A[1]^H)b

	add		rdi,72

	;; rows 1 and 2

	movaps		xmm0,xmm8
	movlps		xmm3,[rdi]		;             a00.i a00.r	<(aa)[1].e[0][0]>
	movhps		xmm3,[rdi+8]		; a01.i a01.r			<(aa)[1].e[0][1]>
	mulps           xmm0,xmm3
	movaps		xmm1,xmm9
	movlps		xmm4,[rdi+24]		;             a10.i a10.r	<(aa)[1].e[1][0]>
	movhps		xmm4,[rdi+32]		; a11.i a11.r			<(aa)[1].e[1][1]>
	mulps           xmm1,xmm4
	movaps		xmm2,xmm10
	movlps		xmm5,[rdi+48]		;             a20.i a20.r	<(aa)[1].e[2][0]>
	movhps		xmm5,[rdi+56]		; a21.i a21.r			<(aa)[1].e[2][1]>
        mulps           xmm2,xmm5
        addps           xmm0,xmm1
        addps           xmm0,xmm2

        ; imaginary part of b
        mulps           xmm3,xmm11
        mulps           xmm4,xmm12
        mulps           xmm5,xmm13
        addps           xmm3,xmm4
        addps           xmm3,xmm5
        xorps           xmm0,[negate]		; <_sse_sgn24>
        shufps          xmm3,xmm3,0xb1
	addps           xmm0,xmm3
	movlps		[rcx],xmm0		; <(cc1)->c[0]>
	movhps		[rcx+8],xmm0		; <(cc1)->c[1]>


	; third row
	movlps		xmm1,[rdi+16]		;             a02.i a02.r	<(aa)[1].e[0][2]>
	unpcklps	xmm1,xmm1		; a02.i a02.i a02.r a02.r
	movaps		xmm6,xmm1
	movlps		xmm2,[rdi+40]		;             a12.i a12.r	<(aa)[1].e[1][2]>
	unpcklps	xmm2,xmm2		; a12.i a12.i a12.r a12.r
	movlhps		xmm6,xmm2		; a12.r a12.r a02.r a02.r
	movhlps		xmm2,xmm1		; a12.i a12.i a02.i a02.i

	mulps		xmm6,xmm14
	mulps		xmm2,xmm14
	shufps		xmm2,xmm2,0xb1
	xorps		xmm2,[negate]		; <_sse_sgn24>
	addps		xmm6,xmm2
	; xmm6 now contains a02*b0 and a12*b1

	movlps		xmm7,[rdi+64]		;             a22.i a22.r	<(aa)[1].e[2][2]>
	unpcklps	xmm7,xmm7		; a22.i a22.i a22.r a22.r
	mulps		xmm7,xmm15
	shufps		xmm7,xmm7,0xb4
	xorps		xmm7,[neg4]		; <_sse_sgn4>
	addps		xmm6,xmm7

	; need to add high bits to low bits in xmm6
	movhlps		xmm5,xmm6
	addps		xmm6,xmm5
	movlps		[rcx+16],xmm6		; <(cc1)->c[2]>


	; *******************************************************************	
	; *******************************************************************	

	;; third product: c2 <- (A[2]^H)b

	add		rdi,72

	;; rows 1 and 2

	movaps		xmm0,xmm8
	movlps		xmm3,[rdi]		;             a00.i a00.r	<(aa)[2].e[0][0]>
	movhps		xmm3,[rdi+8]		; a01.i a01.r			<(aa)[2].e[0][1]>
	mulps           xmm0,xmm3
	movaps		xmm1,xmm9
	movlps		xmm4,[rdi+24]		;             a10.i a10.r	<(aa)[2].e[1][0]>
	movhps		xmm4,[rdi+32]		; a11.i a11.r			<(aa)[2].e[1][1]>
	mulps           xmm1,xmm4
	movaps		xmm2,xmm10
	movlps		xmm5,[rdi+48]		;             a20.i a20.r	<(aa)[2].e[2][0]>
	movhps		xmm5,[rdi+56]		; a21.i a21.r			<(aa)[2].e[2][1]>
        mulps           xmm2,xmm5
        addps           xmm0,xmm1
        addps           xmm0,xmm2

        ; imaginary part of b
        mulps           xmm3,xmm11
        mulps           xmm4,xmm12
        mulps           xmm5,xmm13
        addps           xmm3,xmm4
        addps           xmm3,xmm5
        xorps           xmm0,[negate]		; <_sse_sgn24>
        shufps          xmm3,xmm3,0xb1
	addps           xmm0,xmm3
	movlps		[r8],xmm0		; <(cc2)->c[0]>
	movhps		[r8+8],xmm0		; <(cc2)->c[1]>


	; third row
	movlps		xmm1,[rdi+16]		;             a02.i a02.r	<(aa)[2].e[0][2]>
	unpcklps	xmm1,xmm1		; a02.i a02.i a02.r a02.r
	movaps		xmm6,xmm1
	movlps		xmm2,[rdi+40]		;             a12.i a12.r	<(aa)[2].e[1][2]>
	unpcklps	xmm2,xmm2		; a12.i a12.i a12.r a12.r
	movlhps		xmm6,xmm2		; a12.r a12.r a02.r a02.r
	movhlps		xmm2,xmm1		; a12.i a12.i a02.i a02.i

	mulps		xmm6,xmm14
	mulps		xmm2,xmm14
	shufps		xmm2,xmm2,0xb1
	xorps		xmm2,[negate]		; <_sse_sgn24>
	addps		xmm6,xmm2
	; xmm6 now contains a02*b0 and a12*b1

	movlps		xmm7,[rdi+64]		;             a22.i a22.r	<(aa)[2].e[2][2]>
	unpcklps	xmm7,xmm7		; a22.i a22.i a22.r a22.r
	mulps		xmm7,xmm15
	shufps		xmm7,xmm7,0xb4
	xorps		xmm7,[neg4]		; <_sse_sgn4>
	addps		xmm6,xmm7

	; need to add high bits to low bits in xmm6
	movhlps		xmm5,xmm6
	addps		xmm6,xmm5
	movlps		[r8+16],xmm6		; <(cc2)->c[2]>


	; *******************************************************************	
	; *******************************************************************	

	;; fourth product: c3 <- (A[3]^H)b

	add		rdi,72

	;; rows 1 and 2

	movaps		xmm0,xmm8
	movlps		xmm3,[rdi]		;             a00.i a00.r	<(aa)[3].e[0][0]>
	movhps		xmm3,[rdi+8]		; a01.i a01.r			<(aa)[3].e[0][1]>
	mulps           xmm0,xmm3
	movaps		xmm1,xmm9
	movlps		xmm4,[rdi+24]		;             a10.i a10.r	<(aa)[3].e[1][0]>
	movhps		xmm4,[rdi+32]		; a11.i a11.r			<(aa)[3].e[1][1]>
	mulps           xmm1,xmm4
	movaps		xmm2,xmm10
	movlps		xmm5,[rdi+48]		;             a20.i a20.r	<(aa)[3].e[2][0]>
	movhps		xmm5,[rdi+56]		; a21.i a21.r			<(aa)[3].e[2][1]>
        mulps           xmm2,xmm5
        addps           xmm0,xmm1
        addps           xmm0,xmm2

        ; imaginary part of b
        mulps           xmm3,xmm11
        mulps           xmm4,xmm12
        mulps           xmm5,xmm13
        addps           xmm3,xmm4
        addps           xmm3,xmm5
        xorps           xmm0,[negate]		; <_sse_sgn24>
        shufps          xmm3,xmm3,0xb1
	addps           xmm0,xmm3
	movlps		[r9],xmm0		; <(cc3)->c[0]>
	movhps		[r9+8],xmm0		; <(cc3)->c[1]>


	; third row
	movlps		xmm1,[rdi+16]		;             a02.i a02.r	<(aa)[3].e[0][2]>
	unpcklps	xmm1,xmm1		; a02.i a02.i a02.r a02.r
	movaps		xmm6,xmm1
	movlps		xmm2,[rdi+40]		;             a12.i a12.r	<(aa)[3].e[1][2]>
	unpcklps	xmm2,xmm2		; a12.i a12.i a12.r a12.r
	movlhps		xmm6,xmm2		; a12.r a12.r a02.r a02.r
	movhlps		xmm2,xmm1		; a12.i a12.i a02.i a02.i

	mulps		xmm6,xmm14
	mulps		xmm2,xmm14
	shufps		xmm2,xmm2,0xb1
	xorps		xmm2,[negate]		; <_sse_sgn24>
	addps		xmm6,xmm2
	; xmm6 now contains a02*b0 and a12*b1

	movlps		xmm7,[rdi+64]		;             a22.i a22.r	<(aa)[3].e[2][2]>
	unpcklps	xmm7,xmm7		; a22.i a22.i a22.r a22.r
	mulps		xmm7,xmm15
	shufps		xmm7,xmm7,0xb4
	xorps		xmm7,[neg4]		; <_sse_sgn4>
	addps		xmm6,xmm7

	; need to add high bits to low bits in xmm6
	movhlps		xmm5,xmm6
	addps		xmm6,xmm5
	movlps		[r9+16],xmm6		; <(cc3)->c[2]>

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
	
