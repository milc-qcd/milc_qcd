;
; mult_su3_mat_vec_sum_4dir( su3_matrix *a[4], su3_vector *b0,
;                            su3_vector *b1, su3_vector *b2,
;                            su3_vector *b3, su3_vector *c)
;
; Multiply an array of 4 su3 matrices, a[4], by an array of vectors, b[4],
; summing the result into an su3 vector c.
;
; Steve Whalen
; Cray, Inc.

	bits		64
	
global mult_su3_mat_vec_sum_4dir
mult_su3_mat_vec_sum_4dir:

	;; first product: c <- (A[0])b0

	;; rows 1 and 2

        ; real part of b
	movlps		xmm14,[rsi]		;              b0.i  b0.r	<(bb0)->c[0]>
        movaps          xmm0,xmm14
	unpcklps	xmm0,xmm0		;  b0.i  b0.i  b0.r  b0.r
	movhlps		xmm11,xmm0		;              b0.i  b0.i
	movlhps		xmm0,xmm0		;  b0.r  b0.r  b0.r  b0.r
	movhps		xmm14,[rsi+8]		;  b1.i  b1.r			<(bb0)->c[1]>
        movaps          xmm6,xmm14
	unpckhps	xmm6,xmm6		;  b1.i  b1.i  b1.r  b1.r
	movhlps		xmm12,xmm6		;              b1.i  b1.i
	movlhps		xmm6,xmm6		;  b1.r  b1.r  b1.r  b1.r
	movhps		xmm15,[rsi+16]		;  b2.i  b2.r			<(bb0)->c[2]>
        movaps          xmm7,xmm15
	unpckhps	xmm7,xmm7		;  b2.i  b2.i  b2.r  b2.r
	movhlps		xmm13,xmm7		;              b2.i  b2.i
	movlhps		xmm7,xmm7		;  b2.r  b2.r  b2.r  b2.r

	movlps		xmm1,[rdi]		;             a00.i a00.r	<(aa)[0].e[0][0]>
	movhps		xmm1,[rdi+24]		; a10.i a10.r			<(aa)[0].e[1][0]>
	mulps           xmm0,xmm1
	movlps		xmm9,[rdi+8]		;             a01.i a01.r	<(aa)[0].e[0][1]>
	movhps		xmm9,[rdi+32]		; a11.i a11.r			<(aa)[0].e[1][1]>
	mulps           xmm6,xmm9
	movlps		xmm10,[rdi+16]		;             a02.i a02.r	<(aa)[0].e[0][2]>
	movhps		xmm10,[rdi+40]		; a12.i a12.r			<(aa)[0].e[1][2]>
        mulps           xmm7,xmm10
        addps           xmm0,xmm6
        addps           xmm0,xmm7
	; xmm0 contains part of the answer's first two components

        ; imaginary part of b
        movlhps		xmm11,xmm11		;  b0.i  b0.i  b0.i  b0.i
        movlhps		xmm12,xmm12		;  b1.i  b1.i  b1.i  b1.i
        movlhps		xmm13,xmm13		;  b2.i  b2.i  b2.i  b2.i
        mulps           xmm1,xmm11
        mulps           xmm9,xmm12
        mulps           xmm10,xmm13
        addps           xmm1,xmm9
        addps           xmm1,xmm10
	; xmm1 contains those terms of the first two components that
	; are absent from xmm0. We still need to negate and shuffle xmm1
	; before adding it to xmm0. Rather than doing this now, we'll
	; continue to accumulate results in xmm0 and xmm1, and combine
	; them at the end.

	; third row
	movlps		xmm6,[rdi+48]		;             a20.i a20.r	<(aa)[0].e[2][0]>
	unpcklps	xmm6,xmm6		; a20.i a20.i a20.r a20.r
	movaps		xmm2,xmm6
	movlps		xmm3,[rdi+56]		;             a21.i a21.r	<(aa)[0].e[2][1]>
	unpcklps	xmm3,xmm3		; a21.i a21.i a21.r a21.r
	movlhps		xmm2,xmm3		; a21.r a21.r a20.r a20.r
	movhlps		xmm3,xmm6		; a21.i a21.i a20.i a20.i

	mulps		xmm2,xmm14
	mulps		xmm3,xmm14
	; xmm2 and xmm3 contain term of a20*b0 and a21*b1; xmm3 is waiting
	; for negation and shuffling

	movlps		xmm4,[rdi+64]		;             a22.i a22.r	<(aa)[0].e[2][2]>
	unpcklps	xmm4,xmm4		; a22.i a22.i a22.r a22.r
	movhlps		xmm15,xmm15		;  b2.i  b2.r  b2.i  b2.r
	mulps		xmm4,xmm15
	; xmm4 will also be used for accumulation

	; *******************************************************************	
	; *******************************************************************	

	;; second product: c <- c + (A[1])b1

	add		rdi,72

	;; rows 1 and 2

        ; real part of b
	movlps		xmm14,[rdx]		;              b0.i  b0.r	<(bb1)->c[0]>
        movaps          xmm5,xmm14
	unpcklps	xmm5,xmm5		;  b0.i  b0.i  b0.r  b0.r
	movhlps		xmm11,xmm5		;              b0.i  b0.i
	movlhps		xmm5,xmm5		;  b0.r  b0.r  b0.r  b0.r
	movhps		xmm14,[rdx+8]		;  b1.i  b1.r			<(bb1)->c[1]>
        movaps          xmm6,xmm14
	unpckhps	xmm6,xmm6		;  b1.i  b1.i  b1.r  b1.r
	movhlps		xmm12,xmm6		;              b1.i  b1.i
	movlhps		xmm6,xmm6		;  b1.r  b1.r  b1.r  b1.r
	movhps		xmm15,[rdx+16]		;  b2.i  b2.r			<(bb1)->c[2]>
        movaps          xmm7,xmm15
	unpckhps	xmm7,xmm7		;  b2.i  b2.i  b2.r  b2.r
	movhlps		xmm13,xmm7		;              b2.i  b2.i
	movlhps		xmm7,xmm7		;  b2.r  b2.r  b2.r  b2.r

	movlps		xmm8,[rdi]		;             a00.i a00.r	<(aa)[1].e[0][0]>
	movhps		xmm8,[rdi+24]		; a10.i a10.r			<(aa)[1].e[1][0]>
	mulps           xmm5,xmm8
	movlps		xmm9,[rdi+8]		;             a01.i a01.r	<(aa)[1].e[0][1]>
	movhps		xmm9,[rdi+32]		; a11.i a11.r			<(aa)[1].e[1][1]>
	mulps           xmm6,xmm9
	movlps		xmm10,[rdi+16]		;             a02.i a02.r	<(aa)[1].e[0][2]>
	movhps		xmm10,[rdi+40]		; a12.i a12.r			<(aa)[1].e[1][2]>
        mulps           xmm7,xmm10
        addps           xmm5,xmm6
        addps           xmm5,xmm7
	addps		xmm0,xmm5

        ; imaginary part of b
        movlhps		xmm11,xmm11		;  b0.i  b0.i  b0.i  b0.i
        movlhps		xmm12,xmm12		;  b1.i  b1.i  b1.i  b1.i
        movlhps		xmm13,xmm13		;  b2.i  b2.i  b2.i  b2.i
        mulps           xmm8,xmm11
        mulps           xmm9,xmm12
        mulps           xmm10,xmm13
        addps           xmm8,xmm9
        addps           xmm8,xmm10
	addps		xmm1,xmm8

	; third row
	movlps		xmm6,[rdi+48]		;             a20.i a20.r	<(aa)[1].e[2][0]>
	unpcklps	xmm6,xmm6		; a20.i a20.i a20.r a20.r
	movaps		xmm7,xmm6
	movlps		xmm5,[rdi+56]		;             a21.i a21.r	<(aa)[1].e[2][1]>
	unpcklps	xmm5,xmm5		; a21.i a21.i a21.r a21.r
	movlhps		xmm7,xmm5		; a21.r a21.r a20.r a20.r
	movhlps		xmm5,xmm6		; a21.i a21.i a20.i a20.i

	mulps		xmm7,xmm14
	mulps		xmm5,xmm14
	addps		xmm2,xmm7
	addps		xmm3,xmm5

	movlps		xmm11,[rdi+64]		;             a22.i a22.r	<(aa)[1].e[2][2]>
	unpcklps	xmm11,xmm11		; a22.i a22.i a22.r a22.r
	movhlps		xmm15,xmm15		;  b2.i  b2.r  b2.i  b2.r
	mulps		xmm11,xmm15
	addps		xmm4,xmm11

	; *******************************************************************	
	; *******************************************************************	

	;; third product: c <- c + (A[2])b2

	add		rdi,72

	;; rows 1 and 2

        ; real part of b
	movlps		xmm14,[rcx]		;              b0.i  b0.r	<(bb2)->c[0]>
        movaps          xmm5,xmm14
	unpcklps	xmm5,xmm5		;  b0.i  b0.i  b0.r  b0.r
	movhlps		xmm11,xmm5		;              b0.i  b0.i
	movlhps		xmm5,xmm5		;  b0.r  b0.r  b0.r  b0.r
	movhps		xmm14,[rcx+8]		;  b1.i  b1.r			<(bb2)->c[1]>
        movaps          xmm6,xmm14
	unpckhps	xmm6,xmm6		;  b1.i  b1.i  b1.r  b1.r
	movhlps		xmm12,xmm6		;              b1.i  b1.i
	movlhps		xmm6,xmm6		;  b1.r  b1.r  b1.r  b1.r
	movhps		xmm15,[rcx+16]		;  b2.i  b2.r			<(bb2)->c[2]>
        movaps          xmm7,xmm15
	unpckhps	xmm7,xmm7		;  b2.i  b2.i  b2.r  b2.r
	movhlps		xmm13,xmm7		;              b2.i  b2.i
	movlhps		xmm7,xmm7		;  b2.r  b2.r  b2.r  b2.r

	movlps		xmm8,[rdi]		;             a00.i a00.r	<(aa)[2].e[0][0]>
	movhps		xmm8,[rdi+24]		; a10.i a10.r			<(aa)[2].e[1][0]>
	mulps           xmm5,xmm8
	movlps		xmm9,[rdi+8]		;             a01.i a01.r	<(aa)[2].e[0][1]>
	movhps		xmm9,[rdi+32]		; a11.i a11.r			<(aa)[2].e[1][1]>
	mulps           xmm6,xmm9
	movlps		xmm10,[rdi+16]		;             a02.i a02.r	<(aa)[2].e[0][2]>
	movhps		xmm10,[rdi+40]		; a12.i a12.r			<(aa)[2].e[1][2]>
        mulps           xmm7,xmm10
        addps           xmm5,xmm6
        addps           xmm5,xmm7
	addps		xmm0,xmm5

        ; imaginary part of b
        movlhps		xmm11,xmm11		;  b0.i  b0.i  b0.i  b0.i
        movlhps		xmm12,xmm12		;  b1.i  b1.i  b1.i  b1.i
        movlhps		xmm13,xmm13		;  b2.i  b2.i  b2.i  b2.i
        mulps           xmm8,xmm11
        mulps           xmm9,xmm12
        mulps           xmm10,xmm13
        addps           xmm8,xmm9
        addps           xmm8,xmm10
	addps		xmm1,xmm8

	; third row
	movlps		xmm6,[rdi+48]		;             a20.i a20.r	<(aa)[2].e[2][0]>
	unpcklps	xmm6,xmm6		; a20.i a20.i a20.r a20.r
	movaps		xmm7,xmm6
	movlps		xmm5,[rdi+56]		;             a21.i a21.r	<(aa)[2].e[2][1]>
	unpcklps	xmm5,xmm5		; a21.i a21.i a21.r a21.r
	movlhps		xmm7,xmm5		; a21.r a21.r a20.r a20.r
	movhlps		xmm5,xmm6		; a21.i a21.i a20.i a20.i

	mulps		xmm7,xmm14
	mulps		xmm5,xmm14
	addps		xmm2,xmm7
	addps		xmm3,xmm5

	movlps		xmm11,[rdi+64]		;             a22.i a22.r	<(aa)[2].e[2][2]>
	unpcklps	xmm11,xmm11		; a22.i a22.i a22.r a22.r
	movhlps		xmm15,xmm15		;  b2.i  b2.r  b2.i  b2.r
	mulps		xmm11,xmm15
	addps		xmm4,xmm11

	; *******************************************************************	
	; *******************************************************************	

	;; fourth product: c <- c + (A[3])b3

	add		rdi,72

	;; rows 1 and 2

        ; real part of b
	movlps		xmm14,[r8]		;              b0.i  b0.r	<(bb3)->c[0]>
        movaps          xmm5,xmm14
	unpcklps	xmm5,xmm5		;  b0.i  b0.i  b0.r  b0.r
	movhlps		xmm11,xmm5		;              b0.i  b0.i
	movlhps		xmm5,xmm5		;  b0.r  b0.r  b0.r  b0.r
	movhps		xmm14,[r8+8]		;  b1.i  b1.r			<(bb3)->c[1]>
        movaps          xmm6,xmm14
	unpckhps	xmm6,xmm6		;  b1.i  b1.i  b1.r  b1.r
	movhlps		xmm12,xmm6		;              b1.i  b1.i
	movlhps		xmm6,xmm6		;  b1.r  b1.r  b1.r  b1.r
	movhps		xmm15,[r8+16]		;  b2.i  b2.r			<(bb3)->c[2]>
        movaps          xmm7,xmm15
	unpckhps	xmm7,xmm7		;  b2.i  b2.i  b2.r  b2.r
	movhlps		xmm13,xmm7		;              b2.i  b2.i
	movlhps		xmm7,xmm7		;  b2.r  b2.r  b2.r  b2.r

	movlps		xmm8,[rdi]		;             a00.i a00.r	<(aa)[3].e[0][0]>
	movhps		xmm8,[rdi+24]		; a10.i a10.r			<(aa)[3].e[1][0]>
	mulps           xmm5,xmm8
	movlps		xmm9,[rdi+8]		;             a01.i a01.r	<(aa)[3].e[0][1]>
	movhps		xmm9,[rdi+32]		; a11.i a11.r			<(aa)[3].e[1][1]>
	mulps           xmm6,xmm9
	movlps		xmm10,[rdi+16]		;             a02.i a02.r	<(aa)[3].e[0][2]>
	movhps		xmm10,[rdi+40]		; a12.i a12.r			<(aa)[3].e[1][2]>
        mulps           xmm7,xmm10
        addps           xmm5,xmm6
        addps           xmm5,xmm7
	addps		xmm0,xmm5

        ; imaginary part of b
        movlhps		xmm11,xmm11		;  b0.i  b0.i  b0.i  b0.i
        movlhps		xmm12,xmm12		;  b1.i  b1.i  b1.i  b1.i
        movlhps		xmm13,xmm13		;  b2.i  b2.i  b2.i  b2.i
        mulps           xmm8,xmm11
        mulps           xmm9,xmm12
        mulps           xmm10,xmm13
        addps           xmm8,xmm9
        addps           xmm8,xmm10
	addps		xmm1,xmm8

	; finalize first two components of answer
	xorps		xmm1,[negate]		; <_sse_sgn24>
	shufps		xmm1,xmm1,0xb1
	addps		xmm0,xmm1
	movlps		[r9],xmm0		; <(cc)->c[0]>
	movhps		[r9+8],xmm0		; <(cc)->c[1]>

	; third row
	movlps		xmm6,[rdi+48]		;             a20.i a20.r	<(aa)[3].e[2][0]>
	unpcklps	xmm6,xmm6		; a20.i a20.i a20.r a20.r
	movaps		xmm7,xmm6
	movlps		xmm5,[rdi+56]		;             a21.i a21.r	<(aa)[3].e[2][1]>
	unpcklps	xmm5,xmm5		; a21.i a21.i a21.r a21.r
	movlhps		xmm7,xmm5		; a21.r a21.r a20.r a20.r
	movhlps		xmm5,xmm6		; a21.i a21.i a20.i a20.i

	mulps		xmm7,xmm14
	mulps		xmm5,xmm14
	addps		xmm2,xmm7
	addps		xmm3,xmm5
	xorps		xmm3,[negate]		; <_sse_sgn24>
	shufps		xmm3,xmm3,0xb1
	addps		xmm2,xmm3

	movlps		xmm11,[rdi+64]		;             a22.i a22.r	<(aa)[3].e[2][2]>
	unpcklps	xmm11,xmm11		; a22.i a22.i a22.r a22.r
	movhlps		xmm15,xmm15		;  b2.i  b2.r  b2.i  b2.r
	mulps		xmm11,xmm15
	addps		xmm4,xmm11
	xorps		xmm4,[neg4]		; <_sse_sgn4>
	shufps		xmm4,xmm4,0xb4
	addps		xmm2,xmm4

	; need to add high bits to low bits in xmm2
	movhlps		xmm10,xmm2
	addps		xmm2,xmm10
	movlps		[r9+16],xmm2		; <(cc)->c[2]>

	; *******************************************************************	
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
	
