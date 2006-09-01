;
; mult_adj_su3_mat_hwvec( su3_matrix *a, half_wilson_vector *b, *c)
;
; Steve Whalen
; Cray, Inc.

        bits            64

global mult_adj_su3_mat_hwvec
mult_adj_su3_mat_hwvec:

        ;; perform a 3x3-times-3x2 matrix multiply

	;; first row

	; real part of a
	movlps		xmm7,[rdi]		;               a00.i  a00.r	<(aa)->e[0][0]>
	movaps		xmm0,xmm7
	unpcklps	xmm0,xmm0		; a00.i  a00.i  a00.r  a00.r
	movhlps		xmm3,xmm0		;               a00.i  a00.i
	movlhps		xmm0,xmm0		; a00.r  a00.r  a00.r  a00.r
	movhps		xmm7,[rdi+24]		; a10.i  a10.r			<(aa)->e[1][0]>
	movaps		xmm1,xmm7
	unpckhps	xmm1,xmm1		; a10.i  a10.i  a10.r  a10.r
	movhlps		xmm4,xmm1		;               a10.i  a10.i
	movlhps		xmm1,xmm1		; a10.r  a10.r  a10.r  a10.r
	movlps		xmm10,[rdi+48]		;               a20.i  a20.r	<(aa)->e[2][0]>
	movaps		xmm2,xmm10
	unpcklps	xmm2,xmm2		; a20.i  a20.i  a20.r  a20.r
	movhlps		xmm5,xmm2		;               a20.i  a20.i
	movlhps		xmm2,xmm2		; a20.r  a20.r  a20.r  a20.r

	movlps		xmm12,[rsi]		;               b00.i  b00.r	<(bb)->h[0].c[0]>
	movhps		xmm12,[rsi+24]		; b01.i  b01.r			<(bb)->h[1].c[0]>
	mulps		xmm0,xmm12
	movlps		xmm13,[rsi+8]		;               b10.i  b10.r	<(bb)->h[0].c[1]>
	movhps		xmm13,[rsi+32]		; b11.i  b11.r			<(bb)->h[1].c[1]>
	mulps		xmm1,xmm13

	addps		xmm0,xmm1
	movlps		xmm14,[rsi+16]		;               b20.i  b20.r	<(bb)->h[0].c[2]>
	movhps		xmm14,[rsi+40]		; b21.i  b21.r			<(bb)->h[1].c[2]>
	mulps		xmm2,xmm14
	addps		xmm0,xmm2

	; imaginary part of a
	movlhps		xmm3,xmm3		; a00.i  a00.i  a00.i  a00.i
	movlhps		xmm4,xmm4		; a10.i  a10.i  a10.i  a10.i
	movlhps		xmm5,xmm5		; a20.i  a20.i  a20.i  a20.i
	mulps		xmm3,xmm12
	mulps		xmm4,xmm13
	mulps		xmm5,xmm14
	addps		xmm3,xmm4
	addps		xmm3,xmm5
	shufps		xmm3,xmm3,0xb1
	xorps		xmm3,[negate]		; <_sse_sgn24>
	addps		xmm0,xmm3
	movlps		[rdx],xmm0		; <(cc)->h[0].c[0]>
	movhps		[rdx+24],xmm0		; <(cc)->h[1].c[0]>

        ;; second row

        ; real part of a
	movlps		xmm8,[rdi+8]		;               a01.i  a01.r	<(aa)->e[0][1]>
        movaps          xmm1,xmm8
	unpcklps	xmm1,xmm1		; a01.i  a01.i  a01.r  a01.r
	movhlps		xmm4,xmm1		;               a01.i  a01.i
	movlhps		xmm1,xmm1		; a01.r  a01.r  a01.r  a01.r
	movhps		xmm8,[rdi+32]		; a11.i  a11.r			<(aa)->e[1][1]>
        movaps          xmm2,xmm8
	unpckhps	xmm2,xmm2		; a11.i  a11.i  a11.r  a11.r
	movhlps		xmm5,xmm2		;               a11.i  a11.i
	movlhps		xmm2,xmm2		; a11.r  a11.r  a11.r  a11.r
	movhps		xmm10,[rdi+56]		; a21.i  a21.r			<(aa)->e[2][1]>
	movaps		xmm3,xmm10
	unpckhps	xmm3,xmm3		; a21.i  a21.i  a21.r  a21.r
	movhlps		xmm6,xmm3		;               a21.i  a21.i
	movlhps		xmm3,xmm3		; a21.r  a21.r  a21.r  a21.r

        mulps           xmm1,xmm12
        mulps           xmm2,xmm13
	mulps		xmm3,xmm14
        addps           xmm1,xmm2
        addps           xmm1,xmm3

        ; imaginary part of a
        movlhps		xmm4,xmm4		; a01.i  a01.i  a01.i  a01.i
        movlhps		xmm5,xmm5		; a11.i  a11.i  a11.i  a11.i
	movlhps		xmm6,xmm6		; a21.i  a21.i  a21.i  a21.i
        mulps           xmm4,xmm12
        mulps           xmm5,xmm13
	mulps		xmm6,xmm14
        addps           xmm4,xmm5
        addps           xmm4,xmm6
        shufps          xmm4,xmm4,0xb1
        xorps           xmm4,[negate]		; <_sse_sgn24>
        addps           xmm1,xmm4
	movlps		[rdx+8],xmm1		; <(cc)->h[0].c[1]>
	movhps		[rdx+32],xmm1		; <(cc)->h[1].c[1]>

        ;; third row

        ; real part of a
	movlps		xmm3,[rdi+16]		;               a02.i  a02.r	<(aa)->e[0][2]>
	unpcklps	xmm3,xmm3		; a02.i  a02.i  a02.r  a02.r
	movhlps		xmm6,xmm3		;               a02.i  a02.i
	movlhps		xmm3,xmm3		; a02.r  a02.r  a02.r  a02.r
	movhps		xmm4,[rdi+40]		; a12.i  a12.r			<(aa)->e[1][2]>
	unpckhps	xmm4,xmm4		; a12.i  a12.i  a12.r  a12.r
	movhlps		xmm7,xmm4		;               a12.i  a12.i
	movlhps		xmm4,xmm4		; a12.r  a12.r  a12.r  a12.r
	movlps		xmm11,[rdi+64]		;               a22.i  a22.r	<(aa)->e[2][2]>
        movaps          xmm5,xmm11
	unpcklps	xmm5,xmm5		; a22.i  a22.i  a22.r  a22.r
	movhlps		xmm10,xmm5		;               a22.i  a22.i
	movlhps		xmm5,xmm5		; a22.r  a22.r  a22.r  a22.r

        mulps           xmm3,xmm12
	mulps           xmm4,xmm13
        mulps           xmm5,xmm14
        addps           xmm3,xmm4
        addps           xmm3,xmm5

        ; imaginary part of a
        movlhps		xmm6,xmm6		; a02.i  a02.i  a02.i  a02.i
        movlhps		xmm7,xmm7		; a12.i  a12.i  a12.i  a12.i
        movlhps		xmm10,xmm10		; a22.i  a22.i  a22.i  a22.i
        mulps           xmm12,xmm6
        mulps           xmm13,xmm7
        mulps           xmm14,xmm10
        addps           xmm12,xmm13
        addps           xmm12,xmm14
        shufps          xmm12,xmm12,0xb1
        xorps           xmm12,[negate]		; <_sse_sgn24>
        addps           xmm3,xmm12
	movlps		[rdx+16],xmm3		; <(cc)->h[0].c[2]>
	movhps		[rdx+40],xmm3		; <(cc)->h[1].c[2]>

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
