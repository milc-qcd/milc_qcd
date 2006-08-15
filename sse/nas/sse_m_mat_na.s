;
; mult_su3_na( su3_matrix *a, su3_matrix *b, su3_matrix *c )
;
; Opteron SSE code for 3x3 CGEMM
; C <- AB
; 
; Steve Whalen
; Cray, Inc.

	bits		64

global mult_su3_na
mult_su3_na:

	;; first row, columns 1 and 2

	; real part of a
	movlps		xmm7,[rdi]		;               a00.i  a00.r	<(aa)->e[0][0]>
	movaps		xmm0,xmm7
	unpcklps	xmm0,xmm0		; a00.i  a00.i  a00.r  a00.r
	movhlps		xmm3,xmm0		;               a00.i  a00.i
	movlhps		xmm0,xmm0		; a00.r  a00.r  a00.r  a00.r
	movhps		xmm7,[rdi+8]		; a01.i  a01.r			<(aa)->e[0][1]>
	movaps		xmm1,xmm7
	unpckhps	xmm1,xmm1		; a01.i  a01.i  a01.r  a01.r
	movhlps		xmm4,xmm1		;               a01.i  a01.i
	movlhps		xmm1,xmm1		; a01.r  a01.r  a01.r  a01.r
	movlps		xmm10,[rdi+16]		;               a02.i  a02.r	<(aa)->e[0][2]>
	movaps		xmm2,xmm10
	unpcklps	xmm2,xmm2		; a02.i  a02.i  a02.r  a02.r
	movhlps		xmm5,xmm2		;               a02.i  a02.i
	movlhps		xmm2,xmm2		; a02.r  a02.r  a02.r  a02.r

	movlps		xmm12,[rsi]		;               b00.i  b00.r	<(bb)->e[0][0]>
	movhps		xmm12,[rsi+24]		; b10.i  b10.r			<(bb)->e[1][0]>
	mulps		xmm0,xmm12
	movlps		xmm13,[rsi+8]		;               b01.i  b01.r	<(bb)->e[0][1]>
	movhps		xmm13,[rsi+32]		; b11.i  b11.r			<(bb)->e[1][1]>
	mulps		xmm1,xmm13

	addps		xmm0,xmm1
	movlps		xmm14,[rsi+16]		;               b02.i  b02.r	<(bb)->e[0][2]>
	movhps		xmm14,[rsi+40]		; b12.i  b12.r			<(bb)->e[1][2]>
	mulps		xmm2,xmm14
	addps		xmm0,xmm2

	; imaginary part of a
	movlhps		xmm3,xmm3		; a00.i  a00.i  a00.i  a00.i
	movlhps		xmm4,xmm4		; a01.i  a01.i  a01.i  a01.i
	movlhps		xmm5,xmm5		; a02.i  a02.i  a02.i  a02.i
	mulps		xmm3,xmm12
	mulps		xmm4,xmm13
	mulps		xmm5,xmm14
	addps		xmm3,xmm4
	addps		xmm3,xmm5
	xorps		xmm0,[negate]		; <_sse_sgn24>
	shufps		xmm3,xmm3,0xb1
	addps		xmm0,xmm3
	movlps		[rdx],xmm0		; <(cc)->e[0][0]>
	movhps		[rdx+8],xmm0		; <(cc)->e[0][1]>

        ;; second row, columns 1 and 2

        ; real part of a
	movlps		xmm8,[rdi+24]		;               a10.i  a10.r	<(aa)->e[1][0]>
        movaps          xmm1,xmm8
	unpcklps	xmm1,xmm1		; a10.i  a10.i  a10.r  a10.r
	movhlps		xmm4,xmm1		;               a10.i  a10.i
	movlhps		xmm1,xmm1		; a10.r  a10.r  a10.r  a10.r
	movhps		xmm8,[rdi+32]		; a11.i  a11.r			<(aa)->e[1][1]>
        movaps          xmm2,xmm8
	unpckhps	xmm2,xmm2		; a11.i  a11.i  a11.r  a11.r
	movhlps		xmm5,xmm2		;               a11.i  a11.i
	movlhps		xmm2,xmm2		; a11.r  a11.r  a11.r  a11.r
	movhps		xmm10,[rdi+40]		; a12.i  a12.r			<(aa)->e[1][2]>
	movaps		xmm3,xmm10
	unpckhps	xmm3,xmm3		; a12.i  a12.i  a12.r  a12.r
	movhlps		xmm6,xmm3		;               a12.i  a12.i
	movlhps		xmm3,xmm3		; a12.r  a12.r  a12.r  a12.r

        mulps           xmm1,xmm12
        mulps           xmm2,xmm13
	mulps		xmm3,xmm14
        addps           xmm1,xmm2
        addps           xmm1,xmm3

        ; imaginary part of a
        movlhps		xmm4,xmm4		; a10.i  a10.i  a10.i  a10.i
        movlhps		xmm5,xmm5		; a11.i  a11.i  a11.i  a11.i
	movlhps		xmm6,xmm6		; a12.i  a12.i  a12.i  a12.i
        mulps           xmm4,xmm12
        mulps           xmm5,xmm13
	mulps		xmm6,xmm14
        addps           xmm4,xmm5
        addps           xmm4,xmm6
        xorps           xmm1,[negate]		; <_sse_sgn24>
        shufps          xmm4,xmm4,0xb1
        addps           xmm1,xmm4
	movlps		[rdx+24],xmm1		; <(cc)->e[1][0]>
	movhps		[rdx+32],xmm1		; <(cc)->e[1][1]>

	;; third column, rows 1 and 2
	;; think of this as in the third row of (B^T)(A^T)

	; transpose a
	movaps		xmm5,xmm7		;   *      *    a00.i  a00.r
	movlhps		xmm5,xmm8		; a10.i  a10.r
	movhlps		xmm8,xmm7		; a11.i  a11.r  a01.i  a01.r
  	; xmm10 already contains a02 and a12
	; xmm7 can now be modified

        ; real part of b
	movlps		xmm15,[rsi+48]		;               b20.i  b20.r	<(bb)->e[2][0]>
        movaps          xmm2,xmm15
	unpcklps	xmm2,xmm2		; b20.i  b20.i  b20.r  b20.r
	movhlps		xmm6,xmm2		;               b20.i  b20.i
	movlhps		xmm2,xmm2		; b20.r  b20.r  b20.r  b20.r
	movhps		xmm15,[rsi+56]		; b21.i  b21.r			<(bb)->e[2][1]>
        movaps          xmm3,xmm15
	unpckhps	xmm3,xmm3		; b21.i  b21.i  b21.r  b21.r
	movhlps		xmm7,xmm3		;               b21.i  b21.i
	movlhps		xmm3,xmm3		; b21.r  b21.r  b21.r  b21.r
	movhps		xmm11,[rsi+64]		; b22.i  b22.r			<(bb)->e[2][2]>
        movaps          xmm4,xmm11
	unpckhps	xmm4,xmm4		; b22.i  b22.i  b22.r  b22.r
	movhlps		xmm9,xmm4		;               b22.i  b22.i
	movlhps		xmm4,xmm4		; b22.r  b22.r  b22.r  b22.r

        mulps           xmm2,xmm5
        mulps           xmm3,xmm8
        mulps           xmm4,xmm10
        addps           xmm2,xmm3
        addps           xmm2,xmm4

        ; imaginary part of b
        movlhps		xmm6,xmm6		; b20.i  b20.i  b20.i  b20.i
        movlhps		xmm7,xmm7		; b21.i  b21.i  b21.i  b21.i
        movlhps		xmm9,xmm9		; b22.i  b22.i  b22.i  b22.i
        mulps           xmm6,xmm5
        mulps           xmm7,xmm8
	; xmm8 can now be modified
        mulps           xmm9,xmm10
	; xmm10 can now be modified
        addps           xmm6,xmm7
        addps           xmm6,xmm9
        shufps          xmm6,xmm6,0xb1
        xorps           xmm6,[negate]		; <_sse_sgn24>
	addps           xmm2,xmm6
	movlps		[rdx+16],xmm2		; <(cc)->e[0][2]>
	movhps		[rdx+40],xmm2		; <(cc)->e[1][2]>


	; recap: we are done with registers xmm7, xmm8, xmm10
	; xmm0, xmm1, xmm2 contain the first two rows of the result:
	;
	;      0    1    2
	;   [ --------  --- ]
	; 0 [ | xmm0 |  |x| ]
	;   [ --------  |m| ]
	;   [ --------  |m| ]
	; 1 [ | xmm1 |  |2| ]
	;   [ --------  --- ]
	
        ;; third row, columns 1 and 2

        ; real part of a
	movlps		xmm3,[rdi+48]		;               a20.i  a20.r	<(aa)->e[2][0]>
	unpcklps	xmm3,xmm3		; a20.i  a20.i  a20.r  a20.r
	movhlps		xmm6,xmm3		;               a20.i  a20.i
	movlhps		xmm3,xmm3		; a20.r  a20.r  a20.r  a20.r
	movhps		xmm4,[rdi+56]		; a21.i  a21.r			<(aa)->e[2][1]>
	unpckhps	xmm4,xmm4		; a21.i  a21.i  a21.r  a21.r
	movhlps		xmm7,xmm4		;               a21.i  a21.i
	movlhps		xmm4,xmm4		; a21.r  a21.r  a21.r  a21.r
	movlps		xmm11,[rdi+64]		;               a22.i  a22.r	<(aa)->e[2][2]>
        movaps          xmm5,xmm11
	unpcklps	xmm5,xmm5		; a22.i  a22.i  a22.r  a22.r
	movhlps		xmm10,xmm5		;               a22.i  a22.i
	movlhps		xmm5,xmm5		; a22.r  a22.r  a22.r  a22.r

	; save some values for later
	movaps		xmm9,xmm3
	movlhps		xmm9,xmm4

        mulps           xmm3,xmm12
	mulps           xmm4,xmm13
        mulps           xmm5,xmm14
        addps           xmm3,xmm4
        addps           xmm3,xmm5

        ; imaginary part of a
        movlhps		xmm6,xmm6		; a20.i  a20.i  a20.i  a20.i
        movlhps		xmm7,xmm7		; a21.i  a21.i  a21.i  a21.i
        movlhps		xmm10,xmm10		; a22.i  a22.i  a22.i  a22.i
        mulps           xmm12,xmm6
        mulps           xmm13,xmm7
        mulps           xmm14,xmm10
        addps           xmm12,xmm13
        addps           xmm12,xmm14
        xorps           xmm3,[negate]		; <_sse_sgn24>
        shufps          xmm12,xmm12,0xb1
        addps           xmm3,xmm12
	movlps		[rdx+48],xmm3		; <(cc)->e[2][0]>
	movhps		[rdx+56],xmm3		; <(cc)->e[2][1]>

	; We are now done with registers xmm8, xmm10, xmm12, xmm13, xmm14.
	; xmm0, xmm1, xmm2, xmm3 contain all but the last component of the result:
	;
	;      0    1    2
	;   [ --------  --- ]
	; 0 [ | xmm0 |  |x| ]
	;   [ --------  |m| ]
	;   [ --------  |m| ]
	; 1 [ | xmm1 |  |2| ]
	;   [ --------  --- ]
	;   [ --------      ]
	; 2 [ | xmm3 |      ]
	;   [ --------      ]
	
	; third row, third column
	movlhps		xmm6,xmm7		; a21.i  a21.i  a20.i  a20.i
	mulps		xmm9,xmm15
	mulps		xmm6,xmm15
	xorps		xmm9,[negate]		; <_sse_sgn24>
	shufps		xmm6,xmm6,0xb1
	addps		xmm9,xmm6
	; xmm9 now contains a20*b02 and a21*b12

	movaps		xmm8,xmm11		; b22.i  b22.r  a22.i  a22.r
	unpcklps	xmm8,xmm8		; a22.i  a22.i  a22.r  a22.r
	movhlps		xmm11,xmm11		; b22.i  b22.r  b22.i  b22.r
	mulps		xmm11,xmm8
	xorps		xmm11,[neg2]		; <_sse_sgn2>
	shufps		xmm11,xmm11,0xb4
	addps		xmm9,xmm11

	; need to add high bits to low bits in xmm9
	movhlps		xmm5,xmm9
	addps		xmm9,xmm5
	movlps		[rdx+64],xmm9		; <(cc)->e[2][2]>

	; Done!

	; *******************************************************************	

here:	ret
	
	align		16
negate:	dd		0x00000000
	dd		0x80000000
	dd		0x00000000
	dd		0x80000000

	align		16
neg2:	dd		0x00000000
	dd		0x80000000
	dd		0x00000000
	dd		0x00000000
