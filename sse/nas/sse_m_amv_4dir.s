;
; mult_adj_su3_mat_vec_4dir( su3_matrix *a[4], su3_vector *b, su3_vector *c[4])
;
; Multiply the adjoint of each of four input matrices by an input vector, 
; storing the resulting vectors in an array.
;   

global mult_adj_su3_mat_vec_4dir
mult_adj_su3_mat_vec_4dir:
	push		ebp
	mov		ebp,esp
	push		eax
	push		ebx
	push		ecx
	mov		eax,[ebp+8]			; su3_matrix *a[0]
	mov		ebx,[ebp+12]			; su3_vector *b
	mov		ecx,[ebp+16]			; su3_vector *c[0]

	;  bring in real and imaginary b vector
	movups		xmm0,[ebx]			; c1i,c1r,c0i,c0r  <(bb)->c[0]>
	movaps		xmm1,xmm0
	shufps		xmm1,xmm1,0xB1			; c1r,c1i,c0r,c0i
	movups		xmm2,[ebx+8]			; c2i,c2r,c1i,c1r  <(bb)->c[1]>
	shufps		xmm2,xmm2,0xEB			; c2i,c2r,c2r,c2i
	; xmm0:   c1i, c1r, c0i, c0r
	; xmm1:   c1r, c1i, c0r, c0i
	; xmm2:   c2i, c2r, c2r, c2i
	; *******************************************************************	
	; *******************************************************************	
	; bring in first column of a[0] matrix
	movlps		xmm3,[eax]			;                   <(aa)[0].e[0][0]>
	movhps		xmm3,[eax+24]			; c1i,c1r,c0i,c0r   <(aa)[0].e[1][0]>
	movhps		xmm4,[eax+48]			; c2i,c2r,x,x       <(aa)[0].e[2][0]>
	shufps		xmm4,xmm4,0xEE			; c2i,c2r,c2i,c2r
	; xmm3:   a1i, a1r, a0i, a0r
	; xmm4:   a2i, a2r, a2i, a2r
	movaps		xmm5,xmm3
	mulps		xmm3,xmm0
	mulps		xmm5,xmm1
	mulps		xmm4,xmm2
	; xmm3:  c1i.i, c1r.r, c0i.i, c0r.r
	; xmm4:  c2i.i, c2r.r, c2i.r, c2r.i
	; xmm5:  c1i.r, c1r.i, c0i.r, c0r.i
	movaps		xmm6,xmm5
	shufps		xmm6,xmm3,0x4E
	; xmm6:  c0i.i, c0r.r, c1i.r, c1r.i
	shufps		xmm5,xmm3,0xE4
	; xmm4:  c2i.i, c2r.r, c2i.r, c2r.i
	; xmm5:  c1i.i, c1r.r, c0i.r, c0r.i
	; xmm6:  c0i.i, c0r.r, c1i.r, c1r.i
	addps	xmm4,xmm5
	addps	xmm4,xmm6
	; *******************************************************************	
	; bring in second column of a[0] matrix
	movlps		xmm3,[eax+8]			;		    <(aa)[0].e[0][1]>
	movhps		xmm3,[eax+32]			; c1i,c1r,c0i,c0r   <(aa)[0].e[1][1]>
	movhps		xmm7,[eax+56]			; c2i,c2r,x,x       <(aa)[0].e[2][1]>
	shufps		xmm7,xmm7,0xEE			; c2i,c2r,c2i,c2r
	; xmm3:   a1i, a1r, a0i, a0r
	; xmm7:   a2i, a2r, a2i, a2r
	movaps		xmm5,xmm3
	mulps		xmm3,xmm0
	mulps		xmm5,xmm1
	mulps		xmm7,xmm2
	; xmm3:  c1i.i, c1r.r, c0i.i, c0r.r
	; xmm7:  c2i.i, c2r.r, c2i.r, c2r.i
	; xmm5:  c1i.r, c1r.i, c0i.r, c0r.i
	movaps		xmm6,xmm5
	shufps		xmm6,xmm3,0x4E
	; xmm6:  c0i.i, c0r.r, c1i.r, c1r.i
	shufps		xmm5,xmm3,0xE4
	; xmm7:  c2i.i, c2r.r, c2i.r, c2r.i
	; xmm5:  c1i.i, c1r.r, c0i.r, c0r.i
	; xmm6:  c0i.i, c0r.r, c1i.r, c1r.i
	addps	xmm7,xmm5
	addps	xmm7,xmm6
	; xmm4:  i-sum, r-sum, i.r-sum, r.i-sum  [0]
	; xmm7:  i-sum, r-sum, i.r-sum, r.i-sum  [1]
	movaps		xmm5,xmm4
	shufps		xmm5,xmm7,0x22
	shufps		xmm4,xmm7,0x77
	; xmm4:  i.r-sum[1], i-sum[1], i.r-sum[0], i-sum[0]
	; xmm5:  r.i-sum[1], r-sum[1], r.i-sum[0], r-sum[0]
	xorps		xmm4,[negate]			;                   <_sse_sgn24>
	addps		xmm5,xmm4
	movups		[ecx],xmm5			;		    <(cc)[0].c[0]>
	; *******************************************************************	
	; bring in third column of a[0] matrix
	movlps		xmm3,[eax+16]			;                   <(aa)[0].e[0][2]>
	movhps		xmm3,[eax+40]			; c1i,c1r,c0i,c0r   <(aa)[0].e[1][2]>
	movhps		xmm4,[eax+64]			; c2i,c2r,x,x	    <(aa)[0].e[2][2]>
	shufps		xmm4,xmm4,0xEE			; c2i,c2r,c2i,c2r
	; xmm3:   a1i, a1r, a0i, a0r
	; xmm4:   a2i, a2r, a2i, a2r
	movaps		xmm5,xmm3
	mulps		xmm3,xmm0
	mulps		xmm5,xmm1
	mulps		xmm4,xmm2
	; xmm3:  c1i.i, c1r.r, c0i.i, c0r.r
	; xmm4:  c2i.i, c2r.r, c2i.r, c2r.i
	; xmm5:  c1i.r, c1r.i, c0i.r, c0r.i
	movaps		xmm6,xmm5
	shufps		xmm6,xmm3,0x4E
	; xmm6:  c0i.i, c0r.r, c1i.r, c1r.i
	shufps		xmm5,xmm3,0xE4
	; xmm4:  c2i.i, c2r.r, c2i.r, c2r.i
	; xmm5:  c1i.i, c1r.r, c0i.r, c0r.i
	; xmm6:  c0i.i, c0r.r, c1i.r, c1r.i
	addps	xmm4,xmm5
	addps	xmm4,xmm6
	; xmm4:  i-sum, r-sum, i.r-sum, r.i-sum
	movaps		xmm5,xmm4
	shufps		xmm4,xmm4,0x77
	shufps		xmm5,xmm5,0x22
	; xmm4:  i.r-sum, i-sum, i.r-sum, i-sum
	; xmm5:  r.i-sum, r-sum, r.i-sum, r-sum
	xorps		xmm4,[negate]			;		   <_sse_sgn24>
	addps		xmm5,xmm4
	movhps		[ecx+16],xmm5			;                  <(cc)[0].c[2]>

	; *******************************************************************	
	; *******************************************************************	
	; bring in first column of a[1] matrix
	add		eax,72
	add		ecx,24
	movlps		xmm3,[eax]			;                  <(aa)[1].e[0][0]>
	movhps		xmm3,[eax+24]			; c1i,c1r,c0i,c0r  <(aa)[1].e[1][0]>
	movhps		xmm4,[eax+48]			; c2i,c2r,x,x	   <(aa)[1].e[2][0]>
	shufps		xmm4,xmm4,0xEE			; c2i,c2r,c2i,c2r
	; xmm3:   a1i, a1r, a0i, a0r
	; xmm4:   a2i, a2r, a2i, a2r
	movaps		xmm5,xmm3
	mulps		xmm3,xmm0
	mulps		xmm5,xmm1
	mulps		xmm4,xmm2
	; xmm3:  c1i.i, c1r.r, c0i.i, c0r.r
	; xmm4:  c2i.i, c2r.r, c2i.r, c2r.i
	; xmm5:  c1i.r, c1r.i, c0i.r, c0r.i
	movaps		xmm6,xmm5
	shufps		xmm6,xmm3,0x4E
	; xmm6:  c0i.i, c0r.r, c1i.r, c1r.i
	shufps		xmm5,xmm3,0xE4
	; xmm4:  c2i.i, c2r.r, c2i.r, c2r.i
	; xmm5:  c1i.i, c1r.r, c0i.r, c0r.i
	; xmm6:  c0i.i, c0r.r, c1i.r, c1r.i
	addps	xmm4,xmm5
	addps	xmm4,xmm6
	; *******************************************************************	
	; bring in second column of a[1] matrix
	movlps		xmm3,[eax+8]			;                    <(aa)[1].e[0][1]>
	movhps		xmm3,[eax+32]			; c1i,c1r,c0i,c0r    <(aa)[1].e[1][1]>
	movhps		xmm7,[eax+56]			; c2i,c2r,x,x	     <(aa)[1].e[2][1]>
	shufps		xmm7,xmm7,0xEE			; c2i,c2r,c2i,c2r
	; xmm3:   a1i, a1r, a0i, a0r
	; xmm7:   a2i, a2r, a2i, a2r
	movaps		xmm5,xmm3
	mulps		xmm3,xmm0
	mulps		xmm5,xmm1
	mulps		xmm7,xmm2
	; xmm3:  c1i.i, c1r.r, c0i.i, c0r.r
	; xmm7:  c2i.i, c2r.r, c2i.r, c2r.i
	; xmm5:  c1i.r, c1r.i, c0i.r, c0r.i
	movaps		xmm6,xmm5
	shufps		xmm6,xmm3,0x4E
	; xmm6:  c0i.i, c0r.r, c1i.r, c1r.i
	shufps		xmm5,xmm3,0xE4
	; xmm7:  c2i.i, c2r.r, c2i.r, c2r.i
	; xmm5:  c1i.i, c1r.r, c0i.r, c0r.i
	; xmm6:  c0i.i, c0r.r, c1i.r, c1r.i
	addps	xmm7,xmm5
	addps	xmm7,xmm6
	; xmm4:  i-sum, r-sum, i.r-sum, r.i-sum  [0]
	; xmm7:  i-sum, r-sum, i.r-sum, r.i-sum  [1]
	movaps		xmm5,xmm4
	shufps		xmm5,xmm7,0x22
	shufps		xmm4,xmm7,0x77
	; xmm4:  i.r-sum[1], i-sum[1], i.r-sum[0], i-sum[0]
	; xmm5:  r.i-sum[1], r-sum[1], r.i-sum[0], r-sum[0]
	xorps		xmm4,[negate]			;                     <_sse_sgn24>
	addps		xmm5,xmm4
	movups		[ecx],xmm5			;                     <(cc)[1].c[0]>
	; *******************************************************************	
	; bring in third column of a[1] matrix
	movlps		xmm3,[eax+16]			;                     <(aa)[1].e[0][2]>
	movhps		xmm3,[eax+40]			; c1i,c1r,c0i,c0r     <(aa)[1].e[1][2]>
	movhps		xmm4,[eax+64]			; c2i,c2r,x,x	      <(aa)[1].e[2][2]>
	shufps		xmm4,xmm4,0xEE			; c2i,c2r,c2i,c2r
	; xmm3:   a1i, a1r, a0i, a0r
	; xmm4:   a2i, a2r, a2i, a2r
	movaps		xmm5,xmm3
	mulps		xmm3,xmm0
	mulps		xmm5,xmm1
	mulps		xmm4,xmm2
	; xmm3:  c1i.i, c1r.r, c0i.i, c0r.r
	; xmm4:  c2i.i, c2r.r, c2i.r, c2r.i
	; xmm5:  c1i.r, c1r.i, c0i.r, c0r.i
	movaps		xmm6,xmm5
	shufps		xmm6,xmm3,0x4E
	; xmm6:  c0i.i, c0r.r, c1i.r, c1r.i
	shufps		xmm5,xmm3,0xE4
	; xmm4:  c2i.i, c2r.r, c2i.r, c2r.i
	; xmm5:  c1i.i, c1r.r, c0i.r, c0r.i
	; xmm6:  c0i.i, c0r.r, c1i.r, c1r.i
	addps	xmm4,xmm5
	addps	xmm4,xmm6
	; xmm4:  i-sum, r-sum, i.r-sum, r.i-sum
	movaps		xmm5,xmm4
	shufps		xmm4,xmm4,0x77
	shufps		xmm5,xmm5,0x22
	; xmm4:  i.r-sum, i-sum, i.r-sum, i-sum
	; xmm5:  r.i-sum, r-sum, r.i-sum, r-sum
	xorps		xmm4,[negate]			;                     <_sse_sgn24>
	addps		xmm5,xmm4
	movhps		[ecx+16],xmm5			;                     <(cc)[1].c[2]>

	; *******************************************************************	
	; *******************************************************************	
	; bring in first column of a[2] matrix
	add		eax,72
	add		ecx,24
	movlps		xmm3,[eax]			;                     <(aa)[2].e[0][0]>
	movhps		xmm3,[eax+24]			; c1i,c1r,c0i,c0r     <(aa)[2].e[1][0]>
	movhps		xmm4,[eax+48]			; c2i,c2r,x,x	      <(aa)[2].e[2][0]>
	shufps		xmm4,xmm4,0xEE			; c2i,c2r,c2i,c2r
	; xmm3:   a1i, a1r, a0i, a0r
	; xmm4:   a2i, a2r, a2i, a2r
	movaps		xmm5,xmm3
	mulps		xmm3,xmm0
	mulps		xmm5,xmm1
	mulps		xmm4,xmm2
	; xmm3:  c1i.i, c1r.r, c0i.i, c0r.r
	; xmm4:  c2i.i, c2r.r, c2i.r, c2r.i
	; xmm5:  c1i.r, c1r.i, c0i.r, c0r.i
	movaps		xmm6,xmm5
	shufps		xmm6,xmm3,0x4E
	; xmm6:  c0i.i, c0r.r, c1i.r, c1r.i
	shufps		xmm5,xmm3,0xE4
	; xmm4:  c2i.i, c2r.r, c2i.r, c2r.i
	; xmm5:  c1i.i, c1r.r, c0i.r, c0r.i
	; xmm6:  c0i.i, c0r.r, c1i.r, c1r.i
	addps	xmm4,xmm5
	addps	xmm4,xmm6
	; *******************************************************************	
	; bring in second column of a[1] matrix
	movlps		xmm3,[eax+8]			;                       <(aa)[2].e[0][1]>
	movhps		xmm3,[eax+32]			; c1i,c1r,c0i,c0r	<(aa)[2].e[1][1]>
	movhps		xmm7,[eax+56]			; c2i,c2r,x,x		<(aa)[2].e[2][1]>
	shufps		xmm7,xmm7,0xEE			; c2i,c2r,c2i,c2r
	; xmm3:   a1i, a1r, a0i, a0r
	; xmm7:   a2i, a2r, a2i, a2r
	movaps		xmm5,xmm3
	mulps		xmm3,xmm0
	mulps		xmm5,xmm1
	mulps		xmm7,xmm2
	; xmm3:  c1i.i, c1r.r, c0i.i, c0r.r
	; xmm7:  c2i.i, c2r.r, c2i.r, c2r.i
	; xmm5:  c1i.r, c1r.i, c0i.r, c0r.i
	movaps		xmm6,xmm5
	shufps		xmm6,xmm3,0x4E
	; xmm6:  c0i.i, c0r.r, c1i.r, c1r.i
	shufps		xmm5,xmm3,0xE4
	; xmm7:  c2i.i, c2r.r, c2i.r, c2r.i
	; xmm5:  c1i.i, c1r.r, c0i.r, c0r.i
	; xmm6:  c0i.i, c0r.r, c1i.r, c1r.i
	addps	xmm7,xmm5
	addps	xmm7,xmm6
	; xmm4:  i-sum, r-sum, i.r-sum, r.i-sum  [0]
	; xmm7:  i-sum, r-sum, i.r-sum, r.i-sum  [1]
	movaps		xmm5,xmm4
	shufps		xmm5,xmm7,0x22
	shufps		xmm4,xmm7,0x77
	; xmm4:  i.r-sum[1], i-sum[1], i.r-sum[0], i-sum[0]
	; xmm5:  r.i-sum[1], r-sum[1], r.i-sum[0], r-sum[0]
	xorps		xmm4,[negate]			;                       <_sse_sgn24>
	addps		xmm5,xmm4
	movups		[ecx],xmm5			;			<(cc)[2].c[0]>
	; *******************************************************************	
	; bring in third column of a[1] matrix
	movlps		xmm3,[eax+16]			;			<(aa)[2].e[0][2]>
	movhps		xmm3,[eax+40]			; c1i,c1r,c0i,c0r	<(aa)[2].e[1][2]>
	movhps		xmm4,[eax+64]			; c2i,c2r,x,x		<(aa)[2].e[2][2]>
	shufps		xmm4,xmm4,0xEE			; c2i,c2r,c2i,c2r
	; xmm3:   a1i, a1r, a0i, a0r
	; xmm4:   a2i, a2r, a2i, a2r
	movaps		xmm5,xmm3
	mulps		xmm3,xmm0
	mulps		xmm5,xmm1
	mulps		xmm4,xmm2
	; xmm3:  c1i.i, c1r.r, c0i.i, c0r.r
	; xmm4:  c2i.i, c2r.r, c2i.r, c2r.i
	; xmm5:  c1i.r, c1r.i, c0i.r, c0r.i
	movaps		xmm6,xmm5
	shufps		xmm6,xmm3,0x4E
	; xmm6:  c0i.i, c0r.r, c1i.r, c1r.i
	shufps		xmm5,xmm3,0xE4
	; xmm4:  c2i.i, c2r.r, c2i.r, c2r.i
	; xmm5:  c1i.i, c1r.r, c0i.r, c0r.i
	; xmm6:  c0i.i, c0r.r, c1i.r, c1r.i
	addps	xmm4,xmm5
	addps	xmm4,xmm6
	; xmm4:  i-sum, r-sum, i.r-sum, r.i-sum
	movaps		xmm5,xmm4
	shufps		xmm4,xmm4,0x77
	shufps		xmm5,xmm5,0x22
	; xmm4:  i.r-sum, i-sum, i.r-sum, i-sum
	; xmm5:  r.i-sum, r-sum, r.i-sum, r-sum
	xorps		xmm4,[negate]			;			<_sse_sgn24>
	addps		xmm5,xmm4
	movhps		[ecx+16],xmm5			;			<(cc)[2].c[2]>

	; *******************************************************************	
	; *******************************************************************	
	; bring in first column of a[3] matrix
	add		eax,72
	add		ecx,24
	movlps		xmm3,[eax]			;			<(aa)[3].e[0][0]>
	movhps		xmm3,[eax+24]			; c1i,c1r,c0i,c0r	<(aa)[3].e[1][0]>
	movhps		xmm4,[eax+48]			; c2i,c2r,x,x		<(aa)[3].e[2][0]>
	shufps		xmm4,xmm4,0xEE			; c2i,c2r,c2i,c2r
	; xmm3:   a1i, a1r, a0i, a0r
	; xmm4:   a2i, a2r, a2i, a2r
	movaps		xmm5,xmm3
	mulps		xmm3,xmm0
	mulps		xmm5,xmm1
	mulps		xmm4,xmm2
	; xmm3:  c1i.i, c1r.r, c0i.i, c0r.r
	; xmm4:  c2i.i, c2r.r, c2i.r, c2r.i
	; xmm5:  c1i.r, c1r.i, c0i.r, c0r.i
	movaps		xmm6,xmm5
	shufps		xmm6,xmm3,0x4E
	; xmm6:  c0i.i, c0r.r, c1i.r, c1r.i
	shufps		xmm5,xmm3,0xE4
	; xmm4:  c2i.i, c2r.r, c2i.r, c2r.i
	; xmm5:  c1i.i, c1r.r, c0i.r, c0r.i
	; xmm6:  c0i.i, c0r.r, c1i.r, c1r.i
	addps	xmm4,xmm5
	addps	xmm4,xmm6
	; *******************************************************************	
	; bring in second column of a[1] matrix
	movlps		xmm3,[eax+8]			;			<(aa)[3].e[0][1]>
	movhps		xmm3,[eax+32]			; c1i,c1r,c0i,c0r	<(aa)[3].e[1][1]>
	movhps		xmm7,[eax+56]			; c2i,c2r,x,x		<(aa)[3].e[2][1]>
	shufps		xmm7,xmm7,0xEE			; c2i,c2r,c2i,c2r
	; xmm3:   a1i, a1r, a0i, a0r
	; xmm7:   a2i, a2r, a2i, a2r
	movaps		xmm5,xmm3
	mulps		xmm3,xmm0
	mulps		xmm5,xmm1
	mulps		xmm7,xmm2
	; xmm3:  c1i.i, c1r.r, c0i.i, c0r.r
	; xmm7:  c2i.i, c2r.r, c2i.r, c2r.i
	; xmm5:  c1i.r, c1r.i, c0i.r, c0r.i
	movaps		xmm6,xmm5
	shufps		xmm6,xmm3,0x4E
	; xmm6:  c0i.i, c0r.r, c1i.r, c1r.i
	shufps		xmm5,xmm3,0xE4
	; xmm7:  c2i.i, c2r.r, c2i.r, c2r.i
	; xmm5:  c1i.i, c1r.r, c0i.r, c0r.i
	; xmm6:  c0i.i, c0r.r, c1i.r, c1r.i
	addps	xmm7,xmm5
	addps	xmm7,xmm6
	; xmm4:  i-sum, r-sum, i.r-sum, r.i-sum  [0]
	; xmm7:  i-sum, r-sum, i.r-sum, r.i-sum  [1]
	movaps		xmm5,xmm4
	shufps		xmm5,xmm7,0x22
	shufps		xmm4,xmm7,0x77
	; xmm4:  i.r-sum[1], i-sum[1], i.r-sum[0], i-sum[0]
	; xmm5:  r.i-sum[1], r-sum[1], r.i-sum[0], r-sum[0]
	xorps		xmm4,[negate]			;			<_sse_sgn24>
	addps		xmm5,xmm4
	movups		[ecx],xmm5			;			<(cc)[3].c[0]>
	; *******************************************************************	
	; bring in third column of a[1] matrix
	movlps		xmm3,[eax+16]			;			<(aa)[3].e[0][2]>
	movhps		xmm3,[eax+40]			; c1i,c1r,c0i,c0r	<(aa)[3].e[1][2]>
	movhps		xmm4,[eax+64]			; c2i,c2r,x,x		<(aa)[3].e[2][2]>
	shufps		xmm4,xmm4,0xEE			; c2i,c2r,c2i,c2r
	; xmm3:   a1i, a1r, a0i, a0r
	; xmm4:   a2i, a2r, a2i, a2r
	movaps		xmm5,xmm3
	mulps		xmm3,xmm0
	mulps		xmm5,xmm1
	mulps		xmm4,xmm2
	; xmm3:  c1i.i, c1r.r, c0i.i, c0r.r
	; xmm4:  c2i.i, c2r.r, c2i.r, c2r.i
	; xmm5:  c1i.r, c1r.i, c0i.r, c0r.i
	movaps		xmm6,xmm5
	shufps		xmm6,xmm3,0x4E
	; xmm6:  c0i.i, c0r.r, c1i.r, c1r.i
	shufps		xmm5,xmm3,0xE4
	; xmm4:  c2i.i, c2r.r, c2i.r, c2r.i
	; xmm5:  c1i.i, c1r.r, c0i.r, c0r.i
	; xmm6:  c0i.i, c0r.r, c1i.r, c1r.i
	addps	xmm4,xmm5
	addps	xmm4,xmm6
	; xmm4:  i-sum, r-sum, i.r-sum, r.i-sum
	movaps		xmm5,xmm4
	shufps		xmm4,xmm4,0x77
	shufps		xmm5,xmm5,0x22
	; xmm4:  i.r-sum, i-sum, i.r-sum, i-sum
	; xmm5:  r.i-sum, r-sum, r.i-sum, r-sum
	xorps		xmm4,[negate]			;			<_sse_sgn24>
	addps		xmm5,xmm4
	movhps		[ecx+16],xmm5			;			<(cc)[3].c[2]>

	; *******************************************************************	

here:	pop	ecx
	pop	ebx
	pop	eax
	mov	esp,ebp
	pop	ebp
	ret
	
	align		16
negate:	dd		0x00000000
	dd		0x80000000
	dd		0x00000000
	dd		0x80000000
	
