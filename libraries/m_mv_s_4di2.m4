// mult_su3_mat_vec_sum_4dir( su3_matrix *a, su3_vector *b0,
// su3_vector *b1, su3_vector *b2, su3_vector *b3, su3_vector *c)
// c <- A[0]*B0+A[1]*B1+A[2]*B2+A[3]*B3
// file m_mv_s_4dir.m4, i860 assembler version of m_mv_s_4dir.c
//
// Register usage:
// f8,f9   = b[0]
// f10,f11 = b[1]
// f12,f13 = b[2]
// f14,f15 = a[1][0]
// f16,f17 = a[1][1]
// f18,f19 = a[1][2]
// f20,f21 = c[0]
// f22,f23 = c[1]
// f24,f25 = c[2]
// f26,f27 = a[0][0] and a[2][0]
// f28,f29 = a[0][1] and a[2][1]
// f30,f31 = a[0][2] and a[2][2]
    define(A,r16)	// address of vector of su3_matrices
    define(B0,r17)	// address of vector to be multiplied
    define(B1,r18)	// address of vector to be multiplied
    define(B2,r19)	// address of vector to be multiplied
    define(B3,r20)	// address of vector to be multiplied
    define(C,r21)	// address of result
    define(B,r22)	// temporary address of vector
    define(INC,r23)	// increment for loop
    define(COUNT,r24)	// loop counter

    define(b0,f8)	// complex number = register pair
    define(b0r,f8)	// real part
    define(b0i,f9)	// imag part
    define(b1,f10)
    define(b1r,f10)
    define(b1i,f11)
    define(b2,f12)
    define(b2r,f12)
    define(b2i,f13)

    define(c0,f20)
    define(c0r,f20)
    define(c0i,f21)
    define(c1,f22)
    define(c1r,f22)
    define(c1i,f23)
    define(c2,f24)
    define(c2r,f24)
    define(c2i,f25)

    define(a00,f26)
    define(a00r,f26)
    define(a00i,f27)
    define(a01,f28)
    define(a01r,f28)
    define(a01i,f29)
    define(a02,f30)
    define(a02r,f30)
    define(a02i,f31)

    define(a10,f14)
    define(a10r,f14)
    define(a10i,f15)
    define(a11,f16)
    define(a11r,f16)
    define(a11i,f17)
    define(a12,f18)
    define(a12r,f18)
    define(a12i,f19)

    define(a20,f26)
    define(a20r,f26)
    define(a20i,f27)
    define(a21,f28)
    define(a21r,f28)
    define(a21i,f29)
    define(a22,f30)
    define(a22r,f30)
    define(a22i,f31)

// First accumulate c0 real and imaginary parts and c1 real part,
// then c2 real and imaginary and c1.imag


	.text
//First vector
	.align	8
_mult_su3_mat_vec_sum_4dir:
.align 8
	// start dual mode, start fetching
	// enter zeroes into pipeline
					pfld.d   0(A),f0
					pfld.d   0(B0),f0
					pfld.d   24(A),f0
.align 8
	d.pfadd.ss f0,f0,f0;		pfld.d	8(B0),a00
	d.pfadd.ss f0,f0,f0;		pfld.d	8(A),b0
	d.pfadd.ss f0,f0,f0;		pfld.d	32(A),a10
	// start zero'th column of A times B0 down pipeline
	d.pfmul.ss a00r,b0r,f0;		nop
	d.pfmul.ss a00r,b0i,f0;		pfld.d	16(A),b1
	d.pfmul.ss a10r,b0r,f0;		nop
	// clear part of C
	d.m12apm.ss a00i,b0i,c2r;	pfld.d 	16(B0),a01
	d.m12apm.ss a00i,b0r,c2i;	nop
	d.m12apm.ss a10i,b0i,c1i;	pfld.d	40(A),a11
	// first column of A
	d.m12asm.ss a01r,b1r,f0;	adds -1,r0,INC
	d.m12apm.ss a01r,b1i,f0;	pfld.d	48(A),a02
	d.m12asm.ss a11r,b1r,f0;	or 2,r0,COUNT
	d.m12apm.ss a01i,b1i,f0;	pfld.d	56(A),b2
	d.m12apm.ss a01i,b1r,f0;	bla INC,COUNT,LOOP
	d.m12apm.ss a11i,b1i,f0;	pfld.d	64(A),a12
	// second column of A
LOOP:
	d.m12asm.ss a02r,b2r,f0;	adds 72,A,A
	d.m12apm.ss a02r,b2i,f0;	pfld.d	0(A),a20
	d.m12asm.ss a12r,b2r,f0;	nop
	d.m12apm.ss a02i,b2i,f0;	pfld.d	0(B1),a21
	m12apm.ss a02i,b2r,f0;		nop
	m12apm.ss a12i,b2i,f0;		pfld.d	24(A),a22
	// start multiplies for second half, where we accumulate c[2]
	// real and imaginary and c[1] imaginary.  Sums from first half
	// still going through adder.
	m12asm.ss a20r,b0r,f0
	m12apm.ss a20r,b0i,f0
	m12asm.ss a10r,b0i,f0
	// Now enter zeroes into adder pipe, while results from first
	// half are coming out.
	pfadd.ss	c2r,f0,c0r
	pfadd.ss	c2i,f0,c0i
	pfadd.ss	c1i,f0,c1r
	// continue with multiplies in second half, column 0 of A
	m12apm.ss a20i,b0i,f0
	m12apm.ss a20i,b0r,f0
	m12apm.ss a10i,b0r,f0
	// Row 1 of A
	m12asm.ss a21r,b1r,f0
	m12apm.ss a21r,b1i,f0
	m12apm.ss a11r,b1i,f0
	m12apm.ss a21i,b1i,f0
	m12apm.ss a21i,b1r,f0
.align	8
	d.m12apm.ss a11i,b1r,f0
	// Row 2 of A, start fetching for next vector
	d.m12asm.ss a22r,b2r,f0
	d.m12apm.ss a22r,b2i,f0;	pfld.d	8(B1),a00	
	d.m12apm.ss a12r,b2i,f0;	mov B1,B	
	d.m12apm.ss a22i,b2i,f0;	pfld.d	8(A),b0	
	d.m12apm.ss a22i,b2r,f0;	mov B2,B1	
	d.m12apm.ss a12i,b2r,f0;	pfld.d	32(A),a10	

// Next vector
// First accumulate c0 real and imaginary parts and c1 real part,
// then c2 real and imaginary and c1.imag
					
	// start zero'th column of A times B down pipeline
	// previous vector still going into adder
	d.m12asm.ss a00r,b0r,f0;	mov B3,B2
	d.m12apm.ss a00r,b0i,f0;	pfld.d	16(A),b1
	d.m12apm.ss a10r,b0r,f0;	nop
	// enter C into pipeline, results from last vector coming out
	d.pfadd.ss f0,c0r,c2r;		nop
	d.pfadd.ss f0,c0i,c2i;		nop
	d.pfadd.ss f0,c1r,c1i;		nop
	// continue with multiplications for new vector
	d.m12apm.ss a00i,b0i,f0;	pfld.d 	16(B),a01
	d.m12apm.ss a00i,b0r,f0;	nop
	d.m12apm.ss a10i,b0i,f0;	pfld.d	40(A),a11
	// first column of A
	d.m12asm.ss a01r,b1r,f0;	nop
	d.m12apm.ss a01r,b1i,f0;	pfld.d	48(A),a02
	d.m12asm.ss a11r,b1r,f0;	nop
	d.m12apm.ss a01i,b1i,f0;	pfld.d	56(A),b2
	d.m12apm.ss a01i,b1r,f0;	bla INC,COUNT,LOOP
	d.m12apm.ss a11i,b1i,f0;	pfld.d	64(A),a12
	// go to LOOP from here
	// second column of A
	d.m12asm.ss a02r,b2r,f0;	nop
	d.m12apm.ss a02r,b2i,f0;	pfld.d	64(A),a20
	d.m12asm.ss a12r,b2r,f0;	nop
	d.m12apm.ss a02i,b2i,f0;	pfld.d	0(B1),a21
	m12apm.ss a02i,b2r,f0;		nop
	m12apm.ss a12i,b2i,f0;		pfld.d	64(A),a22
	// start multiplies for second half, where we accumulate c[2]
	// real and imaginary and c[1] imaginary.  Sums from first half
	// still going through adder.
	// end dual mode
	m12asm.ss a20r,b0r,f0
	m12apm.ss a20r,b0i,f0
	m12asm.ss a10r,b0i,f0
	// Now enter C into adder pipe, while results from first
	// half are coming out.
	pfadd.ss	f0,c2r,c0r
	pfadd.ss	f0,c2i,c0i
	pfadd.ss	f0,c1i,c1r
	// continue with multiplies in second half, column 0 of A
	m12apm.ss a20i,b0i,f0
	m12apm.ss a20i,b0r,f0
	m12apm.ss a10i,b0r,f0
	// Row 1 of A
	m12asm.ss a21r,b1r,f0
	m12apm.ss a21r,b1i,f0
	m12apm.ss a11r,b1i,f0
	m12apm.ss a21i,b1i,f0
	m12apm.ss a21i,b1r,f0
	m12apm.ss a11i,b1r,f0
	// Row 2 of A
	m12asm.ss a22r,b2r,f0
	m12apm.ss a22r,b2i,f0
	m12apm.ss a12r,b2i,f0
	m12apm.ss a22i,b2i,f0
	m12apm.ss a22i,b2r,f0
	m12apm.ss a12i,b2r,f0
	// Empty multiplier pipe
	m12asm.ss f0,f0,f0
	m12apm.ss f0,f0,f0
.align 8
	d.m12apm.ss f0,f0,f0
	// empty adder pipe, store results
	d.pfadd.ss f0,f0,c2r
	pfadd.ss f0,f0,c2i;		fst.d	c0,0(C)
	pfadd.ss f0,f0,c1i;		fst.d	c2,16(C)
    bri r1
					fst.d	c1,8(C)

.globl	_mult_su3_mat_vec_sum_4dir
