// mult_su3_nn( su3_matrix *a, su3_matrix *b, su3_matrix *c)
// c <- a*b
// file m_mat_nn.m4, i860 assembler version of m_mat_nn.c
//
// Register usage:
// r19 = increment = -1
// r20 = loop counter
// f8,f9   = b[0][column]
// f10,f11 = b[1][column]
// f12,f13 = b[2][column]
// f14,f15 = a[1][0]
// f16,f17 = a[1][1]
// f18,f19 = a[1][2]
// f20,f21 = c[0][column] and c[2][column]
// f22,f23 = c[1][column]
// f24,f25 = a[0][0] and a[2][0]
// f26,f27 = a[0][1] and a[2][1]
// f28,f29 = a[0][2] and a[2][2]
    define(A,r16)	// address of matrix
    define(B,r17)	// address of vector to be multiplied
    define(C,r18)	// address of result

	// for each column of B
    define(b0,f8)	// complex number = register pair
    define(b0r,f8)	// real part
    define(b0i,f9)	// imag part
    define(b1,f10)
    define(b1r,f10)
    define(b1i,f11)
    define(b2,f12)
    define(b2r,f12)
    define(b2i,f13)

	// for each column of C
    define(c0,f20)
    define(c0r,f20)
    define(c0i,f21)
    define(c1,f22)
    define(c1r,f22)
    define(c1i,f23)
    define(c2,f20)
    define(c2r,f20)
    define(c2i,f21)

    define(a00,f24)
    define(a00r,f24)
    define(a00i,f25)
    define(a01,f26)
    define(a01r,f26)
    define(a01i,f27)
    define(a02,f28)
    define(a02r,f28)
    define(a02i,f29)

    define(a10,f14)
    define(a10r,f14)
    define(a10i,f15)
    define(a11,f16)
    define(a11r,f16)
    define(a11i,f17)
    define(a12,f18)
    define(a12r,f18)
    define(a12i,f19)

    define(a20,f24)
    define(a20r,f24)
    define(a20i,f25)
    define(a21,f26)
    define(a21r,f26)
    define(a21i,f27)
    define(a22,f28)
    define(a22r,f28)
    define(a22i,f29)

// First accumulate c0 real and imaginary parts and c1 real part,
// then c2 real and imaginary and c1.imag


	.text
	.align	8
_mult_su3_nn:
	// loop over columns of B and C
					adds -1,r0,r19
					or 1,r0,r20
.align 8
     	d.pfadd.ss f0,f0,f0;		bla r19,r20,DUMMY
	// start dual mode, start fetching
	// enter zeroes into pipeline
	d.pfadd.ss f0,f0,f0;		fld.d	0(A),a00
DUMMY:
	d.pfadd.ss f0,f0,f0;		fld.d	0(B),b0
	d.pfadd.ss f0,f0,f0;		fld.d	24(A),a10
	// start zero'th column of A times B down pipeline
	d.pfmul.ss a00r,b0r,f0;		nop
	d.pfmul.ss a00r,b0i,f0;		fld.d	24(B),b1
	d.pfmul.ss a10r,b0r,f0;		nop
LOOP:
	d.m12apm.ss a00i,b0i,f0;	fld.d 	8(A),a01
	d.m12apm.ss a00i,b0r,f0;	nop
	d.m12apm.ss a10i,b0i,f0;	fld.d	32(A),a11
	// first column of A
	d.m12asm.ss a01r,b1r,f0;	nop
	d.m12apm.ss a01r,b1i,f0;	fld.d	16(A),a02
	d.m12asm.ss a11r,b1r,f0;	nop
	d.m12apm.ss a01i,b1i,f0;	fld.d	48(B),b2
	d.m12apm.ss a01i,b1r,f0;	adds 8,B,B	// next column
	d.m12apm.ss a11i,b1i,f0;	fld.d	40(A),a12
	// second column of A
	d.m12asm.ss a02r,b2r,f0;	nop
	d.m12apm.ss a02r,b2i,f0;	fld.d	48(A),a20
	d.m12asm.ss a12r,b2r,f0;	nop
	d.m12apm.ss a02i,b2i,f0;	fld.d	56(A),a21
	m12apm.ss a02i,b2r,f0;	nop
	m12apm.ss a12i,b2i,f0;	fld.d	64(A),a22
	// start multiplies for second half, where we accumulate c[2]
	// real and imaginary and c[1] imaginary.  Sums from first half
	// still going through adder.
	// write to c1 is correct for second and third passes, result
	// from first pass will be overwritten
.align 8
	d.m12asm.ss a20r,b0r,f0
	d.m12apm.ss a20r,b0i,f0
	d.m12asm.ss a10r,b0i,f0;	nop
	// Now enter zeroes into adder pipe, while results from first
	// half are coming out.
	d.pfadd.ss	f0,f0,c0r;	fst.d	c1,16(C)
	pfadd.ss	f0,f0,c0i;	nop
	pfadd.ss	f0,f0,c1r;	fst.d	c0,0(C)
	// continue with multiplies in second half, column 0 of A
	// end dual mode
	m12apm.ss a20i,b0i,f0
	m12apm.ss a20i,b0r,f0
	m12apm.ss a10i,b0r,f0
	// Row 1 of A
	m12asm.ss a21r,b1r,f0
	m12apm.ss a21r,b1i,f0
	m12apm.ss a11r,b1i,f0
	m12apm.ss a21i,b1i,f0
	m12apm.ss a21i,b1r,f0
.align 8
	d.m12apm.ss a11i,b1r,f0
	// Row 2 of A
	d.m12asm.ss a22r,b2r,f0
	d.m12apm.ss a22r,b2i,f0;	fld.d 0(A),a00
	d.m12apm.ss a12r,b2i,f0;	nop
	d.m12apm.ss a22i,b2i,f0;	fld.d 0(B),b0
	d.m12apm.ss a22i,b2r,f0;	nop
	d.m12apm.ss a12i,b2r,f0;	fld.d 24(A),a10
	// Empty multiplier pipe
	d.m12asm.ss a00r,b0r,f0;	nop
	d.m12apm.ss a00r,b0i,f0;	fld.d 24(B),b1
	d.m12apm.ss a10r,b0r,f0;	adds 8,C,C	// next column
	// empty adder pipe, store results, 
	d.pfadd.ss f0,f0,c2r;		nop
	d.pfadd.ss f0,f0,c2i;		bla r19,r20,LOOP
	d.pfadd.ss f0,f0,c1i;		fst.d	c2,40(C)
	// the following are executed only for the third column of C
	d.m12apm.ss a00i,b0i,f0;	fld.d 	8(A),a01
	d.m12apm.ss a00i,b0r,f0;	nop
	d.m12apm.ss a10i,b0i,f0;	fld.d	32(A),a11
	// first column of A
	d.m12asm.ss a01r,b1r,f0;	nop
	d.m12apm.ss a01r,b1i,f0;	fld.d	16(A),a02
	d.m12asm.ss a11r,b1r,f0;	nop
	d.m12apm.ss a01i,b1i,f0;	fld.d	48(B),b2
	d.m12apm.ss a01i,b1r,f0;	nop
	d.m12apm.ss a11i,b1i,f0;	fld.d	40(A),a12
	// second column of A
	d.m12asm.ss a02r,b2r,f0;	nop
	d.m12apm.ss a02r,b2i,f0;	fld.d	48(A),a20
	d.m12asm.ss a12r,b2r,f0;	nop
	d.m12apm.ss a02i,b2i,f0;	fld.d	56(A),a21
	m12apm.ss a02i,b2r,f0;	nop
	m12apm.ss a12i,b2i,f0;	fld.d	64(A),a22
	// start multiplies for second half, where we accumulate c[2]
	// real and imaginary and c[1] imaginary.  Sums from first half
	// still going through adder.
	// write to c1 is correct for second and third passes, result
	// from first pass will be overwritten
.align 8
	d.m12asm.ss a20r,b0r,f0
	d.m12apm.ss a20r,b0i,f0
	d.m12asm.ss a10r,b0i,f0;	nop
	// Now enter zeroes into adder pipe, while results from first
	// half are coming out.
	d.pfadd.ss	f0,f0,c0r;	fst.d	c1,16(C)
	pfadd.ss	f0,f0,c0i;	nop
	pfadd.ss	f0,f0,c1r;	fst.d	c0,0(C)
	// continue with multiplies in second half, column 0 of A
	// end dual mode
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
	m12asm.ss a00r,b0r,f0
	m12apm.ss a00r,b0i,f0
	m12apm.ss a10r,b0r,f0
	// empty adder pipe, store results, 
.align 8
	d.pfadd.ss f0,f0,c2r
	pfadd.ss f0,f0,c2i
	pfadd.ss f0,f0,c1i;		fst.d	c2,48(C)
					bri r1
					fst.d	c1,24(C)
.globl	_mult_su3_nn
