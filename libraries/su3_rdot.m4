// float su3_rdot( su3_vector *a,*b );
// return real part of dot product of A and B
// file su3_rdot.m4, i860 assembler version of su3_rdot.c
//
    define(A,r16)	// address of source 1
    define(B,r17)	// address of source 2
    define(sum,f8)	// return value

    define(a0,f8)	// complex number = register pair
    define(a0r,f8)	// real part
    define(a0i,f9)	// imag part
    define(a1,f10)
    define(a1r,f10)
    define(a1i,f11)
    define(a2,f12)
    define(a2r,f12)
    define(a2i,f13)

    define(b0,f14)	// complex number = register pair
    define(b0r,f14)	// real part
    define(b0i,f15)	// imag part
    define(b1,f16)
    define(b1r,f16)
    define(b1i,f17)
    define(b2,f18)
    define(b2r,f18)
    define(b2i,f19)

	.text
	.align	8
_su3_rdot:
				fld.d	0(A),a0
				fld.d	8(A),a1
				fld.d	16(A),a2
				fld.d	0(B),b0
				fld.d	8(B),b1
				fld.d	16(B),b2
	pfmul.ss a0r,b0r,f0
	pfmul.ss a0i,b0i,f0
	pfmul.ss a1r,b1r,f0
	mm12ttpm.ss a1i,b1i,f0	// a0r*b0r into T reg.
	m12tpm.ss a2r,b2r,f0	// start add of a0r*b0r and a0i*b0i
	mm12ttpm.ss a2i,b2i,f0	// a1r*b1r into T reg.
	m12tpm.ss f0,f0,f0	// start add of a1r*b1r and a1i*b1i
				// adder output is now OK
	m12ttpa.ss f0,f0,a0r	// a0r*b0r + a0i*b0i comes out
	m12tpm.ss f0,f0,f0	// start add of a2r*b2r and a2i*b2i

	i2ap1.ss a0r,f0,f0	// start (a0r*b0r+a0i*b0i) + (a1r*b1r+a1i*b1i)
	pfadd.ss f0,f0,f0
	pfadd.ss f0,f0,a0r	// a0r gets (a2r*b2r+a2i*b2i)
//.if FLOATOPTION=X167
	i2ap1.ss a0r,f0,f0
	pfadd.ss f0,f0,f0
	pfadd.ss f0,f0,f0
		bri r1
	pfadd.ss f0,f0,sum
//.else
	//i2ap1.sd a0r,f0,f0
	//pfadd.sd f0,f0,f0
	//pfadd.sd f0,f0,f0
		//bri r1
	//pfadd.sd f0,f0,sum
//.endif
.globl	_su3_rdot
