// float magsq_su3vec( su3_vector *a );
// return |a|^2 = squared magnitude of vector
// file msq_su3vec.m4, i860 assembler version of msq_su3vec.c
//
    define(A,r16)	// address of source 1
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

	.text
	.align	8
_magsq_su3vec:
				fld.d	0(A),a0
				fld.d	8(A),a1
				fld.d	16(A),a2
	pfmul.ss a0r,a0r,f0
	pfmul.ss a0i,a0i,f0
	pfmul.ss a1r,a1r,f0
	mm12ttpm.ss a1i,a1i,f0	// a0r^2 into T reg.
	m12tpm.ss a2r,a2r,f0	// start add of a0r^2 and a0i^2
	mm12ttpm.ss a2i,a2i,f0	// a1r^2 into T reg.
	m12tpm.ss f0,f0,f0	// start add of a1r^2 and a1i^2
				// adder output is now OK
	m12ttpa.ss f0,f0,a0r	// a0r^2 + a0i^2 comes out
	m12tpm.ss f0,f0,f0	// start add of a2r^2 and a2i^2

	i2ap1.ss a0r,f0,f0	// start (a0r^2+a0i^2) + (a1r^2+a1i^2)
	pfadd.ss f0,f0,f0
	pfadd.ss f0,f0,a0r	// a0r gets (a2r^2+a2i^2)
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
.globl _magsq_su3vec
