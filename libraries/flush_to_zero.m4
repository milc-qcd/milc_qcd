// flush_to_zero.m4
// set floating status register so underflow does not cause trap

.text
.align 8
_flush_to_zero:
        ld.c    fsr,r31         // load floating-point status register into r31
        or      1,r31,r31       // set low-order bit to 1
        andnot      32,r31,r31       // set fp trap enable to zero
        bri     r1              // return
        st.c    r31,fsr         // shadow instruction, store r31 to fsr

.globl _flush_to_zero

