#define        v0       $0    /*Integer return value register*/
#define        t0       $1    /*Integer scratch registers (caller saved)*/
#define        t1       $2
#define        t2       $3
#define        t3       $4
#define        t4       $5
#define        t5       $6
#define        t6       $7
#define        t7       $8
#define        s0       $9    /*Integer save registers (callee saved)*/
#define        s1       $10
#define        s2       $11
#define        s3       $12
#define        s4       $13
#define        s5       $14
#define        fp       $15    /*Private frame pointer register*/
#define        a0       $16    /*Integer argument registers*/
#define        a1       $17
#define        a2       $18
#define        a3       $19
#define        a4       $20
#define        a5       $21
#define        t8       $22    /*Scratch registers (continued)*/
#define        t9       $23
#define        t10      $24
#define        t11      $25
#define        ra       $26    /*Return address register*/
#define        t12      $27    /*Scratch registers (continued)*/
#define        at      $28     #reserved for assembler
#define        gp      $29	/*global pointer*/
#define        sp       $30     /*Stack pointer register*/
#define        zero     $31     /*Integer ReadAsZero/Sink register*/

#define        fv0      $f0     /*Floating-point return value register*/
#define        fv1      $f1
#define        fs0      $f2     /*Floating-point save registers (callee saved)*/
#define        fs1      $f3
#define        fs2      $f4
#define        fs3      $f5
#define        fs4      $f6
#define        fs5      $f7
#define        fs6      $f8
#define        fs7      $f9
#define        ft0      $f10    /*Floating-point scratch registers*/
#define        ft1      $f11
#define        ft2      $f12
#define        ft3      $f13
#define        ft4      $f14
#define        ft5      $f15
#define        fa0      $f16    /*Floating-point argument registers*/
#define        fa1      $f17
#define        fa2      $f18
#define        fa3      $f19
#define        fa4      $f20
#define        fa5      $f21
#define        ft6      $f22    /*Floating-point scratch registers (continued)*/
#define        ft7      $f23
#define        ft8      $f24
#define        ft9      $f25
#define        ft10     $f26
#define        ft11     $f27
#define        ft12     $f28
#define        ft13     $f29
#define        ft14     $f30
#define        fzero    $f31    /*Floating-point ReadAsZero/Sink register*/
