/***************************************************
    Include file for the code, that reads 
    propagating heavy-->light form factors
 **************************************************/


#ifndef  READ_HL_FORM_INCLUDE
#define  READ_HL_FORM_INCLUDE 1


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <string.h>

typedef struct {           /* standard complex number declaration for single- */
   Real real;             /* precision complex numbers                       */
   Real imag;
 } complex;


/*** memory location for the form factor operators ******/

#define FORM_WHERE(t,zonk_pt,seq_pt,spect_pt,q_pt,p_pt,oper_pt) \
         (t) + nt*((zonk_pt) + no_zonked*((seq_pt) + no_sequential*((spect_pt) + no_spectator*((q_pt) + no_q_values*((p_pt) + no_p_values*(oper_pt))))))


enum  form_type { HEAVY_TO_HEAVY = 10 , HEAVY_TO_LIGHT  } ;


/**** these flags are used for the binary IO of the two point functions  ****/
enum  two_form_type { HL_2PT_BAG = 20 , HL_REL_2PT  , LL_2PT , HL_LOCAL_SINK } ; 

#define TWOPT_FORM_WHERE(t,zonk_pt,spect_pt,q_pt,oper_pt) \
         (t) + nt*((zonk_pt) + no_zonked*((spect_pt) + no_spectator*((q_pt) + no_q_values*(oper_pt) )))


enum byte_rev_option { do_byte_rev = 10  , do_nothing  } ; 

#define MAX_NAME 300
#define MAX_NO_FILE 120
#define MAX_SELECT_PER_FILE 20

#endif




