#ifndef DEBUG_FORM_INCLUDE
#define DEBUG_FORM_INCLUDE

/*
 *  HEADER FILE FOR DEBUG ROUTINES
 */


/*
 *  This file should NOT be a dependency in the makefile
 *
 *
 */

/* 
   Flags to insert debug code into the production program.
   Because of the name space considerations, the end of each
   flag should be _DD
*/


/*** flag that controls whether the correlators are written out ***/
/*** #define RAW_DUMP_THE_CORR_DD 1  ****/


/** prototypes for the debug functions ******/

void hqet_pion() ;
void wilson_vector_pion(field_offset prop) ;
void zero_wilson_vector(field_offset prop) ;
void loadup_debug(int spin) ;

void  dump_spin_wilson_vector(spin_wilson_vector *in )  ;

void  fake_spin_wilson_src(spin_wilson_vector *in, int color )  ;
void fake_seq_src(int color) ; 


void dump_wilson_vector(field_offset prop) ; 

void mult_gamma5_sink(field_offset ans) ;


void unit_quark_zonked(int color);
void unit_quark_sequential( int color) ;
void gamma_five_sequential(int color) ; 

void dump_quark_sequential(void) ;
void dump_quark_zonked(void) ;

void zero_hqet(void) ;
void unit_hqet(void ) ;

/******* end of the include file for the debug functions ******/

#endif
