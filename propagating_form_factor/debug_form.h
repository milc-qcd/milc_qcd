#ifdef DEBUGDEF
#ifndef DEBUG_FORM_INCLUDE
#define DEBUG_FORM_INCLUDE

/*
 *  HEADER FILE FOR DEBUG ROUTINES
 */


/*
 *  This file should NOT be a dependancy in the makefile
 *
 *
 */

/* 
   Flags to insert debug code into the production program.
   Because of the name space considerations, the end of each
   flag should be _DD
*/


/** prototypes for the debug functions ******/

void  dump_spin_wilson_vector(spin_wilson_vector *in )  ;
void dump_wilson_vector(field_offset prop) ;
void wilson_vector_pion(field_offset prop) ; 
void dump_heavy_smear_func(void) ;


void dump_seq_smear_func(void) ;
/******* end of the include file for the debug functions ******/


#endif
#endif
