#ifndef _MGRID_H
	#define _MGRID_H

	/* This file is for external access */

	// allow the user to specify GRID_Precision and/or GRID_PrecisionInt
	// default to single precision.  This macro is used to convert generic
	// names for types and procedures to precision-specific names
	// This macro can be defined for the entire compilation or for
	// individual compilation units.  It must be defined before including 
	// this file.
	// This logic is the same as the QOP logic.  It is a bit complicated
	// because it allows two ways to define the macro 
	// and it sets a default of single precision.

	#ifndef GRID_Precision
		#ifndef GRID_PrecisionInt
			//#warning "GRID_PrecisionInt NOT Defined"
			#define GRID_PrecisionInt 2
		#else
			//#warning "GRID_PrecisionInt Defined"
		#endif

		#if GRID_PrecisionInt == 1
			#define GRID_Precision 'F'
			#define GRID_PrecisionLetter F
		#elif GRID_PrecisionInt == 2
			#define GRID_Precision 'D'
			#define GRID_PrecisionLetter D
		#else
			#error "bad GRID_PrecisionInt"
		#endif
	#else
		#ifndef GRID_PrecisionInt
			#if GRID_Precision == 'F'
				#define GRID_PrecisionInt 1
				#define GRID_PrecisionLetter F
			#elif GRID_Precision == 'D'
				#define GRID_PrecisionInt 2
				#define GRID_PrecisionLetter D
			#else
				#error "bad GRID_Precision"
			#endif
		#else
			#if GRID_Precision == 'F'
				#if GRID_PrecisionInt != 1
					#error "inconsistent GRID_Precision='F' and GRID_PrecisionInt"
				#endif
				#define GRID_PrecisionLetter F
			#elif GRID_Precision == 'D'
				#if GRID_PrecisionInt != 2
					#error "inconsistent GRID_Precision='D' and GRID_PrecisionInt"
				#endif
				#define GRID_PrecisionLetter D
			#else
				#error "bad GRID_Precision"
			#endif
		#endif
	#endif

	#include "../include/mGrid/mGrid_int.h"
        #include "../include/mGrid/mGrid_f3.h"
        #include "../include/mGrid/mGrid_d3.h"

#endif /* _MGRID_H */
