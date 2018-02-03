#ifndef QPHIXJ_H
#define QPHIXJ_H

/* This file is for external access */

  // allow the user to specify QPHIXJ_Precision and/or QPHIXJ_PrecisionInt
  // default to single precision.  This macro is used to convert generic
  // names for types and procedures to precision-specific names
  // This macro can be defined for the entire compilation or for
  // individual compilation units.  It must be defined before including 
  // this file.
  // This logic is the same as the QOP logic.  It is a bit complicated
  // because it allows two ways to define the macro 
  // and it sets a default of single precision.

#ifndef QPHIXJ_Precision
#  ifndef QPHIXJ_PrecisionInt
//# warning "QPHIXJ_PrecisionInt NOT Defined"
#    define QPHIXJ_PrecisionInt 1
#else
//# warning "QPHIXJ_PrecisionInt Defined"
#  endif
#  if QPHIXJ_PrecisionInt == 1
#    define QPHIXJ_Precision 'F'
#    define QPHIXJ_PrecisionLetter F
#  elif QPHIXJ_PrecisionInt == 2
#    define QPHIXJ_Precision 'D'
#    define QPHIXJ_PrecisionLetter D
#  else
#    error "bad QPHIXJ_PrecisionInt"
#  endif
#else
#  ifndef QPHIXJ_PrecisionInt
#    if QPHIXJ_Precision == 'F'
#      define QPHIXJ_PrecisionInt 1
#      define QPHIXJ_PrecisionLetter F
#    elif QPHIXJ_Precision == 'D'
#      define QPHIXJ_PrecisionInt 2
#      define QPHIXJ_PrecisionLetter D
#    else
#      error "bad QPHIXJ_Precision"
#    endif
#  else
#    if QPHIXJ_Precision == 'F'
#      if QPHIXJ_PrecisionInt != 1
#        error "inconsistent QPHIXJ_Precision='F' and QPHIXJ_PrecisionInt"
#      endif
#      define QPHIXJ_PrecisionLetter F
#    elif QPHIXJ_Precision == 'D'
#      if QPHIXJ_PrecisionInt != 2
#        error "inconsistent QPHIXJ_Precision='D' and QPHIXJ_PrecisionInt"
#      endif
#      define QPHIXJ_PrecisionLetter D
#    else
#      error "bad QPHIXJ_Precision"
#    endif
#  endif
#endif

#include "../include/qphixj/qphixj_int.h"
#include "../include/qphixj/qphixj_f3.h"
#include "../include/qphixj/qphixj_d3.h"


#endif /* QPHIXJ_H */
