/********************** phases.c *********************************/
/* MIMD version 7 */

#include <string.h>
#include "../include/complex.h"
#include "../include/su3.h"

/* Map an ASCIi label to the corresponding integer */
/* Follows gammatypes.h convention for encoding the phase */

static const char *phaselabel[4] = { "1", "i", "-1", "-i" };
static int phase[4] = { 0, 1, 2, 3 };  

int decode_phase(char const *label){
  int i;
  for(i = 0; i < 4; i++){
    if(strcmp(label,phaselabel[i]) == 0)return phase[i];
  }
  return -1;  /* Error condition */
}

/* Multiply b = a times phase using the encoding above */

void mult_c_by_phase(complex *a, complex *b, int ph){

    switch(ph){
    case 0:
      *b =           *a;
      break;
    case 1:
      TIMESPLUSI(    *a, *b);
      break;
    case 2:
      TIMESMINUSONE( *a, *b);
      break;
    case 3:
      TIMESMINUSI(   *a, *b);
    }
}



