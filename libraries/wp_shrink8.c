/*****************  wp_shrink8.c  (in su3.a) ****************************
*									*
* Shrink a wilson vector in eight directions, producing eight		*
*  half_wilson_vectors.							*
* void wp_shrink_8dir(a,b,sign)						*
* wilson_vector *a; half_wilson_vector *b; 				*
* int sign;								*
* B1 <- (1 +- gamma_x)A,, projection					*
*  argument "sign" is sign of gamma matrix.				*
*  See wp_shrink.c for definitions of gamma matrices and eigenvectors.	*
*/
#include "../include/config.h"
#include "../include/complex.h"
#include "../include/su3.h"
#include "../include/dirs.h"

void wp_shrink_8dir( wilson_vector *a, half_wilson_vector *b, int sign) {
    wp_shrink( a,&(b[XUP]),XUP,sign);
    wp_shrink( a,&(b[YUP]),YUP,sign);
    wp_shrink( a,&(b[ZUP]),ZUP,sign);
    wp_shrink( a,&(b[TUP]),TUP,sign);
    wp_shrink( a,&(b[XDOWN]),XDOWN,sign);
    wp_shrink( a,&(b[YDOWN]),YDOWN,sign);
    wp_shrink( a,&(b[ZDOWN]),ZDOWN,sign);
    wp_shrink( a,&(b[TDOWN]),TDOWN,sign);
}
