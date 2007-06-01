#ifndef _GAMMATYPES_H
#define _GAMMATYPES_H
/*********************  gammatypes.h   **********************************
*									*
*  MIMD version 7 							*
*                                                                       *
*	A gamma matrix encoding scheme					*
*                                                                       *
*/

/*  DeGrand-Rossi convention

 gamma(XUP) 
 	    0  0  0  i
            0  0  i  0
            0 -i  0  0
           -i  0  0  0

 gamma(YUP)
 	    0  0  0 -1
            0  0  1  0
            0  1  0  0
           -1  0  0  0

 gamma(ZUP)
 	    0  0  i  0
            0  0  0 -i
           -i  0  0  0
            0  i  0  0

 gamma(TUP)
 	    0  0  1  0
            0  0  0  1
            1  0  0  0
            0  1  0  0

 gamma(FIVE) 
 	    1  0  0  0
            0  1  0  0
            0  0 -1  0
            0  0  0 -1
*/

/* Note: the order of the first five is fixed for compatibility with
   previous versions of mult_by_gamma */

enum gammatype {  GX, GY, GZ, GT, G5, 
		  GYZ, GZX, GXY, GXT, GYT, GZT, 
		  G5X, G5Y, G5Z, G5T, G1,
		  MAXGAMMA };

/* This encoding of gamma matrices works for representations in which
   matrix is just a permutation matrix with a phase change
   corresponding to the fourth roots of 1.  See the convention above.

   The matrix gamma_{r,c} (row index r, column index c) is encoded
   so that gamma.row[r].column = c is the only nonzero column in row r
   and gamma.row[r].phase specifies the phase of that matrix element

   */

typedef struct {
  int column;         /* encodes column with nonzero entry */
  int phase;          /* encodes phase of matrix element: 0,1,2,3 ->
                         1,i,-1,-i */
} gamma_row;

typedef struct {
  gamma_row row[4];
} gamma_matrix_t;

#endif /* _GAMMATYPES_H */
