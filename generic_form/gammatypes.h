#ifndef _GAMMATYPES_H
#define _GAMMATYPES_H
/************* gammatypes.h  **************************/

#ifdef CONTROL
#define EXTERN 
#else
#define EXTERN extern
#endif

/* 

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

enum gammatype {  GX,GY,GZ,GT,G5,GYZ,GZX,GXY,GXT,GYT,GZT,G5X,G5Y,G5Z,G5T,G1,
		MAXGAMMA };

/* This encoding of gamma matrices works for representations in which
   matrix is just a permutation matrix with a phase change
   corresponding to the fourth roots of 1.  See the BJ-Drell convention above.

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
} gamma_matrix;

/* First four gamma matrices are initialized according to the above
   conventions the ones for gamma_x gamma_y gamma_z gamma_t are used
   to generate the rest */

#ifdef CONTROL
 gamma_matrix gamma_mat[MAXGAMMA] = {
    { 3,1 ,  2,1 ,  1,3 , 0,3  },    /* gamma_x     */
    { 3,2 ,  2,0 ,  1,0 , 0,2  },    /* gamma_y     */
    { 2,1 ,  3,3 ,  0,3 , 1,1  },    /* gamma_z     */
    { 2,0 ,  3,0 ,  0,0 , 1,0  }     /* gamma_t     */
};
#else
  EXTERN gamma_matrix gamma_mat[MAXGAMMA] ;
#endif


EXTERN int gamma_initialized ;

typedef struct {
  int gin;         /* Source gamma matrix */
  int gout;        /* Sink gamma matrix */
  int oper;        /* Operator identification */
} gamma_corr;
#endif /* _GAMMATYPES_H */
