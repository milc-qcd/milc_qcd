#define DEPS_SCHUR 1.110223e-16

/* The Schur rotation of a hermitian matrix */
typedef struct 
{
  double c ;
  double_complex s ;
  int p ; 
  int q ;
} h_schur ;

typedef struct 
{
  int N ; 
  double_complex **M ;
} Matrix ;

Matrix AllocateMatrix(int N) ;
void deAllocate(Matrix *A) ;
double FrobeniusNorm(Matrix *A) ;
double OffDiag(Matrix *A) ;
void MatrixMult(Matrix *A, Matrix *B, Matrix *C) ;
void MatrixCopy(Matrix *A, Matrix *B) ;
void ZeroMatrix(Matrix *A) ;
void UnitMatrix(Matrix *A) ;
void HermitianConj(Matrix *A, Matrix *B) ;
void HermitianSchur(Matrix *A, int p, int q, h_schur *schur) ;
void rightMultiplySchur(Matrix *A, h_schur *schur) ;
void leftMultiplySchur(Matrix *A, h_schur *schur) ;
void Jacobi(Matrix *A, Matrix *V, double tolerance) ;
void sort_eigenvectors(Matrix *A, Matrix *V) ;

