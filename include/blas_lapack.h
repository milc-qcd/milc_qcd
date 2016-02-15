#ifndef _BLAS_LAPACK_H
#define _BLAS_LAPACK_H

/* Prototype declaration of blas and lapack routines for inc_eigcg.c */
void zcopy_(int *n, double_complex *x, int *incx, double_complex *y, int *incy);

void zhemm_(char *side, char *uplo, int *m, int *n, double_complex *alpha,
	    double_complex *A, int *ldA, double_complex *B, int *ldB, double_complex *beta,
	    double_complex *C, int *ldC);

void zgemm_(char *tarnsA, char *transB, int *m, int  *n, int *k, double_complex *alpha,
	   double_complex *A, int *ldA, double_complex *B, int *ldB, double_complex *beta,
	   double_complex *C, int *ldC);

void zheev_(char *jobz, char *uplo, int *n, double_complex *A, int *ldA, double *w,
	    double_complex *work, int *lwork, double *rwork, int *info);

void zheevx_(char *jobz, char *range, char *uplo, int *n, double_complex *A, int *ldA,
	     double *vl, double *vu, int *il, int *iu, double *abstol, int *m, double *w,
	     double_complex *Z, int *ldZ, double_complex *work, int *lwork, double *rwork,
	     int *iwork, int *ifail, int *info);

void zgeqrf_(int *m, int *n, double_complex *A, int *ldA, double_complex *tau,
	     double_complex *work, int *lwork, int *info);

void zungqr_(int *m, int *n, int *k, double_complex *A, int *ldA, double_complex *tau,
	     double_complex *work, int *lwork, int *info);

void zpotrf_(char *uplo, int *n, double_complex *A, int *ldA, int *info);

void zpotrs_(char *uplo, int *n, int *nrhs, double_complex *A, int *ldA, double_complex *B,
	     int *ldB, int *info);

#endif /* _BLAS_LAPACK_H */
