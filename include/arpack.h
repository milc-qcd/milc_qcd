/*-------------------------------------
  FORTRAN related header files
  
  A M Abdel-Rehim (ama273@cornell.edu)

  Last updated: 06/14/2016
  -------------------------------------*/

#ifndef _FORTRAN_HEADER_H
#define _FORTRAN_HEADER_H

#include "complex.h"
#if 0
# include "mpi.h"
#endif

//using or not using underscrore for FORTRAN subroutines
#if (defined NOF77UNDERSCORE || defined NOF77_)
#define _AFT(s) s
#else
#define _AFT(s) s ## _
#endif

//ARPACK subroutines

//ARPACK initlog and finilog routines for printing the ARPACK log (same for serial and parallel version)
extern int _AFT(initlog) (int*, char*, int);
extern int _AFT(finilog) (int*);


//ARPACK driver routines (serial version) 
int _AFT(znaupd) (int *ido, char *bmat, int *n, char *which, int *nev, double *tol,
                         complex *resid, int *ncv, complex *v, int *ldv, 
                         int *iparam, int *ipntr, complex *workd, complex *workl, 
                         int *lworkl, double *rwork, int *info);

//                         int *lworkl, double *rwork, int *info, int bmat_size, int which_size );
			
int _AFT(zneupd) (int *comp_evecs, char *howmany, int *select, complex *evals, 
			 complex *v, int *ldv, complex* sigma, complex *workev, 
			 char *bmat, int *n, char *which, int *nev, double *tol, complex *resid, 
                         int *ncv, complex* v1, int *ldv1, int *iparam, int *ipntr, 
                         complex* workd, complex* workl, int *lworkl, double *rwork, int *info) ;

//                         int howmany_size, int bmat_size, int which_size);


int _AFT(mcinitdebug)(int*,int*,int*,int*,int*,int*,int*,int*);


#if 1
//PARPACK routines (parallel version)
int _AFT(pznaupd) (int *comm, int *ido, char *bmat, int *n, char *which, int *nev, double *tol,
                          complex* resid, int *ncv, complex *v, int *ldv, 
                          int *iparam, int *ipntr, complex* workd, complex* workl, 
                          int *lworkl, double *rwork, int *info);

int _AFT(pzneupd) (int *comm, int *comp_evecs, char *howmany, int *select, complex* evals, 
                          complex* v, int *ldv, complex* sigma, complex* workev, 
                          char *bmat, int *n, char *which, int *nev, double *tol, complex* resid, 
                          int *ncv, complex* v1, int *ldv1, int *iparam, int *ipntr, 
                          complex* workd, complex* workl, int *lworkl, double *rwork, int *info
                          );

			 //extern int _AFT(pmcinitdebug)(int*,int*,int*,int*,int*,int*,int*,int*);

#endif

#endif
