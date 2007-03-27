/******************************************************************************
FILE:   instant.c
DESCRIPTION:
   Creates a binary lattice with a classical instanton configuration.

AUTHOR:     Mark Stephenson
DATE:       1-31-97
******************************************************************************/

/*-------------------------------------
   Include Files Directives
-------------------------------------*/
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include <errno.h>
#include "../../include/complex.h"
#include "../../include/su3.h"

/*-------------------------------------
   Global Constant Definitions
-------------------------------------*/
#define SUBLINK 2               /* number of intervals for interpolation */

/*-------------------------------------
   Function Prototypes
-------------------------------------*/
void init_pauli_matrices(su2_matrix *ptau[4]);
void check_unitarity(su2_matrix *psub);
void imbed_su2(su3_matrix *pmat, su2_matrix *psub);
void calc_su2_instanton_link(const int dir, const int xn[4],
                             const Real center[4], const Real rho,
                             su2_matrix *ptau[4], su2_matrix *psu2_link);
void transform(const int dir, const int xn[4],
               const Real center[4], su2_matrix *ptau[4],
               su2_matrix *psu2_link);
void mult_scalar(su2_matrix *pprod, const complex scalar, su2_matrix *pmat);
void mult_su2(su2_matrix *pprod, su2_matrix *pmat1, su2_matrix *pmat2);
void sum_su2(su2_matrix *psum, su2_matrix *pmat);

/******************************************************************************
FUNCTION:
   main
******************************************************************************/
int main(int argc, char *argv[])
{
   FILE *fp;                    /* file pointer, output lattice */
   char fname[FILENAME_MAX];    /* file name, output lattice */
   int version;                 /* code version number */
   int dim[4];                  /* lattice dimensions (order x,y,z,t) */
   Real beta;                  /* inverse coupling */
   Real mass;                  /* quark mass */
   su2_matrix tau[4];           /* Pauli matrix (tau[3] is the unit matrix) */
   su2_matrix *ptau[4];         /* pointer to Pauli matrix */
   int xn[4];                   /* lattice site indices (x,y,z,t) */
   int k;                       /* index over forward directions */
   su3_matrix link;             /* link matrix */
   su2_matrix su2_link;         /* SU(2) submatrix of link */
   Real center[4];             /* center of instanton */
   Real rho;                   /* instanton radius */

   /* get input parameters */
   if (argc == 7)
   {
      strcpy(fname, argv[1]);
      dim[0] = atoi(argv[2]);
      dim[1] = atoi(argv[3]);
      dim[2] = atoi(argv[4]);
      dim[3] = atoi(argv[5]);
      rho = (Real)atof(argv[6]);
   }
   else
   {
      printf("Enter output lattice filename      ");
      scanf("%s", fname);
      printf("Enter lattice dimensions           ");
      scanf("%d %d %d %d", &dim[0], &dim[1], &dim[2], &dim[3]);
      printf("Enter instanton radius             ");
      scanf("%lf", &rho);
   }

   /* open output file */
   if ((fp = fopen(fname, "r")) != NULL)
   {
      fprintf(stderr, "Output file %s already exists.\n", fname);
      fprintf(stderr, "Exiting program %s.\n", argv[0]);
      fclose(fp);
      exit(1);
   }
   if ((fp = fopen(fname, "wb")) == NULL)
   {
      fprintf(stderr, "Unable to open output file %s, ", fname);
      fprintf(stderr, "error %d.\n", errno);
      fprintf(stderr, "Exiting program %s.\n", argv[0]);
      exit(1);
   }

   /* create Pauli and unit matrices */
   for (k = 0; k < 4; k++)
   {
      ptau[k] = &tau[k];
   }
   init_pauli_matrices(ptau);

   /* set lattice parameters */
   version = 59354;
   beta = 6.0;
   mass = 0.0125;

   /* set instanton center */
   for (k = 0; k < 4; k++)
   {
      center[k] = (Real)dim[k] / 2.0 + 0.1;
   }

   /* write lattice file */
   fwrite(&version, sizeof(int), 1, fp);
   fwrite(dim, sizeof(int), 4, fp);
   fwrite(&beta, sizeof(Real), 1, fp);
   fwrite(&mass, sizeof(Real), 1, fp);
   for (xn[3] = 0; xn[3] < dim[3]; xn[3]++)
   {
      for (xn[2] = 0; xn[2] < dim[2]; xn[2]++)
      {
         for (xn[1] = 0; xn[1] < dim[1]; xn[1]++)
         {
            for (xn[0] = 0; xn[0] < dim[0]; xn[0]++)
            {
               for (k = 0; k < 4; k++)
               {
                  calc_su2_instanton_link(k, xn, center,
                                          rho, ptau, &su2_link);
                  transform(k, xn, center, ptau, &su2_link);
                  check_unitarity(&su2_link);
                  imbed_su2(&link, &su2_link);
                  fwrite(&link, sizeof(su3_matrix), 1, fp);
               }
            }
         }
      }
   }

   /* terminate program normally */
   fclose(fp);
   return 0;
}

/******************************************************************************
FUNCTION:
   init_pauli_matrices
DESCRIPTION:
   The tau matrices are set to the Pauli matrices and the identity matrix:
   tau[0] = sigma_1
   tau[1] = sigma_2
   tau[2] = sigma_3
   tau[3] = identity
******************************************************************************/
void init_pauli_matrices(su2_matrix *ptau[4])
{
   int k;                       /* index over forward directions */
   int p;                       /* index over matrix row */
   int q;                       /* index over matrix column */

   for (k =0; k < 4; k++)
   {
      for (p = 0; p < 2; p ++)
      {
         for (q = 0; q < 2; q++)
         {
            ptau[k]->e[p][q].real = 0.0;
            ptau[k]->e[p][q].imag = 0.0;
         }
      }
   }
   ptau[0]->e[0][1].real = 1.0;
   ptau[0]->e[1][0].real = 1.0;
   ptau[1]->e[0][1].imag = -1.0;
   ptau[1]->e[1][0].imag = 1.0;
   ptau[2]->e[0][0].real = 1.0;
   ptau[2]->e[1][1].real = -1.0;
   ptau[3]->e[0][0].real = 1.0;
   ptau[3]->e[1][1].real = 1.0;
}

/******************************************************************************
FUNCTION:
   imbed_su2
******************************************************************************/
void imbed_su2(su3_matrix *pmat, su2_matrix *psub)
{
   int p;                       /* index over matrix row */
   int q;                       /* index over matrix column */

   for (p = 0; p < 2; p++)
   {
      for (q = 0; q < 2; q++)
      {
         pmat->e[p][q].real = psub->e[p][q].real;
         pmat->e[p][q].imag = psub->e[p][q].imag;
      }
   }
   for (p = 0; p < 2; p++)
   {
      pmat->e[p][2].real = 0.0;
      pmat->e[p][2].imag = 0.0;
   }
   for (q = 0; q < 2; q++)
   {
      pmat->e[2][q].real = 0.0;
      pmat->e[2][q].imag = 0.0;
   }
   pmat->e[2][2].real = 1.0;
   pmat->e[2][2].imag = 0.0;
}

/******************************************************************************
FUNCTION:
   calc_su2_instanton_link
******************************************************************************/
void calc_su2_instanton_link(const int dir, const int xn[4],
                             const Real center[4], const Real rho,
                             su2_matrix *ptau[4], su2_matrix *psu2_link)
{
   Real x[4];                  /* position component */
   Real normalization;         /* normalization for unit 3-vector */
   Real shape;                 /* shape factor, averaged over link */
   Real theta;                 /* parameter in SU(2) exponential */
   complex scalar;              /* scalar factor */
   su2_matrix prod;             /* matrix product */
   su2_matrix *pprod;           /* pointer to matrix product */
   int k;                       /* index over forward directions */
   int p;                       /* index over matrix row */
   int q;                       /* index over matrix column */
   int m;                       /* index over link increments */

   /* initialize link to zero */
   for (p = 0; p < 2; p++)
   {
      for (q = 0; q < 2; q++)
      {
         psu2_link->e[p][q].real = 0.0;
         psu2_link->e[p][q].imag = 0.0;
      }
   }

   /* set position components of current site relative to instanton center */
   for (k = 0; k < 4; k++)
   {
      x[k] = (Real)xn[k] - center[k];
   }

   /* initialize shape factor */
   shape = 0.0;
   
   /* initialize pointer to matrix product */
   pprod = &prod;

   /* set instanton link value */
   switch (dir)
   {
   case 0:
      normalization = 1.0 / sqrt(x[1] * x[1] + x[2] * x[2] + x[3] * x[3]);
      for (m = 0; m < SUBLINK; m++)
      {
         shape += 1.0
           / ((x[0] + m / (SUBLINK - 1.0)) * (x[0] + m / (SUBLINK - 1.0))
              + x[1] * x[1] + x[2] * x[2] + x[3] * x[3] + rho * rho);
      }
      shape /= (Real)SUBLINK;
      theta = - shape / normalization;
      scalar.real = 0.0;
      scalar.imag = x[3] * sin(theta) * normalization;
      mult_scalar(pprod, scalar, ptau[0]);
      sum_su2(psu2_link, pprod);
      scalar.real = 0.0;
      scalar.imag = -x[2] * sin(theta) * normalization;
      mult_scalar(pprod, scalar, ptau[1]);
      sum_su2(psu2_link, pprod);
      scalar.real = 0.0;
      scalar.imag = x[1] * sin(theta) * normalization;
      mult_scalar(pprod, scalar, ptau[2]);
      sum_su2(psu2_link, pprod);
      scalar.real = cos(theta);
      scalar.imag = 0.0;
      mult_scalar(pprod, scalar, ptau[3]);
      sum_su2(psu2_link, pprod);
      break;
   case 1:
      normalization = 1.0 / sqrt(x[0] * x[0] + x[2] * x[2] + x[3] * x[3]);
      for (m = 0; m < SUBLINK; m++)
      {
         shape += 1.0
           / (x[0] * x[0]
              + (x[1] + m / (SUBLINK - 1.0)) * (x[1] + m / (SUBLINK - 1.0))
              + x[2] * x[2] + x[3] * x[3] + rho * rho);
      }
      shape /= (Real)SUBLINK;
      theta = - shape / normalization;
      scalar.real = 0.0;
      scalar.imag = x[2] * sin(theta) * normalization;
      mult_scalar(pprod, scalar, ptau[0]);
      sum_su2(psu2_link, pprod);
      scalar.real = 0.0;
      scalar.imag = x[3] * sin(theta) * normalization;
      mult_scalar(pprod, scalar, ptau[1]);
      sum_su2(psu2_link, pprod);
      scalar.real = 0.0;
      scalar.imag = -x[0] * sin(theta) * normalization;
      mult_scalar(pprod, scalar, ptau[2]);
      sum_su2(psu2_link, pprod);
      scalar.real = cos(theta);
      scalar.imag = 0.0;
      mult_scalar(pprod, scalar, ptau[3]);
      sum_su2(psu2_link, pprod);
      break;
   case 2:
      normalization = 1.0 / sqrt(x[0] * x[0] + x[1] * x[1] + x[3] * x[3]);
      for (m = 0; m < SUBLINK; m++)
      {
         shape += 1.0
           / (x[0] * x[0] + x[1] * x[1]
              + (x[2] + m / (SUBLINK - 1.0)) * (x[2] + m / (SUBLINK - 1.0))
              + x[3] * x[3] + rho * rho);
      }
      shape /= (Real)SUBLINK;
      theta = - shape / normalization;
      scalar.real = 0.0;
      scalar.imag = -x[1] * sin(theta) * normalization;
      mult_scalar(pprod, scalar, ptau[0]);
      sum_su2(psu2_link, pprod);
      scalar.real = 0.0;
      scalar.imag = x[0] * sin(theta) * normalization;
      mult_scalar(pprod, scalar, ptau[1]);
      sum_su2(psu2_link, pprod);
      scalar.real = 0.0;
      scalar.imag = x[3] * sin(theta) * normalization;
      mult_scalar(pprod, scalar, ptau[2]);
      sum_su2(psu2_link, pprod);
      scalar.real = cos(theta);
      scalar.imag = 0.0;
      mult_scalar(pprod, scalar, ptau[3]);
      sum_su2(psu2_link, pprod);
      break;
   case 3:
      normalization = 1.0 / sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);
      for (m = 0; m < SUBLINK; m++)
      {
         shape += 1.0 / (x[0] * x[0] + x[1] * x[1] + x[2] * x[2]
                         + (x[3] + m / (SUBLINK - 1.0))
                         * (x[3] + m / (SUBLINK - 1.0)) + rho * rho);
      }
      shape /= (Real)SUBLINK;
      theta = + shape / normalization;
      scalar.real = 0.0;
      scalar.imag = x[0] * sin(theta) * normalization;
      mult_scalar(pprod, scalar, ptau[0]);
      sum_su2(psu2_link, pprod);
      scalar.real = 0.0;
      scalar.imag = x[1] * sin(theta) * normalization;
      mult_scalar(pprod, scalar, ptau[1]);
      sum_su2(psu2_link, pprod);
      scalar.real = 0.0;
      scalar.imag = x[2] * sin(theta) * normalization;
      mult_scalar(pprod, scalar, ptau[2]);
      sum_su2(psu2_link, pprod);
      scalar.real = cos(theta);
      scalar.imag = 0.0;
      mult_scalar(pprod, scalar, ptau[3]);
      sum_su2(psu2_link, pprod);
      break;
   }
}

/******************************************************************************
FUNCTION:
   transform
******************************************************************************/
void transform(const int dir, const int xn[4],
               const Real center[4], su2_matrix *ptau[4],
               su2_matrix *psu2_link)
{
   Real x[4];                  /* position component */
   Real norm_current;          /* normalization factor at current site */
   Real norm_next;             /* normalization factor at next site */
   complex scalar;              /* scalar factor */
   su2_matrix prod;             /* matrix product */
   su2_matrix *pprod;           /* pointer to matrix product */
   su2_matrix gauge;            /* guage transformation matrix */
   su2_matrix *pgauge;          /* pointer to gauge transformation matrix*/
   su2_matrix inv;              /* inverse gauge transformation */
   su2_matrix *pinv;            /* pointer to inverse gauge transformation */
   su2_matrix temp;             /* holding matrix */
   su2_matrix *ptemp;           /* pointer to holding matrix */
   int k;                       /* index over forward directions */
   int p;                       /* index over matrix row */
   int q;                       /* index over matrix column */

   /* set position components of current site relative to instanton center */
   for (k = 0; k < 4; k++)
   {
      x[k] = (Real)xn[k] - center[k];
   }

   /* set normalization factor of current site */
   norm_current = sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2] + x[3] * x[3]);

   /* initialize pointers to matrices */
   pprod = &prod;
   pgauge = &gauge;
   pinv = &inv;
   ptemp = &temp;

   /* set instanton link value */
   switch (dir)
   {
   case 0:
      norm_next = sqrt((x[0] + 1.0) * (x[0] + 1.0) + x[1] * x[1] + x[2] * x[2]
                       + x[3] * x[3]);

      scalar.real = 0.0;
      scalar.imag = (x[0] + 1.0) / norm_next;
      mult_scalar(pgauge, scalar, ptau[0]);
      scalar.real = 0.0;
      scalar.imag = x[1] / norm_next;
      mult_scalar(pprod, scalar, ptau[1]);
      sum_su2(pgauge, pprod);
      scalar.real = 0.0;
      scalar.imag = x[2] / norm_next;
      mult_scalar(pprod, scalar, ptau[2]);
      sum_su2(pgauge, pprod);
      scalar.real = x[3] / norm_next;
      scalar.imag = 0.0;
      mult_scalar(pprod, scalar, ptau[3]);
      sum_su2(pgauge, pprod);

      scalar.real = 0.0;
      scalar.imag = -x[0] / norm_current;
      mult_scalar(pinv, scalar, ptau[0]);
      scalar.real = 0.0;
      scalar.imag = -x[1] / norm_current;
      mult_scalar(pprod, scalar, ptau[1]);
      sum_su2(pinv, pprod);
      scalar.real = 0.0;
      scalar.imag = -x[2] / norm_current;
      mult_scalar(pprod, scalar, ptau[2]);
      sum_su2(pinv, pprod);
      scalar.real = x[3] / norm_current;
      scalar.imag = 0.0;
      mult_scalar(pprod, scalar, ptau[3]);
      sum_su2(pinv, pprod);
      break;
   case 1:
      norm_next = sqrt(x[0] * x[0] + (x[1] + 1.0) * (x[1] + 1.0) + x[2] * x[2]
                       + x[3] * x[3]);

      scalar.real = 0.0;
      scalar.imag = x[0] / norm_next;
      mult_scalar(pgauge, scalar, ptau[0]);
      scalar.real = 0.0;
      scalar.imag = (x[1] + 1.0) / norm_next;
      mult_scalar(pprod, scalar, ptau[1]);
      sum_su2(pgauge, pprod);
      scalar.real = 0.0;
      scalar.imag = x[2] / norm_next;
      mult_scalar(pprod, scalar, ptau[2]);
      sum_su2(pgauge, pprod);
      scalar.real = x[3] / norm_next;
      scalar.imag = 0.0;
      mult_scalar(pprod, scalar, ptau[3]);
      sum_su2(pgauge, pprod);

      scalar.real = 0.0;
      scalar.imag = -x[0] / norm_current;
      mult_scalar(pinv, scalar, ptau[0]);
      scalar.real = 0.0;
      scalar.imag = -x[1] / norm_current;
      mult_scalar(pprod, scalar, ptau[1]);
      sum_su2(pinv, pprod);
      scalar.real = 0.0;
      scalar.imag = -x[2] / norm_current;
      mult_scalar(pprod, scalar, ptau[2]);
      sum_su2(pinv, pprod);
      scalar.real = x[3] / norm_current;
      scalar.imag = 0.0;
      mult_scalar(pprod, scalar, ptau[3]);
      sum_su2(pinv, pprod);
      break;
   case 2:
      norm_next = sqrt(x[0] * x[0] + x[1] * x[1] + (x[2] + 1.0) * (x[2] + 1.0)
                       + x[3] * x[3]);

      scalar.real = 0.0;
      scalar.imag = x[0] / norm_next;
      mult_scalar(pgauge, scalar, ptau[0]);
      scalar.real = 0.0;
      scalar.imag = x[1] / norm_next;
      mult_scalar(pprod, scalar, ptau[1]);
      sum_su2(pgauge, pprod);
      scalar.real = 0.0;
      scalar.imag = (x[2] + 1.0) / norm_next;
      mult_scalar(pprod, scalar, ptau[2]);
      sum_su2(pgauge, pprod);
      scalar.real = x[3] / norm_next;
      scalar.imag = 0.0;
      mult_scalar(pprod, scalar, ptau[3]);
      sum_su2(pgauge, pprod);

      scalar.real = 0.0;
      scalar.imag = -x[0] / norm_current;
      mult_scalar(pinv, scalar, ptau[0]);
      scalar.real = 0.0;
      scalar.imag = -x[1] / norm_current;
      mult_scalar(pprod, scalar, ptau[1]);
      sum_su2(pinv, pprod);
      scalar.real = 0.0;
      scalar.imag = -x[2] / norm_current;
      mult_scalar(pprod, scalar, ptau[2]);
      sum_su2(pinv, pprod);
      scalar.real = x[3] / norm_current;
      scalar.imag = 0.0;
      mult_scalar(pprod, scalar, ptau[3]);
      sum_su2(pinv, pprod);
      break;
   case 3:
      norm_next = sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]
                       + (x[3] + 1.0) * (x[3] + 1.0));

      scalar.real = 0.0;
      scalar.imag = x[0] / norm_next;
      mult_scalar(pgauge, scalar, ptau[0]);
      scalar.real = 0.0;
      scalar.imag = x[1] / norm_next;
      mult_scalar(pprod, scalar, ptau[1]);
      sum_su2(pgauge, pprod);
      scalar.real = 0.0;
      scalar.imag = x[2] / norm_next;
      mult_scalar(pprod, scalar, ptau[2]);
      sum_su2(pgauge, pprod);
      scalar.real = (x[3] + 1.0) / norm_next;
      scalar.imag = 0.0;
      mult_scalar(pprod, scalar, ptau[3]);
      sum_su2(pgauge, pprod);

      scalar.real = 0.0;
      scalar.imag = -x[0] / norm_current;
      mult_scalar(pinv, scalar, ptau[0]);
      scalar.real = 0.0;
      scalar.imag = -x[1] / norm_current;
      mult_scalar(pprod, scalar, ptau[1]);
      sum_su2(pinv, pprod);
      scalar.real = 0.0;
      scalar.imag = -x[2] / norm_current;
      mult_scalar(pprod, scalar, ptau[2]);
      sum_su2(pinv, pprod);
      scalar.real = x[3] / norm_current;
      scalar.imag = 0.0;
      mult_scalar(pprod, scalar, ptau[3]);
      sum_su2(pinv, pprod);
      break;
   }

   for (p = 0; p < 2; p++)
   {
      for (q = 0; q < 2; q++)
      {
         ptemp->e[p][q].real = psu2_link->e[p][q].real;
         ptemp->e[p][q].imag = psu2_link->e[p][q].imag;
      }
   }
   mult_su2(psu2_link, ptemp, pgauge);
   for (p = 0; p < 2; p++)
   {
      for (q = 0; q < 2; q++)
      {
         ptemp->e[p][q].real = psu2_link->e[p][q].real;
         ptemp->e[p][q].imag = psu2_link->e[p][q].imag;
      }
   }
   mult_su2(psu2_link, pinv, ptemp);
}

/******************************************************************************
FUNCTION:
   mult_scalar
******************************************************************************/
void mult_scalar(su2_matrix *pprod, const complex scalar, su2_matrix *pmat)
{
   int p;                       /* index over matrix row */
   int q;                       /* index over matrix column */

   for (p = 0; p < 2; p++)
   {
      for (q = 0; q < 2; q++)
      {
         CMUL(pmat->e[p][q], scalar, pprod->e[p][q])
      }
   }
}

/******************************************************************************
FUNCTION:
   sum_su2
******************************************************************************/
void sum_su2(su2_matrix *psum, su2_matrix *pmat)
{
   int p;                       /* index over matrix row */
   int q;                       /* index over matrix column */

   for (p = 0; p < 2; p++)
   {
      for (q = 0; q < 2; q++)
      {
         psum->e[p][q].real += pmat->e[p][q].real;
         psum->e[p][q].imag += pmat->e[p][q].imag;
      }
   }
}

/******************************************************************************
FUNCTION:
   mult_su2
******************************************************************************/
void mult_su2(su2_matrix *pprod, su2_matrix *pmat1, su2_matrix *pmat2)
{
   int p;                       /* index over matrix row */
   int q;                       /* index over matrix column */
   int k;                       /* summation index */
   complex term;                /* term in row times column */

   for (p = 0; p < 2; p++)
   {
      for (q = 0; q < 2; q++)
      {
         pprod->e[p][q].real = 0.0;
         pprod->e[p][q].imag = 0.0;
         for (k = 0; k < 2; k++)
         {
            CMUL(pmat1->e[p][k], pmat2->e[k][q], term);
            CSUM(pprod->e[p][q], term);
         }
      }
   }
}

/******************************************************************************
FUNCTION:
   check_unitarity
******************************************************************************/
void check_unitarity(su2_matrix *psub)
{
   complex sum;
   complex conj;
   complex prod;
   Real tol;

   tol = 0.0001;
   sum.real = 0.0;
   sum.imag = 0.0;
   CONJG(psub->e[0][0], conj);
   CMUL(psub->e[0][0], conj, prod);
   CSUM(sum, prod);
   CONJG(psub->e[1][0], conj);
   CMUL(psub->e[1][0], conj, prod);
   CSUM(sum, prod);
   if (sum.real < 1.0 - tol || sum.real > 1.0 + tol
       || sum.imag < -tol || sum.imag > tol)
   {
      printf("%g\t%g\n", sum.real, sum.imag);
   }
}
