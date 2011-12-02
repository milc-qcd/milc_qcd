#ifndef _DEFINES_H
#define _DEFINES_H

/* Compiler macros common to all targets in this application */

#ifdef RANDOM
#define SITERAND	/* Use site-based random number generators */
#endif
#define JACOBI_TOL 1.110223e-16

#define GAUGE_FIX_TOL 1.0e-7	/* For gauge fixing */
#define MAX_SPECTRUM_REQUEST 512

#define MAX_MASSES 6

#define HZERO 1
#define HOVERLAP 2
#define HTEST 3 /* for trial action */
#define FN 
/* definitions to switch setup_links */
#define SIMPLE 0
#define HIGHP 1
#define LOWP 2
 
/*number of paths in the hypercube for 1 and gamma-mu paths */
/* #define NLINK 40 */	/* hypercubic action */
#define NLINK 16	/* planar action (use 4 for nearest neighbor!) */

#include "../include/complex.h"
#include "../include/su3.h"

/* does NOT use generic_clover routines! */
typedef struct { complex tr[2][15]; } triangular;
typedef struct { Real di[2][6]; } diagonal;
typedef struct { int chiral; Real coeff; Real mass; 
                 wilson_vector* src; } chiral_src;



/*********************** fat links stuff ****************************/

#ifdef NHYP

/*** number of smearing levels ***/

/* Set the number of smearing levels. Legal values are 1, 2, or 3
   (true NHYP is 3, "stout" is 1).
   NOTE: to save coding, the three alpha_smear parameters must
   always be supplied, either by hard coding below or by reading them from
   the infile.  But if SMEAR_LEVEL=n with n<3, then only the first n values
   of the alpha_smear parameters are used.
*/
#define SMEAR_LEVEL 3

/*** hard coding smearing parameters ***/

/* Comment out the following line
   if you want to read the smearing parameters from the infile.
   Otherwise, they are hard-coded below
*/
#define HARD_CODE_SMEAR
#ifdef HARD_CODE_SMEAR
#ifdef CONTROL
Real alpha_smear[3]={0.75, 0.6, 0.3};
#else
extern Real alpha_smear[3];
#endif /* CONTROL */
#endif /* HARD_CODE_SMEAR */

/*** constants for numerical stability ***/

/* IR regulator:  Q = Omega^dag Omega + IR_STAB
   This slightly changes the definition of the nhyp link. Fine.
   Since we're adding a constant, any derivative of Q is unchanged.
   routines block_nhyp1,2,3() in block_nhyp.c
            Sigma_update1() in force_nhyp.c
*/
#define IR_STAB 1.0e-6

/* calculation of Q^{-1/2}:  bypass R/(-S)^{3/2}=0/0
   routine compute_fhb() in generic_nhyp/nhyp.c
   neighborhood of 0 where we use approximate R and S for u0, u1, p
*/
#define EPS_SQ 1.0e-5

#endif /* NHYP*/
/******************** END fat links stuff ***************************/

#endif /* _DEFINES_H */
#ifndef _FIELD_ALLOC_H
#define _FIELD_ALLOC_H

/* macros to allocate space for fields */

#define FIELD_ALLOC(name,typ)                                \
    name = (typ *)malloc(sites_on_node*sizeof(typ));         \
    if (name==NULL){                                         \
        printf("NODE %d: FIELD_ALLOC failed\n",this_node);   \
        terminate(1);                                        \
    }

#define FIELD_ALLOC_VEC(name,typ,size)                              \
{                                                                   \
    int ifield;                                                     \
    for(ifield=0;ifield<size;ifield++) {                            \
        name[ifield] = (typ *)malloc(sites_on_node*sizeof(typ));    \
        if (name[ifield]==NULL){                                    \
            printf("NODE %d: FIELD_ALLOC_VEC failed\n",this_node);  \
            terminate(1);                                           \
        }                                                           \
    }                                                               \
}

#define FIELD_ALLOC_MAT(name,typ,sizea,sizeb)                              \
{                                                                          \
    int ifield,jfield;                                                     \
    for(ifield=0;ifield<sizea;ifield++)                                    \
    for(jfield=0;jfield<sizeb;jfield++) {                                  \
        name[ifield][jfield] = (typ *)malloc(sites_on_node*sizeof(typ));   \
        if (name[ifield][jfield]==NULL){                                   \
            printf("NODE %d: FIELD_ALLOC_MAT_OFFDIAG failed\n",this_node); \
            terminate(1);                                                  \
        }                                                                  \
    }                                                                      \
}

#define FIELD_ALLOC_MAT_OFFDIAG(name,typ,size)                             \
{                                                                          \
    int ifield,jfield;                                                     \
    for(ifield=0;ifield<size;ifield++)                                     \
    for(jfield=0;jfield<size;jfield++)                                     \
    if(ifield != jfield) {                                                 \
        name[ifield][jfield] = (typ *)malloc(sites_on_node*sizeof(typ));   \
        if (name[ifield][jfield]==NULL){                                   \
            printf("NODE %d: FIELD_ALLOC_MAT_OFFDIAG failed\n",this_node); \
            terminate(1);                                                  \
        }                                                                  \
    }                                                                      \
}

#endif /* _FIELD_ALLOC_H */


/* used in nhyp routines...*/
#ifdef SF
#define FORALLDYNLINKS(i,s,mu) \
    FORALLSITES(i,s) if((mu)==TUP || s->t>0)
#else
#define FORALLDYNLINKS(i,s,mu) \
    FORALLSITES(i,s)
#endif

#ifdef SF
#define CHOOSE_NBR(i,s,nu,chooseme,which_pt) \
  ( ((s)->t==(nt-1) && (nu)==TUP) ?          \
    (char *)&(chooseme)                      \
    :                                        \
    (gen_pt[which_pt][i])                    \
  )
#else
#define CHOOSE_NBR(i,s,nu,chooseme,which_pt) \
  (gen_pt[which_pt][i])
#endif
