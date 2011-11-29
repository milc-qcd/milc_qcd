/******************* make_clov2.c **************************************/

/* MIMD version 7 */
/* version of 1/4/95 by UMH */
/* 1/21/00 combined with Schroedinger functional version - UMH */
/* 4/24/07 Treat Clov_c = 0 case trivially, so we can do Wilson this way */

/* Prepare the "clover" term */

#include "generic_clover_includes.h"
#include "../include/prefetch.h"
#define FETCH_UP 1
#define LOOPEND
#include "../include/loopend.h"

/******************* global_clov *********************************/
/* The global clover term.  Used when we are working with only one
   clover term at a time */

// static clover *global_clov = NULL;

/******************* create_clov ***************************************/
/* Allocate space for a clover term */

clover *create_clov(void){
  
  clover *my_clov = (clover *)malloc(sizeof(clover));
  
  my_clov->clov          = NULL;
  my_clov->clov_diag     = NULL;
  my_clov->clov_raw      = NULL;
  my_clov->clov_diag_raw = NULL;
  my_clov->valid_clov    = 0;
  my_clov->trlogA        = 0.;

  return my_clov;
}


/***************** recompute_clov_parity ************************************/
/* Recomputes the clover field but only for sites with the specified parity */

static void recompute_clov_parity(clover *my_clov, int parity, Real Clov_c)
{
  
  triangular *clov, *clov_raw;
  diagonal *clov_diag, *clov_diag_raw;
  register int i;
  register site *s;
  register int j,k;
  char myname[] = "recompute_clov_parity";

  if(my_clov == NULL){
    printf("%s(%d): Unexpected null clover structure\n", myname, this_node);
  }
  if(my_clov->clov == NULL || my_clov->clov_diag == NULL ||
     my_clov->clov_raw == NULL || my_clov->clov_diag_raw == NULL ){
    printf("%s(%d): Unexpected null clover term\n", myname, this_node);
    terminate(1);
  }

  /* Local defs for convenience */
  clov = my_clov->clov;
  clov_diag = my_clov->clov_diag;
  clov_raw = my_clov->clov_raw;
  clov_diag_raw = my_clov->clov_diag_raw;

  /* Compute clov = 1 + Clov_c * clov_raw */
  FORALLSITESDOMAIN(i,s){
    for(j=0;j<6;j++){
      clov_diag[i].di[0][j] = 1.0 + Clov_c*clov_diag_raw[i].di[0][j];
      clov_diag[i].di[1][j] = 1.0 + Clov_c*clov_diag_raw[i].di[1][j];
    }
    for(k=0;k<15;k++){
      CMULREAL(clov_raw[i].tr[0][k],Clov_c,clov[i].tr[0][k]);
      CMULREAL(clov_raw[i].tr[1][k],Clov_c,clov[i].tr[1][k]);
    }
  }
} /* recompute_clov_parity */

/******************* compute_clov ***************************************/
/* Computes the "raw" clover field as well as the clover field */
/* clov = 1 + clov_raw * Clov_c */

void compute_clov(clover *my_clov, Real Clov_c)
{
  
  triangular *clov_raw;
  diagonal *clov_diag_raw;
  register int i;
  register site *s;
  
  register int j,k,jk,jk2;
  register complex ctmp;
  su3_matrix *f_mn;
  char myname[] = "compute_clov";

  my_clov->Clov_c = Clov_c;

  /* Don't compute if trivial */
  if(Clov_c == 0)return;

  /* Allocate space for clover term if not already allocated */
  if(my_clov->clov == NULL){
    my_clov->clov = (triangular *)malloc(sites_on_node*sizeof(triangular));
    my_clov->valid_clov    = 0;
  }

  if(my_clov->clov_diag == NULL){
    my_clov->clov_diag = (diagonal *)malloc(sites_on_node*sizeof(diagonal));
    my_clov->valid_clov    = 0;
  }

  if(my_clov->clov_raw == NULL){
    my_clov->clov_raw = (triangular *)malloc(sites_on_node*sizeof(triangular));
    my_clov->valid_clov    = 0;
  }

  if(my_clov->clov_diag_raw == NULL){
    my_clov->clov_diag_raw = 
      (diagonal *)malloc(sites_on_node*sizeof(diagonal));
    my_clov->valid_clov    = 0;
  }
    
  if(my_clov->clov == NULL || my_clov->clov_diag == NULL ||
     my_clov->clov_raw == NULL || my_clov->clov_diag_raw == NULL ){
    printf("compute_clov(%d): malloc failed\n",this_node);
    terminate(1);
  }

  /* Local defs for convenience */

  clov_raw = my_clov->clov_raw;
  clov_diag_raw = my_clov->clov_diag_raw;

  /* Allocate temporary */
  f_mn = (su3_matrix *)malloc(sizeof(su3_matrix)*sites_on_node);
  if(f_mn == NULL){
    printf("%s(%d): Can't malloc f_mn\n",myname,this_node);
    terminate(1);
  }
  
  
/* The clover term
*                                         ( X | 0 )
* sum_{mu<nu} sigma_{mu,nu} F_{mu,nu} = ( ----- )
*                                         ( 0 | Y )
*
* is block diagonal in the Weyl basis used, and hermitian. Together with the 1|
* we will store it in a lower complex triangular (block) matrix (without
* diagonal) and a real diagonal. The blocks go over Dirac indices 0,1 and 2,3.
*
*   With
*
*   sigma_{mu,nu} = i/2 [gamma_{mu}, gamma_{nu}]
*
*   and with our gamma matrices (see libraries/mb_gamma.c)
*
*   1234 = XYZT    and   Pauli sigma_1, sigma_2, sigma_3
*
*   sigma_23 =  diag( sigma_1,  sigma_1)
*   sigma_31 =  diag(-sigma_2, -sigma_2)
*   sigma_12 =  diag( sigma_3,  sigma_3)
*   sigma_14 =  diag(-sigma_1,  sigma_1)
*   sigma_24 =  diag( sigma_2, -sigma_2)
*   sigma_34 =  diag(-sigma_3,  sigma_3)
*
*
* Here X = sigma_1 (F_23 - F_14) + sigma_2 (F_13 + F_24) + sigma_3 (F_12 - F_34)
* and  Y = sigma_1 (F_23 + F_14) + sigma_2 (F_13 - F_24) + sigma_3 (F_12 + F_34)
*
* This has to be multiplied by the "clover coefficient" and subtracted from 1|.
*
* Note that F_mn = f_mn/i = -i f_mn!
*
* (The indices above start from 1, in the program they start from 0,
* i.e. subtract -1 from the above) */

/* The clover term R is a hermitian 12 x 12 matrix in color and Dirac spin.
   If we use ab for color and ij for spin, with a = (0,1,2) and
   i = (0,1,2,3), then we store the matrix R(ai,bj) as follows with
   the color index varying most rapidly:

   Upper left block, showing blocks of SU(3) matrices

   [ 1 + ci(f_01 - f_23)                ci(f_12 - f_03) + c(f_02 + f_13) ]
   [ ci(f_12 - f_03) - c(f_02 + f_13)   1 - ci(f_01 - f_23)                ]

   Lower right block, showing blocks of SU(3) matrices

   [ 1 + ci(f_01 + f_23)                 ci(f_12 + f_03) + c(f_02 - f_13)  ]
   [ ci(f_12 + f_03) - c(f_02 - f_13)    1 - ci(f_01 + f_23)               ]

   This hermitian matrix is decomposed into a lower triangular matrix
   and its diagonals.  The diagonals are real.

   clov_diag[site].di[0][0-5] holds the diagonals from the upper left
   clov_diag[site].di[1][0-5] holds the diagonals from the lower right

   clov_diag[site].di[0][a]   = R(a0,a0) = 1 - c(Im f_01 - Im f_23)(a,a)
   clov_diag[site].di[0][a+3] = R(a1,a1) = 1 + c(Im f_01 - Im f_23)(a,a)
   clov_diag[site].di[1][a]   = R(a2,a2) = 1 - c(Im f_01 + Im f_23)(a,a)
   clov_diag[site].di[1][a+3] = R(a3,a3) = 1 + c(Im f_01 + Im f_23)(a,a)

   clov[site].tr[0][jk]:
      Holds the lower triangular elements of the upper left block in
      the pattern

        jk = ...

   [          |        ]      [                      |                     ]
   [  0       |        ]      [  1 + ci(f_01 - f_23) | ci(f_12 - f_03)     ]
   [  1  2    |        ]      [                      |   + c(f_02 + f_13)  ]
   ---------------------   =  ----------------------------------------------
   [  3  4  5 |        ]      [			     |		           ]
   [  6  7  8 |  9     ]      [ ci(f_12 - f_03)      | 1 - ci(f_01 - f_23) ]
   [ 10 11 12 | 13 14  ]      [    - c(f_02 + f_13)  |		           ]

   clov[site].tr[1][jk]:
      Holds the lower triangular elements of the lower right block in
      the pattern

        jk = ...

   [          |        ]      [                      |                     ]
   [  0       |        ]      [  1 + ci(f_01 + f_23) | ci(f_12 + f_03)     ]
   [  1  2    |        ]      [                      |   + c(f_02 - f_13)  ]
   ---------------------   =  ----------------------------------------------
   [  3  4  5 |        ]      [			     |		           ]
   [  6  7  8 |  9     ]      [ ci(f_12 + f_03)      | 1 - ci(f_01 + f_23) ]
   [ 10 11 12 | 13 14  ]      [    - c(f_02 - f_13)  |		           ]
   
   */

  if(!my_clov->valid_clov){
    f_mu_nu(f_mn, 0, 1);
    jk=0;
    for(j=0;j<3;j++){
      FORALLSITESDOMAIN(i,s){
	clov_diag_raw[i].di[0][j]   =  -(f_mn+i)->e[j][j].imag;
	clov_diag_raw[i].di[0][j+3] =   (f_mn+i)->e[j][j].imag;
	clov_diag_raw[i].di[1][j]   =  -(f_mn+i)->e[j][j].imag;
	clov_diag_raw[i].di[1][j+3] =   (f_mn+i)->e[j][j].imag;
      }
      jk2=(j+3)*(j+2)/2+3;
      for(k=0;k<j;k++){
	FORALLSITESDOMAIN(i,s){
	  TIMESPLUSI( (f_mn+i)->e[j][k],  clov_raw[i].tr[0][jk]);
	  TIMESMINUSI( (f_mn+i)->e[j][k], clov_raw[i].tr[0][jk2]);
	  TIMESPLUSI( (f_mn+i)->e[j][k],  clov_raw[i].tr[1][jk]);
	  TIMESMINUSI( (f_mn+i)->e[j][k], clov_raw[i].tr[1][jk2]);
	}
	jk++;
	jk2++;
      }
    }
    
    f_mu_nu(f_mn, 2, 3);
    jk=0;
    for(j=0;j<3;j++){
      FORALLSITESDOMAIN(i,s){
	clov_diag_raw[i].di[0][j] +=	 (f_mn+i)->e[j][j].imag;
	clov_diag_raw[i].di[0][j+3] -= (f_mn+i)->e[j][j].imag;
	clov_diag_raw[i].di[1][j] -=   (f_mn+i)->e[j][j].imag;
	clov_diag_raw[i].di[1][j+3] += (f_mn+i)->e[j][j].imag;
      }
      jk2=(j+3)*(j+2)/2+3;
      for(k=0;k<j;k++){
	FORALLSITESDOMAIN(i,s){
	  TIMESMINUSI( (f_mn+i)->e[j][k], ctmp);
	  CADD( clov_raw[i].tr[0][jk], ctmp, clov_raw[i].tr[0][jk]);
	  CSUB( clov_raw[i].tr[0][jk2], ctmp, clov_raw[i].tr[0][jk2]);
	  CSUB( clov_raw[i].tr[1][jk], ctmp, clov_raw[i].tr[1][jk]);
	  CADD( clov_raw[i].tr[1][jk2], ctmp, clov_raw[i].tr[1][jk2]);
	}
	jk++;
	jk2++;
      }
    }
    
    f_mu_nu(f_mn, 1, 2);
    for(j=0;j<3;j++){
      jk=(j+3)*(j+2)/2;
      for(k=0;k<3;k++){
	FORALLSITESDOMAIN(i,s){
	  TIMESPLUSI( (f_mn+i)->e[j][k], clov_raw[i].tr[0][jk]);
	  TIMESPLUSI( (f_mn+i)->e[j][k], clov_raw[i].tr[1][jk]);
	}
	jk++;
      }
    }
    
    f_mu_nu(f_mn, 0, 3);
    for(j=0;j<3;j++){
      jk=(j+3)*(j+2)/2;
      for(k=0;k<3;k++){
	FORALLSITESDOMAIN(i,s){
	  TIMESMINUSI( (f_mn+i)->e[j][k], ctmp);
	  CADD( clov_raw[i].tr[0][jk], ctmp, clov_raw[i].tr[0][jk]);
	  CSUB( clov_raw[i].tr[1][jk], ctmp, clov_raw[i].tr[1][jk]);
	}
	jk++;
      }
    }
    
    f_mu_nu(f_mn, 0, 2);
    for(j=0;j<3;j++){
      jk=(j+3)*(j+2)/2;
      for(k=0;k<3;k++){
	FORALLSITESDOMAIN(i,s){
	  CSUB( clov_raw[i].tr[0][jk], (f_mn+i)->e[j][k], clov_raw[i].tr[0][jk]);
	  CSUB( clov_raw[i].tr[1][jk], (f_mn+i)->e[j][k], clov_raw[i].tr[1][jk]);
	}
	jk++;
      }
    }
    
    f_mu_nu(f_mn, 1, 3);
    for(j=0;j<3;j++){
      jk=(j+3)*(j+2)/2;
      for(k=0;k<3;k++){
	FORALLSITESDOMAIN(i,s){
	  CSUB( clov_raw[i].tr[0][jk], (f_mn+i)->e[j][k], clov_raw[i].tr[0][jk]);
	  CADD( clov_raw[i].tr[1][jk], (f_mn+i)->e[j][k], clov_raw[i].tr[1][jk]);
	}
	jk++;
      }
    }
    my_clov->valid_clov = 1;
  } /* valid_clov*/

  free(f_mn);

  /* Compute clov = 1 + Clov_c * clov_raw */
  recompute_clov_parity(my_clov, EVENANDODD, Clov_c);

} /* compute_clov */

/******************* compute_clovinv ***********************************/

/* MIMD version 7 */

/* Modifications
   added parity selection CD
   version of 1/3/95 by UMH 
*/

/* Invert the "clover" term on a specified sublattice and return tr(log(A)) */
/* (tr log needed for clover_dynamical) */

void compute_clovinv(clover *my_clov, int parity) {
  
  triangular *clov = my_clov->clov;
  diagonal *clov_diag = my_clov->clov_diag;
  Real Clov_c = my_clov->Clov_c;
  register int i;
  register site *s;
  
  register int b,j,k,kj,l,lk,lj,jl;
  Real f_diag[6];
  complex v1[6],ctmp,sum;
  double trlogA;
  
  if(Clov_c == 0)return;
  if(!my_clov->valid_clov){
    printf("compute_clovinv(%d): Can't invert an invalid clover term\n",
	   this_node);
    terminate(1);
  }

  /* We need to start from a fresh clover term because this module
     changes it */

  recompute_clov_parity(my_clov, parity, Clov_c);

  /* Take the inverse on the odd sublattice for each of the 2 blocks */
  trlogA = (double)0.0;
  FORSOMEPARITYDOMAIN(i,s,parity)for(b=0;b<2;b++){
    
    /* Cholesky decompose. This can be done in place, except that
       we will need a temporary diagonal.
       Uses algorithm 4.2.2 of "Matrix Computations" by Gene H. Golub
       and Charles F. Van Loan (John Hopkins, 1989) */
    for(j=0;j<6;j++){
      
      clov_diag[i].di[b][j] = sqrt((double)(clov_diag[i].di[b][j]));
      f_diag[j] = 1.0 / (clov_diag[i].di[b][j]);
      for(k=j+1;k<6;k++){
	kj=k*(k-1)/2+j;
	CMULREAL(clov[i].tr[b][kj], f_diag[j], clov[i].tr[b][kj]);
      }
      
      for(k=j+1;k<6;k++){
	kj=k*(k-1)/2+j;
	CMUL_J(clov[i].tr[b][kj], clov[i].tr[b][kj], ctmp);
	clov_diag[i].di[b][k] -= ctmp.real;
	for(l=k+1;l<6;l++){
	  lj=l*(l-1)/2+j;
	  lk=l*(l-1)/2+k;
	  CMUL_J(clov[i].tr[b][lj], clov[i].tr[b][kj], ctmp);
	  CSUB(clov[i].tr[b][lk], ctmp, clov[i].tr[b][lk]);
	}
      }
    }
    
    /* Accumulate trlogA */
    for(j=0;j<6;j++) 
      trlogA += (double)2.0 * log((double)(clov_diag[i].di[b][j]));
    
    /* Now use forward and backward substitution to construct inverse */
    for(k=0;k<6;k++){
      for(l=0;l<k;l++) v1[l] = cmplx(0.0, 0.0);
      
      /* Forward substitute */
      v1[k] = cmplx(f_diag[k], 0.0);
      for(l=k+1;l<6;l++){
	sum = cmplx(0.0, 0.0);
	for(j=k;j<l;j++){
	  lj=l*(l-1)/2+j;		    
	  CMUL(clov[i].tr[b][lj], v1[j], ctmp);
	  CSUB(sum, ctmp, sum);
	}
	CMULREAL(sum, f_diag[l], v1[l]);
      }
      
      /* Backward substitute */
      CMULREAL(v1[5], f_diag[5], v1[5]);
      for(l=4;l>=k;l--){
	sum = v1[l];
	for(j=l+1;j<6;j++){
	  jl=j*(j-1)/2+l;
	  CMULJ_(clov[i].tr[b][jl], v1[j], ctmp);
	  CSUB(sum, ctmp, sum);
	}
	CMULREAL(sum, f_diag[l], v1[l]);
      }
      
      /* Overwrite column k */
      clov_diag[i].di[b][k] = v1[k].real;
      for(l=k+1;l<6;l++){
	lk=l*(l-1)/2+k;
	clov[i].tr[b][lk] = v1[l];
      }
    }
  } END_LOOP;
      
  g_doublesum( &trlogA );

  my_clov->trlogA = trlogA;

} /* compute_clovinv */

/******************* copy_site ***********************************/
static void copy_site(
  field_offset src,   /* type wilson_vector RECAST AS wilson_block_vector */
  field_offset dest,  /* type wilson_vector RECAST AS wilson_block_vector */
  int parity
  )
{
  int i;
  site *s;
  FORSOMEPARITYDOMAIN(i,s,parity){
    copy_wvec( (wilson_vector *)F_PT(s,src), (wilson_vector *)F_PT(s,dest) );
  } END_LOOP
}
/******************* mult_this_ldu_site *********************************/

/* version of 12/29/94 by UMH */
/* CD 06/02/98 loop rearrangements, prefetch and macros as suggested by DT */

/* Multiply a Wilson vector (spin,color) with a block-diagonal hermition
   matrix stored as a complex lower triangular matrix (without diagonal)
   and a real diagonal. The blocks are spin indices 0,1 and 2,3. */

/* To simplify the code with the above block structure, introduce
   the somewhat dirty structure, equivalent to a wilson_vector: */
typedef struct { complex b[2][6]; } wilson_block_vector;

/* c += a*b */
#define CMULSUM(a,b,c) { (c).real += (a).real*(b).real - (a).imag*(b).imag; \
				(c).imag += (a).real*(b).imag + (a).imag*(b).real; }
/*  c += conj(a)*b */
#define CMULJ_SUM(a,b,c) { (c).real += (a).real*(b).real + (a).imag*(b).imag; \
				(c).imag += (a).real*(b).imag - (a).imag*(b).real; }

void mult_this_ldu_site(
  clover *my_clov,
  field_offset src,   /* type wilson_vector RECAST AS wilson_block_vector */
  field_offset dest,  /* type wilson_vector RECAST AS wilson_block_vector */
  int parity
  )
{
  triangular *clov = my_clov->clov;
  diagonal *clov_diag = my_clov->clov_diag;

  register int i;
  register site *s;
  register int b,j,k,jk,kj;
  register complex ctmp;
  
  if(my_clov->Clov_c == 0){
    copy_site(src, dest, parity);
    return;
  }

  if(!my_clov->valid_clov){
    printf("mult_this_ldu_site(%d): Can't use an invalid clover term\n",
	   this_node);
    terminate(1);
  }

  FORSOMEPARITYDOMAIN(i,s,parity){
    if( i < loopend-FETCH_UP ){
      prefetch_W( (wilson_vector *)(F_PT((s+FETCH_UP),src)) );
      prefetch_W( (wilson_vector *)(F_PT((s+FETCH_UP),dest)) );
    }
    for(b=0;b<2;b++){
      for(j=0;j<6;j++){
	
	/* diagonal part */
	CMULREAL(((wilson_block_vector *)F_PT(s,src))->b[b][j],
		 clov_diag[i].di[b][j],
		 ctmp );
	
	/* lower triangular part */
	jk=j*(j-1)/2;
	for(k=0;k<j;k++){
	  CMULSUM(clov[i].tr[b][jk],
		  ((wilson_block_vector *)F_PT(s,src))->b[b][k],
		  ctmp );
	  jk++;
	}

	/* upper triangular part */
	for(k=j+1;k<6;k++){
	  kj=k*(k-1)/2+j;
	  CMULJ_SUM(clov[i].tr[b][kj],
		    ((wilson_block_vector *)F_PT(s,src))->b[b][k],
		    ctmp );
	}
	((wilson_block_vector *)F_PT(s,dest))->b[b][j] = ctmp;
      }  /* j loop */
    } /* b loop */
  } END_LOOP /* site loop */
  
} /* mult_this_ldu_site */


/******************* copy_field ***********************************/
static void copy_field( wilson_vector *src, wilson_vector *dest,
			int parity )
{
  int i;
  site *s;
  FORSOMEPARITYDOMAIN(i,s,parity){
    dest[i] = src[i];
  } END_LOOP
}
/******************* mult_this_ldu_field **********************************/

void mult_this_ldu_field(
  clover *my_clov,
  wilson_vector *src,  /* RECAST as wilson_block_vector */
  wilson_vector *dest, /* RECAST as wilson_block_vector */
  int parity
  )
{
  triangular *clov = my_clov->clov;
  diagonal *clov_diag = my_clov->clov_diag;
  register int i;
  register site *s;
  register int b,j,k,jk,kj;
  register complex ctmp;
  wilson_block_vector *srcb = (wilson_block_vector *)src;
  wilson_block_vector *destb = (wilson_block_vector *)dest;
  
  if(my_clov->Clov_c == 0){
    copy_field(src, dest, parity);
    return;
  }

  FORSOMEPARITYDOMAIN(i,s,parity){
    if( i < loopend-FETCH_UP ){
      prefetch_W( &srcb[i+FETCH_UP] );
      prefetch_W( &destb[i+FETCH_UP] );
    }
    for(b=0;b<2;b++){
      for(j=0;j<6;j++){
	
	/* diagonal part */
	CMULREAL( srcb[i].b[b][j], clov_diag[i].di[b][j], ctmp );
	
	/* lower triangular part */
	jk=j*(j-1)/2;
	for(k=0;k<j;k++){
	  CMULSUM(clov[i].tr[b][jk],
		  srcb[i].b[b][k],
		  ctmp );
	  jk++;
	}

	/* upper triangular part */
	for(k=j+1;k<6;k++){
	  kj=k*(k-1)/2+j;
	  CMULJ_SUM(clov[i].tr[b][kj], srcb[i].b[b][k], ctmp );
	}
	destb[i].b[b][j] = ctmp;
      }  /* j loop */
    } /* b loop */
  } END_LOOP /* site loop */
  
} /* mult_this_ldu_field */

/******** tr_sigma_ldu_site *************/
/* MIMD version 7, May 95, MBW */
/* Used in clover_diagonal */
/* Left multiplies a color-dirac matrix stored in "ldu" form
by sigma_mu_nu = gamma_mu gamma_nu where

sigma(XUP,YUP)          sigma(XUP,ZUP)          sigma(XUP,TUP)
        -i  0  0  0              0 -1  0  0              0  i  0  0
         0  i  0  0              1  0  0  0              i  0  0  0
         0  0 -i  0              0  0  0 -1              0  0  0 -i
         0  0  0  i              0  0  1  0              0  0 -i  0

sigma(YUP,ZUP)          sigma(YUP,TUP)          sigma(ZUP,TUP)
         0 -i  0  0              0 -1  0  0              i  0  0  0
        -i  0  0  0              1  0  0  0              0 -i  0  0
         0  0  0 -i              0  0  0  1              0  0 -i  0
         0  0 -i  0              0  0 -1  0              0  0  0  i

and sigma(nu,mu) = -sigma(mu,nu), and sums over the dirac indices.

*/


void tr_sigma_this_ldu_mu_nu_site( clover *my_clov, field_offset mat, 
				   int mu, int nu )
/* triang, diag are input & contain the color-dirac matrix
   mat is output: the resulting su3_matrix  */
{
  triangular *clov;
  diagonal *clov_diag;
  Real Clov_c;
  register site *s;
  register int mm = 0, nn = 0;  /* dummy directions */
  register Real pm = 0;
  register int i,j,k,jk,jk2,kj;
  Real rtmp;
  complex ctmp,ctmp1;
  char myname[] = "tr_sigma_this_ldu_mu_nu_site";

  if(my_clov == NULL){
    printf("%s(%d): NULL clover structure\n",myname, this_node);
  }

  Clov_c = my_clov->Clov_c;

  if(Clov_c == 0){
    FORODDSITESDOMAIN(i,s) {
      clear_su3mat((su3_matrix *)F_PT(s,mat));
    }
    return;
  }
  
  if(!my_clov->valid_clov){
    printf("%s(%d): Can't use an invalid clover term\n", myname, this_node);
    terminate(1);
  }

  clov = my_clov->clov;
  clov_diag = my_clov->clov_diag;

  /* take care of the case mu > nu by flipping them and mult. by -1 */
  if( mu < nu ) {
    mm = mu; nn = nu; pm = 1.0;
  }
  else if( mu > nu) {
    mm = nu; nn = mu; pm = -1.0;
  }
  else
    printf("BAD CALL by %s: mu=%d, nu=%d\n",myname,mu,nu);
  
  switch(mm) {
  case(XUP):
    switch(nn) {
    case(YUP):
      FORODDSITESDOMAIN(i,s) {
	/* diagonal part */
	for(j=0;j<3;j++) {
	  rtmp = clov_diag[i].di[0][j+3]
	    - clov_diag[i].di[0][j];
	  rtmp -= clov_diag[i].di[1][j];
	  rtmp += clov_diag[i].di[1][j+3];
	  ((su3_matrix *)F_PT(s,mat))->e[j][j].real = 0.0;
	  ((su3_matrix *)F_PT(s,mat))->e[j][j].imag = rtmp;
	}
	/* triangular part */
	jk = 0;
	for(j=1;j<3;j++) {
	  jk2 = (j+3)*(j+2)/2 + 3;
	  for(k=0;k<j;k++) {
	    CSUB( clov[i].tr[0][jk2],
		  clov[i].tr[0][jk], ctmp );
	    CSUB( ctmp, clov[i].tr[1][jk],
		  ctmp );
	    CADD( ctmp, clov[i].tr[1][jk2],
		  ctmp );
	    TIMESPLUSI( ctmp, ((su3_matrix *)F_PT(s,mat))->e[j][k] );
	    CONJG( ctmp, ctmp);
	    TIMESPLUSI( ctmp, ((su3_matrix *)F_PT(s,mat))->e[k][j] );
	    jk++; jk2++;
	  }
	}
      }
      break;
    case(ZUP):
      FORODDSITESDOMAIN(i,s) {
	for(j=0;j<3;j++) {
	  jk = (j+3)*(j+2)/2;
	  for(k=0;k<3;k++) {
	    kj = (k+3)*(k+2)/2 + j;
	    CONJG( clov[i].tr[0][kj], ctmp);
	    CSUB( ctmp, clov[i].tr[0][jk],				ctmp );
	    CONJG( clov[i].tr[1][kj], ctmp1);
	    CADD( ctmp, ctmp1, ctmp);
	    CSUB( ctmp, clov[i].tr[1][jk],
		  ((su3_matrix *)F_PT(s,mat))->e[j][k] );
	    jk++;
	  }
	}
      }
      break;
    case(TUP):
      FORODDSITESDOMAIN(i,s) {
	for(j=0;j<3;j++) {
	  jk = (j+3)*(j+2)/2;
	  for(k=0;k<3;k++) {
	    kj = (k+3)*(k+2)/2 + j;
	    CONJG( clov[i].tr[0][kj], ctmp);
	    CADD( ctmp, clov[i].tr[0][jk],
		  ctmp );
	    CONJG( clov[i].tr[1][kj], ctmp1);
	    CSUB( ctmp, ctmp1, ctmp);
	    CSUB( ctmp, clov[i].tr[1][jk],
		  ctmp );
	    TIMESPLUSI( ctmp, ((su3_matrix *)F_PT(s,mat))->e[j][k] );
	    jk++;
	  }
	}
      }
      break;
    default:
      printf("BAD CALL in %s: mu=%d, nu=%d\n",myname,mu,nu);
    }
    break;
  case(YUP):
    switch(nn) {
    case(ZUP):
      FORODDSITESDOMAIN(i,s) {
	for(j=0;j<3;j++) {
	  jk = (j+3)*(j+2)/2;
	  for(k=0;k<3;k++) {
	    kj = (k+3)*(k+2)/2 + j;
	    CONJG( clov[i].tr[0][kj], ctmp);
	    CADD( ctmp, clov[i].tr[0][jk],
		  ctmp );
	    CONJG( clov[i].tr[1][kj], ctmp1);
	    CADD( ctmp, ctmp1, ctmp);
	    CADD( ctmp, clov[i].tr[1][jk],
		  ctmp );
	    TIMESMINUSI( ctmp, ((su3_matrix *)F_PT(s,mat))->e[j][k] );
	    jk++;
	  }
	}
      }
      break;
    case(TUP):
      FORODDSITESDOMAIN(i,s) {
	for(j=0;j<3;j++) {
	  jk = (j+3)*(j+2)/2;
	  for(k=0;k<3;k++) {
	    kj = (k+3)*(k+2)/2 + j;
	    CONJG( clov[i].tr[0][kj], ctmp);
	    CSUB( ctmp, clov[i].tr[0][jk],
		  ctmp );
	    CONJG( clov[i].tr[1][kj], ctmp1);
	    CSUB( ctmp, ctmp1, ctmp);
	    CADD( ctmp, clov[i].tr[1][jk],
		  ((su3_matrix *)F_PT(s,mat))->e[j][k] );
	    jk++;
	  }
	}
      }
      break;
    default:
      printf("BAD CALL in %s: mu=%d, nu=%d\n",myname,mu,nu);
    }
    break;
  case(ZUP):
    switch(nn) {
    case(TUP):
      FORODDSITESDOMAIN(i,s) {
	/* diagonal part */
	for(j=0;j<3;j++) {
	  rtmp = clov_diag[i].di[0][j]
	    - clov_diag[i].di[0][j+3];
	  rtmp -= clov_diag[i].di[1][j];
	  rtmp += clov_diag[i].di[1][j+3];
	  ((su3_matrix *)F_PT(s,mat))->e[j][j].real = 0.0;
	  ((su3_matrix *)F_PT(s,mat))->e[j][j].imag = rtmp;
	}
	/* triangular part */
	jk = 0;
	for(j=1;j<3;j++) {
	  jk2 = (j+3)*(j+2)/2 + 3;
	  for(k=0;k<j;k++) {
	    CSUB( clov[i].tr[0][jk],
		  clov[i].tr[0][jk2], ctmp );
	    CSUB( ctmp, clov[i].tr[1][jk],
		  ctmp );
	    CADD( ctmp, clov[i].tr[1][jk2],
		  ctmp );
	    TIMESPLUSI( ctmp, ((su3_matrix *)F_PT(s,mat))->e[j][k] );
	    CONJG( ctmp, ctmp);
	    TIMESPLUSI( ctmp, ((su3_matrix *)F_PT(s,mat))->e[k][j] );
	    jk++; jk2++;
	  }
	}
      }
      break;
    default:
      printf("BAD CALL in %s: mu=%d, nu=%d\n",myname,mu,nu);
    }
    break;
  default:
    printf("BAD CALL in %s: mu=%d, nu=%d\n",myname,mu,nu);
  }
  
  /* multiply by pm = +/- 1  */
  FORODDSITESDOMAIN(i,s) {
    for(j=0;j<3;j++) for(k=0;k<3;k++)
      CMULREAL( ((su3_matrix *)F_PT(s,mat))->e[j][k], pm,
		((su3_matrix *)F_PT(s,mat))->e[j][k] );
  }
  
} /* tr_sigma_this_ldu_mu_nu_site */


/******************* free_this_clov *********************************/
/* Frees the large allocated members but not the structure itself */

void free_this_clov(clover *my_clov)
{
  
  if(my_clov != NULL){
    if(my_clov->clov != NULL)free(my_clov->clov);
    if(my_clov->clov_diag != NULL)free(my_clov->clov_diag);
    if(my_clov->clov_raw != NULL)free(my_clov->clov_raw);
    if(my_clov->clov_diag_raw != NULL)free(my_clov->clov_diag_raw);
    my_clov->clov = NULL;
    my_clov->clov_diag = NULL;
    my_clov->clov_raw = NULL;
    my_clov->clov_diag_raw = NULL;
  }
} /* free_this_clov */


/******************* destroy_this_clov *********************************/
/* Frees the large allocated members AND the structure itself          */

void destroy_this_clov(clover **my_clov)
{
  
  if(*my_clov != NULL){
    free_this_clov(*my_clov);
    *my_clov = NULL;
  }
} /* destroy_this_clov */


/******************* invalidate_clov ************************************/
/* Mark clover term for refreshing                                      */

void invalidate_this_clov(clover *my_clov){
  if(my_clov != NULL)
    my_clov->valid_clov = 0;
}

// /******************* make_clov ***************************************/
// /* For making only the global clover term */
// 
// void make_clov(Real Clov_c)
// {
//   if(global_clov == NULL){
//     global_clov = create_clov();
//     if(global_clov == NULL)
//       terminate(1);
//   }
//   compute_clov(global_clov, Clov_c);
// }
// 
// /******************* make_clovinv ***********************************/
// /* For inverting only the global clover term */
// 
// double make_clovinv(int parity)
// {
//   compute_clovinv(global_clov, parity);
//   return global_clov->trlogA;
// }
// 
///******************* mult_ldu_site ***********************************/
///* For multiplying only by the global clover term */
//void mult_ldu_site(
//  clover *my_clov;
//  field_offset src,   /* type wilson_vector RECAST AS wilson_block_vector */
//  field_offset dest,  /* type wilson_vector RECAST AS wilson_block_vector */
//  int parity
//  )
//{
//  if(my_clov->Clov_c == 0)
//    copy_site(src, dest, parity);
//  else
//    mult_this_ldu_site(my_clov, src, dest, parity);
//}
//
///******************* mult_ldu_field ***********************************/
///* For multiplying only by the global clover term */
//void mult_ldu_field(
//  wilson_vector *src,
//  wilson_vector *dest,
//  int parity
//  )
//{
//  if(global_clov->Clov_c == 0)
//    copy_field(src, dest, parity);
//  else
//    mult_this_ldu_field(global_clov, src, dest, parity);
//}
//
///******************* free_clov ***************************************/
//
//void free_clov(void)
//{
//  free_this_clov(global_clov);
//  global_clov = NULL;
//
//} /* free_clov */
//
///******************* invalidate_clov ************************************/
///* Mark clover term for refreshing                                      */
//
//void invalidate_clov(void){
//  invalidate_this_clov(global_clov);
//}

