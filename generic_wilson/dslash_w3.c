/*************  dslash_w.c *******************************/
/* MIMD version 6 */

/* 9/06/03 added gathers from temp - CD */

/*
	dslash_w(F_OFFSET(psi),F_OFFSET(mp),isign,l_parity);
Compute SUM_dirs ( 
    ( 1 + isign*gamma[dir] ) * U(x,dir) * src(x+dir)
  + ( 1 - isign*gamma[dir] ) * U_adj(x-dir,dir) * src(x-dir)
)

*/

#include "generic_wilson_includes.h"
/* Temporary work space for dslash_w_on_temp and dslash_w_on_temp_special */ 
static half_wilson_vector *htmp[8] ;

/* Flag indicating if temp is allocated               */
static int temp_not_allocated=1 ;

/* Temporary work space for gauge links */ 
static su3_matrix *t_links;

/* Flag indicating if temp is allocated               */
static int tmp_links_not_set = 1;

void malloc_dslash_temps(){
  int j;

  if(!temp_not_allocated)return;
  for( j = 0; j < 8; j++ ){
    htmp[j] =(half_wilson_vector *)malloc(sites_on_node*sizeof(half_wilson_vector));
    if(htmp[j] == NULL){
      printf("node %d can't malloc htmp[%d]\n",this_node,j);
      terminate(1);
    }
  }
  temp_not_allocated = 0 ;
}

void cleanup_dslash_temps(){
  register int i ;
  if(!temp_not_allocated)
    for(i=0;i<8;i++) {
      free(htmp[i]) ; 
    }
  temp_not_allocated=1 ;
}

void setup_tmp_links(){
  register int i, dir;
  register site *s;

  if(!tmp_links_not_set)return;
  t_links =(su3_matrix *)malloc(sites_on_node*4*sizeof(su3_matrix));
  if(t_links == NULL){
      printf("node %d can't malloc t_links\n",this_node);
      terminate(1);
  }
  FORALLSITES(i,s){
    FORALLUPDIR(dir){
      t_links[4*i+dir] = lattice[i].link[dir];
    }
  }
  tmp_links_not_set = 0 ;
}

void cleanup_tmp_links(){
  register int i ;
  if(!tmp_links_not_set)
      free(t_links) ; 
  tmp_links_not_set=1 ;
}

/* Included here to allow for inlining */
/**************  m_mat_hwvec.c  (in su3.a) ***********************
*									*
* void mult_su3_mat_hwvec(su3_matrix *mat,				*
*	half_wilson_vector *src,*dest)					*
*  multiply a Wilson half-vector by a matrix				*
*  dest  <-  mat*src							*
*/

void mult_su3_mat_hwvec2( su3_matrix *mat, half_wilson_vector *src,
	half_wilson_vector *dest ){

#ifdef NATIVEDOUBLE
  register double a0r,a0i,a1r,a1i,a2r,a2i;
  register double b0r,b0i,b1r,b1i,b2r,b2i;
#else
  register Real a0r,a0i,a1r,a1i,a2r,a2i;
  register Real b0r,b0i,b1r,b1i,b2r,b2i;
#endif
  
/*    mult_su3_mat_vec(mat, &(src->h[0]), &(dest->h[0]) ); */

  a0r=mat->e[0][0].real;    a0i=mat->e[0][0].imag;
  b0r=src->h[0].c[0].real;  b0i=src->h[0].c[0].imag;
  a1r=mat->e[0][1].real;    a1i=mat->e[0][1].imag;
  b1r=src->h[0].c[1].real;  b1i=src->h[0].c[1].imag;
  a2r=mat->e[0][2].real;    a2i=mat->e[0][2].imag;
  b2r=src->h[0].c[2].real;  b2i=src->h[0].c[2].imag;

  dest->h[0].c[0].real = a0r*b0r - a0i*b0i + a1r*b1r - a1i*b1i + a2r*b2r - a2i*b2i;
  dest->h[0].c[0].imag = a0r*b0i + a0i*b0r + a1r*b1i + a1i*b1r + a2r*b2i + a2i*b2r;
  
  a0r=mat->e[1][0].real;    a0i=mat->e[1][0].imag;
  b0r=src->h[0].c[0].real;  b0i=src->h[0].c[0].imag;
  a1r=mat->e[1][1].real;    a1i=mat->e[1][1].imag;
  b1r=src->h[0].c[1].real;  b1i=src->h[0].c[1].imag;
  a2r=mat->e[1][2].real;    a2i=mat->e[1][2].imag;
  b2r=src->h[0].c[2].real;  b2i=src->h[0].c[2].imag;

  dest->h[0].c[1].real = a0r*b0r - a0i*b0i + a1r*b1r - a1i*b1i + a2r*b2r - a2i*b2i;
  dest->h[0].c[1].imag = a0r*b0i + a0i*b0r + a1r*b1i + a1i*b1r + a2r*b2i + a2i*b2r;

  a0r=mat->e[2][0].real;    a0i=mat->e[2][0].imag;
  b0r=src->h[0].c[0].real;  b0i=src->h[0].c[0].imag;
  a1r=mat->e[2][1].real;    a1i=mat->e[2][1].imag;
  b1r=src->h[0].c[1].real;  b1i=src->h[0].c[1].imag;
  a2r=mat->e[2][2].real;    a2i=mat->e[2][2].imag;
  b2r=src->h[0].c[2].real;  b2i=src->h[0].c[2].imag;

  dest->h[0].c[2].real = a0r*b0r - a0i*b0i + a1r*b1r - a1i*b1i + a2r*b2r - a2i*b2i;
  dest->h[0].c[2].imag = a0r*b0i + a0i*b0r + a1r*b1i + a1i*b1r + a2r*b2i + a2i*b2r;

/*    mult_su3_mat_vec(mat, &(src->h[1]), &(dest->h[1]) ); */

  a0r=mat->e[0][0].real;    a0i=mat->e[0][0].imag;
  b0r=src->h[1].c[0].real;  b0i=src->h[1].c[0].imag;
  a1r=mat->e[0][1].real;    a1i=mat->e[0][1].imag;
  b1r=src->h[1].c[1].real;  b1i=src->h[1].c[1].imag;
  a2r=mat->e[0][2].real;    a2i=mat->e[0][2].imag;
  b2r=src->h[1].c[2].real;  b2i=src->h[1].c[2].imag;

  dest->h[1].c[0].real = a0r*b0r - a0i*b0i + a1r*b1r - a1i*b1i + a2r*b2r - a2i*b2i;
  dest->h[1].c[0].imag = a0r*b0i + a0i*b0r + a1r*b1i + a1i*b1r + a2r*b2i + a2i*b2r;
  
  a0r=mat->e[1][0].real;    a0i=mat->e[1][0].imag;
  b0r=src->h[1].c[0].real;  b0i=src->h[1].c[0].imag;
  a1r=mat->e[1][1].real;    a1i=mat->e[1][1].imag;
  b1r=src->h[1].c[1].real;  b1i=src->h[1].c[1].imag;
  a2r=mat->e[1][2].real;    a2i=mat->e[1][2].imag;
  b2r=src->h[1].c[2].real;  b2i=src->h[1].c[2].imag;

  dest->h[1].c[1].real = a0r*b0r - a0i*b0i + a1r*b1r - a1i*b1i + a2r*b2r - a2i*b2i;
  dest->h[1].c[1].imag = a0r*b0i + a0i*b0r + a1r*b1i + a1i*b1r + a2r*b2i + a2i*b2r;

  a0r=mat->e[2][0].real;    a0i=mat->e[2][0].imag;
  b0r=src->h[1].c[0].real;  b0i=src->h[1].c[0].imag;
  a1r=mat->e[2][1].real;    a1i=mat->e[2][1].imag;
  b1r=src->h[1].c[1].real;  b1i=src->h[1].c[1].imag;
  a2r=mat->e[2][2].real;    a2i=mat->e[2][2].imag;
  b2r=src->h[1].c[2].real;  b2i=src->h[1].c[2].imag;

  dest->h[1].c[2].real = a0r*b0r - a0i*b0i + a1r*b1r - a1i*b1i + a2r*b2r - a2i*b2i;
  dest->h[1].c[2].imag = a0r*b0i + a0i*b0r + a1r*b1i + a1i*b1r + a2r*b2i + a2i*b2r;

}


/* Included here to allow for inlining */
/**************  m_amat_hwvec.c  (in su3.a) **********************
*									*
*  void mult_adj_su3_mat_hwvec2( su3_matrix *mat,			*
*	half_wilson_vector *src,*dest )					*
*  multiply a Wilson half-vector by the adjoint of a matrix		*
*/


void mult_adj_su3_mat_hwvec2( su3_matrix *mat,
       half_wilson_vector *src, half_wilson_vector *dest ){

#ifdef NATIVEDOUBLE
  register double a0r,a0i,a1r,a1i,a2r,a2i;
  register double b0r,b0i,b1r,b1i,b2r,b2i;
#else
  register Real a0r,a0i,a1r,a1i,a2r,a2i;
  register Real b0r,b0i,b1r,b1i,b2r,b2i;
#endif

/*    mult_adj_su3_mat_vec(mat, &(src->h[0]), &(dest->h[0]) ); */
  
  a0r=mat->e[0][0].real;   a0i=mat->e[0][0].imag;
  b0r=src->h[0].c[0].real; b0i=src->h[0].c[0].imag;
  a1r=mat->e[1][0].real;   a1i=mat->e[1][0].imag;
  b1r=src->h[0].c[1].real; b1i=src->h[0].c[1].imag;
  a2r=mat->e[2][0].real;   a2i=mat->e[2][0].imag;
  b2r=src->h[0].c[2].real; b2i=src->h[0].c[2].imag;
  
  dest->h[0].c[0].real = a0r*b0r + a0i*b0i + a1r*b1r + a1i*b1i + a2r*b2r + a2i*b2i;
  dest->h[0].c[0].imag = a0r*b0i - a0i*b0r + a1r*b1i - a1i*b1r + a2r*b2i - a2i*b2r;
  
  a0r=mat->e[0][1].real;   a0i=mat->e[0][1].imag;
  b0r=src->h[0].c[0].real; b0i=src->h[0].c[0].imag;  
  a1r=mat->e[1][1].real;   a1i=mat->e[1][1].imag;
  b1r=src->h[0].c[1].real; b1i=src->h[0].c[1].imag;
  a2r=mat->e[2][1].real;   a2i=mat->e[2][1].imag;
  b2r=src->h[0].c[2].real; b2i=src->h[0].c[2].imag;
  
  dest->h[0].c[1].real = a0r*b0r + a0i*b0i + a1r*b1r + a1i*b1i + a2r*b2r + a2i*b2i;
  dest->h[0].c[1].imag = a0r*b0i - a0i*b0r + a1r*b1i - a1i*b1r + a2r*b2i - a2i*b2r;
  
  a0r=mat->e[0][2].real;   a0i=mat->e[0][2].imag;
  b0r=src->h[0].c[0].real; b0i=src->h[0].c[0].imag;
  a1r=mat->e[1][2].real;   a1i=mat->e[1][2].imag;
  b1r=src->h[0].c[1].real; b1i=src->h[0].c[1].imag;
  a2r=mat->e[2][2].real;   a2i=mat->e[2][2].imag;
  b2r=src->h[0].c[2].real; b2i=src->h[0].c[2].imag;
  
  dest->h[0].c[2].real = a0r*b0r + a0i*b0i + a1r*b1r + a1i*b1i + a2r*b2r + a2i*b2i;
  dest->h[0].c[2].imag = a0r*b0i - a0i*b0r + a1r*b1i - a1i*b1r + a2r*b2i - a2i*b2r;


/*    mult_adj_su3_mat_vec(mat, &(src->h[1]), &(dest->h[1]) ); */

  a0r=mat->e[0][0].real;   a0i=mat->e[0][0].imag;
  b0r=src->h[1].c[0].real; b0i=src->h[1].c[0].imag;
  a1r=mat->e[1][0].real;   a1i=mat->e[1][0].imag;
  b1r=src->h[1].c[1].real; b1i=src->h[1].c[1].imag;
  a2r=mat->e[2][0].real;   a2i=mat->e[2][0].imag;
  b2r=src->h[1].c[2].real; b2i=src->h[1].c[2].imag;
  
  dest->h[1].c[0].real = a0r*b0r + a0i*b0i + a1r*b1r + a1i*b1i + a2r*b2r + a2i*b2i;
  dest->h[1].c[0].imag = a0r*b0i - a0i*b0r + a1r*b1i - a1i*b1r + a2r*b2i - a2i*b2r;
  
  a0r=mat->e[0][1].real;   a0i=mat->e[0][1].imag;
  b0r=src->h[1].c[0].real; b0i=src->h[1].c[0].imag;  
  a1r=mat->e[1][1].real;   a1i=mat->e[1][1].imag;
  b1r=src->h[1].c[1].real; b1i=src->h[1].c[1].imag;
  a2r=mat->e[2][1].real;   a2i=mat->e[2][1].imag;
  b2r=src->h[1].c[2].real; b2i=src->h[1].c[2].imag;
  
  dest->h[1].c[1].real = a0r*b0r + a0i*b0i + a1r*b1r + a1i*b1i + a2r*b2r + a2i*b2i;
  dest->h[1].c[1].imag = a0r*b0i - a0i*b0r + a1r*b1i - a1i*b1r + a2r*b2i - a2i*b2r;
  
  a0r=mat->e[0][2].real;   a0i=mat->e[0][2].imag;
  b0r=src->h[1].c[0].real; b0i=src->h[1].c[0].imag;
  a1r=mat->e[1][2].real;   a1i=mat->e[1][2].imag;
  b1r=src->h[1].c[1].real; b1i=src->h[1].c[1].imag;
  a2r=mat->e[2][2].real;   a2i=mat->e[2][2].imag;
  b2r=src->h[1].c[2].real; b2i=src->h[1].c[2].imag;
  
  dest->h[1].c[2].real = a0r*b0r + a0i*b0i + a1r*b1r + a1i*b1i + a2r*b2r + a2i*b2i;
  dest->h[1].c[2].imag = a0r*b0i - a0i*b0r + a1r*b1i - a1i*b1r + a2r*b2i - a2i*b2r;

}

/* Included here to allow for inlining */
/*****************  grow4wvecs.c  (in su3.a) ****************************
*									*
*  If sum=0,								*
*  Grow and add four wilson_vectors 					*
*  If sum=1,								*
*  Grow and sum four wilson_vectors to another wilson_vector		*
* void grow_add_four_wvecs2(a,b1,b2,b3,b4,sign,sum)				*
* wilson_vector *a; half_wilson_vector *b1,*b2,*b3,*b4;			*
* int sign,sum;								*
* A  <-  B1 + B2 + B3 + B4   or						*
* A  <-  A + B1 + B2 + B3 + B4						*
* B1 is expanded using gamma_x, B2 using gamma_y, etc. 			*
*/

/* grow and sum four wilson_vectors */

/* For the RS6000 */

void grow_add_four_wvecs2( wilson_vector *a, half_wilson_vector *b1,
        half_wilson_vector *b2, half_wilson_vector *b3,
	half_wilson_vector *b4, int sign, int sum ){
  int i;
  if(sum==0)
    {
      /* wp_grow( b1,a,XUP,sign); */
      
      /* case XUP: */
      if(sign==PLUS)
	{
	  for(i=0;i<3;i++){
	    a->d[0].c[i]      = b1->h[0].c[i];
	    a->d[1].c[i]      = b1->h[1].c[i];
	    TIMESMINUSI( b1->h[0].c[i], a->d[3].c[i]);
	    TIMESMINUSI( b1->h[1].c[i], a->d[2].c[i]);
	  }
	}
      else
	{
	  /* case XDOWN: */
	  for(i=0;i<3;i++){
	    a->d[0].c[i]      = b1->h[0].c[i];
	    a->d[1].c[i]      = b1->h[1].c[i];
	    TIMESPLUSI( b1->h[0].c[i], a->d[3].c[i]);
	    TIMESPLUSI( b1->h[1].c[i], a->d[2].c[i]);
	  }
	}
    }
  else
    {
      /*wp_grow_add( b1,a,XUP,sign); */
      
      /* case XUP: */
      if(sign==PLUS)
	{
	  for(i=0;i<3;i++){
	    CSUM( a->d[0].c[i], b1->h[0].c[i]);
	    CSUM( a->d[1].c[i], b1->h[1].c[i]);
	    CSUM_TMI( a->d[2].c[i], b1->h[1].c[i] );
	    CSUM_TMI( a->d[3].c[i], b1->h[0].c[i] );
	  }
	}
      else
	{
	  /* case XDOWN: */
	  for(i=0;i<3;i++){
	    CSUM( a->d[0].c[i], b1->h[0].c[i]);
	    CSUM( a->d[1].c[i], b1->h[1].c[i]);
	    CSUM_TPI( a->d[2].c[i], b1->h[1].c[i] );
	    CSUM_TPI( a->d[3].c[i], b1->h[0].c[i] );
	  }
	}
    }
  
  /* wp_grow_add( b2,a,YUP,sign); */
  
  if(sign==PLUS)
    {
      /* case YUP: */
      for(i=0;i<3;i++){
	CSUM( a->d[0].c[i], b2->h[0].c[i]);
	CSUM( a->d[1].c[i], b2->h[1].c[i]);
	CSUM( a->d[2].c[i], b2->h[1].c[i]);
	CSUB( a->d[3].c[i], b2->h[0].c[i], a->d[3].c[i] );
      }
    }
  else
    {
      /* case YDOWN: */
      for(i=0;i<3;i++){
	CSUM( a->d[0].c[i], b2->h[0].c[i]);
	CSUM( a->d[1].c[i], b2->h[1].c[i]);
	CSUB( a->d[2].c[i], b2->h[1].c[i], a->d[2].c[i] );
	CSUM( a->d[3].c[i], b2->h[0].c[i]);
      }
    }
  
  /* wp_grow_add( b3,a,ZUP,sign); */
  
  if(sign==PLUS)
    {
      /* case ZUP: */
      for(i=0;i<3;i++){
	CSUM( a->d[0].c[i], b3->h[0].c[i]);
	CSUM( a->d[1].c[i], b3->h[1].c[i]);
	CSUM_TMI( a->d[2].c[i], b3->h[0].c[i] );
	CSUM_TPI( a->d[3].c[i], b3->h[1].c[i] );
      }
    }
  else
    {
      /* case ZDOWN:*/
      for(i=0;i<3;i++){
	CSUM( a->d[0].c[i], b3->h[0].c[i]);
	CSUM( a->d[1].c[i], b3->h[1].c[i]);
	CSUM_TPI( a->d[2].c[i], b3->h[0].c[i] );
	CSUM_TMI( a->d[3].c[i], b3->h[1].c[i] );
      }
    }
  
  /* wp_grow_add( b4,a,TUP,sign); */
  
  if(sign==PLUS)
    {
      /* case TUP: */
      for(i=0;i<3;i++){
	CSUM( a->d[0].c[i], b4->h[0].c[i]);
	CSUM( a->d[1].c[i], b4->h[1].c[i]);
	CSUM( a->d[2].c[i], b4->h[0].c[i]);
	CSUM( a->d[3].c[i], b4->h[1].c[i]);
      }
    }
  else
    {
      /* case TDOWN: */
      for(i=0;i<3;i++){
	CSUM( a->d[0].c[i], b4->h[0].c[i]);
	CSUM( a->d[1].c[i], b4->h[1].c[i]);
	CSUB( a->d[2].c[i], b4->h[0].c[i], a->d[2].c[i] );
	CSUB( a->d[3].c[i], b4->h[1].c[i], a->d[3].c[i] );
      }
    }
}

/* Included here to allow for inlining */
/*****************  wp_shrink4.c  (in su3.a) ****************************
*									*
* Shrink a wilson vector in four directions, producing four		*
*  half_wilson_vectors.							*
* void wp_shrink_4dir(  wilson_vector *a,  half_wilson_vector *b1,	*
*       half_wilson_vector *b2, half_wilson_vector *b3,			*
*       half_wilson_vector *b4, int sign );				*
* B1 <- (1 +- gamma_x)A,, projection					*
*  argument "sign" is sign of gamma matrix.				*
*  See wp_shrink.c for definitions of gamma matrices and eigenvectors.	*
*/

void wp_shrink_4dir2( wilson_vector *a,  half_wilson_vector *b1,
        half_wilson_vector *b2, half_wilson_vector *b3,
        half_wilson_vector *b4, int sign ){
  register int i; /*color*/
  
/*    wp_shrink( a,b1,XUP,sign); */

  if(sign==PLUS)
    {
    /* case XUP: */
      for(i=0;i<3;i++){
	b1->h[0].c[i].real = a->d[0].c[i].real - a->d[3].c[i].imag;
	b1->h[0].c[i].imag = a->d[0].c[i].imag + a->d[3].c[i].real;
	b1->h[1].c[i].real = a->d[1].c[i].real - a->d[2].c[i].imag;
	b1->h[1].c[i].imag = a->d[1].c[i].imag + a->d[2].c[i].real;
      }
    }
  else
    {
      /* case XDOWN: */
      for(i=0;i<3;i++){
	b1->h[0].c[i].real = a->d[0].c[i].real + a->d[3].c[i].imag;
	b1->h[0].c[i].imag = a->d[0].c[i].imag - a->d[3].c[i].real;
	b1->h[1].c[i].real = a->d[1].c[i].real + a->d[2].c[i].imag;
	b1->h[1].c[i].imag = a->d[1].c[i].imag - a->d[2].c[i].real;
      }
    }
  
  
  /*    wp_shrink( a,b2,YUP,sign); */
  
  if(sign==PLUS)
    {
      /* case YUP: */
      for(i=0;i<3;i++){
	b2->h[0].c[i].real = a->d[0].c[i].real - a->d[3].c[i].real;
	b2->h[0].c[i].imag = a->d[0].c[i].imag - a->d[3].c[i].imag;
	b2->h[1].c[i].real = a->d[1].c[i].real + a->d[2].c[i].real;
	b2->h[1].c[i].imag = a->d[1].c[i].imag + a->d[2].c[i].imag;
      }
      
    }
  else
    {
      /* case YDOWN: */
      for(i=0;i<3;i++){
	b2->h[0].c[i].real = a->d[0].c[i].real + a->d[3].c[i].real;
	b2->h[0].c[i].imag = a->d[0].c[i].imag + a->d[3].c[i].imag;
	b2->h[1].c[i].real = a->d[1].c[i].real - a->d[2].c[i].real;
	b2->h[1].c[i].imag = a->d[1].c[i].imag - a->d[2].c[i].imag;
      }
    }
  
  /*    wp_shrink( a,b3,ZUP,sign); */

  if(sign==PLUS)
    {
      /* case ZUP: */
      for(i=0;i<3;i++){
	b3->h[0].c[i].real = a->d[0].c[i].real - a->d[2].c[i].imag;
	b3->h[0].c[i].imag = a->d[0].c[i].imag + a->d[2].c[i].real;
	b3->h[1].c[i].real = a->d[1].c[i].real + a->d[3].c[i].imag;
	b3->h[1].c[i].imag = a->d[1].c[i].imag - a->d[3].c[i].real;
      }
    }
  else
    {
      /* case ZDOWN: */
      for(i=0;i<3;i++){
	b3->h[0].c[i].real = a->d[0].c[i].real + a->d[2].c[i].imag;
	b3->h[0].c[i].imag = a->d[0].c[i].imag - a->d[2].c[i].real;
	b3->h[1].c[i].real = a->d[1].c[i].real - a->d[3].c[i].imag;
	b3->h[1].c[i].imag = a->d[1].c[i].imag + a->d[3].c[i].real;
      }
      
    }

/*    wp_shrink( a,b4,TUP,sign); */

  if(sign==PLUS)
    {
      /* case TUP: */
      for(i=0;i<3;i++){
	b4->h[0].c[i].real = a->d[0].c[i].real + a->d[2].c[i].real;
	b4->h[0].c[i].imag = a->d[0].c[i].imag + a->d[2].c[i].imag;
	b4->h[1].c[i].real = a->d[1].c[i].real + a->d[3].c[i].real;
	b4->h[1].c[i].imag = a->d[1].c[i].imag + a->d[3].c[i].imag;
      }
    }
  else
    {
      /* case TDOWN: */
      for(i=0;i<3;i++){
	b4->h[0].c[i].real = a->d[0].c[i].real - a->d[2].c[i].real;
	b4->h[0].c[i].imag = a->d[0].c[i].imag - a->d[2].c[i].imag;
	b4->h[1].c[i].real = a->d[1].c[i].real - a->d[3].c[i].real;
	b4->h[1].c[i].imag = a->d[1].c[i].imag - a->d[3].c[i].imag;
      }
    }
}


void dslash_w( field_offset src, field_offset dest, int isign, int parity) {
half_wilson_vector hwvx,hwvy,hwvz,hwvt;

register int i;
register site *s;
register int dir,otherparity;
msg_tag *tag[8];

    switch(parity) {
	case EVEN:      otherparity=ODD; break;
	case ODD:       otherparity=EVEN; break;
	case EVENANDODD:        otherparity=EVENANDODD; break;
    }

#ifdef MAXHTMP
    /* NOTE: We should be defining MAXHTMP in all applications using
       dslash and dslash_w */
    if(MAXHTMP < 8){
      printf("dslash: MAXHTMP must be 8 or more!\n");
      terminate(1);
    }
#endif
    if(N_POINTERS < 8){
      printf("dslash: N_POINTERS must be 8 or more!\n");
      terminate(1);
     }


    /* Take Wilson projection for src displaced in up direction, gather
       it to "our site" */
    FORSOMEPARITY(i,s,otherparity){
        wp_shrink_4dir2( (wilson_vector *)F_PT(s,src), &(s->htmp[XUP]),
	    &(s->htmp[YUP]), &(s->htmp[ZUP]), &(s->htmp[TUP]), isign);
    }
    for( dir=XUP; dir <= TUP; dir++) {
	tag[dir]=start_gather( F_OFFSET(htmp[dir]), sizeof(half_wilson_vector),
	    dir, parity, gen_pt[dir] );
    }

        /* Take Wilson projection for src displaced in down direction,
        multiply it by adjoint link matrix, gather it "up" */
    FORSOMEPARITY(i,s,otherparity){
        wp_shrink_4dir2( (wilson_vector *)F_PT(s,src),
	    &hwvx, &hwvy, &hwvz, &hwvt, -isign);
	mult_adj_su3_mat_hwvec2( &(s->link[XUP]), &hwvx, &(s->htmp[XDOWN]));
	mult_adj_su3_mat_hwvec2( &(s->link[YUP]), &hwvy, &(s->htmp[YDOWN]));
	mult_adj_su3_mat_hwvec2( &(s->link[ZUP]), &hwvz, &(s->htmp[ZDOWN]));
	mult_adj_su3_mat_hwvec2( &(s->link[TUP]), &hwvt, &(s->htmp[TDOWN]));
    }

    for( dir=XUP; dir <= TUP; dir++) {
	tag[OPP_DIR(dir)]=start_gather(F_OFFSET(htmp[OPP_DIR(dir)]), 
		sizeof(half_wilson_vector), OPP_DIR(dir),
		parity, gen_pt[OPP_DIR(dir)] );
    }


	/* Set dest to zero */
        /* Take Wilson projection for src displaced in up direction, gathered,
		multiply it by link matrix, expand it, and add.
		to dest */
    for( dir=XUP; dir <= TUP; dir++) {
	wait_gather(tag[dir]);
    }
    FORSOMEPARITY(i,s,parity){
	mult_su3_mat_hwvec2( &(s->link[XUP]), 
		(half_wilson_vector * )(gen_pt[XUP][i]), &hwvx ); 
	mult_su3_mat_hwvec2( &(s->link[YUP]), 
		(half_wilson_vector * )(gen_pt[YUP][i]), &hwvy ); 
	mult_su3_mat_hwvec2( &(s->link[ZUP]), 
		(half_wilson_vector * )(gen_pt[ZUP][i]), &hwvz ); 
	mult_su3_mat_hwvec2( &(s->link[TUP]), 
		(half_wilson_vector * )(gen_pt[TUP][i]), &hwvt ); 
	grow_add_four_wvecs2( (wilson_vector *)F_PT(s,dest),
	    &hwvx, &hwvy, &hwvz, &hwvt, isign, 0 ); /* "0" is NOSUM */
    }
    for( dir=XUP; dir <= TUP; dir++) {
	cleanup_gather(tag[dir]);
    }

        /* Take Wilson projection for src displaced in down direction,
        expand it, and add to dest */
    for( dir=XUP; dir <= TUP; dir++) {
	wait_gather(tag[OPP_DIR(dir)]);
    }

    FORSOMEPARITY(i,s,parity){
	grow_add_four_wvecs2( (wilson_vector *)F_PT(s,dest),
	    (half_wilson_vector *)(gen_pt[XDOWN][i]),
	    (half_wilson_vector *)(gen_pt[YDOWN][i]),
	    (half_wilson_vector *)(gen_pt[ZDOWN][i]),
	    (half_wilson_vector *)(gen_pt[TDOWN][i]),
	    -isign, 1 );	/* "1" SUMs in current dest */
    }
    for( dir=XUP; dir <= TUP; dir++) {
	cleanup_gather(tag[OPP_DIR(dir)]);
    }

} /* end (of dslash_w() ) */


/* Special dslash for use by congrad.  Uses restart_gather() when
  possible. Last argument is an integer, which will tell if
  gathers have been started.  If is_started=0,use
  start_gather, otherwise use restart_gather.
  Argument "tag" is a vector of a msg_tag *'s to use for
  the gathers.
  The calling program must clean up the gathers! */
void dslash_w_special(field_offset src,field_offset dest,
    int isign,int parity,msg_tag **tag,int is_started)

{
half_wilson_vector hwvx,hwvy,hwvz,hwvt;

register int i;
register site *s;
register int dir,otherparity;

    switch(parity) {
	case EVEN:      otherparity=ODD; break;
	case ODD:       otherparity=EVEN; break;
	case EVENANDODD:        otherparity=EVENANDODD; break;
    }

#ifdef MAXHTMP
    /* NOTE: We should be defining MAXHTMP in all applications using
       dslash and dslash_w */
    if(MAXHTMP < 8){
      printf("dslash_w_special: MAXHTMP must be 8 or more!\n");
      terminate(1);
    }
#endif
    if(N_POINTERS < 8){
      printf("dslash_w_special: N_POINTERS must be 8 or more!\n");
      terminate(1);
     }


    /* Take Wilson projection for src displaced in up direction, gather
       it to "our site" */
    FORSOMEPARITY(i,s,otherparity){
        wp_shrink_4dir2( (wilson_vector *)F_PT(s,src), &(s->htmp[XUP]),
	    &(s->htmp[YUP]), &(s->htmp[ZUP]), &(s->htmp[TUP]), isign);
    }

    for( dir=XUP; dir <= TUP; dir++) {
	if(is_started==0)tag[dir]=start_gather( F_OFFSET(htmp[dir]),
	    sizeof(half_wilson_vector), dir, parity, gen_pt[dir] );
	else restart_gather( F_OFFSET(htmp[dir]),
	    sizeof(half_wilson_vector), dir, parity, gen_pt[dir], 
			     tag[dir] );
    }

        /* Take Wilson projection for src displaced in down direction,
        multiply it by adjoint link matrix, gather it "up" */
    FORSOMEPARITY(i,s,otherparity){
        wp_shrink_4dir2( (wilson_vector *)F_PT(s,src),
	    &hwvx, &hwvy, &hwvz, &hwvt, -isign);
	mult_adj_su3_mat_hwvec2( &(s->link[XUP]), &hwvx, &(s->htmp[XDOWN]));
	mult_adj_su3_mat_hwvec2( &(s->link[YUP]), &hwvy, &(s->htmp[YDOWN]));
	mult_adj_su3_mat_hwvec2( &(s->link[ZUP]), &hwvz, &(s->htmp[ZDOWN]));
	mult_adj_su3_mat_hwvec2( &(s->link[TUP]), &hwvt, &(s->htmp[TDOWN]));
    }

    for( dir=XUP; dir <= TUP; dir++) {
	if(is_started==0)tag[OPP_DIR(dir)]=start_gather(
	    F_OFFSET(htmp[OPP_DIR(dir)]), sizeof(half_wilson_vector),
	    OPP_DIR(dir), parity, gen_pt[OPP_DIR(dir)] );
	else restart_gather(
	    F_OFFSET(htmp[OPP_DIR(dir)]), sizeof(half_wilson_vector),
	    OPP_DIR(dir), parity, gen_pt[OPP_DIR(dir)], 
	    tag[OPP_DIR(dir)] );
    }


	/* Set dest to zero */
        /* Take Wilson projection for src displaced in up direction, gathered,
		multiply it by link matrix, expand it, and add.
		to dest */
    for( dir=XUP; dir <= TUP; dir++) {
	wait_gather(tag[dir]);
    }
    FORSOMEPARITY(i,s,parity){
	mult_su3_mat_hwvec2( &(s->link[XUP]), 
		(half_wilson_vector * )(gen_pt[XUP][i]), &hwvx ); 
	mult_su3_mat_hwvec2( &(s->link[YUP]), 
		(half_wilson_vector * )(gen_pt[YUP][i]), &hwvy ); 
	mult_su3_mat_hwvec2( &(s->link[ZUP]), 
		(half_wilson_vector * )(gen_pt[ZUP][i]), &hwvz ); 
	mult_su3_mat_hwvec2( &(s->link[TUP]), 
		(half_wilson_vector * )(gen_pt[TUP][i]), &hwvt ); 
	grow_add_four_wvecs2( (wilson_vector *)F_PT(s,dest),
	    &hwvx, &hwvy, &hwvz, &hwvt, isign, 0 ); /* "0" is NOSUM */
    }

        /* Take Wilson projection for src displaced in down direction,
        expand it, and add to dest */
    for( dir=XUP; dir <= TUP; dir++) {
	wait_gather(tag[OPP_DIR(dir)]);
    }

    FORSOMEPARITY(i,s,parity){
	grow_add_four_wvecs2( (wilson_vector *)F_PT(s,dest),
	    (half_wilson_vector *)(gen_pt[XDOWN][i]),
	    (half_wilson_vector *)(gen_pt[YDOWN][i]),
	    (half_wilson_vector *)(gen_pt[ZDOWN][i]),
	    (half_wilson_vector *)(gen_pt[TDOWN][i]),
	    -isign, 1 );	/* "1" SUMs in current dest */
    }

} /* end (of dslash_w_special() ) */

void dslash_w_on_temp( wilson_vector *src, wilson_vector *dest, int isign, int parity) {
half_wilson_vector hwvx,hwvy,hwvz,hwvt;

register int i;
register site *s;
register int dir,otherparity;
msg_tag *tag[8];
su3_matrix *linkx,*linky,*linkz,*linkt;

  /* The calling program must clean up the temps! */
  malloc_dslash_temps();
  setup_tmp_links();

    switch(parity) {
	case EVEN:      otherparity=ODD; break;
	case ODD:       otherparity=EVEN; break;
	case EVENANDODD:        otherparity=EVENANDODD; break;
    }

    if(N_POINTERS < 8){
      printf("dslash: N_POINTERS must be 8 or more!\n");
      terminate(1);
     }


    /* Take Wilson projection for src displaced in up direction, gather
       it to "our site" */
    FORSOMEPARITY(i,s,otherparity){
        wp_shrink_4dir2( &(src[i]), &(htmp[XUP][i]),
	    &(htmp[YUP][i]), &(htmp[ZUP][i]), &(htmp[TUP][i]), isign);
    }
    for( dir=XUP; dir <= TUP; dir++) {
	tag[dir]=start_gather_from_temp( htmp[dir], sizeof(half_wilson_vector),
	    dir, parity, gen_pt[dir] );
    }

        /* Take Wilson projection for src displaced in down direction,
        multiply it by adjoint link matrix, gather it "up" */
    FORSOMEPARITY(i,s,otherparity){
        wp_shrink_4dir2( &(src[i]), &hwvx, &hwvy, &hwvz, &hwvt, -isign);

	linkx = &t_links[4*i+XUP];
	linky = &t_links[4*i+YUP];
	linkz = &t_links[4*i+ZUP];
	linkt = &t_links[4*i+TUP];
	mult_adj_su3_mat_hwvec2( linkx, &hwvx, &(htmp[XDOWN][i]));
	mult_adj_su3_mat_hwvec2( linky, &hwvy, &(htmp[YDOWN][i]));
	mult_adj_su3_mat_hwvec2( linkz, &hwvz, &(htmp[ZDOWN][i]));
	mult_adj_su3_mat_hwvec2( linkt, &hwvt, &(htmp[TDOWN][i]));
    }

    for( dir=XUP; dir <= TUP; dir++) {
	tag[OPP_DIR(dir)]=start_gather_from_temp(htmp[OPP_DIR(dir)], 
		sizeof(half_wilson_vector), OPP_DIR(dir),
		parity, gen_pt[OPP_DIR(dir)] );
    }


	/* Set dest to zero */
        /* Take Wilson projection for src displaced in up direction, gathered,
		multiply it by link matrix, expand it, and add.
		to dest */
    for( dir=XUP; dir <= TUP; dir++) {
	wait_gather(tag[dir]);
    }
    FORSOMEPARITY(i,s,parity){
      linkx = &t_links[4*i+XUP];
      linky = &t_links[4*i+YUP];
      linkz = &t_links[4*i+ZUP];
      linkt = &t_links[4*i+TUP];
      mult_su3_mat_hwvec2( linkx, 
			   (half_wilson_vector * )(gen_pt[XUP][i]), &hwvx ); 
      mult_su3_mat_hwvec2( linky, 
			   (half_wilson_vector * )(gen_pt[YUP][i]), &hwvy ); 
      mult_su3_mat_hwvec2( linkz, 
			   (half_wilson_vector * )(gen_pt[ZUP][i]), &hwvz ); 
      mult_su3_mat_hwvec2( linkt, 
			   (half_wilson_vector * )(gen_pt[TUP][i]), &hwvt ); 
      grow_add_four_wvecs2( &(dest[i]),
	    &hwvx, &hwvy, &hwvz, &hwvt, isign, 0 ); /* "0" is NOSUM */
    }
    for( dir=XUP; dir <= TUP; dir++) {
	cleanup_gather(tag[dir]);
    }

        /* Take Wilson projection for src displaced in down direction,
        expand it, and add to dest */
    for( dir=XUP; dir <= TUP; dir++) {
	wait_gather(tag[OPP_DIR(dir)]);
    }

    FORSOMEPARITY(i,s,parity){
	grow_add_four_wvecs2( &(dest[i]),
	    (half_wilson_vector *)(gen_pt[XDOWN][i]),
	    (half_wilson_vector *)(gen_pt[YDOWN][i]),
	    (half_wilson_vector *)(gen_pt[ZDOWN][i]),
	    (half_wilson_vector *)(gen_pt[TDOWN][i]),
	    -isign, 1 );	/* "1" SUMs in current dest */
    }
    for( dir=XUP; dir <= TUP; dir++) {
	cleanup_gather(tag[OPP_DIR(dir)]);
    }

} /* end (of dslash_w_on_temp() ) */


/* Special dslash for use by congrad.  Uses restart_gather() when
  possible. Last argument is an integer, which will tell if
  gathers have been started.  If is_started=0,use
  start_gather, otherwise use restart_gather.
  Argument "tag" is a vector of a msg_tag *'s to use for
  the gathers.
  The calling program must clean up the gathers! */
void dslash_w_on_temp_special( wilson_vector *src, wilson_vector *dest,
    int isign,int parity,msg_tag **tag,int is_started)

{
half_wilson_vector hwvx,hwvy,hwvz,hwvt;

register int i;
register site *s;
register int dir,otherparity;
su3_matrix *linkx,*linky,*linkz,*linkt;
  /* allocate temporary work space only if not already allocated */
  /* The calling program must clean up this space */
  malloc_dslash_temps();
  setup_tmp_links();
  
    switch(parity) {
	case EVEN:      otherparity=ODD; break;
	case ODD:       otherparity=EVEN; break;
	case EVENANDODD:        otherparity=EVENANDODD; break;
    }

    if(N_POINTERS < 8){
      printf("dslash_w_special: N_POINTERS must be 8 or more!\n");
      terminate(1);
     }


    /* Take Wilson projection for src displaced in up direction, gather
       it to "our site" */
    FORSOMEPARITY(i,s,otherparity){
        wp_shrink_4dir2( &(src[i]), &(htmp[XUP][i]),
	    &(htmp[YUP][i]), &(htmp[ZUP][i]), &(htmp[TUP][i]), isign);
    }

    for( dir=XUP; dir <= TUP; dir++) {
	if(is_started==0)tag[dir]=start_gather_from_temp( htmp[dir],
	    sizeof(half_wilson_vector), dir, parity, gen_pt[dir] );
	else restart_gather_from_temp( htmp[dir],
	    sizeof(half_wilson_vector), dir, parity, gen_pt[dir], 
			     tag[dir] );
    }

        /* Take Wilson projection for src displaced in down direction,
        multiply it by adjoint link matrix, gather it "up" */
    FORSOMEPARITY(i,s,otherparity){
      linkx = &t_links[4*i+XUP];
      linky = &t_links[4*i+YUP];
      linkz = &t_links[4*i+ZUP];
      linkt = &t_links[4*i+TUP];
      wp_shrink_4dir2( &(src[i]), &hwvx, &hwvy, &hwvz, &hwvt, -isign);
      mult_adj_su3_mat_hwvec2( linkx, &hwvx, &(htmp[XDOWN][i]));
      mult_adj_su3_mat_hwvec2( linky, &hwvy, &(htmp[YDOWN][i]));
      mult_adj_su3_mat_hwvec2( linkz, &hwvz, &(htmp[ZDOWN][i]));
      mult_adj_su3_mat_hwvec2( linkt, &hwvt, &(htmp[TDOWN][i]));
    }

    for( dir=XUP; dir <= TUP; dir++) {
	if(is_started==0)tag[OPP_DIR(dir)]=start_gather_from_temp(
	    htmp[OPP_DIR(dir)], sizeof(half_wilson_vector),
	    OPP_DIR(dir), parity, gen_pt[OPP_DIR(dir)] );
	else restart_gather_from_temp(
	    htmp[OPP_DIR(dir)], sizeof(half_wilson_vector),
	    OPP_DIR(dir), parity, gen_pt[OPP_DIR(dir)], 
	    tag[OPP_DIR(dir)] );
    }

	/* Set dest to zero */
        /* Take Wilson projection for src displaced in up direction, gathered,
		multiply it by link matrix, expand it, and add.
		to dest */
    for( dir=XUP; dir <= TUP; dir++) {
	wait_gather(tag[dir]);
    }
    FORSOMEPARITY(i,s,parity){
      linkx = &t_links[4*i+XUP];
      linky = &t_links[4*i+YUP];
      linkz = &t_links[4*i+ZUP];
      linkt = &t_links[4*i+TUP];
      mult_su3_mat_hwvec2( linkx, 
			   (half_wilson_vector * )(gen_pt[XUP][i]), &hwvx ); 
      mult_su3_mat_hwvec2( linky, 
			   (half_wilson_vector * )(gen_pt[YUP][i]), &hwvy ); 
      mult_su3_mat_hwvec2( linkz, 
			   (half_wilson_vector * )(gen_pt[ZUP][i]), &hwvz ); 
      mult_su3_mat_hwvec2( linkt, 
			   (half_wilson_vector * )(gen_pt[TUP][i]), &hwvt ); 
      grow_add_four_wvecs2( &(dest[i]),
	    &hwvx, &hwvy, &hwvz, &hwvt, isign, 0 ); /* "0" is NOSUM */
    }

        /* Take Wilson projection for src displaced in down direction,
        expand it, and add to dest */
    for( dir=XUP; dir <= TUP; dir++) {
	wait_gather(tag[OPP_DIR(dir)]);
    }

    FORSOMEPARITY(i,s,parity){
	grow_add_four_wvecs2( &(dest[i]),
	    (half_wilson_vector *)(gen_pt[XDOWN][i]),
	    (half_wilson_vector *)(gen_pt[YDOWN][i]),
	    (half_wilson_vector *)(gen_pt[ZDOWN][i]),
	    (half_wilson_vector *)(gen_pt[TDOWN][i]),
	    -isign, 1 );	/* "1" SUMs in current dest */
    }

} /* end (of dslash_w_on_temp_special() ) */

