/************************* ff_opt.c  *************************************/
/* MIMD version 7 */
#include "generic_ks_includes.h"

#define GOES_FORWARDS(dir) (dir<=TUP)
#define GOES_BACKWARDS(dir) (dir>TUP)
#ifdef XLC
#define cache_touch(A) __dcbt(A)
#else
/*#define cache_touch(A) asm ( "dcbt %0 , %1" : : "r" (0)  , "r" (A) )*/
#define cache_touch(A)
#endif


void mult_adj_su3_fieldlink_lathwvec( su3_matrix *link,
				      half_wilson_vector **src_pt, 
				      half_wilson_vector *dest)
{
  int i;

  half_wilson_vector **b, *c;
  su3_matrix *a;
  Real a0r,a0i,a1r,a1i,a2r,a2i;
  Real b0r,b0i,b1r,b1i,b2r,b2i;
  
  for(i = 0, a = link, b = src_pt, c = dest; 
      i < sites_on_node; 
      i++, b++, c++, a++)
    {
      if(i < sites_on_node - 1){
	su3_matrix *an = a + 1;
	half_wilson_vector **bn = b + 1;

	cache_touch(&(an->e[0][0].real));
	cache_touch(&(an->e[1][1].real));
	cache_touch(&(an->e[2][2].real));
	cache_touch(&((*bn)->h[0].c[0].real));
	cache_touch(&((*bn)->h[1].c[0].real));
      }
#if 0
      mult_adj_su3_mat_hwvec( a, *b, c );
#else
      a0r=a->e[0][0].real;      a0i=a->e[0][0].imag;
      b0r=(*b)->h[0].c[0].real; b0i=(*b)->h[0].c[0].imag;
      a1r=a->e[1][0].real;      a1i=a->e[1][0].imag;
      b1r=(*b)->h[0].c[1].real; b1i=(*b)->h[0].c[1].imag;
      a2r=a->e[2][0].real;      a2i=a->e[2][0].imag;
      b2r=(*b)->h[0].c[2].real; b2i=(*b)->h[0].c[2].imag;
      
      c->h[0].c[0].real = a0r*b0r + a0i*b0i + a1r*b1r 
	+ a1i*b1i + a2r*b2r + a2i*b2i;
      c->h[0].c[0].imag = a0r*b0i - a0i*b0r + a1r*b1i 
	- a1i*b1r + a2r*b2i - a2i*b2r;
      
      a0r=a->e[0][1].real;      a0i=a->e[0][1].imag;
      b0r=(*b)->h[0].c[0].real; b0i=(*b)->h[0].c[0].imag;  
      a1r=a->e[1][1].real;      a1i=a->e[1][1].imag;
      b1r=(*b)->h[0].c[1].real; b1i=(*b)->h[0].c[1].imag;
      a2r=a->e[2][1].real;      a2i=a->e[2][1].imag;
      b2r=(*b)->h[0].c[2].real; b2i=(*b)->h[0].c[2].imag;
      
      c->h[0].c[1].real = a0r*b0r + a0i*b0i + a1r*b1r 
	+ a1i*b1i + a2r*b2r + a2i*b2i;
      c->h[0].c[1].imag = a0r*b0i - a0i*b0r + a1r*b1i 
	- a1i*b1r + a2r*b2i - a2i*b2r;
      
      a0r=a->e[0][2].real;      a0i=a->e[0][2].imag;
      b0r=(*b)->h[0].c[0].real; b0i=(*b)->h[0].c[0].imag;
      a1r=a->e[1][2].real;      a1i=a->e[1][2].imag;
      b1r=(*b)->h[0].c[1].real; b1i=(*b)->h[0].c[1].imag;
      a2r=a->e[2][2].real;      a2i=a->e[2][2].imag;
      b2r=(*b)->h[0].c[2].real; b2i=(*b)->h[0].c[2].imag;
      
      c->h[0].c[2].real = a0r*b0r + a0i*b0i + a1r*b1r 
	+ a1i*b1i + a2r*b2r + a2i*b2i;
      c->h[0].c[2].imag = a0r*b0i - a0i*b0r + a1r*b1i 
	- a1i*b1r + a2r*b2i - a2i*b2r;
      
      a0r=a->e[0][0].real;      a0i=a->e[0][0].imag;
      b0r=(*b)->h[1].c[0].real; b0i=(*b)->h[1].c[0].imag;
      a1r=a->e[1][0].real;      a1i=a->e[1][0].imag;
      b1r=(*b)->h[1].c[1].real; b1i=(*b)->h[1].c[1].imag;
      a2r=a->e[2][0].real;      a2i=a->e[2][0].imag;
      b2r=(*b)->h[1].c[2].real; b2i=(*b)->h[1].c[2].imag;
      
      c->h[1].c[0].real = a0r*b0r + a0i*b0i + a1r*b1r 
	+ a1i*b1i + a2r*b2r + a2i*b2i;
      c->h[1].c[0].imag = a0r*b0i - a0i*b0r + a1r*b1i 
	- a1i*b1r + a2r*b2i - a2i*b2r;
      
      a0r=a->e[0][1].real;      a0i=a->e[0][1].imag;
      b0r=(*b)->h[1].c[0].real; b0i=(*b)->h[1].c[0].imag;  
      a1r=a->e[1][1].real;      a1i=a->e[1][1].imag;
      b1r=(*b)->h[1].c[1].real; b1i=(*b)->h[1].c[1].imag;
      a2r=a->e[2][1].real;      a2i=a->e[2][1].imag;
      b2r=(*b)->h[1].c[2].real; b2i=(*b)->h[1].c[2].imag;
      
      c->h[1].c[1].real = a0r*b0r + a0i*b0i + a1r*b1r 
	+ a1i*b1i + a2r*b2r + a2i*b2i;
      c->h[1].c[1].imag = a0r*b0i - a0i*b0r + a1r*b1i 
	- a1i*b1r + a2r*b2i - a2i*b2r;
      
      a0r=a->e[0][2].real;      a0i=a->e[0][2].imag;
      b0r=(*b)->h[1].c[0].real; b0i=(*b)->h[1].c[0].imag;
      a1r=a->e[1][2].real;      a1i=a->e[1][2].imag;
      b1r=(*b)->h[1].c[1].real; b1i=(*b)->h[1].c[1].imag;
      a2r=a->e[2][2].real;      a2i=a->e[2][2].imag;
      b2r=(*b)->h[1].c[2].real; b2i=(*b)->h[1].c[2].imag;
      
      c->h[1].c[2].real = a0r*b0r + a0i*b0i + a1r*b1r 
	+ a1i*b1i + a2r*b2r + a2i*b2i;
      c->h[1].c[2].imag = a0r*b0i - a0i*b0r + a1r*b1i 
	- a1i*b1r + a2r*b2i - a2i*b2r;
#endif

    }
}

void mult_su3_sitelink_lathwvec( int dir, 
				 half_wilson_vector **src_pt, 
				 half_wilson_vector *dest)
{
  int i;
  site *s;
  half_wilson_vector **b, *c;
  su3_matrix *a;
  Real a0r,a0i,a1r,a1i,a2r,a2i;
  Real b0r,b0i,b1r,b1i,b2r,b2i;

  for(i = 0, s = lattice, b = src_pt, c = dest; 
      i < sites_on_node; 
      i++, s++, b++, c++)
    {
      a = &(s->link[dir]);

      if(i < sites_on_node - 1){
	site *sn = s + 1;
	su3_matrix *an = &(sn->link[dir]);
	half_wilson_vector **bn = b + 1;

	cache_touch(&(an->e[0][0].real));
	cache_touch(&(an->e[1][1].real));
	cache_touch(&(an->e[2][2].real));
	cache_touch(&((*bn)->h[0].c[0].real));
	cache_touch(&((*bn)->h[1].c[0].real));
      }
#if 0
      mult_su3_mat_hwvec( a, *b, c );
#else

      a0r=a->e[0][0].real;       a0i=a->e[0][0].imag;
      b0r=(*b)->h[0].c[0].real;  b0i=(*b)->h[0].c[0].imag;
      a1r=a->e[0][1].real;       a1i=a->e[0][1].imag;
      b1r=(*b)->h[0].c[1].real;  b1i=(*b)->h[0].c[1].imag;
      a2r=a->e[0][2].real;       a2i=a->e[0][2].imag;
      b2r=(*b)->h[0].c[2].real;  b2i=(*b)->h[0].c[2].imag;
      
      c->h[0].c[0].real = a0r*b0r - a0i*b0i + a1r*b1r 
	- a1i*b1i + a2r*b2r - a2i*b2i;
      c->h[0].c[0].imag = a0r*b0i + a0i*b0r + a1r*b1i 
	+ a1i*b1r + a2r*b2i + a2i*b2r;
      
      a0r=a->e[1][0].real;       a0i=a->e[1][0].imag;
      b0r=(*b)->h[0].c[0].real;  b0i=(*b)->h[0].c[0].imag;
      a1r=a->e[1][1].real;       a1i=a->e[1][1].imag;
      b1r=(*b)->h[0].c[1].real;  b1i=(*b)->h[0].c[1].imag;
      a2r=a->e[1][2].real;       a2i=a->e[1][2].imag;
      b2r=(*b)->h[0].c[2].real;  b2i=(*b)->h[0].c[2].imag;
      
      c->h[0].c[1].real = a0r*b0r - a0i*b0i + a1r*b1r 
	- a1i*b1i + a2r*b2r - a2i*b2i;
      c->h[0].c[1].imag = a0r*b0i + a0i*b0r + a1r*b1i 
	+ a1i*b1r + a2r*b2i + a2i*b2r;
      
      a0r=a->e[2][0].real;       a0i=a->e[2][0].imag;
      b0r=(*b)->h[0].c[0].real;  b0i=(*b)->h[0].c[0].imag;
      a1r=a->e[2][1].real;       a1i=a->e[2][1].imag;
      b1r=(*b)->h[0].c[1].real;  b1i=(*b)->h[0].c[1].imag;
      a2r=a->e[2][2].real;       a2i=a->e[2][2].imag;
      b2r=(*b)->h[0].c[2].real;  b2i=(*b)->h[0].c[2].imag;
      
      c->h[0].c[2].real = a0r*b0r - a0i*b0i + a1r*b1r 
	- a1i*b1i + a2r*b2r - a2i*b2i;
      c->h[0].c[2].imag = a0r*b0i + a0i*b0r + a1r*b1i 
	+ a1i*b1r + a2r*b2i + a2i*b2r;
      
      a0r=a->e[0][0].real;       a0i=a->e[0][0].imag;
      b0r=(*b)->h[1].c[0].real;  b0i=(*b)->h[1].c[0].imag;
      a1r=a->e[0][1].real;       a1i=a->e[0][1].imag;
      b1r=(*b)->h[1].c[1].real;  b1i=(*b)->h[1].c[1].imag;
      a2r=a->e[0][2].real;       a2i=a->e[0][2].imag;
      b2r=(*b)->h[1].c[2].real;  b2i=(*b)->h[1].c[2].imag;
      
      c->h[1].c[0].real = a0r*b0r - a0i*b0i + a1r*b1r 
	- a1i*b1i + a2r*b2r - a2i*b2i;
      c->h[1].c[0].imag = a0r*b0i + a0i*b0r + a1r*b1i 
	+ a1i*b1r + a2r*b2i + a2i*b2r;
      
      a0r=a->e[1][0].real;       a0i=a->e[1][0].imag;
      b0r=(*b)->h[1].c[0].real;  b0i=(*b)->h[1].c[0].imag;
      a1r=a->e[1][1].real;       a1i=a->e[1][1].imag;
      b1r=(*b)->h[1].c[1].real;  b1i=(*b)->h[1].c[1].imag;
      a2r=a->e[1][2].real;       a2i=a->e[1][2].imag;
      b2r=(*b)->h[1].c[2].real;  b2i=(*b)->h[1].c[2].imag;
      
      c->h[1].c[1].real = a0r*b0r - a0i*b0i + a1r*b1r 
	- a1i*b1i + a2r*b2r - a2i*b2i;
      c->h[1].c[1].imag = a0r*b0i + a0i*b0r + a1r*b1i 
	+ a1i*b1r + a2r*b2i + a2i*b2r;
      
      a0r=a->e[2][0].real;       a0i=a->e[2][0].imag;
      b0r=(*b)->h[1].c[0].real;  b0i=(*b)->h[1].c[0].imag;
      a1r=a->e[2][1].real;       a1i=a->e[2][1].imag;
      b1r=(*b)->h[1].c[1].real;  b1i=(*b)->h[1].c[1].imag;
      a2r=a->e[2][2].real;       a2i=a->e[2][2].imag;
      b2r=(*b)->h[1].c[2].real;  b2i=(*b)->h[1].c[2].imag;
      
      c->h[1].c[2].real = a0r*b0r - a0i*b0i + a1r*b1r 
	- a1i*b1i + a2r*b2r - a2i*b2i;
      c->h[1].c[2].imag = a0r*b0i + a0i*b0r + a1r*b1i 
	+ a1i*b1r + a2r*b2i + a2i*b2r;
#endif
    }      
}

void scalar_mult_add_lathwvec_proj(su3_matrix *mom, half_wilson_vector *back, 
				   half_wilson_vector *forw, Real coeff[2])
{
  int i,j,k;
  site *s;
  Real tmp_coeff[2];
  su3_matrix tmat, *a;
  half_wilson_vector *b, *f;
  Real tmp0,tmp1;

  for(i = 0, a = mom, b = back, f = forw, s = lattice; 
      i < sites_on_node;
      i++, a++, b++, f++, s++)
  {
      if(i < sites_on_node - 1){
	su3_matrix *an = a + 1;
	half_wilson_vector *bn = b + 1;
	half_wilson_vector *fn = f + 1;

	cache_touch(&(an->e[0][0].real));
	cache_touch(&(an->e[1][1].real));
	cache_touch(&(an->e[2][2].real));
	cache_touch(&(bn->h[0].c[0].real));
	cache_touch(&(bn->h[1].c[0].real));
	cache_touch(&(fn->h[0].c[0].real));
	cache_touch(&(fn->h[1].c[0].real));
      }

    if(s->parity==ODD)
      {	tmp_coeff[0] = -coeff[0] ; tmp_coeff[1] = -coeff[1] ; }
    else
      {	tmp_coeff[0] = coeff[0] ;  tmp_coeff[1] = coeff[1] ;  }

#if SSE_INLINE
    su3_projector(&(b->h[0]), &(f->h[0]), &tmat);
    scalar_mult_add_su3_matrix(a, &tmat, tmp_coeff[0], a );
    su3_projector(&(b->h[1]), &(f->h[1]), &tmat);
    scalar_mult_add_su3_matrix(a, &tmat, tmp_coeff[1], a );
#else
#if 0
    scalar_mult_add_hwvec_proj(a, b, f, tmp_coeff, a );
#else
    for(k=0;k<3;k++)for(j=0;j<3;j++){
	tmp1 = b->h[0].c[k].real * f->h[0].c[j].real;
	tmp0 = b->h[0].c[k].imag * f->h[0].c[j].imag;
	a->e[k][j].real = a->e[k][j].real + (tmp0 + tmp1)*tmp_coeff[0];

	tmp1 = b->h[0].c[k].real * f->h[0].c[j].imag;
	tmp0 = b->h[0].c[k].imag * f->h[0].c[j].real;
	a->e[k][j].imag += (tmp0 - tmp1)*tmp_coeff[0];

	tmp1 = b->h[1].c[k].real * f->h[1].c[j].real;
	tmp0 = b->h[1].c[k].imag * f->h[1].c[j].imag;
	a->e[k][j].real += (tmp0 + tmp1)*tmp_coeff[1];

	tmp1 = b->h[1].c[k].real * f->h[1].c[j].imag;
	tmp0 = b->h[1].c[k].imag * f->h[1].c[j].real;
	a->e[k][j].imag += (tmp0 - tmp1)*tmp_coeff[1];
    }
#endif
#endif
  }
}

void scalar_mult_add_lathwvec(half_wilson_vector *dest, 
			      half_wilson_vector *src, Real s[2])
{
  int i;
  half_wilson_vector *a, *b;

  for(i = 0, a = dest, b = src; 
      i < sites_on_node; 
      i++, a++, b++)
    {
      if(i < sites_on_node - 1){
	half_wilson_vector *bn = b + 1;
	half_wilson_vector *an = a + 1;

	cache_touch(&(an->h[0].c[0].real));
	cache_touch(&(an->h[1].c[0].real));
	cache_touch(&(bn->h[0].c[0].real));
	cache_touch(&(bn->h[1].c[0].real));
      }
#if 0
      scalar_mult_add_su3_vector(&(a->h[0]), &(b->h[0]), s[0],
				 &(a->h[0]));
      scalar_mult_add_su3_vector(&(a->h[1]), &(b->h[1]), s[1],
				 &(a->h[1]));
#else

  a->h[0].c[0].real = a->h[0].c[0].real + s[0]*b->h[0].c[0].real;
  a->h[0].c[0].imag = a->h[0].c[0].imag + s[0]*b->h[0].c[0].imag;
  a->h[0].c[1].real = a->h[0].c[1].real + s[0]*b->h[0].c[1].real;
  a->h[0].c[1].imag = a->h[0].c[1].imag + s[0]*b->h[0].c[1].imag;
  a->h[0].c[2].real = a->h[0].c[2].real + s[0]*b->h[0].c[2].real;
  a->h[0].c[2].imag = a->h[0].c[2].imag + s[0]*b->h[0].c[2].imag;

  a->h[1].c[0].real = a->h[1].c[0].real + s[0]*b->h[1].c[0].real;
  a->h[1].c[0].imag = a->h[1].c[0].imag + s[0]*b->h[1].c[0].imag;
  a->h[1].c[1].real = a->h[1].c[1].real + s[0]*b->h[1].c[1].real;
  a->h[1].c[1].imag = a->h[1].c[1].imag + s[0]*b->h[1].c[1].imag;
  a->h[1].c[2].real = a->h[1].c[2].real + s[0]*b->h[1].c[2].real;
  a->h[1].c[2].imag = a->h[1].c[2].imag + s[0]*b->h[1].c[2].imag;

#endif
    }
}
/* ff_opt.c */
