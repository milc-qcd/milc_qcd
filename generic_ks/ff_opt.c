/************************* ff_opt.c  *************************************/
/* MIMD version 7 */
#include "generic_ks_includes.h"

#define GOES_FORWARDS(dir) (dir<=TUP)
#define GOES_BACKWARDS(dir) (dir>TUP)
#ifdef QCDOC
#ifdef XLC
#define CACHE_TOUCH
#define cache_touch(A) __dcbt(A)
#else
#define CACHE_TOUCH
#define cache_touch(A) asm ( "dcbt %0 , %1" : : "r" (0)  , "r" (A) )
#endif
#else
#undef CACHE_TOUCH
#define cache_touch(A)
#endif

/*#define FF_DEBUG*/

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
#ifdef CACHE_TOUCH
      if(i < sites_on_node - 1){
	su3_matrix *an = a + 1;

	cache_touch(&(an->e[0][0].real));
	cache_touch(&(an->e[1][1].real));
      }
#endif
#ifdef FF_DEBUG
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
      
#ifdef CACHE_TOUCH
      if(i < sites_on_node - 1){
	su3_matrix *an = a + 1;
	half_wilson_vector **bn = b + 1;

	cache_touch(&(an->e[2][2].real));
	cache_touch(&((*bn)->h[0].c[0].real));
	cache_touch(&((*bn)->h[1].c[0].real));
      }
#endif
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
      
#ifdef CACHE_TOUCH
      if(i < sites_on_node - 1){
	half_wilson_vector **bn = b + 1;

	cache_touch(&((*bn)->h[1].c[2].imag));
      }
#endif
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

#ifdef CACHE_TOUCH
      if(i < sites_on_node - 1){
	site *sn = s + 1;
	su3_matrix *an = &(sn->link[dir]);
	half_wilson_vector **bn = b + 1;

	cache_touch(&(an->e[0][0].real));
	cache_touch(&(an->e[1][1].real));
	cache_touch(&(an->e[2][2].real));
      }
#endif
#ifdef FF_DEBUG
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
      
#ifdef CACHE_TOUCH
      if(i < sites_on_node - 1){
	site *sn = s + 1;
	su3_matrix *an = &(sn->link[dir]);
	half_wilson_vector **bn = b + 1;

	cache_touch(&((*bn)->h[0].c[0].real));
	cache_touch(&((*bn)->h[1].c[0].real));
      }
#endif

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

void mult_su3_fieldlink_lathwvec( su3_matrix *link,
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
#ifdef CACHE_TOUCH
      if(i < sites_on_node - 1){
	su3_matrix *an = a + 1;
	half_wilson_vector **bn = b + 1;

	cache_touch(&(an->e[0][0].real));
	cache_touch(&(an->e[1][1].real));
	cache_touch(&(an->e[2][2].real));
      }
#endif
#ifdef FF_DEBUG
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
      
#ifdef CACHE_TOUCH
      if(i < sites_on_node - 1){
	su3_matrix *an = a + 1;
	half_wilson_vector **bn = b + 1;

	cache_touch(&((*bn)->h[0].c[0].real));
	cache_touch(&((*bn)->h[1].c[0].real));
      }
#endif

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

void scalar_mult_add_lathwvec_proj(anti_hermitmat *mom, 
				   half_wilson_vector *back, 
				   half_wilson_vector *forw, Real coeff[2])
{
  int i;
  site *s;
  Real tmp_coeff[2], tmp_coeffd2[2];
#ifdef FF_DEBUG
  su3_matrix tmat, tmat2;
#endif
  anti_hermitmat *a;
  half_wilson_vector *b, *f;
  Real tmp0,tmp1;
  Real trim, tmp;

  for(i = 0, a = mom, b = back, f = forw, s = lattice; 
      i < sites_on_node;
      i++, a++, b++, f++, s++)
  {
    if(i >= even_sites_on_node) /* Odd sites */
      {	tmp_coeff[0] = -coeff[0] ; tmp_coeff[1] = -coeff[1] ; } 
    else /* Even sites */
      {	tmp_coeff[0] = coeff[0] ;  tmp_coeff[1] = coeff[1] ;  }

    tmp_coeffd2[0] = tmp_coeff[0]/2;
    tmp_coeffd2[1] = tmp_coeff[1]/2;

#ifdef CACHE_TOUCH
    if(i < sites_on_node - 1){
      anti_hermitmat *an = a + 1;
      
      cache_touch(&(an->m01.real));
    }
#endif

#ifdef FF_DEBUG
    uncompress_anti_hermitian( a, &tmat2);
    su3_projector(&(b->h[0]), &(f->h[0]), &tmat);
    scalar_mult_add_su3_matrix(&tmat2, &tmat, tmp_coeff[0], &tmat2 );
    su3_projector(&(b->h[1]), &(f->h[1]), &tmat);
    scalar_mult_add_su3_matrix(&tmat2, &tmat, tmp_coeff[1], &tmat2 );
    make_anti_hermitian(&tmat2, a);
    
#else
    trim = 0;

    /* k = 0 */

    tmp1 = b->h[0].c[0].real * f->h[0].c[0].imag;
    tmp0 = b->h[0].c[0].imag * f->h[0].c[0].real;
    tmp = (tmp0 - tmp1)*tmp_coeff[0];
    
    tmp1 = b->h[1].c[0].real * f->h[1].c[0].imag;
    tmp0 = b->h[1].c[0].imag * f->h[1].c[0].real;
    tmp += (tmp0 - tmp1)*tmp_coeff[1];
    a->m00im += tmp;
    trim += tmp;
    
    
    tmp1 = b->h[0].c[0].real * f->h[0].c[1].real;
    tmp0 = b->h[0].c[0].imag * f->h[0].c[1].imag;
    tmp = (tmp0 + tmp1)*tmp_coeffd2[0];
    
    tmp1 = b->h[1].c[0].real * f->h[1].c[1].real;
    tmp0 = b->h[1].c[0].imag * f->h[1].c[1].imag;
    tmp += (tmp0 + tmp1)*tmp_coeffd2[1];
    a->m01.real += tmp;

#ifdef CACHE_TOUCH
    if(i < sites_on_node - 1){
      anti_hermitmat *an = a + 1;
      
      cache_touch(&(an->m22im));
    }
#endif

    tmp1 = b->h[0].c[0].real * f->h[0].c[1].imag;
    tmp0 = b->h[0].c[0].imag * f->h[0].c[1].real;
    tmp = (tmp0 - tmp1)*tmp_coeffd2[0];
    
    tmp1 = b->h[1].c[0].real * f->h[1].c[1].imag;
    tmp0 = b->h[1].c[0].imag * f->h[1].c[1].real;
    tmp += (tmp0 - tmp1)*tmp_coeffd2[1];
    a->m01.imag += tmp;
    
    
    tmp1 = b->h[0].c[0].real * f->h[0].c[2].real;
    tmp0 = b->h[0].c[0].imag * f->h[0].c[2].imag;
    tmp = (tmp0 + tmp1)*tmp_coeffd2[0];
    
    tmp1 = b->h[1].c[0].real * f->h[1].c[2].real;
    tmp0 = b->h[1].c[0].imag * f->h[1].c[2].imag;
    tmp += (tmp0 + tmp1)*tmp_coeffd2[1];
    a->m02.real += tmp;
    
    tmp1 = b->h[0].c[0].real * f->h[0].c[2].imag;
    tmp0 = b->h[0].c[0].imag * f->h[0].c[2].real;
    tmp = (tmp0 - tmp1)*tmp_coeffd2[0];
    
    tmp1 = b->h[1].c[0].real * f->h[1].c[2].imag;
    tmp0 = b->h[1].c[0].imag * f->h[1].c[2].real;
    tmp += (tmp0 - tmp1)*tmp_coeffd2[1];
    a->m02.imag += tmp;
    

#ifdef CACHE_TOUCH
    if(i < sites_on_node - 1){
      half_wilson_vector *bn = b + 1;
      
      cache_touch(&(bn->h[0].c[0].real));
    }
#endif

    /* k = 1 */

    tmp1 = b->h[0].c[1].real * f->h[0].c[0].real;
    tmp0 = b->h[0].c[1].imag * f->h[0].c[0].imag;
    tmp = (tmp0 + tmp1)*tmp_coeffd2[0];
    
    tmp1 = b->h[1].c[1].real * f->h[1].c[0].real;
    tmp0 = b->h[1].c[1].imag * f->h[1].c[0].imag;
    tmp += (tmp0 + tmp1)*tmp_coeffd2[1];
    a->m01.real -= tmp;
    
    tmp1 = b->h[0].c[1].real * f->h[0].c[0].imag;
    tmp0 = b->h[0].c[1].imag * f->h[0].c[0].real;
    tmp = (tmp0 - tmp1)*tmp_coeffd2[0];
    
    tmp1 = b->h[1].c[1].real * f->h[1].c[0].imag;
    tmp0 = b->h[1].c[1].imag * f->h[1].c[0].real;
    tmp += (tmp0 - tmp1)*tmp_coeffd2[1];
    a->m01.imag += tmp;
    
    tmp1 = b->h[0].c[1].real * f->h[0].c[1].imag;
    tmp0 = b->h[0].c[1].imag * f->h[0].c[1].real;
    tmp = (tmp0 - tmp1)*tmp_coeff[0];
    
#ifdef CACHE_TOUCH
    if(i < sites_on_node - 1){
      half_wilson_vector *bn = b + 1;
      
      cache_touch(&(bn->h[1].c[2].imag));
    }
#endif

    tmp1 = b->h[1].c[1].real * f->h[1].c[1].imag;
    tmp0 = b->h[1].c[1].imag * f->h[1].c[1].real;
    tmp += (tmp0 - tmp1)*tmp_coeff[1];
    a->m11im += tmp;
    trim += tmp;
    
    
    tmp1 = b->h[0].c[1].real * f->h[0].c[2].real;
    tmp0 = b->h[0].c[1].imag * f->h[0].c[2].imag;
    tmp = (tmp0 + tmp1)*tmp_coeffd2[0];
    
    tmp1 = b->h[1].c[1].real * f->h[1].c[2].real;
    tmp0 = b->h[1].c[1].imag * f->h[1].c[2].imag;
    tmp += (tmp0 + tmp1)*tmp_coeffd2[1];
    a->m12.real += tmp;
    
    tmp1 = b->h[0].c[1].real * f->h[0].c[2].imag;
    tmp0 = b->h[0].c[1].imag * f->h[0].c[2].real;
    tmp = (tmp0 - tmp1)*tmp_coeffd2[0];
    
    tmp1 = b->h[1].c[1].real * f->h[1].c[2].imag;
    tmp0 = b->h[1].c[1].imag * f->h[1].c[2].real;
    tmp += (tmp0 - tmp1)*tmp_coeffd2[1];
    a->m12.imag += tmp;
    

#ifdef CACHE_TOUCH
    if(i < sites_on_node - 1){
      half_wilson_vector *fn = f + 1;
      
      cache_touch(&(fn->h[0].c[0].real));
    }
#endif

    /* k = 2 */

    tmp1 = b->h[0].c[2].real * f->h[0].c[0].real;
    tmp0 = b->h[0].c[2].imag * f->h[0].c[0].imag;
    tmp = (tmp0 + tmp1)*tmp_coeffd2[0];
    
    tmp1 = b->h[1].c[2].real * f->h[1].c[0].real;
    tmp0 = b->h[1].c[2].imag * f->h[1].c[0].imag;
    tmp += (tmp0 + tmp1)*tmp_coeffd2[1];
    a->m02.real -= tmp;
    
    tmp1 = b->h[0].c[2].real * f->h[0].c[0].imag;
    tmp0 = b->h[0].c[2].imag * f->h[0].c[0].real;
    tmp = (tmp0 - tmp1)*tmp_coeffd2[0];
    
    tmp1 = b->h[1].c[2].real * f->h[1].c[0].imag;
    tmp0 = b->h[1].c[2].imag * f->h[1].c[0].real;
    tmp += (tmp0 - tmp1)*tmp_coeffd2[1];
    a->m02.imag += tmp;
    
    
    tmp1 = b->h[0].c[2].real * f->h[0].c[1].real;
    tmp0 = b->h[0].c[2].imag * f->h[0].c[1].imag;
    tmp = (tmp0 + tmp1)*tmp_coeffd2[0];
    
    tmp1 = b->h[1].c[2].real * f->h[1].c[1].real;
    tmp0 = b->h[1].c[2].imag * f->h[1].c[1].imag;
    tmp += (tmp0 + tmp1)*tmp_coeffd2[1];
    a->m12.real -= tmp;
    
#ifdef CACHE_TOUCH
    if(i < sites_on_node - 1){
      half_wilson_vector *fn = f + 1;
      
      cache_touch(&(fn->h[1].c[2].imag));
    }
#endif

    tmp1 = b->h[0].c[2].real * f->h[0].c[1].imag;
    tmp0 = b->h[0].c[2].imag * f->h[0].c[1].real;
    tmp = (tmp0 - tmp1)*tmp_coeffd2[0];
    
    tmp1 = b->h[1].c[2].real * f->h[1].c[1].imag;
    tmp0 = b->h[1].c[2].imag * f->h[1].c[1].real;
    tmp += (tmp0 - tmp1)*tmp_coeffd2[1];
    a->m12.imag += tmp;
    
    
    tmp1 = b->h[0].c[2].real * f->h[0].c[2].imag;
    tmp0 = b->h[0].c[2].imag * f->h[0].c[2].real;
    tmp = (tmp0 - tmp1)*tmp_coeff[0];
    
    tmp1 = b->h[1].c[2].real * f->h[1].c[2].imag;
    tmp0 = b->h[1].c[2].imag * f->h[1].c[2].real;
    tmp += (tmp0 - tmp1)*tmp_coeff[1];
    a->m22im += tmp;
    trim += tmp;

    trim *= 0.33333333333333333;
    a->m00im -= trim;
    a->m11im -= trim;
    a->m22im -= trim;
#endif
  }
}

void scalar_mult_add_lathwvec_proj_su3mat(su3_matrix *mom, 
					  half_wilson_vector *back, 
					  half_wilson_vector *forw, 
					  Real coeff[2])
{
  int i;
  site *s;
  Real tmp_coeff[2], tmp_coeffd2[2];
#ifdef FF_DEBUG
  su3_matrix tmat;
#endif
  su3_matrix *a;
  half_wilson_vector *b, *f;
  Real tmp0,tmp1;
  Real tmp;

  for(i = 0, a = mom, b = back, f = forw, s = lattice; 
      i < sites_on_node;
      i++, a++, b++, f++, s++)
  {
    if(i >= even_sites_on_node) /* Odd sites */
      {	tmp_coeff[0] = -coeff[0] ; tmp_coeff[1] = -coeff[1] ; }
    else /* Even sites */
      {	tmp_coeff[0] = coeff[0] ;  tmp_coeff[1] = coeff[1] ;  }

    tmp_coeffd2[0] = tmp_coeff[0]/2;
    tmp_coeffd2[1] = tmp_coeff[1]/2;

#ifdef CACHE_TOUCH
    if(i < sites_on_node - 1){
      su3_matrix *an = a + 1;
      
      cache_touch(&(an->e[0][1].real));
    }
#endif

#ifdef FF_DEBUG
    su3_projector(&(b->h[0]), &(f->h[0]), a);
    scalar_mult_add_su3_matrix(a, &tmat, tmp_coeff[0], a );
    su3_projector(&(b->h[1]), &(f->h[1]), &tmat);
    scalar_mult_add_su3_matrix(a, &tmat, tmp_coeff[1], a );
    
#else
    /* k = 0 */

    a->e[0][0].real = 0;

    tmp1 = b->h[0].c[0].real * f->h[0].c[0].imag;
    tmp0 = b->h[0].c[0].imag * f->h[0].c[0].real;
    tmp = (tmp0 - tmp1)*tmp_coeff[0];
    
    tmp1 = b->h[1].c[0].real * f->h[1].c[0].imag;
    tmp0 = b->h[1].c[0].imag * f->h[1].c[0].real;
    tmp += (tmp0 - tmp1)*tmp_coeff[1];
    a->e[0][0].imag += tmp;
    
    
    tmp1 = b->h[0].c[0].real * f->h[0].c[1].real;
    tmp0 = b->h[0].c[0].imag * f->h[0].c[1].imag;
    tmp = (tmp0 + tmp1)*tmp_coeffd2[0];
    
    tmp1 = b->h[1].c[0].real * f->h[1].c[1].real;
    tmp0 = b->h[1].c[0].imag * f->h[1].c[1].imag;
    tmp += (tmp0 + tmp1)*tmp_coeffd2[1];
    a->e[0][1].real += tmp;
    a->e[1][0].real -= tmp;

#ifdef CACHE_TOUCH
    if(i < sites_on_node - 1){
      su3_matrix *an = a + 1;
      
      cache_touch(&(an->e[2][2].imag));
    }
#endif

    tmp1 = b->h[0].c[0].real * f->h[0].c[1].imag;
    tmp0 = b->h[0].c[0].imag * f->h[0].c[1].real;
    tmp = (tmp0 - tmp1)*tmp_coeffd2[0];
    
    tmp1 = b->h[1].c[0].real * f->h[1].c[1].imag;
    tmp0 = b->h[1].c[0].imag * f->h[1].c[1].real;
    tmp += (tmp0 - tmp1)*tmp_coeffd2[1];
    a->e[0][1].imag += tmp;
    a->e[1][0].imag += tmp;
    
    
    tmp1 = b->h[0].c[0].real * f->h[0].c[2].real;
    tmp0 = b->h[0].c[0].imag * f->h[0].c[2].imag;
    tmp = (tmp0 + tmp1)*tmp_coeffd2[0];
    
    tmp1 = b->h[1].c[0].real * f->h[1].c[2].real;
    tmp0 = b->h[1].c[0].imag * f->h[1].c[2].imag;
    tmp += (tmp0 + tmp1)*tmp_coeffd2[1];
    a->e[0][2].real += tmp;
    a->e[2][0].real -= tmp;
    
    tmp1 = b->h[0].c[0].real * f->h[0].c[2].imag;
    tmp0 = b->h[0].c[0].imag * f->h[0].c[2].real;
    tmp = (tmp0 - tmp1)*tmp_coeffd2[0];
    
    tmp1 = b->h[1].c[0].real * f->h[1].c[2].imag;
    tmp0 = b->h[1].c[0].imag * f->h[1].c[2].real;
    tmp += (tmp0 - tmp1)*tmp_coeffd2[1];
    a->e[0][2].imag += tmp;
    a->e[2][0].imag += tmp;
    

#ifdef CACHE_TOUCH
    if(i < sites_on_node - 1){
      half_wilson_vector *bn = b + 1;
      
      cache_touch(&(bn->h[0].c[0].real));
    }
#endif

    /* k = 1 */

    tmp1 = b->h[0].c[1].real * f->h[0].c[0].real;
    tmp0 = b->h[0].c[1].imag * f->h[0].c[0].imag;
    tmp = (tmp0 + tmp1)*tmp_coeffd2[0];
    
    tmp1 = b->h[1].c[1].real * f->h[1].c[0].real;
    tmp0 = b->h[1].c[1].imag * f->h[1].c[0].imag;
    tmp += (tmp0 + tmp1)*tmp_coeffd2[1];
    a->e[0][1].real -= tmp;
    a->e[1][0].real += tmp;
    
    tmp1 = b->h[0].c[1].real * f->h[0].c[0].imag;
    tmp0 = b->h[0].c[1].imag * f->h[0].c[0].real;
    tmp = (tmp0 - tmp1)*tmp_coeffd2[0];
    
    tmp1 = b->h[1].c[1].real * f->h[1].c[0].imag;
    tmp0 = b->h[1].c[1].imag * f->h[1].c[0].real;
    tmp += (tmp0 - tmp1)*tmp_coeffd2[1];
    a->e[0][1].imag += tmp;
    a->e[1][0].imag += tmp;
    
    a->e[1][1].real = 0;
    tmp1 = b->h[0].c[1].real * f->h[0].c[1].imag;
    tmp0 = b->h[0].c[1].imag * f->h[0].c[1].real;
    tmp = (tmp0 - tmp1)*tmp_coeff[0];
    
#ifdef CACHE_TOUCH
    if(i < sites_on_node - 1){
      half_wilson_vector *bn = b + 1;
      
      cache_touch(&(bn->h[1].c[2].imag));
    }
#endif

    tmp1 = b->h[1].c[1].real * f->h[1].c[1].imag;
    tmp0 = b->h[1].c[1].imag * f->h[1].c[1].real;
    tmp += (tmp0 - tmp1)*tmp_coeff[1];
    a->e[1][1].imag += tmp;
    
    
    tmp1 = b->h[0].c[1].real * f->h[0].c[2].real;
    tmp0 = b->h[0].c[1].imag * f->h[0].c[2].imag;
    tmp = (tmp0 + tmp1)*tmp_coeffd2[0];
    
    tmp1 = b->h[1].c[1].real * f->h[1].c[2].real;
    tmp0 = b->h[1].c[1].imag * f->h[1].c[2].imag;
    tmp += (tmp0 + tmp1)*tmp_coeffd2[1];
    a->e[1][2].real += tmp;
    a->e[2][1].real -= tmp;
    
    tmp1 = b->h[0].c[1].real * f->h[0].c[2].imag;
    tmp0 = b->h[0].c[1].imag * f->h[0].c[2].real;
    tmp = (tmp0 - tmp1)*tmp_coeffd2[0];
    
    tmp1 = b->h[1].c[1].real * f->h[1].c[2].imag;
    tmp0 = b->h[1].c[1].imag * f->h[1].c[2].real;
    tmp += (tmp0 - tmp1)*tmp_coeffd2[1];
    a->e[1][2].imag += tmp;
    a->e[2][1].imag += tmp;
    

#ifdef CACHE_TOUCH
    if(i < sites_on_node - 1){
      half_wilson_vector *fn = f + 1;
      
      cache_touch(&(fn->h[0].c[0].real));
    }
#endif

    /* k = 2 */

    tmp1 = b->h[0].c[2].real * f->h[0].c[0].real;
    tmp0 = b->h[0].c[2].imag * f->h[0].c[0].imag;
    tmp = (tmp0 + tmp1)*tmp_coeffd2[0];
    
    tmp1 = b->h[1].c[2].real * f->h[1].c[0].real;
    tmp0 = b->h[1].c[2].imag * f->h[1].c[0].imag;
    tmp += (tmp0 + tmp1)*tmp_coeffd2[1];
    a->e[0][2].real -= tmp;
    a->e[2][0].real += tmp;
    
    tmp1 = b->h[0].c[2].real * f->h[0].c[0].imag;
    tmp0 = b->h[0].c[2].imag * f->h[0].c[0].real;
    tmp = (tmp0 - tmp1)*tmp_coeffd2[0];
    
    tmp1 = b->h[1].c[2].real * f->h[1].c[0].imag;
    tmp0 = b->h[1].c[2].imag * f->h[1].c[0].real;
    tmp += (tmp0 - tmp1)*tmp_coeffd2[1];
    a->e[0][2].imag += tmp;
    a->e[2][0].imag += tmp;
    
    
    tmp1 = b->h[0].c[2].real * f->h[0].c[1].real;
    tmp0 = b->h[0].c[2].imag * f->h[0].c[1].imag;
    tmp = (tmp0 + tmp1)*tmp_coeffd2[0];
    
    tmp1 = b->h[1].c[2].real * f->h[1].c[1].real;
    tmp0 = b->h[1].c[2].imag * f->h[1].c[1].imag;
    tmp += (tmp0 + tmp1)*tmp_coeffd2[1];
    a->e[1][2].real -= tmp;
    a->e[2][1].real += tmp;
    
#ifdef CACHE_TOUCH
    if(i < sites_on_node - 1){
      half_wilson_vector *fn = f + 1;
      
      cache_touch(&(fn->h[1].c[2].imag));
    }
#endif

    tmp1 = b->h[0].c[2].real * f->h[0].c[1].imag;
    tmp0 = b->h[0].c[2].imag * f->h[0].c[1].real;
    tmp = (tmp0 - tmp1)*tmp_coeffd2[0];
    
    tmp1 = b->h[1].c[2].real * f->h[1].c[1].imag;
    tmp0 = b->h[1].c[2].imag * f->h[1].c[1].real;
    tmp += (tmp0 - tmp1)*tmp_coeffd2[1];
    a->e[1][2].imag += tmp;
    a->e[2][1].imag += tmp;
    
    
    a->e[2][2].real = 0;
    tmp1 = b->h[0].c[2].real * f->h[0].c[2].imag;
    tmp0 = b->h[0].c[2].imag * f->h[0].c[2].real;
    tmp = (tmp0 - tmp1)*tmp_coeff[0];
    
    tmp1 = b->h[1].c[2].real * f->h[1].c[2].imag;
    tmp0 = b->h[1].c[2].imag * f->h[1].c[2].real;
    tmp += (tmp0 - tmp1)*tmp_coeff[1];
    a->e[2][2].imag += tmp;

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
#ifdef CACHE_TOUCH
      if(i < sites_on_node - 1){
	half_wilson_vector *bn = b + 1;
	half_wilson_vector *an = a + 1;

	cache_touch(&(an->h[0].c[0].real));
	cache_touch(&(an->h[1].c[2].imag));
      }
#endif
#ifdef FF_DEBUG
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
      
#ifdef CACHE_TOUCH
      if(i < sites_on_node - 1){
	half_wilson_vector *bn = b + 1;
	half_wilson_vector *an = a + 1;
	
	cache_touch(&(bn->h[0].c[0].real));
	cache_touch(&(bn->h[1].c[2].imag));
      }
#endif
      
      a->h[1].c[0].real = a->h[1].c[0].real + s[0]*b->h[1].c[0].real;
      a->h[1].c[0].imag = a->h[1].c[0].imag + s[0]*b->h[1].c[0].imag;
      a->h[1].c[1].real = a->h[1].c[1].real + s[0]*b->h[1].c[1].real;
      a->h[1].c[1].imag = a->h[1].c[1].imag + s[0]*b->h[1].c[1].imag;
      a->h[1].c[2].real = a->h[1].c[2].real + s[0]*b->h[1].c[2].real;
      a->h[1].c[2].imag = a->h[1].c[2].imag + s[0]*b->h[1].c[2].imag;
      
#endif
    }
}

// Routines for type "veclist" = list of su3_vector's

void mult_adj_su3_fieldlink_latveclist( su3_matrix *link,
		veclist **src_pt, veclist *dest, int listlength ) {
	int i,j;

	veclist **b, *c;
	su3_matrix *a;
	Real a0r,a0i,a1r,a1i,a2r,a2i;
	Real b0r,b0i,b1r,b1i,b2r,b2i;

	for(i = 0, a = link, b = src_pt, c = dest; i < sites_on_node; 
			i++, b++, c++, a++){
#ifdef CACHE_TOUCH
	    if(i < sites_on_node - 1){
		su3_matrix *an = a + 1;

		cache_touch(&(an->e[0][0].real));
		cache_touch(&(an->e[1][1].real));
	    }
#endif
#ifdef FF_DEBUG
		mult_adj_su3_mat_hwvec( a, *b, c );
#else
	    for( j=0; j<listlength; j++ ) {
		a0r=a->e[0][0].real;      a0i=a->e[0][0].imag;
		b0r=(*b)->v[j].c[0].real; b0i=(*b)->v[j].c[0].imag;
		a1r=a->e[1][0].real;      a1i=a->e[1][0].imag;
		b1r=(*b)->v[j].c[1].real; b1i=(*b)->v[j].c[1].imag;
		a2r=a->e[2][0].real;      a2i=a->e[2][0].imag;
		b2r=(*b)->v[j].c[2].real; b2i=(*b)->v[j].c[2].imag;

		c->v[j].c[0].real = a0r*b0r + a0i*b0i + a1r*b1r 
			+ a1i*b1i + a2r*b2r + a2i*b2i;
		c->v[j].c[0].imag = a0r*b0i - a0i*b0r + a1r*b1i 
			- a1i*b1r + a2r*b2i - a2i*b2r;

		a0r=a->e[0][1].real;      a0i=a->e[0][1].imag;
		b0r=(*b)->v[j].c[0].real; b0i=(*b)->v[j].c[0].imag;  
		a1r=a->e[1][1].real;      a1i=a->e[1][1].imag;
		b1r=(*b)->v[j].c[1].real; b1i=(*b)->v[j].c[1].imag;
		a2r=a->e[2][1].real;      a2i=a->e[2][1].imag;
		b2r=(*b)->v[j].c[2].real; b2i=(*b)->v[j].c[2].imag;

		c->v[j].c[1].real = a0r*b0r + a0i*b0i + a1r*b1r 
			+ a1i*b1i + a2r*b2r + a2i*b2i;
		c->v[j].c[1].imag = a0r*b0i - a0i*b0r + a1r*b1i 
			- a1i*b1r + a2r*b2i - a2i*b2r;

		a0r=a->e[0][2].real;      a0i=a->e[0][2].imag;
		b0r=(*b)->v[j].c[0].real; b0i=(*b)->v[j].c[0].imag;
		a1r=a->e[1][2].real;      a1i=a->e[1][2].imag;
		b1r=(*b)->v[j].c[1].real; b1i=(*b)->v[j].c[1].imag;
		a2r=a->e[2][2].real;      a2i=a->e[2][2].imag;
		b2r=(*b)->v[j].c[2].real; b2i=(*b)->v[j].c[2].imag;

#ifdef CACHE_TOUCH
		if(i < sites_on_node - 1){
			su3_matrix *an = a + 1;
			veclist **bn = b + 1;

			cache_touch(&(an->e[2][2].real));
			cache_touch(&((*bn)->v[j].c[0].real));
		}
#endif
		c->v[j].c[2].real = a0r*b0r + a0i*b0i + a1r*b1r 
			+ a1i*b1i + a2r*b2r + a2i*b2i;
		c->v[j].c[2].imag = a0r*b0i - a0i*b0r + a1r*b1i 
			- a1i*b1r + a2r*b2i - a2i*b2r;

#endif

	    } //vector list loop
	} //site loop
}


void mult_su3_sitelink_latveclist( int dir, 
	 veclist **src_pt, veclist *dest, int listlength ) {
  int i,j;
  site *s;
  veclist **b, *c;
  su3_matrix *a;
  Real a0r,a0i,a1r,a1i,a2r,a2i;
  Real b0r,b0i,b1r,b1i,b2r,b2i;

  for(i = 0, s = lattice, b = src_pt, c = dest; 
      i < sites_on_node; i++, s++, b++, c++) {
          a = &(s->link[dir]);
    
#ifdef CACHE_TOUCH
      if(i < sites_on_node - 1){
	    site *sn = s + 1;
	    su3_matrix *an = &(sn->link[dir]);
	    veclist **bn = b + 1;
    
	    cache_touch(&(an->e[0][0].real));
	    cache_touch(&(an->e[1][1].real));
	    cache_touch(&(an->e[2][2].real));
      }
#endif
#ifdef FF_DEBUG
          mult_su3_mat_hwvec( a, *b, c );
#else
    
      for( j=0; j<listlength; j++ ) {
          a0r=a->e[0][0].real;       a0i=a->e[0][0].imag;
          b0r=(*b)->v[j].c[0].real;  b0i=(*b)->v[j].c[0].imag;
          a1r=a->e[0][1].real;       a1i=a->e[0][1].imag;
          b1r=(*b)->v[j].c[1].real;  b1i=(*b)->v[j].c[1].imag;
          a2r=a->e[0][2].real;       a2i=a->e[0][2].imag;
          b2r=(*b)->v[j].c[2].real;  b2i=(*b)->v[j].c[2].imag;
          
          c->v[j].c[0].real = a0r*b0r - a0i*b0i + a1r*b1r 
	  - a1i*b1i + a2r*b2r - a2i*b2i;
          c->v[j].c[0].imag = a0r*b0i + a0i*b0r + a1r*b1i 
	  + a1i*b1r + a2r*b2i + a2i*b2r;
          
          a0r=a->e[1][0].real;       a0i=a->e[1][0].imag;
          b0r=(*b)->v[j].c[0].real;  b0i=(*b)->v[j].c[0].imag;
          a1r=a->e[1][1].real;       a1i=a->e[1][1].imag;
          b1r=(*b)->v[j].c[1].real;  b1i=(*b)->v[j].c[1].imag;
          a2r=a->e[1][2].real;       a2i=a->e[1][2].imag;
          b2r=(*b)->v[j].c[2].real;  b2i=(*b)->v[j].c[2].imag;
          
          c->v[j].c[1].real = a0r*b0r - a0i*b0i + a1r*b1r 
	  - a1i*b1i + a2r*b2r - a2i*b2i;
          c->v[j].c[1].imag = a0r*b0i + a0i*b0r + a1r*b1i 
	  + a1i*b1r + a2r*b2i + a2i*b2r;
          
          a0r=a->e[2][0].real;       a0i=a->e[2][0].imag;
          b0r=(*b)->v[j].c[0].real;  b0i=(*b)->v[j].c[0].imag;
          a1r=a->e[2][1].real;       a1i=a->e[2][1].imag;
          b1r=(*b)->v[j].c[1].real;  b1i=(*b)->v[j].c[1].imag;
          a2r=a->e[2][2].real;       a2i=a->e[2][2].imag;
          b2r=(*b)->v[j].c[2].real;  b2i=(*b)->v[j].c[2].imag;
          
          c->v[j].c[2].real = a0r*b0r - a0i*b0i + a1r*b1r 
	  - a1i*b1i + a2r*b2r - a2i*b2i;
          c->v[j].c[2].imag = a0r*b0i + a0i*b0r + a1r*b1i 
	  + a1i*b1r + a2r*b2i + a2i*b2r;
          
#ifdef CACHE_TOUCH
          if(i < sites_on_node - 1){
		site *sn = s + 1;
		su3_matrix *an = &(sn->link[dir]);
		veclist **bn = b + 1;
		cache_touch(&((*bn)->v[j].c[0].real));
          }
#endif

#endif
        } //vector list loop
    } //site loop
}


void scalar_mult_add_latveclist_proj(anti_hermitmat *mom, 
	veclist *back, veclist *forw, Real *coeff, int listlength ) {
  int i,n;
  site *s;
  Real *tmp_coeff, *tmp_coeffd2;
#ifdef FF_DEBUG
  su3_matrix tmat, tmat2;
#endif
  anti_hermitmat *a;
  veclist *b, *f;
  Real tmp0,tmp1,tmp2,tmp3;
  Real trim, tmp;

  tmp_coeff = (Real *)malloc( listlength*sizeof(Real) );
  tmp_coeffd2 = (Real *)malloc( listlength*sizeof(Real) );
  for(i = 0, a = mom, b = back, f = forw, s = lattice; 
      i < sites_on_node; i++, a++, b++, f++, s++) { // site loop

    if(i >= even_sites_on_node) /* Odd sites */
      for( n=0; n<listlength; n++ )tmp_coeff[n] = -coeff[n];
    else /* Even sites */
      for( n=0; n<listlength; n++ )tmp_coeff[n] = coeff[n];

    for( n=0; n<listlength; n++ ) tmp_coeffd2[n] = tmp_coeff[n]/2.0;

#ifdef CACHE_TOUCH
    if(i < sites_on_node - 1){
      anti_hermitmat *an = a + 1;
      cache_touch(&(an->m01.real));
    }
#endif

#ifdef FF_DEBUG
    uncompress_anti_hermitian( a, &tmat2);
    su3_projector(&(b->v[0]), &(f->v[0]), &tmat);
    scalar_mult_add_su3_matrix(&tmat2, &tmat, tmp_coeff[0], &tmat2 );
    su3_projector(&(b->v[1]), &(f->v[1]), &tmat);
    scalar_mult_add_su3_matrix(&tmat2, &tmat, tmp_coeff[1], &tmat2 );
    make_anti_hermitian(&tmat2, a);
    
#else
    trim = 0;

    /* k = 0 */

    for( tmp=0.0,n=0; n<listlength; n++ ){
        tmp1 = b->v[n].c[0].real * f->v[n].c[0].imag;
        tmp0 = b->v[n].c[0].imag * f->v[n].c[0].real;
        tmp += (tmp0 - tmp1)*tmp_coeff[n];
    }
    a->m00im += tmp;
    trim += tmp;
    
    for( tmp2=tmp3=0.0,n=0; n<listlength; n++ ){
        tmp1 = b->v[n].c[0].real * f->v[n].c[1].real;
        tmp0 = b->v[n].c[0].imag * f->v[n].c[1].imag;
        tmp2 += (tmp0 + tmp1)*tmp_coeffd2[n];

        tmp1 = b->v[n].c[0].real * f->v[n].c[1].imag;
        tmp0 = b->v[n].c[0].imag * f->v[n].c[1].real;
        tmp3 += (tmp0 - tmp1)*tmp_coeffd2[n];
    }
    a->m01.real += tmp2;
    a->m01.imag += tmp3;
    
    
    for( tmp2=tmp3=0.0,n=0; n<listlength; n++ ){
        tmp1 = b->v[n].c[0].real * f->v[n].c[2].real;
        tmp0 = b->v[n].c[0].imag * f->v[n].c[2].imag;
        tmp2 += (tmp0 + tmp1)*tmp_coeffd2[n];
    
        tmp1 = b->v[n].c[0].real * f->v[n].c[2].imag;
        tmp0 = b->v[n].c[0].imag * f->v[n].c[2].real;
        tmp3 += (tmp0 - tmp1)*tmp_coeffd2[n];
    }
    a->m02.real += tmp2;
    a->m02.imag += tmp3;
    

#ifdef CACHE_TOUCH
    if(i < sites_on_node - 1){
      anti_hermitmat *an = a + 1;
      veclist *bn = b + 1;
      cache_touch(&(an->m22im));
      cache_touch(&(bn->v[0].c[0].real));
    }
#endif

    /* k = 1 */

    for( tmp2=tmp3=0.0,n=0; n<listlength; n++ ){
        tmp1 = b->v[n].c[1].real * f->v[n].c[0].real;
        tmp0 = b->v[n].c[1].imag * f->v[n].c[0].imag;
        tmp2 += (tmp0 + tmp1)*tmp_coeffd2[n];
    
        tmp1 = b->v[n].c[1].real * f->v[n].c[0].imag;
        tmp0 = b->v[n].c[1].imag * f->v[n].c[0].real;
        tmp3 += (tmp0 - tmp1)*tmp_coeffd2[n];
    }
    a->m01.real -= tmp2;
    a->m01.imag += tmp3;
    
    for( tmp=0.0,n=0; n<listlength; n++ ){
        tmp1 = b->v[n].c[1].real * f->v[n].c[1].imag;
        tmp0 = b->v[n].c[1].imag * f->v[n].c[1].real;
        tmp += (tmp0 - tmp1)*tmp_coeff[n];
    }
    a->m11im += tmp;
    trim += tmp;
    
    
    for( tmp2=tmp3=0.0,n=0; n<listlength; n++ ){
        tmp1 = b->v[n].c[1].real * f->v[n].c[2].real;
        tmp0 = b->v[n].c[1].imag * f->v[n].c[2].imag;
        tmp2 += (tmp0 + tmp1)*tmp_coeffd2[n];
    
        tmp1 = b->v[n].c[1].real * f->v[n].c[2].imag;
        tmp0 = b->v[n].c[1].imag * f->v[n].c[2].real;
        tmp3 += (tmp0 - tmp1)*tmp_coeffd2[n];
    }
    a->m12.real += tmp2;
    a->m12.imag += tmp3;
    

#ifdef CACHE_TOUCH
    if(i < sites_on_node - 1){
      veclist *bn = b + 1;
      veclist *fn = f + 1;
      cache_touch(&(bn->v[1].c[2].imag));
      cache_touch(&(fn->v[0].c[0].real));
    }
    if(i < sites_on_node - 1){
      veclist *fn = f + 1;
      cache_touch(&(fn->v[1].c[2].imag));
    }
#endif

    /* k = 2 */

    for( tmp2=tmp3=0.0,n=0; n<listlength; n++ ){
        tmp1 = b->v[n].c[2].real * f->v[n].c[0].real;
        tmp0 = b->v[n].c[2].imag * f->v[n].c[0].imag;
        tmp2 += (tmp0 + tmp1)*tmp_coeffd2[n];
    
        tmp1 = b->v[n].c[2].real * f->v[n].c[0].imag;
        tmp0 = b->v[n].c[2].imag * f->v[n].c[0].real;
        tmp3 += (tmp0 - tmp1)*tmp_coeffd2[n];
    }
    a->m02.real -= tmp2;
    a->m02.imag += tmp3;
    
    
    for( tmp2=tmp3=0.0,n=0; n<listlength; n++ ){
        tmp1 = b->v[n].c[2].real * f->v[n].c[1].real;
        tmp0 = b->v[n].c[2].imag * f->v[n].c[1].imag;
        tmp2 += (tmp0 + tmp1)*tmp_coeffd2[n];
    
        tmp1 = b->v[n].c[2].real * f->v[n].c[1].imag;
        tmp0 = b->v[n].c[2].imag * f->v[n].c[1].real;
        tmp3 += (tmp0 - tmp1)*tmp_coeffd2[n];
    }
    a->m12.real -= tmp2;
    a->m12.imag += tmp3;
    
    
    for( tmp=0.0,n=0; n<listlength; n++ ){
        tmp1 = b->v[n].c[2].real * f->v[n].c[2].imag;
        tmp0 = b->v[n].c[2].imag * f->v[n].c[2].real;
        tmp += (tmp0 - tmp1)*tmp_coeff[n];
    }
    a->m22im += tmp;
    trim += tmp;

    trim *= 0.33333333333333333;
    a->m00im -= trim;
    a->m11im -= trim;
    a->m22im -= trim;
#endif
  } //site loop
  free( tmp_coeff );
  free( tmp_coeffd2 );
}

void scalar_mult_add_latveclist(veclist *dest, 
	veclist *src, Real *s, int listlength ) {
  int i,j;
  veclist *a, *b;

  for(i = 0, a = dest, b = src; i < sites_on_node; 
      i++, a++, b++) { // site loop
#ifdef CACHE_TOUCH
      if(i < sites_on_node - 1){
	veclist *bn = b + 1;
	veclist *an = a + 1;

	cache_touch(&(an->v[0].c[0].real));
      }
#endif
#ifdef FF_DEBUG
      scalar_mult_add_su3_vector(&(a->v[0]), &(b->v[0]), s[0],
				 &(a->v[0]));
      scalar_mult_add_su3_vector(&(a->v[1]), &(b->v[1]), s[1],
				 &(a->v[1]));
#else

      for( j=0; j<listlength; j++ ){
          a->v[j].c[0].real = a->v[j].c[0].real + s[j]*b->v[j].c[0].real;
          a->v[j].c[0].imag = a->v[j].c[0].imag + s[j]*b->v[j].c[0].imag;
          a->v[j].c[1].real = a->v[j].c[1].real + s[j]*b->v[j].c[1].real;
          a->v[j].c[1].imag = a->v[j].c[1].imag + s[j]*b->v[j].c[1].imag;
          a->v[j].c[2].real = a->v[j].c[2].real + s[j]*b->v[j].c[2].real;
          a->v[j].c[2].imag = a->v[j].c[2].imag + s[j]*b->v[j].c[2].imag;
      }
      
#endif
    } //site loop
}
/* ff_opt.c */
