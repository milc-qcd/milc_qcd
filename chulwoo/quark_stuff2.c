#include <stdio.h>
#include <stdlib.h>
#include "generic_ks_includes.h"
#define NOWHERE -1

#define GOES_FORWARDS(dir) (dir<=TUP)
#define GOES_BACKWARDS(dir) (dir>TUP)
#ifdef XLC
#define cache_touch(A) __dcbt(A)
#define outReal(A) 
#else
#define cache_touch(A) asm ( "dcbt %0 , %1" : : "r" (0)  , "r" (A) )
#endif
#if defined(PPC440QCDOC) || defined(XLC)
#define PRINT(x)
#define PRINT2(x,y)
#else
#define PRINT(x) printf((x))
#define PRINT2(x,y) printf((x),(y))
#endif


void all_mult_adj_su3_mat_vec( int dir, su3_vector * source, su3_vector *result,
int *nei, int local, int nsites, site *lat){
  int i;
  site *s;
  register su3_vector *b, *bn, *c, *cn;
  register su3_matrix *a;
  register Real c0r,c0i,c1r,c1i,c2r,c2i;
  register Real b0r,b0i,b1r,b1i,b2r,b2i;
  register Real a0r,a0i,a1r,a1i,a2r,a2i;
  i=0;
  s = lat;
  b = source;
  c = result;
  if(local)
    while(nei[i]==NOWHERE&&i<nsites) {i++;b++;s++;c++;}
  else
    while(nei[i]!=NOWHERE&&i<nsites) {i++;b++;s++;c++;}
  if (i >=nsites) return;;
  a = &(s->link[dir]);
  cache_touch(&(a->e[0][0].real));
  cache_touch(&(a->e[1][1].real));
  cache_touch(&(a->e[2][2].real));
  cn = c;
  bn = b;
  
  b0r=b->c[0].real; b0i=b->c[0].imag;
  b1r=b->c[1].real; b1i=b->c[1].imag;
  a0r=a->e[0][0].real; a0i=a->e[0][0].imag;
  a1r=a->e[0][1].real; a1i=a->e[0][1].imag;
  a2r=a->e[0][2].real; a2i=a->e[0][2].imag;

  i++;s++;bn++;cn++;
  if(local)
    while(nei[i]==NOWHERE&&i<nsites) {i++;s++;bn++;cn++;}
  else
    while(nei[i]!=NOWHERE&&i<nsites) {i++;s++;bn++;cn++;}

  while(i<nsites){
    b2r=b->c[2].real; b2i=b->c[2].imag;
    c0r = a0r*b0r + a0i*b0i;
    c0i = a0r*b0i - a0i*b0r;
    a0r=a->e[1][0].real;
    a0i=a->e[1][0].imag;
    c1r = a1r*b0r + a1i*b0i;
    c1i = a1r*b0i - a1i*b0r;
    a1r=a->e[1][1].real;
    a1i=a->e[1][1].imag;
    c2r = a2r*b0r + a2i*b0i;
    c2i = a2r*b0i - a2i*b0r;
    a2r=a->e[1][2].real;
    a2i=a->e[1][2].imag;
  
    b0r=bn->c[0].real;    b0i=bn->c[0].imag;
  cache_touch(&((s->link[dir]).e[0][0].real));
  cache_touch(&((s->link[dir]).e[0][2].imag));

    c0r += a0r*b1r + a0i*b1i;
    c0i += a0r*b1i - a0i*b1r;
    a0i=a->e[2][0].imag;
    a0r=a->e[2][0].real;
    c1r += a1r*b1r;
    c1i += a1r*b1i;
    a1r=a->e[2][1].real;
    c2r += a2r*b1r;
    c2i += a2r*b1i;
    a2r=a->e[2][2].real;
    c1r += a1i*b1i;
    c1i -= a1i*b1r;
    a1i=a->e[2][1].imag;
    c2r += a2i*b1i;
    c2i -= a2i*b1r;
    a2i=a->e[2][2].imag;
    b1r=bn->c[1].real;    b1i=bn->c[1].imag;

    a = &(s->link[dir]);
  cache_touch(&(a->e[1][1].real));

    c0r += a0r*b2r;
    c0i += a0r*b2i;
    a0r=a->e[0][0].real;
    c1r += a1r*b2r;
    c1i += a1r*b2i;
    a1r=a->e[0][1].real;
    c2r += a2r*b2r;
    c2i += a2r*b2i;
    a2r=a->e[0][2].real;
    c0r += a0i*b2i;
    c0i -= a0i*b2r;
    a0i=a->e[0][0].imag;
    c1r += a1i*b2i;
    c1i -= a1i*b2r;
    a1i=a->e[0][1].imag;
    c2r += a2i*b2i;
    c2i -= a2i*b2r;
    a2i=a->e[0][2].imag;
  cache_touch(&(a->e[2][2].real));

    c->c[0].real = c0r;
    c->c[0].imag = c0i;
    c->c[1].real = c1r;
    c->c[1].imag = c1i;
    c->c[2].real = c2r;
    c->c[2].imag = c2i;

    c = cn;
    b = bn;
    i++;s++;bn++;cn++;
    if(local)
      while(nei[i]==NOWHERE&&i<nsites) {i++;s++;bn++;cn++;}
    else
      while(nei[i]!=NOWHERE&&i<nsites) {i++;s++;bn++;cn++;}
  }

    b2r=b->c[2].real;    b2i=b->c[2].imag;
    c0r = a0r*b0r;
    c0i = a0r*b0i;
    a0r=a->e[1][0].real;
    c1r = a1r*b0r;
    c1i = a1r*b0i;
    a1r=a->e[1][1].real;
    c2r = a2r*b0r;
    c2i = a2r*b0i;
    a2r=a->e[1][2].real;
  
    c0r += a0i*b0i;
    c0i -= a0i*b0r;
    a0i=a->e[1][0].imag;
    c1r += a1i*b0i;
    c1i -= a1i*b0r;
    a1i=a->e[1][1].imag;
    c2r += a2i*b0i;
    c2i -= a2i*b0r;
    a2i=a->e[1][2].imag;

    c0r += a0r*b1r;
    c0i += a0r*b1i;
    a0r=a->e[2][0].real;
    c1r += a1r*b1r;
    c1i += a1r*b1i;
    a1r=a->e[2][1].real;
    c2r += a2r*b1r;
    c2i += a2r*b1i;
    a2r=a->e[2][2].real;
    c0r += a0i*b1i;
    c0i -= a0i*b1r;
    a0i=a->e[2][0].imag;
    c1r += a1i*b1i;
    c1i -= a1i*b1r;
    a1i=a->e[2][1].imag;
    c2r += a2i*b1i;
    c2i -= a2i*b1r;
    a2i=a->e[2][2].imag;

    c0r += a0r*b2r;
    c0i += a0r*b2i;
    c1r += a1r*b2r;
    c1i += a1r*b2i;
    c2r += a2r*b2r;
    c2i += a2r*b2i;
    c0r += a0i*b2i;
    c0i -= a0i*b2r;
    c1r += a1i*b2i;
    c1i -= a1i*b2r;
    c2r += a2i*b2i;
    c2i -= a2i*b2r;

    c->c[0].real = c0r;
    c->c[0].imag = c0i;
    c->c[1].real = c1r;
    c->c[1].imag = c1i;
    c->c[2].real = c2r;
    c->c[2].imag = c2i;
}

void all_add_force_to_mom(su3_vector *back,su3_vector *forw,int dir,Real coeff, Real onehalf, Real onethird, int neven, int nsites, site *lat){
  register anti_hermitmat *c,*cn;
  register su3_vector *a, *b;
  register Real a0r,a0i,a1r,a1i,a2r,a2i;
  register Real b0r,b0i,b1r,b1i,b2r,b2i;
  register Real c00i,c11i,c22i,c01r,c01i,c02r,c02i,c12r,c12i;
  register Real temp0,temp1,temp2,temp3;
  register Real tmp_coeff, tmp_coeff2;
  register site *s;
  int i;

  if(GOES_BACKWARDS(dir))
  {
    dir = OPP_DIR(dir) ;
    coeff = -coeff ;
  }
  s = lat;
  a = &(back[0]);
  b = &(forw[0]);
  c = &(s->mom[dir]);
  tmp_coeff = coeff;
  tmp_coeff2 = tmp_coeff*onehalf;

  a0r=a->c[0].real; a0i=a->c[0].imag;
  a1r=a->c[1].real; a1i=a->c[1].imag;
  a2r=a->c[2].real; a2i=a->c[2].imag;
  b0r=b->c[0].real; b0i=b->c[0].imag;
  b1r=b->c[1].real; b1i=b->c[1].imag;
  b2r=b->c[2].real; b2i=b->c[2].imag;
  c01r = c->m01.real; c01i = c->m01.imag;
  c00i = c->m00im;
  c11i = c->m11im;
  c22i = c->m22im;

  i=0;
  s++;a++;b++;
  while(i<neven){
    cn = &(s->mom[dir]);
    cache_touch(&(cn->m01.real));
    cache_touch(&(cn->m22im));
    cache_touch(&(a->c[0].real));
    cache_touch(&(b->c[0].real));
    cache_touch(&(a->c[1].imag));
    cache_touch(&(b->c[1].imag));
    c02r = c->m02.real; c02i = c->m02.imag;
    c12r = c->m12.real; c12i = c->m12.imag;
    temp0=a0i*b0r-a0r*b0i;
    temp1=a1i*b1r-a1r*b1i;
    temp2=a2i*b2r-a2r*b2i;
    temp3=(temp0+temp1+temp2)*onethird;
    c00i += tmp_coeff*(temp0-temp3);
    c->m00im = c00i;
    c11i += tmp_coeff*(temp1-temp3);
    c->m11im = c11i;
    c22i += tmp_coeff*(temp2-temp3);
    c->m22im = c22i;
    c00i = cn->m00im;
    c11i = cn->m11im;
    c22i = cn->m22im;

    cache_touch(&(a->c[2].imag));
    cache_touch(&(b->c[2].imag));
    temp0=a0r*b1r+a0i*b1i;
    temp1=a0i*b1r-a0r*b1i;
    temp2=a1r*b0r+a1i*b0i;
    temp3=a1i*b0r-a1r*b0i;
    c01r += tmp_coeff2*(temp0-temp2);
    c01i += tmp_coeff2*(temp1+temp3);
    c->m01.real = c01r;
    c->m01.imag = c01i;
    c01r = cn->m01.real; c01i = cn->m01.imag;

    temp0=a0r*b2r+a0i*b2i;
    temp1=a0i*b2r-a0r*b2i;
    temp2=a2r*b0r+a2i*b0i;
    temp3=a2i*b0r-a2r*b0i;
    a0r=a->c[0].real; a0i=a->c[0].imag;
    b0r=b->c[0].real; b0i=b->c[0].imag;
    c02r += tmp_coeff2*(temp0-temp2);
    c02i += tmp_coeff2*(temp1+temp3);
    c->m02.real = c02r;
    c->m02.imag = c02i;

    temp0=a1r*b2r+a1i*b2i;
    temp1=a1i*b2r-a1r*b2i;
    temp2=a2r*b1r+a2i*b1i;
    temp3=a2i*b1r-a2r*b1i;
    a1r=a->c[1].real; a1i=a->c[1].imag;
    a2r=a->c[2].real; a2i=a->c[2].imag;
    b1r=b->c[1].real; b1i=b->c[1].imag;
    b2r=b->c[2].real; b2i=b->c[2].imag;
    c12r += tmp_coeff2*(temp0-temp2);
    c12i += tmp_coeff2*(temp1+temp3);
    c->m12.real = c12r;
    c->m12.imag = c12i;

#if 0
    PRINT2("i=%d\n",i);
    outReal(lat[i].mom[dir].m00im);
    outReal(lat[i].mom[dir].m11im);
    outReal(lat[i].mom[dir].m22im);
    outReal(lat[i].mom[dir].m01.real);
    outReal(lat[i].mom[dir].m01.imag);
    outReal(lat[i].mom[dir].m02.real);
    outReal(lat[i].mom[dir].m02.imag);
    outReal(lat[i].mom[dir].m12.real);
    outReal(lat[i].mom[dir].m12.imag);
#endif

    c=cn;
    i++;s++;a++;b++;
  }
  tmp_coeff = -coeff;
  tmp_coeff2 = tmp_coeff*onehalf;
  while(i<nsites){
    cache_touch(&(cn->m01.real));
    cache_touch(&(cn->m22im));
    cache_touch(&(a->c[0].real));
    cache_touch(&(b->c[0].real));
    cache_touch(&(a->c[1].imag));
    cache_touch(&(b->c[1].imag));
    cn = &(s->mom[dir]);
    c02r = c->m02.real; c02i = c->m02.imag;
    c12r = c->m12.real; c12i = c->m12.imag;
    temp0=a0i*b0r-a0r*b0i;
    temp1=a1i*b1r-a1r*b1i;
    temp2=a2i*b2r-a2r*b2i;
    temp3=(temp0+temp1+temp2)*onethird;
    c00i += tmp_coeff*(temp0-temp3);
    c->m00im = c00i;
    c00i = cn->m00im;
    c11i += tmp_coeff*(temp1-temp3);
    c->m11im = c11i;
    c11i = cn->m11im;
    c22i += tmp_coeff*(temp2-temp3);
    c->m22im = c22i;
    c22i = cn->m22im;

    cache_touch(&(a->c[2].imag));
    cache_touch(&(b->c[2].imag));
    temp0=a0r*b1r+a0i*b1i;
    temp1=a0i*b1r-a0r*b1i;
    temp2=a1r*b0r+a1i*b0i;
    temp3=a1i*b0r-a1r*b0i;
    c01r += tmp_coeff2*(temp0-temp2);
    c01i += tmp_coeff2*(temp1+temp3);
    c->m01.real = c01r;
    c->m01.imag = c01i;
    c01r = cn->m01.real; c01i = cn->m01.imag;

    temp0=a0r*b2r+a0i*b2i;
    temp1=a0i*b2r-a0r*b2i;
    temp2=a2r*b0r+a2i*b0i;
    temp3=a2i*b0r-a2r*b0i;
    a0r=a->c[0].real; a0i=a->c[0].imag;
    b0r=b->c[0].real; b0i=b->c[0].imag;
    c02r += tmp_coeff2*(temp0-temp2);
    c02i += tmp_coeff2*(temp1+temp3);
    c->m02.real = c02r;
    c->m02.imag = c02i;

    temp0=a1r*b2r+a1i*b2i;
    temp1=a1i*b2r-a1r*b2i;
    temp2=a2r*b1r+a2i*b1i;
    temp3=a2i*b1r-a2r*b1i;
    a1r=a->c[1].real; a1i=a->c[1].imag;
    a2r=a->c[2].real; a2i=a->c[2].imag;
    b1r=b->c[1].real; b1i=b->c[1].imag;
    b2r=b->c[2].real; b2i=b->c[2].imag;
    c12r += tmp_coeff2*(temp0-temp2);
    c12i += tmp_coeff2*(temp1+temp3);
    c->m12.real = c12r;
    c->m12.imag = c12i;

#if 0
    PRINT2("i=%d\n",i);
    outReal(lat[i].mom[dir].m00im);
    outReal(lat[i].mom[dir].m11im);
    outReal(lat[i].mom[dir].m22im);
    outReal(lat[i].mom[dir].m01.real);
    outReal(lat[i].mom[dir].m01.imag);
    outReal(lat[i].mom[dir].m02.real);
    outReal(lat[i].mom[dir].m02.imag);
    outReal(lat[i].mom[dir].m12.real);
    outReal(lat[i].mom[dir].m12.imag);
#endif

    c=cn;
    i++;s++;a++;b++;
  }
    c02r = c->m02.real; c02i = c->m02.imag;
    c12r = c->m12.real; c12i = c->m12.imag;
    temp0=a0i*b0r-a0r*b0i;
    temp1=a1i*b1r-a1r*b1i;
    temp2=a2i*b2r-a2r*b2i;
    temp3=(temp0+temp1+temp2)*onethird;
    c00i += tmp_coeff*(temp0-temp3);
    c11i += tmp_coeff*(temp1-temp3);
    c22i += tmp_coeff*(temp2-temp3);

    temp0=a0r*b1r+a0i*b1i;
    temp1=a0i*b1r-a0r*b1i;
    temp2=a1r*b0r+a1i*b0i;
    temp3=a1i*b0r-a1r*b0i;
    c01r += tmp_coeff2*(temp0-temp2);
    c01i += tmp_coeff2*(temp1+temp3);

    temp0=a0r*b2r+a0i*b2i;
    temp1=a0i*b2r-a0r*b2i;
    temp2=a2r*b0r+a2i*b0i;
    temp3=a2i*b0r-a2r*b0i;
    c02r += tmp_coeff2*(temp0-temp2);
    c02i += tmp_coeff2*(temp1+temp3);

    temp0=a1r*b2r+a1i*b2i;
    temp1=a1i*b2r-a1r*b2i;
    temp2=a2r*b1r+a2i*b1i;
    temp3=a2i*b1r-a2r*b1i;
    c12r += tmp_coeff2*(temp0-temp2);
    c12i += tmp_coeff2*(temp1+temp3);
    c->m00im = c00i;
    c->m11im = c11i;
    c->m22im = c22i;
    c->m01.real = c01r;
    c->m01.imag = c01i;
    c->m02.real = c02r;
    c->m02.imag = c02i;
    c->m12.real = c12r;
    c->m12.imag = c12i;

#if 0
    PRINT2("i=%d\n",i);
    outReal(lat[i].mom[dir].m00im);
    outReal(lat[i].mom[dir].m11im);
    outReal(lat[i].mom[dir].m22im);
    outReal(lat[i].mom[dir].m01.real);
    outReal(lat[i].mom[dir].m01.imag);
    outReal(lat[i].mom[dir].m02.real);
    outReal(lat[i].mom[dir].m02.imag);
    outReal(lat[i].mom[dir].m12.real);
    outReal(lat[i].mom[dir].m12.imag);
#endif
  
#if 0
     PRINT("add_force_to_mom()\n");
  outReal(lat[0].mom[dir].m00im);
  PRINT2("%0.8e\n",lat[0].mom[dir].m00im);
  outReal(lat[0].mom[dir].m11im);
  PRINT2("%0.8e\n",lat[0].mom[dir].m11im);
  outReal(lat[0].mom[dir].m22im);
  PRINT2("%0.8e\n",lat[0].mom[dir].m22im);
  outReal(lat[0].mom[dir].m01.real);
  PRINT2("%0.8e\n",lat[0].mom[dir].m01.real);
  outReal(lat[0].mom[dir].m01.imag);
  PRINT2("%0.8e\n",lat[0].mom[dir].m01.imag);
#if 0
  outReal(lat[0].mom[dir].m02.real);
  PRINT2("%0.8e\n",lat[0].mom[dir].m02.real);
  outReal(lat[0].mom[dir].m02.imag);
  PRINT2("%0.8e\n",lat[0].mom[dir].m02.imag);
  outReal(lat[0].mom[dir].m12.real);
  PRINT2("%0.8e\n",lat[0].mom[dir].m12.real);
  outReal(lat[0].mom[dir].m12.imag);
  PRINT2("%0.8e\n",lat[0].mom[dir].m12.imag);
#endif
#endif

}

void all_mult_su3_mat_vec( int dir, su3_vector ** source, su3_vector *result,
int *nei, int local, int nsites, site *lat){
  int i;
  site *s;
  register su3_vector *b, *bn, *c, *cn;
  register su3_matrix *a;
  register Real c0r,c0i,c1r,c1i,c2r,c2i;
  register Real b0r,b0i,b1r,b1i,b2r,b2i;
  register Real a0r,a0i,a1r,a1i,a2r,a2i;
  i=0;
  s = lat;
  c = result;
  if(local)
    while(nei[i]==NOWHERE&&i<nsites) {i++;s++;c++;}
  else
    while(nei[i]!=NOWHERE&&i<nsites) {i++;s++;c++;}
  if (i >=nsites) return;;
  a = &(s->link[dir]);
  b = (su3_vector *)(source[i]);
  cache_touch(&(a->e[0][0].real));
  cache_touch(&(a->e[1][1].real));
  cache_touch(&(a->e[2][2].real));
  cn = c;
  bn = b;
  
  b0r=b->c[0].real; b0i=b->c[0].imag;
  b1r=b->c[1].real; b1i=b->c[1].imag;
  a0r=a->e[0][0].real; a0i=a->e[0][0].imag;
  a1r=a->e[1][0].real; a1i=a->e[1][0].imag;
  a2r=a->e[2][0].real; a2i=a->e[2][0].imag;

  i++;s++;;cn++;
  if(local)
    while(nei[i]==NOWHERE&&i<nsites) {i++;s++;cn++;}
  else
    while(nei[i]!=NOWHERE&&i<nsites) {i++;s++;cn++;}
  bn = (su3_vector *)(source[i]);

  while(i<nsites){
  
    b2r=b->c[2].real; b2i=b->c[2].imag;
    c0r = a0r*b0r - a0i*b0i;
    c0i = a0r*b0i + a0i*b0r;
    a0r=a->e[0][1].real;
    a0i=a->e[0][1].imag;
    c1r = a1r*b0r - a1i*b0i;
    c1i = a1r*b0i + a1i*b0r;
    a1r=a->e[1][1].real;
    a1i=a->e[1][1].imag;
    c2r = a2r*b0r - a2i*b0i;
    c2i = a2r*b0i + a2i*b0r;
    a2r=a->e[2][1].real;
    a2i=a->e[2][1].imag;
  
    b0r=bn->c[0].real;    b0i=bn->c[0].imag;
  cache_touch(&((s->link[dir]).e[0][0].real));
  cache_touch(&((s->link[dir]).e[1][1].real));
  cache_touch(&((s->link[dir]).e[2][2].real));

    c0r += a0r*b1r - a0i*b1i;
    c0i += a0r*b1i + a0i*b1r;
    a0i=a->e[0][2].imag;
    a0r=a->e[0][2].real;
    c1r += a1r*b1r;
    c1i += a1r*b1i;
    a1r=a->e[1][2].real;
    c2r += a2r*b1r;
    c2i += a2r*b1i;
    a2r=a->e[2][2].real;
    c1r -= a1i*b1i;
    c1i += a1i*b1r;
    a1i=a->e[1][2].imag;
    c2r -= a2i*b1i;
    c2i += a2i*b1r;
    a2i=a->e[2][2].imag;
    b1r=bn->c[1].real;    b1i=bn->c[1].imag;

    a = &(s->link[dir]);

    c0r += a0r*b2r;
    c0i += a0r*b2i;
    a0r=a->e[0][0].real;
    c1r += a1r*b2r;
    c1i += a1r*b2i;
    a1r=a->e[1][0].real;
    c2r += a2r*b2r;
    c2i += a2r*b2i;
    a2r=a->e[2][0].real;
    c0r -= a0i*b2i;
    c0i += a0i*b2r;
    a0i=a->e[0][0].imag;
    c1r -= a1i*b2i;
    c1i += a1i*b2r;
    a1i=a->e[1][0].imag;
    c2r -= a2i*b2i;
    c2i += a2i*b2r;
    a2i=a->e[2][0].imag;

    c->c[0].real = c0r;
    c->c[0].imag = c0i;
    c->c[1].real = c1r;
    c->c[1].imag = c1i;
    c->c[2].real = c2r;
    c->c[2].imag = c2i;
    c = cn;
    b = bn;

    i++;s++;cn++;
    if(local)
      while(nei[i]==NOWHERE&&i<nsites) {i++;s++;cn++;}
    else
      while(nei[i]!=NOWHERE&&i<nsites) {i++;s++;cn++;}
    bn = (su3_vector *)(source[i]);
  }
    b2r=b->c[2].real; b2i=b->c[2].imag;
    c0r = a0r*b0r - a0i*b0i;
    c0i = a0r*b0i + a0i*b0r;
    a0r=a->e[0][1].real;
    a0i=a->e[0][1].imag;
    c1r = a1r*b0r - a1i*b0i;
    c1i = a1r*b0i + a1i*b0r;
    a1r=a->e[1][1].real;
    a1i=a->e[1][1].imag;
    c2r = a2r*b0r - a2i*b0i;
    c2i = a2r*b0i + a2i*b0r;
    a2r=a->e[2][1].real;
    a2i=a->e[2][1].imag;
  
  cache_touch(&((s->link[dir]).e[0][0].real));
  cache_touch(&((s->link[dir]).e[1][1].real));
  cache_touch(&((s->link[dir]).e[2][2].real));

    c0r += a0r*b1r - a0i*b1i;
    c0i += a0r*b1i + a0i*b1r;
    a0i=a->e[0][2].imag;
    a0r=a->e[0][2].real;
    c1r += a1r*b1r;
    c1i += a1r*b1i;
    a1r=a->e[1][2].real;
    c2r += a2r*b1r;
    c2i += a2r*b1i;
    a2r=a->e[2][2].real;
    c1r -= a1i*b1i;
    c1i += a1i*b1r;
    a1i=a->e[1][2].imag;
    c2r -= a2i*b1i;
    c2i += a2i*b1r;
    a2i=a->e[2][2].imag;

    c0r += a0r*b2r;
    c0i += a0r*b2i;
    c1r += a1r*b2r;
    c1i += a1r*b2i;
    c2r += a2r*b2r;
    c2i += a2r*b2i;
    c0r -= a0i*b2i;
    c0i += a0i*b2r;
    c1r -= a1i*b2i;
    c1i += a1i*b2r;
    c2r -= a2i*b2i;
    c2i += a2i*b2r;

    c->c[0].real = c0r;
    c->c[0].imag = c0i;
    c->c[1].real = c1r;
    c->c[1].imag = c1i;
    c->c[2].real = c2r;
    c->c[2].imag = c2i;
}
