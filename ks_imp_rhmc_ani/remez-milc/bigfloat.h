/*

  Mike Clark - 25th May 2005

  bigfloat.h

  Simple C++ wrapper for multiprecision datatype used by AlgRemez
  algorithm

*/

#include <gmp.h>
#include <mpf2mpfr.h>
#include <mpfr.h>

#ifndef INCLUDED_BIGFLOAT_H
#define INCLUDED_BIGFLOAT_H

class bigfloat {

private:

  mpf_t x;

public:

  bigfloat() { mpf_init(x); }
  bigfloat(const bigfloat& y) { mpf_init_set(x, y.x); }
  bigfloat(const unsigned long u) { mpf_init_set_ui(x, u); }
  bigfloat(const long i) { mpf_init_set_si(x, i); }
  bigfloat(const int i) {mpf_init_set_si(x,(long)i);}
  bigfloat(const float d) { mpf_init_set_d(x, (double)d); }
  bigfloat(const double d) { mpf_init_set_d(x, d); }  
  bigfloat(const char *str) { mpf_init_set_str(x, (char*)str, 10); }
  ~bigfloat(void) { mpf_clear(x); }
  operator const double (void) const { return (double)mpf_get_d(x); }
  static void setDefaultPrecision(unsigned long dprec) {
    unsigned long bprec =  (unsigned long)(3.321928094 * (double)dprec);
    mpf_set_default_prec(bprec);
  }

  void setPrecision(unsigned long dprec) {
    unsigned long bprec =  (unsigned long)(3.321928094 * (double)dprec);
    mpf_set_prec(x,bprec);
  }
  
  unsigned long getPrecision(void) const { return mpf_get_prec(x); }

  unsigned long getDefaultPrecision(void) const { return mpf_get_default_prec(); }

  bigfloat& operator=(const bigfloat& y) {
    mpf_set(x, y.x); 
    return *this;
  }

  bigfloat& operator=(const unsigned long y) { 
    mpf_set_ui(x, y);
    return *this; 
  }
  
  bigfloat& operator=(const signed long y) {
    mpf_set_si(x, y); 
    return *this;
  }
  
  bigfloat& operator=(const float y) {
    mpf_set_d(x, (double)y); 
    return *this;
  }

  bigfloat& operator=(const double y) {
    mpf_set_d(x, y); 
    return *this;
  }

  size_t write(void);
  size_t read(void);

  /* Arithmetic Functions */

  bigfloat& operator+=(const bigfloat& y) { return *this = *this + y; }
  bigfloat& operator-=(const bigfloat& y) { return *this = *this - y; }
  bigfloat& operator*=(const bigfloat& y) { return *this = *this * y; }
  bigfloat& operator/=(const bigfloat& y) { return *this = *this / y; }

  friend bigfloat operator+(const bigfloat& x, const bigfloat& y) {
    bigfloat a;
    mpf_add(a.x,x.x,y.x);
    return a;
  }

  friend bigfloat operator+(const bigfloat& x, const unsigned long y) {
    bigfloat a;
    mpf_add_ui(a.x,x.x,y);
    return a;
  }

  friend bigfloat operator-(const bigfloat& x, const bigfloat& y) {
    bigfloat a;
    mpf_sub(a.x,x.x,y.x);
    return a;
  }
  
  friend bigfloat operator-(const unsigned long x, const bigfloat& y) {
    bigfloat a;
    mpf_ui_sub(a.x,x,y.x);
    return a;
  }
  
  friend bigfloat operator-(const bigfloat& x, const unsigned long y) {
    bigfloat a;
    mpf_sub_ui(a.x,x.x,y);
    return a;
  }

  friend bigfloat operator-(const bigfloat& x) {
    bigfloat a;
    mpf_neg(a.x,x.x);
    return a;
  }

  friend bigfloat operator*(const bigfloat& x, const bigfloat& y) {
    bigfloat a;
    mpf_mul(a.x,x.x,y.x);
    return a;
  }

  friend bigfloat operator*(const bigfloat& x, const unsigned long y) {
    bigfloat a;
    mpf_mul_ui(a.x,x.x,y);
    return a;
  }

  friend bigfloat operator/(const bigfloat& x, const bigfloat& y){
    bigfloat a;
    mpf_div(a.x,x.x,y.x);
    return a;
  }

  friend bigfloat operator/(const unsigned long x, const bigfloat& y){
    bigfloat a;
    mpf_ui_div(a.x,x,y.x);
    return a;
  }

  friend bigfloat operator/(const bigfloat& x, const unsigned long y){
    bigfloat a;
    mpf_div_ui(a.x,x.x,y);
    return a;
  }

  friend bigfloat sqrt_bf(const bigfloat& x){
    bigfloat a;
    mpf_sqrt(a.x,x.x);
    return a;
  }

  friend bigfloat sqrt_bf(const unsigned long x){
    bigfloat a;
    mpf_sqrt_ui(a.x,x);
    return a;
  }

  friend bigfloat abs_bf(const bigfloat& x){
    bigfloat a;
    mpf_abs(a.x,x.x);
    return a;
  }

  friend bigfloat pow_bf(const bigfloat& a, long power) {
    bigfloat b;
    mpf_pow_ui(b.x,a.x,power);
    return b;
  }

  friend bigfloat pow_bf(const bigfloat& a, bigfloat &power) {
    bigfloat b;
    mpfr_pow(b.x,a.x,power.x,GMP_RNDN);
    return b;
  }

  friend bigfloat exp_bf(const bigfloat& a) {
    bigfloat b;
    mpfr_exp(b.x,a.x,GMP_RNDN);
    return b;
  }

  /* Comparison Functions */

  friend int operator>(const bigfloat& x, const bigfloat& y) {
    int test;
    test = mpf_cmp(x.x,y.x);
    if (test > 0) return 1;
    else return 0;
  }

  friend int operator<(const bigfloat& x, const bigfloat& y) {
    int test;
    test = mpf_cmp(x.x,y.x);
    if (test < 0) return 1;
    else return 0;
  }

  friend int sgn(const bigfloat&);

  /* Miscellaneous Functions */

  //  friend bigfloat& random(void);
};

/*
class complex_bf {

 private:
  bigfloat x, y;

 public:
  complex_bf() {
    x = 0l;
    y = 0l;
  }

  complex_bf(bigfloat a, bigfloat b) {
    x = a;
    y = b;
  }

  complex_bf(double a) {
    x = a;
    y = 0l;
  }

  complex_bf& operator=(const complex_bf& a) {
    x = a.x;
    y = a.y;
    return *this;
  }

  complex_bf& operator=(const double a) {
    x = a;
    y = 0l;
    return *this;
  }

  complex_bf& operator=(const long a) {
    x = a;
    y = 0l;
    return *this;
  }

  complex_bf conj(const complex_bf &a) {
    complex_bf b(a.x, -a.y);
    return b;
  }

  friend complex_bf operator+(const complex_bf& a, const complex_bf& b) {
    complex_bf c(a.x+b.x, a.y+b.y);
    return c;
  }

  friend complex_bf operator-(const complex_bf& a, const complex_bf& b) {
    complex_bf c(a.x-b.x, a.y-b.y);
    return c;
  }

  friend complex_bf operator*(const complex_bf& a, const complex_bf& b) {
    complex_bf c(a.x*b.x-a.y*b.y, a.x*b.y+a.y*b.x);
    return c;
  }

  friend complex_bf operator*(const long a, const complex_bf& b) {
    complex_bf c(b.x*(bigfloat)a, b.y*(bigfloat)a);
    return c;
  }

  friend complex_bf operator*(const bigfloat &a, const complex_bf& b) {
    complex_bf c(a*b.x, a*b.y);
    return c;
  }

  friend complex_bf operator*(const complex_bf &a, const bigfloat& b) {
    return operator*(b, a);
  }

  friend complex_bf operator/(const complex_bf& a, const complex_bf& b) {
    bigfloat denominator = abs_complex(b);
    complex_bf numerator = a*conj(b);
    return numerator / denominator;
  }

  friend complex_bf operator/(const complex_bf& a, const bigfloat& b) {
    return new complex_bf(a.x/b, a.y/b);
  }

  complex_bf& operator+=(const complex_bf& y) { return *this = *this + y; }
  complex_bf& operator-=(const complex_bf& y) { return *this = *this - y; }
  complex_bf& operator*=(const complex_bf& y) { return *this = *this * y; }
  complex_bf& operator/=(const complex_bf& y) { return *this = *this / y; }

  friend bigfloat abs_complex(const complex_bf& a){
    return a.x*a.x+a.y*a.y;
  }

};
*/
#endif
