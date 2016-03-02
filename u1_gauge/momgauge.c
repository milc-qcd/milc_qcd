/* ************************************************************	*/
/*								*/
/* 			       MOMGAUGE.C	   		*/
/*								*/
/* Generate U(1) fields A_\mu(p) in momentum space: 		*/
/*  S= 1/2 \sum(k) [A0(k)^2 \vec{k'}^2 + \vec{A}^2 k'^2]	*/
/*     k'^2 = \sum [4 (sin k/2)^2] 				*/
/*  \eta(k)_\mu = gaussrand()					*/
/*  A0(k) = [\vec{k'}^2/2]^{-1/2} * \eta{k}			*/
/*  \vec{A}(k) = [k'^2/2]^{-1/2} * \eta{k}			*/
/*								*/
/* 6/5/12 CD Return the FT of the Re A(x) only                  */
/*								*/
/*								*/
/* ************************************************************	*/

#include "include_u1g.h"

void momgauge(complex *u1gf)
{

  complex Am,Ap;
  Real lk[4],mm,mp,r1,r2;
  Real lkssq,lktsq,lksq;
  register site *s;
  register int i,dir;
  int hp;

  /* initialize */
  FORALLSITES(i,s){
      FORALLUPDIR(dir){
	  u1gf[4*i+dir]=cmplx(0.0,0.0);
	 }
    }

  /* A0 */
  FORALLSITES(i,s){
    Real mom[4] = { 2.*PI*s->x/nx, 2.*PI*s->y/ny, 2.*PI*s->z/nz, 2.*PI*s->t/nt };
    
    if(latin[i] == i)hp = 1;
    else hp = 2;


      /* lattice mom: \vec{k}^2, k1,k2,k3 and squares */
      lkssq=0.0;
      for(dir=XUP;dir<=ZUP;dir++){
          lk[dir]=2.0*(Real)sin((double)(mom[dir])/2.0);
	  lkssq+=sqr(lk[dir]);
	 }
      if(lkssq==0.0)				/* \vec{k}==0 */
	{
	if(latin[i]!=junk_id)
	  {
	  u1gf[4*i+TUP]=cmplx(0.0,0.0);
	  }
	} 
      else					/* \vec{k}!=0 */
	{
	if(latin[i]!=junk_id)
	  {
	   u1gf[4*i+TUP] = complex_gaussian_rand_no(&(s->site_prn));
	   CDIVREAL( u1gf[4*i+TUP], sqrt(hp*lkssq/2.0), u1gf[4*i+TUP] );
	  }
	}
     } /* FORALLSITES-ends */

  /* \vec{A} */
  FORALLSITES(i,s){
    Real mom[4] = { 2.*PI*s->x/nx, 2.*PI*s->y/ny, 2.*PI*s->z/nz, 2.*PI*s->t/nt };

    if(latin[i] == i)hp = 1;
    else hp = 2;
      /* lattice mom: k^2, k0,k1,k2,k3 and squares */
      lkssq=lktsq=lksq=0.0;
      for(dir=XUP;dir<=TUP;dir++){
	  lk[dir]=2.0*(Real)sin((double)(mom[dir])/2.0);
	  if(dir!=TUP) lkssq+=sqr(lk[dir]);
	  if(dir==TUP) lktsq+=sqr(lk[dir]);
	 }
      lksq=lkssq+lktsq;

      if(lkssq==0.0 && lktsq==0.0)		/* \vec{k}=0 & k0=0 */
	{
	if(latin[i]!=junk_id)
	  {
	  for(dir=XUP;dir<=ZUP;dir++){
	      u1gf[4*i+dir]=cmplx(0.0,0.0);
	     }
	  }
	}
      else if(lkssq==0.0 && lktsq!=0.0)		/* \vec{k}=0 & k0!=0 */
	{
	if(latin[i]!=junk_id)
	  {
	  for(dir=XUP;dir<=ZUP;dir++){
	      u1gf[4*i+dir] = complex_gaussian_rand_no(&(s->site_prn));
	      CDIVREAL( u1gf[4*i+dir], sqrt(hp*lksq/2.0), u1gf[4*i+dir] );
	     }
	  }
	}
      else if(lkssq!=0.0) 			/* \vec{k}!=0 & k0=any */
	{
	if(lk[YUP]==0.0 && lk[ZUP]==0.0)	/* k1!=0 & k2=k3=0 */
	  {
	  if(latin[i]!=junk_id)
	    {
	    for(dir=YUP;dir<=ZUP;dir++){
		u1gf[4*i+dir] = complex_gaussian_rand_no(&(s->site_prn));
		CDIVREAL( u1gf[4*i+dir], sqrt(hp*lksq/2.0), u1gf[4*i+dir] );
               }
            u1gf[4*i+XUP].real=-(lk[YUP]*u1gf[4*i+YUP].real+
                                lk[ZUP]*u1gf[4*i+ZUP].real)/lk[XUP];
            u1gf[4*i+XUP].imag=-(lk[YUP]*u1gf[4*i+YUP].imag+
                                lk[ZUP]*u1gf[4*i+ZUP].imag)/lk[XUP];
	    }
	  }
	else if(lk[ZUP]==0.0 && lk[XUP]==0.0)	/* k2!=0 & k3=k1=0 */
	  {
	  if(latin[i]!=junk_id)
	    {
	    for(dir=XUP;dir<=ZUP;dir++)if(dir!=YUP){
		u1gf[4*i+dir] = complex_gaussian_rand_no(&(s->site_prn));
		CDIVREAL( u1gf[4*i+dir], sqrt(hp*lksq/2.0), u1gf[4*i+dir] );
               }
            u1gf[4*i+YUP].real=-(lk[ZUP]*u1gf[4*i+ZUP].real+
                                lk[XUP]*u1gf[4*i+XUP].real)/lk[YUP];
            u1gf[4*i+YUP].imag=-(lk[ZUP]*u1gf[4*i+ZUP].imag+
                                lk[XUP]*u1gf[4*i+XUP].imag)/lk[YUP];
	    }
	  }
	else if(lk[XUP]==0.0 && lk[YUP]==0.0)	/* k3!=0 & k1=k2=0 */
	  {
	  if(latin[i]!=junk_id)
	    {
	    for(dir=XUP;dir<=YUP;dir++){
		u1gf[4*i+dir] = complex_gaussian_rand_no(&(s->site_prn));
		CDIVREAL( u1gf[4*i+dir], sqrt(hp*lksq/2.0), u1gf[4*i+dir] );
               }
            u1gf[4*i+ZUP].real=-(lk[XUP]*u1gf[4*i+XUP].real+
                                lk[YUP]*u1gf[4*i+YUP].real)/lk[ZUP];
            u1gf[4*i+ZUP].imag=-(lk[XUP]*u1gf[4*i+XUP].imag+
                                lk[YUP]*u1gf[4*i+YUP].imag)/lk[ZUP];
	    }
	  }
	else if(lk[XUP]!=0.0 && lk[YUP]!=0.0 && lk[ZUP]==0.0)
	  {					/* k1,k2!=0, k3=0 */
	  if(latin[i]!=junk_id)
	    {
	    mm=lksq;
            mp=lksq*(1+((sqr(lk[XUP])+sqr(lk[ZUP]))/sqr(lk[YUP])));
            r1=lk[XUP]/((Real)sqrt(sqr(lk[XUP])+sqr(lk[ZUP])));
            r2=lk[ZUP]/((Real)sqrt(sqr(lk[XUP])+sqr(lk[ZUP])));
	    Am = complex_gaussian_rand_no(&(s->site_prn));
	    CDIVREAL( Am, sqrt(hp*mm/2.0), Am );
	    Ap = complex_gaussian_rand_no(&(s->site_prn));
	    CDIVREAL( Ap, sqrt(hp*mp/2.0), Ap );

            u1gf[4*i+XUP].real=r2*Am.real+r1*Ap.real;
            u1gf[4*i+XUP].imag=r2*Am.imag+r1*Ap.imag;
            u1gf[4*i+ZUP].real=-r1*Am.real+r2*Ap.real;
            u1gf[4*i+ZUP].imag=-r1*Am.imag+r2*Ap.imag;
            u1gf[4*i+YUP].real=-(lk[XUP]*u1gf[4*i+XUP].real+
                                lk[ZUP]*u1gf[4*i+ZUP].real)/lk[YUP];
            u1gf[4*i+YUP].imag=-(lk[XUP]*u1gf[4*i+XUP].imag+
                                lk[ZUP]*u1gf[4*i+ZUP].imag)/lk[YUP];
	    for(dir=XUP;dir<=ZUP;dir++){
	    }
	    }
	  }
	else					/* k3!=0 & k1 &/or k2!=0 */
	  {
	  if(latin[i]!=junk_id)
	    {
	    mm=lksq;
            mp=lksq*(1+((sqr(lk[XUP])+sqr(lk[YUP]))/sqr(lk[ZUP])));
            r1=lk[XUP]/((Real)sqrt(sqr(lk[XUP])+sqr(lk[YUP])));
            r2=lk[YUP]/((Real)sqrt(sqr(lk[XUP])+sqr(lk[YUP])));
	    Am = complex_gaussian_rand_no(&(s->site_prn));
	    CDIVREAL( Am, sqrt(hp*mm/2.0), Am );
	    Ap = complex_gaussian_rand_no(&(s->site_prn));
	    CDIVREAL( Ap, sqrt(hp*mp/2.0), Ap );

            u1gf[4*i+XUP].real=r2*Am.real+r1*Ap.real;
            u1gf[4*i+XUP].imag=r2*Am.imag+r1*Ap.imag;
            u1gf[4*i+YUP].real=-r1*Am.real+r2*Ap.real;
            u1gf[4*i+YUP].imag=-r1*Am.imag+r2*Ap.imag;
            u1gf[4*i+ZUP].real=-(lk[XUP]*u1gf[4*i+XUP].real+
                                lk[YUP]*u1gf[4*i+YUP].real)/lk[ZUP];
            u1gf[4*i+ZUP].imag=-(lk[XUP]*u1gf[4*i+XUP].imag+
                                lk[YUP]*u1gf[4*i+YUP].imag)/lk[ZUP];
	    for(dir=XUP;dir<=ZUP;dir++){
	    }
	    }
	  }
	} /* lkssq!=0 ends */

     } /* FORALLSITES-ends */

  /* Adjust u1gf so we return FT of Re A(k) */
  FORALLSITES(i,s){
    if(i==latin[i]){
      for(dir=XUP;dir<=TUP;dir++){ 
	u1gf[4*i+dir].imag=0.0;
      }
    } else {
      for(dir=XUP;dir<=TUP;dir++){ 
	CMULREAL(u1gf[4*i+dir],2.0,u1gf[4*i+dir]);
      }
    }
  }

} /* end of momgauge() */

/* ************************************************************	*/

