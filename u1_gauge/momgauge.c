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

#define MILC_ORIGINAL 1
#define HATTON_REVISION 2

#define FEYNMAN 1
#define LORENZ 2
#define COULOMB 3

#define QED_L 1
#define QED_TL 2

void show_momgauge_opts(void){
#if GAUGE_METHOD == MILC_ORIGINAL
  node0_printf("MILC original: Coulomb gauge QED_TL scheme\n");
#elif GAUGE_METHOD == HATTON_REVISION
  node0_printf("Hatton proposal: Feynman gauge QED_L scheme\n");
#else
  #if GAUGE == COULOMB
     #if SCHEME == QED_L
  node0_printf("Coulomb gauge QED_L scheme\n");
     #elif SCHEME == QED_TL
  node0_printf("Coulomb gauge QED_TL scheme\n");
     #else
  node0_printf("Coulomb gauge - no zero mode treatment\n");
     #endif
  #elif GAUGE == LORENZ
     #if SCHEME == QED_L
  node0_printf("Lorenz gauge QED_L scheme\n");
     #elif SCHEME == QED_TL
  node0_printf("Lorenz gauge QED_L scheme\n");
     #else
  node0_printf("Lorenz gauge - no zero mode treatment\n");
     #endif
  #else
     #if SCHEME == QED_L
  node0_printf("Feynman gauge QED_L scheme\n");
     #elif SCHEME == QED_TL
  node0_printf("Feynman gauge QED_TL scheme\n");
     #else
  node0_printf("Feynman gauge - no zero mode treatment\n");
     #endif
  #endif
#endif  
}

#if GAUGE_METHOD == MILC_ORIGINAL

/* MILC original version */
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

      /* Insert e^{iak/2} factors */
      complex phase[4] = {ce_itheta(mom[0]/2.),ce_itheta(mom[1]/2.),ce_itheta(mom[2]/2.),ce_itheta(mom[3]/2.)};
      int dir;
      FORALLUPDIR(dir){
	complex cc;
	CMUL(phase[dir],u1gf[4*i+dir], cc);
	u1gf[4*i+dir] = cc;
      }
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

#elif GAUGE_METHOD == HATTON_REVISION

/* Dan Hatton's revision */
void momgauge(complex *u1gf)
{
  
  complex Am,Ap;
  Real lk[4],mm,mp,r1,r2,newlk[4];
  Real lkssq,lktsq,lksq,newlksq;
  register site *s;
  register int i,dir;
  int hp;
  
  /* initialize */
  FORALLSITES(i,s){
    FORALLUPDIR(dir){
      u1gf[4*i+dir]=cmplx(0.0,0.0);
    }
  }
  
  
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
	u1gf[4*i+TUP]=cmplx(0.0,0.0);
	u1gf[4*i+XUP]=cmplx(0.0,0.0);
	u1gf[4*i+YUP]=cmplx(0.0,0.0);
	u1gf[4*i+ZUP]=cmplx(0.0,0.0);
      }
    
    else if(mom[0]==PI && mom[1]==PI && mom[2]==PI && mom[3]==PI)
      {
	
	lkssq=lktsq=lksq=0.0;
	for(dir=XUP;dir<=TUP;dir++){
	  lk[dir]=2.0*(Real)sin((double)(mom[dir])/2.0);
	  if(dir!=TUP) lkssq+=sqr(lk[dir]);
	  if(dir==TUP) lktsq+=sqr(lk[dir]);
	}
	lksq=lkssq+lktsq;
	
        u1gf[4*i+TUP].real=
	  (Real)(gaussian_rand_no(&(s->site_prn))*sqrt(2.0/lksq));
	u1gf[4*i+XUP].real=
	  (Real)(gaussian_rand_no(&(s->site_prn))*sqrt(2.0/lksq));
	u1gf[4*i+YUP].real=
	  (Real)(gaussian_rand_no(&(s->site_prn))*sqrt(2.0/lksq));
	u1gf[4*i+ZUP].real=
	  (Real)(gaussian_rand_no(&(s->site_prn))*sqrt(2.0/lksq));
      }
    
    else	
      {
	
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
	
	
	u1gf[4*i+TUP].real=
	  (Real)(gaussian_rand_no(&(s->site_prn))*sqrt(1.0/lksq));
	u1gf[4*i+TUP].imag=
	  (Real)(gaussian_rand_no(&(s->site_prn))*sqrt(1.0/lksq));
	u1gf[4*i+XUP].real=
	  (Real)(gaussian_rand_no(&(s->site_prn))*sqrt(1.0/lksq));
	u1gf[4*i+XUP].imag=
	  (Real)(gaussian_rand_no(&(s->site_prn))*sqrt(1.0/lksq));
	u1gf[4*i+YUP].real=
	  (Real)(gaussian_rand_no(&(s->site_prn))*sqrt(1.0/lksq));
	u1gf[4*i+YUP].imag=
	  (Real)(gaussian_rand_no(&(s->site_prn))*sqrt(1.0/lksq));
	u1gf[4*i+ZUP].real=
	  (Real)(gaussian_rand_no(&(s->site_prn))*sqrt(1.0/lksq));
	u1gf[4*i+ZUP].imag=
	  (Real)(gaussian_rand_no(&(s->site_prn))*sqrt(1.0/lksq));
	
	
	
      }
    
    if(mom[0]==PI || mom[0]==0)
      {
	if(mom[1]==PI || mom[1]==0)
	  {
	    if(mom[2]==PI || mom[2]==0)
	      {
		if(mom[3]==PI || mom[3]==0)
		  {
		    u1gf[4*i+TUP].real=
		      (Real)(gaussian_rand_no(&(s->site_prn))*sqrt(2.0/lksq));
		    u1gf[4*i+TUP].imag=0.0;
		    u1gf[4*i+XUP].real=
		      (Real)(gaussian_rand_no(&(s->site_prn))*sqrt(2.0/lksq));
		    u1gf[4*i+XUP].imag=0.0;
		    u1gf[4*i+YUP].real=
		      (Real)(gaussian_rand_no(&(s->site_prn))*sqrt(2.0/lksq));
		    u1gf[4*i+YUP].imag=0.0;
		    u1gf[4*i+ZUP].real=
		      (Real)(gaussian_rand_no(&(s->site_prn))*sqrt(2.0/lksq));
		    u1gf[4*i+ZUP].imag=0.0;
		  }
	      }
	  }
      }
    
    
    if(lkssq==0.0)        /* make sure QED_L didn't get messed up */
      {
	u1gf[4*i+TUP]=cmplx(0.0,0.0);
	u1gf[4*i+XUP]=cmplx(0.0,0.0);
	u1gf[4*i+YUP]=cmplx(0.0,0.0);
	u1gf[4*i+ZUP]=cmplx(0.0,0.0);
      } 
    
  } /* FORALLSITES-ends */
  
  /* Extra e^{iak/2} factor */
  
  FORALLSITES(i,s){
    Real mom[4] = { 2.*PI*s->x/nx, 2.*PI*s->y/ny, 2.*PI*s->z/nz, 2.*PI*s->t/nt };
    complex phase[4] = {ce_itheta(mom[0]/2.),ce_itheta(mom[1]/2.),ce_itheta(mom[2]/2.),ce_itheta(mom[3]/2.)};
    complex cc;
    int dir;
    FORALLUPDIR(dir){
      CMUL(phase[dir],u1gf[4*i+dir], cc);
      u1gf[4*i+dir] = cc;
    }
//    u1gf[4*i+TUP] = cmul(&phase[3],&u1gf[4*i+TUP]);
//    u1gf[4*i+XUP] = cmul(&phase[0],&u1gf[4*i+XUP]);
//    u1gf[4*i+YUP] = cmul(&phase[1],&u1gf[4*i+YUP]);
//    u1gf[4*i+ZUP] = cmul(&phase[2],&u1gf[4*i+ZUP]);
  }
  
  /* Arrange components so that Fourier transform is real */
  
  FORALLSITES(i,s){
    int newx,newy,newz,newt;
    newx = ((nx)-s->x)%nx;
    newy = ((ny)-s->y)%ny;
    newz = ((nz)-s->z)%nz;
    newt = ((nt)-s->t)%nt;
    
    Real mom[4] = { 2.*PI*s->x/nx, 2.*PI*s->y/ny, 2.*PI*s->z/nz, 2.*PI*s->t/nt };
    Real newmom[4] = { 2.*PI*newx/nx, 2.*PI*newy/ny, 2.*PI*newz/nz, 2.*PI*newt/nt };
    
    if(s->x==0){newx=0;}
    if(s->x==nx/2){newx=nx/2;}
    if(s->y==0){newy=0;}
    if(s->y==ny/2){newy=ny/2;}
    if(s->z==0){newz=0;}
    if(s->z==nz/2){newz=nz/2;}
    if(s->t==0){newt=0;}
    if(s->t==nt/2){newt=nt/2;}
    
    lksq=0.0;
    newlksq=0.0;
    for(dir=XUP;dir<=TUP;dir++){
      lk[dir]=2.0*(Real)sin((double)(mom[dir])/2.0);
      lksq+=lk[dir]*lk[dir];
      newlk[dir]=2.0*(Real)sin((double)(newmom[dir])/2.0);
      newlksq+=newlk[dir]*newlk[dir];
    }
    
    u1gf[4*i+TUP].real = u1gf[4*node_index(newx,newy,newz,newt )+TUP].real;
    u1gf[4*i+TUP].imag = -u1gf[4*node_index(newx,newy,newz,newt )+TUP].imag;
    u1gf[4*i+XUP].real = u1gf[4*node_index(newx,newy,newz,newt )+XUP].real;
    u1gf[4*i+XUP].imag = -u1gf[4*node_index(newx,newy,newz,newt )+XUP].imag;
    u1gf[4*i+YUP].real = u1gf[4*node_index(newx,newy,newz,newt )+YUP].real;
    u1gf[4*i+YUP].imag = -u1gf[4*node_index(newx,newy,newz,newt )+YUP].imag;
    u1gf[4*i+ZUP].real = u1gf[4*node_index(newx,newy,newz,newt )+ZUP].real;
    u1gf[4*i+ZUP].imag = -u1gf[4*node_index(newx,newy,newz,newt )+ZUP].imag;
    
  }

#if 0
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
#endif

} /* end of momgauge() */


#else

/* ************************************************************	*/

/* Transform to Coulomb gauge */

/* The transformation is

     A_mu(k) = A_mu(k) - c \hat k_mu

   where

     c = vec(k) vec(A)/ vec(k)^2

   Note that we also need to deal with the case \vec(k) = 0!
   See QED_L or QED_TL below

*/
static void u1_Coulomb(complex *u1gfsite, Real *lk, Real lkssq){
  complex cdot = cmplx(0.0,0.0);
  complex cc;

  if(lkssq == 0.0)return;
  int dir;
  FORALLUPDIRBUT(TUP,dir){
    CMULREAL( u1gfsite[dir], lk[dir]/lkssq, cc );
    CSUM( cdot, cc );
  }
  FORALLUPDIR(dir){
    CMULREAL( cdot, lk[dir], cc );
    CSUB( u1gfsite[dir], cc, u1gfsite[dir] );
  }
}

static void u1_Lorenz(complex *u1gfsite, Real *lk, Real lksq){
  /* Project out the longitudinal mode */
  complex cdot = cmplx(0.0,0.0);
  complex cc;

  if(lksq == 0.0)return;
  int dir;
  FORALLUPDIR(dir){
    CMULREAL( u1gfsite[dir], lk[dir]/lksq, cc );
    CSUM( cdot, cc );
  }
  FORALLUPDIR(dir){
    CMULREAL( cdot, lk[dir], cc );
    CSUB( u1gfsite[dir], cc, u1gfsite[dir] );
  }
}


/* Start with Feynman gauge and then fix to the desired gauge.  Note
   that with Lorenz-gauge fixing we get the same ensemble of vector
   potentials as with the quenched non-compact U(1) action.
 */

/* The Lorenz gauge action in momentum space is

   S = A_mu(k) (k^2 delta_mu,nu - k_mu k_nu) A^*_nu(k)
     = \sum_i (A_mu eps^i_mu)^2 k^2
  
   where eps^i_mu is a normalized polarization vector orthogonal to k_mu.

   We drop the k = 0 mode as well as the longitudinal mode A_mu ~ k_mu.

   Then we generate A_mu by projecting out the gauge mode from a Gaussian random
   complex four-vector \eta_nu
    
       A_mu(k) = (delta_mu,nu - k_mu k_nu/k^2) \eta_nu/\sqrt(k^2)

   where 
 
      <\eta_nu \eta^_mu> = delta_mu,nu ,

   Without the projection, the vector potential is in "Feynman gauge".
      
*/

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

  FORALLSITES(i,s){
    Real mom[4] = { 2.*PI*s->x/nx, 2.*PI*s->y/ny, 2.*PI*s->z/nz, 2.*PI*s->t/nt };
    
    if(latin[i] == i)hp = 1;
    else hp = 2;
    
    /* lattice mom: kx,ky,kz,kt and squares */
    lkssq=0.0; lksq = 0.0;
    for(dir=XUP;dir<=TUP;dir++){
      lk[dir]=2.0*(Real)sin((double)(mom[dir])/2.0);
      lksq+=sqr(lk[dir]);
      if(dir != TUP)lkssq += sqr(lk[dir]);
    }

    if(lksq == 0.0){
      /* No constant mode */
      FORALLUPDIR(dir){
	u1gf[4*i+dir]=cmplx(0.0,0.0);
      }
      continue;
    }

    /* Skip the repeated components */
    if(latin[i]==junk_id)continue;
    
    /* Start with a random complex four-vector */
    FORALLUPDIR(dir){
      u1gf[4*i+dir] = complex_gaussian_rand_no(&(s->site_prn));
      CDIVREAL( u1gf[4*i+dir], sqrt(hp*lksq/2.0), u1gf[4*i+dir] );
    }
    
    /* Defaults to Feynman gauge */
    
#if GAUGE == COULOMB
    /* Fix Coulomb gauge */
    u1_Coulomb( u1gf + 4*i, lk, lkssq );
    
#elif GAUGE == LORENZ      
    /* Fix Lorenz gauge */
    u1_Lorenz( u1gf + 4*i, lk, lksq );
#endif
    
#if SCHEME == QED_L
    /* QED_L condition */
    if(lkssq == 0.0){
      /* All components are zero */
      int dir;
      FORALLUPDIR(dir){
	u1gf[4*i+dir] = cmplx(0.0, 0.0);
      }
    }
    
#elif SCHEME == QED_TL
    /* QED_TL condition */
    if(lkssq == 0.0){
      /* Time component is zero regargless of lk[TUP] */
      u1gf[TUP] = cmplx(0.0, 0.0);
      int dir;
      /* Spatial components are zero only if lk[TUP] is zero */
      if(lk[TUP] == 0.0){
	FORALLUPDIRBUT(TUP,dir){
	  u1gf[4*i+dir] = cmplx(0.0, 0.0);
	}
      }
    }
#endif    
	
    /* Insert e^{iak/2} factors */
    complex phase[4] = {ce_itheta(mom[0]/2.),ce_itheta(mom[1]/2.),ce_itheta(mom[2]/2.),ce_itheta(mom[3]/2.)};
    int dir;
    FORALLUPDIR(dir){
      complex cc;
      CMUL(phase[dir],u1gf[4*i+dir], cc);
      u1gf[4*i+dir] = cc;
    }
  } /* FORALLSITES-ends */
  
  /* Adjust u1gf so we return FT of Re A(k) */
  FORALLSITES(i,s){
    if(i==latin[i]){
      int dir;
      FORALLUPDIR(dir){
	u1gf[4*i+dir].imag=0.0;
      }
    } else {
      int dir;
      FORALLUPDIR(dir){
	CMULREAL(u1gf[4*i+dir],2.0,u1gf[4*i+dir]);
      }
    }
  }
} /* end of momgauge() */

#endif
