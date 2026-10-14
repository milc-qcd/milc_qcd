/******************** u1_gauge_field.c **********************************/
/* MIMD version 7 */

/* Self-contained generation of a stochastic U(1) vector potential
   A_mu(x).

   In quantum electrodynamics the momentum-space distribution of the
   photon field A_mu is Gaussian, so the field can be generated
   independently at each momentum (no Markov chain is needed) and then
   transformed to coordinate space with a fast Fourier transform.

   Usage:

     initialize_u1_gauge_field(U1_COULOMB, U1_QED_L);   // once
     Real *A = create_u1_A_field();                     // caller owns A
     generate_u1_gauge_field(A);                        // as often as wanted
     ...
     destroy_u1_A_field(A);
     finalize_u1_gauge_field();

   The caller allocates and frees the A field.  It must have room for
   4*sites_on_node Reals, indexed as A[4*i+dir], which is what
   create_u1_A_field() returns.

   The gauge choice and the treatment of the vec(k) = 0 modes are
   arguments to initialize_u1_gauge_field(), so no compilation macros
   are needed to call this routine.  The choice of algorithm
   (GAUGE_METHOD) is still a compile-time macro; when it is left
   undefined the general algorithm is used, and only that algorithm
   honors the gauge and scheme arguments.

   This code supports only scalar running.

   Originally u1_gauge/momgauge.c and the set_mom_table() routine in
   u1_gauge/setup.c.  Collected here, with the former globals latin and
   junk_id made file-static, so that any application can generate a
   U(1) field without defining them.

   6/5/12 CD Return the FT of the Re A(x) only
*/

#include "generic_u1_includes.h"
#include "../include/generic_u1.h"
#include "../include/random.h"

/* The momentum-space generation visits only half of the momenta and
   builds the other half by reflection, so that the Fourier transform
   of the result is real.  mom_status[i], indexed by the LOCAL site
   index i, says which role the momentum carried by that site plays.

   This replaces the old "latin" table, which stored, at the
   whole-lattice (single-process) index of the site, the index of its
   parity-inverted partner, and junk_id for the reflected copies.  Two
   distinct indices were conflated there: the whole-lattice index, used
   to decide which member of the +-k pair is generated, and the local
   index, used to address the table.  They agree only on one process,
   which is why the old table was wrong under MPP.  Keeping the roles
   in an explicit code lets the two indices be used for their separate
   purposes. */

#define U1_MOM_SKIP  0   /* reflected copy; not generated independently */
#define U1_MOM_PAIR  1   /* representative of a +k,-k pair             */
#define U1_MOM_SELF  2   /* self-conjugate momentum, k = -k            */

static char *mom_status = NULL;

/* Run-time options, set by initialize_u1_gauge_field() */
static u1_gauge_fix_t gauge_fix = U1_FEYNMAN;
static u1_zero_mode_t zero_mode_scheme = U1_QED_NONE;
static int is_initialized = 0;

static Real sqr(Real val)
{
  return(val*val);
}

/*------------------------------------------------------------------*/
/* Added to support MPP gauge field generation that gives the same
   random field as the original single-process field generation
   code.  This is the single-process node_index function.

   It must use the whole lattice volume, not sites_on_node: the two are
   equal only when there is a single process, and the value is used to
   order the whole lattice.  It is used only to choose, consistently
   across processes, which member of each +-k pair is generated; the
   table itself is addressed by the local node_index(). */

static size_t scalar_node_index(const int coords[], int dim, int size[])
{
  int d;
  size_t rank = coords[dim-1];
  int sum = coords[dim-1];

  for(d = dim-2; d >= 0; d--){
    rank = rank * size[d] + coords[d];
    sum += coords[d];
  }

  if(sum % 2 == 0)return rank/2;
  else return (rank + volume)/2;

}

/* ************************************************************	*/

/* setup momentum table */
static void set_mom_table(void)
{

  int x,y,z,t;
  int i,dir;
  size_t j,k;
  int dimn[4] = {nx, ny, nz, nt};
  site *s;

  /* initialization */

  mom_status=(char *)malloc(sites_on_node*sizeof(char));
  if(mom_status == NULL){
    node0_printf("set_mom_table: No room for the momentum table\n");
    terminate(1);
  }

  FORALLSITES(i,s){
    mom_status[i] = U1_MOM_SKIP;
  }

  /* This algorithm is needed to achieve backward compatibility with
     the single-process code */

  for(x=0;x<nx;x++)for(y=0;y<ny;y++)for(z=0;z<nz;z++)for(t=0;t<nt;t++){
	  int st[4] = { x, y, z, t };
	  int stin[4];
	  
	  /* The parity-inverted coordinate */
	  FORALLUPDIR(dir){
	    stin[dir] = (dimn[dir] - st[dir]) % dimn[dir];
	  }
	    
	  /* Of the pair (st, stin), the one with the smaller
	     whole-lattice index is the one we generate.  The Gaussian
	     numbers are drawn from that site's own site_prn, which is
	     seeded from the layout-independent lattice[i].index, so
	     making this choice the same on any number of processes is
	     what reproduces the single-process field. */

	  j = scalar_node_index(stin, 4, dimn);
	  k = scalar_node_index(st, 4, dimn);
	  if(j <= k){
	    /* The site we generate is stin.  Record its role on the
	       node that owns it, at its local index. */
	    if(node_number(stin[XUP], stin[YUP], stin[ZUP], stin[TUP]) == this_node){
	      mom_status[node_index(stin[XUP], stin[YUP], stin[ZUP], stin[TUP])]
		= (j == k) ? U1_MOM_SELF : U1_MOM_PAIR;
	    }
	  }
	}

} /* end of set_mom_table() */

/* ************************************************************	*/

void show_u1_gauge_field_opts(u1_gauge_fix_t gauge, u1_zero_mode_t scheme){
#if GAUGE_METHOD == U1_MILC_ORIGINAL
  node0_printf("MILC original: Coulomb gauge QED_TL scheme\n");
#elif GAUGE_METHOD == U1_HATTON_REVISION
  node0_printf("Hatton proposal: Feynman gauge QED_L scheme\n");
#else
  const char *gauge_name;
  const char *scheme_name;

  switch(gauge){
  case U1_COULOMB: gauge_name = "Coulomb gauge"; break;
  case U1_LORENZ:  gauge_name = "Lorenz gauge";  break;
  default:         gauge_name = "Feynman gauge"; break;
  }

  switch(scheme){
  case U1_QED_L:  scheme_name = " QED_L scheme";  break;
  case U1_QED_TL: scheme_name = " QED_TL scheme"; break;
  default:        scheme_name = " - no zero mode treatment"; break;
  }

  node0_printf("%s%s\n", gauge_name, scheme_name);
#endif
}

/* ************************************************************	*/
/* Generate A_mu(k) in momentum space:				*/
/*  S= 1/2 \sum(k) [A0(k)^2 \vec{k'}^2 + \vec{A}^2 k'^2]	*/
/*     k'^2 = \sum [4 (sin k/2)^2] 				*/
/*  \eta(k)_\mu = gaussrand()					*/
/*  A0(k) = [\vec{k'}^2/2]^{-1/2} * \eta{k}			*/
/*  \vec{A}(k) = [k'^2/2]^{-1/2} * \eta{k}			*/
/* ************************************************************	*/

#if GAUGE_METHOD == U1_MILC_ORIGINAL

/* MILC original version */
static void momgauge(complex *u1gf)
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
    
    if(mom_status[i] == U1_MOM_SELF)hp = 1;
    else hp = 2;


      /* lattice mom: \vec{k}^2, k1,k2,k3 and squares */
      lkssq=0.0;
      for(dir=XUP;dir<=ZUP;dir++){
          lk[dir]=2.0*(Real)sin((double)(mom[dir])/2.0);
	  lkssq+=sqr(lk[dir]);
	 }
      if(lkssq==0.0)				/* \vec{k}==0 */
	{
	if(mom_status[i]!=U1_MOM_SKIP)
	  {
	  u1gf[4*i+TUP]=cmplx(0.0,0.0);
	  }
	} 
      else					/* \vec{k}!=0 */
	{
	if(mom_status[i]!=U1_MOM_SKIP)
	  {
	   u1gf[4*i+TUP] = complex_gaussian_rand_no(&(s->site_prn));
	   CDIVREAL( u1gf[4*i+TUP], sqrt(hp*lkssq/2.0), u1gf[4*i+TUP] );
	  }
	}
     } /* FORALLSITES-ends */

  /* \vec{A} */
  FORALLSITES(i,s){
    Real mom[4] = { 2.*PI*s->x/nx, 2.*PI*s->y/ny, 2.*PI*s->z/nz, 2.*PI*s->t/nt };

    if(mom_status[i] == U1_MOM_SELF)hp = 1;
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
	if(mom_status[i]!=U1_MOM_SKIP)
	  {
	  for(dir=XUP;dir<=ZUP;dir++){
	      u1gf[4*i+dir]=cmplx(0.0,0.0);
	     }
	  }
	}
      else if(lkssq==0.0 && lktsq!=0.0)		/* \vec{k}=0 & k0!=0 */
	{
	if(mom_status[i]!=U1_MOM_SKIP)
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
	  if(mom_status[i]!=U1_MOM_SKIP)
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
	  if(mom_status[i]!=U1_MOM_SKIP)
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
	  if(mom_status[i]!=U1_MOM_SKIP)
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
	  if(mom_status[i]!=U1_MOM_SKIP)
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
	    }
	  }
	else					/* k3!=0 & k1 &/or k2!=0 */
	  {
	  if(mom_status[i]!=U1_MOM_SKIP)
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
    if(mom_status[i]==U1_MOM_SELF){
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

#elif GAUGE_METHOD == U1_HATTON_REVISION

/* Dan Hatton's revision */

/* NOTE: this variant is single-process only, independently of the
   momentum table.  It builds the reflected half of the momenta with
   u1gf[4*node_index(newx,newy,newz,newt)+dir], but node_index() is a
   local index and the reflected site generally lives on another
   process; it also reads that element out of the same array it is
   overwriting in place.  Running it on more than one rank gives
   neither the single-process answer nor a correct one. */
static void momgauge(complex *u1gf)
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
    
    if(mom_status[i] == U1_MOM_SELF)hp = 1;
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
	
	if(mom_status[i] == U1_MOM_SELF)hp = 1;
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

static void momgauge(complex *u1gf)
{

  Real lk[4];
  Real lkssq,lksq;
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
    
    if(mom_status[i] == U1_MOM_SELF)hp = 1;
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
    if(mom_status[i]==U1_MOM_SKIP)continue;
    
    /* Start with a random complex four-vector */
    FORALLUPDIR(dir){
      u1gf[4*i+dir] = complex_gaussian_rand_no(&(s->site_prn));
      CDIVREAL( u1gf[4*i+dir], sqrt(hp*lksq/2.0), u1gf[4*i+dir] );
    }
    
    /* Defaults to Feynman gauge */
    
    if(gauge_fix == U1_COULOMB){
      /* Fix Coulomb gauge */
      u1_Coulomb( u1gf + 4*i, lk, lkssq );
      
    } else if(gauge_fix == U1_LORENZ){
      /* Fix Lorenz gauge */
      u1_Lorenz( u1gf + 4*i, lk, lksq );
    }
    
    if(zero_mode_scheme == U1_QED_L){
      /* QED_L condition */
      if(lkssq == 0.0){
	/* All components are zero */
	int dir;
	FORALLUPDIR(dir){
	  u1gf[4*i+dir] = cmplx(0.0, 0.0);
	}
      }
      
    } else if(zero_mode_scheme == U1_QED_TL){
      /* QED_TL condition */
      if(lkssq == 0.0){
	u1gf[TUP] = cmplx(0.0, 0.0);
	int dir;
	/* Spatial components are zero only if lk[TUP] is zero */
	if(lk[TUP] == 0.0){
	  FORALLUPDIRBUT(TUP,dir){
	    u1gf[4*i+dir] = cmplx(0.0, 0.0);
	  }
	}
      }
    }
	
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
    if(mom_status[i]==U1_MOM_SELF){
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

/* ************************************************************	*/
/* Public entry points						*/
/* ************************************************************	*/

/* Prepare the momentum table and the Fourier transform.  Must be
   called once, after make_lattice(), before generate_u1_gauge_field().
   The tables do not depend on the random numbers, so this is done only
   once no matter how many fields are generated.

   NOTE: the Fourier transform setup is global state shared with the
   rest of the application.  We configure it here for a full
   four-dimensional transform.  setup_restrict_fourier() may be called
   repeatedly without leaking, so a caller that needs the same full 4D
   transform elsewhere is unaffected. */

void initialize_u1_gauge_field(u1_gauge_fix_t gauge, u1_zero_mode_t scheme){

  int key[4], slc[4];
  int dir;

  if(is_initialized)finalize_u1_gauge_field();

  gauge_fix = gauge;
  zero_mode_scheme = scheme;

  /* Momentum table */
  set_mom_table();

  /* Fourier transform over all four directions */
  FORALLUPDIR(dir){
    key[dir] = 1;
    slc[dir] = 0;
  }
  setup_restrict_fourier(key, slc);

  is_initialized = 1;

} /* end of initialize_u1_gauge_field() */

/* ************************************************************	*/

/* Generate one stochastic U(1) vector potential.

   The caller provides u1_A with room for 4*sites_on_node Reals,
   indexed as u1_A[4*i+dir], and owns the storage.  See
   create_u1_A_field() and destroy_u1_A_field() in u1link.c. */

void generate_u1_gauge_field(Real *u1_A){

  complex *u1gf;
  Real vnorm;
  register site *s;
  register int i,dir;

  if(!is_initialized){
    node0_printf("generate_u1_gauge_field: must call initialize_u1_gauge_field first\n");
    terminate(1);
  }

  if(u1_A == NULL){
    node0_printf("generate_u1_gauge_field: null u1_A\n");
    terminate(1);
  }

  /* Generate the gauge config in momentum space */

  u1gf = (complex *)malloc(sites_on_node*4*sizeof(complex));
  if(u1gf == NULL){
    node0_printf("generate_u1_gauge_field: No room for the momentum-space field\n");
    terminate(1);
  }

  momgauge(u1gf);

  /* FFT the momentum space gauge config to coordinate space */

  restrict_fourier_field(u1gf, 4*sizeof(complex), BACKWARDS);

  vnorm=(Real)sqrt((double)volume);

  FORALLSITES(i,s){
    FORALLUPDIR(dir){
      u1_A[4*i+dir] = u1gf[4*i+dir].real/vnorm;
    }
  }

  free(u1gf);

} /* end of generate_u1_gauge_field() */

/* ************************************************************	*/

/* Release the momentum table.  The Fourier transform setup is global
   state that we do not own, so it is left in place. */

void finalize_u1_gauge_field(void){

  if(mom_status != NULL){
    free(mom_status);
    mom_status = NULL;
  }
  is_initialized = 0;

} /* end of finalize_u1_gauge_field() */

/* ************************************************************	*/
