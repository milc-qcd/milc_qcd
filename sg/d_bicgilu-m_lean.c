/******* d_bicgilu-m_w.c - BiCGstab-M-ILU for  Wilson fermions ****/
/* MIMD version 4 */
/* BJ 01/16/97    */
/* See paper in Edinbugh Lattice Proceedings: hep-lat/9708029 */

/* Memory requirements:
   7 + N_kappa * 1.5 wilson_vectors
*/

#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <su3.h>
#include LATDEF
#include <comdefs.h>
double dtime,dclock();


#define SQUINT_FACTOR 1.2   /* see below! */

/* some lazy macros */

#define W_PT(s,pt)            ((wilson_vector *)F_PT(s,pt))
#define COPY(sr, ds, par)  FORSOMEPARITY(i,s,par) \
                                copy_wvec(W_PT(s,sr), W_PT(s,ds))
#define CLEAR(sr, par)       FORSOMEPARITY(i,s,par) \
                                clear_wvec(W_PT(s,sr))

/* this routine solves the even or odd part of the system. */
/* convergence checks and right-preconditioning are done in the 
   driver routine below */

int h_cgilu_m_w(src,dest,MaxCG,RsdCG,size_r,Kappa, num_kappa, mem, pari, k_sol)
field_offset src,   
             *dest, /* array [0..num_kappa-1] of solutions */
                    
              *mem; /* auxiliary vectors: (num_kappa+6) wilson-vectors */
int MaxCG;
float *size_r,      /* residual of the main system */
      *RsdCG;       /* array of desired residuals */
float *Kappa;       /* array of kappas, no specific order necessary */
int num_kappa;      /* # of kappas */
int pari;           /* solve even or odd half of the preconditioned system */
int k_sol;          /* index of the main system */
{
    int N_iter;         /* general loop variables */
    register int i;
    register site *s;
    register float masssq;
    int k, i1,i2;

    float size_src;     /* residues */

    double_complex pinm1[MAX_KAP],    /* normalisation parameters for the */
                   pin[MAX_KAP],      /* shifted systems */
                   pinp1[MAX_KAP],
                   rhon[MAX_KAP],
                   rhonp1[MAX_KAP],
                   cgalpha,           /* BiCG parameters */
                   cgbeta,
                   cgbetanm1,
                   cgchi,             /* MR parameter    */
                   cgalphas[MAX_KAP], /* BiCG shifted parameters */
                   cgbetas[MAX_KAP],  
                   cgchis[MAX_KAP];   /* MR shifted parameters */

    int            converged[MAX_KAP];/* 1 if system has converged */

    double_complex omega, phi,        /* further BiCGstab parameters */
                   deln, delnp1,
                   psi;

    double_complex z1,z2,z3,z4,       /* auxiliary variables */
                   zone, zzero;
    double         d1, d2;
    complex        c1, c2, c3, c4, c5;
    float          f1;

    double         sigma[MAX_KAP];    /* shift = 1/kappa_0^2 - 1/kappa^2 */
		   
                                      /* we rename the fields (convenience) */
    field_offset   y,                 
                   r,
                   rn,                /* r may contain r_n+1 */
                   ss[MAX_KAP],        
                   ms,                /* M*s */
                   w,
                   mw,
                   tmp1;              /* temporary vector */

#ifndef EVENFIRST /* Bomb the compile */
    DIE  /*You must use EVENFIRST with this routine*/
#endif

    if(even_sites_on_node!=odd_sites_on_node){
	printf("Need same number of even and odd sites on each node\n");
	terminate(1);
    }

    /* rename auxiliary vectors, else the code is unreadable */

    y  = src;   
    r  = mem[0];
    ms = mem[1];
    w  = mem[2];
    mw = mem[3];
    rn = mem[4];
    tmp1 = mem[5];
    for (k = 0; k < (num_kappa+1)/2;k++)
      {

	/* since the ss vectors only have to be even or odd parts of fields
           we can use the same trick as in d_bicgilu_lean.c to save memory */

	if (pari == EVEN)
	  {
	    ss[k] = mem[6+k];
	    ss[k+(num_kappa+1)/2] = ss[k] + even_sites_on_node*sizeof(site);
	  }
	else
	  {
	    ss[k] = mem[6+k] - even_sites_on_node*sizeof(site);
	    ss[k+(num_kappa+1)/2] = mem[6+k];
	  }
      }

    /* initialize vectors */

    COPY(src, r, pari);
    for (k=0;k < num_kappa;k++)
      {
	COPY(r, ss[k], pari);
	CLEAR(dest[k], pari);
      }

    /* COPY(r, y); see assignment of y */
    
    /* initialize variables */

    zzero.real = 0.0; zzero.imag = 0.0;
    zone.real = 1.0; zone.imag = 0.0;

    d1 = 0.0;
    FORSOMEPARITY(i, s, pari)
      d1 += magsq_wvec(W_PT(s, r));
    g_doublesum(&d1);

    CMULREAL(zone, d1, deln);      

    size_src    = (float) sqrt(d1);
    *size_r = 1.0;

/* if (this_node == 0) printf("Beginning inversion, size_src = %g\n",size_src); */
    
    for (k = 0;k < num_kappa;k++)
      sigma[k] = - (float)(1.0/((double)Kappa[k]*(double)Kappa[k])-
			   1.0/((double)Kappa[k_sol]*(double)Kappa[k_sol]));
    
    masssq = - (float)(1.0/((double)Kappa[k_sol]*(double)Kappa[k_sol]));

    cgalpha = zzero;
    cgbetanm1 = zone;
    for (k=0;k < num_kappa;k++)
      {
	pinm1[k] = zone;
	pin[k]   = zone;
	rhon[k]  = zone;

	converged[k] = 0;
      }
    for( N_iter = 0; N_iter < MaxCG && RsdCG[k_sol]  < *size_r; 
        N_iter = N_iter + 1) {


      /* calculate ms = MPRE * ss[k_sol], beta = deltan/<y, ms> */
      if (pari == EVEN) {
	dslash_w_site(ss[k_sol], tmp1, PLUS, ODD);
	dslash_w_site(tmp1, ms, PLUS, EVEN);
      } else {
	dslash_w_site(ss[k_sol], tmp1, PLUS, EVEN);
	dslash_w_site(tmp1, ms, PLUS, ODD);
      }	
      z1 = zzero;
      FORSOMEPARITY(i, s, pari) {
	scalar_mult_add_wvec(W_PT(s, ms), W_PT(s, ss[k_sol]), masssq,
			     W_PT(s, ms));
	c1 = wvec_dot(W_PT(s, y), W_PT(s, ms));
	CSUM(z1, c1);
      }
      g_dcomplexsum(&z1);
      CDIV(deln, z1, cgbeta);
      CNEGATE(cgbeta, cgbeta);

      /* calculate shifted parameters pinp1, betas */

      for (k = 0;k < num_kappa;k++)
	if (k != k_sol && !converged[k]) {
	  CSUB(pinm1[k], pin[k], z1);
	  CMUL(z1, cgalpha, z2);
	  CMUL(z2, cgbeta, z1);
	  CMULREAL(cgbeta, sigma[k], z2);
	  CSUB(zone, z2, z2);
	  CMUL(z2, cgbetanm1, z3);
	  CMUL(z3, pinm1[k], z2);
	  CADD(z1, z2, z1);
	  CDIV(pin[k], z1, z2);
	  CMUL(z2, pinm1[k], z1);
	  CMUL(z1, cgbetanm1, pinp1[k]);
	  CMUL(pinp1[k], cgbeta, z1);
	  CDIV(z1, pin[k], cgbetas[k]);
	}

      /* w = r + beta * ms */

      CMULREAL(cgbeta, 1.0, c1);
      FORSOMEPARITY(i,s,pari) 
	c_scalar_mult_add_wvec( W_PT(s, r), W_PT(s, ms), &c1, W_PT(s, w));
      
      /* mw = MPRE * w, chi = <mw, w>/|mw|^2 */

      if (pari == EVEN) {
	dslash_w_site(w, mw, PLUS, ODD);
	dslash_w_site(mw, mw, PLUS, EVEN);
      } else {
	dslash_w_site(w, mw, PLUS, EVEN);
	dslash_w_site(mw, mw, PLUS, ODD);
      }

      z1 = zzero;
      d1 = 0.0;
      FORSOMEPARITY(i,s,pari) {
	scalar_mult_add_wvec( W_PT(s, mw), W_PT(s,w), masssq, W_PT(s, mw));
	d1 += (double)magsq_wvec(W_PT(s,mw));
	c1 = wvec_dot(W_PT(s,mw), W_PT(s, w));
	CSUM(z1, c1);
      }
      g_dcomplexsum(&z1);
      g_doublesum(&d1);
      CDIVREAL(z1, d1, cgchi);
	
      /* calculate shifted parameters chis and rhonp1 */

      for (k=0;k < num_kappa;k++) 
	if (k != k_sol && !converged[k]) {
	CMULREAL(cgchi, sigma[k], z1);
	CADD(zone,z1,z1);
	CDIV(zone,z1,z2);
	CMUL(cgchi,z2,cgchis[k]);
	CMUL(rhon[k],z2,rhonp1[k]);
      }

      /* rn+1 = wn - chi * mw , deltanp1 = <y, rnp1>, 
         dest = dest - beta * s + chi * w */
      z1 = zzero;
      d1 = 0.0;
      CMULREAL(cgchi, -1.0, c1);
      CMULREAL(cgchi,  1.0, c2);
      CMULREAL(cgbeta, -1.0, c3);
      FORSOMEPARITY(i,s,pari) {
	copy_wvec(W_PT(s, r), W_PT(s, rn));
	c_scalar_mult_add_wvec(W_PT(s,w), W_PT(s, mw), &c1, W_PT(s,r));
	c4 = wvec_dot (W_PT(s,y), W_PT(s,r));
	CSUM(z1, c4);
	d1 += (double)magsq_wvec(W_PT(s,r));
	c_scalar_mult_add_wvec(W_PT(s,dest[k_sol]), W_PT(s, ss[k_sol]), &c3, 
			       W_PT(s,dest[k_sol]));
	c_scalar_mult_add_wvec(W_PT(s,dest[k_sol]), W_PT(s, w), &c2,
			       W_PT(s,dest[k_sol]));
      }
      g_dcomplexsum(&z1);
      g_doublesum(&d1);
      *size_r = sqrt(d1)/size_src;
      delnp1 = z1;
      /* update the shifted solutions */

      for (k=0;k < num_kappa;k++)
	if (k != k_sol && !converged[k]) {
	  CNEGATE(cgbetas[k], c1);
	  CMUL(cgchis[k], rhon[k], z1);
	  CMUL(z1, pinp1[k], c2);
	  FORSOMEPARITY(i,s,pari) {
	    c_scalar_mult_add_wvec(W_PT(s, dest[k]), W_PT(s, ss[k]), &c1,
				   W_PT(s, dest[k]));
	    c_scalar_mult_add_wvec(W_PT(s, dest[k]), W_PT(s, w), &c2,
				   W_PT(s, dest[k]));
	  }
	}

      /* calculate alpha */

      CMUL(deln, cgchi, z1);
      CDIV(delnp1, z1, cgalpha);
      z1 = cgalpha;
      CMUL(z1, cgbeta, cgalpha);
      CNEGATE(cgalpha, cgalpha);
      /* calculate shifted parameters */

      for (k = 0;k < num_kappa;k++) 
	if (k != k_sol && !converged[k]) {
	  CMUL(pinp1[k], cgbetas[k], z1);
	  CMUL(z1, cgalpha, z2);
	  CDIV(z2, pin[k], z1);
	  CDIV(z1, cgbeta, cgalphas[k]);
	}

      /* s = wn + alpha*(s - chi*ms) */

      CNEGATE(cgchi, c1);
      CMULREAL(cgalpha, 1.0, c2);
      FORSOMEPARITY(i,s,pari) {
	c_scalar_mult_add_wvec(W_PT(s, ss[k_sol]), W_PT(s, ms), &c1, 
			       W_PT(s, ss[k_sol]));
	c_scalar_mult_add_wvec(W_PT(s, r), W_PT(s,ss[k_sol]), &c2, 
			       W_PT(s, ss[k_sol]));
      }

      /* now update s for the shifted systems */

      for (k=0;k < num_kappa;k++)
	if (k != k_sol && !converged[k]) {
	  CMUL(pin[k], rhon[k], z1);
	  CNEGATE(z1, c1);
	  CMUL(pinp1[k], rhon[k], c2);
	  CDIV(cgchis[k], cgbetas[k], c3);
	  CNEGATE(c3,c3);
	  CMUL(pinp1[k], rhonp1[k], c4);
	  CMULREAL(cgalphas[k], 1.0, c5);

	  FORSOMEPARITY(i,s,pari) {
	    clear_wvec(W_PT(s, tmp1));
	  }
	  FORSOMEPARITY(i,s,pari) {

	    c_scalar_mult_add_wvec(W_PT(s,tmp1), W_PT(s, w), &c2, W_PT(s, ms));
	    c_scalar_mult_add_wvec(W_PT(s, ms), W_PT(s,rn), &c1, W_PT(s, ms));

	    c_scalar_mult_add_wvec(W_PT(s,ss[k]),W_PT(s, ms), &c3,W_PT(s,ss[k]));
	    c_scalar_mult_add_wvec(W_PT(s,tmp1),W_PT(s,ss[k]),&c5,W_PT(s,ss[k]));
	    c_scalar_mult_add_wvec(W_PT(s,ss[k]),W_PT(s,r), &c4, W_PT(s,ss[k]));
	  }
	}
      /* n -> n+1 */

      deln = delnp1;
      cgbetanm1 = cgbeta;
      for (k=0;k < num_kappa;k++)
	if (k != k_sol && !converged[k]) {
	  pinm1[k] = pin[k];
	  pin[k] = pinp1[k];
	  rhon[k] = rhonp1[k];
	  CMUL(pin[k], rhon[k], z1);
	  if (cabs(&z1)*(*size_r) < RsdCG[k] / 5.0)
	    converged[k] = 1;    /* in case that the real residue is more than
				    a factor 5 (maybe less?) greater than the 
				    estimated residue, we dont expect the 
				    system to converge any more anyways */
	}
/*      if (this_node == 0)
	printf("Iter : %d Resd : %g\n", N_iter, *size_r); fflush(stdout); */

    }  

    /* adjust normalisation */

    for (k = 0; k < num_kappa;k++)
      {
	f1 =  - (float)1.0/((double)Kappa[k]*(double)Kappa[k]);
	FORSOMEPARITY(i,s,pari)
	  scalar_mult_wvec(W_PT(s,dest[k]), f1, W_PT(s, dest[k]));
      }

    return N_iter;
}

/* the driver routine. performs convergence checks and preconditioning */
/* note: since we abuse the quark-propagator storage as auxiliary storage,
   we prefer to provide this storage as a parameter to avoid strange 
   side-effects. */

int cgilu_m_w(src,dest,MaxCG,RsdCG,size_r,Kappa, num_kappa, mem, sourcekind)
field_offset src,   
             *dest, /* array [0..num_kappa-1] of solutions */
                    
              *mem; /* auxiliary vectors: (num_kappa+x) wilson-vectors 
		       NOTE: mem[0] MUST NOT BE USED BY cgilu_w !!!!! */
int MaxCG;
float *size_r,      /* array of residuals !!! */
      *RsdCG;       /* array of desired residuals */
float *Kappa;       /* array of kappas, no specific order necessary */
int num_kappa;      /* # of kappas */
int sourcekind;     /* EVEN = source only on even sites, ODD only on odd,
                       EVENANDODD = general source */
{
  double dt;        /* time variable */
  int k,i;
  site *s;
  int N_iter, k_sol;
  float rsd;
  double d1;
  float size_src;

  /* select the main system to solve */

  for (k_sol = 0,k = 1;k < num_kappa;k++) /* we solve the system with the */
    if (Kappa[k] > Kappa[k_sol])          /* largest kappa */
      k_sol = k;

  /* start timing */

  dt = - dclock();

  
  if (sourcekind == EVEN)
    {
      N_iter = h_cgilu_m_w(src, dest, MaxCG, RsdCG, &rsd, Kappa, num_kappa, mem, 
			     EVEN, k_sol);
      for (k=0;k < num_kappa;k++)
	CLEAR(dest[k], ODD);
    }
  else if (sourcekind == ODD)
    {
      N_iter = h_cgilu_m_w(src, dest, MaxCG, RsdCG, &rsd, Kappa, num_kappa, mem, 
			     ODD, k_sol);
      for (k=0;k < num_kappa;k++)
	CLEAR(dest[k], EVEN);
    }
  else
    {
      N_iter = h_cgilu_m_w(src, dest, MaxCG, RsdCG,&rsd, Kappa, num_kappa, mem, 
			     EVEN, k_sol);
      N_iter +=  h_cgilu_m_w(src, dest, MaxCG, RsdCG, &rsd, Kappa, num_kappa, mem, 
			       ODD, k_sol);
    }
  
  /* apply preconditioning, multiply results with 1+kappa Dslash */

  for (k=0;k < num_kappa;k++)
    {
      dslash_w_site(dest[k], mem[0], PLUS, EVENANDODD);
      FORALLSITES(i,s) 
	scalar_mult_add_wvec(W_PT(s,dest[k]), W_PT(s,mem[0]), Kappa[k], 
			     W_PT(s,dest[k]));
    }

  /* calculate source size */
  
  d1 = 0.0;
  FORALLSITES(i,s)
    d1 += magsq_wvec(W_PT(s,src));
  g_doublesum(&d1);
  size_src = (float)sqrt(d1);

  /* now we can check for convergence of each single system */
  
  for (k = 0;k < num_kappa;k++)
    {
      dslash_w_site(dest[k], mem[0], PLUS, EVENANDODD);
      d1 = 0.0;
      FORALLSITES(i,s) {
	scalar_mult_add_wvec(W_PT(s,dest[k]), W_PT(s,mem[0]), -Kappa[k],
			     W_PT(s,mem[0]));
	sub_wilson_vector(W_PT(s,src), W_PT(s, mem[0]), W_PT(s,mem[0]));
	d1 += (double) magsq_wvec(W_PT(s, mem[0]));
      }
      g_doublesum(&d1);
      size_r[k] = sqrt(d1)/size_src;

      /* check for convergence. squint_factor takes into account that the
	 true residual can be somewhat different from the iterated residual.
	 we will not try to improve the main system */

      if (size_r[k] > RsdCG[k] * SQUINT_FACTOR && k != k_sol)
	{
	  if (this_node == 0)
	    {
	      printf("CGILU-M: Kappa %d not converged: %g, calling BiCGstab\n",
		     k, size_r[k]); fflush(stdout);
	    }
	  /* cgilu_w destroys the source, so save it */
	  COPY(src, mem[0], EVENANDODD);
	  N_iter += 
	    cgilu_w(mem[0], dest[k], 20, RsdCG[k], &(size_r[k]), 1, Kappa[k]);
	}
    }

  dt += dclock();

  if (this_node == 0)
    {
      if (N_iter != 0) 
	printf("CGILU-M: time = %e = %e/site-iter\n",
	       dt, dt/((N_iter+2*num_kappa)*volume));
      else
	printf("CGILU-M: NO ITERATIONS TAKEN\n");
      fflush(stdout);
    }

  for (k = 0;k < num_kappa;k++)
    if (size_r[k] > RsdCG[k] * SQUINT_FACTOR)
      if (this_node == 0)
	{
	  printf("CGILU-M: Kappa %d Not Converged\n", k); fflush(stdout);
	}

  return N_iter;
}
