/********** w_nrbaryon.c *************/
/* MIMD version 6 */
/* C. DeTar 4/6/97 */

/* Construct baryon propagators for the baryon octet
   based on a nonrelativistic treatment of the spin states.

   For the nucleon, xi, and sigma we can use a diquark propagator
   together with a quark propagator.  All three states are
   calculated as a generalized "nucleon" with two like quark flavors
   and one distinct flavor.  The propagator for src_nr1 is
   the distinct flavor and the diquark src_di2 is used for
   the two like quark flavors.

   For the lambda, we use three quark propagators. In this case
   the propagator for src_nr1 is treated as the s quark
   and the propagator for src_nr2 is used for both u and d quarks.
*/

#include "h_dibaryon_includes.h"
#include <string.h>

void make_nucleon_wf(twoplusoneqkwf wf[],int *terms, Real *norm)
{
  /* See diquarprop.c for translation from diquark to quark color-spin
     indices */

  wf[0].dq =  5 ; wf[0].cs = 3 ; wf[0].wt =  2 ;
  wf[1].dq =  8 ; wf[1].cs = 0 ; wf[1].wt = -1 ;
  wf[2].dq =  6 ; wf[2].cs = 2 ; wf[2].wt =  1 ;
  wf[3].dq = 10 ; wf[3].cs = 0 ; wf[3].wt =  1 ;
  wf[4].dq =  9 ; wf[4].cs = 1 ; wf[4].wt = -1 ;
  wf[5].dq =  0 ; wf[5].cs = 5 ; wf[5].wt =  2 ;
  wf[6].dq =  1 ; wf[6].cs = 4 ; wf[6].wt = -2 ;
  wf[7].dq =  3 ; wf[7].cs = 2 ; wf[7].wt = -1 ;
  wf[8].dq =  4 ; wf[8].cs = 1 ; wf[8].wt =  1 ;

  *terms = 9;
  *norm = 3*sqrt(2.);
}

void make_lambda_wf(threeqkwf wf[], int *terms, Real *norm)
{

  wf[ 0].cs[0] = 1; wf[ 0].cs[1] = 5; wf[ 0].cs[2] = 0; wf[ 0].wt = -1; 
  wf[ 1].cs[0] = 1; wf[ 1].cs[1] = 3; wf[ 1].cs[2] = 2; wf[ 1].wt =  1;
  wf[ 2].cs[0] = 2; wf[ 2].cs[1] = 4; wf[ 2].cs[2] = 0; wf[ 2].wt =  1;
  wf[ 3].cs[0] = 2; wf[ 3].cs[1] = 3; wf[ 3].cs[2] = 1; wf[ 3].wt = -1;
  wf[ 4].cs[0] = 0; wf[ 4].cs[1] = 4; wf[ 4].cs[2] = 2; wf[ 4].wt = -1;
  wf[ 5].cs[0] = 0; wf[ 5].cs[1] = 5; wf[ 5].cs[2] = 1; wf[ 5].wt =  1;
  wf[ 6].cs[0] = 4; wf[ 6].cs[1] = 2; wf[ 6].cs[2] = 0; wf[ 6].wt =  1;
  wf[ 7].cs[0] = 4; wf[ 7].cs[1] = 0; wf[ 7].cs[2] = 2; wf[ 7].wt = -1;
  wf[ 8].cs[0] = 5; wf[ 8].cs[1] = 1; wf[ 8].cs[2] = 0; wf[ 8].wt = -1;
  wf[ 9].cs[0] = 5; wf[ 9].cs[1] = 0; wf[ 9].cs[2] = 1; wf[ 9].wt =  1;
  wf[10].cs[0] = 3; wf[10].cs[1] = 1; wf[10].cs[2] = 2; wf[10].wt =  1;
  wf[11].cs[0] = 3; wf[11].cs[1] = 2; wf[11].cs[2] = 1; wf[11].wt = -1;

  *terms = 12;
  *norm = 2*sqrt(3.);
}

void w_nrbaryon(field_offset src_1, 
		field_offset src_2, field_offset src_2di, 
		propagator prop[]) 
{

  int i,t;
  register site *s;

  twoplusoneqkwf nucleon[9];
  threeqkwf lambda[12];

  int iterm,jterm,nucleon_terms,lambda_terms;
  int csuui,csuuj,csdi,csdj;
  int csui,csuj,cssi,cssj;
  int cui,cuj,cdi,cdj,csi,csj;
  int sui,suj,sdi,sdj,ssi,ssj;
  int wti,wtj;
  complex ctemp;
  complex *prop_tmp;
  Real norm,nucleon_norm,lambda_norm;

  prop_tmp = (complex *)malloc(nt*sizeof(complex));


  /* NUCLEON, i.e. a q1 q1 q2 state in the spin 1/2 octet */

  make_nucleon_wf(nucleon,&nucleon_terms,&nucleon_norm);

  for(t=0;t<nt;t++)
    {
      prop_tmp[t].real = 0.;
      prop_tmp[t].imag = 0.;
    }
  
  for(iterm = 0; iterm < nucleon_terms; iterm++)
    {
      /* We pretend the state is uud, but the
	 same wf works for the sigma and xi */
      /* So the diquark is uu */
      csuui = nucleon[iterm].dq;
      /* and the third quark is d */
      csdi  = nucleon[iterm].cs;
      /* Decode color and spin for the d quark */
      cdi = csdi % 3; sdi = csdi/3;

      for(jterm = 0; jterm < nucleon_terms; jterm++)
	{
	  csuuj = nucleon[jterm].dq;
	  csdj  = nucleon[jterm].cs;
	  wti   = nucleon[iterm].wt;
	  wtj   = nucleon[jterm].wt;
	  cdj = csdj % 3; sdj = csdj/3;
	  
	  FORALLSITES(i,s)
	    {
	      /* Multiply the diquark propagator by the quark propagator */
	      CMUL(((dipauli_propagator *)F_PT(s,src_2di))->q[csuui].q[csuuj],
		   ((pauli_propagator *)
		   F_PT(s,src_1))->c[cdi].p[sdi].p[sdj].c[cdj],
		   ctemp);
	      CMULREAL(ctemp,(Real)wti*wtj,ctemp);
	      CSUM(prop_tmp[s->t],ctemp);
	    }
	}
    }
  
  norm = 1./(volume*nucleon_norm*nucleon_norm);
  
  for(t=0;t<nt;t++)
    {
      g_complexsum(&prop_tmp[t]);
      prop[0].c[t] = prop_tmp[t];
      CMULREAL(prop[0].c[t],norm,prop[0].c[t]);
    }
  
  /* Label for this channel pair */
  strcpy(prop[0].label,"NUC");
  
  /* Mark it done */
  prop[0].done = 1;

  /* LAMBDA */

  make_lambda_wf(lambda,&lambda_terms,&lambda_norm);

  for(t=0;t<nt;t++)
    {
      prop_tmp[t].real = 0.;
      prop_tmp[t].imag = 0.;
    }
  
  for(iterm = 0; iterm < lambda_terms; iterm++)
    {
      csui = lambda[iterm].cs[0];
      csdi = lambda[iterm].cs[1];
      cssi = lambda[iterm].cs[2];
      wti   = lambda[iterm].wt;
      /* Decode color and spin labels */
      cui = csui % 3; sui = csui/3;
      cdi = csdi % 3; sdi = csdi/3;
      csi = cssi % 3; ssi = cssi/3;
      
      for(jterm = 0; jterm < lambda_terms; jterm++)
	{
	  csuj = lambda[jterm].cs[0];
	  csdj = lambda[jterm].cs[1];
	  cssj = lambda[jterm].cs[2];
	  wtj   = lambda[jterm].wt;
	  
	  cuj = csuj % 3; suj = csuj/3;
	  cdj = csdj % 3; sdj = csdj/3;
	  csj = cssj % 3; ssj = cssj/3;
	  
	  FORALLSITES(i,s)
	    {
	
	      /* Multiply three quark propagators */
	      CMUL(((pauli_propagator *)
		   F_PT(s,src_2))->c[cui].p[sui].p[suj].c[cuj],
		   ((pauli_propagator *)
		   F_PT(s,src_2))->c[cdi].p[sdi].p[sdj].c[cdj],
		   ctemp);
	      
	      CMUL(ctemp,
		   ((pauli_propagator *)
		   F_PT(s,src_1))->c[csi].p[ssi].p[ssj].c[csj],
		   ctemp);
	      
	      CMULREAL(ctemp,(Real)wti*wtj,ctemp);
	      CSUM(prop_tmp[s->t],ctemp);
	    }
	}
    }

  norm = 1./(volume*lambda_norm*lambda_norm);
  
  for(t=0;t<nt;t++)
    {
      g_complexsum(&prop_tmp[t]);
      prop[1].c[t] = prop_tmp[t];
      CMULREAL(prop[1].c[t],norm,prop[1].c[t]);
    }
  
  /* Label for this channel pair */
  strcpy(prop[1].label,"LAM");

  /* Mark it done */
  prop[1].done = 1;

  free(prop_tmp);

}  /* w_nrbaryon */

