/********** w_hdibaryon.c *************/
/* MIMD version 6 */
/* 4/5/97 C. DeTar */

/* H di-baryon propagator */

/* Construction of propagators in flavor uuddss channels connected
   to the H di-baryon 

   It is assumed that the u and d quarks are degenerate.

   The spin states are taken to be nonrelativistic.

   Terms in the channel wavefunctions are enumerated in terms of
   the color and spin of diquark pairs:  Each quark has 6 color-spin
   degrees of freedom, but due to antisymmetry among a flavor pair,
   only 15 combinations of color and spin of a quark pair need be
   enumerated.  These 15 diquark color-spin indices form the basis
   of the diquark propagators.  Three such propagators are then
   combined according to the weights of the terms in the wavefunctions.

   Three channels are considered with all quarks in the same spatial
   orbital: lambda-lambda, nucleon-xi, and sigma-sigma

 */

#include "h_dibaryon_includes.h"
#include <string.h>

#define MAXCHANNEL 3
#define MAXWF 200
#define MAXLABEL 5

void w_hdibaryon(field_offset src_u, field_offset src_s, propagator prop[])
{

  diqkwf channel_wf[MAXCHANNEL][MAXWF];
  int iterm,jterm,channel_terms[MAXCHANNEL];
  Real channel_norm[MAXCHANNEL];
  char channel_label[MAXCHANNEL][MAXLABEL];
  int ichan,jchan,nchannel;
  int csui,csuj,csdi,csdj,cssi,cssj,wti,wtj;
  complex ctemp;

  int i,t;
  register site *s;
  
  complex *prop_tmp;
  Real norm;
  
  prop_tmp = (complex *)malloc(nt*sizeof(complex));
  
  /* Define the channel wave functions */
  make_channel_wfs(channel_wf, channel_terms, channel_label,  
		   channel_norm, &nchannel);

  /* Compute the channel propagators */
  for(ichan = 0; ichan < nchannel; ichan++)
    for(jchan = 0; jchan < nchannel; jchan++)
      {
	for(t=0;t<nt;t++)
	  {
	    prop_tmp[t].real = 0.;
	    prop_tmp[t].imag = 0.;
	  }

	for(iterm = 0; iterm < channel_terms[ichan]; iterm++)
	  for(jterm = 0; jterm < channel_terms[jchan]; jterm++)
	    {
	      csui = channel_wf[ichan][iterm].dq[0];
	      csuj = channel_wf[jchan][jterm].dq[0];
	      csdi = channel_wf[ichan][iterm].dq[1];
	      csdj = channel_wf[jchan][jterm].dq[1];
	      cssi = channel_wf[ichan][iterm].dq[2];
	      cssj = channel_wf[jchan][jterm].dq[2];
	      wti  = channel_wf[ichan][iterm].wt;
	      wtj  = channel_wf[jchan][jterm].wt;
		
	      FORALLSITES(i,s)
		{
		  CMUL(((dipauli_propagator *)F_PT(s,src_u))->q[csui].q[csuj],
		       ((dipauli_propagator *)F_PT(s,src_u))->q[csdi].q[csdj],
		       ctemp);
		  CMUL(((dipauli_propagator *)F_PT(s,src_s))->q[cssi].q[cssj],
		       ctemp,ctemp);
		  CMULREAL(ctemp,(Real)wti*wtj,ctemp);
		  CSUM(prop_tmp[s->t],ctemp);
		}
	    }

	norm = 1./(volume*channel_norm[ichan]*channel_norm[jchan]);
	
	for(t=0;t<nt;t++)
	  {
	    g_complexsum(&prop_tmp[t]);
	    prop[ichan + nchannel*jchan].c[t] = prop_tmp[t];
	    CMULREAL(prop[ichan + nchannel*jchan].c[t],
		     norm,
		     prop[ichan + nchannel*jchan].c[t]);
	  }

	/* Label for this channel pair */
	strcpy(prop[ichan + nchannel*jchan].label,channel_label[ichan]);
	strcat(prop[ichan + nchannel*jchan].label,channel_label[jchan]);

	/* Mark propagator completed */
	prop[ichan + nchannel*jchan].done = 1;
      }

  free(prop_tmp);

} /* w_hdibaryon */
