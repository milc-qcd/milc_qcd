/****************** bar_corr11.c ************************************/
/* MIMD version 6 */
/* NOT MAINTAINED.  TEST BEFORE USE! */
/* C. DeTar 20 Nov 1995 Increased NRANDMAX and QRAND for 12^3 x 6 studies
/* C. DeTar 15 Jan 1995 using LU preconditioning for propagators */
/* C. DeTar 23 Mar 1994 rearranged output format */
/* C. DeTar 22 Mar 1994 measures both nonfuzzed correlation in place of slice */
/* C. DeTar 18 Jan 1994 generalized fft to allow dimension p*n^k */
/* C. DeTar 29 Jul 1993 modified to do fuzzy ploop correlations */
/* C. DeTar  2 Apr 1992 modified to do psi-bar-psi as well */
/* C. DeTar 14 Feb 1992 modified to do line density as well as slice */
/* C. DeTar 29 Sept 1991 with Fourier transform and symmetry folding by Doug Toussaint */
/* Evaluate baryon density, Polyakov loop correlations */
/* This version uses a Fourier transform and a random source */
/* This version also measures the single-color-trace density-density 
   correlation */

/* It is assumed that ploop_staple has just been called and 
   the fuzzy Polyakov loop link matrix products are in tempmat1 on even sites 
   in the first two time slices. Note that these need to be normalized 
   by multiplying by (1 - alpha)^N_t.  This is done at the end when the 
   correlations are finally normalized */

/* This procedure computes eight correlations:
   (1) the Polyakov loop/ Polyakov loop correlation (1 value)
   (2) the correlation between a static test quark and a fuzzed test quark
   and the dynamical baryon density, measured as an average along the time
   direction (line density). (2 values)   
   (3) three terms contributing to the density-density correlation,
   for line/line correlations  (3 values)

   Correlation (1) includes the KS phases and correlates the imaginary parts
   of the Polyakov loop.

   Correlation (2) is defined in terms of two quantities, 
   the static Polyakov loop P_stat (fuzzed and nonfuzzed) and the 
   "dynamical" Polyakov loop P_dyn.  The static loop is the usual
   one, except it is defined here including the KS phases (and 
   antiperiodic b.c. phase).  The dynamical one is defined as 
   P_dyn(r) = Im <q bar(r,t+1) eta U_t q(r,t)> for any choice of t
   (here we average over t for
   the line density).  eta is the KS phase.  The Im takes account
   of the contribution of both quarks and antiquarks to the density.  
   Thus this expectation value involves computing the quark
   propagator from t to t + 1 and then multiplying by the link
   matrix U.
   
   The correlation computed here is the convolution of f_s(r) = Im P_stat(r) 
   and f_d(r) = Im P_dyn(r).  The convolution is obtained from the 
   Fourier transforms

        bpcorr(k) = - f_s^*(k) f_d(k)/2

   Its Fourier transform back to r is reported as a list of values as a 
   function of the offset vector r between the operators.  To conserve space 
   the computed values are averaged over symmetry-related offsets.

   Since computing f_d(k) would require knowing Im P_dyn(r) at every point, 
   which is too expensive, we use the random vector trick and compute
   
        f_d(k) = sum_rs e^(ikr) X^*(s) Im Tr [M^{-1}(s,t+1;r,t) eta U_t] X(r)

   where X(r) is a Gaussian random vector.  Averaging over X gives the
   desired result.  This routine repeats the calculation for NRAND X's
   for each configuration.

   Correlation (3) has three contributions: 
   (a) a single point, single trace term S11
   (b) a two point, single trace term S21
   (c) a two point, two trace term S22
   The first of these is related to Re P_dyn and contributes to the correlation 
   only at zero separation.  The second is essentially 
   a hadron propagator from a baryon density source to a baryon density sink
   and the third is a correlation between Im P_dyn and Im P_dyn.

   Output pattern

   BCOR0 bbs11 bb0s11
                   has two values, giving the term S11 for the 
                   slice-density/slice-density correlation
   BCOR1 x y z LPAVG LFAVG LLS22 LLS21 
                   line density/nonfuzz correlation, 
		   line density/ploop fuzz correlation,
		   term S22 of line density/line density,
		   term S21 of line density/line density
   BCOR2 x y z RRAVG IIAVG FFAVG GGAVG 
                   Re ploop nonfuzz/Re ploop nonfuzz correlation, 
		   Im ploop nonfuzz/Im ploop nonfuzz correlation
                   Re ploop fuzz/Re ploop fuzz correlation, 
                   Im ploop fuzz/Im ploop fuzz correlation, 
   BCOR3 x y z SPAVG SFAVG SSAVG 
                   pbp - Polyakov loop nonfuzz correlation
                   pbp - Polyakov loop fuzz correlation
                   pbp - pbp correlation
  
   */ 

#include "ks_dyn_includes.h"
#include "defines.h"

#define TRUE 1
#define FALSE 0
#define restrict rstrict /* C-90 T3D cludge */
/* This is where we stash some complex numbers needed in the computation */
#define PLPRE tempmat2.e[0][0]   /* Re part of Polyakov loop */
#define PLPIM tempmat2.e[0][1]   /* Polyakov loop */
#define FLPRE tempmat2.e[0][2]   /* Re part of fuzzy Polyakov loop */
#define FLPIM tempmat2.e[1][0]   /* Im part of fuzzy Polyakov loop */
#define TEMP1 tempmat2.e[1][1]   /* complex */
#define TEMP2 tempmat2.e[1][2]   /* complex */
#define PBP   tempmat2.e[2][0]   /* Psi-bar-psi (single site) */
#define BBS21 tempmat2.e[2][1]
#define DEBUG tempmat2.e[2][2]
#define LPAVG avg[0]   /* Baryon line density/nonfuzz P loop corrln average */
#define LFAVG avg[1]   /* Baryon line density/Polyakov loop corrln average */
#define LLS22 avg[2]   /* 2 point 2 trace contribution to line density/line density corrln */
#define LLS21 avg[3]   /* 2 point, 1 trace contribution to line density/line density corrln */
#define RRAVG avg[4]   /* nonfuzz Re P - nonfuzz Re P correlation */
#define IIAVG avg[5]   /* nonfuzz Im P - nonfuzz Im P correlation */
#define FFAVG avg[6]   /* fuzz Re P- fuzz Re P correlation */
#define GGAVG avg[7]   /* fuzz Im P- fuzz Im P correlation */
#define SPAVG avg[8]   /* pbp - nonfuzz Polyakov loop correlation */
#define SFAVG avg[9]   /* pbp -fuzz ploop correlation */
#define SSAVG avg[10]  /* pbp - pbp correlation */

#define NRAND NBPRAND  /* Min no of random sources. NBPRAND is defined in lattice.h */
#define NRANDMAX NBPRAND*32 /* Maximum number of random sources */
#define QRAND 32      /* Number of random sources for psi-bar-psi measurement */
#define MAXVAR 40.      /* Maximum variance tolerated in integrated LF correlation */

  
/* Global variables */
  
/* Time slices for calculation */
int t0,t0p1;
/* index for gathering Polyakov loop variable from time slice 1*/
int sl_swap_dir;
/* Gathers for symmetry transforms.) */ 
int xy_dir,xz_dir,yz_dir,minus_dir[3];

/************ slice_swap *************************************************/
/* Mapping for gather to swap values between even and odd time slices  */
void slice_swap(int x,int y,int z,int t,int *arg,
		int forw_back,int *xpt,int *ypt,int *zpt,int *tpt)
{
  *xpt = x;  *ypt = y;  *zpt = z;
  /* Swaps t = 0 with t = 1, t = 2 with t = 3, etc. */
  if(t % 2 ==0)*tpt = (t + 1) % nt;
  else *tpt = t - 1;
} /* slice_swap */

/************ set_slice_swap *********************************************/
/* Set up gathers for swapping values between even and odd time slices */
void set_slice_swap()
{
  int dummy;
  sl_swap_dir = make_gather( slice_swap, &dummy, NO_INVERSE, 
			  ALLOW_EVEN_ODD, SCRAMBLE_PARITY);
}  /* set_slice_swap */

/************ slice_switch_map *********************************************/
/* From Doug's ks wave function code */
/* Function that defines switching two coordinates (arg tells which two) */
/* Map only on time slice t = t0p1 */
void slice_switch_map(int x,int y,int z,int t,int *arg,
		      int fb,int *xp,int *yp,int *zp,int *tp)
{
  int tt,d[3];
  d[XUP]=x; d[YUP]=y; d[ZUP]=z;
  tt = d[arg[0]];
  d[arg[0]]=d[arg[1]];
  d[arg[1]]=tt;
  *tp = t;
  if(t==t0p1)
    {
      *xp=d[XUP]; *yp=d[YUP]; *zp=d[ZUP];
    }
  else
    {
      *xp = x; *yp = y; *zp = z;
    }
} /* slice_switch_map */

/************ slice_minus_map *********************************************/
/* From Doug's ks wave function code */
/* Function that changes "sign" of selected coordinates. (arg tells which) */
/* Map only on time slice t = t0p1 */
/* if arg[i]==1, coord[i] -> dim[i]-coord[i], except 0 -> 0 */
void slice_minus_map(int x,int y,int z,int t,int *arg,
		     int fb,int *xp,int *yp,int *zp,int *tp)
{
  *tp = t;
  if(t==t0p1)
    {
      if(arg[XUP]==1){ if(x==0)*xp=0; else *xp = nx-x; } else *xp = x;
      if(arg[YUP]==1){ if(y==0)*yp=0; else *yp = ny-y; } else *yp = y;
      if(arg[ZUP]==1){ if(z==0)*zp=0; else *zp = nz-z; } else *zp = z;
    }
  else
    {
      *xp = x; *yp = y; *zp = z;
    }
} /* slice_minus_map */

/************ set_sl_sym_gath *********************************************/
/* From Doug's ks wave function code */
void set_sl_sym_gath()
{
  int i,key[3],arg[2];

  if(nx==ny){  /* interchange x and y coordinates */
    arg[0]=XUP; arg[1]=YUP;
    xy_dir = make_gather( slice_switch_map, arg, OWN_INVERSE,
			 NO_EVEN_ODD, SAME_PARITY);
  }
  if(nx==nz){  /* interchange x and z coordinates */
    arg[0]=XUP; arg[1]=ZUP;
    xz_dir = make_gather( slice_switch_map, arg, OWN_INVERSE,
			 NO_EVEN_ODD, SAME_PARITY);
  }
  if(ny==nz){  /* interchange y and z coordinates */
    arg[0]=YUP; arg[1]=ZUP;
    yz_dir = make_gather( slice_switch_map, arg, OWN_INVERSE,
			 NO_EVEN_ODD, SAME_PARITY);
  }
  for(i=XUP;i<TUP;i++){
    key[XUP]=0; key[YUP]=0; key[ZUP]=0;
    key[i]=1;
    minus_dir[i] = make_gather( slice_minus_map, key, OWN_INVERSE,
			       NO_EVEN_ODD, SAME_PARITY);
  }
} /* set_sl_sym_gath */

/************ grslice  *********************************************/
/* construct a gaussian random vector on time slice t, put it in g_rand  */
void grslice(int t) {
  register int i,j;
  register site *s;
  
  FORALLSITES(i,s){
    if(s->t!=t)continue;
    for(j=0;j<3;j++){
#ifdef SITERAND
      s->g_rand.c[j].real = gaussian_rand_no(&(s->site_prn));
      s->g_rand.c[j].imag = gaussian_rand_no(&(s->site_prn));
#else
      s->g_rand.c[j].real = gaussian_rand_no(&node_prn);
      s->g_rand.c[j].imag = gaussian_rand_no(&node_prn);
#endif
    }
  }
}/* grslice */

/************ grfull  *********************************************/
/* construct a gaussian random vector on the full lattice, put it in g_rand  */
void grfull() {
  register int i,j;
  register site *s;
  
  FORALLSITES(i,s){
    for(j=0;j<3;j++){
#ifdef SITERAND
      s->g_rand.c[j].real = gaussian_rand_no(&(s->site_prn));
      s->g_rand.c[j].imag = gaussian_rand_no(&(s->site_prn));
#else
      s->g_rand.c[j].real = gaussian_rand_no(&node_prn);
      s->g_rand.c[j].imag = gaussian_rand_no(&node_prn);
#endif
    }
  }
}/* grfull */

/************ slice_sym_combine *********************************************/
/* Fetch four contiguous complex values of AVG at symmetry related sites 
   gathered in direction "dir".  Average them if avg = 1 and leave result in
   avg[j].  Otherwise, leave gathered values in dest */
int slice_sym_combine(int dir,field_offset dest,int avg) 
{
  msg_tag *tag;
  register int i,j;
  register site *s;

  tag = start_gather_site( F_OFFSET(avg[0]), NBPAVRG*sizeof(complex), dir,
		     EVENANDODD, gen_pt[0]);
  wait_gather(tag);
  FORALLSITES(i,s)
    {
      if(s->t!=t0p1)continue;
      for (j=0; j<NBPAVRG; j++)
	{
	  ((complex *)(F_PT(s,dest)))[j] = ((complex *)gen_pt[0][i])[j];
	}
    }
  
  if(avg==1)
    {
      FORALLSITES(i,s)
	{ 
	  if(s->t!=t0p1)continue;
	  for(j = 0; j < NBPAVRG; j++)
	    {
	      CSUM( s->avg[j], ((complex *)F_PT(s,dest))[j]);
	      CMULREAL( s->avg[j], 0.5, s->avg[j]);
	    }
	}
    }
  
  cleanup_gather(tag);
}

/************ write_corr *********************************************/
write_corr()
{
  /* Doug Toussaint's */
  /* this subroutine writes the correlations */

  int currentnode,newnode;
  int i,j,l,x,y,z,t;
  complex lbuf[NBPAVRG];
  int node0=0;
  
  g_sync();
  currentnode=0;
  
  /* write the elements */
  for(z=0;z<=nz/2;z++)for(y=0;y<=ny/2;y++)for(x=0;x<=nx/2;x++) 
    {
      
      /* Skip elements that are redundant due to symmetries */
      if(nx==ny && y<x )continue;
      if(nx==nz && z<x )continue;
      if(ny==nz && z<y )continue;
      
      newnode=node_number(x,y,z,t0p1);
      if(newnode != currentnode) {	/* switch to another node */
	g_sync();
	currentnode=newnode;
      }
      
      if(this_node==0) {
	if(currentnode==0) {
	  l=node_index(x,y,z,t0p1);
	  for (j = 0; j < NBPAVRG; j++)
	    lbuf[j] = lattice[l].avg[j];
	}
	else {
	  get_field((char *)&lbuf,NBPAVRG*sizeof(complex),currentnode);
	}
	
	/* See defines for avg[] above for the order of the 
	   output correlation values */
	if( (printf("BCOR1 %d %d %d\t%.7e\t%.7e\t%.7e\t%.7e\n",
		    x,y,z,
		    (double)lbuf[0].real, (double)lbuf[1].real, 
		    (double)lbuf[2].real, (double)lbuf[3].real)== EOF)) {
	  printf("bar_corr: Write error\n"); 
	  terminate(1);
	} 
	if( (printf("BCOR2 %d %d %d\t%.7e\t%.7e\t%.7e\t%.7e\n",
		    x,y,z,
		    (double)lbuf[4].real, (double)lbuf[5].real,
		    (double)lbuf[6].real, (double)lbuf[7].real)== EOF)) {
	  printf("bar_corr: Write error\n"); 
	  terminate(1);
	} 
	if( (printf("BCOR3 %d %d %d\t%.7e\t%.7e\t%.7e\n",
		    x,y,z,
		    (double)lbuf[8].real, 
		    (double)lbuf[9].real,
		    (double)lbuf[10].real)== EOF)) {
	  printf("bar_corr: Write error\n"); 
	  terminate(1);
	} 
      }
      else {	/* for nodes other than 0 */
	if(this_node==currentnode) {
	  if(this_node==0) printf("this node is zero\n");
	  l=node_index(x,y,z,t0p1);
	  for (j = 0; j < NBPAVRG; j++)
	    lbuf[j] = lattice[l].avg[j];
	  send_field((char *)&lbuf,NBPAVRG*sizeof(complex),node0);
	}
      }
    }
  g_sync();
  if(this_node==0) fflush(stdout);

}

/************ line_avg *******************************************/

void line_avg(field_offset source,field_offset dest,field_offset temp1)
     /* Sums complex lattice values along time direction in source field 
	and puts total in time slice t0p1 in dest field */
     /* Uses the field temp1 for work space */
     /* Assumes t0 = 0 and t0p1 = 1 */
     /* Patterned after ploop */
{
  int t;
  register int i;
  register site *st;
  msg_tag *tag;
  int d[4];
  
  /* First add "source" value on every even site to value above it */
  
  tag=start_gather_site( source, sizeof(complex),
		   TUP, EVEN, gen_pt[0] );
  wait_gather(tag);
  FOREVENSITES(i,st)
    {
      CADD((*(complex *)(F_PT(st,source))),
	   (*(complex *)gen_pt[0][i]),
	   (*(complex *)(F_PT(st,temp1))));
    }
  cleanup_gather(tag);
  
  /* Next total all paired results, putting the sums on even sites in the
     first two time slices */
  
  d[XUP] = d[YUP] = d[ZUP] = 0;
  for(t=2;t<nt;t+=2)
    {
      d[TUP] = t;	/* distance from which to gather */
      tag=start_general_gather_site( temp1, sizeof(complex),
			       d, EVEN, gen_pt[0] );
      wait_general_gather(tag);
      FOREVENSITES(i,st)
	{
	  if( st->t > 1 )continue;  /* compute only on first two slices */
	  CSUM((*(complex *)(F_PT(st,temp1))),
	       (*(complex *)gen_pt[0][i]));
	  /* We overwrite temp1 on the first two time slices,
	     leaving the others undisturbed so we can still gather
	     them. */
	}
      cleanup_general_gather(tag);
    }
  
  /* Finally move totals in temp1 to the time 1 slice and average them */

  /* Move result in temp1 up to time slice 1*/
  
  tag = start_gather_site( temp1 , sizeof(complex), 
		     sl_swap_dir, ODD, gen_pt[0]);
  wait_gather(tag);

  FORODDSITES(i,st)
    {
      if( st->t != 1)continue;
      *(complex *)(F_PT(st,temp1)) = *(complex *)gen_pt[0][i];
    }

  FORALLSITES(i,st)
    {
      if( st->t != 1)continue;
      CMULREAL((*(complex *)(F_PT(st,temp1))),
	       1./((Real)nt),
	       (*(complex *)(F_PT(st,dest))));
    }

  cleanup_gather(tag);
}  


/************ sngl_trace *******************************************/
int sngl_trace(Real *bb0s11)
     /* For the density-density correlation (two-point, single trace term) */
     /* (Also computes the single slice, single trace term S11 in bb0s11) */
     /* Computes propagator between baryon density source and sink */
     /* There are two terms.  One correlates quark with quark and antiquark */
     /* with antiquark, the other quark with antiquark. */
{
  Real mass_x2;
  int cgn,icol,jcol,this_cgn;
  int x,y,z,t;
  register int i;
  register site *st;
  su3_vector pp;
  complex cc;
  Real finalrsq;
  msg_tag *tag;
  int d[4];

  mass_x2 = 2.*mass;
  cgn=0;

  /* Single site, single trace term, computed at the origin */
  *bb0s11 = 0;

  for(icol=0; icol<3; icol++)
    {
      /* First calculate propagator from point source at 0,0,0,1 parallel
	 transported back to 0,0,0,0 */

      /* initialize phi and xxx */
      clear_latvec( F_OFFSET(phi), EVENANDODD);
      clear_latvec( F_OFFSET(xxx), EVENANDODD);
      if( node_number(0,0,0,0) == mynode() )
	{
	  i=node_index(0,0,0,0);
	  for(jcol=0;jcol<3;jcol++)
	    /* Parallel transport unit color vector icol from 0,0,0,1*/
	    lattice[i].phi.c[jcol] = lattice[i].link[TUP].e[jcol][icol];
	}

      g_sync();

      /* do a C.G. (source in phi, result in xxx) */
      /* Preconditioning in this case is trivial and not needed,
	 since the source is only on one site.  
	 However, by not preconditioning as elsewhere,
	 we are giving the near zero modes more relative weight
	 compared with the other propagators in this module.
	 So we divide by the square of the quark mass 
	 in both congrad calls in this procedure to
	 give a comparably stringent stopping criterion 
	 throughout */
      load_ferm_links(&fn_links);
      this_cgn = ks_congrad(F_OFFSET(phi),F_OFFSET(xxx),mass,
			    niter, rsqprop/(mass*mass), PRECISION, 
			    EVEN,&finalrsq, &fn_links);
      cgn += this_cgn;
      if(this_node==0)printf("bar_corr: ss1 cgn %d\n",this_cgn);
      /* Multiply by -Madjoint */
      /* Notice that this product includes an overall minus sign */
      /* which will be compensated by propmat2 */
      load_ferm_links(&fn_links);
      dslash_site( F_OFFSET(xxx), F_OFFSET(ttt), ODD, &fn_links);
      scalar_mult_latvec( F_OFFSET(xxx), -mass_x2, F_OFFSET(ttt), EVEN);
      
      /* Result goes into propmat */
      copy_latvec( F_OFFSET(ttt), F_OFFSET(propmat[icol]), EVENANDODD);

      /* Compute contribution to single-site single trace term */
      /* This term measures only at the origin.  */
      if(node_number(0,0,0,1) == mynode())
	{
	  i = node_index(0,0,0,1);
	  *bb0s11 += lattice[i].propmat[icol].c[icol].real/2.;
	}

      /* Next calculate propagator from point source at 0,0,0,1 */

      /* reinitialize phi and xxx */
      clear_latvec( F_OFFSET(phi), EVENANDODD);
      clear_latvec( F_OFFSET(xxx), EVENANDODD);
      if( node_number(0,0,0,1) == mynode() )
	{
	  i=node_index(0,0,0,1);
	  lattice[i].phi.c[icol].real = 1.;
	}

      g_sync();

      /* do a C.G. (source in phi, result in xxx) now ODD source*/
      load_ferm_links(&fn_links);
      this_cgn = ks_congrad(F_OFFSET(phi),F_OFFSET(xxx),mass,
			    niter, rsqprop/(mass*mass), PRECISION, 
			    ODD, &finalrsq, &fn_links);
      cgn += this_cgn;
      if(this_node==0)printf("bar_corr: ss1: cgn %d\n",this_cgn);
      /* Multiply by -Madjoint */
      /* Again this product includes an overall minus sign which 
	 will cancel the sign in propmat */
      dslash_site( F_OFFSET(xxx), F_OFFSET(ttt), EVEN, &fn_links);
      scalar_mult_latvec( F_OFFSET(xxx), -mass_x2, F_OFFSET(ttt), ODD);
      
      /* This time result goes into propmat2 */
      copy_latvec( F_OFFSET(ttt), F_OFFSET(propmat2[icol]), EVENANDODD);
      
    } /* end loop on icol */
  
  /* Parallel transport propmat2 forward one time link */
  
  FORALLSITES(i,st)
    {
      /* Multiply propmat2 by adjoint of gauge link */
      /* Put product in tempvec[0-2] */
      for(icol=0;icol<3;icol++)
	{
	  mult_adj_su3_mat_vec(&lattice[i].link[TUP],
			       &lattice[i].propmat2[icol], &pp);
	  su3vec_copy( &pp, &lattice[i].tempvec[icol]);
	}
    }
      
  tag = start_gather_site( F_OFFSET(tempvec[0]), 3*sizeof(su3_vector), 
		     TDOWN, EVENANDODD, gen_pt[0] );
  wait_gather(tag);

  FORALLSITES(i,st)
    {
      lattice[i].BBS21 = cmplx(0.,0.);

      /* Combine with propmat */
      
      for(icol=0;icol<3;icol++)
	{
	  cc = su3_dot( (su3_vector *)gen_pt[0][i] + icol, 
		       &lattice[i].propmat[icol] );
	  /* propmat2 time reversal for propagation from r=0,t=1 to r=r,t=t */
	  if((st->x + st->y + st->z + st->t) %2 == 0)
	    lattice[i].BBS21.real -= cc.real;
	  else lattice[i].BBS21.real += cc.real;
	}
    }

  cleanup_gather(tag);
  g_sync();

  /* Parallel transport propmat forward one time link */
  
  FORALLSITES(i,st)
    {

      /* Multiply propmat by adjoint of gauge link */
      /* Put product in tempvec[0-2] */
      for(icol=0;icol<3;icol++)
	{
	  mult_adj_su3_mat_vec(&lattice[i].link[TUP],
			       &lattice[i].propmat[icol], &pp);
	  su3vec_copy( &pp, &lattice[i].tempvec[icol]);
	}
    }
      
  
  tag = start_gather_site( F_OFFSET(tempvec[0]), 3*sizeof(su3_vector), 
		     TDOWN, EVENANDODD, gen_pt[0] );
  wait_gather(tag);
  
  FORALLSITES(i,st)
    {
      /* Combine with propmat2 */
      
      for(icol=0;icol<3;icol++)
	{
	  cc = su3_dot( (su3_vector *)gen_pt[0][i] + icol, 
		       &lattice[i].propmat2[icol] );
	  /* propmat2 time reversal for propagation from r=0,t=1 to r=r,t=t */
	  if((st->x + st->y + st->z + st->t) %2 == 0)
	    lattice[i].BBS21.real += cc.real;
	  else lattice[i].BBS21.real -= cc.real;
	}
      /* Normalization */
      lattice[i].BBS21.real /= 2.;
    }

  cleanup_gather(tag);
  g_sync();

  /* Average line results in BBS21 and put in LLS21 */

  line_avg(F_OFFSET(BBS21),F_OFFSET(LLS21),F_OFFSET(TEMP1));

  return cgn;
}
/************ bar_corr *********************************************/
int bar_corr()
{
  register int i;
  int parity;
  register site *st;
  Real mass_x2;
  Real finalrsq;
  register int icol,j,cgn,this_cgn;
  complex plp;
  int x,y,z,t;
  Real vol,vol2;
  msg_tag* tag;
  complex* pt;
  complex bb, cc;
  int key[4],restrict[4];
  Real bbs11,bb0s11;
  complex bdensum;
  complex lls22,lls22cross;
  Real corrf,corrf0,corrf02,corrf0avg,corrf0var;
  Real fact;
  complex err;

  int irand,nrand,jrand,nbrand;
  static int setup = 0;

  /* These values MUST BE 0 and 1 */
  t0 = 0;
  t0p1 = 1;
  
  /* Do setup only once for any run */
  if(setup==0)
    {
      setup = 1;

      set_slice_swap();
      set_sl_sym_gath();

      key[XUP] = 1;
      key[YUP] = 1;
      key[ZUP] = 1;
      key[TUP] = 2;
      restrict[XUP] = 0;
      restrict[YUP] = 0;
      restrict[ZUP] = 0;
      restrict[TUP] = t0p1;

      setup_restrict_fourier(key,restrict);
    }

  /* Clear averages */
  FORALLSITES(i,st)
    {
      for(j=0;j<NBPAVRG;j++)st->avg[j] = cmplx(99.,99.);
    }

  /* Polyakov loop and fuzzy loops are computed in ploop and ploop_staple 
     and results are stored in even sites on slices 0 and 1 */

  /* Debug */
  if(node_number(0,0,0,t0) == mynode())
    {
      i = node_index(0,0,0,t0);
      printf("PLP 0 0 0 %f %f\n",
	     lattice[i].ploop_fuzz.real,lattice[i].ploop_fuzz.imag);
      fflush(stdout);
    }

  /* Move result in ploop up to time slice t0p1*/
  
  tag = start_gather_site( F_OFFSET(ploop), sizeof(complex), 
		     sl_swap_dir, ODD, gen_pt[0]);
  wait_gather(tag);

  /* Store gathered real and imaginary parts separately on odd sites */

  FORODDSITES(i,st)
    {
      if( st->t != t0p1)continue;
      st->PLPRE.real = ((complex *)gen_pt[0][i])->real;
      st->PLPRE.imag = 0.;
      st->PLPIM.real = ((complex *)gen_pt[0][i])->imag;
      st->PLPIM.imag = 0.;
    }

  cleanup_gather(tag);


  /* Move result in ploop up to time slice t0p1*/
  
  tag = start_gather_site( F_OFFSET(ploop_fuzz), sizeof(complex), 
		     sl_swap_dir, ODD, gen_pt[0]);
  
  wait_gather(tag);

  /* Store gathered real and imaginary parts separately on odd sites */

  FORODDSITES(i,st)
    {
      if( st->t != t0p1)continue;
      st->FLPRE.real = ((complex *)gen_pt[0][i])->real;
      st->FLPRE.imag = 0.;
      st->FLPIM.real = ((complex *)gen_pt[0][i])->imag;
      st->FLPIM.imag = 0.;
    }

  /* Separate real and imag parts on even sites for Fourier transform */

  FOREVENSITES(i,st)
    {
      if( st->t != t0p1)continue;
      st->PLPRE = cmplx(st->ploop.real,0.);
      st->PLPIM = cmplx(st->ploop.imag,0.);
      st->FLPRE = cmplx(st->ploop_fuzz.real,0.);
      st->FLPIM = cmplx(st->ploop_fuzz.imag,0.);
    }

  cleanup_gather(tag);
  
  /* Do Fourier transform of plp values */
  /* Uses TEMP1 and TEMP2 for working space */

  restrict_fourier_site( F_OFFSET(PLPRE), sizeof(complex), FORWARDS);

  restrict_fourier_site( F_OFFSET(PLPIM), sizeof(complex), FORWARDS);
	  
  restrict_fourier_site( F_OFFSET(FLPRE), sizeof(complex), FORWARDS);

  restrict_fourier_site( F_OFFSET(FLPIM), sizeof(complex), FORWARDS);
	  
  /* Compute ploop - ploop correlations */

  FORALLSITES(i,st)
    {
      /* All results are now on time slice t0p1 */
      if( st->t != t0p1)continue;
      CMULJ_(st->PLPRE, st->PLPRE, st->RRAVG);
      CMULJ_(st->PLPIM, st->PLPIM, st->IIAVG);
      CMULJ_(st->FLPRE, st->FLPRE, st->FFAVG);
      CMULJ_(st->FLPIM, st->FLPIM, st->GGAVG);
      /* Normalization of FT's */
      vol = nx*ny*nz;
      vol2 = vol*vol;
      CMULREAL(st->RRAVG,1./vol2,st->RRAVG);
      CMULREAL(st->IIAVG,1./vol2,st->IIAVG);
      CMULREAL(st->FFAVG,1./vol2,st->FFAVG);
      CMULREAL(st->GGAVG,1./vol2,st->GGAVG);
    }

  mass_x2 = 2.*mass;
  cgn=0;

  /* Compute psi-bar-psi (single site definition) */

  /* Clear sums of correlations */
  FORALLSITES(i,st)
    {
      st->SPAVG = cmplx(0.,0.);
      st->SFAVG = cmplx(0.,0.);
      st->SSAVG = cmplx(0.,0.);
    }

  for(irand=0; irand<QRAND; irand++)
    {
      /* Construct random source on full lattice */
      grfull();
      
      /* Compute quark propagator from this source */

      /* initialize ttt, xxx, and phi */
      clear_latvec( F_OFFSET(ttt), EVENANDODD);
      clear_latvec( F_OFFSET(xxx), EVENANDODD);
      FORALLSITES(i,st)
	{
	  st->phi = st->g_rand;
	}

      /* Compute quark propagator from random source with LU preconditioning */
      /* Preconditioning step */
      /* -phi_e' <- -2ma phi_e + Dslash_eo phi_o */
      dslash_site( F_OFFSET(phi), F_OFFSET(ttt), EVEN, &fn_links);
      scalar_mult_add_latvec( F_OFFSET(ttt), F_OFFSET(phi), 
			     -mass_x2, F_OFFSET(phi), EVEN);
      /* Invert on even sites only */
      /* -x_e = -(4m^2a^2 - D_eo D_oe)^(-1) phi_e' */
      load_ferm_links(&fn_links);
      this_cgn = ks_congrad(F_OFFSET(phi),F_OFFSET(xxx),mass,
			    niter, rsqprop, PRECISION, EVEN, &finalrsq, &fn_links);
      cgn += this_cgn;
      if(this_node==0)printf("bar_corr: pbp irand %2d cgn %d\n",irand,this_cgn);
      /* Even sites are now OK, except for a minus sign */
      /* Multiply by 1/2ma L for odd sites */
      /* -x_o <- [Dslash_eo (-x_e) + phi_o]/2ma */
      dslash_site( F_OFFSET(xxx), F_OFFSET(ttt), ODD, &fn_links);
      FORODDSITES(i,st){
	add_su3_vector( &(st->ttt), &(st->phi), &(st->xxx));
	scalar_mult_su3_vector( &(st->xxx), -1./(mass_x2), &(st->xxx));
      }
      /* The resulting propagator includes an overall minus sign */
      /* This sign affects the correlation SPCORR which 
	 must be taken into account in the post analysis*/
      
      /* Result of inversion from random source appears in xxx */

      /* Project back upon random source to get psi-bar-psi */
      FORALLSITES(i,st)
	{
	  st->PBP = su3_dot(&(st->g_rand), &(st->xxx));
	}


      /* Average results and put on time slice t0p1 */
      line_avg(F_OFFSET(PBP),F_OFFSET(PBP),F_OFFSET(TEMP1));
      

      /* Debug */
      if(node_number(0,0,0,t0p1) == mynode())
	{
	  i = node_index(0,0,0,t0p1);
	  printf("PBP0 %.8e\n",-lattice[i].PBP.real/(nx*ny*nz));
	  fflush(stdout);
	}

      /* Do Fourier transform on time slice t0p1 */
      restrict_fourier_site( F_OFFSET(PBP), sizeof(complex), FORWARDS);

      /* Debug */
      if(node_number(0,0,0,t0p1) == mynode())
	{
	  i = node_index(0,0,0,t0p1);
	  printf("PBP %.8e QRAND %d\n",
		 -lattice[i].PBP.real/(nx*ny*nz),QRAND);
	  fflush(stdout);
	}

      /* Compute pbp correlations with Polyakov loop and pbp itself */
      
      FORALLSITES(i,st)
	{
	  if( st->t != t0p1)continue;
	  CMULJ_(st->PLPRE, st->PBP, bb);
	  CSUM(st->SPAVG, bb);
	  CMULJ_(st->FLPRE, st->PBP, bb);
	  CSUM(st->SFAVG, bb);
	  CMULJ_(st->PBP, st->PBP, bb);
	  CSUM(st->SSAVG, bb);
	}
    }

  /* Normalization of correlations */
  vol = nx*ny*nz;
  vol2 = vol*vol;

  FORALLSITES(i,st)
    {
      CMULREAL(st->SPAVG,1./(vol2*QRAND),st->SPAVG);
      CMULREAL(st->SFAVG,1./(vol2*QRAND),st->SFAVG);
      CMULREAL(st->SSAVG,1./(vol2*QRAND),st->SSAVG);
    }

  /* Clear sum for 1 trace - single site contribution */
  bbs11 = 0.;

  /* Computation of baryon number density and related correlations */

  FORALLSITES(i,st)
    {
      st->bdensum = cmplx(0.,0.);
      st->LLS22 = cmplx(0.,0.);
    }

  /* For testing for convergence */
  corrf0 = corrf02 = 0.;

  /* Loop over up to NRANDMAX random sources */
  /* We deal out NRAND sources, do the FT and accumulate the result */
  /* If the variance of the integrated Fuzzy loop/ density correlation is not */
  /* good enough, we deal out another NRAND sources, and keep going until the */
  /* variance is acceptable or we reach NRANDMAX, whichever is sooner */


  for(irand=0; irand<NRANDMAX; irand++)
    {
      nbrand = irand % NRAND;  /* Where we store the result of this round */

      grfull();

      /* Parallel transport source back to time slice t */

      tag = start_gather_site( F_OFFSET(g_rand), sizeof(su3_vector), 
			 TUP, EVENANDODD, gen_pt[0] );
      /* Compute quark propagator from random source on time slice t */
      
      /* initialize ttt and xxx */
      clear_latvec( F_OFFSET(ttt), EVENANDODD);
      clear_latvec( F_OFFSET(xxx), EVENANDODD);
      clear_latvec( F_OFFSET(phi), EVENANDODD);
      
      wait_gather(tag);

      /* Multiply by the forward link and put result in phi */

      for(x=0; x<nx; x++)for(y=0; y<ny; y++)for(z=0; z<nz; z++)for(t=0; t<nt; t++)
	if( node_number(x,y,z,t) == mynode() )
	  { 
	    i=node_index(x,y,z,t);
	    mult_su3_mat_vec(&lattice[i].link[TUP],
			     (su3_vector *)gen_pt[0][i], &lattice[i].phi);
	  }

      cleanup_gather(tag);
      
      g_sync();
      

      /* Compute quark propagator from random source with LU preconditioning */
      /* Preconditioning step */
      /* -phi_e' <- -2ma phi_e + Dslash_eo phi_o */
      dslash_site( F_OFFSET(phi), F_OFFSET(ttt), EVEN, &fn_links);
      scalar_mult_add_latvec( F_OFFSET(ttt), F_OFFSET(phi), 
			     -mass_x2, F_OFFSET(phi), EVEN);
      /* Invert on even sites only */
      /* -x_e = -(4m^2a^2 - D_eo D_oe)^(-1) phi_e' */
      load_ferm_links(&fn_links);
      this_cgn = ks_congrad(F_OFFSET(phi),F_OFFSET(xxx),mass,
			    niter, rsqprop, PRECISION, EVEN, &finalrsq,
			    &fn_links);
      cgn += this_cgn;
      if(this_node==0)printf("bar_corr: bdens irand %2d cgn %d\n",irand,this_cgn);
      /* Even sites are now OK, except for a minus sign */
      /* Multiply by 1/2ma L for odd sites */
      /* -x_o <- [Dslash_eo (-x_e) + phi_o]/2ma */
      dslash_site( F_OFFSET(xxx), F_OFFSET(ttt), ODD, &fn_links);
      FORODDSITES(i,st){
	add_su3_vector( &(st->ttt), &(st->phi), &(st->xxx));
	scalar_mult_su3_vector( &(st->xxx), -1./(mass_x2), &(st->xxx));
      }
      /* The resulting propagator includes an overall minus sign */
      /* This sign affects both the S11 term and the ploop correlation term which 
	 must be taken into account in the post analysis*/
      
      /* Result of inversion from random source appears in xxx */

      /* Project result back onto random source */
      /* Stash result in bdens */

      for(x=0; x<nx; x++)for(y=0; y<ny; y++)for(z=0; z<nz; z++)for(t=0; t<nt; t++)
	if(node_number(x,y,z,t) == mynode())
	  {
	    i = node_index(x,y,z,t);
	    cc = su3_dot(&lattice[i].g_rand, &lattice[i].xxx);
	    lattice[i].bdens[nbrand].real =  -cc.imag/2.;
	    lattice[i].bdens[nbrand].imag =  0.;
/*	      printf("DEB %d %d %d %d %d %e\n",nbrand,x,y,z,t,-cc.imag/2.); */

	      bbs11 += cc.real/2.;
	    }
       /* Average results on line and put in bdens on time slice t0p1 */
      
      line_avg(F_OFFSET(bdens[nbrand]),F_OFFSET(bdens[nbrand]),F_OFFSET(TEMP1));
      
      FORALLSITES(i,st)
	{
	  if((st->t != 1) || (st->x != 0) || (st->y != 0) || (st->z != 0))continue;
	  printf("DEBL %d %d %d %d %e %e\n",
		 st->x,st->y,st->z,irand,
		 st->bdens[nbrand].real,st->bdens[nbrand].imag);
	  fflush(stdout);
	}
      
      /* Do Fourier transform of baryon density estimator */
      /* Uses TEMP1 for working space */
      /* If conserving memory were not an object, then we could wait
	 and do the FT on the entire set at once outside the nbrand loop */
      
      restrict_fourier_site(F_OFFSET(bdens[nbrand]), sizeof(complex), FORWARDS);

      nrand = irand + 1;
      if((nrand % NRAND) == 0)
	{
	  
	  /* Calculate correlations and accumulate results for this round */
	  corrf0var = 0.;

 	  FORALLSITES(i,st)
	    {
	      /* All results are now on time slice t0p1 */
	      if( st->t != t0p1)continue;
	      
	      /* Total up baryon density results for all random sources */
	      
	      bdensum = cmplx(0.,0.);

	      for(jrand = 0; jrand < NRAND; jrand++)
		{
		  CSUM(bdensum,st->bdens[jrand]);

		  /* Diagnostic */
		  if((st->x == 0) && (st->y == 0) &&(st->z == 0) && 
		     (node_number(0,0,0,st->t) == mynode()))
		    {
		      corrf = st->FLPIM.real*st->bdens[jrand].real;
		      corrf0 += corrf;
		      corrf02 += corrf*corrf;
		      printf("DIAG0 %d %f\n",jrand,corrf);
		    }
		}

	      /* Compute and report current variance */
	      if((st->x == 0) && (st->y == 0) &&(st->z == 0) &&
		 (node_number(0,0,0,st->t) == mynode()))
		{
		  corrf0avg = corrf0/nrand;
		  corrf0var = (corrf02/nrand) - 
		    (corrf0/nrand)*(corrf0/nrand);
		  corrf0var = corrf0var/nrand;
		  
		  printf("BPVAR %d %f %f %f\n",
			 nrand,corrf0avg,sqrt(corrf0var),MAXVAR);
		}

	      /* Compute density-density correlation */

	      CMULJ_(bdensum, bdensum, lls22);
	      /* Remove unwanted diagonal terms from lls22:*/
	      for(jrand = 0; jrand < NRAND; jrand++)
		{
		  CMULJ_(st->bdens[jrand], st->bdens[jrand], bb);
		  CSUB(lls22, bb, lls22);
		}
	      /* Compute cross term from previous blocks.
	         These are of the form 2*Re(bi^* * bj) 
		 where bi is the baryon density sum from the ith block 
		 of NRAND random sources */
	      /* At this point st->bdensum contains the cumulative
		 sum of bdensum from previous blocks */
	      CMULJ_(bdensum,st->bdensum,lls22cross);
	      
	      /* Cumulative baryon density sum is kept in the site structure */
	      /* local bdensum is only for this block */
	      CSUM(st->bdensum,bdensum);
	      /* LFAVG and LPAVG are cumulative over the entire sample */
	      CMULJ_(st->bdensum, st->FLPIM, st->LFAVG);
	      CMULJ_(st->bdensum, st->PLPIM, st->LPAVG);

	      vol = nx*ny*nz;
	      vol2 = vol*vol;
	      /* Normalization of FT's */
	      CMULREAL(st->LPAVG,1./(nrand*vol2),st->LPAVG);
	      CMULREAL(st->LFAVG,1./(nrand*vol2),st->LFAVG);
	      /* Undo previous normalization of LLS22 */
	      CMULREAL(st->LLS22,(Real)(nrand-NRAND)*(nrand-NRAND-1),
		       st->LLS22);
	      /* Volume factors for FFT and our sign convention */
	      CMULREAL(lls22,-1./vol2,lls22);
	      CMULREAL(lls22cross,-1./vol2,lls22cross);
	      CSUM(st->LLS22,lls22);
	      st->LLS22.real += 2.*lls22cross.real;
	      /* New normalization of LLS22 */
	      CMULREAL(st->LLS22,1./(nrand*(nrand-1)),st->LLS22);
	    }

	  /* Decide whether to quit now or not */
	  /* Only one node should have the nonzero value of corrf0var */
	  g_floatsum( &corrf0var );
	  if(corrf0var < MAXVAR)break;
	}	  
    } /* irand loop */
  
  /* Do inverse Fourier transform of all NBPAVRG correlations */
  /* Note that two of the terms in this FT are irrelevant,
     since they have yet to be computed in "sngl_trace" below */
  

  restrict_fourier_site(F_OFFSET(avg[0]), NBPAVRG*sizeof(complex), BACKWARDS);
  
  g_sync();

  /* Calculate single-color-trace density-density correlation (point-to-point
     hadron propagator) */

  cgn += sngl_trace(&bb0s11);

  /* Average 1 trace, single site contributions */
  /* bbs11 is mesured as a space-time average and an average over random sources.  */
  /* bb0s11 is measured only at origin. */

  g_floatsum( &bbs11 );
  g_floatsum( &bb0s11);
  if(mynode()==0)printf("BCOR0 %.8e %.8e\n",
			bbs11/(nrand*vol*nt),bb0s11);
  
  /* Combine results at displacements that are related by symmetry transforms */
  
  if(nx==ny) 
    {
      if(nx==nz)
	{
	  /* nx = ny = nz:  Fetch but do not combine*/
	  slice_sym_combine(xy_dir,F_OFFSET(gath1[0]),0);
	  slice_sym_combine(xz_dir,F_OFFSET(gath2[0]),0);
	  
	  /* Combine and normalize four values */
	  FORALLSITES(i,st)
	    {
	      if( st->t != t0p1)continue;
	      for (j = 0; j < NBPAVRG; j++)
		{
		  CSUM( st->avg[j], st->gath1[j]);
		  CSUM( st->avg[j], st->gath2[j]);
		  CMULREAL( st->avg[j], 0.3333333333333, st->avg[j]);
		}
	    }
	  
	}
      /* nx = ny != nz Combine only x<->y permutation */
      else slice_sym_combine(xy_dir,F_OFFSET(gath1[0]),1);
    }
  /* nx = nz != ny Combine only x<->z permutation*/
  else if(nx==nz) slice_sym_combine(xz_dir,F_OFFSET(gath1[0]),1);
  
  /* Now combine with y <-> z permutation */
  if(ny==nz) slice_sym_combine(yz_dir,F_OFFSET(gath1[0]),1); 
  
  /* Finally, combine with reflections */
  for(i=XUP;i<TUP;i++)
    {
      slice_sym_combine(minus_dir[i],F_OFFSET(gath1[0]),1);
    }
  
  /* Report results */

  write_corr();

  return cgn;

} /* bar_corr11.c */
