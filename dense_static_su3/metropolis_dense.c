/**************************** metropolis_dense.c ****************************/
/* MIMD version 6 */

/************************** metropolis_dense.c *******************************/
/* Metropolis updating for SU3 pure gauge */
/* MIMD version 4 */
/* J. Hetrick and D. Toussaint June 1995 */
/* update with "almost quenched" approximation for nonzero density 6/1/95 */
/* monte_space() does spatial links, monte_time() does temporal */


#include "su3_dense_includes.h"

/* Generic definitions - could be useful elsewhere */
#define FORSPACEUPDIR(dir) for(dir=XUP; dir<=ZUP; dir++)
#define OPP_PAR(parity) (0x03 ^ parity) /* Switches EVEN and ODD. 
                                           Nulls EVENANDODD*/
#define CMULINC(a,b,c) { (c).real += (a).real*(b).real - (a).imag*(b).imag; \
                      (c).imag += (a).real*(b).imag + (a).imag*(b).real; }
#define CMULDEC(a,b,c) { (c).real -= (a).real*(b).real - (a).imag*(b).imag; \
                      (c).imag -= (a).real*(b).imag + (a).imag*(b).real; }

Real make_change(su3_matrix *mat, site *st, Real scale);
void base_measurements(int time);
complex numerator(Real c, complex T);
complex det_special(Real c, complex T);

/* det. phase, P. loop, and total density before this processor starts
   updating this parity and slice */
Real base_phase;
complex base_ploop, base_density;

void monte_space(int NumStp) {
  int Nhit;
  int parity;
  Real scale;		/* limits size of change matrix */
  Real theta,theta_rms;	/* rms angle of change matrix */
  su3_matrix change;	/* proposed change in link */
  su3_matrix newlink;	/* change * oldlink */
  su3_matrix tmat;
  int dir, i;
  register site *st;
  register int c;
  int accept, reject;	/* number of accepts and rejects */
  Real oldaction, newaction;

  accept = reject = 0;
  theta_rms = 0.0;
  scale = 1.1;
  for(parity=ODD;parity<=EVEN;parity++) {
    FORSPACEUPDIR(dir) {
      /* compute the gauge force */
      dsdu_qhb(dir,parity); 
      /* now for the Metropolis updating */
      FORSOMEPARITY(i,st,parity) {
	/* generate random SU(3) matrix */
	/* scale < 2/sqrt(3), so vector magnitude < 1 */

        for( Nhit = 0 ; Nhit < NumStp; Nhit++) { 
	  theta = make_change(&change,st,scale); theta_rms += theta*theta;
	  mult_su3_nn( &change, &(st->link[dir]), &newlink );
	  
	  /* compute old action and new action */
	  oldaction=(0.333333*beta)*realtrace_su3( &(st->link[dir]), 
						       &(st->staple) );
	  newaction=(0.333333*beta)*realtrace_su3( &newlink, 
						       &(st->staple) );
	  /* accept or reject */
	  if( newaction > oldaction ){
	    st->link[dir]=newlink;
	    accept++;
	  }
	  else{ /* trace decreased */
	    if( myrand(&(st->site_prn)) < exp( newaction-oldaction ) ){
	      st->link[dir]=newlink;
	      accept++;
	    }
	    else{
	      reject++;
	    }
	  }
        } /* Nhit */
      } /*   st */
    } /*  direction */
  } /* parity */

  /* diagnostics */
/**/
  printf("monte_space: accept= %d reject= %d fraction= %.2f  theta_rms= %.3f\n",
    accept,reject,accept/(Real)(accept+reject),
    sqrt(theta_rms/(accept+reject)) );/**/

} /* monte_space */


/* time direction update includes det(Polyakov+C) in its action */

void monte_time(int NumStp) {
  int Nhit;
  int time,parity;
  Real scale;		/* limits size of change matrix */
  Real theta,theta_rms;	/* rms angle of change matrix */
  su3_matrix change;	/* proposed change in link */
  su3_matrix newlink;	/* change * oldlink */
  su3_matrix tmat;	/* scratch */
  register int i;
  register site *st;
  register int c;
  int accept, reject;	/* number of accepts and rejects */
  int mcount;		/* number of measurements */
  Real oldaction,newaction;
  Real oldphase,newphase,phase;
  complex phasefactor;
  complex oldtrace,newtrace,trace;
  complex oldzdet,newzdet,zdet;
  complex oldznum,newznum,znum;
  complex ploop;
  Real mag_ploop;
  complex cc;
  complex olddensity,newdensity,density;
  complex phase_aver,ploop_aver,mag_ploop_aver,density_aver;

  scale = 1.1;
  mcount = accept = reject = 0;
  theta_rms = 0.0;
  phase_aver = ploop_aver = mag_ploop_aver = density_aver = cmplx(0.0,0.0);

  for(parity=ODD;parity<=EVEN;parity++) {
    /* compute the gauge force */
    dsdu_qhb(TUP,parity); 
    /* update time slices separately, because they are linked by P. loops */
    for(time=0;time<nt;time++){
      /* evaluate "almost Polyakov loops", total P. loop, total phase */
      /* see extended comment at bottom */
      ploop_less_slice(time,parity);
      ploop_less_slice(time,OPP_PAR(parity)); /* temporary? */
      base_measurements(time);
/*printf("base:  %e    %e  %e     %e  %e\n",base_phase,base_ploop.real,
base_ploop.imag, base_density.real, base_density.imag);*/
      phase = base_phase;
      ploop = base_ploop;
      density = base_density;
      phasefactor = ce_itheta(phase);

      /* now for the Metropolis updating */
      FORSOMEPARITY(i,st,parity)if(st->t==time) {
	/* generate random SU(3) matrix */
	for( Nhit = 0 ; Nhit < NumStp; Nhit++) {
	  theta = make_change(&change,st,scale); theta_rms += theta*theta;
	  mult_su3_nn( &change, &(st->link[TUP]), &newlink );

	  /* compute old action and new action */

	  oldaction=(0.333333*beta)*realtrace_su3(&(st->link[TUP]),
		  			    	  &(st->staple));
	  /* finish P. loop,  add C*identity matrix to matrix */
	  mult_su3_nn( &(st->link[TUP]), &(st->ploop_t), &tmat );
	  oldtrace = trace_su3( &tmat );
	  /*for(c=0; c<3; ++c) tmat.e[c][c].real += C;*/
	  oldzdet = det_special(C,oldtrace);
	  oldphase = carg( &oldzdet );
	  oldznum = numerator(C, oldtrace);
	  CDIV(oldznum, oldzdet, olddensity);
	  oldaction += log(cabs(&oldzdet));

	  newaction=(0.333333*beta)*realtrace_su3( &newlink, &(st->staple) );
	  mult_su3_nn( &newlink, &(st->ploop_t), &tmat );
	  newtrace = trace_su3( &tmat );
	  newzdet = det_special(C,newtrace);
	  newphase = carg( &newzdet );
	  newznum = numerator(C, newtrace);
	  CDIV(newznum, newzdet, newdensity);
/*printf("DEN:  %e  %e  %e  %e  %e  %e\n",
newznum.real,newznum.imag,newzdet.real,newzdet.imag,newdensity.real,newdensity.imag);*/
	  newaction += log(cabs(&newzdet));

	  /* accept or reject */
	  if( newaction > oldaction ){
	    st->link[TUP]=newlink;
	    accept++;
	    phase += newphase-oldphase;
	    phasefactor = ce_itheta(phase);
	    CSUM(ploop,newtrace); CSUB(ploop,oldtrace,ploop);
	    CSUM(density,newdensity); CSUB(density,olddensity,density);
	  }
	  else{ /* trace decreased */
	    if( myrand(&(st->site_prn)) < exp( newaction-oldaction ) ){
	      st->link[TUP]=newlink;
	      accept++;
	      phase += newphase-oldphase;
	      phasefactor = ce_itheta(phase);
	      CSUM(ploop,newtrace); CSUB(ploop,oldtrace,ploop);
	      CSUM(density,newdensity); CSUB(density,olddensity,density);
	    }
	    else{
	      reject++;
	    }
	  }

	  /* accumulate measurements */
	  CSUM(phase_aver,phasefactor);
	  CMULINC(ploop,phasefactor,ploop_aver);
	  CMULINC(density,phasefactor,density_aver);
	  mag_ploop = cabs(&ploop);
	  CMULREAL(phasefactor,mag_ploop,cc);
	  CSUM(mag_ploop_aver,cc);
	  mcount++;
/*DEBUG*/
/*printf("MT:  %d  %d  %d  %d, accept = %d\n",st->x,st->y,st->z,st->t,accept);
printf("action:  %.3e   %.3e\n",oldaction,newaction);
printf("phase:  %.3e  %.3e  %.3e\n",oldphase,newphase,phase);
printf("ploop:  %.3e  %.3e    %.3e  %.3e    %.3e  %.3e\n",
oldtrace.real,oldtrace.imag,newtrace.real,newtrace.imag,ploop.real,ploop.imag);
printf("density: %.3e  %.3e    %.3e  %.3e    %.3e  %.3e\n",
olddensity.real,olddensity.imag,newdensity.real,newdensity.imag,density.real,density.imag);
printf("\n");*/

	} /* Nhit */
      } /*   st */

    } /* time slice */
  } /* parity */

  /* normalize and dump averages */
  CDIVREAL(phase_aver,mcount,phase_aver);
  CDIVREAL(ploop_aver,mcount*nx*ny*nz,ploop_aver);
  CDIVREAL(mag_ploop_aver,mcount*nx*ny*nz,mag_ploop_aver);
  CDIVREAL(density_aver,mcount*nx*ny*nz,density_aver);
/**
  if(this_node==0)printf(
	"MTMES:  %d     %e  %e     %e  %e     %e  %e     %e  %e\n",
	mcount,
	phase_aver.real, phase_aver.imag,
	density_aver.real, density_aver.imag,
	ploop_aver.real, ploop_aver.imag,
	mag_ploop_aver.real, mag_ploop_aver.imag );
**/

  /* diagnostics */
/**/
  printf("monte_time:  accept= %d reject= %d fraction= %.2f  theta_rms= %.3f\n",
    accept,reject,accept/(Real)(accept+reject),
    sqrt(theta_rms/(accept+reject)) );
/**/

} /* monte_time */



/* Calculate product of timelike links for all time slices except "time",
   for sites where parity at "time" is "parity".
   Put the result in the matrix "ploop_t"
*/
void ploop_less_slice(int time,int parity){
  int l_time;	/* time at which we are multiplying */
  int c_time;   /* l_time%nt - use this one in accessing sites */
  int l_parity; 
	/* parity of sites at l_time where parity at "time" is "parity" */
  register int i;
  register site *s;
  msg_tag *tag;

  if(parity==EVENANDODD){
    if(this_node==0)printf("Bad parity in ploop_less_slice()\n");
    terminate(0);
  }

  FORSOMEPARITY( i,s,OPP_PAR(parity) )if(s->t==(time-1+nt)%nt){
    s->ploop_t = s->link[TUP];
  }

  for( l_time = time + nt-2; l_time > time; l_time--){
    c_time = l_time >= nt ? l_time-nt : l_time;
    if( (l_time-time)%2 ==0){l_parity=parity;}
    else {l_parity = OPP_PAR(parity);}

    /* gather current product from slice above */
    tag = start_gather_site( F_OFFSET(ploop_t), sizeof(su3_matrix), TUP,
      l_parity, gen_pt[0]);
    wait_gather(tag);
    FORSOMEPARITY( i,s,l_parity )if(s->t==c_time){
      mult_su3_nn( &(s->link[TUP]), (su3_matrix *)gen_pt[0][i],
	&(s->ploop_t) );
	/* since only on one time slice, don't have to worry if
	   gen_pt points to ploop_t */
    cleanup_gather(tag);
    } /* end loop over sites at time l_time */
  } /* end loop over l_times */

  /* finish by bringing result to desired time slice */
  tag = start_gather_site( F_OFFSET(ploop_t), sizeof(su3_matrix), TUP,
    parity, gen_pt[0]);
  wait_gather(tag);
  FORSOMEPARITY( i,s,parity )if(s->t==time)
    s->ploop_t = *(su3_matrix *)gen_pt[0][i];
  cleanup_gather(tag);
}


/* Make change matrix for Metropolis step.  Choose random generator
  and rotation angle */
Real make_change(su3_matrix *mat, site *st, Real scale){
  register int ia,ib,gen;
  register Real theta,c,s;

  /* load identity into "mat" */
  for(ia=0; ia<3; ia++) {
    for(ib=0; ib<3; ib++){
      mat->e[ia][ib].real = mat->e[ia][ib].imag = 0.0;
    }
    mat->e[ia][ia].real = 1.0;
  }

  /* random (symmetric) change */
  /* angles range from +- theta/4 to theta/2 */
  do{ theta = scale * (myrand(&(st->site_prn)) - 0.5); }
	while( fabs(theta) < 0.25*scale );
  c = cos(theta); s = sin(theta);
  /* random generator, range 1 through 8 */
  gen = ((int)(8.0*myrand(&(st->site_prn))))+1;

  switch(gen){
    case 1:
      mat->e[0][0].real = mat->e[1][1].real = c;
      mat->e[0][1].imag = mat->e[1][0].imag = s;
      break;
    case 2:
      mat->e[0][0].real = mat->e[1][1].real = c;
      mat->e[0][1].real = s;
      mat->e[1][0].real = -s;
      break;
    case 3:
      mat->e[0][0].real = mat->e[1][1].real = c;
      mat->e[0][0].imag = s;
      mat->e[1][1].imag = -s;
      break;
    case 4:
      mat->e[0][0].real = mat->e[2][2].real = c;
      mat->e[0][2].imag = mat->e[2][0].imag = s;
      break;
    case 5:
      mat->e[0][0].real = mat->e[2][2].real = c;
      mat->e[0][2].real = s;
      mat->e[2][0].real = -s;
      break;
    case 6:
      mat->e[1][1].real = mat->e[2][2].real = c;
      mat->e[1][2].imag = mat->e[2][1].imag = s;
      break;
    case 7:
      mat->e[1][1].real = mat->e[2][2].real = c;
      mat->e[1][2].real = s;
      mat->e[2][1].real = -s;
      break;
    case 8:
      mat->e[0][0].real = mat->e[1][1].real = cos(theta/sqrt(3.));
      mat->e[0][0].imag = mat->e[1][1].imag = sin(theta/sqrt(3.));
      mat->e[2][2].real = cos(-2.*theta/sqrt(3.));
      mat->e[2][2].imag = sin(-2.*theta/sqrt(3.));
      break;
  }
/*printf("make_change: theta = %e\n",theta);*/
  return(theta);
}

/* Numerator in expression for density at a site */
complex numerator(Real c, complex T) {
   complex z0, z1, z2;

   z0.real=T.real; z0.imag= -T.imag; /* z0 = conjg(T) */
   CMULREAL(T, c*c, z1);
   CMULREAL(z0, 2.0*c, z2);
   CSUM(z1,z2);
   z1.real += 3.0;

   return(z1);
}

/* determinant of "unitary + c*unit" matrix, T is trace of matrix */
/* det = C^3 + C^2*T + C*T^* + 1 */
complex det_special(Real c, complex T) {
   complex z0, z1, z2;

   z0.real=T.real; z0.imag= -T.imag; /* z0 = conjg(T) */
   CMULREAL(T, c*c, z1);
   CMULREAL(z0, c, z2);
   CSUM(z1,z2);
   z1.real += 1.0+c*c*c;

   return(z1);
}

/* Measure polyakov loop and density to get started.  Assumes that
   ploop_less_slice(time) has just been called for both parities. */
void base_measurements(int time){
   register int i;
   register site *s;
   complex ztr, zcof, znum, zdet, TC, zd, density, zphase;
   complex zplp, plp;
   Real locphase, phase;

   /* First make T (= timelike P-loop) from s->ploop_t 
      T stored in s->tempmat1
      (assumes ploop_less_slice was recently called...) 
   */

   phase = 0.;
   density = plp = cmplx(0.0, 0.0);
   FORALLSITES(i,s) {
      if(s->t != time) continue;
      mult_su3_nn(&(s->link[TUP]), &(s->ploop_t), &(s->tempmat1));

      ztr = trace_su3(&(s->tempmat1));
      CSUM(plp, ztr);

      zdet = det_special(C,ztr);
      znum = numerator(C, ztr);
      CDIV(znum, zdet, zd);
      CSUM(density, zd);

      locphase = carg(&zdet);
      phase += locphase;

   }
   g_floatsum( &phase );
   g_complexsum( &density );
   g_complexsum( &plp );
   base_phase = phase;
   base_ploop = plp;
   base_density = density;
}

/* EXTENDED COMMENT
  The time direction link updating also averages the phase and a few
observables, since everything must be weighted by the phase factor and
we need really good statistics.
  The logic is a little complicated on parallel machines.  Each node is
doing some subset of the links with the current parity and time slice,
but it doesn't know what the other nodes are doing to their links.  It
doesn't need to for the updating of the link.  However, to weight things
by the total phase one would think that you need to know how the other
nodes are changing their contributions to the phase.  However, we may
imagine that the node we are on is working, while all the other nodes
are turned off.  In this case we can keep track of the phase, Polyakov
loop, etc. just from recording the cumulative changes we make on our
node.  Thus we are generating a legitimate set of configurations.
Now imagine returning to the starting point, and letting the next node
change its links, with the rest of the lattice fixed.  Again, this
is a set of equilibrium configurations.  Repeat this process for
each node, always starting at the same configuration.  
  However, at the end we take as the starting point for the next time
slice updates the configuration obtained by taking the changes made by
all the nodes.  Since we have checkerboarded carefully and do one
time slice at a time, this is also a legal equilibrium configuration.
*/
