/********************* density.c ****************************/
/* MIMD version 6 */
/*			    -*- Mode: C -*-
// File: density.c
// Created: Thu Jun 15 1995
// Author: J. E. Hetrick <hetrick@corelli.physics.Arizona.EDU>
//
// Description: measure singlet and non-singlet SU(2) densities at
//              finite chemical potential.
//
//              Uses ploop_t, "P-loop minus time slice", which is defined
//              in metropolis_dense.c. Thus a call to either update_dense(),
//              monte_time(), or explicitly ploop_less_slice() must be made
//              so that s->ploop_t is loaded. Puts Tr P into s->ploop on 
//              time slice t=0 en passent.
// 7/21/95 Add probabilities for 0,1,2,3 quarks at a site
*/

#include "su3_dense_includes.h"

#define CMULINC(a,b,c) { (c).real += (a).real*(b).real - (a).imag*(b).imag; \
                      (c).imag += (a).real*(b).imag + (a).imag*(b).real; }

#define CMULDEC(a,b,c) { (c).real -= (a).real*(b).real - (a).imag*(b).imag; \
                      (c).imag -= (a).real*(b).imag + (a).imag*(b).real; }



complex numer(Real c, complex T, complex cof) {
   complex z1, z2;

   CMULREAL(T, c*c, z1);
   CMULREAL(cof, 2.*c, z2);
   CSUM(z1,z2);
   z1.real += 3.;

   return(z1);
}


/* do measurements: load density, ploop, etc. and phases onto lattice */
void measure() {
   register int i,j,k, c;
   register site *s;
   int dx,dy,dz;	/* separation for correlated observables */
   int dir;		/* direction of separation */
   msg_tag *tag;
   register complex cc,dd;	/*scratch*/
   complex ztr, zcof, znum, zdet, TC, zd, density, zphase;
   complex p[4]; /* probabilities of n quarks at a site */
   complex np[4]; /* probabilities at neighbor site */
   complex pp[4][4]; /* joint probabilities of n here and m there */
   complex zplp, plp;
   Real locphase, phase;


   /* First make T (= timelike P-loop) from s->ploop_t 
      T stored in s->tempmat1
   */
   ploop_less_slice(nt-1,EVEN);
   ploop_less_slice(nt-1,ODD);

   phase = 0.;
   density = plp = cmplx(0.0, 0.0);
   for(j=0;j<4;j++){
	p[j]=cmplx(0.0,0.0);
	for(k=0;k<4;k++)pp[j][k]=cmplx(0.0,0.0);
   }
   FORALLSITES(i,s) {
      if(s->t != nt-1) continue;
      mult_su3_nn(&(s->link[TUP]), &(s->ploop_t), &(s->tempmat1));

      zplp = trace_su3(&(s->tempmat1));
      CSUM(plp, zplp);

      ztr = trace_su3(&(s->tempmat1));
      CONJG(ztr, zcof);

      for(c=0; c<3; ++c) s->tempmat1.e[c][c].real += C;
      zdet = det_su3(&(s->tempmat1));
      znum = numer(C, ztr, zcof);
      CDIV(znum, zdet, zd);
      CSUM(density, zd);

      /* store n_quark probabilities at this site in lattice variable
	qprob[], accumulate sum over lattice in p[] */
      cc = cmplx(C*C*C,0.0); CDIV(cc,zdet,s->qprob[0]); CSUM(p[0],s->qprob[0]);
      CMULREAL(ztr,C*C,cc); CDIV(cc,zdet,s->qprob[1]); CSUM(p[1],s->qprob[1]);
      CMULREAL(zcof,C,cc); CDIV(cc,zdet,s->qprob[2]); CSUM(p[2],s->qprob[2]);
      cc = cmplx(1.0,0.0); CDIV(cc,zdet,s->qprob[3]); CSUM(p[3],s->qprob[3]);

      locphase = carg(&zdet);
      phase += locphase;

   }
   g_floatsum( &phase );
   g_complexsum( &density );
   g_complexsum( &p[0] );
   g_complexsum( &p[1] );
   g_complexsum( &p[2] );
   g_complexsum( &p[3] );
   g_complexsum( &plp );
   CDIVREAL(density,(Real)(nx*ny*nz),density);
   CDIVREAL(p[0],(Real)(nx*ny*nz),p[0]);
   CDIVREAL(p[1],(Real)(nx*ny*nz),p[1]);
   CDIVREAL(p[2],(Real)(nx*ny*nz),p[2]);
   CDIVREAL(p[3],(Real)(nx*ny*nz),p[3]);
   CDIVREAL(plp,(Real)(nx*ny*nz),plp);

   zphase = ce_itheta(phase);
   if(this_node == 0) {
      printf("ZMES\t%e\t%e\t%e\t%e\t%e\t%e\n", zphase.real, zphase.imag, 
	                               density.real, density.imag,
	                               plp.real, plp.imag);
      printf("PMES\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n",
				p[0].real, p[0].imag, p[1].real, p[1].imag,
				p[2].real, p[2].imag, p[3].real, p[3].imag );
   }

#ifdef PPCORR
   dx=1; dy=0; dz=0;	/* Temporary - right now we just do nearest neighbor */
   for(dir=XUP;dir<=ZUP;dir++){
      tag = start_gather_site( F_OFFSET(qprob[0]), 4*sizeof(complex), dir,
	   EVENANDODD, gen_pt[0] );
      wait_gather(tag);
      FORALLSITES(i,s)if(s->t==nt-1){
        for(j=0;j<4;j++)for(k=0;k<4;k++){
	   CMUL( (s->qprob)[j],((complex *)gen_pt[0][i])[k],cc);
           CSUM(pp[j][k],cc);
        }
      }
      cleanup_gather(tag);
   }

   /* density correlation format:
	PP dx dy dz n1 n2 real imag */
   for(j=0;j<4;j++)for(k=0;k<4;k++){
     g_complexsum( &pp[j][k] );
     CDIVREAL(pp[j][k],(Real)(3*nx*ny*nz),pp[j][k]);
     if(this_node==0)
       printf("PP %d %d %d   %d %d   %e   %e\n",dx,dy,dz,j,k,
	  pp[j][k].real,pp[j][k].imag);
   }
#endif /*PPCORR*/
}
