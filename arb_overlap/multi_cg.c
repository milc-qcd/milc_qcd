/******* multi_cg.c - conjugate gradient for SU3/fermions ****/
/* MIMD version 7 */
/* Calls hdelta0.c for fermion matrix */

/* version of 22 Feb 00, updated 1 April 2002  */

/* This is an implementation of a multi-mass conjugate gradient
algorithm a la B. Jegerlehner. The number of masses, num_masses,
will be set in lattice.h, and so will the array of mass values,
shift0[num_masses].  */

/* The source vector is in "src", and the answer (initial guess is zero)
   in "dest". 
   "psim[jmass]" are working vectors for the conjugate gradient.
   rsqmin = desired rsq, quit when we reach rsq = rsqmin*source_norm.

In this version of the code, all the scalars are real because M=D^\dagger D

All the inverter stuff usually put in control() is here, just for neatness.


Conventions:
massless case
H=R0(gamma-5 + epsilon(h(-R0)))
H_0^2 = R)^2(2+gamma-5 eps + eps gamma-5)

Massive:
D(m)=(1-m/(2R0))D(0) +m
H(m)^2 = C(H^2+S)
C=(1-m*m/(4R0*R0))) S=m^2/C

Multimass inverter finds psi = (H^2+S)^{-1} chi
Before subtraction
D^{-1} = D(m)^\dagger (H(m))^{-2} = [(1-m/(2R0))H(0) gamma-5 + m] psi/C


Then subtraction for psibar-psi, etc...

Dtilde chi = (D^{-1} chi - chi/(2R0))/(1-m/(2R0))


NOTE: This code uses trial vectors of H(0)^2 in the projection, with
the two chirality states packed in to the upper and lower halves of
eigVec (since eigenvectors of H(0)^2 can be taken to be chiral...)

Eigenmodes treated explicitly:

The  propagator matrix mode by mode, is equal to

   | a -b|
   | b  a|

with mu = mass/(2*R0), eps = lambda/(2*R0)
and lambda the eigenmode of H(0) = sqrt(eigValH0)
a = mu*(1-eps*eps)/(2*R0)/(mu*mu + eps*eps*(1.0-mu*mu))
b = eps*sqrt(1-eps*eps)/(2*R0)/(mu*mu + eps*eps*(1.0-mu*mu))



Chiral modes are
 |1/m 0|   -> |1/m 0|
 |0   0|      |0 1/m|
and we can work with the identity since the wrong-chirality
part of the state is zeroed.

This change was made, so that no multiplications of Dslash onto low
modes are ever performed.


*/


#include "arb_ov_includes.h"

int congrad_xxx(
    field_offset src,   /* type wilson_vector (where source is to be created)*/
    Real cgmass, /* unused here*/
    int source_chirality /* chirality sector for inversion (NOT USED)  */
    )
{
register int i;
register site *s;
int j,k, avs_iters, avm_iters,status,flag;
int MaxCG;
int ksource, spin,color,my_chirality,chb,che,chbo,cheo,ii,jj;
Real *RsdCG;
Real size_r,one_minus_m,r02inv;

wilson_vector **psim;

void setup_multi();

w_prop_file *fp_out_w[MAX_MASSES];       /* For propagator files */
w_prop_file *fp_in_w[MAX_MASSES];        /* For propagator files */
w_prop_file *h0_out_w[MAX_MASSES];       /* For intermediate propagator files */


#ifdef EIGO
wilson_vector wproj;
complex ctmp,cd,*cproj;

int l;
int icount, ivec;
int *chiral_check;
Real cdp, cdm;
Real *ca, *cb;
Real eps, mu, denom;
#endif

double source_norm;

RsdCG=resid;
MaxCG=niter;
avs_iters=0;
r02inv= -0.5/R0;

#ifdef MINN
  do_minn=1;
#endif

    setup_multi();

#ifdef EIGO
  if(Nvecs_hov != 0)cproj = (complex *)malloc(Nvecs_hov*sizeof(complex));
  /* check chirality of your modes (to identify zero modes) */
  if(Nvecs_hov != 0)chiral_check= (int *)malloc(Nvecs_hov*sizeof(int));
  for(j=0;j<Nvecs_hov;j++){
    cdp=0.0;
    cdm=0.0;
    FORALLSITES(i,s){
      for(l=0;l<2;l++)for(k=0;k<3;k++){
        cdp += cabs_sq(&(eigVec[j][i].d[l].c[k]));
      }
      for(l=2;l<4;l++)for(k=0;k<3;k++){
        cdm += cabs_sq(&(eigVec[j][i].d[l].c[k]));
      }
    }
    g_floatsum(&cdp);
    g_floatsum(&cdm);

    if(cdm< 1.e-6 && cdp >1.e-6)
      chiral_check[j] =1;
    else if (cdm >1.e-6 && cdp < 1.e-6)
      chiral_check[j] = -1;
    else if (cdm >1.e-6 && cdp > 1.e-6)
      chiral_check[j] =0;
    else{
      node0_printf("eigVec0[%d] is a null vector!\n",j);
      exit(1);
    }
  }
    /* the  mode  propagator matrix */
  /* I am stupid--how to do this in a 2-d array?? */
  if(Nvecs_hov != 0){
    ca= (Real *)malloc(num_masses*Nvecs_hov*sizeof(Real));
    cb= (Real *)malloc(num_masses*Nvecs_hov*sizeof(Real));
  }

  /* initialize the coefficients of the propagator matrix for modes */

  for(k=0;k<num_masses;k++)for(ivec=0;ivec<Nvecs_hov;ivec++){
    icount=Nvecs_hov*k + ivec;

    if(chiral_check[ivec]==0){
      mu=mass[k]/(2.0*R0);
      eps= sqrt(eigVal[ivec])/(2.0*R0);
      denom= (mu*mu+eps*eps*(1.0-mu*mu))*2.0*R0;
      ca[icount]= mu*(1.0-eps*eps)/denom;
      cb[icount]= eps*sqrt(1.0-eps*eps)/denom;
    }
    else{
      ca[icount]= 1.0/mass[k];
      cb[icount]= 0.0;
    }
    node0_printf("mass %e mode %d %d %e %e\n",mass[k],ivec,
                 chiral_check[ivec],ca[icount],cb[icount]);
  }
#endif


    /* open the prop files */

    for(k=0;k<num_masses;k++){
      fp_in_w[k]  = r_open_wprop(startflag_w[k], startfile_w[k]);
      fp_out_w[k] = w_open_wprop(saveflag_w[k],  savefile_w[k], wqs.type);
#ifdef H0INV
      h0_out_w[k] = w_open_wprop(saveflag_w3[k],  savefile_w3[k], wqs.type);
#endif
    }

  for(ksource = 0; ksource < wqs.nsource; ksource++){
    spin = convert_ksource_to_spin(ksource);
    color = convert_ksource_to_color(ksource);

//                /* Loop over source spins */
//    for(spin=0;spin<4;spin++){
//            /* Loop over source colors */
//    for(color=0;color<3;color++){

node0_printf("Propagator color %d spin %d\n",color,spin);
if(startflag_w[0] == FRESH){flag=0;}
else{
      /* check if there's a propagator already there--Do for all masses */
      flag=1;
      for(k=0;k<num_masses && flag==1 ;k++){
#ifdef IOTIME
      status = reload_wprop_sc_to_site( startflag_w[k], fp_in_w[k],
                                   &wqs, spin, color, F_OFFSET(psi),1);
#else
      status = reload_wprop_sc_to_site( startflag_w[k], fp_in_w[k],
                               &wqs, spin, color, F_OFFSET(psi),0);
#endif
      if(status != 0){
	node0_printf("congrad_outer_p: computing prop\n");
	/*
	reload_wprop_sc_to_site( FRESH, fp_in_w[k],
                               &wqs, spin, color, F_OFFSET(psi),0);
			       */
	flag = 0;
      }
      else{ /* status = 1--put the propagator in the new output file
so all the elements are in one place. This will fail if 
the propagator generation did not write the same number of elements
for each mass value propagator */
#ifdef IOTIME
                    save_wprop_sc_from_site( saveflag_w[k],fp_out_w[k],
                                    &wqs, spin,color,F_OFFSET(psi),1);
#else
                    save_wprop_sc_from_site( saveflag_w[k],fp_out_w[k],
                                    &wqs, spin,color,F_OFFSET(psi),0);
#endif
      }
      } /* k loop */
} /*startflag_w != FRESH */

      if(flag==0){  /* proceed to inversion */
      if(spin<2){my_chirality=1;chb=0;che=2;chbo=2;cheo=4;}
      else {my_chirality= -1;chb=2,che=4;chbo=0;cheo=2;}
      chirality_flag=my_chirality;

      /* Make source */

             /* Complete the source structure */

      /* NEEDS FIXING!! */
//            wqs.color = color;
//            wqs.spin = spin;

            /* For wilson_info */
            wqstmp = wqs;
	    //	    status = w_source_site(src,&wqs);
	    status = wv_source_site(src,&wqs);

	    /* check original source size... */
	    source_norm=0.0;
	    FORALLSITES(i,s){
	      source_norm += (double)magsq_wvec(((wilson_vector *)F_PT(s,src))  );
	    }
	    g_doublesum( &source_norm );

  if(this_node==0){
    printf("Original: source_norm = %e\n",source_norm);
    fflush(stdout);
  } 



	  FORALLSITES(i,s) copy_wvec((wilson_vector *)F_PT(s,src),&(s->chi0));
#ifdef EIGO
      /* project out the eigenvectors from the source */
node0_printf("removing %d modes from source\n",Nvecs_hov);
	  for(j=0;j<Nvecs_hov;j++){
	    cd=cmplx(0.0,0.0);
            FORALLSITES(i,s){
	      /* wproj will hold the chiral projections--
	       recall we have ``packed'' two chiralities into eigVec */
	      clear_wvec(&wproj);
	      for(ii=chb;ii<che;ii++)for(jj=0;jj<3;jj++){
		wproj.d[ii].c[jj]=eigVec[j][i].d[ii].c[jj];
	      }
	      ctmp =  wvec_dot( &(wproj),(wilson_vector *)F_PT(s,src));
	      CSUM(cd,ctmp);
	    }
	    g_complexsum(&cd);
	    cproj[j]=cd;
node0_printf("projector %d %e %e\n",j,cproj[j].real,cproj[j].imag);

	    CMULREAL(cd,-1.0,cd);

	    FORALLSITES(i,s){
	      clear_wvec(&wproj);
	      for(ii=chb;ii<che;ii++)for(jj=0;jj<3;jj++){
		wproj.d[ii].c[jj]=eigVec[j][i].d[ii].c[jj];
	      }
	      c_scalar_mult_add_wvec(&(s->chi0), &(wproj),
                             &cd, &(s->chi0) );
	    }
	  }
#endif


	  psim = (wilson_vector **)malloc(num_masses*sizeof(wilson_vector*));
	  for(i=0;i<num_masses;i++)
	    psim[i]=
	      (wilson_vector*)malloc(sites_on_node*sizeof(wilson_vector));

    FORALLSITES(i,s) {
                for(j=0;j<num_masses;j++){
                  clear_wvec(&(psim[j][i]));
                }}

	  node0_printf("inverting source color %d spin %d chirality %d\n",
color,spin,my_chirality);

	  size_r=1.e8;
	  avm_iters= -1;
	  if(MaxCG != 0)
	  avm_iters = congrad_multi_o( F_OFFSET(chi0),psim,
                                      MaxCG,RsdCG,&size_r,0); 

	  node0_printf("cg iters %d residue %e\n",avm_iters,size_r);
	  avs_iters += avm_iters;

          /* save the inverse of H(0)^2+s^2  for restarting */
#ifdef H0INV
	  for(j=0;j<num_masses;j++){
	    FORALLSITES(i,s) {
	      copy_wvec(&(psim[j][i]),&(s->psi));
/* test 
if(s->y == 0 && s->z == 0 && s->t == 0) dump_wvec(&(s->psi));
*/
	    }
#ifdef IOTIME
	    save_wprop_sc_from_site( saveflag_w3[j],h0_out_w[j],
			     &wqs, spin,color,F_OFFSET(psi),1);
#else
	    save_wprop_sc_from_site( saveflag_w3[j],h0_out_w[j],
			     &wqs, spin,color,F_OFFSET(psi),0);
#endif
	  }
#endif



#ifdef DEBUG
	    /*
if(this_node==0)printf("outer cg iters %d\n",avm_iters);
if(this_node==0)printf("outer cg res %e\n",size_r);
*/
#endif



    /* psim[k]= inverses of (H(0)^2+S) from non-eigenmode part.
For this part ONLY, compute the numerator, to get D^{-1} */
	  for(k=0;k<num_masses;k++){
	    one_minus_m=1.0- 0.5*mass[k]/R0;
	    FORALLSITES(i,s){
	      mult_by_gamma(&(psim[k][i]),&(s->psi),GAMMAFIVE);
	    }
	    hoverlap(F_OFFSET(psi),F_OFFSET(pm0));

	    FORALLSITES(i,s){
	      scalar_mult_wvec(&(s->pm0),one_minus_m,&(s->pm0));
	      scalar_mult_add_wvec(&(s->pm0),&(psim[k][i]),
				   mass[k],&(s->psi));
              scalar_mult_wvec(&(s->psi),coeff0[k],&(s->psi));
	    }

	    /* and the subtracted, rescaled propagator for true psibar-psi--
recall, we are only inverting in the subspace of missing modes! */
	    one_minus_m= 1.0/one_minus_m;
	    FORALLSITES(i,s){
	      /*
	      scalar_mult_add_wvec(&(s->psi),(wilson_vector *)F_PT(s,src),
				   r02inv,&(s->psi));
				   */
	      scalar_mult_add_wvec(&(s->psi),&(s->chi0),
				   r02inv,&(s->psi));
	      scalar_mult_wvec(&(s->psi),one_minus_m,&(s->psi));

	    }


#ifdef EIGO
      /* add the ''projector term'' back into dest...
	 Recall that eigValH0 holds the (packed) eigenvalues of H(0)^2*/

        /* form the wilson_matrix (propagator) from the modes */


/* cproj[ivec] is the overlap of the source with the ivec-th mode,
and we just have to remember what our chirality is... */
      for(ivec=0;ivec<Nvecs_hov;ivec++){
        icount = Nvecs_hov*k + ivec;

	    FORALLSITES(i,s){
	      clear_wvec(&wproj);

	      for(ii=chb;ii<che;ii++)for(jj=0;jj<3;jj++){
		wproj.d[ii].c[jj]=eigVec[ivec][i].d[ii].c[jj];
		CMULREAL(wproj.d[ii].c[jj],ca[icount],wproj.d[ii].c[jj]);
	      }

	      for(ii=chbo;ii<cheo;ii++)for(jj=0;jj<3;jj++){
		wproj.d[ii].c[jj]=eigVec[ivec][i].d[ii].c[jj];
		CMULREAL(wproj.d[ii].c[jj],cb[icount],wproj.d[ii].c[jj]);
		if(chbo==0) 
		  CMULREAL(wproj.d[ii].c[jj], -1.0,wproj.d[ii].c[jj]);
	      }
	      c_scalar_mult_add_wvec(&(s->psi), &(wproj),
                             &cproj[ivec], &(s->psi) );
	    }
 
      } /* modes */
#endif



	      /** DEBUGGING
	    FORALLSITES(i,s){
             printf("%d %d %d %d\n",s->x,s->y,s->z,s->t);
                dump_wvec(&(s->psi));
	    }
	    **/		


#ifdef IOTIME
                    save_wprop_sc_from_site( saveflag_w[k],fp_out_w[k],
                                    &wqs, spin,color,F_OFFSET(psi),1);
#else
                    save_wprop_sc_from_site( saveflag_w[k],fp_out_w[k],
                                    &wqs, spin,color,F_OFFSET(psi),0);
#endif

	    /*print out */

	  } /* k loop--num_masses */



	  for(i=0;i<num_masses;i++) free(psim[i]);
	  free(psim) ;
      }/* flag for inversion */
            }/* ksource (spins, colors)  */

	    /* close the prop files */
for(k=0;k<num_masses;k++){
  r_close_wprop(startflag_w[k],fp_in_w[k]);
  w_close_wprop(saveflag_w[k],fp_out_w[k]);
#ifdef H0INV
      w_close_prop(saveflag_w3[k],h0_out_w[k]);
#endif
}
#ifdef EIGO
 if(Nvecs_hov != 0){  
   free(cproj);
   free(chiral_check);
   free(ca);
   free(cb);
 }
#endif



return(avs_iters);

}



void setup_multi()
{
int i;


shift0=(Real *)malloc(num_masses*sizeof(Real));
coeff0=(Real *)malloc(num_masses*sizeof(Real));

for(i=0;i<num_masses;i++){
  coeff0[i]=1.0-0.25*mass[i]*mass[i]/R0/R0;
  shift0[i]=mass[i]*mass[i]/coeff0[i];
  coeff0[i]= 1.0/coeff0[i];


node0_printf("%d mass %e shift0 %e coeff0 %e\n",i,mass[i],
shift0[i],coeff0[i]);
}

}
