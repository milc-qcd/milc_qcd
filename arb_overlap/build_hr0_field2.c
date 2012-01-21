/************* build_hr0_field2.c ****************/
/* MIMD version 7 */

#include "arb_ov_includes.h"

void build_hr0( int flag, int precflag )/*  0, fresh start, 1--restart with stored vectors,*/
{

    register int i,j;
    register site *s;
    int total_R_iters ;
    Real re,im;
    complex cc;
    int nn0;
    Real eigenval_tol_save;

    /* code for good h(-r0) eigenmodes */
    int failed,ifailed,nextra;
    double *MyeigVal0;
    wilson_vector **MyeigVec0;
    wilson_vector *tmpvec;
    Real ftmp;

    Real invp,invp5;



    tmpvec=(wilson_vector*)malloc(sites_on_node*sizeof(wilson_vector));
    nn0=-ndelta0;

    if (precflag==HIGHP) 
    {
	eigenval_tol=eigenval_tol_high;
	error_decr=error_decr_high;
    } else {
	eigenval_tol=eigenval_tol_low;
	error_decr=error_decr_low;
    }


    node0_printf("\n");
    eigenval_tol_save=eigenval_tol;


    build_params(-R0);
    make_clov1();
    kind_of_h0=HZERO;
    if(Maxr0Iter != 0){

      /*
     FORALLSITES(i,s){
       clear_wvec(&(eigVec0[0][i]));
       if(s->x==0 && s->y==0 && s->z==0 && s->t==0) eigVec0[0][i].d[0].c[0].real=1.0;
     }
    delta0_field(eigVec0[0],eigVec0[1],PLUS);
     FORALLSITES(i,s){
       if(s->x ==0 && s ->y==0 && s->z==0){
         node0_printf("DT %d %d %d %d\n",s->x,s->y,s->z,s->t);
         dump_wvec(&(eigVec0[1][i]));
       }
     }
      */


	/* initialize the CG vectors */
	if(flag==0)
	{
	    if(this_node==0) printf("random initial vectors for h(m0)\n");
	    /* Initiallize all the eigenvectors to a random vector */
	    for(j=0;j<Nvecs_h0r0;j++){
		grsource_dirac(EVENANDODD);
		FORSOMEPARITY(i,s,EVENANDODD){
		    copy_wvec(&(s->g_rand),&(eigVec0[j][i]));
		}
		eigVal0[j]=1.0e+16;
	    }
	}
	else 
	{

	    /* eigenvectors and eigenvalues of h(-r0) already exist, but we need eigenvalues of h^2 */
	    for(j=0;j<Nvecs_h0r0;j++){
		eigVal0[j] = eigVal0[j]*eigVal0[j];
	    }
	}

	failed=0;
	nextra=0;
	MyeigVec0 = (wilson_vector**)malloc((Nvecs_h0r0+5)*sizeof(wilson_vector*));
	MyeigVal0= (double *)malloc((Nvecs_h0r0+5)*sizeof(double));
	for (i=0;i<Nvecs_h0r0;i++) {
	    MyeigVec0[i]=eigVec0[i];
	    MyeigVal0[i]= eigVal0[i];
	}

	for (i=Nvecs_h0r0;i<Nvecs_h0r0+5;i++)
	    MyeigVec0[i]=(wilson_vector*)malloc(sites_on_node*sizeof(wilson_vector));



	for(ifailed=0;ifailed< 5;ifailed++){

	    if(ifailed > 0){
		node0_printf("\nRETRY %d\n",ifailed);
		MyeigVec0[Nvecs_h0r0+ifailed-1]=
		    (wilson_vector*)malloc(sites_on_node*sizeof(wilson_vector));
		grsource_dirac(EVENANDODD);
		FORSOMEPARITY(i,s,EVENANDODD){
		    copy_wvec(&(s->g_rand),&(MyeigVec0[Nvecs_h0r0+ifailed-1][i]));
		}
		MyeigVal0[Nvecs_h0r0+ifailed-1]=1.0e+16;
		/*
		eigenval_tol/=10.;
		*/
		eigenval_tol /= 2.0;
		if (eigenval_tol <1.e-14)eigenval_tol=1.e-14;
	    }
	    node0_printf("calling hr0 with Maxr0Iter %d\n",Maxr0Iter);

	    total_R_iters=Kalkreuter(MyeigVec0, MyeigVal0, eigenval_tol, 
		    error_decr, Nvecs_h0r0+ifailed, Maxr0Iter, Restart, 
		    Kiters, EVENANDODD) ;
	    if(this_node==0)printf("total Rayleigh iters = %d\n",total_R_iters);
	    fflush(stdout);




	    for(j=0;j<Nvecs_h0r0+ifailed;j++){
		node0_printf("F3HSQ %d %e\n",j,MyeigVal0[j]);
	    }
	    fflush(stdout);

	    /* now we do tests on eigenvectors */
	    failed=0;
	    for(j=0;j<Nvecs_h0r0;j++){


		/* expectation value of action r=g5 dslash chi*/
		delta0_field(eigVec0[j],tmpvec,PLUS);
		FORALLSITES(i,s){
		    mult_by_gamma(&tmpvec[i],&tmpvec[i],GAMMAFIVE);
		}
		re=im=0.0;
		FORALLSITES(i,s){
		    cc = wvec_dot( &(eigVec0[j][i]), &(tmpvec[i]) );
		    re += cc.real ;
		    im += cc.imag ;
		}
		g_floatsum(&re);
		g_floatsum(&im);

		eigVal0[j]=re;


		FORALLSITES(i,s)
		    scalar_mult_add_wvec(&(tmpvec[i]),&eigVec0[j][i],
			    (Real)-eigVal0[j],&(tmpvec[i]));
		ftmp=vectornorm(tmpvec);



		if(this_node==0) printf("F3MEX %d %e %e %e ",
			j,(double)eigVal0[j], (double)(im),ftmp);


		if( fabs(ftmp) >= eigenvec_quality) {
		    failed++;
		    node0_printf("F");
		}
		node0_printf("\n");
	    }
	    if(failed==0) break;
	} /* failed */

	fflush(stdout);

	for (i=Nvecs_h0r0;i<Nvecs_h0r0+5;i++) free(MyeigVec0[i]);
	free(MyeigVec0);
	free(MyeigVal0);

	node0_printf("\n");
	node0_printf("START_EIGVALS of h(-r0): re im  IPR invp5\n");

	for(j=0;j<Nvecs_h0r0;j++){
	    /* expectation value of action r=g5 dslash chi*/
	    delta0_field(eigVec0[j],tmpvec,PLUS);
	    FORALLSITES(i,s){
		mult_by_gamma(&tmpvec[i],&tmpvec[i],GAMMAFIVE);
	    }
	    re=im=0.0;
	    FORALLSITES(i,s){
		cc = wvec_dot( &(eigVec0[j][i]), &(tmpvec[i]) );
		re += cc.real ;
		im += cc.imag ;
	    }
	    g_floatsum(&re);
	    g_floatsum(&im);

	    eigVal0[j]=re;

	    /* inverse participation for eigenmodes h(-r0) */
	    invp5=0.0;
	    FORALLSITES(i,s){
		mult_by_gamma(&(eigVec0[j][i]),&tmpvec[i],GAMMAFIVE);
		cc = wvec_dot( &(eigVec0[j][i]), &tmpvec[i]);
		invp5 += cc.real*cc.real + cc.imag*cc.imag;
	    }
	    g_floatsum(&invp5);
	    invp=0.0;
	    FORALLSITES(i,s){
		cc = wvec_dot( &(eigVec0[j][i]), &(eigVec0[j][i]) );
		invp += cc.real*cc.real + cc.imag*cc.imag;
	    }
	    g_floatsum(&invp);
	    invp *= (Real)volume;
	    invp5 *= (Real)volume;

	    if(this_node==0) printf("F3MES %d %e %e   %e %e\n",
		    j,(double)eigVal0[j],
		    (double)(im),invp,invp5);
	    fflush(stdout);
	    /*
	    if( fabs(eigVal0[j]*eigVal0[j]-MyeigVal0[j]) >= 1.e-5) failed++;
	    */
	}

    } /* Maxr0Iter != 0 */

    else{
	/* if we're not iterating, need to compute the eigenvalue of h(-r0) */
	node0_printf("\n");
	node0_printf("START_EIGVALS of h(-r0): re im  IPR invp5\n");

	for(j=0;j<Nvecs_h0r0;j++){
	    /* expectation value of action r=g5 dslash chi*/
	    delta0_field(eigVec0[j],tmpvec,PLUS);
	    FORALLSITES(i,s){
		mult_by_gamma(&tmpvec[i],&tmpvec[i],GAMMAFIVE);
	    }
	    re=im=0.0;
	    FORALLSITES(i,s){
		cc = wvec_dot( &(eigVec0[j][i]), &(tmpvec[i]) );
		re += cc.real ;
		im += cc.imag ;
	    }
	    g_floatsum(&re);
	    g_floatsum(&im);

	    eigVal0[j]=re;
	    FORALLSITES(i,s)scalar_mult_add_wvec(&(tmpvec[i]),&eigVec0[j][i],(Real)-eigVal0[j],&(tmpvec[i]));
	    ftmp=vectornorm(tmpvec);
	    node0_printf("QUALITY OF EIGENMODE %i : %e\n",j,ftmp);
	    /* inverse participation for eigenmodes h(-r0) */
	    invp5=0.0;
	    FORALLSITES(i,s){
		mult_by_gamma(&(eigVec0[j][i]),&tmpvec[i],GAMMAFIVE);
		cc = wvec_dot( &(eigVec0[j][i]), &tmpvec[i]);
		invp5 += cc.real*cc.real + cc.imag*cc.imag;
	    }
	    g_floatsum(&invp5);
	    invp=0.0;
	    FORALLSITES(i,s){
		cc = wvec_dot( &(eigVec0[j][i]), &(eigVec0[j][i]) );
		invp += cc.real*cc.real + cc.imag*cc.imag;
	    }
	    g_floatsum(&invp);
	    invp *= (Real)volume;
	    invp5 *= (Real)volume;

	    if(this_node==0) printf("F3MES %d %e %e   %e %e\n",
		    j,(double)eigVal0[j],
		    (double)(im),invp,invp5);


	    fflush(stdout);
	}
    } /* else */

    node0_printf("END_EIGVALS of h(-r0)\n");
    node0_printf("\n");
    node0_printf("mxv operations for eigenvecs %ld\n",ndelta0+nn0);
    /*
       eigenval_tol=eigenval_tol_save;
       resid_inner=resid_inner_save;

*/
    free(tmpvec);
    fflush(stdout);
    eigenval_tol=eigenval_tol_save;
}

