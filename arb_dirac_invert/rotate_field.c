/******* rotate_field.c -   FP fermions ****/
/* MIMD version 6 */

/*
dest= rotate*src


*/



#include "arb_dirac_inv_includes.h"
#include <string.h>

/* FP operators extrapolation */



static Real r_fit_params[9][3]={
{1.047189e-01, 3.287866e-01, -6.392501e-02},
{-3.407017e+00, -6.776410e-01, -4.341875e-02},
{-4.230266e+00, -9.707742e-01, -7.609335e-02},
{-4.949531e+00, -1.086064e+00, -1.146840e-01},
{-5.612031e+00, -1.143715e+00, -1.509096e-01},
{-2.474636e+00, -9.150334e-01, -1.917872e-02},
{-4.238351e+00, -9.947357e-01, -8.460870e-02},
{-5.506068e+00, -9.809474e-01, -1.387683e-01},
{-6.621816e+00, -8.942916e-01, -1.969786e-01}
};








void build_rot_params(Real mass_0, Real *lam, Real *rh)
{
int i,j,k;
Real ans;
Real p[3];
void fpoly(Real x,Real *p, int j);


/* lam and rh terms */
	fpoly(mass_0,p,3);

	ans=0.0;
	for(k=0;k<3;k++) ans += p[k]*r_fit_params[0][k];
	lam[0]= ans;


	for(j=1;j<5;j++) {

	ans=0.0;
	for(k=0;k<3;k++) ans += p[k]*r_fit_params[j][k];
	lam[j]= exp(ans);

	ans=0.0;
	for(k=0;k<3;k++) ans += p[k]*r_fit_params[j+4][k];
	rh[j]= exp(ans);


	}


}




void rotate_field(field_offset src,field_offset dest, Real mass_0)
{
	register int i;
	register site *s;
	msg_tag *tag[2];
	wilson_vector wvec1,wvec2;
	int ipath,n[4],k,mu,isign;
	Real lam[5],rh[5];
	void build_rot_params(Real mass_0, Real *lam, Real *rh);


	isign=1;



	build_rot_params(mass_0,lam,rh);



	/* begin with the origin */
	FORALLSITES(i,s) {
/*
                dump_wvec((wilson_vector *)F_PT(s,src));
*/

		scalar_mult_wvec((wilson_vector *)F_PT(s,src), lam[0],
		    (wilson_vector *)F_PT(s,dest));
	}




	/* now loop over the offsets and collect them */
	for(ipath=0;ipath<off_max;ipath++){

		k=label[ipath];
/*  printf("ipath=%d k=%d\n",ipath,k); */

		if(rh[k] != 0.0)if(lam[k]!= 0.0){
/*
printf("ipath=%d k=%d  %e   %e \n",ipath,k,rh[k],lam[k]); 
*/

			for(mu=0;mu<4;mu++) n[mu] = -offset[ipath][mu];


        /* Take  src displaced in up direction, gather it to "our site" */
 /*printf("ahead %d %d %d %d\n",offset[ipath][0], offset[ipath][1],
offset[ipath][2],offset[ipath][3]);  */

			tag[0] = start_general_gather_site( src,
			    sizeof(wilson_vector), offset[ipath],
				 EVENANDODD, gen_pt[0] );


        /* Take  src displaced in up direction, gathered,
                multiply it by link matrix,  and add to dest */
			wait_general_gather(tag[0]);
/* printf("mult 1\n"); */

			FORALLSITES(i,s) {
			mult_mat_wilson_vec( &(s->blocked_link[ipath]), 
				(wilson_vector * )(gen_pt[0][i]), &wvec1 ); 
			scalar_mult_add_wvec((wilson_vector *)F_PT(s,dest),
			    &wvec1,lam[k], 
			    (wilson_vector *)F_PT(s,dest));

		/* multiply by the four Dirac matrices ahead*/
			for(mu=0;mu<4;mu++)if(n[mu]!=0){
			mult_by_gamma(&wvec1,&wvec2,mu);

			scalar_mult_add_wvec((wilson_vector *)F_PT(s,dest),
			    &wvec2, -isign*n[mu]*rh[k], 
			    (wilson_vector *)F_PT(s,dest));
			} /* dirac */

				}

		cleanup_general_gather(tag[0]);

        /* Take  src displaced in down direction,
        multiply it by adjoint link matrix, gather it "up" */


			FORALLSITES(i,s){
			  mult_adj_mat_wilson_vec(&(s->blocked_link[ipath]),
			  (wilson_vector *)F_PT(s,src), &(s->htmp[1]));
			}

/* printf("behind %d %d %d %d\n",n[0],n[1],n[2],n[3]); */
			tag[1] = start_general_gather_site( F_OFFSET(htmp[1]),
			    sizeof(wilson_vector), n, EVENANDODD, gen_pt[1] );
			wait_general_gather(tag[1]);
/* printf("mult 2\n"); */

	/* multiply by the scalar term behind*/

			FORALLSITES(i,s) {
			scalar_mult_add_wvec((wilson_vector *)F_PT(s,dest),
			    (wilson_vector *)gen_pt[1][i],lam[k], 
			    (wilson_vector *)F_PT(s,dest));
				}

		/* multiply by the four Dirac matrices behind*/
			for(mu=0;mu<4;mu++)if(n[mu]!=0){
			  FORALLSITES(i,s) {
			mult_by_gamma((wilson_vector *)gen_pt[1][i],&wvec1,mu);

			scalar_mult_add_wvec((wilson_vector *)F_PT(s,dest),
			    &wvec1,isign*n[mu]*rh[k], 
			    (wilson_vector *)F_PT(s,dest));
					}
			} /* dirac */


		cleanup_general_gather(tag[1]);
		} /* if (rh and lam) */


	} /* offsets */


}
/* rotate_field.c */


