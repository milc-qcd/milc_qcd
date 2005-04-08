/******** hax.c *************/
/* MIMD version 6 */
#include "arb_dirac_inv_includes.h"
void mult_by_gamma_l( spin_wilson_vector *src, spin_wilson_vector *dest, 
		     int dir);
void mult_by_gamma_r( spin_wilson_vector *src, spin_wilson_vector *dest, 
		     int dir);

/*  ``hyp'' axial current for hypercube */

void hax(field_offset src1,field_offset src2,Real prop[MAX_P][MAX_NT]) 
{
/* pion and rho with pointlike sink, p=0 */

register int i,x,y,z,t;
register site *s; 

int my_t,k;
int cf, sf, ci, si;
complex g1,g2;
int j;
Real co;

wilson_vector wvec1,wvec2,wvec3,wvec4;

spin_wilson_vector localmat;       /* temporary storage for quark */

spin_wilson_vector antiquark;      /* temporary storage for antiquark */
spin_wilson_vector quark;      /* temporary storage for quark */

/* CVC stuff */
msg_tag *tag[2];
int ipath,n[4],mu,nu;
/* data for paths: [k][npath[k]][along each path] */
int dir[5][24][4],sign[5][24][4],length[5][24],npath[5];
Real cc[5][24];

int ioffset, iflag,l1,l2,keep,ahead;
int mid_offset[4],end_offset[4],n_offset[4],dd;
Real pathsign;

int ix,iy,iz,my_x,my_y,my_z;

/* for nonzero momentum */
Real cx,cy,cz,cxy,cyz,cxz,c111;




nu=TUP;{




	FORALLSITES(i,s){
		s->cvc[nu] = cmplx(0.0,0.0);
		}


        /* now loop over the offsets and collect them */
        for(ioffset=0;ioffset<off_max;ioffset++){
                k=label[ioffset];

/* this generates the paths in the action */
        path_stuff(iflag,ioffset,k,npath,dir,sign,length,cc);


/* now loop over every path */
	for(ipath=0;ipath<npath[k];ipath++){
	keep=0;

/* zero mid and set end to the end of the path */
          for(l2=0;l2<4;l2++) mid_offset[l2]=0;
          for(l2=0;l2<4;l2++) end_offset[l2]=offset[ioffset][l2];

/* examine the path for the index nu */
          for(l1=0;l1<length[k][ipath];l1++){
		if(dir[k][ipath][l1]==nu){
		keep=1;
		pathsign=cc[k][ipath];
		  if(sign[k][ipath][l1]<0)
			{
			pathsign= -cc[k][ipath];
			dd=dir[k][ipath][l1];mid_offset[dd]-=1;
			}
		break;
		}
		else{dd=dir[k][ipath][l1];mid_offset[dd]+=sign[k][ipath][l1];}
        }


	if(keep){

/*  multiply links along path (uses tempmat2, returns tempmat1) */

                  path(dir[k][ipath],sign[k][ipath],length[k][ipath]);
/* path.c puts the end of the path at location ``zero.''
We must gather so that the beginning of the path (psibar end) is at zero */
                tag[0]=start_general_gather_site(F_OFFSET(tempmat1),
		sizeof(su3_matrix),
                offset[ioffset], EVENANDODD, gen_pt[0] );

                wait_general_gather(tag[0]);

                FORALLSITES(i,s){
                        su3mat_copy((su3_matrix *)(gen_pt[0][i]),
                                &(s->tempmat2));
		}
                FORALLSITES(i,s){
                        su3mat_copy(&(s->tempmat2),
                                &(s->tempmat1));
		}
                cleanup_general_gather(tag[0]);

/* we have figured out a contribution to the current. The antiquark propagator
ends on ``our site'' while the quark is at site end_offset, and the current
is at site mid_offset. We have to run through this loop twice, for the
``ahead'' and ``behind'' offsets of the action */

		for(ahead=0;ahead<2;ahead++){

                if(rho[k]!= 0.0){


/* for the behind case, gather the path to our site, too
(ordinarily, we would multiply and then gather, but too much storage
needed...) */
	     if(ahead==1){
		pathsign *= -1.0;
		for(l2=0;l2<4;l2++){mid_offset[l2]*= -1;end_offset[l2]*= -1;}

		tag[0]=start_general_gather_site(F_OFFSET(tempmat1),
		sizeof(su3_matrix),
                end_offset, EVENANDODD, gen_pt[0] );

                wait_general_gather(tag[0]);

                FORALLSITES(i,s){
                        su3mat_copy((su3_matrix *)(gen_pt[0][i]),
                                &(s->tempmat2));
		}
                cleanup_general_gather(tag[0]);
                 FORALLSITES(i,s){
                        su3_adjoint(&(s->tempmat2),
                                &(s->tempmat1));
		}

	     }


/* we gather the quark sitting on site end_offset to our site */

                  tag[0] = start_general_gather_site( src2,
                  sizeof(spin_wilson_vector), end_offset,
                                 EVENANDODD, gen_pt[0] );

                        wait_general_gather(tag[0]);


				

/* next we multiply up all the traces for this particular combination 
of offsets (yes, it is expensive) */

    FORALLSITES(i,s) {
	s->cvct[nu]= cmplx(0.0,0.0);

        /*first, assemble the antiquark */

        /*  antiquark = c.c. of quark propagator */
        for(si=0;si<4;si++)
        for(sf=0;sf<4;sf++)
        for(cf=0;cf<3;cf++){
            CONJG(((spin_wilson_vector *)F_PT(s,src1))->d[si].d[sf].c[cf],
                antiquark.d[si].d[sf].c[cf]);
         }


        /* left multiply antiquark by source gamma matrices,
           beginning with gamma_5 for quark -> antiquark */
        mult_by_gamma_l( &(antiquark), &(localmat), GAMMAFIVE);


        /* right dirac multiplication by gamma-5 (finishing up antiquark) */
        mult_by_gamma_r( &(localmat),&(antiquark), GAMMAFIVE);


        /* dirac multiplication by the source gamma matrices (on left) */
        mult_by_gamma_l( &(antiquark), &(localmat), GAMMAFIVE);
        /* dirac multiplication by the sink gamma-5 matrix (on right) */
        mult_by_gamma_r(&(localmat), &(antiquark),  GAMMAFIVE);


/* now we can  multiply the quark by the link and by the
dirac matrices for the current operator */

        quark= *(spin_wilson_vector *)(gen_pt[0][i]);


        for(si=0;si<4;si++){

            mult_mat_wilson_vec(&(s->tempmat1),&(quark.d[si]), 
			&wvec1);


		clear_wvec(&wvec4);

		for(mu=0;mu<4;mu++) n[mu] = -end_offset[mu];

                /* multiply by the four Dirac matrices*/
                        for(mu=0;mu<4;mu++)if(n[mu]!=0){
                        mult_by_gamma(&wvec1,&wvec2,mu);
                        scalar_mult_add_wvec(&wvec4,
                            &wvec2, -pathsign*n[mu]*rho[k], 
                            &wvec4);
                        } /* dirac */



	/* trace over propagators to form the (offset) cvc */
	    for(sf=0;sf<4;sf++)
	      for(cf=0;cf<3;cf++) {
		  g1 = wvec4.d[sf].c[cf];
                  g2 = antiquark.d[si].d[sf].c[cf];
		  s->cvct[nu].real += (g1.real*g2.real - g1.imag*g2.imag);
		  s->cvct[nu].imag += (g1.real*g2.imag + g1.real*g2.imag);

			}
	} /* si  */
    } /* sites */

                cleanup_general_gather(tag[0]);

/* we have computed the current from one contribution, but it is on site
mid_offset. We must gather up from below and add to current*/

		for(l2=0;l2<4;l2++){n_offset[l2]= -mid_offset[l2];}

                tag[0]=start_general_gather_site(F_OFFSET(cvct[nu]),
		sizeof(complex),
                n_offset, EVENANDODD, gen_pt[0] );

                wait_general_gather(tag[0]);
                FORALLSITES(i,s){CSUM(s->cvc[nu], *(complex *)(gen_pt[0][i]));}
                cleanup_general_gather(tag[0]);



		} /* rho[k] and lambda[k] */

	} /* ahead and behind */

	} /* keep */

} /* ipath */

} /* ioffset */


/* print out time averaged observable */

    FORALLSITES(i,s) {
        my_t = s->t;


                for(j=0;j<3;j++){
                cz=cos(2.0*PI/(Real)nz*(Real)(s->z)*(Real)j);
                cx=cos(2.0*PI/(Real)nx*(Real)(s->x)*(Real)j);
                cy=cos(2.0*PI/(Real)ny*(Real)(s->y)*(Real)j);


                prop[j][my_t] += (s->cvc[nu].real)*
                                (cx+cy+cz)/3.0;
                }
                cxy=cos(2.0*PI/(Real)nz*(Real)(s->x +s->y));
                cxz=cos(2.0*PI/(Real)nz*(Real)(s->x +s->z));
                cyz=cos(2.0*PI/(Real)nz*(Real)(s->y +s->z));
                c111=cos(2.0*PI/(Real)nz*(Real)(s->x +s->y + s->z));
                prop[3][my_t] += (s->cvc[nu].real)*(cxy+cyz+cxz)/3.0;
                prop[4][my_t] += (s->cvc[nu].real)*c111;
/* to print out the correlator 
if(my_t== nt/4){
            my_x = ((s->x)+nx) % nx;
            ix = (my_x < (nx-my_x)) ?  my_x :  (my_x-nx);
            my_y = ((s->y)+ny) % ny;
            iy = (my_y < (ny-my_y)) ?  my_y :  (my_y-ny);
            my_z = ((s->z)+nz) % nz;
            iz = (my_z < (nz-my_z)) ?  my_z :  (my_z-nz);
printf("VECTOR %d %d %d %e\n",ix,iy,iz,s->cvc[nu].real);
}
*/
	}


} /* nu */

/* now let's check current conservation: j_mu(n)-j_mu(n-mu) 



    FORALLSITES(i,s) {
	s->cvct[0]=0.0;
	}

	for(nu=0;nu<4;nu++){
        tag[0]=start_gather_site( F_OFFSET(cvc[nu]), sizeof(Real),
            OPP_DIR(nu), EVENANDODD, gen_pt[0] );
                wait_gather(tag[0]);

	FORALLSITES(i,s){
		s->cvct[0] += (s->cvc[nu]- *(Real *)(gen_pt[0][i]));}
	} 

    FORALLSITES(i,s) {
	printf("%d %d %d %d %e %e %e %e : %e\n",
 s->x, s->y, s->z,s->t,s->cvc[0],s->cvc[1],s->cvc[2],s->cvc[3],s->cvct[0]);	}
 cleanup_gather(tag[0]);
*/


} /* spectrum */

