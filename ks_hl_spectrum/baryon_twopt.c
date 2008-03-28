/***************** baryon_twpt.c ********************************************/
/* MIMD version 7 */
/* Heechang Na 2006/7 */

/* Calculate Two point function of baryon with two staggered light quarks, and
   one wilson heavy quark. The 2-pt for the operator O5 is;
   
    (2)                   -i*p*x   
   C      (p,t) = SUM    e       (-4)     epsilon     epsilon         
    alpha,beta     x                             a b c       a' b' c'

                  a a'        b b'      c c'
               * G    (x,0)  G   (x,0) G    (x,0)
                  1 chi       2 chi      H  alpha,beta

  where a,b, and c are color indices, and mu, nu are spin indices.
  Since G_1, and G_2 are single component staggered propagators, 
  they don't have spin indices.
  G_H is wilson heavy quark propagator, which has spin and color index. 


  for the operator Omu is;


    (2)                   -i*p*x    x_mu
   C      (p,t) = SUM    e       (-4)     epsilon     epsilon         
     mu,nu         x                             a b c       a' b' c'

                  a a'        b b'      c c'
               * G    (x,0)  G   (x,0) G    (x,0)
                  1 chi       2 chi      H  mu,nu

  where mu, mu are the Lorentz indices for the spin directions, and when
  mu is not equal to nu, 2-pt is zero, because the cancellation of the copy
  indices. So, we calculate only for mu = nu. 


  for the operator Omu_HH which is for doubly heavy baryons;


    (2)                   -i*p*x    
   C      (p,t) = SUM    e         epsilon     epsilon         
     mu,nu         x                      a b c       a' b' c'

                               a a'                    b b'         
               *tr[gamma_nu C G   (x,0)  C  gamma_mu  G    (x,0) ]  
                               H                       H            

                           c c'
               * Omega(x) G   (x,0) 
                           chi

  where C is the charge conjugation operator, and Omega is the spin diagonalize
  matrix which is multiplied on the single component staggered quark, G_chi, 
  in order to convert to naive quark from staggered quark.  
  -Heechang 3/20/2008
*/                         


#include "ks_hl_spectrum_includes.h"


int spin_trace_two_clover_double(wilson_propagator *quark1,int mu, int nu, double_complex tr[3][3][3][3])
{

	wilson_propagator heavy_quark1, heavy_quark2;
	wilson_propagator tmp1, tmp2;
	double_complex temp, temp1;
    double_complex temp_a, temp_b;	
	int si,sf,ss,a,b,ap,bp,j;

	/* Find out : Trace = tr[ ( g_nu g_2 g_0 G^T ) ( g_0 g_2 g_mu G ) ] 
	   Note that the trace is only for spin. -Heechang*/

	heavy_quark1 = * quark1;
	heavy_quark2 = * quark1;

	/*Transpose of the first quark*/
	for(si=0;si<4;si++)for(sf=0;sf<4;sf++)for(a=0;a<3;a++)for(b=0;b<3;b++){
		tmp1.c[a].d[si].d[sf].c[b] = heavy_quark1.c[a].d[sf].d[si].c[b];
	}
	heavy_quark1 = tmp1;

	/*( g_nu g_2 g_0 G^T ) */
	for(a=0;a<3;a++){
		mult_swv_by_gamma_l( &(heavy_quark1.c[a]), &(tmp1.c[a]), TUP);
		heavy_quark1.c[a] = tmp1.c[a];
	} 
	for(a=0;a<3;a++){
		mult_swv_by_gamma_l( &(heavy_quark1.c[a]), &(tmp1.c[a]), YUP);
		heavy_quark1.c[a] = tmp1.c[a];
	}
	for(a=0;a<3;a++){
		mult_swv_by_gamma_l( &(heavy_quark1.c[a]), &(tmp1.c[a]), nu);
		heavy_quark1.c[a] = tmp1.c[a];
	}

	/*( g_0 g_2 g_mu G ) */
	for(a=0;a<3;a++){
		mult_swv_by_gamma_l( &(heavy_quark2.c[a]), &(tmp2.c[a]), mu);
		heavy_quark2.c[a] = tmp2.c[a];
	}
	for(a=0;a<3;a++){
		mult_swv_by_gamma_l( &(heavy_quark2.c[a]), &(tmp2.c[a]), YUP);
		heavy_quark2.c[a] = tmp2.c[a];
	}
	for(a=0;a<3;a++){
		mult_swv_by_gamma_l( &(heavy_quark2.c[a]), &(tmp2.c[a]), TUP);
		heavy_quark2.c[a] = tmp2.c[a];
	}

	/*Do product and trace */

	/* Simple version to product*/
	/*
	for(a=0;a<3;a++)for(b=0;b<3;b++)for(ap=0;ap<3;ap++)for(bp=0;bp<3;bp++){
		tr[a][ap][b][bp].real = 0; tr[a][ap][b][bp].imag = 0;

		for(si=0;si<4;si++){
			temp1.real = 0; temp1.imag = 0;

			for(ss=0;ss<4;ss++){
				temp_a.real = (double) heavy_quark1.c[a].d[si].d[ss].c[ap].real;
				temp_a.imag = (double) heavy_quark1.c[a].d[si].d[ss].c[ap].imag;
				temp_b.real = (double)  heavy_quark2.c[b].d[ss].d[si].c[bp].real;
				temp_b.imag = (double)  heavy_quark2.c[b].d[ss].d[si].c[bp].imag;
				CMUL(temp_a,temp_b,temp);
				CSUM(temp1 , temp);
			}

			CSUM(tr[a][ap][b][bp], temp1);
		}
	}           */

	/*However, we know which color combinations will be non-zero, so
	  it's fast version */
	for(a=0;a<3;a++){ ap = a;
		for(j=-1;j<=1;j=j+2){
			b = bp = (a+j+3)%3;

			tr[a][ap][b][bp].real = 0; tr[a][ap][b][bp].imag = 0;

			for(si=0;si<4;si++){
				temp1.real = 0; temp1.imag = 0;

				for(ss=0;ss<4;ss++){
					temp_a.real = (double) heavy_quark1.c[a].d[si].d[ss].c[ap].real;
					temp_a.imag = (double) heavy_quark1.c[a].d[si].d[ss].c[ap].imag;
					temp_b.real = (double)  heavy_quark2.c[b].d[ss].d[si].c[bp].real;
					temp_b.imag = (double)  heavy_quark2.c[b].d[ss].d[si].c[bp].imag;
					CMUL(temp_a,temp_b,temp);
					CSUM(temp1 , temp);
				}

				CSUM(tr[a][ap][b][bp], temp1);
			}
		}
	}

	for(a=0;a<3;a++){ ap = a;
		for(j=-1;j<=1;j=j+2){
			b = (a+j+3)%3;
			bp = (a-j+3)%3;
			tr[a][ap][b][bp].real = 0; tr[a][ap][b][bp].imag = 0;

			for(si=0;si<4;si++){
				temp1.real = 0; temp1.imag = 0;

				for(ss=0;ss<4;ss++){
					temp_a.real = (double) heavy_quark1.c[a].d[si].d[ss].c[ap].real;
					temp_a.imag = (double) heavy_quark1.c[a].d[si].d[ss].c[ap].imag;
					temp_b.real = (double)  heavy_quark2.c[b].d[ss].d[si].c[bp].real;
					temp_b.imag = (double)  heavy_quark2.c[b].d[ss].d[si].c[bp].imag;
					CMUL(temp_a,temp_b,temp);
					CSUM(temp1 , temp);
				}

				CSUM(tr[a][ap][b][bp], temp1);
			}
		}
	}

	for(a=0;a<3;a++){ bp = a;
		for(j=-1;j<=1;j=j+2){
			b = ap = (a+j+3)%3;

			tr[a][ap][b][bp].real = 0; tr[a][ap][b][bp].imag = 0;

			for(si=0;si<4;si++){
				temp1.real = 0; temp1.imag = 0;

				for(ss=0;ss<4;ss++){
					temp_a.real = (double) heavy_quark1.c[a].d[si].d[ss].c[ap].real;
					temp_a.imag = (double) heavy_quark1.c[a].d[si].d[ss].c[ap].imag;
					temp_b.real = (double)  heavy_quark2.c[b].d[ss].d[si].c[bp].real;
					temp_b.imag = (double)  heavy_quark2.c[b].d[ss].d[si].c[bp].imag;
					CMUL(temp_a,temp_b,temp);
					CSUM(temp1 , temp);
				}

				CSUM(tr[a][ap][b][bp], temp1);
			}
		}
	}


	for(a=0;a<3;a++){ bp = a;
		for(j=-1;j<=1;j=j+2){
			b  = (a+j+3)%3;
			ap = (a-j+3)%3;
			tr[a][ap][b][bp].real = 0; tr[a][ap][b][bp].imag = 0;

			for(si=0;si<4;si++){
				temp1.real = 0; temp1.imag = 0;

				for(ss=0;ss<4;ss++){
					temp_a.real = (double) heavy_quark1.c[a].d[si].d[ss].c[ap].real;
					temp_a.imag = (double) heavy_quark1.c[a].d[si].d[ss].c[ap].imag;
					temp_b.real = (double)  heavy_quark2.c[b].d[ss].d[si].c[bp].real;
					temp_b.imag = (double)  heavy_quark2.c[b].d[ss].d[si].c[bp].imag;
					CMUL(temp_a,temp_b,temp);
					CSUM(temp1 , temp);
				}

				CSUM(tr[a][ap][b][bp], temp1);
			}
		}
	}

	for(a=0;a<3;a++){ 
		for(j=-1;j<=1;j=j+2){
			b = bp = (a+j+3)%3;
			ap = (a-j+3)%3;
			tr[a][ap][b][bp].real = 0; tr[a][ap][b][bp].imag = 0;

			for(si=0;si<4;si++){
				temp1.real = 0; temp1.imag = 0;

				for(ss=0;ss<4;ss++){
					temp_a.real = (double) heavy_quark1.c[a].d[si].d[ss].c[ap].real;
					temp_a.imag = (double) heavy_quark1.c[a].d[si].d[ss].c[ap].imag;
					temp_b.real = (double)  heavy_quark2.c[b].d[ss].d[si].c[bp].real;
					temp_b.imag = (double)  heavy_quark2.c[b].d[ss].d[si].c[bp].imag;
					CMUL(temp_a,temp_b,temp);
					CSUM(temp1 , temp);
				}

				CSUM(tr[a][ap][b][bp], temp1);
			}
		}
	}

	for(a=0;a<3;a++){ 
		for(j=-1;j<=1;j=j+2){
			b = ap = (a+j+3)%3;
			bp = (a-j+3)%3;
			tr[a][ap][b][bp].real = 0; tr[a][ap][b][bp].imag = 0;

			for(si=0;si<4;si++){
				temp1.real = 0; temp1.imag = 0;

				for(ss=0;ss<4;ss++){
					temp_a.real = (double) heavy_quark1.c[a].d[si].d[ss].c[ap].real;
					temp_a.imag = (double) heavy_quark1.c[a].d[si].d[ss].c[ap].imag;
					temp_b.real = (double)  heavy_quark2.c[b].d[ss].d[si].c[bp].real;
					temp_b.imag = (double)  heavy_quark2.c[b].d[ss].d[si].c[bp].imag;
					CMUL(temp_a,temp_b,temp);
					CSUM(temp1 , temp);
				}

				CSUM(tr[a][ap][b][bp], temp1);
			}
		}
	}


	return 0;
}



/*for Omu operator in doubly heavy baryons*/
int ks_baryon_2point_Omu_HH(field_offset ks_1, field_offset heavy_quark, double_complex *propagator[4][4][3][3])
{
	int i, j, t, my_x, my_y, my_z,  mu, nu;
	site *s;
	double_complex epx, epx1;
	double pi, mom[3];
	int eps;
	int a, b, c, ap, bp, cp; /*color indices. bp -> b prime */
	int si,sf; /*dirac indices for heavy quark propagator */

	su3_matrix *ks_1st_quark;
	wilson_propagator *quark1;
	wilson_propagator light, tmp;
	double_complex tr[3][3][3][3]; /* The trace has four color indices */
	double_complex prop_matrix, prop_tot[4][4], prop_all;
	complex zero;

		int p[3] ={0,0,0};

	pi = 4.0 * atan( 1.);
	mom[0] = -2.*pi/(float)nx;
	mom[1] = -2.*pi/(float)ny;
	mom[2] = -2.*pi/(float)nz;

	zero.real = 0.0; zero.imag = 0.0;

	FORALLSITES(i,s){
		t = s->t;
		my_x = s->x;
		my_y = s->y;
		my_z = s->z;

		/* asign each field */
		ks_1st_quark = (su3_matrix *)F_PT(s, ks_1);
		quark1 = (wilson_propagator *)F_PT(s, heavy_quark);


		/* We need light quark with spin indices, i.e. multiply by Omega matrix */
		/* First, just let me construct a wilson_propagator with unit matrix for dirac 
		indices. Heechang*/
		for(a=0;a<3;a++)for(b=0;b<3;b++)for(si=0;si<4;si++)for(sf=0;sf<4;sf++){
			if(si == sf)light.c[a].d[si].d[sf].c[b] = ks_1st_quark->e[a][b];
			else light.c[a].d[si].d[sf].c[b] =  zero;
		}
		/*Now, multiply by Omega matrix */
		if((my_z % 2) == 1)
			for(a=0;a<3;a++){
				mult_swv_by_gamma_l( &(light.c[a]), &(tmp.c[a]), ZUP);
				light.c[a] = tmp.c[a]; 
			}  

		if((my_y % 2) == 1)
			for(a=0;a<3;a++){
				mult_swv_by_gamma_l( &(light.c[a]), &(tmp.c[a]), YUP);
				light.c[a] = tmp.c[a];
			}  

		if((my_x % 2) == 1)
			for(a=0;a<3;a++){
				mult_swv_by_gamma_l( &(light.c[a]), &(tmp.c[a]), XUP);      
				light.c[a] = tmp.c[a]; 
			}  

		if((t % 2) == 1)
			for(a=0;a<3;a++){
				mult_swv_by_gamma_l( &(light.c[a]), &(tmp.c[a]), TUP);
				light.c[a] = tmp.c[a]; 
			}

		/* Start the mu and nu loops */
		for(mu=0;mu<3;mu++)for(nu=0;nu<3;nu++){

			/*initialize propagator buffer*/
			/*accumulate at propagator */
			for(si=0;si<4;si++)for(sf=0;sf<4;sf++){
				prop_tot[si][sf].real = 0.0;
				prop_tot[si][sf].imag = 0.0;
			}



			/*Calculate the Trace part*/
			/* Trace = tr[ ( g_nu g_2 g_0 G^T ) ( g_0 g_2 g_mu G ) ] */

			spin_trace_two_clover_double(quark1,mu,nu,tr);


			/*product of three propagators with summation of colors */
			for(a=0;a<3;a++){ ap = a;
				for(j=-1;j<=1;j=j+2){
					b = bp = (a+j+3)%3;
					c = cp = (a-j+3)%3;
					eps = 1;

					CMULREAL(tr[a][ap][b][bp],eps,tr[a][ap][b][bp]);
					for(si=0;si<4;si++)for(sf=0;sf<4;sf++){

						CMUL(tr[a][ap][b][bp],light.c[c].d[si].d[sf].c[cp],prop_matrix);

						CSUM(prop_tot[si][sf],prop_matrix);
					}

				}
			}

			for(a=0;a<3;a++){ ap = a;
				for(j=-1;j<=1;j=j+2){
					b = cp = (a+j+3)%3;
					c = bp = (a-j+3)%3;
					eps = -1;

					CMULREAL(tr[a][ap][b][bp],eps,tr[a][ap][b][bp]);
					for(si=0;si<4;si++)for(sf=0;sf<4;sf++){

						CMUL(tr[a][ap][b][bp],light.c[c].d[si].d[sf].c[cp],prop_matrix);

						CSUM(prop_tot[si][sf],prop_matrix);

					}

				}
			}

			for(a=0;a<3;a++){ bp = a;
				for(j=-1;j<=1;j=j+2){
					b = ap = (a+j+3)%3;
					c = cp = (a-j+3)%3;
					eps = -1;

					CMULREAL(tr[a][ap][b][bp],eps,tr[a][ap][b][bp]);
					for(si=0;si<4;si++)for(sf=0;sf<4;sf++){

						CMUL(tr[a][ap][b][bp],light.c[c].d[si].d[sf].c[cp],prop_matrix);

						CSUM(prop_tot[si][sf],prop_matrix);

					}

				}
			}

			for(a=0;a<3;a++){ bp = a;
				for(j=-1;j<=1;j=j+2){
					b = cp = (a+j+3)%3;
					c = ap = (a-j+3)%3;
					eps = 1;

					CMULREAL(tr[a][ap][b][bp],eps,tr[a][ap][b][bp]);
					for(si=0;si<4;si++)for(sf=0;sf<4;sf++){

						CMUL(tr[a][ap][b][bp],light.c[c].d[si].d[sf].c[cp],prop_matrix);

						CSUM(prop_tot[si][sf],prop_matrix);

					}

				}
			}


			for(a=0;a<3;a++){ cp = a;
				for(j=-1;j<=1;j=j+2){
					b = bp = (a+j+3)%3;
					c = ap = (a-j+3)%3;
					eps = -1;

					CMULREAL(tr[a][ap][b][bp],eps,tr[a][ap][b][bp]);
					for(si=0;si<4;si++)for(sf=0;sf<4;sf++){

						CMUL(tr[a][ap][b][bp],light.c[c].d[si].d[sf].c[cp],prop_matrix);

						CSUM(prop_tot[si][sf],prop_matrix);

					}

				}
			}

			for(a=0;a<3;a++){ cp = a;
				for(j=-1;j<=1;j=j+2){
					b = ap = (a+j+3)%3;
					c = bp = (a-j+3)%3;
					eps = 1;

					CMULREAL(tr[a][ap][b][bp],eps,tr[a][ap][b][bp]);
					for(si=0;si<4;si++)for(sf=0;sf<4;sf++){

						CMUL(tr[a][ap][b][bp],light.c[c].d[si].d[sf].c[cp],prop_matrix);

						CSUM(prop_tot[si][sf],prop_matrix);

					}

				}
			}


			/* exp(i*vector_p.vector_x) */
			epx1.real = 0.0;
			epx1.imag = (mom[0]*(double)p[0]*(double)my_x + mom[1]*(double)p[1]*(double)my_y + mom[2]*(double)p[2]*(double)my_z);
			epx = dcexp(&epx1);


			for(si=0;si<4;si++)for(sf=0;sf<4;sf++){
				CMUL(epx, prop_tot[si][sf], prop_all);
				CSUM(propagator[si][sf][mu][nu][t],prop_all);
			}
		}/*mu and nu loops*/
	}/*foralllattice loop*/

	return 0;
}



/*for the Omu operator in singly heavy baryons */
int ks_baryon_2point_Omu(field_offset ks_1, field_offset ks_2, field_offset heavy_quark, double_complex *propagator[4][4][3][3])
{
        int i, j, t, my_x, my_y, my_z, mu, my, nu;
	site *s;
	double_complex epx, epx1, EP[3][3];
	double pi, mom[3];
	int eps;
	int a, b, c, ap, bp, cp; /*color indices. bp -> b prime */
	int si,sf; /*dirac indices for heavy quark propagator */

	su3_matrix *ks_1st_quark;
	su3_matrix *ks_2nd_quark;
	wilson_propagator *quark;
	double_complex prop;
	double_complex prop_matrix[4][4], prop_tot[4][4], prop_all[4][4][3][3];

		int p[3] ={0,0,0};

	pi = 4.0 * atan( 1.);
	mom[0] = -2.*pi/(double)nx;
	mom[1] = -2.*pi/(double)ny;
	mom[2] = -2.*pi/(double)nz;

	FORALLSITES(i,s){
		t = s->t;
		my_x = s->x;
		my_y = s->y;
		my_z = s->z;

		/*initialize propagator buffer*/
		/*accumulate at propagator */
		for(si=0;si<4;si++)for(sf=0;sf<4;sf++){
			prop_tot[si][sf].real = 0.0;
			prop_tot[si][sf].imag = 0.0;
		}
		/* asign each field */
		ks_1st_quark = (su3_matrix *)F_PT(s, ks_1);
		ks_2nd_quark = (su3_matrix *)F_PT(s, ks_2);
		quark = (wilson_propagator *)F_PT(s, heavy_quark);

		/*product of three propagators with sumation of colors */
		for(a=0;a<3;a++){ ap = a;
			for(j=-1;j<=1;j=j+2){
				b = bp = (a+j+3)%3;
				c = cp = (a-j+3)%3;
				eps = 1;

				CMUL(ks_1st_quark->e[a][ap],ks_2nd_quark->e[b][bp],prop);
				prop.real *= eps;
				prop.imag *= eps;
				for(si=0;si<4;si++)for(sf=0;sf<4;sf++){

					CMUL(prop,quark->c[c].d[si].d[sf].c[cp],prop_matrix[si][sf]);

					prop_tot[si][sf].real += prop_matrix[si][sf].real;
					prop_tot[si][sf].imag += prop_matrix[si][sf].imag;
				}

			}
		}

		for(a=0;a<3;a++){ ap = a;
			for(j=-1;j<=1;j=j+2){
				b = cp = (a+j+3)%3;
				c = bp = (a-j+3)%3;
				eps = -1;

				CMUL(ks_1st_quark->e[a][ap],ks_2nd_quark->e[b][bp],prop);
				prop.real *= eps;
				prop.imag *= eps;
				for(si=0;si<4;si++)for(sf=0;sf<4;sf++){

					CMUL(prop,quark->c[c].d[si].d[sf].c[cp],prop_matrix[si][sf]);

					prop_tot[si][sf].real += prop_matrix[si][sf].real;
					prop_tot[si][sf].imag += prop_matrix[si][sf].imag;
				}

			}
		}

		for(a=0;a<3;a++){ bp = a;
			for(j=-1;j<=1;j=j+2){
				b = ap = (a+j+3)%3;
				c = cp = (a-j+3)%3;
				eps = -1;

				CMUL(ks_1st_quark->e[a][ap],ks_2nd_quark->e[b][bp],prop);
				prop.real *= eps;
				prop.imag *= eps;
				for(si=0;si<4;si++)for(sf=0;sf<4;sf++){

					CMUL(prop,quark->c[c].d[si].d[sf].c[cp],prop_matrix[si][sf]);

					prop_tot[si][sf].real += prop_matrix[si][sf].real;
					prop_tot[si][sf].imag += prop_matrix[si][sf].imag;
				}

			}
		}

		for(a=0;a<3;a++){ bp = a;
			for(j=-1;j<=1;j=j+2){
				b = cp = (a+j+3)%3;
				c = ap = (a-j+3)%3;
				eps = 1;

				CMUL(ks_1st_quark->e[a][ap],ks_2nd_quark->e[b][bp],prop);
				prop.real *= eps;
				prop.imag *= eps;
				for(si=0;si<4;si++)for(sf=0;sf<4;sf++){

					CMUL(prop,quark->c[c].d[si].d[sf].c[cp],prop_matrix[si][sf]);

					prop_tot[si][sf].real += prop_matrix[si][sf].real;
					prop_tot[si][sf].imag += prop_matrix[si][sf].imag;
				}

			}
		}


		for(a=0;a<3;a++){ cp = a;
			for(j=-1;j<=1;j=j+2){
				b = bp = (a+j+3)%3;
				c = ap = (a-j+3)%3;
				eps = -1;

				CMUL(ks_1st_quark->e[a][ap],ks_2nd_quark->e[b][bp],prop);
				prop.real *= eps;
				prop.imag *= eps;
				for(si=0;si<4;si++)for(sf=0;sf<4;sf++){

					CMUL(prop,quark->c[c].d[si].d[sf].c[cp],prop_matrix[si][sf]);

					prop_tot[si][sf].real += prop_matrix[si][sf].real;
					prop_tot[si][sf].imag += prop_matrix[si][sf].imag;
				}

			}
		}

		for(a=0;a<3;a++){ cp = a;
			for(j=-1;j<=1;j=j+2){
				b = ap = (a+j+3)%3;
				c = bp = (a-j+3)%3;
				eps = 1;

				CMUL(ks_1st_quark->e[a][ap],ks_2nd_quark->e[b][bp],prop);
				prop.real *= eps;
				prop.imag *= eps;
				for(si=0;si<4;si++)for(sf=0;sf<4;sf++){

					CMUL(prop,quark->c[c].d[si].d[sf].c[cp],prop_matrix[si][sf]);

					prop_tot[si][sf].real += prop_matrix[si][sf].real;
					prop_tot[si][sf].imag += prop_matrix[si][sf].imag;
				}

			}
		}


		/* exp(i*vector_p.vector_x) */
		epx1.real = 0.0;
		epx1.imag = (mom[0]*(double)p[0]*(double)my_x + mom[1]*(double)p[1]*(double)my_y
				+ mom[2]*(double)p[2]*(double)my_z);
		epx = dcexp(&epx1);

		/* 4*(-1)^(x_mu) factor, mu is only for x,y, and z, so 0 to 2 */
		for(mu=0;mu<3;mu++){
			nu = mu;
			if(mu==0) my = my_x;
			if(mu==1) my = my_y;
			if(mu==2) my = my_z;
			if((my)%2 == 1) {
				EP[mu][nu].real = -4.0*epx.real;
				EP[mu][nu].imag = -4.0*epx.imag;
			}
			else{
				EP[mu][nu].real = 4.0*epx.real;
				EP[mu][nu].imag = 4.0*epx.imag;
			}/*In case of Omu */


			for(si=0;si<4;si++)for(sf=0;sf<4;sf++){
				CMUL(EP[mu][nu], prop_tot[si][sf], prop_all[si][sf][mu][nu]);
				propagator[si][sf][mu][nu][t].real += prop_all[si][sf][mu][nu].real;
				propagator[si][sf][mu][nu][t].imag += prop_all[si][sf][mu][nu].imag;
			}
		}
	}/*foralllattice loop*/

	return 0;
}


/*for O5 operator */
int ks_baryon_2point_O5(field_offset ks_1, field_offset ks_2, field_offset heavy_quark, double_complex *propagator[4][4])
{
	int i, j, t, my_x, my_y, my_z;
	site *s;
	double_complex epx, epx1;
	double pi, mom[3];
	int eps;
	int a, b, c, ap, bp, cp; /*color indices. bp -> b prime */
	int si,sf; /*dirac indices for heavy quark propagator */

	su3_matrix *ks_1st_quark;
	su3_matrix *ks_2nd_quark;
	wilson_propagator *quark;
	double_complex prop;
	double_complex prop_matrix[4][4], prop_tot[4][4];

		int p[3] ={0,0,0};

	pi = 4.0 * atan( 1.);
	mom[0] = -2.*pi/(double)nx;
	mom[1] = -2.*pi/(double)ny;
	mom[2] = -2.*pi/(double)nz;

	FORALLSITES(i,s){
		t = s->t;
		my_x = s->x;
		my_y = s->y;
		my_z = s->z;

		/*initialize propagator buffer*/
		/*accumulate at propagator */
		for(si=0;si<4;si++)for(sf=0;sf<4;sf++){
			prop_tot[si][sf].real = 0.0;
			prop_tot[si][sf].imag = 0.0;
		}
		/* asign each field */
		ks_1st_quark = (su3_matrix *)F_PT(s, ks_1);
		ks_2nd_quark = (su3_matrix *)F_PT(s, ks_2);
		quark = (wilson_propagator *)F_PT(s, heavy_quark);

		/*product of three propagators with sumation of colors */
		for(a=0;a<3;a++){ ap = a;
			for(j=-1;j<=1;j=j+2){
				b = bp = (a+j+3)%3;
				c = cp = (a-j+3)%3;
				eps = 1;

				CMUL(ks_1st_quark->e[a][ap],ks_2nd_quark->e[b][bp],prop);
				prop.real *= eps;
				prop.imag *= eps;
				for(si=0;si<4;si++)for(sf=0;sf<4;sf++){

					CMUL(prop,quark->c[c].d[si].d[sf].c[cp],prop_matrix[si][sf]);

					prop_tot[si][sf].real += prop_matrix[si][sf].real;
					prop_tot[si][sf].imag += prop_matrix[si][sf].imag;
				}

			}
		}

		for(a=0;a<3;a++){ ap = a;
			for(j=-1;j<=1;j=j+2){
				b = cp = (a+j+3)%3;
				c = bp = (a-j+3)%3;
				eps = -1;

				CMUL(ks_1st_quark->e[a][ap],ks_2nd_quark->e[b][bp],prop);
				prop.real *= eps;
				prop.imag *= eps;
				for(si=0;si<4;si++)for(sf=0;sf<4;sf++){

					CMUL(prop,quark->c[c].d[si].d[sf].c[cp],prop_matrix[si][sf]);

					prop_tot[si][sf].real += prop_matrix[si][sf].real;
					prop_tot[si][sf].imag += prop_matrix[si][sf].imag;
				}

			}
		}

		for(a=0;a<3;a++){ bp = a;
			for(j=-1;j<=1;j=j+2){
				b = ap = (a+j+3)%3;
				c = cp = (a-j+3)%3;
				eps = -1;

				CMUL(ks_1st_quark->e[a][ap],ks_2nd_quark->e[b][bp],prop);
				prop.real *= eps;
				prop.imag *= eps;
				for(si=0;si<4;si++)for(sf=0;sf<4;sf++){

					CMUL(prop,quark->c[c].d[si].d[sf].c[cp],prop_matrix[si][sf]);

					prop_tot[si][sf].real += prop_matrix[si][sf].real;
					prop_tot[si][sf].imag += prop_matrix[si][sf].imag;
				}

			}
		}

		for(a=0;a<3;a++){ bp = a;
			for(j=-1;j<=1;j=j+2){
				b = cp = (a+j+3)%3;
				c = ap = (a-j+3)%3;
				eps = 1;

				CMUL(ks_1st_quark->e[a][ap],ks_2nd_quark->e[b][bp],prop);
				prop.real *= eps;
				prop.imag *= eps;
				for(si=0;si<4;si++)for(sf=0;sf<4;sf++){

					CMUL(prop,quark->c[c].d[si].d[sf].c[cp],prop_matrix[si][sf]);

					prop_tot[si][sf].real += prop_matrix[si][sf].real;
					prop_tot[si][sf].imag += prop_matrix[si][sf].imag;
				}

			}
		}


		for(a=0;a<3;a++){ cp = a;
			for(j=-1;j<=1;j=j+2){
				b = bp = (a+j+3)%3;
				c = ap = (a-j+3)%3;
				eps = -1;

				CMUL(ks_1st_quark->e[a][ap],ks_2nd_quark->e[b][bp],prop);
				prop.real *= eps;
				prop.imag *= eps;
				for(si=0;si<4;si++)for(sf=0;sf<4;sf++){

					CMUL(prop,quark->c[c].d[si].d[sf].c[cp],prop_matrix[si][sf]);

					prop_tot[si][sf].real += prop_matrix[si][sf].real;
					prop_tot[si][sf].imag += prop_matrix[si][sf].imag;
				}

			}
		}

		for(a=0;a<3;a++){ cp = a;
			for(j=-1;j<=1;j=j+2){
				b = ap = (a+j+3)%3;
				c = bp = (a-j+3)%3;
				eps = 1;

				CMUL(ks_1st_quark->e[a][ap],ks_2nd_quark->e[b][bp],prop);
				prop.real *= eps;
				prop.imag *= eps;
				for(si=0;si<4;si++)for(sf=0;sf<4;sf++){

					CMUL(prop,quark->c[c].d[si].d[sf].c[cp],prop_matrix[si][sf]);

					prop_tot[si][sf].real += prop_matrix[si][sf].real;
					prop_tot[si][sf].imag += prop_matrix[si][sf].imag;
				}

			}
		}


		/* exp(i*vector_p.vector_x) */
		epx1.real = 0.0;
		epx1.imag = (mom[0]*(double)p[0]*(double)my_x + mom[1]*(double)p[1]*(double)my_y
				+ mom[2]*(double)p[2]*(double)my_z);
		epx = dcexp(&epx1);

		/* 4*(-1)^(...) factor */
		/* epx has information of exponential and 4*(-1)^(..) */

		if((1)%2 == 1) {
			epx.real *= 4.0;
			epx.imag *= 4.0;
		}
		else{
			epx.real *= -4.0;
			epx.imag *= -4.0;
		}   /*In case of O_5 */


		/*if((my_x)%2 == 1) {
		  epx.real *= -4.0;
		  epx.imag *= -4.0;
		  }
		  else{
		  epx.real *= 4.0;
		  epx.imag *= 4.0;
		  }*/ /*In case of O_mu */


		for(si=0;si<4;si++)for(sf=0;sf<4;sf++){
			CMUL(epx, prop_tot[si][sf], prop_tot[si][sf]);
			propagator[si][sf][t].real += prop_tot[si][sf].real;
			propagator[si][sf][t].imag += prop_tot[si][sf].imag;
		}
	}/*foralllattice loop*/

	return 0;
}



