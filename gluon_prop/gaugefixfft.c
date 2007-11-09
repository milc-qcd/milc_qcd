/************************** gaugefixfft.c *******************************/
/* Fix Coulomb or Lorentz gauge using Fourier acceleration */
/* Uses double precision global sums */
/* MIMD version 7 */
/* U.M. Heller 12-29-00 */

/* Prototype...

   void gaugefixfft(int gauge_dir, Real accel_param, int max_gauge_iter,
	      Real gauge_fix_tol,
	      int nvector, field_offset vector_offset[], int vector_parity[],
	      int nantiherm, field_offset antiherm_offset[], 
	      int antiherm_parity[] )
   -------------------------------------------------------------------

   NOTE: We assume that the FFT setup was already done!

   -------------------------------------------------------------------

   NOTE: For staggered fermion applications, it is necessary to remove
   the KS phases from the gauge links before calling this procedure.
   See "rephase" in setup.c.

   -------------------------------------------------------------------
   EXAMPLE:  Fixing only the link matrices to Coulomb gauge

   gaugefixfft(TUP,(Real)0.065,500,(Real)1.0e-7,
	       0,NULL,NULL,0,NULL,NULL);

   -------------------------------------------------------------------
   EXAMPLE:  Fixing Coulomb gauge with respect to the y direction
      in the staggered fermion scheme and simultaneously transforming
      the pseudofermion fields and gauge-momenta involved in updating:

   int nvector = 3;
   field_offset vector_offset[3] = { F_OFFSET(g_rand), F_OFFSET(phi), 
        F_OFFSET(xxx) };
   int vector_parity[3] = { EVENANDODD, EVEN, EVEN };
   int nantiherm = 4;
   field_offset antiherm_offset[4] = { F_OFFSET(mom[0]), F_OFFSET(mom[1]),
       F_OFFSET(mom[2]), F_OFFSET(mom[3]) };
   field_offset antiherm_parity[4] = { EVENANDODD, EVENANDODD, EVENANDODD,
       EVENANDODD }

   rephase( OFF );
   gaugefixfft(YUP,(Real)1.8,500,(Real)2.0e-6,
       nvector,vector_offset,vector_parity,
       nantiherm,antiherm_offset,antiherm_parity);
   rephase( ON );

   -------------------------------------------------------------------

   gauge_dir     specifies the direction of the "time"-like hyperplane
                 for the purposes of defining Coulomb or Lorentz gauge 
      TUP    for evaluating propagators in the time-like direction
      ZUP    for screening lengths.
      8      for Lorentz gauge
   accel_param	   Parameter for Fourier acceleration
   max_gauge_iter  Maximum number of iterations 
   gauge_fix_tol   Stop if change is less than this 
   */

#include "gluon_prop_includes.h"
#define REUNIT_INTERVAL 20

/* Scratch space */

su3_matrix *diffmatp;	/* malloced diffmat pointer */
su3_matrix *tmpmat1p;	/* malloced tmpmat1 pointer */
su3_matrix *tmpmat2p;	/* malloced tmpmat2 pointer */
Real *p_rat;		/* malloced momentum ratio for Fourier accelaration */
su3_matrix *gt_matrix;	/* malloced gauge transformation matrix, if
			   gauge transformation on vectors, etc. are needed */
int fft_vol;

/* Prototypes */

double dmu_amu(int gauge_dir);
void gf_fft_step(int gauge_dir, double *dmu_amu_norm, Real accel_param,
		 int save_gt);
void gaugefixsetup(int gauge_dir, int save_gt);


double dmu_amu(int gauge_dir)
{

/* Accumulates d_mu A_mu in diffmat and computes the norm */

register int dir, i;
register site *s;
msg_tag *mtag[6];
double damu;
su3_matrix *m1;
anti_hermitmat ahtmp;

    /* Start gathers of downward links */
    FORALLUPDIRBUT(gauge_dir,dir){
	mtag[dir] = start_gather_site( F_OFFSET(link[dir]), sizeof(su3_matrix),
				       OPP_DIR(dir), EVENANDODD, gen_pt[dir] );
    }

    /* Clear diffmat */
    FORALLSITES(i,s){
	clear_su3mat(&diffmatp[i]);
    }
  
    /* Add upward link contributions */
    FORALLSITES(i,s)
    {
	FORALLUPDIRBUT(gauge_dir,dir){
	    add_su3_matrix( &diffmatp[i], &(s->link[dir]), &diffmatp[i]);
	}
    }

    FORALLUPDIRBUT(gauge_dir,dir){
	wait_gather(mtag[dir]);
    }

    /* Subtract downward link contributions, make traceless anti-hermitian
       and accumulate norm. Note we store it as SU(3) matrix for
       calls to FFT */
    damu = 0.0;
    FORALLSITES(i,s)
    {
	m1 = &diffmatp[i];
	FORALLUPDIRBUT(gauge_dir,dir){
	    sub_su3_matrix( m1, (su3_matrix *)gen_pt[dir][i], m1);
	}
	make_anti_hermitian( m1, &ahtmp);
	uncompress_anti_hermitian( &ahtmp, m1);
#ifndef IMP_GFIX
	damu += (double)realtrace_su3( m1, m1);
#endif
    }

#ifdef IMP_GFIX
    /* Clear tmpmat1, for the 2-link term */
    FORALLSITES(i,s){
	clear_su3mat( &tmpmat1p[i]);
    }

    FORALLUPDIRBUT(gauge_dir,dir){
	/* Double link U(mu,x-mu) U(mu,x) */
	FORALLSITES(i,s)
	{
	    mult_su3_nn( (su3_matrix *)gen_pt[dir][i],  &(s->link[dir]),
			 &tmpmat2p[i]);
	}

	mtag[4] = start_gather_field( tmpmat2p, sizeof(su3_matrix),
				     dir, EVENANDODD, gen_pt[4] );
	mtag[5] = start_gather_field( tmpmat2p, sizeof(su3_matrix),
				     OPP_DIR(dir), EVENANDODD, gen_pt[5] );
	wait_gather(mtag[4]);
	wait_gather(mtag[5]);

	FORALLSITES(i,s)
	{
	    add_su3_matrix( &tmpmat1p[i], (su3_matrix *)gen_pt[4][i],
			    &tmpmat1p[i]);
	    sub_su3_matrix( &tmpmat1p[i], (su3_matrix *)gen_pt[5][i],
			    &tmpmat1p[i]);
	}

	cleanup_gather(mtag[4]);
	cleanup_gather(mtag[5]);
    }

    /* Now compute the improved d_mu A_mu, and its norm^2 */
    ftmp1 = 4.0/3.0;
    ftmp2 = - 1.0/(12.0*u0);
    FORALLSITES(i,s)
    {
	m1 = &diffmatp[i];
	scalar_mult_su3_matrix( m1, ftmp1, m1);
	scalar_mult_add_su3_matrix( m1, &tmpmat1p[i], ftmp2, m1);
	make_anti_hermitian( m1, &ahtmp);
	uncompress_anti_hermitian( &ahtmp, m1);
	damu += (double)realtrace_su3( m1, m1);
    }
#endif

    FORALLUPDIRBUT(gauge_dir,dir){
	cleanup_gather(mtag[dir]);
    }

    g_doublesum( &damu);

    /* Normalize, with factor "no. of active dirs" * (Nc^2-1) / 2 */
    if (gauge_dir > TUP){
	damu /= (double)(16*volume);
    }
    else{
	damu /= (double)(12*volume);
    }

    return(sqrt(damu));
} /* dmu_amu */


void gf_fft_step(int gauge_dir, double *dmu_amu_norm, Real accel_param,
		 int save_gt)
{
/* Carry out one iteration in the gauge-fixing process */

msg_tag *mtag[4];
register int dir, i, j;
register site *s;
su3_matrix *matp1, *matp2, tmat1, tmat2;
Real ftmp1, ftmp2;
int err;

    /* Compute d_mu A_mu(x) in diffmat and its norm */
    *dmu_amu_norm = dmu_amu(gauge_dir);

    /* Fourier transform d_mu A_mu(x) */
    g_sync();
    restrict_fourier_field((complex *)diffmatp, sizeof(su3_matrix), FORWARDS);

    /* Multiply with (accel_param * p_rat / fft_vol) */
    ftmp1 = accel_param / (Real)fft_vol;
    FORALLSITES(i,s){
	ftmp2 = ftmp1 * p_rat[i];
	scalar_mult_su3_matrix( &diffmatp[i], ftmp2, &diffmatp[i]);
    }

    /* Fourier transform back */
    g_sync();
    restrict_fourier_field((complex *)diffmatp, sizeof(su3_matrix), BACKWARDS);

    /* Now exponentiate, using 6-th order expansion */
    FORALLSITES(i,s){
	matp1 = &diffmatp[i];
	matp2 = &tmpmat1p[i];
	clear_su3mat( matp2);
	for(j=0; j<3; j++) matp2->e[j][j].real = 1.0;
	scalar_mult_add_su3_matrix( matp2, matp1, 0.16666666, &tmat2);
	mult_su3_nn( matp1, &tmat2, &tmat1);
	scalar_mult_add_su3_matrix( matp2, &tmat1, 0.2, &tmat2);
	mult_su3_nn( matp1, &tmat2, &tmat1);
	scalar_mult_add_su3_matrix( matp2, &tmat1, 0.25, &tmat2);
	mult_su3_nn( matp1, &tmat2, &tmat1);
	scalar_mult_add_su3_matrix( matp2, &tmat1, 0.33333333, &tmat2);
	mult_su3_nn( matp1, &tmat2, &tmat1);
	scalar_mult_add_su3_matrix( matp2, &tmat1, 0.5, &tmat2);
	mult_su3_nn( matp1, &tmat2, &tmat1);
	add_su3_matrix( matp2, &tmat1, &tmat2);
	su3mat_copy( &tmat2, matp2);
	/* reunitarize */
	err =  reunit_su3( matp2);
	if(err > 0) printf("Node %d site %d: G nut unitary, err = %d\n",
			    this_node, i, err);
    }

    /* The gauge matrices are now in tmpmat1p */
    /* Do the gauge transformation of the links */
    g_sync();
    FORALLUPDIR(dir)
	mtag[dir] = start_gather_field( tmpmat1p, sizeof(su3_matrix),
				  dir, EVENANDODD, gen_pt[dir] );

    FORALLUPDIR(dir){
	/* First multiply with the gauge matrices on site */
	FORALLSITES(i,s){
	    mult_su3_nn( &tmpmat1p[i], &(s->link[dir]), &tmpmat2p[i]);
	}

	/* Then multiply with forward gauge matrices */
	wait_gather(mtag[dir]);
	FORALLSITES(i,s){
	    mult_su3_na( &tmpmat2p[i], (su3_matrix *)gen_pt[dir][i],
			 &(s->link[dir]));
	}
	cleanup_gather(mtag[dir]);
    }

    if (save_gt == 1){
	/* Also multiply with existing gauge transfromation matrices */
	FORALLSITES(i,s){
	    matp1 = &gt_matrix[i];
	    mult_su3_nn( &tmpmat1p[i], matp1, &tmat1);
	    su3mat_copy( &tmat1, matp1);
	}

    }

} /* gf_fft_step */


void gaugefixsetup(int gauge_dir, int save_gt)
{
register int dir, i, j, pmu;
register site *s;
Real pix, piy, piz, pit;
Real sin_pmu, sum_p2, sum_p2_max;

    diffmatp = (su3_matrix *)malloc(sizeof(su3_matrix)*sites_on_node);
    if(diffmatp == NULL){
	node0_printf("gaugefix: Can't malloc diffmat\n");
	fflush(stdout);terminate(1);
    }
    tmpmat1p = (su3_matrix *)malloc(sizeof(su3_matrix)*sites_on_node);
    if(tmpmat1p == NULL){
	node0_printf("gaugefix: Can't malloc tmpmat1\n");
	fflush(stdout);terminate(1);
    }
    tmpmat2p = (su3_matrix *)malloc(sizeof(su3_matrix)*sites_on_node);
    if(tmpmat2p == NULL){
	node0_printf("gaugefix: Can't malloc tmpmat2\n");
	fflush(stdout);terminate(1);
    }

    fft_vol = 1;
    FORALLUPDIRBUT(gauge_dir,dir){
	switch(dir){
	    case XUP: fft_vol *= nx; break;
	    case YUP: fft_vol *= ny; break;
	    case ZUP: fft_vol *= nz; break;
	    case TUP: fft_vol *= nt; break;
	}
    }

    pix = PI / (Real)nx;
    piy = PI / (Real)ny;
    piz = PI / (Real)nz;
    pit = PI / (Real)nt;

    /* Allocate space for p_rat = (hat)p^2_max / (hat)p^2 */
    p_rat = (Real *)malloc(sizeof(Real)*sites_on_node);
    if(p_rat == NULL){
	node0_printf("gaugefixfft: Can't malloc p_rat\n");
	fflush(stdout);terminate(1);
    }

    sum_p2_max = 0.0;
    /* Now compute sum_mu sin^2(p_mu/2) */
    FORALLSITES(i,s){
	sum_p2 = 0.0;
	FORALLUPDIRBUT(gauge_dir,dir){
	    switch(dir){
		case XUP:
		    pmu = s->x;
		    sin_pmu = sin((double)(pmu*pix));
		    sum_p2 += sin_pmu * sin_pmu;
		    break;

		case YUP:
		    pmu = s->y;
		    sin_pmu = sin((double)(pmu*piy));
		    sum_p2 += sin_pmu * sin_pmu;
		    break;

		case ZUP:
		    pmu = s->z;
		    sin_pmu = sin((double)(pmu*piz));
		    sum_p2 += sin_pmu * sin_pmu;
		    break;

		case TUP:
		    pmu = s->t;
		    sin_pmu = sin((double)(pmu*pit));
		    sum_p2 += sin_pmu * sin_pmu;
		    break;

                default: printf("BOTCH: bad direction\n"); exit(1);
            }
        }

	p_rat[i] = sum_p2;
	if (sum_p2 > sum_p2_max) sum_p2_max = sum_p2;
    }
    g_floatmax( &sum_p2_max);

    if (gauge_dir > TUP ){
	FORALLSITES(i,s){
	    if( s->x == 0 && s->y == 0 && s->z == 0 && s->t == 0 ){
		p_rat[i] = 0.0;
	    }
	    else{
		p_rat[i] = sum_p2_max / p_rat[i];
	    }
	}
    }
    else{
	switch(gauge_dir){
	    case XUP:
		FORALLSITES(i,s){
		    if( s->y == 0 && s->z == 0 && s->t == 0 ){
			p_rat[i] = 0.0;
		    }
		    else{
			p_rat[i] = sum_p2_max / p_rat[i];
		    }
		}
		break;

	    case YUP:
		FORALLSITES(i,s){
		    if( s->x == 0 && s->z == 0 && s->t == 0 ){
			p_rat[i] = 0.0;
		    }
		    else{
			p_rat[i] = sum_p2_max / p_rat[i];
		    }
		}
		break;

	    case ZUP:
		FORALLSITES(i,s){
		    if( s->x == 0 && s->y == 0 && s->t == 0 ){
			p_rat[i] = 0.0;
		    }
		    else{
			p_rat[i] = sum_p2_max / p_rat[i];
		    }
		}
		break;

	    case TUP:
		FORALLSITES(i,s){
		    if( s->x == 0 && s->y == 0 && s->z == 0 ){
			p_rat[i] = 0.0;
		    }
		    else{
			p_rat[i] = sum_p2_max / p_rat[i];
		    }
		}
		break;

	    default: printf("BOTCH: bad direction\n"); exit(1);
	}
    }

    /* Allocate space for gauge transformation, if it is needed */
    if (save_gt == 1){
	gt_matrix = (su3_matrix *)malloc(sizeof(su3_matrix)*sites_on_node);
	if(gt_matrix == NULL){
	    node0_printf("gaugefixfft: Can't malloc gt_matrix\n");
	    fflush(stdout);terminate(1);
	}

	/* Initialize it to unity */
	FORALLSITES(i,s){
	    clear_su3mat( &gt_matrix[i]);
	    for(j=0; j<3; j++) gt_matrix[i].e[j][j].real = 1.0;
	}
    }
 
} /* gaugefixsetup */

void gaugefixfft_combo(int gauge_dir, Real accel_param, int max_gauge_iter,
		       Real gauge_fix_tol, int nvector,
		       field_offset vector_offset[], int vector_parity[],
		       int nantiherm, field_offset antiherm_offset[], 
		       int antiherm_parity[] )
{
register int i;
register site *s;
int gauge_iter, save_gt, j;
double curr_damu;
su3_vector vtmp;
su3_matrix htmp1, htmp2;

    if(nvector > 0 || nantiherm > 0){
	save_gt = 1;
    }
    else{
	save_gt = 0;
    }

    /* Set up work space */
    gaugefixsetup(gauge_dir, save_gt);

    /* Do at most max_gauge_iter iterations, but stop after the first step */
    /* if |d_mu A_mu| is smaller than gauge_fix_tol */

    for (gauge_iter=0; gauge_iter < max_gauge_iter; gauge_iter++)
    {
	gf_fft_step(gauge_dir, &curr_damu, accel_param, save_gt);

if(gauge_iter==0) node0_printf("Starting |d_mu A_mu| = %e\n", curr_damu);

	if (curr_damu < gauge_fix_tol) break;

	/* Reunitarize when iteration count is a multiple of REUNIT_INTERVAL */
	if((gauge_iter % REUNIT_INTERVAL) == (REUNIT_INTERVAL - 1)){
	    node0_printf("step %d: |d_mu A_mu| = %e\n", gauge_iter, curr_damu);
	    reunitarize();
	    fflush(stdout);
	}
    }
    /* Reunitarize at the end, unless we just did it in the loop */
    if((gauge_iter % REUNIT_INTERVAL) != 0)
	reunitarize();

    /* Transform vectors and gauge momenta if requested */
    for(j = 0; j < nvector; j++)
	FORSOMEPARITY(i,s,vector_parity[j]){
	    mult_su3_mat_vec( &gt_matrix[i],
		(su3_vector *)F_PT(s,vector_offset[j]), &vtmp);
	    su3vec_copy( &vtmp, (su3_vector *)F_PT(s,vector_offset[j]));
	}

    for(j = 0; j < nantiherm; j++)
	FORSOMEPARITY(i,s,antiherm_parity[j]){
	    uncompress_anti_hermitian(
		(anti_hermitmat *)F_PT(s,antiherm_offset[j]), &htmp1);
	    mult_su3_nn( &gt_matrix[i], &htmp1, &htmp2);
	    mult_su3_na( &htmp2, &gt_matrix[i], &htmp1);
	    make_anti_hermitian( &htmp1,
		(anti_hermitmat *)F_PT(s,antiherm_offset[j]));
	}


    /* Free workspace */
    free(diffmatp);
    free(tmpmat1p);
    free(tmpmat2p);
    free(p_rat);
    if(save_gt == 1) free(gt_matrix);
  
    if(this_node==0)
	printf("GFIX: Ended at step %d. |d_mu A_mu| = %e\n",
	       gauge_iter, curr_damu);
}

/* Abbreviated form for fixing only gauge field */

void gaugefixfft(int gauge_dir, Real accel_param, int max_gauge_iter,
		 Real gauge_fix_tol)
{
    gaugefixfft_combo(gauge_dir, accel_param, max_gauge_iter, gauge_fix_tol,
		      0,NULL,NULL,0,NULL,NULL);
}
