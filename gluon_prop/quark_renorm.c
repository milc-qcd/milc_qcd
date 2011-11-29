/******** quark_renorm.c *********/
/* MIMD version 7 */
/* Initial version by UMH 4/1/06 */
/* Converted from version 6 to version 7 by UMH 8/12/06 */

/* Computes the momentum space quark propagators for 16 point sources,
   needed for the computation of the renormalization constant of bilinear
   quark operators for the given (gauge fixed) gauge field configuration */

/* This is a based on a "stripped down" version of quark_prop.c without
   reading or writting of real space propagators */

#include "gluon_prop_includes.h"
#include <errno.h>

#undef MAX_BUF_LENGTH
#define MAX_BUF_LENGTH 4096

void get_qprop_data(int ix, int iy, int iz, int it, Real *r1,
                    Real *r2, Real *r3, int *currentnode);
void write_mom_ksprop(FILE *fp, int xi, int j, field_offset src);

int quark_renorm( void ) {
register int i, dir;
register site *s;
int j, cgn;
Real mass_x2, finalrsq;
Real pix, piy, piz, pit;
Real sin_pmu, q_mu, prop_a, prop_b, z_fac, m_func, ftmp = 0;
Real r1, r2, r3;
int pmu, px, py, pz, pt;
int pxn, pyn, pzn, ptn;
int currentnode;
int j1, jm2, k, dirs[4];
msg_tag *mtag[2];
int j_mass;
su3_vector **psim = NULL;
int xi, j2, j3, j4, parity;
  imp_ferm_links_t **fn;

int multiflag;
FILE *fp_mom_ks[MAX_NUM_MASS];	/* for writing mom propagator files */
char filename[50];
 int prec = PRECISION;   /* Make internal precision for CG the same as
			    the prevailing precision */

    pix = 2.*PI / (Real)nx;
    piy = 2.*PI / (Real)ny;
    piz = 2.*PI / (Real)nz;
    pit = 2.*PI / (Real)nt;

    cgn = 0;

    if( num_mass == 1){
	multiflag = 0;
    }
    else{
	multiflag = 1;
	psim = (su3_vector **)malloc(num_mass*sizeof(su3_vector *));
	for(j=0; j<num_mass; j++){
	    psim[j] = (su3_vector *)malloc(sites_on_node*sizeof(su3_vector));
	}
    }

    /* Open the momentum propagator files */
    if(this_node == 0){
	for(j=0; j<num_mass; j++){
	    sprintf(filename,"mom_pt_prop.m_%d",j);
	    fp_mom_ks[j] = fopen(filename, "ab");
	    if(fp_mom_ks[j] == NULL){
		printf("quark_renorm: Node %d can't open file %s, error %d\n",
		       this_node,filename,errno);fflush(stdout);
		terminate(1);
	    }
	}
    }

    rephase( ON );	/* Turn staggered phases on */

    /* Create fat and long links */
    restore_fermion_links_from_site(fn_links, PRECISION);
    fn = get_fm_links(fn_links);
    /* If we want HISQ support, we need the Naik term epsilon index
        here for now we use fn[0] only, which has no epsilon
        correction. */

    /* Loop over the 16 source points */
    for(xi=0; xi<16; xi++){

	/* Initialize color trace of the propagator */
	FORALLSITES(i,s){
	    for(j=0; j<num_mass; j++){
		s->trace_prop[j].real = 0.0;
		s->trace_prop[j].imag = 0.0;
	    }
	}

	j1 = xi%2;
	k = xi/2;
	j2 = k%2;
	k /= 2;
	j3 = k%2;
	k /= 2;
	j4 = k%2;
	parity = (j1+j2+j3+j4)%2;
/*	dirs[XUP] = j1;
	dirs[YUP] = j2;
	dirs[ZUP] = j3;
	dirs[TUP] = j4; */

	/* Loop over colors of source vector */
	for(j=0; j<3; j++){

	    /* initialize the source in phi */
	    FORALLSITES(i,s){
		clearvec( &(s->phi));
	    }

	    /* Point source at site xi in the hypercube at origin */
	    if( node_number(j1,j2,j3,j4) == this_node ){
		i=node_index(j1,j2,j3,j4);
		lattice[i].phi.c[j].real = -1.0;
	    }

	    if( multiflag == 0){
		for(j_mass=0; j_mass<num_mass; j_mass++){
		    FORALLSITES(i,s){
			clearvec( &(s->xxx1));
		    }

		    if(parity == 0){
			/* do a C.G. (source in phi, result in xxx1) */
		      cgn += ks_congrad( F_OFFSET(phi), F_OFFSET(xxx1),
					 mass[j_mass], niter, nrestart, 
					 rsqprop,  PRECISION, 
					 EVEN, &finalrsq, fn[0]);
		      /* Multiply by -Madjoint */
		      dslash_site( F_OFFSET(xxx1), F_OFFSET(ttt), ODD,
				   fn[0]);
		      mass_x2 = 2.*mass[j_mass];
		      FOREVENSITES(i,s){
			scalar_mult_su3_vector( &(s->xxx1), -mass_x2,
						&(s->ttt));
		      }
		    }
		    else{
			/* do a C.G. (source in phi, result in xxx1) */
			cgn += ks_congrad( F_OFFSET(phi), F_OFFSET(xxx1),
					   mass[j_mass], niter, nrestart, 
					   rsqprop, PRECISION, ODD, &finalrsq,
					   fn[0]);
			/* Multiply by -Madjoint */
			dslash_site( F_OFFSET(xxx1), F_OFFSET(ttt), EVEN,
				     fn[0]);
			mass_x2 = 2.*mass[j_mass];
			FORODDSITES(i,s){
			    scalar_mult_su3_vector( &(s->xxx1), -mass_x2,
						    &(s->ttt));
			}
		    }

		    /* Get contribution to color trace of the propagator */
		    FORALLSITES(i,s){
			CSUM( s->trace_prop[j_mass], s->ttt.c[j]);
		    }

		    /* Fourier transform and save */
		    restrict_fourier_site(F_OFFSET(ttt), sizeof(su3_vector),
					  FORWARDS);
		    write_mom_ksprop(fp_mom_ks[j_mass], xi, j, F_OFFSET(ttt));
		}
	    }
	    else{
	      /* Copy pointers for fermion links.  Should be based on Naik epsilon indices */
	      imp_ferm_links_t **fn_multi = 
		(imp_ferm_links_t **)malloc(sizeof(imp_ferm_links_t *)*num_mass);
	      for(j_mass = 0; j_mass < num_mass; j_mass++)
		/* Here we should use the Naik epsilon index appropriate to each fermion */
		fn_multi[j_mass] = fn[0];
	      if(parity == 0){
		/* do a multi-cg */
		cgn += ks_multicg_mass_site( F_OFFSET(phi), psim, mass, num_mass,
					     niter, rsqprop, prec, EVEN, &finalrsq,
					     fn_multi);
		/* Multiply by -Madjoint */
		for(j_mass=0; j_mass<num_mass; j_mass++){
		  FORALLSITES(i,s){
		    su3vec_copy( &(psim[j_mass][i]), &(s->xxx1));
		  }
		  dslash_site( F_OFFSET(xxx1), F_OFFSET(ttt), ODD,
			       fn_multi[j_mass]);
		  mass_x2 = 2.*mass[j_mass];
		  FOREVENSITES(i,s){
		    scalar_mult_su3_vector( &(s->xxx1), -mass_x2,
					    &(s->ttt));
		  }
		  
		  /* Get contribution to color trace of the propagator */
		  FORALLSITES(i,s){
		    CSUM( s->trace_prop[j_mass], s->ttt.c[j]);
		  }
		  
		  /* Fourier transform and save */
		  restrict_fourier_site(F_OFFSET(ttt), 
					sizeof(su3_vector),
					FORWARDS);
		  write_mom_ksprop(fp_mom_ks[j_mass], xi, j, F_OFFSET(ttt));
		}
	      }
	      else{
		/* do a multi-cg */
		cgn += ks_multicg_mass_site( F_OFFSET(phi), psim, mass, num_mass,
					     niter, rsqprop, prec, ODD, &finalrsq,
					     fn_multi);
		/* Multiply by -Madjoint */
		for(j_mass=0; j_mass<num_mass; j_mass++){
		  FORALLSITES(i,s){
		    su3vec_copy( &(psim[j_mass][i]), &(s->xxx1));
		  }
		  dslash_site( F_OFFSET(xxx1), F_OFFSET(ttt), EVEN,
			       fn_multi[j_mass]);
		  mass_x2 = 2.*mass[j_mass];
		  FORODDSITES(i,s){
		    scalar_mult_su3_vector( &(s->xxx1), -mass_x2,
					    &(s->ttt));
		  }
		  
		  /* Get contribution to color trace of the propagator */
		  FORALLSITES(i,s){
		    CSUM( s->trace_prop[j_mass], s->ttt.c[j]);
		  }
		  
		  /* Fourier transform and save */
		  restrict_fourier_site(F_OFFSET(ttt),
					sizeof(su3_vector),
					FORWARDS);
		  write_mom_ksprop(fp_mom_ks[j_mass], xi, j, F_OFFSET(ttt));
		}
	      } /* source parity */
	      free(fn_multi);
	    } /* multiflag */
	} /* colors */

	if ( xi > 0 ) continue;

/*	if ( xi > 0 ){
	    g_sync();
	    mtag[0] = start_general_gather_site( F_OFFSET(trace_prop[0]),
			num_mass*sizeof(complex), dirs, EVENANDODD, gen_pt[0]);
	    wait_general_gather(mtag[0]);
	    for(j_mass=0; j_mass<num_mass; j_mass++){
		FORALLSITES(i,s){
		    s->ttt.c[0] = ((complex *)gen_pt[0][i])[j_mass];
		}
		FORALLSITES(i,s){
		    s->trace_prop[j_mass] = s->ttt.c[0];
		}
	    }
	    cleanup_general_gather(mtag[0]);
	} */

	g_sync();
	/* Now Fourier transform */
	for(j_mass=0; j_mass<num_mass; j_mass+=3){
	    if( (j_mass+2) < num_mass ){
		j = 3;
	    }
	    else{
		j = num_mass - j_mass;
	    }
	    restrict_fourier_site(F_OFFSET(trace_prop[j_mass]),
				  j*sizeof(complex), FORWARDS);
	}

	for(j_mass=0; j_mass<num_mass; j_mass++){
	  /* Now multiply accumulate sum_mu q_mu Im(prop) in qprop[0] */
	  /* Also accumulate sum_mu q^2_mu, where q_mu = sin(p_mu) for 1-link, */
	  /* ie. standard, an q_mu = sin(p_mu)*(1+sin^2(p_mu)/6) for improved */
	  /* fermions with Naik term */
	  FORALLSITES(i,s){
	    for(dir=XUP; dir<=TUP; dir++)
	    {
	      switch(dir){
		case XUP:
		  pmu = s->x;
#ifdef NAIK
		  sin_pmu = sin((double)(pmu*pix));
		  q_mu = sin_pmu * (1.0 + sin_pmu*sin_pmu/6.0);
#else
		  q_mu = sin((double)(pmu*pix));
#endif
		  s->tempfloat = q_mu * q_mu;
		  ftmp = q_mu;
		  break;

		case YUP:
		  pmu = s->y;
#ifdef NAIK
		  sin_pmu = sin((double)(pmu*piy));
		  q_mu = sin_pmu * (1.0 + sin_pmu*sin_pmu/6.0);
#else
		  q_mu = sin((double)(pmu*piy));
#endif
		  s->tempfloat += q_mu * q_mu;
		  ftmp += q_mu;
		  break;

		case ZUP:
		  pmu = s->z;
#ifdef NAIK
		  sin_pmu = sin((double)(pmu*piz));
		  q_mu = sin_pmu * (1.0 + sin_pmu*sin_pmu/6.0);
#else
		  q_mu = sin((double)(pmu*piz));
#endif
		  s->tempfloat += q_mu * q_mu;
		  ftmp += q_mu;
		  break;

		case TUP:
		  pmu = s->t;
#ifdef NAIK
		  sin_pmu = sin((double)(pmu*pit));
		  q_mu = sin_pmu * (1.0 + sin_pmu*sin_pmu/6.0);
#else
		  q_mu = sin((double)(pmu*pit));
#endif
		  s->tempfloat += q_mu * q_mu;
		  ftmp += q_mu;
		  break;

		default: printf("BOTCH: bad direction\n"); exit(1);
	      }

	    }

	    s->qprop[0] = -ftmp * s->trace_prop[j_mass].imag;
	    s->qprop[1] = s->trace_prop[j_mass].real;
	  }

	  /* Now add the contributions from the 16 parts of the Brillouin zone */
	  for(j=0; j<16; j++){
	    jm2 = j%2;

	    if (j > 0) wait_general_gather(mtag[jm2]);

	    if (j < 15){
	      j1 = j + 1;
	      k = j1%2;
	      dirs[XUP] = k * (nx/2);
	      j1 /= 2;
	      dirs[YUP] = (j1%2) * (ny/2);
	      j1 /= 2;
	      dirs[ZUP] = (j1%2) * (nz/2);
	      j1 /= 2;
	      dirs[TUP] = (j1%2) * (nt/2);
	      mtag[k] = start_general_gather_site( F_OFFSET(qprop[0]),
				2*sizeof(Real), dirs, EVENANDODD, gen_pt[k]);
	    }

	    if (j == 0){
	      FORALLSITES(i,s){
		s->sum_qprop[0] = s->qprop[0];
		s->sum_qprop[1] = s->qprop[1];
	      }
	    }
	    else{
	      FORALLSITES(i,s){
		s->sum_qprop[0] += ((Real *)gen_pt[jm2][i])[0];
		s->sum_qprop[1] += ((Real *)gen_pt[jm2][i])[1];
	      }
	    }
	  }

	  cleanup_general_gather(mtag[0]);
	  cleanup_general_gather(mtag[1]);

	  /* Finally node 0 averages over permutations of spatial momentum
	     components and writes the quark propagator */
	  if(this_node==0){
	    printf("START_QUARK_PROP\n");
	    printf("mass_no %d mass %f\n", j_mass, (double)mass[j_mass]);
	  }
	  g_sync();
	  currentnode=0;

	  for(pt=0;pt<nt/4;pt++)for(px=0;px<=nx/4;px++)if(pt!=0 || px!=0){
	    for(py=0;py<=px;py++)for(pz=0;pz<=py;pz++){

	      if(pt>0){
		ptn = nt/2 - pt;
	      }
	      else
		ptn = 0;
	      if(px>0){
		pxn = nx/2 - px;
	      }
	      else
		pxn = 0;
	      if(py>0){
		pyn = ny/2 - py;
	      }
	      else
		pyn = 0;
	      if(pz>0){
		pzn = nz/2 - pz;
	      }
	      else
		pzn = 0;

	      get_qprop_data(px, py, pz, pt, &prop_a, &prop_b, &sin_pmu,
			     &currentnode);
	      get_qprop_data(pxn, py, pz, pt, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(px, pyn, pz, pt, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(px, py, pzn, pt, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(px, py, pz, ptn, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(pxn, pyn, pz, pt, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(pxn, py, pzn, pt, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(pxn, py, pz, ptn, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(px, pyn, pzn, pt, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(px, pyn, pz, ptn, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(px, py, pzn, ptn, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(pxn, pyn, pzn, pt, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(pxn, pyn, pz, ptn, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(pxn, py, pzn, ptn, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(px, pyn, pzn, ptn, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(pxn, pyn, pzn, ptn, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;

	      /* 2nd permutation */
	      get_qprop_data(py, px, pz, pt, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(pyn, px, pz, pt, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(py, pxn, pz, pt, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(py, px, pzn, pt, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(py, px, pz, ptn, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(pyn, pxn, pz, pt, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(pyn, px, pzn, pt, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(pyn, px, pz, ptn, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(py, pxn, pzn, pt, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(py, pxn, pz, ptn, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(py, px, pzn, ptn, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(pyn, pxn, pzn, pt, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(pyn, pxn, pz, ptn, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(pyn, px, pzn, ptn, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(py, pxn, pzn, ptn, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(pyn, pxn, pzn, ptn, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;

	      /* 3rd permutation */
	      get_qprop_data(pz, py, px, pt, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(pzn, py, px, pt, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(pz, pyn, px, pt, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(pz, py, pxn, pt, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(pz, py, px, ptn, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(pzn, pyn, px, pt, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(pzn, py, pxn, pt, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(pzn, py, px, ptn, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(pz, pyn, pxn, pt, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(pz, pyn, px, ptn, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(pz, py, pxn, ptn, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(pzn, pyn, pxn, pt, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(pzn, pyn, px, ptn, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(pzn, py, pxn, ptn, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(pz, pyn, pxn, ptn, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(pzn, pyn, pxn, ptn, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;

	      /* 4th permutation */
	      get_qprop_data(px, pz, py, pt, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(pxn, pz, py, pt, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(px, pzn, py, pt, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(px, pz, pyn, pt, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(px, pz, py, ptn, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(pxn, pzn, py, pt, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(pxn, pz, pyn, pt, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(pxn, pz, py, ptn, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(px, pzn, pyn, pt, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(px, pzn, py, ptn, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(px, pz, pyn, ptn, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(pxn, pzn, pyn, pt, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(pxn, pzn, py, ptn, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(pxn, pz, pyn, ptn, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(px, pzn, pyn, ptn, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(pxn, pzn, pyn, ptn, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;

	      /* 5th permutation */
	      get_qprop_data(py, pz, px, pt, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(pyn, pz, px, pt, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(py, pzn, px, pt, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(py, pz, pxn, pt, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(py, pz, px, ptn, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(pyn, pzn, px, pt, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(pyn, pz, pxn, pt, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(pyn, pz, px, ptn, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(py, pzn, pxn, pt, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(py, pzn, px, ptn, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(py, pz, pxn, ptn, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(pyn, pzn, pxn, pt, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(pyn, pzn, px, ptn, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(pyn, pz, pxn, ptn, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(py, pzn, pxn, ptn, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(pyn, pzn, pxn, ptn, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;

	      /* 6th permutation */
	      get_qprop_data(pz, px, py, pt, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(pzn, px, py, pt, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(pz, pxn, py, pt, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(pz, px, pyn, pt, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(pz, px, py, ptn, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(pzn, pxn, py, pt, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(pzn, px, pyn, pt, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(pzn, px, py, ptn, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(pz, pxn, pyn, pt, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(pz, pxn, py, ptn, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(pz, px, pyn, ptn, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(pzn, pxn, pyn, pt, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(pzn, pxn, py, ptn, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(pzn, px, pyn, ptn, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(pz, pxn, pyn, ptn, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;
	      get_qprop_data(pzn, pxn, pyn, ptn, &r1, &r2, &r3, &currentnode);
	      prop_a += r1;
	      prop_b += r2;

	      if(this_node==0){
		/* factor 16*6*16*3*0.5 (perm,2^4,Nc,MILC convention for dslash) */
		prop_a /= (2304.0);
		prop_b /= (2304.0);
		if(sin_pmu > 1.0e-6){
		  z_fac = prop_a; 
		  prop_a /= sin_pmu;
		  m_func = prop_b / prop_a;
		  z_fac *= prop_a;
		  z_fac += prop_b * prop_b;
		  z_fac /= prop_a;
		}
		else{
		  m_func = -mass[j_mass];
		  z_fac = -1.0;
		}
		printf(" %e %e  %e %e\n", (double)prop_a, (double)prop_b,
		       (double)z_fac, (double)m_func);
	      }
	    }
	  }
	  if(this_node==0) printf("END_QUARK_PROP\n");
	}	/* j_mass */
    }	/* xi */

    rephase( OFF );	/* Turn staggered phases off */

    if( num_mass > 1){
	for(j=0; j<num_mass; j++) free(psim[j]);
	free(psim);
    }

    /* Close the momentum propagator files */
    if(this_node == 0){
	for(j=0; j<num_mass; j++){
	    fclose(fp_mom_ks[j]);
	}
    }

    return(cgn);
} /* quark_prop */

void get_qprop_data(int ix, int iy, int iz, int it, Real *r1,
                    Real *r2, Real *r3, int *currentnode){

int i, newnode;
struct {
  Real f1, f2, f3;
} msg;

  newnode = node_number(ix,iy,iz,it);
  if(newnode != *currentnode){	/* switch to another node */
    /* tell newnode it's OK to send */
    if( this_node==0 && newnode!=0 )send_field((char *)&i,4,newnode);
    if( this_node==newnode && newnode!=0 )get_field((char *)&i,4,0);
    *currentnode = newnode;
  }

  if(this_node==0){
    if(*currentnode==0){
      i = node_index(ix,iy,iz,it);
      *r1 = lattice[i].sum_qprop[0];
      *r2 = lattice[i].sum_qprop[1];
      *r3 = lattice[i].tempfloat;
    }
    else{
      get_field((char *)&msg,3*sizeof(Real),*currentnode);
      *r1 = msg.f1;
      *r2 = msg.f2;
      *r3 = msg.f3;
    }
  }
  else{	/* for nodes other than 0 */
    if(this_node==*currentnode){
      i = node_index(ix,iy,iz,it);
      msg.f1 = lattice[i].sum_qprop[0];
      msg.f2 = lattice[i].sum_qprop[1];
      msg.f3 = lattice[i].tempfloat;
      send_field((char *)&msg,3*sizeof(Real),0);
    }
  }
}

void write_mom_ksprop(FILE *fp, int xi, int color, field_offset src){

int i,j,x,y,z,t;
int currentnode,newnode;
int buf_length;
su3_vector *pbuf = NULL;

struct {
  su3_vector ksv;
} msg;

  if(this_node==0){
    pbuf = (su3_vector *)malloc(MAX_BUF_LENGTH*sizeof(su3_vector));
    if(pbuf == NULL){
      printf("write_mom_ksprop: Node 0 can't malloc pbuf\n");
      fflush(stdout); terminate(1);
    }
    if( (int)fwrite(&xi,sizeof(int),1,fp) != 1){
      printf("write_mom_ksprop: Node %d prop write error %d\n",this_node,errno);
      fflush(stdout);
      terminate(1);
    }
    if( (int)fwrite(&color,sizeof(int),1,fp) != 1){
      printf("write_mom_ksprop: Node %d prop write error %d\n",this_node,errno);
      fflush(stdout);
      terminate(1);
    }
  }

  /* Buffered algorithm for writing fields in serial order */

  g_sync();
  currentnode=0;
 
  buf_length = 0;
 
  for(j=0,t=0;t<nt;t++)for(z=0;z<nz;z++)for(y=0;y<ny;y++)for(x=0;x<nx;x++,j++){
    newnode=node_number(x,y,z,t);
    if(newnode != currentnode){       /* switch to another node */
      /* Send a few bytes of garbage to tell newnode it's OK to send */
      if( this_node==0 && newnode!=0 )
	send_field((char *)&msg,sizeof(msg),newnode);
      if( this_node==newnode && newnode!=0 )
	get_field((char *)&msg,sizeof(msg),0);
      currentnode=newnode;
    }
                                                                                
    if(this_node==0){
      if(currentnode==0){
	i=node_index(x,y,z,t);
	msg.ksv = *((su3_vector *)F_PT(&lattice[i],src));
      }
      else{
	get_field((char *)&msg, sizeof(msg), currentnode);
      }
                                                                                
      pbuf[buf_length] = msg.ksv;
      buf_length++;
           
      if( (buf_length == MAX_BUF_LENGTH) || (j == volume-1)){
	/* write out buffer */
	if( (int)fwrite(pbuf,sizeof(su3_vector),buf_length,fp) != buf_length){
	  printf("write_mom_ksprop: Node %d prop write error %d\n",
		 this_node,errno);
       	  fflush(stdout);
       	  terminate(1);
       	}
       	buf_length = 0;           /* start again after write */
      }
    }
    else{ /* for nodes other than 0 */
      if(this_node==currentnode){
	i=node_index(x,y,z,t);
	/* Copy data into send buffer and send to node 0 with padding */
	msg.ksv = *((su3_vector *)F_PT(&lattice[i],src));
	send_field((char *)&msg, sizeof(msg),0);
      }
    }
  } /*close x,y,z,t loops */
   
  g_sync();
  if(this_node==0){
    free(pbuf);
  }
}
