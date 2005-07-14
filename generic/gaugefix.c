/************************** gaugefix.c *******************************/
/* Fix Coulomb or Lorentz gauge by doing successive SU(2) gauge hits */
/* Uses double precision global sums */
/* MIMD version 7 */
/* C. DeTar 10-22-90 */
/* T. DeGrand 1993 */
/* U.M. Heller 8-31-95 */
/* C. DeTar 10-11-97 converted to generic */

/* Prototype...

   void gaugefix(int gauge_dir,Real relax_boost,int max_gauge_iter,
	      Real gauge_fix_tol, field_offset diffmat, field_offset sumvec,
	      int nvector, field_offset vector_offset[], int vector_parity[],
	      int nantiherm, field_offset antiherm_offset[], 
	      int antiherm_parity[] )
   -------------------------------------------------------------------

   NOTE: For staggered fermion applications, it is necessary to remove
   the KS phases from the gauge links before calling this procedure.
   See "rephase" in setup.c.

   -------------------------------------------------------------------
   EXAMPLE:  Fixing only the link matrices to Coulomb gauge with scratch
     space in mp (su3_matrix) and chi (su3_vector):

   gaugefix(TUP,(Real)1.5,500,(Real)1.0e-7,
		     F_OFFSET(mp),F_OFFSET(chi),0,NULL,NULL,0,NULL,NULL);

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
   gaugefix(YUP,(Real)1.8,500,(Real)2.0e-6,
       F_OFFSET(tempmat1),F_OFFSET(tempvec[0]),
       nvector,vector_offset,vector_parity,
       nantiherm,antiherm_offset,antiherm_parity);
   rephase( ON );

   -------------------------------------------------------------------

   gauge_dir     specifies the direction of the "time"-like hyperplane
                 for the purposes of defining Coulomb or Lorentz gauge 
      TUP    for evaluating propagators in the time-like direction
      ZUP    for screening lengths.
      8      for Lorentz gauge
   relax_boost	   Overrelaxation parameter 
   max_gauge_iter  Maximum number of iterations 
   gauge_fix_tol   Stop if change is less than this 
   diffmat         Scratch space for an su3 matrix
   sumvec          Scratch space for an su3 vector
   NOTE: if diffmat or sumvec are negative, gaugefix mallocs its own
   scratch space. */

#include "generic_includes.h"

/* Generic definitions - could be useful elsewhere */

/*    CDIF(a,b)         a -= b						      */
								/*  a -= b    */
#define CDIF(a,b) { (a).real -= (b).real; (a).imag -= (b).imag; }

/* Scratch space */

su3_matrix *diffmatp;                       /* malloced diffmat pointer */
su3_vector *sumvecp;                        /* malloced sumvec pointer */
field_offset diffmat_offset,sumvec_offset;  /* field offsets */

void mult_su2_mat_vec_elem_n(su2_matrix *u,complex *x0,complex *x1)
{
  /* Multiplies the complex column spinor (x0, x1) by the SU(2) matrix u */
  /* and puts the result in (x0,x1).  */
  /* Thus x <- u * x          */
  /* C. DeTar 3 Oct 1990 */
  
  complex z0, z1, t0, t1;

  t0 = *x0; t1 = *x1;

  CMUL(u->e[0][0], t0, z0);
  CMUL(u->e[0][1], t1, z1);
  CADD(z0, z1, *x0);
  CMUL(u->e[1][0], t0, z0);
  CMUL(u->e[1][1], t1, z1);
  CADD(z0, z1, *x1);

} /* mult_su2_mat_vec_elem_n */

void mult_su2_mat_vec_elem_a(su2_matrix *u,complex *x0,complex *x1)
{
  /* Multiplies the complex row spinor (x0, x1) by the adjoint of the */
  /* SU(2) matrix u and puts the result in (x0,x1).  */
  /* Thus x <-  x * u-adj       */
  /* C. DeTar 3 Oct 1990 */
  
  complex z0, z1, t0, t1;

  t0 = *x0; t1 = *x1;

  CMUL_J(t0, u->e[0][0], z0);
  CMUL_J(t1, u->e[0][1], z1);
  CADD(z0, z1, *x0);
  CMUL_J(t0, u->e[1][0], z0);
  CMUL_J(t1, u->e[1][1], z1);
  CADD(z0, z1, *x1);

} /* mult_su2_mat_vec_elem_a */

void dumpsu2(su2_matrix *u)
{
  int i,j;
  for(i=0;i<2;i++){
    for(j=0;j<2;j++)printf("(%.2e,%.2e)\t",
			  (double)u->e[i][j].real,(double)u->e[i][j].imag);
    printf("\n");
  }
  printf("\n");
}

void left_su2_hit_n(su2_matrix *u,int p,int q,su3_matrix *link)
{
  /* link <- u * link */
  /* The 0 row of the SU(2) matrix u matches row p of the SU(3) matrix */
  /* The 1 row of the SU(2) matrix u matches row q of the SU(3) matrix */
  /* C. DeTar 18 Oct 1990 */

  register int m;

  for (m = 0; m < 3; m++)
    mult_su2_mat_vec_elem_n(u, &(link->e[p][m]), &(link->e[q][m]));

} /* left_su2_hit_n */

void right_su2_hit_a(su2_matrix *u,int p,int q,su3_matrix *link)
{
  /* link <-  link * u adj */
  /* The 0 column of u-adjoint matches column p of the SU(3) matrix */
  /* The 1 column of u-adjoint matches column q of the SU(3) matrix */
  /* C. DeTar 18 Oct 1990 */

  register int m;

  for (m = 0; m < 3; m++)
    mult_su2_mat_vec_elem_a(u, &(link->e[m][p]), &(link->e[m][q]));

} /*right_su2_hit_a */

void accum_gauge_hit(int gauge_dir,int parity)
{

/* Accumulates sums and differences of link matrices for determining optimum */
/* hit for gauge fixing */
/* Differences are kept in diffmat and the diagonal elements of the sums */
/* in sumvec  */

  register int j;
  register su3_matrix *m1,*m2;
  register int dir,i;
  register site *s;

  /* Clear sumvec and diffmat */

  FORSOMEPARITY(i,s,parity)
    {
      if(diffmat_offset >= 0)
	clear_su3mat((su3_matrix *)F_PT(s,diffmat_offset));
      else
	clear_su3mat(&diffmatp[i]);
      if(sumvec_offset >= 0)
	clearvec((su3_vector *)F_PT(s,sumvec_offset));
      else
	clearvec(&sumvecp[i]);
    }
  
  /* Subtract upward link contributions */

  FORSOMEPARITY(i,s,parity)
    {
      FORALLUPDIRBUT(gauge_dir,dir)
	{
	  /* Upward link matrix */
	  m1 = &(s->link[dir]);
	  if(diffmat_offset >= 0)
	    sub_su3_matrix((su3_matrix *)F_PT(s,diffmat_offset), 
			   m1, (su3_matrix *)F_PT(s,diffmat_offset)); 
	  else
	    sub_su3_matrix( &diffmatp[i], m1, &diffmatp[i]); 

	  if(sumvec_offset >= 0)
	    {
	      for(j=0;j<3;j++)CSUM( ((su3_vector *)F_PT(s,sumvec_offset))->c[j],
				    m1->e[j][j]);
	    }
	  else
	    {
	      for(j=0;j<3;j++)CSUM( sumvecp[i].c[j],m1->e[j][j]);
	    }
	}
    }

  /* Add downward link contributions */

  FORSOMEPARITY(i,s,parity)
    {
      FORALLUPDIRBUT(gauge_dir,dir)
	{
	  /* Downward link matrix */
	  m2 = (su3_matrix *)gen_pt[dir][i];

	  if(diffmat_offset >= 0)
	    add_su3_matrix((su3_matrix *)F_PT(s,diffmat_offset), m2, 
			   (su3_matrix *)F_PT(s,diffmat_offset));
	  else
	    add_su3_matrix( &diffmatp[i], m2, &diffmatp[i]);
	    
	  if(sumvec_offset >= 0)
	    {
	      for(j=0;j<3;j++)CSUM( ((su3_vector *)F_PT(s,sumvec_offset))->c[j],
				    m2->e[j][j]);
	    }
	  else
	    {
	      for(j=0;j<3;j++)CSUM( sumvecp[i].c[j], m2->e[j][j]);
	    }

	  /* Add diagonal elements to sumvec  */
	}
    }
} /* accum_gauge_hit */


void do_hit(int gauge_dir, int parity, int p, int q, Real relax_boost,
	      int nvector, field_offset vector_offset[], int vector_parity[],
	      int nantiherm, field_offset antiherm_offset[], 
	      int antiherm_parity[] )
{
  /* Do optimum SU(2) gauge hit for p, q subspace */

  Real a0,a1,a2,a3,asq,a0sq,x,r,xdr;
  register int dir,i,j;
  register site *s;
  su2_matrix u;
  su3_matrix htemp;

  /* Accumulate sums for determining optimum gauge hit */

  accum_gauge_hit(gauge_dir,parity);

  FORSOMEPARITY(i,s,parity)
    {
      /* The SU(2) hit matrix is represented as a0 + i * Sum j (sigma j * aj)*/
      /* The locally optimum unnormalized components a0, aj are determined */
      /* from the current link in direction dir and the link downlink */
      /* in the same direction on the neighbor in the direction opposite dir */
      /* The expression is */
      /* a0 = Sum dir Tr Re 1       * (downlink dir + link dir) */
      /* aj = Sum dir Tr Im sigma j * (downlink dir - link dir)  j = 1,2, 3 */
      /*   where 1, sigma j are unit and Pauli matrices on the p,q subspace */
 /*     
      a0 =  s->sumvec.c[p].real + s->sumvec.c[q].real;
      a1 =  s->diffmat.e[q][p].imag + s->diffmat.e[p][q].imag;
      a2 = -s->diffmat.e[q][p].real + s->diffmat.e[p][q].real;
      a3 =  s->diffmat.e[p][p].imag - s->diffmat.e[q][q].imag;
*/  
      if(sumvec_offset >= 0)
	a0 =     ((su3_vector *)F_PT(s,sumvec_offset))->c[p].real +
	  ((su3_vector *)F_PT(s,sumvec_offset))->c[q].real;
      else
	a0 =     sumvecp[i].c[p].real +  sumvecp[i].c[q].real;
	
      if(diffmat_offset >= 0)
	{
	  a1 =     ((su3_matrix *)F_PT(s,diffmat_offset))->e[q][p].imag +
	    ((su3_matrix *)F_PT(s,diffmat_offset))->e[p][q].imag;
	  a2 =    -((su3_matrix *)F_PT(s,diffmat_offset))->e[q][p].real +
	    ((su3_matrix *)F_PT(s,diffmat_offset))->e[p][q].real;
	  a3 =     ((su3_matrix *)F_PT(s,diffmat_offset))->e[p][p].imag -
	    ((su3_matrix *)F_PT(s,diffmat_offset))->e[q][q].imag;
	}
      else
	{
	  a1 =     diffmatp[i].e[q][p].imag + diffmatp[i].e[p][q].imag;
	  a2 =    -diffmatp[i].e[q][p].real + diffmatp[i].e[p][q].real;
	  a3 =     diffmatp[i].e[p][p].imag - diffmatp[i].e[q][q].imag;
	}

      /* Over-relaxation boost */

      /* This algorithm is designed to give little change for large |a| */
      /* and to scale up the gauge transformation by a factor of relax_boost*/
      /* for small |a| */

      asq = a1*a1 + a2*a2 + a3*a3;
      a0sq = a0*a0;
      x = (relax_boost*a0sq + asq)/(a0sq + asq);
      r = sqrt((double)(a0sq + x*x*asq));
      xdr = x/r;
      /* Normalize and boost */
      a0 = a0/r; a1 = a1*xdr; a2 = a2*xdr; a3 = a3*xdr;

      /* Elements of SU(2) matrix */

      u.e[0][0] = cmplx( a0, a3);
      u.e[0][1] = cmplx( a2, a1);
      u.e[1][0] = cmplx(-a2, a1);
      u.e[1][1] = cmplx( a0,-a3);

      
      /* Do SU(2) hit on all upward links */

      FORALLUPDIR(dir)
	left_su2_hit_n(&u,p,q,&(s->link[dir]));
      
      /* Do SU(2) hit on all downward links */
      
      FORALLUPDIR(dir)
	right_su2_hit_a(&u,p,q,(su3_matrix *)gen_pt[dir][i]);
      
      /* Transform vectors and gauge momentum if requested */

      for(j = 0; j < nvector; j++)
	  
	/* Do SU(2) hit on specified su3 vector for specified parity */
	  
	/* vector <- u * vector */
	if(vector_parity[j] == EVENANDODD || vector_parity[j] == parity)
	  mult_su2_mat_vec_elem_n(&u, 
		 &((su3_vector *)F_PT(s,vector_offset[j]))->c[p],
		 &((su3_vector *)F_PT(s,vector_offset[j]))->c[q]);

      /* Transform antihermitian matrices if requested */

      for(j = 0; j < nantiherm; j++)
	/* antiherm <- u * antiherm * u^dagger */
	if(antiherm_parity[j] == EVENANDODD || antiherm_parity[j] == parity)
	  {
	    uncompress_anti_hermitian( 
		     (anti_hermitmat *)F_PT(s,antiherm_offset[j]), &htemp);
	    /* If the next 2 steps prove too time consuming, */
	    /* they can be simplified algebraically, and sped up by ~2 */
	    left_su2_hit_n(&u,p,q,&htemp);
	    right_su2_hit_a(&u,p,q,&htemp);
	    make_anti_hermitian( &htemp, 
		     (anti_hermitmat *)F_PT(s,antiherm_offset[j]));
	  }
    }
  
  /* Exit with modified downward links left in communications buffer */
} /* do_hit */

double get_gauge_fix_action(int gauge_dir,int parity)
{
  /* Adds up the gauge fixing action for sites of given parity */
  /* Returns average over these sites */
  /* The average is normalized to a maximum of 1 when all */
  /* links are unit matrices */

  register int dir,i,ndir;
  register site *s;
  register su3_matrix *m1, *m2;
  double gauge_fix_action;
  complex trace;

  gauge_fix_action = 0.0;
  
  FORSOMEPARITY(i,s,parity)
    {
      FORALLUPDIRBUT(gauge_dir,dir)
	{
	  m1 = &(s->link[dir]);
	  m2 = (su3_matrix *)gen_pt[dir][i];

	  trace = trace_su3(m1);
	  gauge_fix_action += (double)trace.real;

	  trace = trace_su3(m2); 
 	  gauge_fix_action += (double)trace.real;
	}
    }

  /* Count number of terms to average */
  ndir = 0; FORALLUPDIRBUT(gauge_dir,dir)ndir++;
  
  /* Sum over all sites of this parity */
  g_doublesum( &gauge_fix_action);
  
  /* Average is normalized to max of 1/2 on sites of one parity */
  return(gauge_fix_action /((double)(6*ndir*nx*ny*nz*nt)));
} /* get_gauge_fix_action */

void gaugefixstep(int gauge_dir,double *av_gauge_fix_action,Real relax_boost,
	      int nvector, field_offset vector_offset[], int vector_parity[],
	      int nantiherm, field_offset antiherm_offset[], 
	      int antiherm_parity[] )
{
  /* Carry out one iteration in the gauge-fixing process */

  int parity;
  msg_tag *mtag[8];
  Real gauge_fix_action;
  register int dir,i;
  register site *s;

  /* Alternate parity to prevent interactions during gauge transformation */
  *av_gauge_fix_action = 0.;
  g_sync();
  fflush(stdout);
  
  for(parity = ODD; parity <= EVEN; parity++)
    {
      /* Start gathers of downward links */
      
      FORALLUPDIR(dir)
	{
	  mtag[dir] = start_gather_site( F_OFFSET(link[dir]), sizeof(su3_matrix),
			   OPP_DIR(dir), parity, gen_pt[dir] );
	}
      
      /* Wait for gathers */
      
      FORALLUPDIR(dir)
         {
	  wait_gather(mtag[dir]);
	 }

      /* Total gauge fixing action for sites of this parity: Before */
      gauge_fix_action = get_gauge_fix_action(gauge_dir,parity);

      /* Do optimum gauge hit on various subspaces */

      do_hit(gauge_dir,parity,0,1, relax_boost,
		   nvector, vector_offset, vector_parity,
		   nantiherm, antiherm_offset, antiherm_parity);
      do_hit(gauge_dir,parity,1,2, relax_boost,
		   nvector, vector_offset, vector_parity,
		   nantiherm, antiherm_offset, antiherm_parity);
      do_hit(gauge_dir,parity,2,0, relax_boost,
		   nvector, vector_offset, vector_parity,
		   nantiherm, antiherm_offset, antiherm_parity);

      /* Total gauge fixing action for sites of this parity: After */
      gauge_fix_action = get_gauge_fix_action(gauge_dir,parity);
      
      *av_gauge_fix_action += gauge_fix_action;

      /* Scatter downward link matrices by gathering to sites of */
      /* opposite parity */

      FORALLUPDIR(dir)
	{
	  /* Synchronize before scattering to be sure the new modified link */
	  /* matrices are all ready to be scattered and diffmat is not */
	  /* overwritten before it is used */
	  g_sync();

	  /* First copy modified link for this dir */
	  /* from comm buffer or node to diffmat */

	  FORSOMEPARITY(i,s,parity)
	    {
	      if(diffmat_offset >= 0)
		su3mat_copy((su3_matrix *)(gen_pt[dir][i]),(su3_matrix *)F_PT(s,diffmat_offset));
	      else
		su3mat_copy((su3_matrix *)(gen_pt[dir][i]), &diffmatp[i]);
	    }
	  
	  /* Now we are finished with gen_pt[dir] */
	  cleanup_gather(mtag[dir]);
      
	  /* Synchronize to make sure the previous copy happens before the */
	  /* subsequent gather below  */
	  g_sync();

	  /* Gather diffmat onto sites of opposite parity */
	  if(diffmat_offset >= 0)
	    mtag[dir] = start_gather_site( diffmat_offset, sizeof(su3_matrix),
				      dir, OPP_PAR(parity), gen_pt[dir] );
	  else
	    mtag[dir] = start_gather_field( diffmatp, sizeof(su3_matrix),
				      dir, OPP_PAR(parity), gen_pt[dir] );

	  wait_gather(mtag[dir]);

         /* Copy modified matrices into proper location */

         FORSOMEPARITY(i,s,OPP_PAR(parity))
	      su3mat_copy((su3_matrix *)(gen_pt[dir][i]),&(s->link[dir]));

	  cleanup_gather(mtag[dir]);
	}

    }
} /* gaugefixstep */

void gaugefixscratch(field_offset diffmat, field_offset sumvec)
{
  diffmat_offset = diffmat;
  diffmatp = NULL;
  if(diffmat_offset < 0)
    {
      diffmatp = (su3_matrix *)malloc(sizeof(su3_matrix)*sites_on_node);
      if(diffmatp == NULL)
	{
	  node0_printf("gaugefix: Can't malloc diffmat\n");
	  fflush(stdout);terminate(1);
	}
    }

  sumvec_offset = sumvec;
  sumvecp = NULL;
  if(sumvec_offset < 0)
    {
      sumvecp = (su3_vector *)malloc(sizeof(su3_vector)*sites_on_node);
      if(sumvecp == NULL)
	{
	  node0_printf("gaugefix: Can't malloc sumvec\n");
	  fflush(stdout);terminate(1);
	}
    }
} /* gaugefixscratch */

void gaugefix(int gauge_dir,Real relax_boost,int max_gauge_iter,
	      Real gauge_fix_tol, field_offset diffmat, field_offset sumvec,
	      int nvector, field_offset vector_offset[], int vector_parity[],
	      int nantiherm, field_offset antiherm_offset[], 
	      int antiherm_parity[] )
{
  int gauge_iter;
  double current_av, old_av, del_av;

  /* We require at least 8 gen_pt values for gauge fixing */
  if(N_POINTERS < 8)
    {
      printf("gaugefix: N_POINTERS must be at least %d.  Fix the code.\n",
	     N_POINTERS);
      fflush(stdout); terminate(1);
    }
    

  /* Set up work space */
  gaugefixscratch(diffmat,sumvec);

  /* Do at most max_gauge_iter iterations, but stop after the second step if */
  /* the change in the avg gauge fixing action is smaller than gauge_fix_tol */

  for (gauge_iter=0; gauge_iter < max_gauge_iter; gauge_iter++)
    {
      gaugefixstep(gauge_dir,&current_av,relax_boost,
		   nvector, vector_offset, vector_parity,
		   nantiherm, antiherm_offset, antiherm_parity);

      if(gauge_iter != 0)
	{
	  del_av = current_av - old_av;
	  if (fabs(del_av) < gauge_fix_tol) break;
	}
      old_av = current_av;
    }
  /* Free workspace */
  if(diffmat_offset < 0)free(diffmatp);
  if(sumvec_offset < 0)free(sumvecp);
  
  if(this_node==0)
    printf("GFIX WITHOUT REUNITARIZATION: Ended at step %d. Av gf action %.3e, delta %.3e\n",
	   gauge_iter,(double)current_av,(double)del_av);
}
