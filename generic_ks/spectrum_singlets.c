/******** spectrum_singlets.c *************/
/* MIMD version 7*/
/* DT 2/04
   point source/sink propagators for eta_prime, sigma
   Use unit random vectors for point sources

   Any antiperiodic boundary conditions (time direction usually)
   are absorbed in the link matrices (ie rephase(ON) )

    Hopefully, the result is gauge invariant when averaged over
    sources -- NOT gauge invariant for one source.

   Kogut-Susskind fermions  -- this version for "fat plus Naik"
   or general "even plus odd" quark actions.

*/

#include "generic_ks_includes.h"

enum prop_name {
    SIGMA_DISC,
    ETA_L_DISC,		/* local (in time) */
    ETA_NL_DISC,	/* nonlocal (in time) */
    SIGMA_SIGMA_CONN,
    ETA_L_ETA_L_CONN,	/* local-local (in time) */
    ETA_L_ETA_NL_CONN,	/* local-nonlocal (in time) */
    ETA_NL_ETA_L_CONN,	/* nonlocal-local (in time) */
    ETA_NL_ETA_NL_CONN,	/* nonlocal-nonlocal (in time) */
    NPROPS
};
#define NSOURCEVECS 1
#define SMEAR 0
#define STAPLE_WEIGHT 0.25
#define APE_NHIT 100
#define APE_TOL 1e-5
#define mat_invert mat_invert_uml

int test_converge(int t_source);
void mult_sigma(su3_vector *src, su3_vector *dest, su3_vector *temp );
void mult_eta_l(su3_vector *src, su3_vector *dest, su3_vector *temp, su3_matrix *smearlink );
void mult_eta_nl(su3_vector *src, su3_vector *dest, su3_vector *temp );
static su3_matrix *smear_links();
void sym_shift_smear(int dir, su3_vector * src, su3_vector * dest,
		     su3_matrix *smearlink);
void zeta_shift_smear(int n, int *d, su3_vector * src, su3_vector * dest,
		      su3_matrix *smearlink);
void eta_shift_smear(int n, int *d, su3_vector *src, su3_vector *dest,
		     su3_matrix *smearlink);
void mult_flavor_vector_smear(int mu, su3_vector *src, su3_vector *dest, su3_matrix *smearlink );
void mult_flavor_tensor_smear(int mu, int nu, su3_vector *src, su3_vector *dest, su3_matrix *smearlink );
void mult_flavor_pseudovector_smear(int mu, su3_vector * src, su3_vector * dest, su3_matrix *smearlink );
void mult_flavor_pseudoscalar_smear(su3_vector * src, su3_vector * dest,
				    su3_matrix *smearlink);
void mult_spin_vector_smear(int mu, su3_vector *src, su3_vector *dest, su3_matrix *smearlink );
void mult_spin_tensor_smear(int mu, int nu, su3_vector *src, su3_vector *dest, su3_matrix *smearlink );
void mult_spin_pseudovector_smear(int mu, su3_vector *src, su3_vector *dest, su3_matrix *smearlink );
void mult_spin_pseudoscalar_smear(su3_vector *src, su3_vector *dest, su3_matrix *smearlink );
void rephase_smear(su3_matrix *smearlink );

su3_vector *temp, *R, *R1,*R2,*R3, *X,*X_sigma,*X_eta_l,*X_eta_nl;

#ifdef ONEMASS
#define RSRC phi
#else
#define RSRC phi1
#endif

int spectrum_singlets( Real mass, Real tol, field_offset temp_offset,
		       ferm_links_t *fn){
  /* arguments are mass, tolerance for inverter check,
   * temporary lattice su3_vector.
     return C.G. iteration number */

  int cgn;
  register int i,j;
  register site* s;
  register complex cc,czero;
  register int t_source;
  int sourcevec,color;	/* color for source */
  Real x;
  complex **props;	/* arrays of propagators */
  complex cc1,cc2,cc3,cc4,cc5;
  su3_matrix *smearlink;

  cgn=0; /* number of CG iterations */
  czero.real = czero.imag = 0.0;

  /* allocate space for quark propagators */
  temp = (su3_vector *)malloc(sites_on_node*sizeof(su3_vector));
  /* allocate space for "R" and "X" vectors.  "R..." is source multiplied
   * by meson operator, "X..." are inversion solutions */
  R = (su3_vector *)malloc(sites_on_node*sizeof(su3_vector));
  R1 = (su3_vector *)malloc(sites_on_node*sizeof(su3_vector));
  R2 = (su3_vector *)malloc(sites_on_node*sizeof(su3_vector));
  R3 = (su3_vector *)malloc(sites_on_node*sizeof(su3_vector));
  X = (su3_vector *)malloc(sites_on_node*sizeof(su3_vector));
  X_sigma = (su3_vector *)malloc(sites_on_node*sizeof(su3_vector));
  X_eta_l = (su3_vector *)malloc(sites_on_node*sizeof(su3_vector));
  X_eta_nl = (su3_vector *)malloc(sites_on_node*sizeof(su3_vector));

  /* allocate space for meson propagators */
  props = (complex **)malloc( NPROPS*sizeof(complex *) );
  props[0] = (complex *)malloc( NPROPS*nt*sizeof(complex) );
  for(i=1; i<NPROPS; i++) props[i] = props[i-1] + nt;

  /* set propagators to zero */
  for(i=0;i<NPROPS; i++)for(j=0;j<nt;j++){
    props[i][j]=czero;
  }


  /* construct APE smeared links */
  /* Note we need to get KS phases into smeared links, so use rephase_smear()
   * (later in this file) */
  smearlink = smear_links();
  rephase_smear(smearlink);

  for(sourcevec=0;sourcevec<NSOURCEVECS;sourcevec++){

    /* Construct R, R^\sigma and R^\eta as full vectors.  Take slices
     *        of them as sources for mat_invert */
    FORALLSITES(i,s){
	for(color=0;color<3;color++){
	    R[i].c[color].real  = gaussian_rand_no(&(s->site_prn));
	    R[i].c[color].imag  = gaussian_rand_no(&(s->site_prn));
	}
	x = 2.0/sqrt( magsq_su3vec( &(R[i]) ) );
	/* "2" because we invert 2m + 2 Dslash */
	scalar_mult_su3_vector( &(R[i]), x, &(R[i]) );
    }

    /* loop over "source" time slice */
    for(t_source=0; t_source<nt; t_source ++){

      /* Use mat_invert to construct X, source on one timeslice */
      FORALLSITES(i,s){
	  if(s->t==t_source){ s->RSRC = R[i]; }
	  else clearvec( &(s->RSRC) );
      }
      cgn += mat_invert( F_OFFSET(RSRC), F_OFFSET(g_rand), temp_offset, 
			 mass, PRECISION, fn  );
      FORALLSITES(i,s){ X[i]=s->g_rand; }
  
      /* disconnected parts (Eqs 24 and 29)  */
      /* Note R is full vector, only want one time slice */
      mult_sigma(  X, R1, temp );
      mult_eta_l(  X, R2, temp, smearlink );
      mult_eta_nl( X, R3, temp );
      FORALLSITES(i,s){
        if(s->t==t_source){
  	  cc = su3_dot( &(R[i]), &(R1[i]) );
	  CSUM(props[SIGMA_DISC][t_source],cc);
  	  cc = su3_dot( &(R[i]), &(R2[i]) );
	  CSUM(props[ETA_L_DISC][t_source],cc);
  	  cc = su3_dot( &(R[i]), &(R3[i]) );
	  CSUM(props[ETA_NL_DISC][t_source],cc);
        }
      }

      if(sourcevec==0){ // Only use one source vector for connected  part
        /* connected parts (Eqs 22 and 28) */
        /* Use mat_invert to construct  X^\sigma and X^\eta_l and X^\eta_nl */
        /* These are antiquark propagators, so need (-1)^X at each end */
        FORALLSITES(i,s){
  	    if(s->t==t_source) R1[i] = R[i];
  	    else clearvec( &(R1[i]) );
        }
  
        mult_sigma(  R1, R2, temp );
        FORALLSITES(i,s){
  	  s->RSRC = R2[i];
          if(s->parity==ODD)scalar_mult_su3_vector( &(s->RSRC), -1.0, &(s->RSRC) );
        }

        cgn += mat_invert( F_OFFSET(RSRC), F_OFFSET(g_rand), temp_offset, 
			   mass, PRECISION, fn  );

        FORALLSITES(i,s){
  	  X_sigma[i]=s->g_rand;
          if(s->parity==ODD)scalar_mult_su3_vector( &(X_sigma[i]), -1.0, &(X_sigma[i]) );
        }
    
        mult_eta_l(  R1, R2, temp, smearlink );
        FORALLSITES(i,s){
  	s->RSRC = R2[i];
          if(s->parity==ODD)scalar_mult_su3_vector( &(s->RSRC), -1.0, &(s->RSRC) );
        }

        cgn += mat_invert( F_OFFSET(RSRC), F_OFFSET(g_rand), temp_offset, 
			   mass, PRECISION, fn );

        FORALLSITES(i,s){
  	  X_eta_l[i]=s->g_rand;
          if(s->parity==ODD)scalar_mult_su3_vector( &(X_eta_l[i]), -1.0, &(X_eta_l[i]) );
        }
    
        // mult_eta_nl is not written yet
             FORALLSITES(i,s)clearvec( &(X_eta_nl[i]) ); //TEMPORARY
        //mult_eta_nl(  R1, R2, temp );
        //FORALLSITES(i,s){
  	//s->RSRC = R2[i];
          //if(s->parity==ODD)scalar_mult_su3_vector( &(s->RSRC), -1.0, &(s->RSRC) );
        //}
        //cgn += mat_invert( F_OFFSET(RSRC), F_OFFSET(g_rand), 
	//     temp_offset, mass, PRECISION, fn );
        //FORALLSITES(i,s){
  	//X_eta_nl[i]=s->g_rand;
          //if(s->parity==ODD)scalar_mult_su3_vector( &(X_eta_nl[i]), -1.0, &(X_eta_nl[i]) );
        //}
    
        /* store X multiplied by "sigma", xxx to store an
    	"eta_multiplied" X */
        mult_sigma( X, R1, temp ); /* no longer need R1 */
        mult_eta_l( X, R2, temp, smearlink );
        mult_eta_nl( X, R3, temp );
        FORALLSITES(i,s){
          cc1 = su3_dot( &(X_sigma[i]), &(R1[i]) );
          cc2 = su3_dot( &(X_eta_l[i]), &(R2[i]) );
          cc3 = su3_dot( &(X_eta_l[i]), &(R3[i]) );
          cc4 = su3_dot( &(X_eta_nl[i]), &(R2[i]) );
          cc5 = su3_dot( &(X_eta_nl[i]), &(R3[i]) );
          CSUM( props[SIGMA_SIGMA_CONN]  [(s->t+nt-t_source)%nt], cc1 );
          CSUM( props[ETA_L_ETA_L_CONN]  [(s->t+nt-t_source)%nt], cc2 );
          CSUM( props[ETA_L_ETA_NL_CONN] [(s->t+nt-t_source)%nt], cc3 );
          CSUM( props[ETA_NL_ETA_L_CONN] [(s->t+nt-t_source)%nt], cc4 );
          CSUM( props[ETA_NL_ETA_NL_CONN][(s->t+nt-t_source)%nt], cc5 );
        }
      } // end if(sourcevec==0)

    } /* end loop on t_source */
  } /* sourcevec */


  /* Sum propagator arrays over nodes */
  /* print out propagators */
  for(i=0;i<NPROPS;i++){
      g_veccomplexsum( props[i] , nt );
  }
//IS NORMALIZATION SAME FOR CONNECTED AND DISCONNECTED?
  for(j=0;j<nt;j++){
    CDIVREAL(props[SIGMA_DISC][j],(Real)NSOURCEVECS*nx*ny*nz,props[SIGMA_DISC][j]);
    CDIVREAL(props[ETA_L_DISC][j],(Real)NSOURCEVECS*nx*ny*nz,props[ETA_L_DISC][j]);
    CDIVREAL(props[ETA_NL_DISC][j],(Real)NSOURCEVECS*nx*ny*nz,props[ETA_NL_DISC][j]);
    CDIVREAL(props[SIGMA_SIGMA_CONN][j],(Real)nx*ny*nz,props[SIGMA_SIGMA_CONN][j]);
    CDIVREAL(props[ETA_L_ETA_L_CONN][j],(Real)nx*ny*nz,props[ETA_L_ETA_L_CONN][j]);
    CDIVREAL(props[ETA_L_ETA_NL_CONN][j],(Real)nx*ny*nz,props[ETA_L_ETA_NL_CONN][j]);
    CDIVREAL(props[ETA_NL_ETA_L_CONN][j],(Real)nx*ny*nz,props[ETA_NL_ETA_L_CONN][j]);
    CDIVREAL(props[ETA_NL_ETA_NL_CONN][j],(Real)nx*ny*nz,props[ETA_NL_ETA_NL_CONN][j]);
  }

  if(this_node==0){

    printf("STARTPROP\n");
    printf("MASSES:  %.5e   %.5e\n",mass,mass);
    printf("SOURCE: SIGMA\n");
    printf("SINKS: SIGMA_CONN\n");
    for(j=0;j<nt;j++){
      printf("%d %e %e\n",j,
      props[SIGMA_SIGMA_CONN][j].real,
      props[SIGMA_SIGMA_CONN][j].imag);
    }
    printf("ENDPROP\n");

    printf("STARTPROP\n");
    printf("MASSES:  %.5e   %.5e\n",mass,mass);
    printf("SOURCE: ETA_L\n");
    printf("SINKS: ETA_L_CONN  ETA_NL_CONN\n");
    for(j=0;j<nt;j++){
      printf("%d %e %e %e %e\n",j,
      props[ETA_L_ETA_L_CONN][j].real,
      props[ETA_L_ETA_L_CONN][j].imag,
      props[ETA_L_ETA_NL_CONN][j].real,
      props[ETA_L_ETA_NL_CONN][j].imag);
    }
    printf("ENDPROP\n");

    //printf("STARTPROP\n");
    //printf("MASSES:  %.5e   %.5e\n",mass,mass);
    //printf("SOURCE: ETA_NL\n");
    //printf("SINKS: ETA_L_CONN ETA_NL_CONN\n");
    //for(j=0;j<nt;j++){
      //printf("%d %e %e %e %e\n",j,
      //props[ETA_NL_ETA_L_CONN][j].real,
      //props[ETA_NL_ETA_L_CONN][j].imag,
      //props[ETA_NL_ETA_NL_CONN][j].real,
      //props[ETA_NL_ETA_NL_CONN][j].imag);
    //}
    //printf("ENDPROP\n");

    /* disconnected parts, by time slice */
    printf("STARTDISC\n");
    printf("MASSES:  %.5e   %.5e\n",mass,mass);
    printf("SOURCE: RANDOMWALL\n");
    printf("SINKS: SIGMA_DISC  ETA_L_DISC  ETA_NL_DISC\n");
    for(j=0;j<nt;j++){
      printf("%d %e %e %e %e %e %e\n",j,
      props[SIGMA_DISC][j].real,
      props[SIGMA_DISC][j].imag,
      props[ETA_L_DISC][j].real,
      props[ETA_L_DISC][j].imag,
      props[ETA_NL_DISC][j].real,
      props[ETA_NL_DISC][j].imag);
    }
    printf("ENDDISC\n");
  } // end node 0 

  free(props[0]); free(props); free(temp);
  free(R); free(R1); free(R2); free(R3);
  free(X); free(X_sigma); free(X_eta_l); free(X_eta_nl);
  
  return(cgn);
} /* spectrum_singlets */



/* "Multiply by" taste singlet 0++ operator.  */
/* mostly this is a dummy for testing - the operator is "one" */
void mult_sigma(su3_vector *src, su3_vector *dest, su3_vector *temp ) {
  register int i;
  register site *s;
  FORALLSITES(i,s){ dest[i]=src[i]; }
}

/* "Multiply by" the three link taste singlet pion operator.  */
void mult_eta_l(su3_vector *src, su3_vector *dest, su3_vector *temp, su3_matrix *smearlink ) {
  register int i;
  register site *s;
  int c ;
  /* All permutation with appropriate sign */
  struct {
    int d[3];
    Real sign ;
  } p[6]={{{0,1,2},+1.0/6.0},
	  {{1,2,0},+1.0/6.0},
	  {{2,0,1},+1.0/6.0},
	  {{0,2,1},-1.0/6.0},
	  {{1,0,2},-1.0/6.0},
	  {{2,1,0},-1.0/6.0}}; /* The factor of 6 accounts for the *
				* multiplicity of the permutations */
  
  /*clean up dest */
  FORALLSITES(i,s){
    clearvec( &(dest[i]) );
  }    
  for(c=0;c<6;c++)
    {
      //zeta_shift_smear(3,p[c].d,src,temp,smearlink) ; // Eric Gregory 3/1/05
      eta_shift_smear(3,p[c].d,src,temp,smearlink) ;
      FORALLSITES(i,s){
	scalar_mult_sum_su3_vector( &(dest[i]),
				   &(temp[i]), p[c].sign );
      }
    }
}

/* "Multiply by" the four link  nonlocal-time taste singlet pion operator.  */
void mult_eta_nl(su3_vector *src, su3_vector *dest, su3_vector *temp ) {
  register int i;
  register site *s;
  //NOT WRITTEN YET, set to zero
  FORALLSITES(i,s) clearvec( &(dest[i]) );
}


static su3_matrix *smear_links()
{
  su3_matrix *smearlink, *smearlink0;
  int i,j, dir;
  site *s;

  rephase(OFF);
  if(this_node==0)printf("SMEARING IS ON, level %d\n",SMEAR);
  smearlink = (su3_matrix *)malloc(sites_on_node*4*sizeof(su3_matrix));
  if(smearlink == NULL){
    printf("spectrum_singlets: No room for smeared links\n");
    terminate(1);
  }
  smearlink0 = (su3_matrix *)malloc(sites_on_node*4*sizeof(su3_matrix));
  if(smearlink0 == NULL){
    printf("spectrum_singlets: No room for smeared links\n");
    terminate(1);
  }
  if(SMEAR>=1){
    node0_printf( "APE Smearing: staple_weight = %e\n", STAPLE_WEIGHT);
    FORALLSITES(i,s)for(dir=XUP;dir<=TUP;dir++){
      smearlink0[4*i+dir] = s->link[dir];
    }
    ape_smear_field( smearlink0, smearlink, STAPLE_WEIGHT, u0, 1, APE_NHIT, APE_TOL);
  }
  else FORALLSITES(i,s)for(dir=XUP;dir<=TUP;dir++){
    FORALLSITES(i,s)for(dir=XUP;dir<=TUP;dir++){
      smearlink[4*i+dir] = s->link[dir];
    }
  }
  for(j=2;j<=SMEAR;j++){
    FORALLSITES(i,s)for(dir=XUP;dir<=TUP;dir++){
      smearlink0[4*i+dir] = smearlink[4*i+dir];
    }
    ape_smear_field( smearlink0, smearlink, STAPLE_WEIGHT, u0, 1, APE_NHIT, APE_TOL);
  }
  free(smearlink0);
  rephase(ON);

  return smearlink;
}

/******************************** flavor_ops_smear ***************************/
/* This is "flavor_ops.c", except that we use the APE smeared links for
 * parallel transport.  Also, field offsets have been replace by temp vectors */
/* MIMD version 7*/
/* Implementation of the flavor (\Xi_\mu) operators. 
   See Golterman & Smit Nulc. Phys. B245 (1984) 61.
   They are used for constructing the non-local pion sources
   and in the measurement of all the componenets of \bar{\psi}\psi
*/

/* This is the definition of the epsilon symbol */
static struct {
  int d[4];
  Real sign ;
} eps[24]={
  {{0,1,2,3},+1.0},
  {{0,3,1,2},+1.0},
  {{0,2,3,1},+1.0},
  {{0,3,2,1},-1.0},
  {{0,1,3,2},-1.0},
  {{0,2,1,3},-1.0},

  {{1,0,2,3},-1.0},
  {{1,3,0,2},-1.0},
  {{1,2,3,0},-1.0},
  {{1,3,2,0},+1.0},
  {{1,0,3,2},+1.0},
  {{1,2,0,3},+1.0},

  {{2,1,0,3},-1.0},
  {{2,3,1,0},-1.0},
  {{2,0,3,1},-1.0},
  {{2,3,0,1},+1.0},
  {{2,1,3,0},+1.0},
  {{2,0,1,3},+1.0},

  {{3,1,2,0},-1.0},
  {{3,0,1,2},-1.0},
  {{3,2,0,1},-1.0},
  {{3,0,2,1},+1.0},
  {{3,1,0,2},+1.0},
  {{3,2,1,0},+1.0}
};

/* Apply the symmetric shift opperator in direction "dir" *
 * This is the explicit version                           *
 * Covariant shifts are used                              */
void sym_shift_smear(int dir, su3_vector * src, su3_vector * dest,
		     su3_matrix *smearlink) {
  register int i ;
  register site *s ;
  msg_tag *tag[2];
  su3_vector *tvec;
  
  tvec = (su3_vector *)malloc( sites_on_node*sizeof(su3_vector) );

  tag[0] = start_gather_field( src, sizeof(su3_vector), dir, EVENANDODD ,gen_pt[0] );
  FORALLSITES(i,s)
    {
      mult_adj_su3_mat_vec( &smearlink[4*i+dir], &(src[i]),
			    &(tvec[i]) ) ;
    }
  tag[1] = start_gather_field(tvec, sizeof(su3_vector), OPP_DIR(dir), 
				  EVENANDODD ,gen_pt[1] );
  wait_gather(tag[0]);
  FORALLSITES(i,s)
    {
    mult_su3_mat_vec( &smearlink[4*i+dir], (su3_vector *)gen_pt[0][i], 
		      &(dest[i]) ) ;    
    }

  wait_gather(tag[1]);
  FORALLSITES(i,s)
    {
      add_su3_vector( &(dest[i]), (su3_vector *)gen_pt[1][i], &(dest[i]) ) ;    
    }
  /* Now divide by 2 eq. (4.2b) of Golterman's Meson paper*/
 FORALLSITES(i,s)
   {
     scalar_mult_su3_vector( &(dest[i]), 0.5, &(dest[i]) );
   }
  for(i=0;i<2;i++) cleanup_gather(tag[i]) ;
  free(tvec);
}

/* It applies the symmetric shift with directions                       *
 * stored in the array d. Each shift is multiplied by \zeta_k           *
 * n is the number of shifts                                            *
 * This is the E_\mu(x,y)=\Xi_\mu operator defined by Golterman.        *
 * Nucl. Phys. B245  eq.3.5 and eq. 4.2b                                */
void zeta_shift_smear(int n, int *d, su3_vector * src, su3_vector * dest,
		      su3_matrix *smearlink) {
  register int i,c ;
  short coords[4];
  register site *s;
  register Real sign ;
  su3_vector *ltemp;
  ltemp = (su3_vector *)malloc(sites_on_node*sizeof(su3_vector));

  for(c=0;c<n;c++)
    {  
      /* Do the shift in d[c] */ 
      if(c==0){
	/* first time from source */
	sym_shift_smear( d[c], src, ltemp, smearlink );
      }
      else{
	/* second time from dest */
	sym_shift_smear( d[c], dest, ltemp, smearlink );
      }
      /* Multiply by \zeta_d[c]. Because the phases are               *
       * on we multiply by \zeta * \eta = \epsilon * (-1)^coord[d[c]] */
      FORALLSITES(i,s){
	coords[XUP] = s->x; coords[YUP] = s->y;
	coords[ZUP] = s->z; coords[TUP] = s->t;
	/* The \epsilon */
	if (s->parity==EVEN) 
	  sign =  1.0 ;
	else
	  sign = -1.0 ;
	/* And the (-1)^coord[d[c]] */
	if(coords[d[c]]%2==1) sign=-sign ;
	scalar_mult_su3_vector( &(ltemp[i]), sign, &(dest[i]) );
      }
    }
    free(ltemp);
}

/* It applies the symmetric shift with directions                       *
 * stored in the array d. Each shift is multiplied by \eta_k            *
 * In fact since \eta_k are already absorbed into the U matrices do     *
 * no action is needed. Just do the symmetric shift.                    *
 * n is the number of shifts                                            */
void eta_shift_smear(int n, int *d, su3_vector * src, su3_vector * dest,
		     su3_matrix *smearlink) {
  int c ;
  su3_vector *temp0, *temp1, *tmp ;
  temp0 = (su3_vector *)malloc(sites_on_node*sizeof(su3_vector));
  temp1 = (su3_vector *)malloc(sites_on_node*sizeof(su3_vector));
  
  if(n==1)
    {
      sym_shift_smear(d[0], src, dest, smearlink);
    }
  else
    {
      sym_shift_smear(d[0], src, temp0, smearlink );
      for(c=1;c<n-1;c++)
	{  
	  sym_shift_smear(d[c], temp0, temp1, smearlink);
	  /* switch the pointers */
	  tmp = temp0   ;
	  temp0 = temp1 ;
	  temp1 = tmp   ;
	}
      /* do the last shift */
      sym_shift_smear(d[n-1], temp0, dest, smearlink );
    }
    free(temp0);
    free(temp1);
}

/* Multiply by the Xi_mu flavor operator */
void mult_flavor_vector_smear(int mu, su3_vector * src, su3_vector * dest, su3_matrix *smearlink ) {
  int d[1] ;
  
  d[0] = mu ;
  zeta_shift_smear(1,d,src,dest,smearlink) ;
}

/* Multiply by the 1/2(Xi_mu Xi_nu - Xi_nu Xi_mu)flavor operator */
void mult_flavor_tensor_smear(int mu, int nu, su3_vector * src, su3_vector * dest, su3_matrix *smearlink ) { 
  register int i;
  register site *s;
  int d[2];
  su3_vector *ltemp;
  ltemp = (su3_vector *)malloc(sites_on_node*sizeof(su3_vector));

  d[0] = mu ; d[1]=nu ;
  zeta_shift_smear(2, d,src,dest,smearlink) ;
 
  d[0] = nu ; d[1]=mu ;
  zeta_shift_smear(2,d,src,ltemp,smearlink) ;

  
  FORALLSITES(i,s){
    scalar_mult_sum_su3_vector( &(dest[i]), &(ltemp[i]), -1.0);
    scalar_mult_su3_vector( &(dest[i]), 0.5, &(dest[i]) );
  }
  free(ltemp);
}

/* Multiply by the Xi_mu Xi_5 flavor operator */
void mult_flavor_pseudovector_smear(int mu, su3_vector * src, su3_vector * dest, su3_matrix *smearlink ) {
  register int i;
  register site *s;
  int p ; 
  su3_vector *ltemp;
  ltemp=(su3_vector *)malloc(sites_on_node*sizeof(su3_vector));

  /*clean up dest */
  FORALLSITES(i,s){
    clearvec( &(dest[i]) );
  }   
   for(p=0;p<24;p++)
      if(eps[p].d[0]==mu)
	{
	  zeta_shift_smear( 3, &eps[p].d[1], src, &(ltemp[i]),smearlink ) ;
	  /* Multiply the extra 1/6 needed by the definition    *
	   * of the operator (number of permutations)           */
	  FORALLSITES(i,s){
	    scalar_mult_sum_su3_vector( &(dest[i]), &(ltemp[i]), eps[p].sign/6.0);
	  }
	}
   free(ltemp);
}

/* Multiply by the Xi_5 flavor operator */
void mult_flavor_pseudoscalar_smear(su3_vector * src, su3_vector * dest,
				    su3_matrix *smearlink) { 
  register int i;
  register site *s;
  int p ; 
  su3_vector *ltemp;
  ltemp=(su3_vector *)malloc(sites_on_node*sizeof(su3_vector));
  
  /*clean up dest */
  FORALLSITES(i,s){
    clearvec( &(dest[i]) );
  }   
  for(p=0;p<24;p++) {
       zeta_shift_smear( 4, eps[p].d, src, &(ltemp[i]),smearlink ) ;
       /*  Multiply the the extra 1/24 needed by the            *
	* definition of the operator (number of permutations)   */
       FORALLSITES(i,s){
	 scalar_mult_sum_su3_vector( &(dest[i]), &(ltemp[i]), eps[p].sign/24.0);
       }
     } 
  free(ltemp);
}

/* Multiply by the Gamma_mu spin operator */
void mult_spin_vector_smear(int mu, su3_vector * src, su3_vector * dest, su3_matrix *smearlink ) {
  int d[1] ;
  
  d[0] = mu ;
  eta_shift_smear(1,d,src,dest,smearlink) ;
}

/* Multiply by the 1/2(Gamma_mu Gamma_nu - Gamma_nu gamma_mu) spin operator */
void mult_spin_tensor_smear(int mu, int nu, su3_vector * src, su3_vector * dest, su3_matrix *smearlink ) { 
  register int i;
  register site *s;
  int d[2];
  su3_vector *ltemp;
  ltemp = (su3_vector *)malloc(sites_on_node*sizeof(su3_vector));

  d[0] = mu ; d[1]=nu ;
  eta_shift_smear( 2, d, src, dest, smearlink) ;
 
  d[0] = nu ; d[1]=mu ;
  eta_shift_smear( 2, d, src, ltemp, smearlink ) ;

  FORALLSITES(i,s){
    scalar_mult_sum_su3_vector( &(dest[i]), &(ltemp[i]), -1.0);
    scalar_mult_su3_vector( &(dest[i]), 0.5, &(dest[i]) );
  }
  free(ltemp);
}

/* Multiply by the Gamma_mu Gamma_5 spin operator */
void mult_spin_pseudovector_smear(int mu, su3_vector * src, su3_vector * dest, 
				  su3_matrix *smearlink ) {
  register int i;
  register site *s;
  int p ; 
  su3_vector *ltemp;
  ltemp = (su3_vector *)malloc(sites_on_node*sizeof(su3_vector));

  /*clean up dest */
  FORALLSITES(i,s){
    clearvec( &(dest[i]) );
  }   
   for(p=0;p<24;p++)
      if(eps[p].d[0]==mu)
	{
	  eta_shift_smear( 3, &eps[p].d[1], src, &(ltemp[i]), smearlink ) ;
	  /* Multiply the extra 1/6 needed by the definition    *
	     of the operator (number of permutations)           */
	  FORALLSITES(i,s){
	    scalar_mult_sum_su3_vector( &(dest[i]), &(ltemp[i]), eps[p].sign/6.0);
	  }
	}
   free(ltemp);
}

/* Multiply by the Gamma_5 spin operator */
void mult_spin_pseudoscalar_smear(su3_vector * src, su3_vector * dest, su3_matrix *smearlink ) { 
  register int i;
  register site *s;
  int p ; 
  su3_vector *ltemp;
  ltemp = (su3_vector *)malloc(sites_on_node*sizeof(su3_vector));
  
  /*clean up dest */
  FORALLSITES(i,s){
    clearvec( &(dest[i]) );
  }   
  for(p=0;p<24;p++)
     {
       eta_shift_smear( 4, eps[p].d, src, &(ltemp[i]), smearlink ) ;
       /*  Multiply the the extra 1/24 needed by the            *
	* definition of the operator (number of permutations)   */
       FORALLSITES(i,s){
	 scalar_mult_sum_su3_vector( &(dest[i]), &(ltemp[i]), eps[p].sign/24.0);
       }
     } 
  free(ltemp);
}

/* Put KS phases into smeared links */
void rephase_smear(su3_matrix *smearlink ){
register int i,j,k,dir;
register site *s;
    FORALLSITES(i,s){
	for(dir=XUP;dir<=TUP;dir++){
	    for(j=0;j<3;j++)for(k=0;k<3;k++){
		smearlink[4*i+dir].e[j][k].real *= s->phase[dir];
		smearlink[4*i+dir].e[j][k].imag *= s->phase[dir];
	    }
	}
    }
} /* rephase_smear.c */
