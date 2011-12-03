/******** gluon_prop.c *********/
/* MIMD version 7 */
/* Initial version by UMH 12/21/00 */

/* Computes the gluon propagator in momentum space for the given
   (gauge fixed) gauge field configuration */

#include "gluon_prop_includes.h"

void gluon_prop( void ) {
register int i,dir;
register int pmu;
register site *s;
anti_hermitmat ahtmp;
Real pix, piy, piz, pit;
Real sin_pmu, prop_s = 0., prop_l = 0., ftmp1, ftmp2;
complex ctmp;
su3_matrix mat;
struct {
  Real f1, f2;
} msg;
double trace, dmuAmu;
int px, py, pz, pt;
int currentnode,newnode;

    pix = PI / (Real)nx;
    piy = PI / (Real)ny;
    piz = PI / (Real)nz;
    pit = PI / (Real)nt;

    trace = 0.0;
    /* Make A_mu as anti-hermition traceless part of U_mu */
    /* But store as SU(3) matrix for call to FFT */
    for(dir=XUP; dir<=TUP; dir++)
    {
	FORALLSITES(i,s){
	    trace += (double)(trace_su3( &(s->link[dir]))).real;
	    make_anti_hermitian( &(s->link[dir]), &ahtmp);
	    uncompress_anti_hermitian( &ahtmp, &(s->a_mu[dir]));
	}

	g_sync();
	/* Now Fourier transform */
	restrict_fourier_site(F_OFFSET(a_mu[dir]),
			      sizeof(su3_matrix), FORWARDS);
    }
    g_doublesum( &trace);
    trace /= (double)(12*volume);

    /* Correct for mid-link definition of A_mu(p), multiply A_mu(p) with */
    /* sin(p_mu/2) and add. Also accumulate sum_mu sin^2(p_mu/2) */
    for(dir=XUP; dir<=TUP; dir++)
    {
	FORALLSITES(i,s){
	    switch(dir){
		case XUP:
		    pmu = s->x;
		    sin_pmu = sin((double)(pmu*pix));
		    s->tempfloat = sin_pmu * sin_pmu;
		    ctmp.real = cos((double)(pmu*pix));
		    ctmp.imag = -sin_pmu;
		    c_scalar_mult_su3mat( &(s->a_mu[dir]), &ctmp, &mat);
		    su3mat_copy( &mat, &(s->a_mu[dir]));
		    scalar_mult_su3_matrix( &(s->a_mu[dir]), sin_pmu,
			&(s->tempmat1));
		    break;

		case YUP:
		    pmu = s->y;
		    sin_pmu = sin((double)(pmu*piy));
		    s->tempfloat += sin_pmu * sin_pmu;
		    ctmp.real = cos((double)(pmu*piy));
		    ctmp.imag = -sin_pmu;
		    c_scalar_mult_su3mat( &(s->a_mu[dir]), &ctmp, &mat);
		    su3mat_copy( &mat, &(s->a_mu[dir]));
		    scalar_mult_add_su3_matrix( &(s->tempmat1), 
			&(s->a_mu[dir]), sin_pmu, &(s->tempmat1));
		    break;

		case ZUP:
		    pmu = s->z;
		    sin_pmu = sin((double)(pmu*piz));
		    s->tempfloat += sin_pmu * sin_pmu;
		    ctmp.real = cos((double)(pmu*piz));
		    ctmp.imag = -sin_pmu;
		    c_scalar_mult_su3mat( &(s->a_mu[dir]), &ctmp, &mat);
		    su3mat_copy( &mat, &(s->a_mu[dir]));
		    scalar_mult_add_su3_matrix( &(s->tempmat1), 
			&(s->a_mu[dir]), sin_pmu, &(s->tempmat1));
		    break;

		case TUP:
		    pmu = s->t;
		    sin_pmu = sin((double)(pmu*pit));
		    s->tempfloat += sin_pmu * sin_pmu;
		    ctmp.real = cos((double)(pmu*pit));
		    ctmp.imag = -sin_pmu;
		    c_scalar_mult_su3mat( &(s->a_mu[dir]), &ctmp, &mat);
		    su3mat_copy( &mat, &(s->a_mu[dir]));
		    scalar_mult_add_su3_matrix( &(s->tempmat1), 
			&(s->a_mu[dir]), sin_pmu, &(s->tempmat1));
		    break;

		default: printf("BOTCH: bad direction\n"); exit(1);
	    }

	}
    }

    dmuAmu = 0.0;
    /* Now get the longitudinal part and then the scalar part of
       the gluon propagator */
    FORALLSITES(i,s){
	/* The longitudinal part [*(N^2-1)] */
	if( s->x == 0 && s->y == 0 && s->z == 0 && s->t == 0 ){
	    ftmp1 = 0.0;
	}
	else{
	    ftmp1 = realtrace_su3( &(s->tempmat1), &(s->tempmat1));
	    dmuAmu += (double)ftmp1;
	    ftmp1 /= s->tempfloat;
	}
	s->long_prop = ftmp1 / (Real)(12*volume);

	/* The scalar part */
	ftmp2 = realtrace_su3( &(s->a_mu[XUP]), &(s->a_mu[XUP]));
        for(dir=YUP; dir<=TUP; dir++)
	    ftmp2 += realtrace_su3( &(s->a_mu[dir]), &(s->a_mu[dir]));
	/* With normalization from (d-1)*(N^2-1)/2 */
#if 1
	s->scalar_prop = (ftmp2 - ftmp1) / (Real)(12*volume);
#else
	s->scalar_prop = ftmp2 / (Real)(12*volume);
#endif
    }
    g_doublesum( &dmuAmu);
    dmuAmu /= (double)(volume/4);
    dmuAmu /= (double)(3*volume);

    /* Finally node 0 averages over permutations of spatial momentum
       components and writes the gluon propagator */
    if(this_node==0)printf("START_GLUON_PROP\n");
    g_sync();
    currentnode=0;

    for(pt=0;pt<=nt/2;pt++){
      for(px=0;px<=nx/2;px++)for(py=0;py<=px;py++)for(pz=0;pz<=py;pz++){

	newnode = node_number(px,py,pz,pt);
	if(newnode != currentnode){	/* switch to another node */
	  /* tell newnode it's OK to send */
	  if( this_node==0 && newnode!=0 )send_field((char *)&ftmp2,4,newnode);
	  if( this_node==newnode && newnode!=0 )get_field((char *)&ftmp2,4,0);
	  currentnode = newnode;
	}

	if(this_node==0){
	  if(currentnode==0){
	    i = node_index(px,py,pz,pt);
	    prop_s = lattice[i].scalar_prop;
	    prop_l = lattice[i].long_prop;
	  }
	  else{
	    get_field((char *)&msg,sizeof(msg),currentnode);
	    prop_s = msg.f1;
	    prop_l = msg.f2;
	  }
	}
	else{	/* for nodes other than 0 */
	  if(this_node==currentnode){
	    i = node_index(px,py,pz,pt);
	    msg.f1 = lattice[i].scalar_prop;
	    msg.f2 = lattice[i].long_prop;
	    send_field((char *)&msg,sizeof(msg),0);
	  }
	}

	/* 2nd permutation */
	newnode = node_number(py,px,pz,pt);
	if(newnode != currentnode){	/* switch to another node */
	  /* tell newnode it's OK to send */
	  if( this_node==0 && newnode!=0 )send_field((char *)&ftmp2,4,newnode);
	  if( this_node==newnode && newnode!=0 )get_field((char *)&ftmp2,4,0);
	  currentnode = newnode;
	}

	if(this_node==0){
	  if(currentnode==0){
	    i = node_index(py,px,pz,pt);
	    prop_s += lattice[i].scalar_prop;
	    prop_l += lattice[i].long_prop;
	  }
	  else{
	    get_field((char *)&msg,sizeof(msg),currentnode);
	    prop_s += msg.f1;
	    prop_l += msg.f2;
	  }
	}
	else{	/* for nodes other than 0 */
	  if(this_node==currentnode){
	    i = node_index(py,px,pz,pt);
	    msg.f1 = lattice[i].scalar_prop;
	    msg.f2 = lattice[i].long_prop;
	    send_field((char *)&msg,sizeof(msg),0);
	  }
	}

	/* 3rd permutation */
	newnode = node_number(pz,py,px,pt);
	if(newnode != currentnode){	/* switch to another node */
	  /* tell newnode it's OK to send */
	  if( this_node==0 && newnode!=0 )send_field((char *)&ftmp2,4,newnode);
	  if( this_node==newnode && newnode!=0 )get_field((char *)&ftmp2,4,0);
	  currentnode = newnode;
	}

	if(this_node==0){
	  if(currentnode==0){
	    i = node_index(pz,py,px,pt);
	    prop_s += lattice[i].scalar_prop;
	    prop_l += lattice[i].long_prop;
	  }
	  else{
	    get_field((char *)&msg,sizeof(msg),currentnode);
	    prop_s += msg.f1;
	    prop_l += msg.f2;
	  }
	}
	else{	/* for nodes other than 0 */
	  if(this_node==currentnode){
	    i = node_index(pz,py,px,pt);
	    msg.f1 = lattice[i].scalar_prop;
	    msg.f2 = lattice[i].long_prop;
	    send_field((char *)&msg,sizeof(msg),0);
	  }
	}

	/* 4th permutation */
	newnode = node_number(px,pz,py,pt);
	if(newnode != currentnode){	/* switch to another node */
	  /* tell newnode it's OK to send */
	  if( this_node==0 && newnode!=0 )send_field((char *)&ftmp2,4,newnode);
	  if( this_node==newnode && newnode!=0 )get_field((char *)&ftmp2,4,0);
	  currentnode = newnode;
	}

	if(this_node==0){
	  if(currentnode==0){
	    i = node_index(px,pz,py,pt);
	    prop_s += lattice[i].scalar_prop;
	    prop_l += lattice[i].long_prop;
	  }
	  else{
	    get_field((char *)&msg,sizeof(msg),currentnode);
	    prop_s += msg.f1;
	    prop_l += msg.f2;
	  }
	}
	else{	/* for nodes other than 0 */
	  if(this_node==currentnode){
	    i = node_index(px,pz,py,pt);
	    msg.f1 = lattice[i].scalar_prop;
	    msg.f2 = lattice[i].long_prop;
	    send_field((char *)&msg,sizeof(msg),0);
	  }
	}

	/* 5th permutation */
	newnode = node_number(py,pz,px,pt);
	if(newnode != currentnode){	/* switch to another node */
	  /* tell newnode it's OK to send */
	  if( this_node==0 && newnode!=0 )send_field((char *)&ftmp2,4,newnode);
	  if( this_node==newnode && newnode!=0 )get_field((char *)&ftmp2,4,0);
	  currentnode = newnode;
	}

	if(this_node==0){
	  if(currentnode==0){
	    i = node_index(py,pz,px,pt);
	    prop_s += lattice[i].scalar_prop;
	    prop_l += lattice[i].long_prop;
	  }
	  else{
	    get_field((char *)&msg,sizeof(msg),currentnode);
	    prop_s += msg.f1;
	    prop_l += msg.f2;
	  }
	}
	else{	/* for nodes other than 0 */
	  if(this_node==currentnode){
	    i = node_index(py,pz,px,pt);
	    msg.f1 = lattice[i].scalar_prop;
	    msg.f2 = lattice[i].long_prop;
	    send_field((char *)&msg,sizeof(msg),0);
	  }
	}

	/* 6th permutation */
	newnode = node_number(pz,px,py,pt);
	if(newnode != currentnode){	/* switch to another node */
	  /* tell newnode it's OK to send */
	  if( this_node==0 && newnode!=0 )send_field((char *)&ftmp2,4,newnode);
	  if( this_node==newnode && newnode!=0 )get_field((char *)&ftmp2,4,0);
	  currentnode = newnode;
	}

	if(this_node==0){
	  if(currentnode==0){
	    i = node_index(pz,px,py,pt);
	    prop_s += lattice[i].scalar_prop;
	    prop_l += lattice[i].long_prop;
	  }
	  else{
	    get_field((char *)&msg,sizeof(msg),currentnode);
	    prop_s += msg.f1;
	    prop_l += msg.f2;
	  }
	}
	else{	/* for nodes other than 0 */
	  if(this_node==currentnode){
	    i = node_index(pz,px,py,pt);
	    msg.f1 = lattice[i].scalar_prop;
	    msg.f2 = lattice[i].long_prop;
	    send_field((char *)&msg,sizeof(msg),0);
	  }
	}

	if(this_node==0){
	  prop_s /= (6.0*trace*trace);
	  prop_l /= (6.0*trace*trace);
	  printf(" %e %e\n", (double)prop_s, (double)prop_l);
	}
      }
    }
    if(this_node==0){
	printf("END_GLUON_PROP\n");
	printf("TRACE_U %e\n", trace);
	printf("DMU_AMU %e\n", sqrt(dmuAmu));
    }

} /* gluon_prop */
