/******** nl_spectrum.c *************/
/* MIMD version 7 */
/* C.E.D.  7/91  first verion
* D.T. modified 5/97 for improved action program, changed source
*   slice logic, use gen_pt[8-15] instead of gen_pt2[].
* For measuring propagation IN THE T DIRECTION
* With correct t offset in E_MES_PRO.  Now reported as E_PI_PRO. C.D.7/17/91
* modified 2/99 for new output format, DT
*
* Spectrum for Kogut-Susskind nonlocal hadrons, including delta
* pion, and a few local ones for checking.  Uses E and O wall sources
* This version does arbitrary number of wall sources 
* This version DOES NOT FIX THE GAUGE.  MUST BE DONE BEFORE CALLING!!
*/

#include "generic_ks_includes.h"
#include "../include/dslash_ks_redefine.h"

#define NL_PI_DIR ZUP  /* Defines direction for non-local pi propagator */
                       /* This direction depends on the choice of KS */
                       /* phase! */


short get_coord(int i,int dir)
     /* Returns the lattice coordinate at site i in the direction dir */
{
  switch (dir)
    {
    case XUP:   return(lattice[i].x);
    case YUP:   return(lattice[i].y);
    case ZUP:   return(lattice[i].z);
    case TUP:   return(lattice[i].t);
    case XDOWN: return(lattice[i].x);
    case YDOWN: return(lattice[i].y);
    case ZDOWN: return(lattice[i].z);
    case TDOWN: return(lattice[i].t);
    }
  return -1;  /* Error */
}

void accum_delta_prop(int i,int c[3], field_offset tempmat1)
     /* Load a determinant for the calculation of the delta propagator */
     /* Each column of the determinant is the symmetrically shifted   */
     /* "q" propagator.  Each column has a differenct shift direction */
     /* Each column also corresponds to a different source wall color */
     /* as specified by the color vector c                            */
     /* Result is put in tempmat1 */
{
  int j=-1,dir;
  su3_matrix *af, *ab;
  register su3_matrix *tmp1;
  
  
  FORALLUPDIRBUT(TUP,dir)
    {
      /* For timeward propagation, this step is pedantic */
      /* since we could simply take j = dir */
      /* but for spacelike propagation, it is necessary */
      switch (dir)
	{
	case XUP: j = 0; break;
	case YUP: j = 1; break;
	case ZUP: j = 2;
	}

      /* gen_pt points to the "q" propagator on the next neighbor    */
      af = (su3_matrix *)gen_pt[dir][i];
      ab = (su3_matrix *)gen_pt[dir+4][i];
      tmp1 = (su3_matrix *)F_PT(&lattice[i],tempmat1);
      su3vec_copy( (su3_vector *)(af->e[c[j]]),
		   (su3_vector *)(tmp1->e[j]) );
      add_su3_vector( (su3_vector *)(ab->e[c[j]]),
		     (su3_vector *)(tmp1->e[j]),
		     (su3_vector *)(tmp1->e[j]) );

    }
}

void delta_prop(int c0,int c1,int c2,int perm,int t,double *delprop,
		field_offset tempmat1)
     /* Calculates the contribution to the delta propagator at time t */
     /* With source colors c0, c1, c2 */
     /* perm specifies the sign of the contribution to be calculated */
{
  register int i;
  register complex cc;
  int x,y,z;
  int c[3];

  c[0] = c0;
  c[1] = c1;
  c[2] = c2;

  for(x=0;x<nx;x+=2)for(y=0;y<ny;y+=2)for(z=0;z<nz;z+=2) {
      if( node_number(x,y,z,t) != mynode() )continue;
      i=node_index(x,y,z,t);

      /*  Calculate for the cube origin only */

      accum_delta_prop(i,c,tempmat1);

      cc = det_su3( (su3_matrix *)F_PT(&lattice[i],tempmat1) );

      if(perm>0)*delprop += cc.real;
      else      *delprop -= cc.real;  

    }
}

void accum_nl_meson_prop(int i,int dir,
			 field_offset destq, field_offset desto)
     /* Apply symmetric shift to "q" and "o" propagators */
     /* for all source wall colors */
     /* at site i in direction dir */
     /* destq <- D_dir q; desto <- D_dir o */
     /* destq and desto must be of size su3_matrix */
{
  su3_matrix *af, *ab;
  register su3_matrix *dstq, *dsto;
  
  /* gen_pt[0+dir] points to the "q" propagator       */
  /* destq <- D_dir q                       */
  af = (su3_matrix *)gen_pt[dir][i];
  ab = (su3_matrix *)gen_pt[dir+4][i];
  dstq = (su3_matrix *)F_PT(&lattice[i],destq);
  su3mat_copy( af, dstq );
  add_su3_matrix( ab, dstq, dstq);

  /* gen_pt[8+dir] points to the "o" propagator      */
  /* desto <- D_dir o                       */
  af = (su3_matrix *)gen_pt[8+dir][i];
  ab = (su3_matrix *)gen_pt[8+dir+4][i];
  dsto = (su3_matrix *)F_PT(&lattice[i],desto);
  su3mat_copy( af, dsto);
  add_su3_matrix( ab, dsto, dsto);

}

void nl_meson_prop (int t,double *nlpiprop,double *nlpi2prop,
		    double *ckpiprop,double *ckpi2prop,
		    field_offset tempmat1, field_offset tempmat2)
     /* Calculates non-local pion propagator pi_3 and pi_3 tilde */
     /* and local pion propagators for a check */
     /* tempmat1 and tempmat2 are scratch space of size su3_matrix each */
{
  int x,y,z,i,icol;
  int coord;
  register su3_matrix *tmp1, *tmp2;

  complex cc;
  /* Sum over all x,y z              */
  for(x=0;x<nx;x++)for(y=0;y<ny;y++)for(z=0;z<nz;z++) {
    if( node_number(x,y,z,t) != mynode() )continue;
    i=node_index(x,y,z,t);
    
    coord = get_coord(i,NL_PI_DIR);
    if(coord%2 == 0)accum_nl_meson_prop(i,NL_PI_DIR,tempmat1,tempmat2);
    
    for(icol=0;icol<3;icol++) {
      
      /* Calculate non-local propagator only on even coordinate */
      
      if(coord%2 ==0) {
	/* propmat contains "q" and tempmat1 contains "Dq" */
	/* q^adj D q  */
	tmp1 = (su3_matrix *)F_PT(&lattice[i],tempmat1);
	cc = su3_dot( &lattice[i].propmat[icol],
		      (su3_vector *)(tmp1->e[icol]) );
	*nlpiprop += cc.real;
	
	/* propmat2 contains "o" and tempmat2 contains "Do" */
	/* q^adj D q - o^adj D o  */
	tmp2 = (su3_matrix *)F_PT(&lattice[i],tempmat2);
	cc = su3_dot( &(lattice[i].propmat2[icol]),
		      (su3_vector *)(tmp2->e[icol]) );
	*nlpiprop -= cc.real;
	
	/* o^adj D q          */
	cc = su3_dot( &lattice[i].propmat2[icol],
		      (su3_vector *)(tmp1->e[icol]) );
	
	if((x+y+z)%2 == 0)*nlpi2prop += cc.real;
	else              *nlpi2prop -= cc.real;
	
	/* o^adj D q - q^adj D o          */
	cc = su3_dot( &(lattice[i].propmat[icol]),
		      (su3_vector *)(tmp2 ->e[icol]) );
	
	if((x+y+z)%2 == 0)*nlpi2prop -= cc.real;
	else              *nlpi2prop += cc.real;
	
      }
      
      /* Calculate local check for all z */
      
      /* q^adj q  */
      cc = su3_dot( &lattice[i].propmat[icol],
		    &lattice[i].propmat[icol] );
      *ckpiprop += cc.real;
      
      /* q^adj q + o^adj o          */
      cc = su3_dot( &lattice[i].propmat2[icol],
		    &lattice[i].propmat2[icol] );
      *ckpiprop += cc.real;
      
      
      /* q^adj o  */
      cc = su3_dot( &lattice[i].propmat[icol],
		    &lattice[i].propmat2[icol] );
      
      if( (x+y+z)%2==0)*ckpi2prop += cc.real;
      else	           *ckpi2prop -= cc.real;
      
    }
  }
}


#include <assert.h>

int nl_spectrum( Real vmass, field_offset tempvec1, field_offset tempvec2,
		 field_offset tempmat1, field_offset tempmat2,
		 ferm_links_t *fn) { 
  /* return the C.G. iteration number */
  double *piprop,*pi2prop,*rhoprop,*rho2prop,*barprop;
  double *nlpiprop,*nlpi2prop,*ckpiprop,*ckpi2prop;
  double *delprop,*ckbarprop;
  Real vmass_x2;
  site* s;
  register complex cc;
  Real finalrsq;
  register int i,x,y,z,t,icol,cgn;
  register int t_source,t_off;
  int dir,isrc;
  register su3_matrix *tmp1;
  double proptmp;

  msg_tag *mtag[16];

  /* Required gen_pt dimension */
  assert(N_POINTERS >= 16);

  piprop = (double *)malloc( nt*sizeof(double) );
  pi2prop = (double *)malloc( nt*sizeof(double) );
  rhoprop = (double *)malloc( nt*sizeof(double) );
  rho2prop = (double *)malloc( nt*sizeof(double) );
  barprop = (double *)malloc( nt*sizeof(double) );
  nlpiprop = (double *)malloc( nt*sizeof(double) );
  nlpi2prop = (double *)malloc( nt*sizeof(double) );
  ckpiprop = (double *)malloc( nt*sizeof(double) );
  ckpi2prop = (double *)malloc( nt*sizeof(double) );
  delprop = (double *)malloc( nt*sizeof(double) );
  ckbarprop = (double *)malloc( nt*sizeof(double) );

  for( t=0; t<nt; t++ ){
    piprop[t]=0.0; pi2prop[t]=0.0; rhoprop[t]=0.0; rho2prop[t]=0.0;
    nlpiprop[t]=0.0; nlpi2prop[t]=0.0;
    ckpiprop[t]=0.0; ckpi2prop[t]=0.0;
    barprop[t]=0.0; delprop[t]=0.0; ckbarprop[t]=0.0;
  }

  vmass_x2 = 2.*vmass;
  cgn=0;

  /* Unlike spectrum.c, here we calculate only with wall sources */
  for(t_source=source_start, isrc=0; t_source<2*nt && isrc < n_sources;
        ++isrc, t_source += source_inc ) {
      
    if(this_node==0)printf("nl_spectrum(): source time = %d\n",t_source);
      /* Only work for even source slices */
      if( t_source%2 != 0 ){
	printf("DUMMY:  Use even time slices for nl_spectrum()\n");
	terminate(0);
      }

      /* Compute propagator from even wall sites */
      /* Sources are normalized to 1/8 to make them comparable to */
      /* propagators from a wall with ones on the cube origin. */
      /* Put result in propmat */
      
      for(icol=0; icol<3; icol++) {
	  
	  /* initialize tempvec1 and tempvec2 */
	  clear_latvec( tempvec1, EVEN);
	  clear_latvec( tempvec2, EVEN);
	  
	  for(x=0;x<nx;x++)for(y=0;y<ny;y++)for(z=0;z<nz;z++) {
	      if((x+y+z) % 2 == 0) {
		  if( node_number(x,y,z,t_source) != mynode() )continue;
		  i=node_index(x,y,z,t_source);
		  ((su3_vector *)(F_PT(&lattice[i],tempvec1)))->c[icol].real = 
		    -0.25;
		}
	    }
	  
	  /* do a C.G.: source tempvec1, result tempvec2 */
	  cgn += ks_congrad( tempvec1, tempvec2, vmass,
	    niter, nrestart, rsqprop, PRECISION, EVEN, &finalrsq, fn);
	  /* Multiply by -Madjoint */
	  dslash_site( tempvec2, tempvec2, ODD, fn);
	  scalar_mult_latvec( tempvec2, -vmass_x2, tempvec2, EVEN);
	  
	  /* fill the hadron matrix */
	  copy_latvec( tempvec2, F_OFFSET(propmat[icol]), EVENANDODD);
	} /* end loop on icol */
      
      
      /* Compute propagator from odd wall sites */
      /* Put result in propmat2 */
      
      for(icol=0; icol<3; icol++) {
	  
	  /* initialize tempvec1 and tempvec2 */
	  clear_latvec( tempvec1, ODD);
	  clear_latvec( tempvec2, ODD);
	  for(x=0;x<nx;x++)for(y=0;y<ny;y++)for(z=0;z<nz;z++) {
	      if((x+y+z) % 2 == 1) {
		  if( node_number(x,y,z,t_source) != mynode() )continue;
		  i=node_index(x,y,z,t_source);
		  ((su3_vector *)(F_PT(&lattice[i],tempvec1)))->c[icol].real = 
		    -0.25;
		}
	    }
	  
	  /* do a C.G.: source tempvec1, result tempvec2 */
	  cgn += ks_congrad( tempvec1, tempvec2, vmass,
	    niter, nrestart, rsqprop, PRECISION, ODD, &finalrsq, fn);
	  /* Multiply by -Madjoint */
	  dslash_site( tempvec2, tempvec2, EVEN, fn);
	  scalar_mult_latvec( tempvec2, -vmass_x2, tempvec2, ODD);
	  
	  /* fill the hadron matrix */
	  copy_latvec( tempvec2, F_OFFSET(propmat2[icol]), EVENANDODD);
	} /* end loop on icol */
      
      /* cgn now gives the sum for both inversions */
      
      
      /* measure the meson propagator for the E wall source */
      for(t=0; t<nt; t++) {
	  /* define the time value offset t from t_source */
	  t_off = (t+t_source)%nt;
	  
	  for(x=0;x<nx;x++)for(y=0;y<ny;y++)for(z=0;z<nz;z++)
	    for(icol=0;icol<3;icol++) {
		if( node_number(x,y,z,t_off) != mynode() )continue;
		i=node_index(x,y,z,t_off);
		cc = su3_dot( &lattice[i].propmat[icol],
			     &lattice[i].propmat[icol] );
		
		piprop[t] += cc.real;
                /* (rhoprop and rho2prop are not generated by this source) */
		
		if( (x+y+z)%2==0)pi2prop[t] += cc.real;
		else	     pi2prop[t] -= cc.real;
		
	      }
	  
	} /* nt-loop */
      
      /* measure the baryon propagator for the E wall source */
      for(t=0; t<nt; t++) {
	  /* define the time value offset t from t_source */
	  t_off = (t+t_source)%nt;

	  proptmp = 0;
	  for(x=0;x<nx;x+=2)for(y=0;y<ny;y+=2)for(z=0;z<nz;z+=2) {
	      if( node_number(x,y,z,t_off) != mynode() )continue;
	      i=node_index(x,y,z,t_off);
	      cc = det_su3( (su3_matrix *)(lattice[i].propmat) );
	      proptmp += cc.real;
	  }
	  
	  /* must get sign right.  This looks to see if we have
	     wrapped around the lattice.  "t" is the distance
	     from the source to the measurement, so we are
	     trying to find out if t_source+t is greater than
	     or equal to nt. */
	  /*if( (((t+t_source)/nt-t_source/nt)%2) == 1 )barprop[t] *= -1.0;*/
	  /* change sign because antiperiodic b.c.  sink point
	     should really be in a copy of the lattice */

	  if( (((t+t_source)/nt-t_source/nt)%2) == 0 )
	    barprop[t] += proptmp;
	  else
	    barprop[t] -= proptmp;

	} /* nt-loop */
      
      
      /* Measure nonlocal (and some local for checking) propagators    */
      /* These propagators include the delta and some nonlocal mesons  */
      /* The method for the delta is described in M.F.L. Golterman and */
      /* J. Smit, Nucl. Phys. B 255, 328 (1985)                        */
      /* Equation (6.3) defines the sink operator for the delta        */
      /* The method for the mesons is described in M.F.L. Golterman    */
      /* Nucl. Phys. B 273, 663 (1986)                                 */
      /* The treatment of the source wall is described in Gupta,       */
      /* Guralnik, Kilcup, and Sharpe, (GGKS) NSF-ITP-90-172 (1990)    */
      /* Phys.Rev.D43:2003-2026,1991  */
      /* To get the delta propagator, we take the "q" propagator       */
      /* matrices for each of the wall colors and antisymmetrize over  */
      /* wall color as well as s                    */
      
      /* First construct the "q" and "o" propagators                   */
      /* Put q = E + O in propmat and o = E - O in propmat2 */
      
      FORALLSITES(i,s) {
	tmp1 = (su3_matrix *)F_PT(s,tempmat1);
	for(icol=0; icol<3; icol++) {
	  add_su3_vector (&(s->propmat[icol]), &(s->propmat2[icol]), 
			  (su3_vector *)(tmp1->e[icol]) );
	  sub_su3_vector (&(s->propmat[icol]), &(s->propmat2[icol]), 
			  &(s->propmat2[icol]) );
	  su3vec_copy( (su3_vector *)(tmp1->e[icol]),
		       &(s->propmat[icol]) );
	}
      }
      
      
      
      /* Next gather the propagators in preparation for calculating   */
      /* shifted propagators Dq and Do                                */
      
      FORALLUPDIRBUT(TUP,dir) {
	  /* Start bringing "q" = propmat from forward sites    */
	  
	  mtag[dir] = start_gather_site(F_OFFSET(propmat[0]), 
		   sizeof(su3_matrix), dir, EVENANDODD, gen_pt[dir]);
	  
	  /* Start bringing "q" from backward neighbors       */
	  
	  mtag[dir+4] = start_gather_site(F_OFFSET(propmat[0]), 
		   sizeof(su3_matrix), OPP_DIR(dir), EVENANDODD,
		   gen_pt[dir+4]);
	  wait_gather(mtag[dir]);
	  wait_gather(mtag[dir+4]);
	  
	  /* Start bringing "o" = propmat2 from forward sites   */
	  
	  mtag[8+dir] = start_gather_site(F_OFFSET(propmat2[0]), 
		   sizeof(su3_matrix), dir, EVENANDODD, gen_pt[8+dir]);
	  
	      /* Start bringing "o" from backward neighbors       */
	  
	  mtag[8+dir+4] = start_gather_site(F_OFFSET(propmat2[0]), 
		    sizeof(su3_matrix), OPP_DIR(dir), EVENANDODD,
		    gen_pt[8+dir+4]);
	  wait_gather(mtag[8+dir]);
	  wait_gather(mtag[8+dir+4]);
	  
	}
      
      
      /* Calculate and dump delta propagator */
      for(t=0; t<nt; t++) {
	  /* define the time value offset t from t_source */
	  t_off = (t+t_source)%nt;
	  
	  /* Calculate contribution for each permutation of source color */
	  proptmp = 0;
	  delta_prop (0,1,2, 1, t_off, &proptmp, tempmat1);
	  delta_prop (1,2,0, 1, t_off, &proptmp, tempmat1);
	  delta_prop (2,0,1, 1, t_off, &proptmp, tempmat1);
	  delta_prop (1,0,2,-1, t_off, &proptmp, tempmat1);
	  delta_prop (0,2,1,-1, t_off, &proptmp, tempmat1);
	  delta_prop (2,1,0,-1, t_off, &proptmp, tempmat1);
	  
	  if( (((t+t_source)/nt-t_source/nt)%2) == 0 )
	    delprop[t] += proptmp;
	  else
	    delprop[t] -= proptmp;

	  /*if( (((t+t_source)/nt-t_source/nt)%2) == 1 ) delprop[t] *= -1;*/
	} /* nt-loop */
      
      /* Calculate the "q" source nucleon as a check */
      
      /* Calculate and dump nucleon check propagator */
      for(t=0; t<nt; t++) {
	  /* define the time value offset t from t_source */
	  t_off = (t+t_source)%nt;
	  
	  proptmp = 0;
	  for(x=0;x<nx;x+=2)for(y=0;y<ny;y+=2)for(z=0;z<nz;z+=2) {
	      if( node_number(x,y,z,t_off) != mynode() )continue;
	      i=node_index(x,y,z,t_off);
	      /* The q propagator is in propmat */
	      cc = det_su3( (su3_matrix *)(lattice[i].propmat) );
	      proptmp += cc.real;
	    }

	  if( (((t+t_source)/nt-t_source/nt)%2) == 0 )
	    ckbarprop[t] += proptmp;
	  else
	    ckbarprop[t] -= proptmp;

	  /*if( (((t+t_source)/nt-t_source/nt)%2) == 1 )ckbarprop[t]*= -1.0;*/
	  /* change sign because antiperiodic b.c.  sink point
	     should really be in a copy of the lattice */
	} /* nt-loop */
      
      /* Calculate nonlocal meson propagators and local check */
      for(t=0; t<nt; t++) {
	  /* Calculate two nonlocal pion propagators */
	  /* These are pi_1 and pi_1 tilde of Gupta et al */
	  /* Also calculate two local propagators as a check */
	  
	  /* define the time value offset t from t_source */
	  t_off = (t+t_source)%nt;
	  
	  nl_meson_prop(t_off,&nlpiprop[t],&nlpi2prop[t],&ckpiprop[t],
	     &ckpi2prop[t],tempmat1,tempmat2);
	  
	} /* nt-loop */
      
      /* Clean up gathers */
      FORALLUPDIRBUT(TUP,dir) {
	  cleanup_gather(mtag[dir]);
	  cleanup_gather(mtag[dir+4]);
	  cleanup_gather(mtag[8+dir]);
	  cleanup_gather(mtag[8+dir+4]);
	}
      
    } /* end loop on t_source */

  g_vecdoublesum( piprop, nt );
  g_vecdoublesum( pi2prop, nt );
  g_vecdoublesum( rhoprop, nt );
  g_vecdoublesum( rho2prop, nt );
  g_vecdoublesum( barprop, nt );
  g_vecdoublesum( nlpiprop, nt );
  g_vecdoublesum( nlpi2prop, nt );
  g_vecdoublesum( ckpiprop, nt );
  g_vecdoublesum( ckpi2prop, nt );
  g_vecdoublesum( delprop, nt );
  g_vecdoublesum( ckbarprop, nt );

  /* dump propagators */
  if( this_node==0 ){
    printf("STARTPROP\n");
    printf("MASSES:   %e  %e\n",vmass,vmass);
    printf("SOURCE: EVEN_WALL\n");
    printf("SINKS: PION_PS PION_SC\n");
    for(t=0;t<nt;t++) printf("%d %e 0.0 %e 0.0\n",t,
      piprop[t]/n_sources, pi2prop[t]/n_sources);
    printf("ENDPROP\n");

    printf("STARTPROP\n");
    printf("MASSES:  %e   %e  %e\n",vmass,vmass,vmass);
    printf("SOURCE: EVEN_WALL\n");
    printf("SINKS: NUCLEON\n");
    for(t=0;t<nt;t++) printf("%d %e 0.0\n",t, barprop[t]/n_sources);
    printf("ENDPROP\n");

    printf("STARTPROP\n");
    printf("MASSES:   %e  %e\n",vmass,vmass);
    printf("SOURCE: EVENANDODD_WALL\n");
    printf("SINKS: PION_PS PION_SC PION_i5 PION_ij\n");
    for(t=0;t<nt;t++) printf("%d %e 0.0 %e 0.0 %e 0.0 %e 0.0\n",t,
      ckpiprop[t]/n_sources, ckpi2prop[t]/n_sources,
      nlpiprop[t]/n_sources, nlpi2prop[t]/n_sources);
    printf("ENDPROP\n");

    printf("STARTPROP\n");
    printf("MASSES:  %e   %e  %e\n",vmass,vmass,vmass);
    printf("SOURCE: EVENANDODD_WALL\n");
    printf("SINKS: NUCLEON DELTA\n");
    for(t=0;t<nt;t++) printf("%d %e 0.0 %e 0.0\n",t, 
	ckbarprop[t]/n_sources, delprop[t]/n_sources );
    printf("ENDPROP\n");

  }

  free( piprop  ); free( pi2prop  ); free( rhoprop  ); free( rho2prop  );
  free( nlpiprop  ); free( nlpi2prop  );
  free( ckpiprop  ); free( ckpi2prop  );
  free( barprop  ); free( delprop  ); free( ckbarprop  );
  
  return(cgn);
} /* nl_spectrum */
