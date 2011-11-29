/******** wf_pt_t.c *************/
/* MIMD version 6 */
/* NOT MAINTAINED.  TEST BEFORE USE! */

/* Wave functions along principal axis with parallel transport.
   Slight modification of spectrum.c */
/* This version does two wall sources */
/* The wall sources are on slices 0 and nt/2.  */
/* For measuring propagation in the t direction */

/* Spectrum for Kogut-Susskind pointlike hadrons, wall source */
#include "ks_dyn_includes.h"

#define TMAX 32 /* Maximum nt for which this works */
int wf_pt_t() /* return the C.G. iteration number */
{
  complex mesprop[2][2][2][TMAX];
  complex barprop[TMAX];
  Real mass_x2;
  register complex cc;
  Real finalrsq;
  register int i,j,k,x,y,z,t,icol,cgn;
  int isite,dz,t_source,t_off;
  site *s;
  msg_tag *tag;
  su3_vector tmat[3];

  
  if( nt > TMAX && this_node==0){printf("nt TOO BIG!\n"); terminate(0);}
  /* Fix TUP Coulomb gauge - gauge links only*/
  rephase( OFF );
  gaugefix(TUP,(Real)1.8,500,(Real)GAUGE_FIX_TOL);
  rephase( ON );
  
  mass_x2 = 2.0*mass;
  cgn=0;
  for(t_source=0; t_source<nt; t_source += (nt/2) )
    {
	  for(icol=0; icol<3; icol++) {
		      
	      /* initialize ttt and xxx */
	      clear_latvec( F_OFFSET(phi), EVENANDODD);
              clear_latvec( F_OFFSET(xxx), EVENANDODD);

	      /* wall source */
              for(x=0;x<nx;x+=2)for(y=0;y<ny;y+=2)for(z=0;z<nz;z+=2) {
                     if( node_number(x,y,z,t_source) != mynode())continue;
                     i=node_index(x,y,z,t_source);
                     lattice[i].phi.c[icol].real = -1.0;
		}

	      /* do a C.G. (source in phi, result in xxx) */
	      load_ferm_links(&fn_links);
	      cgn += ks_congrad(F_OFFSET(phi),F_OFFSET(xxx),mass,
				niter, rsqprop, PRECISION, EVEN, &finalrsq,
				&fn_links);
	      /* Multiply by -Madjoint */
	      dslash_site( F_OFFSET(xxx), F_OFFSET(ttt), ODD, &fn_links);
	      scalar_mult_latvec( F_OFFSET(xxx), -mass_x2, F_OFFSET(ttt), EVEN);
	      
	      /* fill the hadron matrix */
	      copy_latvec( F_OFFSET(ttt), F_OFFSET(propmat[icol]), EVENANDODD);
	    } /* end loop on icol */
	  
	/* Make copy of propmat to parallel transport */
	FORALLSITES(isite,s)for(icol=0;icol<3;icol++){
	  s->propmat2[icol] = s->propmat[icol];
	}
	/* Loop over displacements */
	for(dz=0;dz<=nx/2;dz+=2){

    /* make pion and rho wave functions in mesprop[i][j][k][t] */
    /* index     sign            particle
        i,j,k
        0,0,0   +1              pi (PS)
        1,0,0   (-1)^x          rho-b1 (VT), gamma_x polarization
        0,1,0   (-1)^y          rho-b1 (VT), gamma_y polarization
        0,0,1   (-1)^z          rho-b1 (VT), gamma_z polarization
        1,1,0   (-1)^(x+y)      rho-a1 (PV), gamma_z polarization
        1,0,1   (-1)^(x+z)      rho-a1 (PV), gamma_y polarization
        0,1,1   (-1)^(y+z)      rho-a1 (PV), gamma_x polarization
        1,1,1   (-1)^(x+y+z)    pi_2, sigma (SC)
    */
	  /* clear propagators */
	  for(i=0;i<2;i++)for(j=0;j<2;j++)for(k=0;k<2;k++)for(t=0;t<nt;t++){
	    mesprop[i][j][k][t]= cmplx(0.0,0.0);
	  }
	  for(t=0;t<nt;t++)barprop[t]= cmplx(0.0,0.0);

	  /* measure the meson propagator */
	  for(t=0; t<nt; t++) {
	      /* define the time value offset t from t_source */
	      t_off = (t+t_source)%nt;
	      
	      for(x=0;x<nx;x++)for(y=0;y<ny;y++)for(z=0;z<nz;z++){
		for(icol=0;icol<3;icol++) {
		    if( node_number(x,y,z,t_off) != mynode() )continue;
		    isite=node_index(x,y,z,t_off);
		    cc = su3_dot( &lattice[isite].propmat[icol],
				 &lattice[isite].propmat2[icol] );
		    
		    {CADD( mesprop[0][0][0][t], cc, mesprop[0][0][0][t] );}
		    
		    if( x%2==0){CADD( mesprop[1][0][0][t], cc,
			mesprop[1][0][0][t] );}
		    else       {CSUB( mesprop[1][0][0][t], cc,
			mesprop[1][0][0][t] );}
		    if( y%2==0){CADD( mesprop[0][1][0][t], cc,
			mesprop[0][1][0][t] );}
		    else       {CSUB( mesprop[0][1][0][t], cc,
			mesprop[0][1][0][t] );}
		    if( z%2==0){CADD( mesprop[0][0][1][t], cc,
			mesprop[0][0][1][t] );}
		    else       {CSUB( mesprop[0][0][1][t], cc,
			mesprop[0][0][1][t] );}
		    
		    if( (x+y)%2==0){CADD(mesprop[1][1][0][t],cc,
			mesprop[1][1][0][t]);}
		    else           {CSUB(mesprop[1][1][0][t],cc,
			mesprop[1][1][0][t]);}
		    if( (y+z)%2==0){CADD(mesprop[0][1][1][t],cc,
			mesprop[0][1][1][t]);}
		    else           {CSUB(mesprop[0][1][1][t],cc,
			mesprop[0][1][1][t]);}
		    if( (x+z)%2==0){CADD(mesprop[1][0][1][t],cc,
			mesprop[1][0][1][t]);}
		    else           {CSUB(mesprop[1][0][1][t],cc,
			mesprop[1][0][1][t]);}
		    
		   if((x+y+z)%2==0){CADD(mesprop[1][1][1][t],cc,
			mesprop[1][1][1][t]);}
		   else            {CSUB(mesprop[1][1][1][t],cc,
			mesprop[1][1][1][t]);}
		    
		}
	      }
	      
	      /* measure the baryon propagator */
	      for(x=0;x<nx;x+=2)for(y=0;y<ny;y+=2)for(z=0;z<nz;z+=2) {
		  if( node_number(x,y,z,t_off) != mynode() )continue;
		  i=node_index(x,y,z,t_off);
		  for(icol=0;icol<3;icol++){
		    for(j=0;j<3;j++){
		 	if(j==icol)tmat[j]=lattice[i].propmat2[j];
			else       tmat[j]=lattice[i].propmat[j];
		    }
		    cc = det_su3( tmat );
		    CSUM(barprop[t],cc);
		  }
		}
		      
	  } /* nt-loop */

	  for(t=0; t<nt; t++) {
	      /* sum propagators */
	      for(i=0;i<2;i++)for(j=0;j<2;j++)for(k=0;k<2;k++){
	        g_complexsum( &(mesprop[i][j][k][t]) );
	      }
	      g_complexsum( &(barprop[t]) );
	  }

	  if( this_node==0) for(t=0; t<=nt/2; t++) {
	      /* dump meson propagators */
	      /* Since displacement is along z axis, average gamma_x
		and gamma_y mesons  (i <-> j) */
	      for(i=0;i<2;i++)for(j=0;j<2;j++)for(k=0;k<2;k++){
		if( i>j ) continue;
	        printf("PMF:  %d   %d %d %d\t%.7e\t%.7e\n",
		  t,i,j,dz+k,
		  mesprop[i][j][k][t].real+mesprop[i][j][k][(nt-t)%nt].real+
		  mesprop[j][i][k][t].real+mesprop[j][i][k][(nt-t)%nt].real,
		  mesprop[i][j][k][t].imag+mesprop[i][j][k][(nt-t)%nt].imag+
		  mesprop[j][i][k][t].imag+mesprop[j][i][k][(nt-t)%nt].imag);
	      }
	      
	      /* dump baryon propagators */
		  /* must get sign right.  This looks to see if we have
			wrapped around the lattice.  "t" is the distance
			from the source to the measurement, so we are
			trying to find out if t_source+t is greater than
			or equal to nt. */
	      cc = cmplx(0.0,0.0);

	      if( t+t_source <= nt ){
		CADD(cc, barprop[t], cc);
	      }
	      else{
	        CSUB(cc, barprop[t], cc);
	      }

	      /**
	      if( (nt-t)+t_source >= nt ){
		if( (nt-t)%2 == 1 ){ CADD(cc, barprop[(nt-t)%nt], cc);}
		else               { CSUB(cc, barprop[(nt-t)%nt], cc);}
	      }
	      else{
		if( (nt-t)%2 == 1 ){ CSUB(cc, barprop[(nt-t)%nt], cc);}
		else               { CADD(cc, barprop[(nt-t)%nt], cc);}
	      }
	      **/

	      printf("PBF:  %d   %d %d %d\t%.7e\t%.7e\n", t,0,0,dz,
		      cc.real,cc.imag);
	    } /* nt-loop */
	  

	  /* parallel transport propmat2 by two units */
/* Temporary: use CG vectors for temps.  Have to be careful not
to overwrite propmat2 before cleanup_gather, because we are still
using it */
	  tag=start_gather_site( F_OFFSET(propmat2[0]),3*sizeof(su3_vector),ZUP,
	    EVENANDODD,gen_pt[ZUP]);
	  wait_gather(tag);
	  FORALLSITES(i,s){
	    mult_su3_mat_vec( &(s->link[ZUP]),
	      ((su3_vector *)gen_pt[ZUP][i])+0, &(s->xxx) );
	    mult_su3_mat_vec( &(s->link[ZUP]),
	      ((su3_vector *)gen_pt[ZUP][i])+1, &(s->ttt) );
	    mult_su3_mat_vec( &(s->link[ZUP]),
	      ((su3_vector *)gen_pt[ZUP][i])+2, &(s->phi) );
	  }
	  cleanup_gather(tag);
	  FORALLSITES(i,s){
	    s->propmat2[0] = s->xxx;
	    s->propmat2[1] = s->ttt;
	    s->propmat2[2] = s->phi;
	  }
	  tag=start_gather_site( F_OFFSET(propmat2[0]),3*sizeof(su3_vector),ZUP,
	    EVENANDODD,gen_pt[ZUP]);
	  wait_gather(tag);
	  FORALLSITES(i,s){
	    mult_su3_mat_vec( &(s->link[ZUP]),
	      ((su3_vector *)gen_pt[ZUP][i])+0, &(s->xxx) );
	    mult_su3_mat_vec( &(s->link[ZUP]),
	      ((su3_vector *)gen_pt[ZUP][i])+1, &(s->ttt) );
	    mult_su3_mat_vec( &(s->link[ZUP]),
	      ((su3_vector *)gen_pt[ZUP][i])+2, &(s->phi) );
	  }
	  cleanup_gather(tag);
	  FORALLSITES(i,s){
	    s->propmat2[0] = s->xxx;
	    s->propmat2[1] = s->ttt;
	    s->propmat2[2] = s->phi;
	  }

	}	/* end loop on displacements */
  
    } /* end loop on t_source */
  return(cgn);
} /* wf_pt_t */

