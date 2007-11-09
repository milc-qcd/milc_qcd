/*************************** spectrum_fzw.c *********************************/
/* MIMD version 6                                                           */
/* This version uses gathers to get the neighbors                           */
/*                                                                          */
/* Version 4.0  08/06/05 by Ziwen Fu     For his Ph.D Thesis                */
/* This version save the pbp to file also                                   */
/*                                                                          */
/* ------------------------------------------------------------------------ */
/* Warning: 1: This version DOES NOT fix the gauge -- you should do that    */
/*             before calling it.                                           */
/*          2: This is a version for temporal g_link.                       */
/*                                                                          */
/* ------------------------------------------------------------------------ */
/*  In this version. T=0, It is real 0                                      */
/*                                                                          */
/* ------------------------------------------------------------------------ */
/* Purpose: To calculate the connected part and its                         *
**          standard error of the corralator by direct method               */
/****************************************************************************/

#include "generic_ks_includes.h"
#include "../include/dslash_ks_redefine.h"


/* ////////////////////////////////////////////////////////////////////// */
/*                           MAIN PROCEDURE                               */
/* ---------------------------------------------------------------------- */
int spectrum_fzw( Real vmass, field_offset temp1, field_offset temp2,
		  ferm_links_t *fn){ 
  
  double *prop5,   *Prop5; /* Goldstone */
  double *prop000, *prop111, *prop1_11, *prop11_1, *prop1_1_1;
  double *prop200, *prop100, *prop110,  *prop1_10, *prop101, *prop10_1;
  double *Prop000, *Prop111, *Prop1_11, *Prop11_1, *Prop1_1_1;
  double *Prop200, *Prop100, *Prop110,  *Prop1_10, *Prop101, *Prop10_1;
  double *rho_pv_prop, *rho_vt_prop, *Rho_pv_prop, *Rho_vt_prop ;
  Real    vmass_x2, finalrsq;
  register double  tt, cc, c1;
  register int Sign,sign,icol,cgn, x,y,z,t, t_src, t_off, xs, ys, zs;
  register site *s = NULL;
  complex  CCC;
  int      i, I, N , nsrc, *SCT;
  int      mysource_start,mysource_inc,myn_sources;
    
  vmass_x2 = 2.0*vmass;
  cgn     = 0;
  c1      = (2.0*PI/(double)nx);
  N = 3;   /* 2^6 = 64 = nt */
  
  I = nt*sizeof(double);
  prop5     = (double *)malloc(I);   Prop5     = (double *)malloc(I);   
 
  prop000   = (double *)malloc(I);   prop111   = (double *)malloc(I);
  prop1_11  = (double *)malloc(I);   prop11_1  = (double *)malloc(I);
  prop1_1_1 = (double *)malloc(I);   prop100   = (double *)malloc(I);
  prop110   = (double *)malloc(I);   prop1_10  = (double *)malloc(I);
  prop101   = (double *)malloc(I);   prop10_1  = (double *)malloc(I);
  prop200   = (double *)malloc(I);
 
  Prop000   = (double *)malloc(I);   Prop111   = (double *)malloc(I);
  Prop1_11  = (double *)malloc(I);   Prop11_1  = (double *)malloc(I);
  Prop1_1_1 = (double *)malloc(I);   Prop100   = (double *)malloc(I);
  Prop110   = (double *)malloc(I);   Prop1_10  = (double *)malloc(I);
  Prop101   = (double *)malloc(I);   Prop10_1  = (double *)malloc(I);
  Prop200   = (double *)malloc(I);
  
  SCT  = (int *)malloc( nt*sizeof(int) ); 
  
  rho_vt_prop = (double *)malloc(I);  /* "rho" */
  rho_pv_prop = (double *)malloc(I);  /* "rho2" */
  Rho_vt_prop = (double *)malloc(I);  /* "rho" */
  Rho_pv_prop = (double *)malloc(I);  /* "rho2" */

  for( t=0; t<nt; t++){
    prop000[t]=0.0;    prop111[t]=0.0;    prop1_1_1[t]=0.0;
    prop1_11[t]=0.0;   prop11_1[t]=0.0;   prop100[t]=0.0;
    prop110[t]=0.0;    prop1_10[t]=0.0;   prop101[t]=0.0;
    prop10_1[t]=0.0;   prop200[t]=0.0;    SCT[t] = 0;  prop5[t]=0.0;
    rho_vt_prop[t]=0;  Rho_vt_prop[t]=0;
    rho_pv_prop[t]=0;  Rho_pv_prop[t]=0;
  }

  mysource_start = 0;  mysource_inc = 2*nt;   myn_sources = nt;
  
  /* -------------------------------------------------------------------- */
  xs = 0;   ys =0;   zs =0;  nsrc=0;

  
  for(i=0; i<=N; i++){
    
    mysource_inc /= 2;
    if( mysource_inc <=0 ) exit(0);
    if( this_node  ==0 ) printf("\n\nSource_inc= %d\n", mysource_inc);
   
    for(t_src=mysource_start;  t_src<myn_sources; t_src+=mysource_inc ) {
      if(this_node==0)printf("spectrum_fzw(): source time = %d\n",t_src);

      if ( SCT[ t_src ] == 1 )   continue;
      else SCT[ t_src ] = 1;
           
      nsrc++;
      
      if(this_node==0) printf("spectrum(): source time = %d\n", t_src);
 
      for(icol=0; icol<3; icol++) {
	
      /* initialize temp1 and temp2 */
      clear_latvec( temp1, EVENANDODD);
      clear_latvec( temp2,    EVENANDODD);
      
      if( node_number(xs,ys,zs,t_src) == mynode() ) {       
	 ((su3_vector *)(F_PT(&lattice[ node_index(xs,ys,zs,t_src)],temp1)))->c[icol].real=-1.0;
      }
      
      /* do a C.G. (source in temp1, result in temp2) */
      if( (t_src + xs + ys + zs)%2 == 0 ) {     
           cgn +=  ks_congrad(temp1, temp2,vmass,
                              niter, nrestart, rsqprop, PRECISION, 
			      EVEN, &finalrsq, fn);
           dslash_site( temp2, F_OFFSET(ttt), ODD, fn);
           scalar_mult_latvec( temp2, -vmass_x2, F_OFFSET(ttt), EVEN);
      }else {      
           cgn +=  ks_congrad(temp1, temp2,vmass,
			      niter, nrestart, rsqprop, PRECISION, 
			      ODD, &finalrsq, fn);
           dslash_site( temp2, F_OFFSET(ttt), EVEN, fn);
           scalar_mult_latvec( temp2, -vmass_x2, F_OFFSET(ttt), ODD);
       }	  
	    copy_latvec( F_OFFSET(ttt), F_OFFSET(propmat[icol]), EVENANDODD);   
       }/* end of icol */

      /* --------------------------------------------------------- */
      /* --------------------------------------------------------- */
      /* measure the meson propagator */
      for(t=0; t<nt; t++) for(x=0;x<nx;x++)for(y=0;y<ny;y++)for(z=0;z<nz;z++){

	if( node_number(x,y,z,t) != mynode() )continue;
	I=node_index(x,y,z,t);
	
	   for(icol=0;icol<3;icol++) {
	   CCC = su3_dot( &lattice[I].propmat[icol],
			  &lattice[I].propmat[icol] );
	    cc = CCC.real; 

	   /* define the time value offset t from t_source */
	   t_off = ( t - t_src + nt)%nt;		 
           Sign  = ( x+ y + z+ t - xs-ys-zs-t_src)%2 ; /*Sign factor*/
		 
	   prop5[t_off] += cc;  /* Goldstone */

	   if(Sign==0)  prop000[t_off] += cc; /* Momentum (0,0,0) */
	   else         prop000[t_off] -= cc;

	   tt = cos( c1*( x-xs ) );         /* Momentum (1,0,0) */
	   if(Sign==0)  prop100[t_off] += tt*cc;
	   else         prop100[t_off] -= tt*cc;

	   tt   = cos( c1*( x-xs + y-ys) ); /* Momentum (1,1,0) */
	   if(Sign==0)  prop110[t_off]  += tt*cc;
	   else         prop110[t_off]  -= tt*cc;
 
	   tt =  cos( c1*(x-xs - y+ys) );  /* Momentum (1,-1,0) */
	   if(Sign==0)  prop1_10[t_off]  += tt*cc;
	   else         prop1_10[t_off]  -= tt*cc;

	   tt = cos( c1*(x-xs + z-zs) );  /* Momentum (1,0,1) */
	   if(Sign==0)  prop101[t_off]  += tt*cc;
	   else         prop101[t_off]  -= tt*cc;

	   tt = cos( c1*( x-xs - z+zs) );   /* Momentum (1,0,-1) */
	   if(Sign==0)  prop10_1[t_off]  += tt*cc;
	   else         prop10_1[t_off]  -= tt*cc;
		
	   tt = cos( 2.0*c1*(s->x-xs) );    /* Momentum (2,0,0) */
	   if(Sign==0)  prop200[t_off]  += tt*cc;
	   else         prop200[t_off]  -= tt*cc;

	   tt = cos( c1*( (x-xs) + (y-ys) + (z-zs) ) ); 
	   if(Sign==0)  prop111[t_off]  += tt*cc;
	   else         prop111[t_off]  -= tt*cc;

	   tt = cos( c1*( (x-xs) - (y-ys) + (z-zs) ) ); /*(1,-1,1) */
	   if(Sign==0)  prop1_11[t_off]  += tt*cc;
	   else         prop1_11[t_off]  -= tt*cc;

	   tt = cos( c1*( (x-xs) + (y-ys) - (z-zs) ) ); /*(1,1,-1)*/
	   if(Sign==0)  prop11_1[t_off]  += tt*cc;
	   else         prop11_1[t_off]  -= tt*cc;

	   tt = cos( c1*( (x-xs)-(y-ys) - (z-zs) ) ); /*(1,-1,-1)*/
	   if(Sign==0)  prop1_1_1[t_off]  += tt*cc;
	   else         prop1_1_1[t_off]  -= tt*cc;

	   if( (x+y)%2==0) rho_pv_prop[t] += cc;
	   else	                 rho_pv_prop[t] -= cc;
	   if( (y+z)%2==0) rho_pv_prop[t] += cc;
	   else	                 rho_pv_prop[t] -= cc;
	   if( (z+x)%2==0) rho_pv_prop[t] += cc;
	   else	                 rho_pv_prop[t] -= cc;
		    
	   if( (x%2)==0 ) rho_vt_prop[t] += cc;
	   else              rho_vt_prop[t] -= cc;
	   if( (y%2)==0 ) rho_vt_prop[t] += cc;
	   else              rho_vt_prop[t] -= cc;
	   if( (z%2)==0 ) rho_vt_prop[t] += cc;
	   else              rho_vt_prop[t] -= cc;
		 
	 } /* icol */
       } /* nt-loop */
     } /* end loop on t_source */
  fflush(stdout);
  /* dump the propagators */
   for( t=0; t<nt; t++){
    Prop000[t]=prop000[t];     Prop100[t]=prop100[t];   Prop200[t]=prop200[t];
    Prop111[t]=prop111[t];     Prop1_11[t]=prop1_11[t]; Prop11_1[t]=prop11_1[t];
    Prop1_1_1[t]=prop1_1_1[t]; Prop110[t]=prop110[t];   Prop1_10[t]=prop1_10[t];
    Prop101[t]=prop101[t];     Prop10_1[t]=prop10_1[t]; Prop5[t]=prop5[t];
    Rho_vt_prop[t] = rho_vt_prop[t];
    Rho_pv_prop[t] = rho_pv_prop[t];
  }

  g_vecdoublesum( Prop000,   nt );    g_vecdoublesum( Prop100,    nt );
  g_vecdoublesum( Prop200,   nt );    g_vecdoublesum( Prop111,    nt );
  g_vecdoublesum( Prop1_11,  nt );    g_vecdoublesum( Prop11_1,   nt );
  g_vecdoublesum( Prop1_1_1, nt );    g_vecdoublesum( Prop110,    nt );
  g_vecdoublesum( Prop1_10,  nt );    g_vecdoublesum( Prop101,    nt );
  g_vecdoublesum( Prop10_1,  nt );    g_vecdoublesum( Prop5,      nt );
  g_vecdoublesum( Rho_vt_prop,nt);    g_vecdoublesum( Rho_pv_prop,nt);

  /* ++++++++++++++++++++++++++ Output the data ++++++++++++++++++++++++++++ */
  if( this_node==0 ){
    printf("STARTPROP\n");
    printf("MASSES:  %e   %e\n", vmass, vmass);
    printf("SOURCE: POINT\n\n\n");
    
    printf("PROPS::   T   PION_SC000 \t  PION_SC100 \t PION_SC200 \t Goldstone\n");
    for(t=0; t<nt; t++)  printf("Prop0_%d   %d  %e \t %e \t %e \t %e\n",
	i,t,Prop000[t]/nsrc, Prop100[t]/nsrc, Prop200[t]/nsrc, Prop5[t]/nsrc);

    printf("\n\nPROPS::   T   PION_SC111 \t PION_SC1_11 \t PION_SC11_1 \t PION_SC1_1_1\n");
    for(t=0; t<nt; t++)  printf("Prop1_%d   %d   %e\t%e\t%e\t%e\n", i, t,
	Prop111[t]/nsrc, Prop1_11[t]/nsrc, Prop11_1[t]/nsrc,Prop1_1_1[t]/nsrc);

    printf("\n\nPROPS::   T   PION_SC110 \t PION_SC1_10 \t PION_SC101 \t PION_SC10_1 \n");
    for(t=0; t<nt; t++)  printf("Prop2_%d   %d   %e\t%e\t%e\t%e\n", i, t,
	Prop110[t]/nsrc, Prop1_10[t]/nsrc, Prop101[t]/nsrc, Prop10_1[t]/nsrc);

    printf("\n\nPROPS::   T   RHO_VT \t\t  RHO_PV \n");
    for(t=0; t<nt; t++)  printf("Prop3_%d   %d   %e\t%e \n", i, t,
	Rho_vt_prop[t]/nsrc, Rho_pv_prop[t]/nsrc);

    /* MILC convention */
    printf( "\nData:: >>>>>MILC convention<<<< \n");
   printf("PROPS::   T   PION_SC000 \t  PION_SC100 \t PION_SC200 \t Goldstone\n");
   for(t=0; t<nt; t++) {
      if( (t%2)==0 )  sign = 1;
      else            sign = -1;

      printf("prop0_%d   %d  %e \t %e \t %e \t %e\n",
	i, t, sign*Prop000[t]/nsrc, sign*Prop100[t]/nsrc,
	      sign*Prop200[t]/nsrc, sign*Prop5[t]/nsrc );
    }  
    printf("\n\nPROPS::   T   PION_SC111 \t PION_SC1_11 \t PION_SC11_1 \t PION_SC1_1_1\n");
    for(t=0; t<nt; t++)  {
      if( (t%2)==0 )  sign = 1;
      else            sign = -1;printf("prop1_%d   %d   %e\t%e\t%e\t%e\n", i, t,
	Prop111[t]/nsrc, Prop1_11[t]/nsrc, Prop11_1[t]/nsrc,Prop1_1_1[t]/nsrc);
    } 
    printf("\n\nPROPS::   T   PION_SC110 \t PION_SC1_10 \t PION_SC101 \t PION_SC10_1 \n");
    for(t=0; t<nt; t++)  {
      if( (t%2)==0 )  sign = 1;
      else            sign = -1;
      printf("prop2_%d   %d   %e\t%e\t%e\t%e\n", i, t,
	sign*Prop110[t]/nsrc, sign*Prop1_10[t]/nsrc,
	sign*Prop101[t]/nsrc, sign*Prop10_1[t]/nsrc);
    }    
    printf("ENDPROP\n");
  }
  fflush(stdout);  
  }
  fflush(stdout);
  return cgn;
}
