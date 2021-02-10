/****** plane_avg_plc.c  -- ******************/
/* Plane averaged Polyakov loop correlator
* MIMD version 7
*
* New Strategy:
* Gather entire local lattice blocks of the fields and
* interpret their entries by adding to the corresponding
* site indices.
*
* PLEASE BE AWARE THAT THERE MIGHT BE A CONFLICT!
* Introduced modified multi-distance gathers in the scheme
* of the 3n gathers for the naik term (-> setup.c). This
* is done by replacing the 3n gathers for the naik term with
* fixed nzyxt for transporting the full local lattice.
*
* Evaluate in different spatial directions, to check rotational
* invariance.  Gauge fix to Coulomb gauge to do spatial segments.
*
*/

/*
*
* Iteration 1)
*
* Modified these things in "setup.c":
* 1) Get the size of the local lattices
* 2) Make gathers for displacements in local lattice
*
* These things are done in this file (plane_avg_plc.c):
* 3) Compute the Polyakov lines in buffer wline
*    on even sites in first two time slices
*
* 4) Trace the Polyakov lines and compute plane averages
*
* 5) Contract the Polyakov loops
*
* 6) Cleanup
*/

#include "../generic/generic_includes.h"	/* definitions files and prototypes */
#include "./hvy_qpot_includes.h"
#include "../include/openmp_defs.h"
#include <assert.h>

#ifdef PLANE_AVERAGED_PLC

#ifdef OMP
#define ATOMIC \
_Pragma( STR(omp atomic) ) 
#else
#define ATOMIC
#endif

void plane_averaged_plc( su3_matrix *links ) {
    register int i,j,mu;
    int disp[4]={0,0,0,0};
    int ngather=0;
    int ns[3]={nc[XUP],nc[YUP],nc[ZUP]};
    int maxns,minns, maxnc;
    site *s=NULL;

    msg_tag *mtag0=NULL;
    su3_matrix *owline=NULL,buffer;
    complex ctr1, cprod;
    double_complex **plp=NULL;
    double_complex **plc=NULL;
    double_complex *plf=NULL;
    double_complex ppmu={0,0};
    double_complex ppa={0,0};
    double_complex pp2={0,0};
    double_complex ppf={0,0};
    double iv2[3]={1./(nc[YUP]*nc[ZUP]),1./(nc[XUP]*nc[ZUP]),1./(nc[XUP]*nc[YUP])},inz[3]={1./nc[XUP],1./nc[YUP],1./nc[ZUP]};
    char myname[] = "plane_avg_plc";

    /*******************************************************
     *
     * Start of header:
     *
     * Text output concerning the geometry and its consistency.
     * In particular, the correlation directions 'xc[mu]' are
     * permuted from the original lattice directions 'mu'
     * to permit static quark correlations with arbitrary
     * correlation time direction and arbitrary anisotropic
     * direction.
     *
     ******************************************************/

    node0_printf("%s with %d times smeared links\n",myname,tot_smear);
#ifdef DEBUG
    node0_printf(" [min,max](CT)[%d] = [%d,%d], \
                \n max(C0)[%x] = %d, max(C1)[%x] = %d, max(C2)[%x] = %d\n",
                 cor_dir, min_ct,maxc[TUP], xc[XUP],maxc[XUP] ,xc[YUP],maxc[YUP] ,xc[ZUP],maxc[ZUP] );
#endif
    assert ( maxc[XUP] <= nc[XUP]/2 && maxc[YUP] <= nc[YUP]/2 && maxc[ZUP] <= nc[ZUP]/2 );
    assert ( min_ct <= maxc[TUP] && maxc[TUP] <= nc[TUP] );
    assert ( maxc[TUP] == nc[TUP] );

    maxnc=(maxc[XUP]>maxc[ZUP]?maxc[XUP]>maxc[YUP]?maxc[XUP]:maxc[YUP]:maxc[ZUP]);
    maxns=(nc[XUP]>nc[ZUP]?nc[XUP]>nc[YUP]?nc[XUP]:nc[YUP]:nc[ZUP]);
    minns=(nc[XUP]<nc[ZUP]?nc[XUP]<nc[YUP]?nc[XUP]:nc[YUP]:nc[ZUP]);
//    node0_printf("####################################\n\n");

    /****************************************************************
     *
     * End of header
     * Setup of buffers
     *
     ***************************************************************/

    owline = hqp_alloc_su3mat_buffer(1);
    /* averaged correlations of PL */
    plc = (double_complex **)malloc(3*sizeof(double_complex*));
    if(plc == NULL){
      printf("%s(%d): Can't malloc temporary\n",myname,this_node);
      terminate(1);
    }
    for ( mu=XUP; mu<=ZUP; mu++) {
      plc[ mu] = (double_complex *)malloc(maxns*sizeof(double_complex));
      memset (plc[ mu],0,maxns*sizeof(double_complex));
    }
    /* directed axial distribution of PLP */
    plp = (double_complex **)malloc(3*sizeof(double_complex*));
    if(plp == NULL){
      printf("%s(%d): Can't malloc temporary\n",myname,this_node);
      terminate(1);
    }
    for ( mu=XUP; mu<=ZUP; mu++) {
      plp[ mu] = (double_complex *)malloc(maxns*sizeof(double_complex));
      memset (plp[ mu],0,maxns*sizeof(double_complex));
    }
    FORALLFIELDSITES_OMP(i, default(shared) ){
      su3mat_copy( links+4*i+xc[TUP], owline+i );
    } END_LOOP_OMP;
    /* distribution of PLP fluctutations */
    plf = (double_complex *)malloc((volume)*sizeof(double_complex));
    if(plf == NULL){
      printf("%s(%d): Can't malloc temporary\n",myname,this_node);
      terminate(1);
    }
    memset (plf,0,(volume)*sizeof(double_complex));

    /****************************************************************
     *
     * End of setup
     * build Wilson lines and plane-summed Polyakov loop distribution
     *
     ***************************************************************/

    /* Start the Wilson lines on even sites */
    mtag0 = start_gather_field( owline, sizeof(su3_matrix), xc[TUP], EVEN, gen_pt[0] );
    wait_gather(mtag0);
    ngather++;
    FOREVENFIELDSITES_OMP(i, private(buffer) ){
      mult_su3_nn( links+4*i+xc[TUP], (su3_matrix *)gen_pt[0][i], &buffer );
      su3mat_copy( &buffer, owline+i );
    } END_LOOP_OMP;
    cleanup_gather( mtag0 );
    ngather++;

    /* Gather the Wilson line fragments */
    for(j=2;j<nc[TUP];j+=2){
      disp[xc[TUP]] = j;     /* distance from which to gather */
      mtag0 = start_general_gather_field( owline, sizeof(su3_matrix),disp, EVEN, gen_pt[0] );
      wait_general_gather(mtag0);
      ngather++;
      FORSOMEPARITY_OMP(i,s,EVEN,private(j,disp,buffer) ) {
        if ( site_coord(s,xc[TUP]) > 1 ) continue;
        mult_su3_nn( owline+i, (su3_matrix *)gen_pt[0][i], &buffer );
        su3mat_copy( &buffer, owline+i );
      } END_LOOP_OMP;
      cleanup_general_gather(mtag0);
    }

    /* Trace the Polyakov loop and sum up the xy,yz and zx slices */
    FORSOMEPARITY_OMP(i,s, EVEN, private(ctr1) ) {
      if ( site_coord(s,xc[TUP]) <= 1 ) {
//printf("i %d xyzt %d %d %d %d\n",i,site_coord(s,xc[XUP]),site_coord(s,xc[YUP]),site_coord(s,xc[ZUP]),site_coord(s,xc[TUP]));
        ctr1 = trace_su3( owline+i );
        plf[i].real = (double)(ctr1.real);
        plf[i].imag = (double)(ctr1.imag);
        ATOMIC
        plp[XUP][site_coord(s,xc[XUP])].real += plf[i].real;
        ATOMIC
        plp[XUP][site_coord(s,xc[XUP])].imag += plf[i].imag;
        ATOMIC
        plp[YUP][site_coord(s,xc[YUP])].real += plf[i].real;
        ATOMIC
        plp[YUP][site_coord(s,xc[YUP])].imag += plf[i].imag;
        ATOMIC
        plp[ZUP][site_coord(s,xc[ZUP])].real += plf[i].real;
        ATOMIC
        plp[ZUP][site_coord(s,xc[ZUP])].imag += plf[i].imag;
      }
    } END_LOOP_OMP;

    /****************************************************************
     *
     * build correlators and averages
     *
     ***************************************************************/

    for ( mu=XUP; mu<=ZUP; mu++) {
      ppmu.real = 0;
      ppmu.imag = 0;
      g_vecdcomplexsum( plp[ mu],maxns );
      for (j=0;j<ns[ mu];j++) CMULREAL(plp[ mu][j],iv2[ mu],plp[ mu][j]);

      for (i=0;i<ns[ mu];i++) {
        ppmu.real += plp[ mu][i].real*inz[ mu]/3.;
        ppmu.imag += plp[ mu][i].imag*inz[ mu]/3.;
        pp2.real += ( plp[ mu][i].real*plp[ mu][i].real + plp[ mu][i].imag*plp[ mu][i].imag )*inz[ mu]/3.;
      }
      ppa.real += ppmu.real;
      ppa.imag += ppmu.imag;
    }
    for ( mu=XUP; mu<=ZUP; mu++) {
      /* first take out the average */
      for (j=0;j<ns[ mu];j++) {
        plp[ mu][j].real -= ppa.real;
        plp[ mu][j].imag -= ppa.imag;
      }
      /* then look at correlations */
      for (j=0;j<=maxc[mu];j++) {
        for (i=0;i<ns[ mu];i++) {
          CMUL_J(plp[ mu][i],plp[ mu][(i+j+ns[ mu])%(ns[ mu])],cprod);
          CMULREAL(cprod,inz[ mu]/3.,cprod);
          plc[mu][j].real += (double)cprod.real;
          plc[mu][j].imag += (double)cprod.imag;
        }
      }
    }

    g_vecdcomplexsum( plf,volume );
    double ppfr=ppf.real, ppfi=ppf.imag;
    FORSOMEPARITY_OMP(i,s,EVEN, default(shared) reduction(+:ppfr,ppfi) ) {
      if ( site_coord(s,xc[TUP]) > 1 ) continue;

      plf[i].real -= (double)(ppa.real);
      plf[i].imag -= (double)(ppa.imag);

      ppfr += plf[i].real*plf[i].real/(nc[XUP]*nc[YUP]*nc[ZUP]);
      ppfi += plf[i].imag*plf[i].imag/(nc[XUP]*nc[YUP]*nc[ZUP]);
    } END_LOOP_OMP;
    ppf.real=ppfr; ppf.imag=ppfi;

#ifdef SMEARING
          node0_printf("P_LOOP_%d:\t%e\t%e\n", tot_smear, ppa.real, ppa.imag );
#else
          node0_printf("P_LOOP:\t%e\t%e\n", ppa.real, ppa.imag );
#endif
    for ( mu=XUP; mu<=ZUP; mu++) {
      if ( maxc[mu]>0 ) {
#ifdef SMEARING
        for (j=0;j<=maxc[mu];j++) node0_printf("PAP_CORR_%d: [%d] %d %d %d %d \t%.6e %.6e\n",tot_smear, mu, 0, 0, j, nc[TUP], plc[mu][j].real,plc[mu][j].imag);
#else
        for (j=0;j<=maxc[mu];j++) node0_printf("PAP_CORR: [%d] %d %d %d %d \t%.6e %.6e\n", mu, 0, 0, j, nc[TUP], plc[mu][j].real,plc[mu][j].imag );
#endif
      }
    }

    /* Cleanup */
    free(owline); owline=NULL;
    for ( mu=XUP; mu<=ZUP; mu++) {
      free(plp[ mu]); plp[ mu]=NULL;
    }
    free(plp); plp=NULL;
    for ( mu=XUP; mu<=ZUP; mu++) {
      free(plc[ mu]); plc[ mu]=NULL;
    }
    free(plc); plc=NULL;
#ifdef DEBUG
  node0_printf("NUMBER OF GATHERS = %d\n",ngather);
#endif
}

#endif
