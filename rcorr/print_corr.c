/***************** print_oorr.c *****************************************/

/* Extrapolate to inifinite block size and print results */

/* MIMD version 7 */

/* 04/05/15 C. DeTar */

#include "rcorr_includes.h"

#define FORALLSITESLIN(x,y,z,t,a,b,c,d) for(x=0;x<a;x++)for(y=0;y<b;y++)for(z=0;z<c;z++)for(t=0;t<d;t++)

static FILE *
open_corr_file(void){
  FILE *fp;
  int jflav;

  if(this_node != 0)
    return NULL;

  fp = fopen(param.corrfile,"a");
  if(fp == NULL){
    printf("open_corr_file: ERROR. Can't open %s\n",
	   param.corrfile);
    return NULL;
  }

  /* Write YAML header */
  fprintf(fp,"---\n");
  fprintf(fp,"JobID:                        %s\n",param.job_id);
  fprintf(fp,"date:                         \"%s UTC\"\n",utc_date_time);
  fprintf(fp,"lattice_size:                 %d,%d,%d,%d\n", nx, ny, nz, nt);
  fprintf(fp,"random_sources_sloppy:        %d\n", param.nrand_sloppy);
  fprintf(fp,"random_sources_diff:          %d\n", param.nrand_diff);
  fprintf(fp,"files_sloppy:                 [\n");
  for(jflav = 0; jflav < param.nflav; jflav++){
    fprintf(fp,"{ charge: %g , mass: %g , file: %s }", 
	    param.charges[jflav], param.mass[jflav], param.fname_sloppy[jflav]);
    if(jflav < param.nflav-1)fprintf(fp,",");
    fprintf(fp,"\n");
  }
  fprintf(fp,"]\n");

  fprintf(fp,"files_diff:                 [\n");
  for(jflav = 0; jflav < param.nflav; jflav++){
    fprintf(fp,"{ charge: %g , mass: %g , file: %s }", 
	    param.charges[jflav], param.mass[jflav], param.fname_diff[jflav]);
    if(jflav < param.nflav-1)fprintf(fp,",");
    fprintf(fp,"\n");
  }
  fprintf(fp,"]\n");

  fprintf(fp,"...\n");
  return fp;
}

static void
close_corr_file(FILE *fp){
  if(this_node != 0)return;
  if(fp != NULL)fclose(fp);
}

static int 
wrap(int x, int nx){
  return (x < nx/2) ? x : x - nx;
}

static double 
radius(int x, int y, int z, int t){
  int k;
  double r;

  k = wrap(x,nx);
  r = k*k;
  k = wrap(y,ny);
  r += k*k;
  k = wrap(z,nz);
  r += k*k;
  k = wrap(t,nt);
  r += k*k;

  return sqrt(r);
}

/* The bin selection is supposed to make it possible to
   distinguish precisely the radii smaller than about 3
   and then give gradually increasing bin widths with
   increasing r */

/* Binning algorithm. Convert radius to bin count */
static int 
r2bin(double r){
  int k;
  k = 75.*log((r + 15.1)/15.);
  if(k >= MAXBIN)k = MAXBIN - 1;
  return k;
}

#if 0
/* Binning algorithm.  Convert bin count to radius */
static double 
bin2r(int k){
  return 15.*exp(k/75.)-15.1;
}
#endif

void
print_result(Real *q[], Real *q2[], int nblock, int block_size[]){
  double corrvsr[MAXBIN][nblock];
  double corrvsr2[MAXBIN][nblock];
  int nvsr[MAXBIN];
  double rvsr[MAXBIN];
  double r, myq, dmyq, qtot[nblock];
  double myqb[nblock], myqb2[nblock];
  double bsinv[nblock], sd[nblock];
  double m, dm, chisq;
  int mult;
  int totmult[nblock];
  int i, ib;
  int x,y,z,t;
  FILE *fp;

  fp = open_corr_file();

  for(ib = 0; ib < nblock; ib++){
    bsinv[ib] = 1./((double) block_size[ib]);
  }

  /* Unbinned output for r <= RMAX */

  for(ib = 0; ib < nblock; ib++){
    totmult[ib] = 0;
    qtot[ib] = 0.;
  }
  
  for(x = 0; x <= nx/2; x++)
    for(y = x; y <= ny/2; y++)
      for(z = y; z <= nz/2; z++)
	for(t = 0; t <= nt/2; t++){
	  mult = 0;
	  for(ib = 0; ib < nblock; ib++){
	    myqb[ib] = 0.;
	    myqb2[ib] = 0.;
	  }
	  r = radius(x,y,z,t);
	  if(node_number(x,y,z,t) == this_node){
	    if(r <= RMAX){
	      /* Compute multiplicity */
	      mult = 6*16;
	      if(x == 0 || x == nx/2){mult /= 2;}
	      if(y == 0 || y == ny/2){mult /= 2;}
	      if(z == 0 || z == nz/2){mult /= 2;}
	      if(t == 0 || t == nt/2){mult /= 2;}
	      if((x == y) && (x == z)){mult /= 6;}
	      else if(x == y){mult /= 2;}
	      else if(x == z){mult /= 2;}
	      else if(y == z){mult /= 2;}
	    }
	    for(ib = 0; ib < nblock; ib++){
	      myqb[ib] = q[ib][node_index(x,y,z,t)];
	      myqb2[ib] = q2[ib][node_index(x,y,z,t)];
	    }
	  }

	  /* Collect values on all nodes by broadcasting the values 
	     at each node to all other nodes (only the values at node 0 
	     needed for printing them) */
	  g_intsum(&mult);
	  g_vecdoublesum(myqb, nblock);
	  g_vecdoublesum(myqb2, nblock);

	  if(this_node == 0)
	    //	    if(r <= RMAX && x%2==0 && y%2==0 && z%2==0 && t%2==0){
	    if(r <= RMAX && (x+y+z+t)%2 == 0){
	      /* Get stdev of mean over symmetry-related sites */
	      /* Assumes statistical independence  - not assured */
	      for(ib = 0; ib < nblock; ib++)
		sd[ib] = sqrt(myqb2[ib]/mult);
	      chisq = linearlsq(&m, &dm, &myq, &dmyq, bsinv, myqb, sd, nblock);

	      /* DEBUG */
	      for(ib = 0; ib < nblock; ib++){
		fprintf(fp, "%.8e %.8e ", myqb[ib], sd[ib]);
	      }
	      fprintf(fp, "\n");
	      
	      fprintf(fp, "%5d %7.3f %15.8e %15.8e %7.3f %2d %2d %2d %2d \n",
		      mult,r,myq,dmyq,chisq,x,y,z,t);
	      totmult[ib] += mult;
	      qtot[ib] += mult*myq;
	    }
	}

  /* Binned output for r > RMAX */
  /* We have output only for stride 2 */
  
  for(i = 0; i < MAXBIN; i++){
    nvsr[i] = 0;
    for(ib = 0; ib < nblock; ib++){
      corrvsr[i][ib] = 0.;
      corrvsr2[i][ib] = 0.;
    }
    rvsr[i] = 0.;
  }

  for(x = 0; x < nx; x++)
    for(y = 0; y < ny; y++)
      for(z = 0; z < nz; z++)
	for(t = 0; t < nt; t++){
	  if(node_number(x,y,z,t) == this_node){
	    r = radius(x,y,z,t);
	    if(r > RMAX){
	      i = r2bin(r);
	      nvsr[i]++;
	      rvsr[i] += r;
	      for(ib = 0; ib < nblock; ib++){
		corrvsr[i][ib] += q[ib][node_index(x,y,z,t)];
		corrvsr2[i][ib] += q2[ib][node_index(x,y,z,t)];
	      }
	    }
	  }
	}  

  /* Collect values on all nodes by broadcasting the values
     at each node to all other nodes (only node 0 needed) */
  g_vecdoublesum(rvsr, MAXBIN);
  for(i = 0; i < MAXBIN; i++){
    g_intsum(&nvsr[i]);
    g_vecdoublesum(&corrvsr[i][0], nblock);
    g_vecdoublesum(&corrvsr2[i][0], nblock);
    if(nvsr[i] != 0)
      rvsr[i] /= nvsr[i];
  }

  if(this_node == 0)
    for(i = 0; i < MAXBIN; i++){
      if(nvsr[i] != 0 && rvsr[i] > RMAX){
	for(ib = 0; ib < nblock; ib++){
	  corrvsr[i][ib] /= nvsr[i];
	  /* Get stdev of mean over binned sites */
	  /* Assumes statistical independence  - not assured */
	  sd[ib] = sqrt(corrvsr2[i][ib])/nvsr[i];
	}
	chisq = linearlsq(&m, &dm, &myq, &dmyq, bsinv, &corrvsr[i][0], sd, nblock);

	/* DEBUG */
	for(ib = 0; ib < nblock; ib++){
	  fprintf(fp, "%.8e %.8e ", corrvsr[i][ib], sd[ib]);
	}
	fprintf(fp, "\n");
	
	fprintf(fp, "%5d %7.3f %15.8e %15.8e %7.3f  bin\n",
		nvsr[i], rvsr[i], myq, dmyq, chisq);
	totmult[ib] += nvsr[i];
	for(ib = 0; ib < nblock; ib++){
	  qtot[ib] += nvsr[i]*corrvsr[i][ib];
	}
      }
    }

  if(this_node == 0){
    close_corr_file(fp);
    for(ib = 0; ib < nblock; ib++)
      node0_printf("qtot2[%d] = %e  mult = %d\n", ib, qtot[ib], totmult[ib]);
  }
  
} /* print_corr.c */

void
print_result_time(Real *q[], Real *q2[], int nblock, int block_size[]){
  double myq, dmyq;
  double myqb[nblock], myqb2[nblock];
  double bsinv[nblock];
  double m, dm, chisq;
  int t, ib;
  FILE *fp;

  fp = open_corr_file();

  for(ib = 0; ib < nblock; ib++)
    bsinv[ib] = 1./((double) block_size[ib]);

  for(t=0;t<nt;t++){
    for(ib=0;ib<nblock;ib++){
      myqb[ib]=q[ib][t];
      myqb2[ib]=sqrt(q2[ib][t]);
    }
    if(this_node == 0){
      chisq = linearlsq(&m, &dm, &myq, &dmyq, bsinv, myqb, myqb2, nblock);

      /* DEBUG */
      for(ib = 0; ib < nblock; ib++){
	fprintf(fp, "%.8e %.8e ", myqb[ib], myqb2[ib]);
      }
      fprintf(fp, "\n");
      
      fprintf(fp, "%2d %15.8e %15.8e %7.3f \n",t,myq,dmyq,chisq);
    }
  }

  
  if(this_node == 0)
    close_corr_file(fp);
    
} /* print_corr.c */

void
print_result_time2(Real *q[], Real *q2[], int nblock, int block_size[]){
  double corrvsr[MAXBIN][nblock];
  double corrvsr2[MAXBIN][nblock];
  int nvsr[MAXBIN];
  double rvsr[MAXBIN];
  double r, myq, dmyq, qtot[nblock];
  double myqb[nblock], myqb2[nblock];
  double bsinv[nblock], sd[nblock];
  double m, dm, chisq;
  double *c_t = (double *) malloc(sizeof(double)*nt);
  double *c2_t= (double *) malloc(sizeof(double)*nt);
  int mult;
  int totmult[nblock];
  int i, ib;
  int x,y,z,t;
  FILE *fp;

  fp = open_corr_file();

  for(ib = 0; ib < nblock; ib++){
    bsinv[ib] = 1./((double) block_size[ib]);
    qtot[ib] = 0.;
  }
  
  for(t=0;t<nt;t++){
    c_t[t] = 0;
    c2_t[t]=0;
  }
  
  FORALLSITESLIN(x,y,z,t,nx,ny,nz,nt){
    for(ib = 0; ib < nblock; ib++){
      myqb[ib] = 0.;
      myqb2[ib] = 0.;
    }
    if(node_number(x,y,z,t) == this_node){
      for(ib = 0; ib < nblock; ib++){
	myqb[ib] = q[ib][node_index(x,y,z,t)];
	myqb2[ib] = q2[ib][node_index(x,y,z,t)];
      }
    }
    
    /* Collect values on all nodes by broadcasting the values 
       at each node to all other nodes (only the values at node 0 
       needed for printing them) */
    g_vecdoublesum(myqb, nblock);
    g_vecdoublesum(myqb2, nblock);

    if(this_node==0){
      for(ib = 0; ib < nblock; ib++)
	sd[ib] = sqrt(myqb2[ib]);
      chisq = linearlsq(&m, &dm, &myq, &dmyq, bsinv, myqb, sd, nblock);
      node0_printf("%e %e %d %d %d %d\n",myq, dmyq,x,y,z,t);
      c_t[t]+=myq/3;
      c2_t[t]+=dmyq*dmyq/9;
    }
  }
  
  if(this_node == 0){
    for(t=0;t<nt;t++)
      fprintf(fp, "%2d %15.8e %15.8e %7.3f \n",t,c_t[t],sqrt(c2_t[t]),0);
  }

  if(this_node == 0){
    close_corr_file(fp);
    for(ib = 0; ib < nblock; ib++)
      node0_printf("qtot2[%d] = %e  mult = %d\n", ib, qtot[ib], totmult[ib]);
  }
  
} /* print_corr.c */
