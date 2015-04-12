/***************** print_oorr.c *****************************************/

/* Symmetrize the current-current correlator over hypercubic group transformations */

/* MIMD version 7 */

/* 04/05/15 C. DeTar */

#include "rcorr_includes.h"

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
  fprintf(fp,"random_sources:               %d\n", param.nrand);
  fprintf(fp,"files:                        [\n");
  for(jflav = 0; jflav < param.nflav; jflav++){
    fprintf(fp,"{ charge: %g , mass: %g , file: %s }", 
	    param.charges[jflav], param.mass[jflav], param.fname[jflav]);
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
print_result(Real *q, int nrand){
  double rhovsr[MAXBIN];
  int nvsr[MAXBIN];
  double rvsr[MAXBIN];
  double r, myq, qtot;
  int mult, totmult;
  int i;
  int x,y,z,t;
  FILE *fp;

  fp = open_corr_file();

  for(i = 0; i < MAXBIN; i++){
    nvsr[i] = 0;
    rhovsr[i] = 0.;
    rvsr[i] = 0.;
  }

  /* Unbinned output for r <= RMAX */
  
  totmult = 0;
  qtot = 0.;
  
  for(x = 0; x <= nx/2; x++)
    for(y = x; y <= ny/2; y++)
      for(z = y; z <= nz/2; z++)
	for(t = 0; t <= nt/2; t++){
	  mult = 0;
	  myq = 0;
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
	    myq = q[node_index(x,y,z,t)];
	  }

	  g_intsum(&mult);
	  g_doublesum(&myq);

	  if(this_node == 0)
	    if(r <= RMAX){
	      fprintf(fp, "%5d %7.3f %15.8e %2d %2d %2d %2d \n",
		      mult,r,myq,x,y,z,t);
	      totmult += mult;
	      qtot += mult*myq;
	    }
	}

  /* Binned output for r > RMAX */
  
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
	      rhovsr[i] += q[node_index(x,y,z,t)];
	    }
	  }
	}  

  /* Normalization and output */
  
  for(i = 0; i < MAXBIN; i++){
    g_intsum(&nvsr[i]);
    g_doublesum(&rvsr[i]);
    g_doublesum(&rhovsr[i]);
    if(nvsr[i] != 0)
      rvsr[i] /= nvsr[i];
  }

  if(this_node == 0)
    for(i = 0; i < MAXBIN; i++){
      if(nvsr[i] != 0 && rvsr[i] > RMAX){
	rhovsr[i] /= nvsr[i];
	fprintf(fp, "%5d %7.3f %15.8e   bin\n",
		nvsr[i], rvsr[i], rhovsr[i]);
	totmult += nvsr[i];
	qtot += nvsr[i]*rhovsr[i];
      }
    }

  if(this_node == 0){
    close_corr_file(fp);
    node0_printf("qtot2 = %e  mult = %d\n", qtot, totmult);
  }

} /* print_corr.c */


