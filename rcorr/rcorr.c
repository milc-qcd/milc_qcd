/***************** rcorr.c *****************************************/

/* Measure the correlator from the given densities                   */

/* MIMD version 7 */

/* 04/03/15 C. DeTar (after topo_rcorr.c) */

#include "rcorr_includes.h"
#define FORALLSITESLIN(x,y,z,t,a,b,c,d) for(x=0;x<a;x++)for(y=0;y<b;y++)for(z=0;z<c;z++)for(t=0;t<d;t++)
#define COMPSUM(mu,mu0,muf) for(mu=mu0;mu<muf;mu++)
#if defined(COMPX)
#define MU0 0
#define MUF 1
#elif defined(COMPY)
#define MU0 1
#define MUF 2
#elif defined(COMPZ)
#define MU0 2
#define MUF 3
#elif defined(COPMT)
#define MU0 3
#define MUF 4
#else
#define MU0 0
#define MUF NMU
#endif

static void 
mulreal_c_field( complex *c, double x, int count )
{
  int i, j;
  FORALLFIELDSITES(i){
    for(j = 0; j < count; j++){
      CMULREAL(c[count*i+j], x, c[count*i+j]);
    }
  }
}

static void 
sum_c_array_field( complex *dest, complex *src, int count ){
  int i, j;

  FORALLFIELDSITES(i){
    for(j = 0; j < count; j++){
      CSUM(dest[count*i+j], src[count*i+j]);
    }
  }
}

static void 
dot_corr( complex *dest, complex *src, int count )
{
  int i, j;

  FORALLFIELDSITES(i){
    dest[i].real = 0.;
    for(j = 0; j < count; j++)
      dest[i].real += 
	(src[i*count+j].real*src[i*count+j].real + src[i*count+j].imag*src[i*count+j].imag)/volume/volume; 
    /* factor of squared inverse volume comes from normalization of Fourier decomposition of 
       current in momentum space.
       Note: Backward FFTW does not normalize the input array. 
    */
    dest[i].imag = 0.;
  }
}

static void 
sum_field_c2r( Real *dest, complex *src )
{
  int i;
  FORALLFIELDSITES(i){
    dest[i] += src[i].real;
  }
}

static void 
sum_sq_field_c2r( Real *dest, complex *src )
{
  int i;
  FORALLFIELDSITES(i){
    dest[i] += src[i].real*src[i].real;
  }
}

static void 
copy_mul_c2r( Real *dest, complex *src, double x )
{
  int i;
  FORALLFIELDSITES(i){
    dest[i] = src[i].real * x;
  }
}

static void
mean_var_r_field(Real *q, Real *q2, int n)
{
  int i;
  Real d;

  if(n > 1){
    FORALLFIELDSITES(i){
      q2[i] = q2[i]/n;
      q[i] = q[i]/n;
      d = q2[i] - q[i]*q[i];
      if(d < 0.)d = 0.;
      q2[i] = d/(n-1);
    }
  } else {
    clear_r_field(q2);
  }
}

static complex
cmul_add(complex *x, complex *y, int n){
  complex z = (complex) {0.0,0.0};
  for(int i=0;i<n;i++){
    Real a=x[i].real,b=x[i].imag,c=y[i].real,d=y[i].imag;
    z.real += a*c-b*d;
    z.imag += a*d+b*c;
  }
  return z;
}

static void
print_corr(complex *out, int x, int y, int z, int t){
  double msg = 0;
  if(this_node == node_number(x%nx,y%ny,z%nz,t%nt)) msg = (double) out[node_index(x%nx,y%ny,z%nz,t%nt)].real;
  g_doublesum(&msg);
  node0_printf("%.16e\n",msg);fflush(stdout);
}

static double
weighted_lin_extp(double * xs, Real *ys, Real *yds, int len){
  double sum_e=0,sum_x=0,sum_x2=0,sum_y=0,sum_xy=0,inv_sdm2=0;
  for(int i=0;i<len;i++){
    inv_sdm2=1.0/yds[i];
    sum_e+=inv_sdm2;
    sum_x+=xs[i]*inv_sdm2;
    sum_x2+=xs[i]*xs[i]*inv_sdm2;
    sum_y+=ys[i]*inv_sdm2;
    sum_xy+=xs[i]*ys[i]*inv_sdm2;
  }
  return (sum_x2*sum_y-sum_xy*sum_x)/(sum_e*sum_x2-sum_x*sum_x);
}

static double
avg_over_cube(Real *corr, int cx, int cy, int cz, int ct, int csize){
  double rtmp=0;
  int dx,dy,dz,dt,x,y,z,t;
  FORALLSITESLIN(dx,dy,dz,dt,csize,csize,csize,csize){
    x=dx+cx*csize;y=dy+cy*csize;z=dz+cz*csize;t=dt+ct*csize;
    if(this_node==node_number(x,y,z,t)) rtmp+=(double)corr[node_index(x,y,z,t)];
  }
  g_doublesum(&rtmp);
  return rtmp/pow(csize,4);
}

static int
is_spatial_symm_related(int x, int y, int z, int t, int rx, int ry, int rz, int rt){
  /* Assume: The representive point (rx,ry,rz) is located in the first octant. */
  int perm[3]={x,y,z},repr[3]={rx,ry,rz},dim[3]={nx,ny,nz};
  int tmp=0,flag=0;

  if((t-rt)%nt != 0) return 0;
  for(int i=0; i<6; i++){
    tmp=perm[(i+1)%3];
    perm[(i+1)%3]=perm[i%3];
    perm[i%3]=tmp;
    for(int j=0; j<3; j++) flag += (perm[j]-repr[j])%dim[j];
    if(!flag) return 1;
  }
  return 0;
}

static int 
r2(int x, int y, int z, int t){
  return x*x+y*y+z*z+t*t;
}

/******************************************************************************/
/* Input fields are have NMU values per site (one for each current component)
   There is one such field per random number group */
/* Output field has one real value per site, one such field for each
   blocking size */

void
rcorr(Real *qblock[], Real *q2block[], 
      complex *qin_sloppy[], int nrand_sloppy, 
      complex *qin_diff[], int nrand_diff,
      int nblock, int block_size[]){
  complex *qtmp;
  complex *qcorr;
  Real *q;
  int i, jrand, krand;
  
  /* Create blocked averages of the corrected result, 
     and compute the correlation within each block */
  complex *outb = create_c_field();
  complex *outf = create_c_field();
#ifdef BF
#define out outb
#else
#define out outf
#endif
  qtmp  = create_c_array_field(NMU);
  qcorr = create_c_array_field(NMU);
#ifdef CORR_VAR_DISP
  /* initialization */
  /* We compute the disc. corr. values of the given displacements
     at each lattice site.  
     To eliminate the correlation among corr. values at nearby sites, 
     we average the corr. values over a cube of the given sizes.
   */
  int cx,cy,cz,ct,ibc;
  int n_csizes = 3, csizes[3]={1,4,8}; // number of sizes & sizes of a cubic box
  int n_disp=3;
  Real *corr_avg[n_disp];
  Real *corr_sq[n_disp];
  for(int disp=0;disp<n_disp;disp++){
    corr_avg[disp]=create_r_array_field(nblock);
    corr_sq[disp]=create_r_array_field(nblock);
  }
#endif

  for(int ib = 0; ib < nblock; ib++){
    clear_r_field(qblock[ib]);
    clear_r_field(q2block[ib]);
    int bss = block_size[ib];        /* Size of one block for sloppy random sources */
    int nblocks = nrand_sloppy/bss;  /* Number of blocks */
    int bsfs[nblocks];               /* Size of each block for fine random sources */
    int nsamp = 0;                   /* Counter for block number */
    int fine_start = 0;              /* The starting index of a block of fine random solves */
    for(int ibf=0; ibf<nblocks; ibf++) bsfs[ibf] = ibf<nrand_diff%nblocks ? nrand_diff/nblocks+1: nrand_diff/nblocks;

    for(jrand = 0; jrand < nrand_sloppy; jrand += bss){
      
      /* Average qin_diff, the differnece between precise and sloppy. */
      /* Result in qcorr */
      clear_c_array_field(qcorr, NMU);      

      /* Add up the results in qin_diff */
      for(krand = 0; krand < bsfs[nsamp]; krand++)
	sum_c_array_field(qcorr, qin_diff[fine_start+krand], NMU);
      
      /* Average the accumulated results to get bias correction */
      if(bsfs[nsamp] > 0)
	mulreal_c_field(qcorr, 1./((double) bsfs[nsamp]), NMU);
      
      /* Compute average of current density for this block */
      clear_c_array_field(qtmp, NMU);
      for(int kb = 0; kb < bss; kb++)
	sum_c_array_field(qtmp, qin_sloppy[jrand+kb], NMU);
      mulreal_c_field(qtmp, 1./((double) bss), NMU);

      /* Correct the sloppy results for this block by adding the
	 difference between precise and sloppy*/
      sum_c_array_field(qtmp, qcorr, NMU);
      
      /* Compute correlation function for this block */
      msg_tag *tag;
      complex ctmp;
      int rx,ry,rz,rt,d[4],disp=0;
      FORALLSITESLIN(rx,ry,rz,rt,nx,ny,nz,nt){
#ifdef BF
	ctmp = (complex) {0.0,0.0};
	d[XUP]=rx;d[YUP]=ry;d[ZUP]=rz;d[TUP]=rt;
	tag = start_general_gather_field(qtmp,sizeof(complex)*NMU,d,EVENANDODD,gen_pt[0]);
	wait_general_gather(tag);
	FORALLFIELDSITES(i){
	  CSUM(ctmp,cmul_add(qtmp+i*NMU,(complex *) gen_pt[0][i],NMU));
	}
	g_complexsum(&ctmp);
	CDIVREAL(ctmp,volume,ctmp);
	if(this_node == node_number(rx,ry,rz,rt)) outb[node_index(rx,ry,rz,rt)]=ctmp;
#endif

#ifdef CORR_VAR_DISP
	/* Find the statistics of correlation values over the blocks of sites */
	if((rx==1 && ry==1 && rz==0 && rt==0) || (rx==2 && ry==0 && rz==0 && rt==0) || (rx==2 && ry==2 && rz==0 && rt==0)){
#ifndef BF
	  ctmp = (complex) {0.0,0.0};
	  d[XUP]=rx;d[YUP]=ry;d[ZUP]=rz;d[TUP]=rt;
	  tag = start_general_gather_field(qtmp,sizeof(complex)*NMU,d,EVENANDODD,gen_pt[0]);
	  wait_general_gather(tag);
#endif
	  complex *j2;
	  FORALLFIELDSITES(i){
	    j2=(complex *) gen_pt[0][i];
	    for(int mu =0; mu<NMU;mu++) corr_avg[disp][i*nblock+ib]+=(Real) qtmp[i*NMU+mu].real*((j2+mu)->real);
	    corr_sq[disp][i*nblock+ib]+=corr_avg[disp][i*nblock+ib]*corr_avg[disp][i*nblock+ib];
	  }
	  disp++;
#ifndef BF
	  cleanup_general_gather(tag);
#endif
	}
#endif
#ifdef BF
	cleanup_general_gather(tag);
#endif
      }

#ifdef FFT
      /* The forward FT is done separately for each current component */
      restrict_fourier_field((complex *)qtmp, NMU*sizeof(complex), FORWARDS);
      
      /* Square and sum over components to get the correlator */
      dot_corr(outf, qtmp, NMU);

      /* Consistency check */
      double qtot = 0.;
      FORALLFIELDSITES(i){
	qtot += out[i].real;
      }
      g_doublesum(&qtot);
      node0_printf("qtot[%d][%d] = %g\n",ib, jrand, qtot);fflush(stdout);
      
      /* The backward FT: out = C(-r) = C(r) */
      restrict_fourier_field(outf, sizeof(complex), BACKWARDS);
#endif
      /* Accumulate the result for this block size */
      sum_field_c2r(qblock[ib], out);

      /* Accumulate the squares of the results */
      sum_sq_field_c2r(q2block[ib], out);

      fine_start+=bsfs[nsamp];
      nsamp++;
    } /* jrand: loop over blocks*/

    /* Compute the mean and the variance of the mean */
    mean_var_r_field(qblock[ib], q2block[ib], nsamp);

#ifdef CORR_VAR_DISP
    /* Accumulate the result of statistics for the given displacements */
    for(int disp=0;disp<n_disp;disp++){
      Real d;
      FORALLFIELDSITES(i){
	corr_avg[disp][i*nblock+ib]/=nsamp;
	corr_sq[disp][i*nblock+ib]/=nsamp;
	d=corr_sq[disp][i*nblock+ib]-corr_avg[disp][i*nblock+ib]*corr_avg[disp][i*nblock+ib];
	if(d<0.) d=0.;
	corr_sq[disp][i*nblock+ib]=d/(nsamp-1);
      }
    }
#endif
  } /* ib: loop over block sizes */

#ifdef DEBUG
  int perm[4]={0,2,2,0},perm0[4]={0,2,2,0};
  int tmp=0,rx,ry,rz,rt,sum;
  double x_inv[nblock],ys[nblock],yds[nblock],y;
  for(int ib=0;ib<nblock;ib++) x_inv[ib]=1/((double) block_size[ib]);
    
  for(i=0; i<6; i++){
    sum=0;
    for(int j=0;j<3;j++) sum+=(perm[j]-perm0[j])*(perm[j]-perm0[j]);
    if(!(perm[i%3]==perm[(i+1)%3] || (i!=0 && sum==0))){
      tmp=perm[(i+1)%3];
      perm[(i+1)%3]=perm[i%3];
      perm[i%3]=tmp;
      rt=perm[3];
      for(int ix=0; ix<(perm[0]?2:1);ix++){
	for(int iy=0; iy<(perm[1]?2:1);iy++){
	  for(int iz=0; iz<(perm[2]?2:1);iz++){
	    y=0;
	    rx=(nx+(ix?-1:1)*perm[0])%nx;
	    ry=(ny+(iy?-1:1)*perm[1])%ny;
	    rz=(nz+(iz?-1:1)*perm[2])%nz;
	    if(this_node == node_number(rx,ry,rz,rt)){
	      for(int ib=0;ib<nblock;ib++){
		ys[ib]=(double) qblock[ib][node_index(rx,ry,rz,rt)];
		yds[ib]=(double) q2block[ib][node_index(rx,ry,rz,rt)];
	      }
	      y=weighted_lin_extp(x_inv,ys,yds,nblock);
	    }
	    g_doublesum(&y);
	    node0_printf("C(%d,%d,%d,%d)= %.16e\n",rx,ry,rz,rt,y);
	  }
	}
      }
    }
  }
#endif

#ifdef CORR_VAR_DISP
  Real *corr=create_r_field();
  for(int disp=0;disp<n_disp;disp++){
    clear_r_field(corr);
    double x=0,c_mean=0,c_var=0,x_inv[nblock],xm;
    for(int ib=0;ib<nblock;ib++) x_inv[ib]=1/((double) block_size[ib]);
    FORSOMEFIELDPARITY(i,EVEN){
      x=weighted_lin_extp(x_inv,corr_avg[disp]+i*nblock,corr_sq[disp]+i*nblock,nblock);
      corr[i]=x;
      c_mean+=x;
      c_var+=x*x;
      node0_printf("cube_size=1&disp=%d) %.16e at node 0\n",disp,x);
      if(this_node != 0) send_field((char*)  &x,sizeof(double),0);
      if(this_node == 0){
	for(int nd=1;nd<number_of_nodes;nd++){
	  get_field((char*) &xm,sizeof(double),nd);
	  printf("cube_size=1&disp=%d) %.16e at node=%d\n",disp,xm,nd);
	}
      }
      g_sync();
    }
    g_doublesum(&c_mean);
    g_doublesum(&c_var);
    node0_printf("cube_size=%d&disp=%d) mean: %.16e sd: %.16e\n", csizes[0],disp,c_mean/volume,sqrt((c_var-c_mean*c_mean/volume)/volume));fflush(stdout);
    for(int ibc=1;ibc<n_csizes;ibc++){
      c_mean=0; c_var=0;
      FORALLSITESLIN(cx,cy,cz,ct,nx/csizes[ibc],ny/csizes[ibc],nz/csizes[ibc],nt/csizes[ibc]){
	x=avg_over_cube(corr,cx,cy,cz,ct,csizes[ibc]);
	c_mean+=x;
	c_var+=x*x;
	node0_printf("cube_size=%d&disp=%d) %.16e at (%d,%d,%d,%d)\n",csizes[ibc],disp,x,cx,cy,cz,ct);
      }
      int ncubes = (int) volume/pow(csizes[ibc],4);
      node0_printf("cube_size=%d&disp=%d) mean: %.16e sd: %.16e\n", csizes[ibc],disp,c_mean/ncubes,sqrt((c_var-c_mean*c_mean/ncubes)/ncubes));fflush(stdout);
    }
  }
  destroy_r_field(corr);
#endif

  destroy_c_field(qtmp);
  destroy_c_field(qcorr);
  destroy_c_field(outb);
  destroy_c_field(outf);
#ifdef CORR_VAR_DISP
  for(int disp=0;disp<n_disp;disp++){
    destroy_r_array_field(corr_avg[disp],nblock);
    destroy_r_array_field(corr_sq[disp],nblock);
  }
#endif

} /* rcorr.c */

/******************************************************************************/
/* Input fields are have NMU values per site (one for each current component)
   There is one such field per random number group */
/* Output field has one real value per time-slice, one such field for each
   blocking size */

void
rcorr_time(Real *qblock[], Real *q2block[], 
	   complex *qin_sloppy[], int nrand_sloppy, 
	   complex *qin_diff[], int nrand_diff,
	   int nblock, int block_size[]){
  complex *qtmp;
  complex *qcorr;
  int i, x, y, z, t, mu, jrand, krand;
  
  /* Create blocked averages of the corrected result, 
     and compute the correlation within each block */
  Real *out = malloc(sizeof(Real)*nt);
  Real *j_t = malloc(sizeof(Real)*nt*NMU);
  if(out == NULL || j_t == NULL){
    node0_printf("Not enough memory space.\n");
    exit(1);
  }


  qtmp = create_c_array_field(NMU);
  qcorr = create_c_array_field(NMU);
  
  for(int ib = 0; ib < nblock; ib++){
    int bss = block_size[ib];        /* Size of one block for sloppy random sources */
    int nblocks = nrand_sloppy/bss;  /* Number of blocks */
    int bsfs[nblocks];               /* Size of each block for fine random sources */
    int nsamp = 0;                   /* Counter for block number */
    int fine_start = 0;              /* The starting index of a block of fine random solves */
    for(int ibf=0; ibf<nblocks; ibf++) bsfs[ibf] = ibf<nrand_diff%nblocks ? nrand_diff/nblocks+1: nrand_diff/nblocks;

    for(jrand = 0; jrand < nrand_sloppy; jrand += bss){
      
      /* Average qin_diff, the differnece between precise and sloppy. */
      /* Result in qcorr */
      clear_c_array_field(qcorr, NMU);      
      
      /* Add up the results in qin_diff */
      for(krand = 0; krand < bsfs[nsamp]; krand++)
	sum_c_array_field(qcorr, qin_diff[fine_start+krand], NMU);
      
      /* Average the accumulated results to get bias correction */
      if(bsfs[nsamp] > 0)
	mulreal_c_field(qcorr, 1./((double) bsfs[nsamp]), NMU);
      
      /* Compute average of current density for this block */
      clear_c_array_field(qtmp, NMU);
      for(int kb = 0; kb < bss; kb++)
	sum_c_array_field(qtmp, qin_sloppy[jrand+kb], NMU);
      mulreal_c_field(qtmp, 1./((double) bss), NMU);
      
      /* Correct the sloppy results for this block by adding the
	 difference between precise and sloppy*/
      sum_c_array_field(qtmp, qcorr, NMU);
      
      /* Compute time-to-time correlation function for this block */
      for(t=0;t<nt*NMU;t++) j_t[t]=0;
      FORALLSITESLIN(x,y,z,t,nx,ny,nz,nt){
	if(this_node == node_number(x,y,z,t)){
	  for(mu=0;mu<NMU;mu++) j_t[t*NMU+mu] += (Real) qtmp[node_index(x,y,z,t)*NMU+mu].real;
	}
      }
      g_vecfloatsum(j_t,nt*NMU);
      for(t=0;t<nt*NMU;t++) j_t[t]/=(Real)nx*ny*nz; //node0_printf("j_t[%d]_0= %f\n",t,j_t[t]);}

      for(t=0;t<nt;t++) out[t]=0;
      for(int dt=0;dt<nt;dt++){
	for(t=0;t<nt;t++)
	  for(mu=0;mu<NMU-1;mu++) out[dt]+=j_t[((t+dt)%nt)*NMU+mu]*j_t[t*NMU+mu];
	  //COMPSUM(mu,MU0,MUF) out[dt]+=j_t[((t+dt)%nt)*NMU+mu]*j_t[t*NMU+mu];
	out[dt]/=nt;
      }
      

      /* Accumulate the result for this block size */
      for(t=0;t<nt;t++){
        qblock[ib][t]+=out[t];
	q2block[ib][t]+=out[t]*out[t];
      }

      fine_start+=bsfs[nsamp];
      nsamp++;
    } /* jrand: loop over blocks*/

    /* Compute the mean and the variance of the mean */
    Real d;
    for(t=0;t<nt;t++){
      if(nsamp>1){
        qblock[ib][t]/=nsamp;
	q2block[ib][t]/=nsamp;
	d = q2block[ib][t] - qblock[ib][t]*qblock[ib][t];
	if(d < 0.)d = 0.;
	q2block[ib][t] = d/(nsamp-1);
      }
      else{
        q2block[ib][t]=0;
      }
    }
  } /* ib: loop over block sizes */

  destroy_c_field(qtmp);
  destroy_c_field(qcorr);
  free(out);
  free(j_t);

}

/******************************************************************************/
/* Input fields are have NMU values per site (one for each current component)
   There is one such field per random number group */
/* Output field has one real value per time-slice, one such field for each
   blocking size */

void
rcorr_t2tfrmpt2pt(Real *qblock[], Real *q2block[], 
      complex *qin_sloppy[], int nrand_sloppy, 
      complex *qin_diff[], int nrand_diff,
      int nblock, int block_size[]){
  complex *qtmp;
  complex *qcorr;
  Real *q = malloc(sizeof(Real)*nt);
  int i, jrand, krand;
  int x, y, z, t, dx, dy, dz, dt, h = 2*RMAX+1;
  int *multp = malloc(sizeof(int)*h*h*h);

#if 0
  FORALLSITESLIN(dx,dy,dz,dt,h,h,h,1){
    multp[dx+h*dy+h*h*dz]=0;
    if(r2(dx-RMAX,dy-RMAX,dz-RMAX,0)<RMAX*RMAX+1) 
      FORALLSITESLIN(x,y,z,t,h,h,h,1){
        if((r2(x-RMAX,y-RMAX,z-RMAX,0)<=RMAX*RMAX) && (r2(x+dx-2*RMAX,y+dy-2*RMAX,z+dz-2*RMAX,0)<=RMAX*RMAX))
          multp[dx+h*dy+h*h*dz]++;
	//node0_printf("%d\n",multp[dx+h*dy+h*h*dz]);fflush(stdout);
      }
  }
#endif
  
  /* Create blocked averages of the corrected result, 
     and compute the correlation within each block */
  complex *out = create_c_field();
  qtmp  = create_c_array_field(NMU);
  qcorr = create_c_array_field(NMU);

  for(int ib = 0; ib < nblock; ib++){
    clear_r_field(qblock[ib]);
    clear_r_field(q2block[ib]);
    int bss = block_size[ib];        /* Size of one block for sloppy random sources */
    int nblocks = nrand_sloppy/bss;  /* Number of blocks */
    int bsfs[nblocks];               /* Size of each block for fine random sources */
    int nsamp = 0;                   /* Counter for block number */
    int fine_start = 0;              /* The starting index of a block of fine random solves */
    for(int ibf=0; ibf<nblocks; ibf++) bsfs[ibf] = ibf<nrand_diff%nblocks ? nrand_diff/nblocks+1: nrand_diff/nblocks;

    for(jrand = 0; jrand < nrand_sloppy; jrand += bss){
      
      /* Average qin_diff, the differnece between precise and sloppy. */
      /* Result in qcorr */
      clear_c_array_field(qcorr, NMU);      

      /* Add up the results in qin_diff */
      for(krand = 0; krand < bsfs[nsamp]; krand++)
	sum_c_array_field(qcorr, qin_diff[fine_start+krand], NMU);
      
      /* Average the accumulated results to get bias correction */
      if(bsfs[nsamp] > 0)
	mulreal_c_field(qcorr, 1./((double) bsfs[nsamp]), NMU);
      
      /* Compute average of current density for this block */
      clear_c_array_field(qtmp, NMU);
      for(int kb = 0; kb < bss; kb++)
	sum_c_array_field(qtmp, qin_sloppy[jrand+kb], NMU);
      mulreal_c_field(qtmp, 1./((double) bss), NMU);

      /* Correct the sloppy results for this block by adding the
	 difference between precise and sloppy*/
      sum_c_array_field(qtmp, qcorr, NMU);
      
      /* Compute correlation function for this block */

      /* The forward FT is done separately for each current component */
      restrict_fourier_field((complex *)qtmp, NMU*sizeof(complex), FORWARDS);
      
      /* Square and sum over components to get the correlator */
      dot_corr(out, qtmp, NMU);

      /* Consistency check */
      double qtot = 0.;
      //int i;
      FORALLFIELDSITES(i){
	qtot += out[i].real;
      }
      g_doublesum(&qtot);
      node0_printf("qtot[%d][%d] = %g\n",ib, jrand, qtot);fflush(stdout);
      
      /* The backward FT: out = C(-r) = C(r) */
      restrict_fourier_field(out, sizeof(complex), BACKWARDS);

      /* Compute the restricted sum of time-slice correlator for this block */
      for(t=0;t<nt;t++) q[t]=0;
      FORALLSITESLIN(dx,dy,dz,t,h,h,h,nt){
        if(this_node==node_number((dx-RMAX+nx)%nx,(dy-RMAX+ny)%ny,(dz-RMAX+nz)%nz,t))
          //q[t]+=multp[dx+h*dy+h*h*dz]*out[node_index((dx-RMAX+nx)%nx,(dy-RMAX+ny)%ny,(dz-RMAX+nz)%nz,t)].real/pow(RMAX,6);
	  q[t]+=out[node_index((dx-RMAX+nx)%nx,(dy-RMAX+ny)%ny,(dz-RMAX+nz)%nz,t)].real/(RMAX*RMAX*RMAX);
      }
      g_vecfloatsum(q,nt);

      for(t=0;t<nt;t++){
	/* Accumulate the result for this block size */
	qblock[ib][t]+=q[t];
	/* Accumulate the squares of the results */
	q2block[ib][t]+=q[t]*q[t];
      }

      fine_start+=bsfs[nsamp];
      nsamp++;
    } /* jrand: loop over blocks*/

    /* Compute the mean and the variance of the mean */
    Real d;
    for(t=0;t<nt;t++){
      if(nsamp>1){
        qblock[ib][t]/=nsamp;
        q2block[ib][t]/=nsamp;
        d = q2block[ib][t] - qblock[ib][t]*qblock[ib][t];
        if(d < 0.)d = 0.;
        q2block[ib][t] = d/(nsamp-1);
      }
      else{
        q2block[ib][t]=0;
      }
    }

  } /* ib: loop over block sizes */

  free(multp); multp=NULL;
  free(q); q=NULL;
  destroy_c_field(qtmp);
  destroy_c_field(qcorr);
  destroy_c_field(out);

} /* rcorr.c */

/******************************************************************************/
/* Input fields are have NMU values per site (one for each current component)
   There is one such field per random number group */
/* This function computes the extrapolated current densities to infinite 
   blocking size and print the reuslts into the file */

void
compute_current_density_and_print(complex *qin_sloppy[], int nrand_sloppy,
				  complex *qin_diff[], int nrand_diff,
				  char *filename ){

  char myname[] = "compute_current_density_and_print";

  int jrand,i;
  FILE *fp;

  fp = fopen(filename,"w");
  if((fp) == NULL){
    node0_printf("%s: Failed to open %s\n", myname, filename);
    exit(1);
  }

  /* Average qin_diff, the differnece between precise and sloppy. */
  /* Result in qcurr */
  complex *qcurr = create_c_array_field(NMU);

  /* Add up the results in qin_diff */
  for(jrand = 0; jrand < nrand_diff; jrand++)
    sum_c_array_field(qcurr, qin_diff[jrand], NMU);

  /* Average the accumulated results and multiply it by a factor nrand_sloppy for later calculation */
  if(nrand_diff > 0)
    mulreal_c_field(qcurr, ((double) nrand_sloppy)/((double) nrand_diff), NMU);
  
  /* Add up the results in qin_sloppy */
  for(jrand = 0; jrand < nrand_sloppy; jrand++)
    sum_c_array_field(qcurr, qin_sloppy[jrand], NMU);

  /* Average the accumulated results */
  if(nrand_sloppy > 0)
    mulreal_c_field(qcurr, 1./((double) nrand_sloppy), NMU);

  int x,y,z,t;
  double msg[NMU];
  FORALLSITESLIN(x,y,z,t,nx,ny,nz,nt){
    if((x+y+z+t)%2==0){
      for(int mu=0;mu<NMU;mu++) msg[mu]=0;
      if(this_node==node_number(x,y,z,t))
	for(int mu=0;mu<NMU;mu++) msg[mu]= (double) qcurr[NMU*node_index(x,y,z,t)+mu].real;
      g_vecdoublesum(msg,NMU);
      if(this_node==0) for(int mu=0;mu<NMU;mu++) fprintf(fp,"j(%d,%d,%d,%d)_%d= %.16e\n",x,y,z,t,mu,msg[mu]);
      g_sync();
    }
  }

  if(this_node==0 && fclose(fp) == EOF) node0_printf("closing vcj file, %s, failed\n",filename);
}
