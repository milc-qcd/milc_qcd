/***************** control_w.c *****************************************/

/* Main procedure for quenched SU3 Wilson fermions 			*/
/* MIMD version 4 */

/* This version computes propagators for Wilson fermions on a
 supplied background field config */

/* Modifications ...

   8/15/96 Made scratch file name a variable C.D. 
   8/10/96 Installed new propagator IO and added timing C.D. 

   2/03/97 Use multiple mass algorithm. MAX_KAP == 6 is used 
           look for MAX_KAP in comments */

#define CONTROL
#include <stdio.h>
#include <math.h>
#include <complex.h>
#include <su3.h>
#include <io_lat.h>
#include <io_wprop.h>
#include "lattice_w.h" 	/* global variables for lattice fields */
#include <comdefs.h>	/* definitions and variables for communications */
#include <string.h>
#define POINT 1
#define WALL 2

/* Comment these out if you want to suppress detailed timing */
#define IOTIME
#define PRTIME

EXTERN  gauge_header start_lat_hdr;     /* Input gauge field header */

/* Object with two Dirac and two color indices. A given element
   of a "wilson_propagator" is accessed by
   object.c[color1].d[spin1].d[spin2].c[color2].real , etc.
   As alway, "d" denotes a Dirac index and "c" a color index.
   "1" refers to the source, "2" to the sink.
*/

main(argc,argv)
{
int readin();
int meascount,todo;
int prompt;
float avm_iters,avs_iters;

double starttime,endtime,dclock();
double dtime;

int MaxCG;
float size_r[MAX_KAP],RsdCG;

register int i;
register site *s;

int spin,color,j,k,t;
int flag;
int ci,si,sf,cf;
int num_prop;
float space_vol;
int cgilu_w(), setup_w();
void gaugefix();
void w_source(), w_meson(), w_baryon(), w_sink();

#define MAX_AUX 9 /* MAX_KAP * 0.5 + 6 */
                  /* see below, do not change this value */

field_offset dest[MAX_KAP],  /* these are the fields we need for */
             mem[MAX_AUX];     /* the multiple masses */


int key[4];
#ifdef RESTRICT_FFT
#define restrict rstrict /* C-90 T3D cludge */
void restrict_fourier(), setup_restrict_fourier();
int restrict[4];
#else
void fourier(), setup_fourier();
#endif

float norm_fac[10];

static char *mes_kind[10] = {"PION","PS505","PS055","PS0505",
		"RHO33","RHO0303","SCALAR","SCALA0","PV35","B12"};
static char *bar_kind[4] = {"PROTON","PROTON0","DELTA","DELTA0"};

complex *pmes_prop[MAX_KAP][10];
complex *smes_prop[MAX_KAP][10];
complex *bar_prop[MAX_KAP][4];

FILE  *fp_in[MAX_KAP],*fp_out[MAX_KAP];

w_prop_file *fp_in_w[MAX_KAP];        /* For reading binary propagator files */
w_prop_file *fp_out_w[MAX_KAP];       /* For writing binary propagator files */
gauge_file *fp_out_g;

FILE *w_ascii_w_i(); void w_ascii_w(), w_ascii_w_f();
w_prop_file *w_binary_w_i(); void w_binary_w(), w_binary_w_f();
FILE *r_ascii_w_i(); void r_ascii_w(), r_ascii_w_f();
w_prop_file *r_binary_w_i(); void r_binary_w(), r_binary_w_f();

w_prop_file *w_parallel_w_i(); void w_parallel_w(); void w_parallel_w_f();
w_prop_file *w_parallel_w_o(); void w_parallel_w_c();
w_prop_file *r_parallel_w_i(); void r_parallel_w(); void r_parallel_w_f();
w_prop_file *r_parallel_w_o(); void r_parallel_w_c();

void save_ascii(), save_old_binary();

gauge_file *w_parallel_i(char *);
void w_parallel(gauge_file *);
void w_parallel_f(gauge_file *);
gauge_file *w_binary_i(char *);
void w_binary(gauge_file *);
void w_binary_f(gauge_file *);

char scratch_file[MAX_KAP][80];
FILE *fb_scr[MAX_KAP];

FILE *open_w_fast_quarkdump(), *reopen_w_fast_quarkdump(), *open_r_fast_quarkdump();
void fast_quarkdump(), fast_quarkread();
void close_w_fast_quarkdump(), close_r_fast_quarkdump();

    initialize_machine(argc,argv);

    g_sync();
    /* set up */
    prompt = setup_w();
    /* loop over input sets */

    /* Set up Fourier transform */
    key[XUP] = 1;
    key[YUP] = 1;
    key[ZUP] = 1;
    key[TUP] = 0;
#ifdef RESTRICT_FFT
    setup_restrict_fourier(key, restrict);
#else
    setup_fourier(key);
#endif
/* here we assume that MAX_KAP is 6! */
/* note: the vector of mem[0] must not be used by the d_bicgilu_lean.c 
         routine */

    for (k = 0;k < 6;k++)
      {
	dest[k] = F_OFFSET(quark_propagator.c[k/4].d[(k % 4)]);
	mem[k]  = F_OFFSET(quark_propagator.c[(k+6)/4].d[((k+6)%4)]);
      }

/* we got 3 more vectors to assign */

   mem[6] = F_OFFSET(sss);
   mem[7] = F_OFFSET(tmp);
   mem[8] = F_OFFSET(mp);

    while( readin(prompt) == 0)
    {

	starttime=dclock();
	MaxCG=niter;

	avm_iters=0.0;
	meascount=0;

	for(num_prop=0;num_prop<10;num_prop++)
	for(i=0;i<num_kap;i++){
	    pmes_prop[i][num_prop] = (complex *)malloc(nt*sizeof(complex));
	    smes_prop[i][num_prop] = (complex *)malloc(nt*sizeof(complex));
	    for(t=0;t<nt;t++){
		pmes_prop[i][num_prop][t] = cmplx(0.0,0.0); 
		smes_prop[i][num_prop][t] = cmplx(0.0,0.0); 
	    }
	}

	for(num_prop=0;num_prop<4;num_prop++)
	for(i=0;i<num_kap;i++){
	    bar_prop[i][num_prop] = (complex *)malloc(nt*sizeof(complex));
	    for(t=0;t<nt;t++)bar_prop[i][num_prop][t] = cmplx(0.0,0.0); 
	}

	if( fixflag == COULOMB_GAUGE_FIX)
	  {
	    if(this_node == 0) 
	      printf("We are gauge fixing to Coulomb gauge\n");
#ifdef IOTIME
	    dtime = -dclock();
#endif
	    gaugefix(TUP,(float)1.5,500,(float)1.0e-7);
#ifdef IOTIME
	    dtime += dclock();
	    if(this_node==0)printf("Time to gauge fix = %e\n",dtime);
#endif
	    /* Append information to gauge field description */
	    strncpy(start_lat_hdr.gauge_field.descript,
		    " (Coulomb gauge)",
		    sizeof(start_lat_hdr.gauge_field.descript));
	  }
	else
	  if(this_node == 0)printf("COULOMB GAUGE FIXING SKIPPED.\n");

        /* save lattice if requested */
        if( saveflag == SAVE_ASCII){
#ifdef IOTIME
	  dtime = -dclock();
#endif
	  save_ascii(savefile,
		     start_lat_hdr.gauge_field.param[0],
		     start_lat_hdr.gauge_field.param[1]);
#ifdef IOTIME
            dtime += dclock();
            if(this_node==0){
                printf("Time for saving lattice = %e seconds\n",dtime);
            }
            fflush(stdout);
#endif
        }
        else if( saveflag == SAVE_BINARY){
#ifdef IOTIME
            dtime = -dclock();
#endif
	    fp_out_g = w_binary_i(savefile);
	    w_binary(fp_out_g);
	    w_binary_f(fp_out_g);
#ifdef IOTIME
            dtime += dclock();
            if(this_node==0)printf("Time to save = %e\n",dtime);
#endif
        }
	else if( saveflag == SAVE_OLD_BINARY){
#ifdef IOTIME
	  dtime = -dclock();
#endif
	  save_old_binary(savefile,
			  start_lat_hdr.gauge_field.param[0],
			  start_lat_hdr.gauge_field.param[1]);
#ifdef IOTIME
	  dtime += dclock();
	  if(this_node==0)printf("Time to save = %e\n",dtime);
#endif
	  }
        else if( saveflag == SAVE_PARALLEL){
#ifdef IOTIME
            dtime = -dclock();
#endif
	    fp_out_g = w_parallel_i(savefile);
	    w_parallel(fp_out_g);
	    w_parallel_f(fp_out_g);
#ifdef IOTIME
            dtime += dclock();
            if(this_node==0)printf("Time to save = %e\n",dtime);
#endif
        }
	if(this_node==0)printf("END OF HEADER\n");

	if(this_node==0) printf("num_kap = %d\n", num_kap);


	/* Loop over kappas */


	for(k=0;k<num_kap;k++){

	  kappa=kap[k];
	  source_r0=r0[k];
	  RsdCG=resid[k];
	  if(this_node==0)printf("Kappa=%e r0=%e residue=%e\n",
		(double)kappa,(double)source_r0,(double)RsdCG);

	  if (source_r0 != r0[0])
	    if (this_node == 0) 
	      {
		printf("SOURCE SIZE MUST BE EQUAL FOR ALL KAPPA!\n");
		terminate(1);
	      }

	  /* open files for wilson propagators */
	  
	  num_prop=12;	/* Size of one propagator in "wilson_vector" units */
	  
	  if( startflag_w[k] == RELOAD_ASCII)
	    fp_in[k] = r_ascii_w_i(startfile_w[k]);
	  
	  else if( startflag_w[k] == RELOAD_BINARY)
	    fp_in_w[k] = r_binary_w_i(startfile_w[k]);
	  
	  else if( startflag_w[k] == RELOAD_PARALLEL)
	    {
	      /* Read header and close temporarily */
	      fp_in_w[k] = r_parallel_w_i(startfile_w[k]);
	      r_parallel_w_c(fp_in_w[k]);
	    }
	  
	  if( saveflag_w[k] == SAVE_ASCII)
	    fp_out[k] = w_ascii_w_i(savefile_w[k]);
	  
	  else if( saveflag_w[k] == SAVE_BINARY)
	    fp_out_w[k] = w_binary_w_i(savefile_w[k]);
	  
	  else if( saveflag_w[k] == SAVE_PARALLEL)
	    {
	      /* Write header and close temporarily */
	      fp_out_w[k] = w_parallel_w_i(savefile_w[k]);
	      w_parallel_w_c(fp_out_w[k]);
	    }


	  /* have to use scratch file here */

	    sprintf(scratch_file[k],"%s_%02d",scratchstem_w,k);
	    width=r0[0];	/* For write in opening the file */
	    fb_scr[k] = open_w_fast_quarkdump(scratch_file[k],num_prop);
	    close_w_fast_quarkdump(fb_scr[k],scratch_file[k]);
	}

	/* Loop over source colors */
	for(color=0;color<3;color++){

	  /* Loop over source spins */
	  for(spin=0;spin<4;spin++){

	    meascount ++;
	    if(this_node==0)printf("color=%d spin=%d\n",color,spin);
	    w_source(F_OFFSET(chi),color,spin,wallflag,0,0,0,0,r0[0]);

	    flag = 1;      /* Saves one multiplication in cgilu */

	    /* load psi if requested */
	    if( startflag_w[k] == RELOAD_ASCII)
	      { /*r_ascii_w(fp_in[k],spin,color,F_OFFSET(psi));*/
		if(this_node == 0) printf("RELOAD NOT POSSIBLE\n");
	      }
	    else if( startflag_w[k] == RELOAD_BINARY)
	      {
/*
#ifdef IOTIME
			dtime = -dclock();
#endif
			r_binary_w(fp_in_w[k],spin,color,F_OFFSET(psi)); 
#ifdef IOTIME
			dtime += dclock();
			if(this_node==0) printf("Time for serial binary read %e\n",dtime);
#endif
*/
		if(this_node == 0) printf("BINARY RELOAD NOT POSSIBLE\n");
	      }

	    /* site member mp is being used for scratch storage */
	    else if( startflag_w[k] == RELOAD_PARALLEL)
	      {
			/* Reopen, read, and close temporarily */
/*
#ifdef IOTIME
			dtime = -dclock();
#endif
			r_parallel_w_o(fp_in_w[k]);
			r_parallel_w(fp_in_w[k],spin,color,F_OFFSET(psi),F_OFFSET(mp));
			r_parallel_w_c(fp_in_w[k]);
#ifdef IOTIME
			dtime += dclock();
			if(this_node==0) printf("Time for parallel read %e\n",dtime);
#endif
*/
		if (this_node == 0) printf("PARALLEL RELOAD NOT POSSIBLE\n");
	      }
	    
	    avs_iters = (float)cgilu_m_w(F_OFFSET(chi),dest,
				MaxCG,resid,size_r,kap, num_kap,mem,
					 EVENANDODD);
	    avm_iters += avs_iters;
	    for (k=0;k < num_kap;k++)
	      if(this_node==0)
		printf("size_r= %e, iters= %e\n",(double)size_r[k],
		       (double)avs_iters);

	    for (k=0;k < num_kap;k++)
	      {

	    /* Write solutions to scratch disk */
#ifdef IOTIME
		    dtime = -dclock();
#endif
		    fb_scr[k] = reopen_w_fast_quarkdump(scratch_file[k]);
		    fast_quarkdump(fb_scr[k],spin,color,dest[k]);
		    close_w_fast_quarkdump(fb_scr[k],scratch_file[k]);
#ifdef IOTIME
		    dtime += dclock();
		    if(this_node==0) printf("Time for fast_quarkdump %e\n",dtime);
#endif

		    /* Write psi to scratch disk */
		    if( saveflag_w[k] == SAVE_ASCII)
		      {
#ifdef IOTIME
			dtime = -dclock();
#endif
			w_ascii_w(fp_out[k],spin,color,dest[k]);
#ifdef IOTIME
			dtime += dclock();
			if(this_node==0)printf("Time to save a Wilson prop= %e\n",dtime);
#endif
		      }

		    if( saveflag_w[k] == SAVE_BINARY) 
		      {
#ifdef IOTIME
			dtime = -dclock();
#endif
			w_binary_w(fp_out_w[k],spin,color,dest[k]);
#ifdef IOTIME
			dtime += dclock();
			if(this_node==0)printf("Time to save a Wilson prop= %e\n",dtime);
#endif
		      }

		    if( saveflag_w[k] == SAVE_PARALLEL) 
		      {
			/* Reopen, write, and close temporarily */
#ifdef IOTIME
			dtime = -dclock();
#endif
			w_parallel_w_o(fp_out_w[k]);
			w_parallel_w(fp_out_w[k],spin,color,dest[k]);
			w_parallel_w_c(fp_out_w[k]);
#ifdef IOTIME
			dtime += dclock();
			if(this_node==0)printf("Time to save a Wilson prop= %e\n",dtime);
#endif
		      }
		  } /* loop over kappas */

	  } /* source spins */
	} /* source colors */
	
	for (k=0;k < num_kap;k++) 
	  {
	    /* Close scratch files */
	    close_w_fast_quarkdump(fb_scr[k],scratch_file[k]); 
	    if(this_node==0)printf("Saved binary wilson_vector in file  %s\n",
				   scratch_file[k]);
	    
	    /* close files for wilson propagators */
	  if( startflag_w[k] == RELOAD_ASCII) 
	    r_ascii_w_f(fp_in[k],startfile_w[k]);
	    
	    if( startflag_w[k] == RELOAD_BINARY) 
	      r_binary_w_f(fp_in_w[k]); 
	    
	    if( startflag_w[k] == RELOAD_PARALLEL) 
	      r_parallel_w_f(fp_in_w[k]); 
	    
	    if( saveflag_w[k] == SAVE_ASCII) 
	      w_ascii_w_f(fp_out[k],savefile_w[k]);
	    
	    if( saveflag_w[k] == SAVE_BINARY) 
	      w_binary_w_f(fp_out_w[k]); 
	    
	    if( saveflag_w[k] == SAVE_PARALLEL) 
	      w_parallel_w_f(fp_out_w[k]); 
	    
	    
      } /* kappas */
	
	
	for(k=0;k<num_kap;k++) {
	  
	  /* Read the propagator from the scratch file */
	  kappa=kap[k];
	  width=r0[0];
	  fb_scr[k] = open_r_fast_quarkdump(scratch_file[k],num_prop);
	  
#ifdef IOTIME
	    dtime = -dclock();
#endif
	    for(color=0;color<3;color++) for(spin=0;spin<4;spin++){
		fast_quarkread(fb_scr[k], spin, color,
		    F_OFFSET(quark_propagator.c[color].d[spin])); 
	    }

#ifdef IOTIME
	    dtime += dclock();
	    if(this_node==0) 
	      {
		printf("Time to read 12 spin,color combinations %e\n",dtime);
		fflush(stdout);
	      }
#endif

	    close_r_fast_quarkdump(fb_scr[k],scratch_file[k]); 

	    if(this_node==0)printf("Closed fast quark dump file\n");fflush(stdout);

      /*******************************************************************/
      /* FROM HERE ON, EVERYTHING IS THE SAME AS IN control_w.c          */
      /*******************************************************************/

	
    /* spectrum */
#ifdef PRTIME
	    dtime = -dclock();
#endif
	    w_baryon(F_OFFSET(quark_propagator), F_OFFSET(quark_propagator),
		F_OFFSET(quark_propagator), bar_prop[k]);
#ifdef PRTIME
	    dtime += dclock();
	    if(this_node==0)
	      {
		printf("Time for diagonal baryons %e\n",dtime);
		fflush(stdout);
	      }
#endif

#ifdef PRTIME
	    dtime = -dclock();
#endif
	    for(color=0;color<3;color++){
		w_meson(F_OFFSET(quark_propagator.c[color]),
		    F_OFFSET(quark_propagator.c[color]), pmes_prop[k]);
	    }
#ifdef PRTIME
	    dtime += dclock();
	    if(this_node==0)
	      {
		printf("Time for diagonal mesons %e\n",dtime);
		fflush(stdout);
	      }
#endif
#ifdef PRTIME
	    dtime = -dclock();
#endif

	    /* Now convolute the quark propagator with a Gaussian for
	       the smeared mesons. This is done with FFT's */

	    /* fft quark_propagator (in place)--use mp as working space */
	    for(color=0;color<3;color++)for(spin=0;spin<4;spin++){
#ifdef RESTRICT_FFT
		/* Use htmp as a second working space */
		restrict_fourier(F_OFFSET(quark_propagator.c[color].d[spin]),
		    F_OFFSET(mp), F_OFFSET(htmp),
		    sizeof(wilson_vector), FORWARDS);
#else
		fourier(F_OFFSET(quark_propagator.c[color].d[spin]),
		    F_OFFSET(mp), sizeof(wilson_vector), FORWARDS);
#endif
	    }

#ifdef PRTIME
	    dtime += dclock();
	    if(this_node==0)
	      {
		printf("Time for 12 FFT's of quark %d %e\n",k,dtime);
		fflush(stdout);
	      }
#endif
	    /* Use chi, for spin=color=0 for the sink wave function */
	    FORALLSITES(i,s) clear_wvec( &(s->chi));

	    spin=0;color=0;
	    w_sink(F_OFFSET(chi),color,spin,WALL,0,0,0,r0[k]);

	    /* We want chi(-k)* -- the complex conjugate of FFT of the
	       complex conjugate of the quark sink. */
	    FORALLSITES(i,s){
		CONJG(s->chi.d[spin].c[color], s->chi.d[spin].c[color]);
	    }
#ifdef PRTIME
	    dtime = -dclock();
#endif
#ifdef RESTRICT_FFT
	    restrict_fourier(F_OFFSET(chi.d[spin].c[color]),
		F_OFFSET(mp.d[spin].c[color]),
		F_OFFSET(mp.d[1].c[1]), sizeof(complex), FORWARDS);
#else
	    fourier(F_OFFSET(chi.d[spin].c[color]),
		F_OFFSET(mp.d[spin].c[color]), sizeof(complex), FORWARDS);
#endif
	    FORALLSITES(i,s){
		CONJG(s->chi.d[spin].c[color], s->chi.d[spin].c[color]);
	    }

	    /* Now multiply quark by sink wave function */
	    FORALLSITES(i,s)
	    for(ci=0;ci<3;ci++)for(si=0;si<4;si++)
	    for(sf=0;sf<4;sf++)for(cf=0;cf<3;cf++){
		CMUL(s->quark_propagator.c[ci].d[si].d[sf].c[cf],
		    s->chi.d[spin].c[color],
		    s->quark_propagator.c[ci].d[si].d[sf].c[cf]);
	    }

#ifdef PRTIME
	    dtime += dclock();
	    if(this_node==0)
	      {
		printf("Time for FFT of chi and multiply %d %e\n",k,dtime);
		fflush(stdout);
	      }
#endif
	    /* Note: we do not have to FFT the quark propagator back,
	       as long as we don't compute smeared baryons,
	       since the meson construction works the same in momentum space!
	       We do need to devide by an additional (space) volume
	       factor, though! */

	    for(color=0;color<3;color++){
		w_meson(F_OFFSET(quark_propagator.c[color]),
		    F_OFFSET(quark_propagator.c[color]), smes_prop[k]);
	    }

	
	} /* kappas */

	space_vol = (float)(nx*ny*nz);
	for(num_prop=0;num_prop<10;num_prop++) norm_fac[num_prop] = space_vol;
	norm_fac[4] *= 3.0;
	norm_fac[5] *= 3.0;
	norm_fac[8] *= 3.0;
	norm_fac[9] *= 3.0;

	/* print meson propagators */
	for(num_prop=0;num_prop<10;num_prop++)
	for(i=0;i<num_kap;i++){
	    for(t=0; t<nt; t++){
		g_floatsum( &pmes_prop[i][num_prop][t].real );
		pmes_prop[i][num_prop][t].real  /= norm_fac[num_prop];
		g_floatsum( &pmes_prop[i][num_prop][t].imag );
		pmes_prop[i][num_prop][t].imag  /= norm_fac[num_prop];
		if(this_node == 0)
		    printf("POINT%s %d %d  %e %e\n",mes_kind[num_prop],i,t,
			(double)pmes_prop[i][num_prop][t].real,
			(double)pmes_prop[i][num_prop][t].imag);
	    }
	    for(t=0; t<nt; t++){
		g_floatsum( &smes_prop[i][num_prop][t].real );
		smes_prop[i][num_prop][t].real  /=
		    (space_vol*norm_fac[num_prop]);
		g_floatsum( &smes_prop[i][num_prop][t].imag );
		smes_prop[i][num_prop][t].imag  /=
		    (space_vol*norm_fac[num_prop]);
		if(this_node == 0)
		    printf("SMEAR%s %d %d  %e %e\n",mes_kind[num_prop],i,t,
			(double)smes_prop[i][num_prop][t].real,
			(double)smes_prop[i][num_prop][t].imag);
	    }
	}

	/* print baryon propagators */
	for(num_prop=0;num_prop<4;num_prop++)
	for(i=0;i<num_kap;i++){
	    for(t=0; t<nt; t++){
		g_floatsum( &bar_prop[i][num_prop][t].real );
		bar_prop[i][num_prop][t].real  /= space_vol;
		g_floatsum( &bar_prop[i][num_prop][t].imag );
		bar_prop[i][num_prop][t].imag  /= space_vol;
		if(this_node == 0)
		    printf("POINT%s %d %d  %e %e\n",bar_kind[num_prop],i,t,
			(double)bar_prop[i][num_prop][t].real,
			(double)bar_prop[i][num_prop][t].imag);
	    }
	}

	for(num_prop=0;num_prop<10;num_prop++)
	for(i=0;i<num_kap;i++){
	    free(pmes_prop[i][num_prop]);
	    free(smes_prop[i][num_prop]);
	}

	for(num_prop=0;num_prop<4;num_prop++)
	for(i=0;i<num_kap;i++){
	    free(bar_prop[i][num_prop]);
	}

	if(this_node==0)printf("RUNNING COMPLETED\n");
	if(meascount>0){
	    if(this_node==0)printf("total cg iters for measurement= %e\n",
		(double)avm_iters);
	    if(this_node==0)printf("cg iters for measurement= %e\n",
		(double)avm_iters/(double)meascount);
	}

	endtime=dclock();
	if(this_node==0){
	    printf("Time = %e seconds\n",(double)(endtime-starttime));
	    printf("total_iters = %d\n",total_iters);
	}
	fflush(stdout);

    }
}
#undef POINT
#undef WALL
