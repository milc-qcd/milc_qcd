/***************** smear.c *********************************************/

/* MIMD version 7 */

#include "pw_nr_meson_includes.h"
#include <sys/types.h>

static complex **phi_src, **phi_snk;   /* Source and sink wave functions */
static complex *phi;  /* Temporary smearing function */

static Real cartesian_Y(int dir, int x, int y, int z){
  int n[3] = {nx, ny, nz};
  int r[3];
  Real rmag;

  if(x > nx-x) r[0] = x-nx; else r[0] = x;
  if(y > ny-y) r[1] = y-ny; else r[1] = y;
  if(z > nz-z) r[2] = z-nz; else r[2] = z;

  rmag = sqrt(r[0]*r[0] + r[1]*r[1] + r[2]*r[2]);

  if(r[dir] == n[dir]/2 || r[dir] == 0)  /* Reflection symmetric */
    //  if(r[dir] == 0) 
    return 0;
  else
    return r[dir]/rmag;
}

static void copy_swave_to_complex_swave(complex *dst, wilson_vector *src, 
					int spin, int color)
{
  int i;
  site *s;

  FORALLSITES(i,s){
    dst[i] = src[i].d[spin].c[color];
  }
}
  

static void swave_to_pwave(complex *pw, complex *sw, int dir,
			   wilson_quark_source *wqs){
  int i;
  site *s;
  Real xs;
  int rx, ry, rz;

  FORALLSITES(i,s){

    /* Shift according to origin of wave */
    rx = (s->x - wqs->x0 + nx) % nx;
    ry = (s->y - wqs->y0 + ny) % ny;
    rz = (s->z - wqs->z0 + nz) % nz;
    
    xs = cartesian_Y(dir, rx, ry, rz);
    CMULREAL(sw[i], xs, pw[i]);
  }
}

/* Load and allocate space for smearing functions */

void load_smearing(wilson_quark_source source_wqs[], 
		   wilson_quark_source sink_wqs[], int n){
  int k;
  wilson_vector *chi; /* Temporary field for reading source */
  
  /* Space for smearing functions */
  phi_src = (complex **)malloc(sizeof(complex *)*n);
  phi_snk = (complex **)malloc(sizeof(complex *)*n);
  if(phi_src == NULL || phi_snk == NULL){
    printf("load_smearing(%d): Out of memory\n",this_node);
    terminate(1);
  }

  /* Space for temporary smearing functions */
  phi = (complex *)malloc(sizeof(complex)*sites_on_node);
  if(phi == NULL){
    printf("load_smearing(%d): Out of memory\n",this_node);
    terminate(1);
  }

  chi = (wilson_vector *)malloc(sizeof(wilson_vector)*sites_on_node);
  if(chi == NULL){
    printf("load_smearing(%d): Out of memory\n",this_node);
    terminate(1);
  }

  for(k=0; k<n; k++){

    phi_src[k] = (complex *)malloc(sizeof(complex)*sites_on_node);
    phi_snk[k] = (complex *)malloc(sizeof(complex)*sites_on_node);
    if(phi_src[k] == NULL || phi_snk[k] == NULL){
      printf("load_smearing(%d): Out of memory\n",this_node);
      terminate(1);
    }

    /* Load source function to chi */
    w_source_field(chi, &source_wqs[k]);

    /* Save it */
    copy_swave_to_complex_swave(phi_src[k], chi, source_wqs[k].spin,
				source_wqs[k].color);

    /* Switch to taking source from stored field from now on */
    source_wqs[k].type = COMPLEX_FIELD_STORE;
    source_wqs[k].c_src = phi_src[k];

    /* Load sink function */
    w_sink_field(phi_snk[k], &sink_wqs[k]);

    /* Switch to taking sink from stored field from now on */
    sink_wqs[k].type = COMPLEX_FIELD_STORE;
    sink_wqs[k].c_src = phi_snk[k];
  }

  free(chi);
}


void free_smearing(int n){
  int k;

  if(phi_snk != NULL)
    for(k = 0; k < n; k++)
      if(phi_snk[k] != NULL)free(phi_snk[k]);

  if(phi_src != NULL)
    for(k = 0; k < n; k++)
      if(phi_src[k] != NULL)free(phi_src[k]);

  if(phi != NULL)free(phi);

  phi_snk = NULL; phi_src = NULL; phi = NULL;
}

/* Construct smeared source by converting to P-wave */

void make_pwave_source(wilson_quark_source wqs[], int ksrc, int dir){

  /* Make P-wave source in phi*/
  swave_to_pwave(phi, phi_src[ksrc], dir, &wqs[ksrc]);
  
  /* Point source to phi */
  wqs[ksrc].type = COMPLEX_FIELD_STORE;
  wqs[ksrc].c_src = phi;
}
  

/* smear in dir2 quark prop. with source in dir1. the input
   quark_propagator is assumed to be already fft'd */

void smear_quark(int ksnk, int dir1, int dir2){
  int i, color, color1,spin, spin1;
  site *s;
  complex pr_tmp, pr_smear;
  int key[4];
  int slice[4];
  
  
  key[XUP] = 1;
  key[YUP] = 1;
  key[ZUP] = 1;
  key[TUP] = 0;
  
  setup_restrict_fourier(key, slice);

  // add Ylm to preloaded source wave function in chi. Result in w_copy.

  swave_to_pwave(phi, phi_snk[ksnk], dir2, &sink_wqs[ksnk]);

  restrict_fourier_field(phi, sizeof(complex), FORWARDS);
  
  FORALLSITES(i,s){
    for(color=0;color<3;color++)
      for(spin=0;spin<2;spin++)   /* Note 2 Pauli spins here */
	for(color1=0;color1<3;color1++)
	  for(spin1=0;spin1<2;spin1++){
	    
	    pr_tmp = quark_prop[i].up.c[color].d[spin].d[spin1].c[color1];
	    CMUL(pr_tmp,phi[i],pr_smear);
	    quark_prop_smear[i].up.c[color].d[spin].d[spin1].c[color1] =
	      pr_smear;

	    pr_tmp = quark_prop[i].dn.c[color].d[spin].d[spin1].c[color1];
	    CMUL(pr_tmp,phi[i],pr_smear);
	    quark_prop_smear[i].dn.c[color].d[spin].d[spin1].c[color1] =
	      pr_smear;
	  }
  }
}

