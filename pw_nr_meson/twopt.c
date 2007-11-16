/***************** twopt.c ****************************************/

/* Compute the meson propagators */
/* MIMD version 7 */

#include "pw_nr_meson_includes.h"
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#define MAXT 96

static complex prop_smear[MAXT][12];
static complex raw_prop[MAXT][3][3][4][4];

static complex meson_nr(pauli_propagator * antiquark, 
			pauli_propagator * quark, 
			int i1, int i2, int j1, int j2)
{
  
  complex x, y, z;
  
  complex result;
  int c0,c1;

  result.real = 0.;
  result.imag = 0.;

  for(c1=0; c1<3; c1++) for(c0=0; c0<3; c0++){
    x = antiquark -> c[c1].d[i1].d[j1].c[c0];
    y = quark -> c[c1].d[i2].d[j2].c[c0];
    CMULJ_(x,y,z);
    /* debug */
    result.real += z.real;
    result.imag += z.imag;
  }
  return(result);
}


/* Caution: elements of mes are not set initially to zero here! 
   In mes averages of propagators are accumulated!  */

static void pw_nr_meson(  pauli_propagator *q1, pauli_propagator *q2 , complex mes[4][4])
{

  complex m_00_x_00, m_00_x_11, m_01_x_01, m_10_x_10, m_01_x_10,
    m_10_x_01, m_11_x_00, m_11_x_11, m_00_x_10, m_01_x_00, m_10_x_00, m_11_x_10,
    m_11_x_01, m_01_x_11, m_10_x_11, m_00_x_01;
  int i1, i2, j1, j2;

/* make only the upper components of the pion and the rho's */

/* 00_x_00  */
i1 = 0;   i2 = 0;    j1 = 0;   j2 = 0;
m_00_x_00 = meson_nr(q1, q2, i1, i2, j1, j2);

/* 00_x_11  */
i1 = 0;   i2 = 0;    j1 = 1;   j2 = 1;
m_00_x_11 = meson_nr(q1, q2, i1, i2, j1, j2);

/* 01_x_01  */
i1 = 0;   i2 = 1;    j1 = 0;   j2 = 1;
m_01_x_01 = meson_nr(q1, q2, i1, i2, j1, j2);

/* 10_x_10  */
i1 = 1;   i2 = 0;    j1 = 1;   j2 = 0;
m_10_x_10 = meson_nr(q1, q2, i1, i2, j1, j2);

/* 01_x_10  */
i1 = 0;   i2 = 1;    j1 = 1;   j2 = 0;
m_01_x_10 = meson_nr(q1, q2, i1, i2, j1, j2);

/* 10_x_01  */
i1 = 1;   i2 = 0;    j1 = 0;   j2 = 1;
m_10_x_01 = meson_nr(q1, q2, i1, i2, j1, j2);

/* 11_x_00  */
i1 = 1;   i2 = 1;    j1 = 0;   j2 = 0;
m_11_x_00 = meson_nr(q1, q2, i1, i2, j1, j2);

/* 11_x_11  */
i1 = 1;   i2 = 1;    j1 = 1;   j2 = 1;
m_11_x_11 = meson_nr(q1, q2, i1, i2, j1, j2);


// below are the mes[sigma_i-source transposed] [sigma_j-sink] combinations. The sigma's are the canopy ones.
//sigma=0 is I

/* A = antiquark; Q = quark */
/* m_i1j1_x_i2j2 = A^*_i1i2 Q_j1j2 */
/* j = src; k = snk */
/* mes[j][k] = Tr(A^* sig_k Q^T sig_j) = sum m_i1j1_x_i2j2 c^jk_i1j1_x_i2j2 */


mes[0][0].real += m_00_x_00.real + m_11_x_11.real + m_11_x_00.real 
                     + m_00_x_11.real;
mes[0][0].imag += m_00_x_00.imag + m_11_x_11.imag + m_11_x_00.imag
                     + m_00_x_11.imag;


mes [1][1].real += m_10_x_10.real + m_01_x_01.real + m_10_x_01.real
                     + m_01_x_10.real;
mes [1][1].imag += m_10_x_10.imag + m_01_x_01.imag + m_10_x_01.imag
                     + m_01_x_10.imag;


mes [3][3].real += m_00_x_00.real + m_11_x_11.real - m_11_x_00.real 
                     - m_00_x_11.real;
mes [3][3].imag += m_00_x_00.imag + m_11_x_11.imag - m_11_x_00.imag
                     - m_00_x_11.imag;


 mes [2][2].real +=  m_10_x_10.real + m_01_x_01.real - m_10_x_01.real 
                     - m_01_x_10.real;
 mes [2][2].imag +=  m_10_x_10.imag + m_01_x_01.imag - m_10_x_01.imag
                     - m_01_x_10.imag;

/*  make off - diagonal spin combinations here */

/* 00_x_10  */
i1 = 0;   i2 = 0;    j1 = 1;   j2 = 0;
m_00_x_10 = meson_nr(q1, q2, i1, i2, j1, j2);

/* 00_x_01  */
i1 = 0;   i2 = 0;    j1 = 0;   j2 = 1;
m_00_x_01 = meson_nr(q1, q2, i1, i2, j1, j2);

/* 01_x_00  */
i1 = 0;   i2 = 1;    j1 = 0;   j2 = 0;
m_01_x_00 = meson_nr(q1, q2, i1, i2, j1, j2);

/* 10_x_00  */
i1 = 1;   i2 = 0;    j1 = 0;   j2 = 0;
m_10_x_00 = meson_nr(q1, q2, i1, i2, j1, j2);

/* 11_x_10  */
i1 = 1;   i2 = 1;    j1 = 1;   j2 = 0;
m_11_x_10 = meson_nr(q1, q2, i1, i2, j1, j2);

/* 11_x_01  */
i1 = 1;   i2 = 1;    j1 = 0;   j2 = 1;
m_11_x_01 = meson_nr(q1, q2, i1, i2, j1, j2);

/* 01_x_11  */
i1 = 0;   i2 = 1;    j1 = 1;   j2 = 1;
m_01_x_11 = meson_nr(q1, q2, i1, i2, j1, j2);

/* 10_x_11  */
 i1 = 1;   i2 = 0;    j1 = 1;   j2 = 1;
m_10_x_11 = meson_nr(q1, q2, i1, i2, j1, j2);

 mes [1][2].real += -1*(- m_01_x_01.imag + m_01_x_10.imag - m_10_x_01.imag   // sing changed for Im
                       + m_10_x_10.imag);
 mes [1][2].imag +=  -1*( m_01_x_01.real - m_01_x_10.real + m_10_x_01.real
                       - m_10_x_10.real);


 mes [1][3].real +=  m_01_x_00.real - m_01_x_11.real + m_10_x_00.real //sign changed
                       - m_10_x_11.real;
 mes [1][3].imag +=  m_01_x_00.imag - m_01_x_11.imag  + m_10_x_00.imag
                       - m_10_x_11.imag;


 mes [2][1].real +=   -1*(m_01_x_01.imag + m_01_x_10.imag - m_10_x_01.imag // sing changed for Im
                       - m_10_x_10.imag);
 mes [2][1].imag += -(- m_01_x_01.real - m_01_x_10.real + m_10_x_01.real
                       + m_10_x_10.real);


 mes [2][3].real +=   -(m_01_x_00.imag - m_01_x_11.imag - m_10_x_00.imag // sign changed for Im
                       + m_10_x_11.imag);
 mes [2][3].imag += -(- m_01_x_00.real + m_01_x_11.real + m_10_x_00.real
                       - m_10_x_11.real);
 

 mes [3][1].real +=   m_00_x_10.real + m_00_x_01.real - m_11_x_01.real // OK
                       - m_11_x_10.real;
 mes [3][1].imag +=   m_00_x_10.imag + m_00_x_01.imag - m_11_x_01.imag
                       - m_11_x_10.imag;


 mes [3][2].real += -(- m_00_x_01.imag + m_00_x_10.imag + m_11_x_01.imag // sign changed for Im
                      - m_11_x_10.imag);
 mes [3][2].imag +=  -( m_00_x_01.real - m_00_x_10.real - m_11_x_01.real
                       + m_11_x_10.real);

  
}  /* pw_nr_meson */

// in raw_prop dim [4][4] is the sigma structures and [3][3] is the polarization structure (directions)


/* rp[dir_src][dir_snk][pol_src][pol_snk] */
static void assemble_prop(complex rp[3][3][4][4], complex *pw_temp)
{
                                           /* a0 (\chi_c0)   */
  pw_temp[0].real += - rp[0][0][1][1].real - rp[1][1][2][2].real 
                    - rp[2][2][3][3].real
                    - rp[0][2][1][3].real - rp[2][0][3][1].real
                    - rp[1][2][2][3].real - rp[2][1][3][2].real
                    - rp[0][1][1][2].real - rp[1][0][2][1].real;

  pw_temp[0].imag += - rp[0][0][1][1].imag - rp[1][1][2][2].imag 
                    - rp[2][2][3][3].imag
                    - rp[0][2][1][3].imag - rp[2][0][3][1].imag
                    - rp[1][2][2][3].imag - rp[2][1][3][2].imag
                    - rp[0][1][1][2].imag - rp[1][0][2][1].imag;

                                           /*  b1 (h_c)     */
//node0_printf("before b1x %e\n",pw_temp[1].real);
  pw_temp[1].real += - rp[0][0][0][0].real;
  pw_temp[1].imag += - rp[0][0][0][0].imag;
//node0_printf("after b1x %e\n",pw_temp[1].real);
  pw_temp[2].real += - rp[1][1][0][0].real;
  pw_temp[2].imag += - rp[1][1][0][0].imag;

  pw_temp[3].real += - rp[2][2][0][0].real;
  pw_temp[3].imag += - rp[2][2][0][0].imag;

                                           /*  a1_x (\chi_c1) */
  pw_temp[4].real += - rp[2][2][2][2].real - rp[1][1][3][3].real
                    + rp[2][1][2][3].real + rp[1][2][3][2].real;
  pw_temp[4].imag += - rp[2][2][2][2].imag - rp[1][1][3][3].imag
                    + rp[2][1][2][3].imag + rp[1][2][3][2].imag;

                                           /*  a1_y    */
  pw_temp[5].real += - rp[2][2][1][1].real - rp[0][0][3][3].real
                    + rp[2][0][1][3].real + rp[0][2][3][1].real;
  pw_temp[5].imag += - rp[2][2][1][1].imag - rp[0][0][3][3].imag
                    + rp[2][0][1][3].imag + rp[0][2][3][1].imag;

                                           /*  a1_z    */
  pw_temp[6].real += - rp[0][0][2][2].real - rp[1][1][1][1].real
                    + rp[0][1][2][1].real + rp[1][0][1][2].real;
  pw_temp[6].imag += - rp[0][0][2][2].imag - rp[1][1][1][1].imag
                    + rp[0][1][2][1].imag + rp[1][0][1][2].imag;

                                           /*  a2_x (T2) (\chi_c2) */
  pw_temp[7].real += - rp[2][2][2][2].real - rp[1][1][3][3].real
                    - rp[2][1][2][3].real - rp[1][2][3][2].real;
  pw_temp[7].imag += - rp[2][2][2][2].imag - rp[1][1][3][3].imag
                    - rp[2][1][2][3].imag - rp[1][2][3][2].imag;

                                           /*  a2_y (T2)    */
  pw_temp[8].real += - rp[2][2][1][1].real - rp[0][0][3][3].real
                    - rp[2][0][1][3].real - rp[0][2][3][1].real;
  pw_temp[8].imag += - rp[2][2][1][1].imag - rp[0][0][3][3].imag
                    - rp[2][0][1][3].imag - rp[0][2][3][1].imag;

                                           /*  a2_z (T2)   */
  pw_temp[9].real += - rp[0][0][2][2].real - rp[1][1][1][1].real
                    - rp[0][1][2][1].real - rp[1][0][1][2].real;
  pw_temp[9].imag += - rp[0][0][2][2].imag - rp[1][1][1][1].imag
                    - rp[0][1][2][1].imag - rp[1][0][1][2].imag;

                                           /*  a2_x (E2)    */
  pw_temp[10].real += - rp[0][0][1][1].real - rp[1][1][2][2].real
                     + rp[0][1][1][2].real + rp[1][0][2][1].real;
  pw_temp[10].imag += - rp[0][0][1][1].imag - rp[1][1][2][2].imag
                     + rp[0][1][1][2].imag + rp[1][0][2][1].imag;

                                           /*  a2_y (E2)    */
  pw_temp[11].real += - rp[2][2][3][3].real - rp[1][1][2][2].real
                     + rp[1][2][2][3].real + rp[2][1][3][2].real;
  pw_temp[11].imag += - rp[2][2][3][3].imag - rp[1][1][2][2].imag
                     + rp[1][2][2][3].imag + rp[2][1][3][2].imag;

}


static void init_raw_prop(void){

  memset((void *)raw_prop, 0, sizeof(raw_prop));

}


/* Caution: elements of propagator are not set initially to zero here!
   Space summation done here */

void all_pw_prop(int dir1, int dir2)
{

  int i;
  int my_t;

  pauli_propagator *quark;
  pauli_propagator *antiquark;
  site *s;

  FORALLSITES(i,s){
    my_t = s->t;

    antiquark = &antiquark_prop[i].up;
    quark     = &quark_prop_smear[i].up;
    pw_nr_meson(antiquark, quark, raw_prop[my_t][dir1][dir2]);
    
#if 0   /* Debug */
    if(s->t==0&&s->x==0&&s->y==0&&s->z==0){
      int j,k;
      for(j=0;j<4;j++)
	for(k=0;k<4;k++){
	  node0_printf("[%d][%d][%d][%d] %e\n", 
		       dir1,dir2,j,k, raw_prop[my_t][dir1][dir2][j][k].real);
	}      
    }
#endif

    antiquark = &antiquark_prop[i].dn;
    quark     = &quark_prop_smear[i].dn;
    pw_nr_meson(antiquark, quark, raw_prop[my_t][dir1][dir2]);
    
  }
}

void assemble_pw_prop(void){
  int t;

  for(t = 0; t < nt; t++){
    assemble_prop(raw_prop[t], prop_smear[t]);
    g_veccomplexsum( prop_smear[t], 12 );
  }
}

/* calculates averages of propagators over polarizations and prints them out
   Expects that prop[nt][12] has been already summed over nodes !!!!*/

void av_pw_prop_out(void)
{
  complex prop_av[6][nt];
  int n,t;
  FILE * kind[6];
  struct stat fs;

  char *trace_kind[6]=
   {"a0","b1", "a1",
    "a2_t2", "a2_e", "a2"};

  char *trace_filename[6] = 
    {a0_file, b1_file, a1_file, a2_t2_file, a2_e_file, a2_file};

  if(this_node==0){

    for(n=0; n<6; n++) {
      if(stat(trace_filename[n], &fs)){
	/* If the file doesn't exist, open it for writing and
	   write the header */
	kind[n] = fopen(trace_filename[n],"w");
	fprintf(kind[n], "# %s\n",trace_kind[n]);
      } else {
	/* Otherwise, append to the file */
	kind[n] = fopen(trace_filename[n],"a");
      }
    }

    for(n=0; n<6; n++) {
      for(t=0;t<nt;t++){
	prop_av[n][t].real = 0.0;
	prop_av[n][t].imag = 0.0;
      }
    }
    /*    a0      */
    for(t=0;t<nt;t++){
      
      prop_av[0][t].real = prop_smear[t][0].real;
      prop_av[0][t].imag = prop_smear[t][0].imag;
      
      /*    b1      */
      for(n=1; n<4; n++) {
	prop_av[1][t].real += prop_smear[t][n].real;
	prop_av[1][t].imag += prop_smear[t][n].imag;
      }
      
      /*    a1      */
      for(n=4; n<7; n++) {
	prop_av[2][t].real += prop_smear[t][n].real;
	prop_av[2][t].imag += prop_smear[t][n].imag;
      }
      
      /*    a2(T2)  */
      for(n=7; n<10; n++) {
	prop_av[3][t].real += prop_smear[t][n].real;
	prop_av[3][t].imag += prop_smear[t][n].imag;
      }
      /*    a2(E)   */
      for(n=10; n<12; n++) {
	prop_av[4][t].real += prop_smear[t][n].real;
	prop_av[4][t].imag += prop_smear[t][n].imag;
      }
      
      /*    a2      */
      for(n=7; n<12; n++) {
	prop_av[5][t].real += prop_smear[t][n].real;
	prop_av[5][t].imag += prop_smear[t][n].imag;
      }
      
    }

    for(n=0; n<6; n++) {
      fprintf(kind[n],"# configuration %06d done\n",sequence);
      fprintf(kind[n],"# nx ny nz = %d %d %d\n", nx, ny,nz);
      fprintf(kind[n],"# nt = %d\n",nt);
      for(t=0;t<nt;t++)
	fprintf(kind[n],"%.8e %.8e\n",prop_av[n][t].real,prop_av[n][t].imag);
      fclose(kind[n]);
    }
  }
}

void init_pw_prop(void){
  int t, i;

  if(nt > MAXT){
    node0_printf("init_pw_prop: recompile with MAXT > %d\n",nt);
    terminate(1);
  }

  for(t=0;t<nt;t++)
    {
      for(i=0;i<12;i++){
	prop_smear[t][i].real = 0.0;
	prop_smear[t][i].imag = 0.0; 
      }
    }

  init_raw_prop();
}

void print_pw_prop(int ksrc, int ksnk){
  int t, i;
  char *trace_kind[12]=
   {"a0_","b1_x_", "b1_y_", "b1_z_", "a1_x_", "a1_y_", "a1_z_",
    "a2_t2x_","a2_t2y_","a2_t2z_", "a2_ex_", "a2_ey_"}; 

  if(this_node != 0)return;

  printf("BEGIN\n");

  for(i=0;i<12;i++)
    {
      printf("\n\n %s%s_%s\n_________________________________\n",
	     trace_kind[i],source_wf_label[ksrc],
	     sink_wf_label[ksnk]);
      for(t=0;t<nt;t++)
	printf("%d %e %e\n", t,
	       prop_smear[t][i].real, prop_smear[t][i].imag);
    }
  
  printf("\n\n");
}

