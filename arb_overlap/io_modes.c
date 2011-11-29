/************* io_modes.c ******************************/
/* MIMD version 7 */

/* wrappers to do serial i/o for eigenmodes */
#include "arb_ov_includes.h"

#define MAX_MODES 4

void read_eigen(wilson_vector **eigVec,int Nmodes, int in_flag, char *in_file)
{
w_prop_file *fp_in_w[MAX_MODES];       /* For eigenmode files */
char in_file_append[MAX_MODES][MAXFILENAME];

int spin=0,color=0,j,jmax,jmark,status;
register site *s;
register int i;

if(Nmodes > 12*MAX_MODES){
i=12*MAX_MODES;
node0_printf("Can only read %d eigenmodes! You wanted %d\n",i,Nmodes);
exit(1);
}
 jmax=Nmodes/(12) + 1;
node0_printf("Opening %d files to read eigenmodes\n",jmax);

/* create the input files and open them */
for(j=0;j<jmax;j++){
  sprintf(in_file_append[j],"%s%d",in_file,j);
  fp_in_w[j]  = r_open_wprop(in_flag, in_file_append[j]);
}

/* loop over all the modes and read them in  by (12)'s */
  for(j=0;j<Nmodes;j++){
    jmark=j/(12);
    if(12*jmark == j){spin=color=0;} /*reset Wilson counters */

    /** TEST **/
node0_printf("reading mode %d with spin %d color %d file %d\n",j,spin,color,jmark);
#ifdef IOTIME
                    status = reload_wprop_sc_to_site( in_flag, fp_in_w[jmark], 
                                      &wqs, spin, color, F_OFFSET(psi),1);
#else
                    status = reload_wprop_sc_to_site( in_flag, fp_in_w[jmark], 
                                      &wqs, spin, color, F_OFFSET(psi),0);
#endif
                  FORALLSITES(i,s){eigVec[j][i]=s->psi;}
spin++;
if((spin %4) == 0){spin=0;color++;}
  }
  /* and close all the files */
node0_printf("preparing to close files\n");
for(j=0;j<jmax;j++)r_close_wprop(in_flag,fp_in_w[j]);
node0_printf("closed files\n");
fflush(stdout);
}

void write_eigen(wilson_vector **eigVec,int Nmodes, int out_flag,
 char *out_file)
{
w_prop_file *fp_out_w[MAX_MODES];       /* For eigenmode files */
char out_file_append[MAX_MODES][MAXFILENAME];
int spin=0,color=0,j,jmax,jmark;
register site *s;
register int i;

if(Nmodes > 12*MAX_MODES){
i=12*MAX_MODES;
node0_printf("Can only write %d eigenmodes! You wanted %d\n",i,Nmodes);
exit(1);
}


 jmax=Nmodes/(12) + 1;
node0_printf("Opening %d files to write %d eigenmodes\n",jmax,Nmodes);

/* create the input files and open them */
for(j=0;j<jmax;j++){
  sprintf(out_file_append[j],"%s%d",out_file,j);
  fp_out_w[j]  = w_open_wprop(out_flag, out_file_append[j], wqs.type);
}


/* loop over all the modes and read them in  by (12)'s */
  for(j=0;j<Nmodes;j++){
    jmark=j/(12);
    if(12*jmark == j){spin=color=0;} /*reset Wilson counters */

    /** TEST **/
node0_printf("writing mode %d with spin %d color %d file %d\n",j,spin,color,jmark);
                  FORALLSITES(i,s){s->psi=eigVec[j][i];}
#ifdef IOTIME
                    save_wprop_sc_from_site( out_flag,fp_out_w[jmark],
                                    &wqs, spin,color,F_OFFSET(psi),1);
#else
                    save_wprop_sc_from_site( out_flag,fp_out_w[jmark],
                                    &wqs, spin,color,F_OFFSET(psi),0);
#endif
spin++;
if((spin %4) == 0){spin=0;color++;}
  }
          /* close files for wilson propagators */
for(j=0;j<jmax;j++)w_close_wprop(out_flag,fp_out_w[j]);



}



void read_eigenval(double *eigVal,int Nmodes, char *in_file)
{
FILE *fp_in_w;       /* For eigenmode files */
char in_file_append[MAXFILENAME];
char extra[]="eigVal";
Real xxx;

int j;


/* create the input files and open them */

strcpy(in_file_append,in_file);
strcat(in_file_append,extra);

if(this_node==0){
  fp_in_w  = fopen(in_file_append,"r");

/* loop over all the modes and read them in  */
  for(j=0;j<Nmodes;j++){
    fscanf(fp_in_w,"%e",&xxx);
    eigVal[j]=(double)xxx;
    printf("reading eigenval %d %e\n",j,eigVal[j]);
  }
fclose(fp_in_w);
}
fflush(stdout);
}
void write_eigenval(double *eigVal,int Nmodes, char *out_file)
{
FILE *fp_out_w;       /* For eigenmode files */
char out_file_append[MAXFILENAME];
char extra[]="eigVal";

int j;


/* create the input files and open them */

strcpy(out_file_append,out_file);
strcat(out_file_append,extra);

if(this_node==0){
  fp_out_w  = fopen(out_file_append,"w");

/* loop over all the modes and write them out */
  for(j=0;j<Nmodes;j++){
    fprintf(fp_out_w,"%e\n",(double)eigVal[j]);
    printf("writing eigenval %d %e\n",j,eigVal[j]);
  }
fclose(fp_out_w);
}
}
