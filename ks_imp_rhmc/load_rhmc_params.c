/********************** load_rhmc_params.c ****************************/
/* MIMD version 7 */
/* Read RHMC parameters from file */
/* Created 10/24/2006 C. DeTar
 */

#include "ks_imp_includes.h"
#define IF_OK if(status==0)

static int read_broadcast_ratfunc(FILE *fp, 
				  char tag[], params_ratfunc *rf)
{
  int status = 0;
  int i;
  int prompt = 0;  /* We don't use prompts for this file */
  char prompt_order[16];
  char prompt_res[16];
  char prompt_pole[16];
  char prompt_y[16];
  char prompt_z[16];
  char prompt_m[16];

  /* Build prompt strings based on tags */
  sprintf(prompt_y,"y_%s",tag);
  sprintf(prompt_z,"z_%s",tag);
  sprintf(prompt_m,"m_%s",tag);
  sprintf(prompt_order,"order_%s",tag);
  sprintf(prompt_res,"res_%s",tag);
  sprintf(prompt_pole,"pole_%s",tag);

  /* Read and broadcast fractional powers and masses for this term */
  if(mynode() == 0){
    IF_OK status += get_vi(fp, prompt, prompt_y, rf->y, NMASS);
    IF_OK status += get_vi(fp, prompt, prompt_z, rf->z, NMASS);
    IF_OK status += get_vf(fp, prompt, prompt_m, rf->m, NMASS);
  }

  broadcast_bytes((char *)rf->y, sizeof(Real)*NMASS);
  broadcast_bytes((char *)rf->z, sizeof(Real)*NMASS);
  broadcast_bytes((char *)rf->m, sizeof(Real)*NMASS);

  /* Node 0 gets order and broadcasts */
  if(mynode()==0){
    rf->order = 0;
    IF_OK status += get_i(fp, prompt, prompt_order,&(rf->order));
  }

  /* Bail if a problem thus far */
  broadcast_bytes((char *)&status, sizeof(int));
  if(status != 0)return status;

  broadcast_bytes((char *)&(rf->order),sizeof(int));
  
  /* Print interpretation of parameters */
  if(mynode() == 0){
    printf("Loading order %d rational function approximation for %s:\n",
	   rf->order, tag);
    printf("f(x) = (x+4*%f^2)^(%d/%d) (x+4*%f^2)^(%d/%d)\n",
	   rf->m[0], rf->y[0], rf->z[0], rf->m[1], rf->y[1], rf->z[1]);
    printf("       (x+4*%f^2)^(%d/%d) (x+4*%f^2)^(%d/%d)\n",
	   rf->m[2], rf->y[2], rf->z[2], rf->m[3], rf->y[3], rf->z[3]);
  }

  /* All nodes create storage for the A and B arrays */
  /* NB: the dimension is one more than the order */
  rf->res = (Real *)malloc(sizeof(Real)*(rf->order+1));
  if(rf->res == NULL)status = 1;
  rf->pole = (Real *)malloc(sizeof(Real)*(rf->order+1));
  if(rf->pole == NULL)status = 1;
  /* Poll status and bail if any node failed */
  g_intsum(&status);
  if(status != 0)return 1;
  
  /* Node 0 reads values and broadcasts them */
  if(mynode()==0){
    IF_OK for(i = 0; i <= rf->order; i++){
      IF_OK status += get_f(fp, prompt, prompt_res,rf->res+i);
    }
    IF_OK for(i = 0; i <= rf->order; i++){
      IF_OK status += get_f(fp, prompt, prompt_pole,rf->pole+i);
    }
  }
  broadcast_bytes((char *)rf->res,sizeof(Real)*(rf->order+1));
  broadcast_bytes((char *)rf->pole,sizeof(Real)*(rf->order+1));

  /* Final status check */
  broadcast_bytes((char *)&status, sizeof(int));
  return status;
}

/* Returns an array of rhmc parameter structures, one for each
   pseudofermion field and the number of such structures */
   
params_rhmc *load_rhmc_params(char filename[], int n_pseudo)
{
  FILE *fp = NULL;
  char myname[] = "load_rhmc_params";
  int status = 0;
  int prompt = 0;   /* We don't do prompting for these parameters */
  int i;
  int my_n_pseudo;
  params_rhmc *p;
  
  /* Node zero starts reading */
  if(mynode() == 0){
    fp = fopen(filename,"r");
    if(fp == NULL){
      printf("%s: Can't open %s\n",myname,filename);
      status = 1;
    }
    
    /* Read the number of pseudofermion fields and check */
    my_n_pseudo = 0;
    IF_OK status += get_i(fp,prompt,"n_pseudo",&my_n_pseudo);
    IF_OK if(my_n_pseudo != n_pseudo){
      printf("%s: Wrong value of n_pseudo.  Expected %d\n",myname,n_pseudo);
      status += 1;
    }
  }

  /* Bail here if error */
  broadcast_bytes((char *)&status, sizeof(int));
  if(status > 0)return NULL;
  
  /* Make space for the parameters */
  p = (params_rhmc *)malloc(n_pseudo*sizeof(params_rhmc));
  if(p == NULL){
    status = 1;
  }
  
  /* Poll nodes for any errors and bail if so */
  g_intsum(&status);  /* We don't have an intsum */
  if(status != 0)return NULL;

  /* Read rational function parameters for each pseudofermion field */
  
  for(i = 0; i < n_pseudo; i++){
    node0_printf("Loading rational function parameters for phi field %d\n",i);
    if(mynode()==0) { 
      IF_OK status += get_f(fp, prompt, "naik_term_epsilon",&p[i].naik_term_epsilon);
    }
    broadcast_bytes((char *)&(p[i].naik_term_epsilon), sizeof(Real));
    IF_OK status += read_broadcast_ratfunc(fp,"MD",&p[i].MD);
    IF_OK status += read_broadcast_ratfunc(fp,"GR",&p[i].GR);
    IF_OK status += read_broadcast_ratfunc(fp,"FA",&p[i].FA);
  }

  /* Poll nodes for any errors */
  g_intsum(&status);  /* We don't have an intsum */
  
  /* Bail out here if there is a problem */
  if(status != 0)return NULL;

  return p;
}

