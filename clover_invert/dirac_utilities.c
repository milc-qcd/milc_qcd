/***************** dirac_utilities.c *************************************/

/* Diagonal clover spectrum procedures */
/* MIMD version 7 */

#include "cl_inv_includes.h"

/*--------------------------------------------------------------------*/
double start_timing(void){
  double dtime;

#ifdef PRTIME
  dtime = -dclock();
#else
  dtime = 0;
#endif
  return dtime;
}

/*--------------------------------------------------------------------*/
void print_timing(double dtime, char *str){

#ifdef PRTIME
  dtime += dclock();
  node0_printf("Time for %s %e\n",str, dtime);  fflush(stdout);
#endif

}

/*--------------------------------------------------------------------*/
spin_wilson_vector *create_swv_field(void){
  spin_wilson_vector *swv;
  
  swv = (spin_wilson_vector *)
    malloc(sizeof(spin_wilson_vector)*sites_on_node);
  
  if(swv == NULL){
    printf("create_swv_field(%d): No room for temporary\n",this_node);
    terminate(1);
  }
  return swv;
}

/*--------------------------------------------------------------------*/
spin_wilson_vector *extract_swv_from_wp(wilson_prop_field wp, int color){
  return wp[color];
}

/*--------------------------------------------------------------------*/
wilson_prop_field create_wp_field(void){
  wilson_prop_field wp;
  int color;
  
  wp = (wilson_prop_field) malloc(4*sizeof(spin_wilson_vector *));
  if(wp == NULL){
    printf("create_wp_field(%d): No room for temporary\n",this_node);
    terminate(1);
  }

  for(color= 0; color < 3; color++)
    wp[color] = create_swv_field();
  
  return wp;
}

/*--------------------------------------------------------------------*/
wilson_vector *create_wv_field(void){
  wilson_vector *wv;
  
  wv = (wilson_vector *)
    malloc(sizeof(wilson_vector)*sites_on_node);
  
  if(wv == NULL){
    printf("create_wv_field(%d): No room for temporary\n",this_node);
    terminate(1);
  }
  return wv;
}

/*--------------------------------------------------------------------*/
void copy_wv_from_wp(wilson_vector *wv, wilson_prop_field wp, 
		     int color, int spin){
  int i;
  site *s;
  
  FORALLSITES(i,s){
    wv[i] = wp[color][i].d[spin];
  }
}

/*--------------------------------------------------------------------*/
void copy_wp_from_wv(wilson_prop_field wp, wilson_vector *wv, 
		     int color, int spin){
  int i;
  site *s;
  
  FORALLSITES(i,s){
    wp[color][i].d[spin] = wv[i];
  }
}

/*--------------------------------------------------------------------*/
void destroy_swv_field(spin_wilson_vector *swv){
  free(swv);
}

/*--------------------------------------------------------------------*/
void destroy_wv_field(wilson_vector *wv){
  free(wv);
}

/*--------------------------------------------------------------------*/
void destroy_wp_field(wilson_prop_field wp){
  int color;

  if(wp == NULL)return;
  for(color = 0; color < 3; color++){
    destroy_swv_field(wp[color]);
  }
  free(wp);
}
